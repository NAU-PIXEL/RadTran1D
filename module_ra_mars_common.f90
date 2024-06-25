!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_mars_common

  USE module_wrf_error
  USE module_nrutils, only: linear_interpolate
  use module_mars24
  use netcdf
  use module_model_constants, only: pi2
  use module_planet_utilities, only: read_wrf_profile_file, &
                                check_ncdf, dimension_size
  IMPLICIT NONE

  PRIVATE

public :: newton20, dust_distrib_fixed, dust_distribution, oxford_dust, &
     mcd_mgs, mcd_mgsx2, mcd_viking, mcs_dust, init_mcd_mgs, &
     sw_aerosol_scatter, lw_aerosol_heat, demiss3, &
         inject_mcd_mgs, delta_mcd_mgs, da_mcd_mgs, diagnose_tau_2d, &
         mars_msr_fit, chicagodust, get_lower_pbl_ir_heat, cp_mars, &
         dust_distrib_truly_fixed, forced_opac_profile, fop_init, &
         prescribed_dust_storm_init

  PUBLIC :: dust_tes_limb
  PUBLIC :: dust_tes_limb_init

  PUBLIC :: dust_montabone,dust_montabone_init
  INTEGER, SAVE :: n_ls_tes
  INTEGER, SAVE :: n_times_tes
  INTEGER, SAVE :: n_altitudes_tes
  INTEGER, SAVE :: n_latitudes_tes
  INTEGER, SAVE :: n_longitudes_tes
  INTEGER, SAVE :: n_aerosol_types = 2   ! This does *not* include co2 ice, which is done elsewhere

  REAL, SAVE :: dseason_avg_tes
  REAL, SAVE :: dtime_tes, time_offset_tes
  REAL, SAVE :: dlat_tes
  REAL, SAVE :: dlon_tes

  REAL, SAVE, ALLOCATABLE, DIMENSION(:)           :: seasons_tes
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)           :: local_times_tes
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)           :: altitudes_tes
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)           :: latitudes_tes
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)           :: longitudes_tes
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)           :: ls_tes
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: tes_avg_tau
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: tes_avg_tau_err

!montabone dust
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: montabone_cdod610
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) ::     montabone_longitude, &
                                               montabone_latitude, &
                                               montabone_ls
! imposed dust storms from file
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: storm_cdod610
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) ::     storm_longitude, &
                                               storm_latitude, &
                                               storm_ls

  REAL, SAVE    :: montabone_old_ls = 0.
  INTEGER, SAVE :: montabone_old_year = 0
  INTEGER, SAVE :: montabone_rec_len

!end montabone dust

!forced opacity profile
  INTEGER, SAVE :: fop_n_heights, fop_n_times
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: fop_pres, fop_opac
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)   :: fop_time
!end forced opacity profile

CONTAINS

!-----------------------------------------------------------------------
!+ Function generates specific heat values (cp) to be used in flux     !
!    calculations.  The gas mixture assumes a modern-day composition   !
!    of CO2 and N2.  Coefficients are derived from cp/R vs T data in   !
!    Hilsenrath (1955).                                                !
!-----------------------------------------------------------------------
      FUNCTION cp_mars(t)

#ifdef mpas
      USE mpas_atmphys_constants, only : co2_mixing_ratio, cp
#else
      USE module_model_constants, ONLY: co2_mixing_ratio, cp
#endif

      IMPLICIT NONE

      logical, parameter :: fixed_cp = .true.

      REAL(KIND(0.d0)), INTENT(IN   ) ::       t
      
      real(kind(0.d0)) :: cp_mars

      if(fixed_cp) then
        cp_mars = cp
      else
        cp_mars=(((3.47D-07*t*t*t)-(1.269D-03*t*t)+(1.688*t)+443.1)* &  ! First term is for CO2, second is
              co2_mixing_ratio)+(1059.*(1.-co2_mixing_ratio))         !    for N2.
      endif

      END FUNCTION

!====================================================================
  SUBROUTINE newton20( RTHRATEN, t, znu, pi3d,       &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte)
!--------------------------------------------------------------------
!   newton20 is an old subroutine that calculates the
!   Newtonian damping of temperature.
!   The 20 is in reference to the vertical levels in the GFDL Skyhi
!   setup, and the the subsequent "tuning" of parameters (decay constants,
!   height dependence, etc.) to that vertical structure.
!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------
    INTEGER ,       INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte

    REAL , DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(INOUT) ::  RTHRATEN
    REAL , DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(IN   ) ::  t,pi3d
    REAL , DIMENSION( kms:kme ) :: znu

    INTEGER, PARAMETER :: nlevels = 5
    REAL, PARAMETER :: t_naught = 140.     ! Kelvin
    REAL, PARAMETER :: damp_time = 5.e4    ! seconds
    REAL, PARAMETER :: eta_naught = 9.447012e-06   ! "q" level (like sigma)
    INTEGER :: i, j, k
    REAL :: rnx, tbar, vfac, deltat

    ! Use these variables to choose the type of newtonian damping
    ! If both are false, default is point to zonal average
    LOGICAL :: zonal_to_t0 = .FALSE.
    LOGICAL :: point_to_t0 = .FALSE.

    ! This subroutine is meaningless unless we have full access to all grid
    ! points in the E/W direction.  Check for that now.
    IF ((its /= ids) .OR. (ite /= ide)) THEN
       WRITE ( wrf_err_message , * ) 'module_damping: damptop: (its /= ids) or (ite /= ide)',its,ids,ite,ide
       CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
    END IF
    

    rnx = 1./REAL(ide-ids)  ! Remember that ide has extra point for u wind...
                            ! so that this is (ide-1) - ids + 1 = ide-ids

    DO j=jts,jte  ! We don't care if we are averaging an extra row or not
    DO k=kde-nlevels,kde-1
       vfac = ( (LOG(znu(k))/LOG(eta_naught))**6 )/damp_time
       tbar = SUM(t(ids:ide-1,k,j))*rnx
       DO i=its,ite
          IF (zonal_to_t0) THEN
             deltat = tbar - t_naught
          ELSE IF (point_to_t0) THEN
             deltat = t(i,k,j) - t_naught
          ELSE
             deltat = t(i,k,j) - tbar
          END IF
          RTHRATEN(i,k,j)=RTHRATEN(i,k,j)-vfac*deltat/pi3D(i,k,j)
       END DO
    END DO
    END DO

  END SUBROUTINE newton20

   subroutine diagnose_tau_2d(tau_3d, psurf, tau_2d, &
                                ids,ide, jds,jde, kds,kde,         &
                                ims,ime, jms,jme, kms,kme,         &
                                its,ite, jts,jte, kts,kte          )

     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN  ) ::    &
                                                        tau_3d

     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT  ) ::    &
                                                        tau_2d

     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN  ) ::    &
                                                        psurf

        tau_2d(its:ite,jts:jte) = tau_3d(its:ite,kts,jts:jte) * ( 700. / psurf(its:ite,jts:jte) )
        
        return
   
   end subroutine

!====================================================================
   SUBROUTINE dust_distrib_fixed(pph,pp,opt_depth,a,optdpth_array, &
                                 ids,ide, jds,jde, kds,kde,        &
                                 ims,ime, jms,jme, kms,kme,        &
                                 its,ite, jts,jte, kts,kte        )
!---------------------------------------------------------------------
! Based on formula of Conrath (1975)
! Evaluate optical depth at given pressure levels
! Assume pressure increases with index number.
! The parameter a serves to pick at a height above which
! the dust mixing ratio decreases stongly.
! For a=0.01, this is at about 35 km.
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL,  INTENT(IN   ) :: opt_depth, a
!
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::    &
                                                        optdpth_array

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                   pp, &
                                                                  pph
     ! pph is our full sigma level, pp is our half sigma level
 
     INTEGER:: i,j,k
     REAL, DIMENSION(kts:kte+1) :: opt
     REAL :: tau, rnorm, pbar, pbarinv

     pbar = 758.8782           ! reference pressure level (optical depth 
                               ! refers to a mass of this pressure)
     pbarinv= 1.0/pbar

     ! Array indices are reversed in WRF compared to MM5

     DO j=jts,jte
     DO i=its,ite
        ! kte+1 is the top border (full-eta level) of the model
        opt(kte+1)= 0.0
        DO k= kte, kts, -1
           tau=exp( a*(1.0 -  (pbar/pp(i,k,j))**(70./55.) ) )
           opt(k)= opt(k+1) + tau*( pph(i,k,j)-pph(i,k+1,j) )
        END DO
        rnorm= opt_depth * pph(i,kts,j)*pbarinv / opt(kts)
        DO k= kts, kte+1
           optdpth_array(i,k,j)= opt(k)*rnorm
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE dust_distrib_fixed

   !====================================================================
   SUBROUTINE dust_distrib_truly_fixed(pph,pp,opt_depth,a,optdpth_array, &
                                 ids,ide, jds,jde, kds,kde,        &
                                 ims,ime, jms,jme, kms,kme,        &
                                 its,ite, jts,jte, kts,kte        )

! The "truly" means that opt_depth is always the optical depth at the
! actual surface - this will give weird results in a GCM with topography
! but is useful as a proper control in 1D or flat topo idealized sims
!---------------------------------------------------------------------
! Based on formula of Conrath (1975)
! Evaluate optical depth at given pressure levels
! Assume pressure increases with index number.
! The parameter a serves to pick at a height above which
! the dust mixing ratio decreases stongly.
! For a=0.01, this is at about 35 km.
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL,  INTENT(IN   ) :: opt_depth, a
!
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::    &
                                                        optdpth_array

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                   pp, &
                                                                  pph
     ! pph is our full sigma level, pp is our half sigma level
 
     INTEGER:: i,j,k
     REAL, DIMENSION(kts:kte+1) :: opt
     REAL :: tau, rnorm

     ! assumes pph(i,kts,j) = psfc(i,j)

     DO j=jts,jte
     DO i=its,ite
        ! kte+1 is the top border (full-eta level) of the model
        opt(kte+1)= 0.0
        DO k= kte, kts, -1
           tau=exp( a*(1.0 -  (pph(i,kts,j)/pp(i,k,j))**(70./55.) ) )
           opt(k)= opt(k+1) + tau*( pph(i,k,j)-pph(i,k+1,j) )
        END DO
        rnorm= opt_depth / opt(kts)
        DO k= kts, kte+1
           optdpth_array(i,k,j)= opt(k)*rnorm
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE dust_distrib_truly_fixed

   !====================================================================
   SUBROUTINE forced_opac_profile(p_in, optdpth_array, looper,     &
                                 ids,ide, jds,jde, kds,kde,        &
                                 ims,ime, jms,jme, kms,kme,        &
                                 its,ite, jts,jte, kts,kte        )

! this routine forces the dust opacity to be that read in from a horizontal
! averge profile file, but with time and vertical variation
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     INTEGER, INTENT(IN)  :: looper ! 0 = match time exactly to the second
                                    ! 1 = day looper, match local time to second
                                    ! 2 = year looper, match ls
!
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::    &
                                                        optdpth_array

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 p_in
 
     INTEGER:: i,j,k, it

     ! put code in here to find right it, or instead do bilinear interp

     DO j=jts,jte
     DO i=its,ite
        optdpth_array(i,kts:kte,j)=0.
        DO k= kts, kte+1
           optdpth_array(i,k,j)=linear_interpolate(x0=p_in(i,k,j), &
                   x=fop_pres(:,it),y=fop_opac(:,it),              &
                   log_x=.true.,                                   &
                   out_of_range_use_nearest_edge=.true.)
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE forced_opac_profile

   SUBROUTINE fop_init(filename)

   implicit none

   character(len=*), intent(in) :: filename

   CALL read_wrf_profile_file(filename=filename,       &
                              memo="forced opacity",   &
                              n_times=fop_n_times,     &
                              n_heights=fop_n_heights, &
                              pres=fop_pres,           &
                              opac=fop_opac,           &
                              time=fop_time            )

   END SUBROUTINE fop_init

!====================================================================
!!! LJS
   SUBROUTINE chicagodust(pph,pp,xlat,xlong,                  &
                           dust_prs, dust_vsc,                &
                           opt_depth,optdpth_array,           &
                           ids,ide, jds,jde, kds,kde,         &
                           ims,ime, jms,jme, kms,kme,         &
                           its,ite, jts,jte, kts,kte          )
!---------------------------------------------------------------------
! Scheme used to test cloud formation and heating in past climates
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::         ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, jms:jme ),   INTENT(IN   ) :: xlat, xlong

     REAL, INTENT(IN   ) :: opt_depth, dust_prs, dust_vsc

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::    &
                                                           optdpth_array
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                 pp, pph
     INTEGER:: i,j,k
     REAL, DIMENSION(kts:kte+1) :: opt
     REAL :: tau, rnorm, pref, vscale_fac, optdepth, od_inc

     if((dust_prs > -1.) .and. (dust_vsc > -1.)) then  ! dust_prs and _vsc come from namelist
       pref = dust_prs
       vscale_fac = dust_vsc
     else
     ! if they were set to -1. then use average surface pressure over whole domain
       pref = 0.
       DO j=jds,jde
         DO i=ids,ide
           pref = pref + pph(i,1,j)/(jde-jds+1)/(ide-ids+1)
         ENDDO
       ENDDO
       vscale_fac = 0.03
     endif
     
     DO j=jts,jte
     DO i=its,ite
     
        ! Calculate vertical opacity profile
        opt(kte+1) = 0.0
        DO k=kte,kts,-1
           tau = exp(vscale_fac*(1.0-pref/pp(i,k,j)))
           opt(k) = opt(k+1) + tau*(pph(i,k,j)-pph(i,k+1,j))
        END DO
        optdepth = opt_depth
        
!        ! Alter optical depth (specified in namelist) so it decreases from
!        ! full value at +/-50 deg lat to 1/10 of its value at the poles
!        if (abs(xlat(i,j)) .gt. 50.) then
!          od_inc = (1./900.)*(abs(xlat(i,j))-50.)**2
!          optdepth = opt_depth*(1 - 0.9*od_inc)
!        endif
        
        ! Normalize for surface pressure
        rnorm = optdepth*pph(i,kts,j)/pref/opt(kts)
        DO k=kts,kte+1
           optdpth_array(i,k,j) = max(opt(k)*rnorm,1e-10)
        END DO
        
     END DO
     END DO
     
     RETURN
   END SUBROUTINE chicagodust

!====================================================================
   SUBROUTINE dust_distribution(qst01,qst02,p,optdpth_array,       &
                                ids,ide, jds,jde, kds,kde,         &
                                ims,ime, jms,jme, kms,kme,         &
                                its,ite, jts,jte, kts,kte          )
!---------------------------------------------------------------------
! evaluate optical depth at given pressure levels give the
! chemical tracer array
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                qst01, &
                                                                qst02, &
                                                                    p

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                        optdpth_array

     ! Local variables
     INTEGER:: i,j,k
     REAL :: dust, wt1, wt2

     wt1 = 0.25
     wt2 = 1.0

     ! Remember that the order of K's is reversed between WRF and MM5.
     ! Arrays will be reversed to MM5 order in the actual radiation scheme
     ! so keep them in WRF order here
     DO j=jts,jte
     DO i=its,ite
        optdpth_array(i,kte+1,j)= 0.
        DO k=kte, kts, -1
           dust = wt1*qst01(i,k,j) + wt2*qst02(i,k,j)
           dust = max( 1.e-8, dust ) !CEN 25 MAY 2007 - using this lower lim
           optdpth_array(i,k,j)= optdpth_array(i,k+1,j) + &
                                 dust*(p(i,k,j)-p(i,k+1,j))
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE dust_distribution
!====================================================================
!CEN aug 2011 - adding in rad ac dust the 'oxford' way with only one psize:
   SUBROUTINE oxford_dust(dust1,p,dust_array,g,ddpp,rhod, qext,    &
                                ids,ide, jds,jde, kds,kde,         &
                                ims,ime, jms,jme, kms,kme,         &
                                its,ite, jts,jte, kts,kte          )
!---------------------------------------------------------------------
! evaluate optical depth at given pressure levels give the
! chemical tracer array
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
                                                                dust1, &
                                                                    p

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                           dust_array

     REAL, INTENT(IN) :: g, ddpp, qext, rhod

     ! Local variables
     INTEGER:: i,j,k
     REAL, parameter :: dustmin=1.e-8 !CEN 25 MAY 2007 - using this lower lim
     REAL :: aerosol, qextrhor

     qextrhor=0.75*qext/(rhod*(0.5*ddpp)) !!ddpp is diameter, want radius.

     ! Remember that the order of K's is reversed between WRF and MM5.
     ! Arrays will be reversed to MM5 order in the actual radiation scheme
     ! so keep them in WRF order here
     DO j=jts,jte
     DO i=its,ite
        dust_array(i,kte+1,j)= 0.
        DO k=kte, kts, -1
! aerosol(ig,l)=qextrhor*(pq(ig,l)+1.e-8)*(pplev(ig,l)-pplev(ig,l+1))/g
! where qextrhor=(3./4.)*qext/(1250.*ddpp)
! We want to ignore negative dust amounts in radiation scheme
           aerosol = max(dustmin,qextrhor*dust1(i,k,j)*(p(i,k,j)-p(i,k+1,j))/g)
           dust_array(i,k,j)= dust_array(i,k+1,j) + aerosol
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE oxford_dust
!====================================================================

   subroutine inject_mcd_mgs(optdpth_array,tau_od, &
                      l_s,glat,glon,p,             &
                      ids,ide, jds,jde, kds,kde,   &
                      ims,ime, jms,jme, kms,kme,   &
                      its,ite, jts,jte, kts,kte,   &
                      julian, julday, gmt,         &
                      secs_into_sol, dt            )
!---------------------------------------------------------------------
! calculate dust injection rate necessary to match mcd dust
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: p
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: glat, glon
     REAL, INTENT(IN   ) :: l_s
                                            
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                        optdpth_array

     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT  ) ::    &
                                                        tau_od

     real, intent(in) :: julian, dt
     integer, intent(in) :: julday
     real, intent(in) :: gmt, secs_into_sol
!local
    REAL, dimension( ims:ime, jms:jme) :: delta

    real(kind(0d0)) :: dls_prev, j2000, j2000_ott
    real(kind(0d0)) :: djulian_prev, djulian, tmp
    real :: ls_prev, delta_ls
    real :: pi
    pi = ACOS(-1.)
        !dt in minutes

      !tmp is the number of minutes into a day, negative for the day before.
    tmp = mod((gmt*60. + (secs_into_sol/60.)),1440.) - dt
    !again, julday is the start day, not the truncated current day
!bad code    djulian_prev = mod(gmt*60 + xtime - dt, 1440.)/1440.0d0 + (julday-1)
!okay    djulian_prev = (gmt*60 + xtime - dt)/1440.0d0 + (julday-1)
    djulian_prev = mod((gmt*60 + (secs_into_sol/60.) - dt), 1440.0)/1440.0d0 + (floor(julian))
    if (tmp < 0) then
       !subtract 1 because I have gone to the previous day
       djulian_prev=djulian_prev-1
    endif

    call mars24_j2000_ott_from_Mars_Solar_Date(djulian_prev, j2000_ott)
    call mars24_mars_ls(j2000_ott, dls_prev)
    ls_prev = dls_prev
       delta_ls = l_s - ls_prev
    if(ls_prev .gt. l_s) then
       !we wrapped
       delta_ls = delta_ls + 360.
    endif

    call delta_mcd_mgs(tau_od, delta,l_s,glat,glon,     &
         ids,ide, jds,jde,         &
         ims,ime, jms,jme,         &
         its,ite, jts,jte          )
         
    tau_od = tau_od + delta * delta_ls * pi/180.

      where(tau_od < 0) tau_od = 0.0

    call da_mcd_mgs(optdpth_array, tau_od, l_s, glat, glon,p, &
                      ids,ide, jds,jde, kds,kde,         &
                      ims,ime, jms,jme, kms,kme,         &
                      its,ite, jts,jte, kts,kte          )


  end subroutine inject_mcd_mgs

  !====================================================================
  SUBROUTINE delta_mcd_mgs(tau_od, delta, l_s,glat,glon,     &
       ids,ide, jds,jde,         &
       ims,ime, jms,jme,         &
       its,ite, jts,jte          )
    !---------------------------------------------------------------------
    ! calculate delta dust amount at given pressure levels using formula
    ! given in Mars Climate Database for an "MGS scenario"
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !---------------------------------------------------------------------
    INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, &
                                           ims,ime, jms,jme, &
                                           its,ite, jts,jte

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: glat, glon
    REAL, INTENT(IN   ) :: l_s

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT  ) ::    &
         delta

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT  ) ::    &
      tau_od
         

    ! Local variables
    INTEGER:: i,j,k
    REAL :: pi, gls, zls
    REAL :: taueq, tauS, tauN
    REAL :: topdust, tauref, zp

    real :: delta_zls, delta_taueq,  &
         delta_tauS, delta_tauN, delta_topdust,  &
         delta_tempval, delta_res, delta_zp, delta_tauref

    real :: omin, omin_log, tempval

     pi = ACOS(-1.)
     gls = l_s*pi/180.
     zls = SIN(gls-2.76)
     taueq = 0.2 + (0.5-0.2)*(COS(0.5*(gls-4.363))**14)
     tauS  = 0.1 + (0.5-0.1)*(COS(0.5*(gls-4.363))**14)
     tauN  = 0.1

     delta_zls   = cos(gls - 2.76)
     delta_taueq = -(0.5-0.2)*14*0.5*(cos(0.5*(gls-4.363))**13)*sin(0.5*(gls-4.363))
     delta_tauS  = -(0.5-0.1)*14*0.5*(cos(0.5*(gls-4.363))**13)*sin(0.5*(gls-4.363))
     delta_tauN  = 0.0

     DO j=jts,jte
     DO i=its,ite
        IF (glat(i,j) >= 0.) THEN
           ! Northern Hemisphere
           tauref = tauN + (taueq-tauN)*0.5*(1+TANH(4.5-glat(i,j)*18./pi))
           delta_tauref = delta_tauN + 0.5*(delta_taueq-delta_tauN)*(1+TANH(4.5-glat(i,j)*18./pi))
        ELSE
           ! Southern Hemisphere
           tauref = tauS + (taueq-tauS)*0.5*(1+TANH(4.5+glat(i,j)*18./pi))
           delta_tauref = delta_tauS + 0.5*(delta_taueq-delta_tauS)*(1+TANH(4.5+glat(i,j)*18./pi))
        END IF
        
        delta(i,j) = delta_tauref

     END DO
     END DO

     RETURN
   END SUBROUTINE delta_mcd_mgs


  !====================================================================
  SUBROUTINE init_mcd_mgs(tau_od, l_s,glat,glon,     &
       ids,ide, jds,jde,         &
       ims,ime, jms,jme,         &
       its,ite, jts,jte          )
    !---------------------------------------------------------------------
    ! Evaluate optical depth at given pressure levels using formula
    ! given in Mars Climate Database for an "MGS scenario"
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !---------------------------------------------------------------------
    INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, &
         ims,ime, jms,jme, &
         its,ite, jts,jte

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: glat, glon
    REAL, INTENT(IN   ) :: l_s

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT  ) ::    &
      tau_od
         

    ! Local variables
    INTEGER:: i,j,k
    REAL :: pi, gls, zls
    REAL :: taueq, tauS, tauN
    REAL :: topdust, tauref, zp

    real :: omin, omin_log, tempval

     pi = ACOS(-1.)
     gls = l_s*pi/180.
     zls = SIN(gls-2.76)
     taueq = 0.2 + (0.5-0.2)*(COS(0.5*(gls-4.363))**14)
     tauS  = 0.1 + (0.5-0.1)*(COS(0.5*(gls-4.363))**14)
     tauN  = 0.1

     DO j=jts,jte
     DO i=its,ite

        IF (glat(i,j) >= 0.) THEN
           ! Northern Hemisphere
!     glat is in degrees here
           tauref = tauN + (taueq-tauN)*0.5*(1+TANH(4.5-glat(i,j)*0.1))
        ELSE
           ! Southern Hemisphere
           tauref = tauS + (taueq-tauS)*0.5*(1+TANH(4.5+glat(i,j)*0.1))
        END IF
        tau_od(i,j) = tauref

     END DO
     END DO

     RETURN
END SUBROUTINE init_mcd_mgs


!====================================================================
   SUBROUTINE da_mcd_mgs(optdpth_array,tau_od, l_s,glat,glon,p,     &
                      ids,ide, jds,jde, kds,kde,         &
                      ims,ime, jms,jme, kms,kme,         &
                      its,ite, jts,jte, kts,kte          )
!---------------------------------------------------------------------
! Evaluate optical depth at given pressure levels using formula
! given in Mars Climate Database for an "MGS scenario"
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: p
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: glat, glon
     REAL, INTENT(IN   ) :: l_s

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                        optdpth_array
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN  ) ::    &
                                                        tau_od

     ! Local variables
     INTEGER:: i,j,k
     REAL :: pi, gls, zls
     REAL :: taueq, tauS, tauN
     REAL :: topdust, tauref, zp

     pi = ACOS(-1.)
     gls = l_s*pi/180.
     zls = SIN(gls-2.76)


     DO j=jts,jte
     DO i=its,ite
        optdpth_array(i,kte+1,j) = 0.
        topdust = 60.+18.*zls - (32.+18.*zls)*(SIN(glat(i,j))**4) &
                  - 8.*zls*(SIN(glat(i,j))**5)
        tauref = tau_od(i,j)
        DO k=kts,kte
           zp=(700./p(i,k,j))**(70./topdust)
           optdpth_array(i,k,j)= (tauref/700.) * p(i,k,j) * &
                                 MAX( EXP(.007*(1.-MAX(zp,1.))) , 1.e-3 )
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE da_mcd_mgs


!====================================================================
   SUBROUTINE mcd_mgs(optdpth_array,l_s,glat,glon,p,     &
                      optical_depth_multiplier,          &
                      ids,ide, jds,jde, kds,kde,         &
                      ims,ime, jms,jme, kms,kme,         &
                      its,ite, jts,jte, kts,kte          )
!---------------------------------------------------------------------
! Evaluate optical depth at given pressure levels using formula
! given in Mars Climate Database for an "MGS scenario"
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: p
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: glat, glon
     REAL, INTENT(IN   ) :: optical_depth_multiplier ! coopted 'optical depth' namelist parameter
     REAL, INTENT(IN   ) :: l_s
                                            
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                        optdpth_array

     ! Local variables
     INTEGER:: i,j,k
     REAL :: pi, gls, zls
     REAL :: taueq, tauS, tauN
     REAL :: topdust, tauref, zp

     pi = ACOS(-1.)
     gls = l_s*pi/180.
     zls = SIN(gls-2.76)
     taueq = 0.2 + (0.5-0.2)*(COS(0.5*(gls-4.363))**14)
     tauS  = 0.1 + (0.5-0.1)*(COS(0.5*(gls-4.363))**14)
     tauN  = 0.1

     DO j=jts,jte
     DO i=its,ite
        optdpth_array(i,kte+1,j) = 0.
        topdust = 60.+18.*zls - (32.+18.*zls)*(SIN(glat(i,j))**4) &
                  - 8.*zls*(SIN(glat(i,j))**5)
        IF (glat(i,j) >= 0.) THEN
           ! Northern Hemisphere
           tauref = tauN + (taueq-tauN)*0.5*(1+TANH(4.5-glat(i,j)*18./pi))
        ELSE
           ! Southern Hemisphere
           tauref = tauS + (taueq-tauS)*0.5*(1+TANH(4.5+glat(i,j)*18./pi))
        END IF
        DO k=kts,kte
           zp=(700./p(i,k,j))**(70./topdust)
           optdpth_array(i,k,j)= optical_depth_multiplier * &
                                 (tauref/700.) * p(i,k,j) * &
                                 MAX( EXP(.007*(1.-MAX(zp,1.))) , 1.e-3 )
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE mcd_mgs


!====================================================================
   SUBROUTINE mcd_mgsx2(optdpth_array,l_s,glat,glon,p,   & !CEN May 27 2010 - added 
                      ids,ide, jds,jde, kds,kde,         &
                      ims,ime, jms,jme, kms,kme,         &
                      its,ite, jts,jte, kts,kte          )
!---------------------------------------------------------------------
! Evaluate optical depth at given pressure levels using formula
! given in Mars Climate Database for an "MGS scenario"

! HOWEVER, worried about this definition being for IR (not VIS as assumed)
! am doubling amounts to look at result
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: p
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: glat, glon
     REAL, INTENT(IN   ) :: l_s
                                            
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                        optdpth_array

     ! Local variables
     INTEGER:: i,j,k
     REAL :: pi, gls, zls
     REAL :: taueq, tauS, tauN
     REAL :: topdust, tauref, zp

     pi = ACOS(-1.)
     gls = l_s*pi/180.
     zls = SIN(gls-2.76)
     taueq = 0.2 + (0.5-0.2)*(COS(0.5*(gls-4.363))**14)
     tauS  = 0.1 + (0.5-0.1)*(COS(0.5*(gls-4.363))**14)
     tauN  = 0.1

     DO j=jts,jte
     DO i=its,ite
        optdpth_array(i,kte+1,j) = 0.
        topdust = 60.+18.*zls - (32.+18.*zls)*(SIN(glat(i,j))**4) &
                  - 8.*zls*(SIN(glat(i,j))**5)
        IF (glat(i,j) >= 0.) THEN
           ! Northern Hemisphere
           tauref = tauN + (taueq-tauN)*0.5*(1+TANH(4.5-glat(i,j)*18./pi))
        ELSE
           ! Southern Hemisphere
           tauref = tauS + (taueq-tauS)*0.5*(1+TANH(4.5+glat(i,j)*18./pi))
        END IF
        DO k=kts,kte
           zp=(700./p(i,k,j))**(70./topdust)
           optdpth_array(i,k,j)= 2.*(tauref/700.) * p(i,k,j) * &
                                 MAX( EXP(.007*(1.-MAX(zp,1.))) , 1.e-3 )
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE mcd_mgsx2


!====================================================================
   SUBROUTINE mcs_dust(optdpth_array,l_s,glat,glon,p,     &
                      ids,ide, jds,jde, kds,kde,         &
                      ims,ime, jms,jme, kms,kme,         &
                      its,ite, jts,jte, kts,kte          )
!---------------------------------------------------------------------
! calculate dust amount using a very simplistic "fit" to the MCS
! "High Altitude Tropical Dust Maximum"
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: p
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: glat, glon
     REAL, INTENT(IN   ) :: l_s

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                        optdpth_array

     ! Local variables
     INTEGER:: i,j,k
     REAL :: pi, gls, zls
     REAL :: taueq, tauS, tauN
     REAL :: topdust, zp
     
     real :: pmax, log_pmax, pref, &
          log_pref, frac, inv_frac, scale, &
          q_con, phi, log_pval
     REAL, DIMENSION( ims:ime, jms:jme ) ::  tau_od

     pi = ACOS(-1.)

     gls = l_s*pi/180.
     zls = SIN(gls-2.76)

     call init_mcd_mgs(tau_od, l_s, glat, glon, &
       ids,ide, jds,jde,         &
       ims,ime, jms,jme,         &
       its,ite, jts,jte          )

     pmax = 40.
     log_pmax = log(pmax)
     pref = 700.
     log_pref=log(pref)
     frac = 0.5
     inv_frac=1./frac

     DO j=jts,jte
     DO i=its,ite
        optdpth_array(i,kte+1,j) = 0.
        topdust = 60.+18.*zls - (32.+18.*zls)*(SIN(glat(i,j))**4) &
                  - 8.*zls*(SIN(glat(i,j))**5)
        
        DO k=kts,kte
           zp=(pref/p(i,k,j))**(70./topdust)
           q_con =  MAX(EXP(.007*(1.-MAX(zp,1.))), 1e-3)
           log_pval=log_pmax
           if (p(i,k,j) >  pmax) log_pval = log(p(i,k,j))

           if(p(i,k,j) < pref) then
              phi = pi - (pi/2.) * (log_pmax-log_pval)/(log_pmax-log_pref)
              scale = (1-frac)*cos(phi)**2 + frac
              scale = scale*inv_frac
           else 
              scale = 1.0
           endif
           
           optdpth_array(i,k,j)= (tau_od(i,j)/700.) * p(i,k,j) * scale * q_con
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE mcs_dust

!====================================================================
   SUBROUTINE mcd_viking(optdpth_array,l_s,glat,glon,p,     &
                         ids,ide, jds,jde, kds,kde,         &
                         ims,ime, jms,jme, kms,kme,         &
                         its,ite, jts,jte, kts,kte          )
!---------------------------------------------------------------------
! Evaluate optical depth at given pressure levels using formula
! given in Mars Climate Database for a "Viking scenario"
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: p
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: glat, glon
     REAL, INTENT(IN   ) :: l_s

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT  ) ::    &
                                                        optdpth_array

     ! Local variables
     INTEGER:: i,j,k
     REAL :: pi, gls, zls
     REAL :: topdust, tauref, zp

     pi = ACOS(-1.)
     gls = l_s*pi/180.
     zls = SIN(gls-2.76)
     tauref = 0.7 + 0.3*(COS(gls+1.3962634))

     DO j=jts,jte
     DO i=its,ite
        optdpth_array(i,kte+1,j) = 0.
        topdust = 60.+18.*zls - 22.*(SIN(glat(i,j))**2)
        DO k=kts,kte
           zp=(700./p(i,k,j))**(70./topdust)
           optdpth_array(i,k,j)= (tauref/700.) * p(i,k,j) * &
                                 MAX( EXP(.007*(1.-MAX(zp,1.))) , 1.e-3 )
        END DO
     END DO
     END DO

     RETURN
   END SUBROUTINE mcd_viking

!====================================================================
   SUBROUTINE dust_tes_limb(optdpth_array,cloud_array,            &
                            l_s, xlat, xlong, gmt, secs_into_sol, &
                            radfrq, z_at_w,opt_dpth_mult,         &
                            ids,ide, jds,jde, kds,kde,            &
                            ims,ime, jms,jme, kms,kme,            &
                            its,ite, jts,jte, kts,kte,            &
                            annual_average                        )
!---------------------------------------------------------------------
! Evaluate optical depth at a given altitude from TES limb observation
! database
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z_at_w
     REAL, DIMENSION( ims:ime, jms:jme ),          INTENT(IN   ) :: xlat, xlong
     REAL,                                         INTENT(IN   ) :: l_s
     REAL,                                         INTENT(IN   ) :: gmt
     REAL,                                         INTENT(IN   ) :: secs_into_sol
     REAL,                                         INTENT(IN   ) :: radfrq
     REAL,                                         INTENT(IN   ) :: opt_dpth_mult

     LOGICAL, OPTIONAL,                            INTENT(IN   ) :: &
          annual_average

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::    &
          optdpth_array, cloud_array

     ! Local variables
     LOGICAL :: use_annual_averages, use_daily_averages, use_zonal_averages
     INTEGER :: i,j,k,l
     REAL :: elapsed_minutes, local_time, day
     REAL :: day_tau, night_tau, tau
     REAL :: sine_amplitude, sine_offset, pi

     INTEGER :: ilatitude, ilongitude, ilocaltime, iyear

     REAL, DIMENSION(n_altitudes_tes) :: interpolation_vector

     pi = ACOS(-1.)

     use_daily_averages = .FALSE.
     use_zonal_averages = .FALSE.

     use_annual_averages = .FALSE.
     IF (PRESENT(annual_average)) use_annual_averages = annual_average

     IF (use_annual_averages) THEN

        ! XTIME and RADFRQ should be in fractions of a sol ("minutes") so no 
        ! adjustment is necessary here.  Use RADFRQ*0.5 to calculate
        ! climatology from the middle of the radiation timestep.
        ! GMT is set from "simulation_start_[hour,minute,second]" in
        ! module_domain, and so should correctly reflect start times that
        ! start at times other than midnight.

        elapsed_minutes = GMT + (secs_into_sol/60.) + RADFRQ*0.5

     END IF

     optdpth_array(:,:,:) = 0.
     cloud_array(:,:,:) = 0.

     DO j=jts,jte
        DO i=its,ite

           local_time = GMT + MODULO(elapsed_minutes,1440.)/60. + &
                        XLONG(i,j)/15.
           local_time = MODULO(local_time,24.)
           ! local_time must be a number between 0 and 24 ... MODULO lets
           ! negative numbers slip through.
           IF (local_time < 0.) local_time = local_time+24.

           ! Full trilinear interpolation (with linear extrapolation) in
           ! all 3 dimensions (horizontal and vertical) and a sinusoidal
           ! time of day fit (peak and trough at TES day/night measurements)

           DO k=kts,kte+1
              aerosol_types: DO l=1,n_aerosol_types
                 tau = 0.
                 IF (use_annual_averages) THEN
                    night_tau = &
                         LINEAR_INTERPOLATE(xlong(i,j),     &
                                            xlat(i,j),      &
                                            z_at_w(i,k,j),  &
                                            l_s,            &
                                            longitudes_tes, &
                                            latitudes_tes,  &
                                            altitudes_tes,  &
                                            ls_tes,         &
                                            tes_avg_tau(:,:,:,:,1,l))
                    day_tau = &
                         LINEAR_INTERPOLATE(xlong(i,j),     &
                                            xlat(i,j),      &
                                            z_at_w(i,k,j),  &
                                            l_s,            &
                                            longitudes_tes, &
                                            latitudes_tes,  &
                                            altitudes_tes,  &
                                            ls_tes,         &
                                            tes_avg_tau(:,:,:,:,2,l))
                    IF (use_daily_averages) THEN
                       tau = 0.5*(day_tau + night_tau)
                    ELSE
                       sine_amplitude = 0.5 * (day_tau - night_tau)
                       sine_offset    = 0.5 * (day_tau + night_tau)
                       tau = sine_amplitude * SIN(pi*(local_time-8.)/12.) + &
                            sine_offset
                    END IF
                 END IF
                 SELECT CASE (l)
                 CASE (1) ! Dust aerosols
                    ! The following will only work when all longitudes
                    ! are on one tile
                    ! ADT: I can fix this, but I'm not sure it's worth
                    ! the effort...
                    IF (use_daily_averages.AND.use_zonal_averages) THEN
                       optdpth_array(i,k,j) = optdpth_array(i,k,j) + &
                                              tau/float(ide-ids)
                    ELSE
                       optdpth_array(i,k,j) = tau
                    END IF
                 CASE (2) ! Ice aerosols
                    cloud_array(i,k,j)   = tau
                 END SELECT
              END DO aerosol_types

              IF(k > kts) THEN
                IF(optdpth_array(i,k,j) > optdpth_array(i,k-1,j)) THEN
                  optdpth_array(i,k,j) = optdpth_array(i,k-1,j)
                  ! WRITE(0,*) 'alert: dust opt depth inv',i,j,k,&
                  !            optdpth_array(i,k,j),optdpth_array(i,k-1,j)
                ENDIF
                IF(cloud_array(i,k,j) > cloud_array(i,k-1,j)) THEN
                  cloud_array(i,k,j) = cloud_array(i,k-1,j)
                ENDIF
              ENDIF

           END DO ! altitude loop

        END DO ! lon loop
     END DO    ! lat loop

     WHERE (optdpth_array < 0.) optdpth_array = 0.
     WHERE (cloud_array   < 0.)   cloud_array = 0.

     optdpth_array = optdpth_array * opt_dpth_mult

     RETURN
   END SUBROUTINE dust_tes_limb

!====================================================================
  SUBROUTINE SW_AEROSOL_SCATTER(MU0, ALBEDO, OPTD_ARRAY, HEAT, GRNFLX, &
                                sfc_dir,sfc_diff,kts, kte)
!--------------------------------------------------------------------
!
! Calculate heating in shortwave from dust scattering
!
!--------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------

    INTEGER, INTENT(IN) :: kts,kte
    REAL, INTENT(IN) :: MU0, ALBEDO
    REAL, INTENT(OUT) :: GRNFLX,sfc_dir,sfc_diff
    REAL, DIMENSION( kts:kte+1 ), INTENT(IN   ) :: OPTD_ARRAY 
    REAL, DIMENSION( kts:kte ),   INTENT(OUT  ) :: HEAT

! Local Variables
    INTEGER :: k

    REAL, DIMENSION( kts:kte+1 ) ::                         NETFLX, &
                                                            DIRECT, &
                                                               RUP, &
                                                              RUP1, &
                                                               TDN, &
                                                              RDN1

    REAL, DIMENSION( kts:kte ) ::                            SSCAT, &
                                                             FSCAT, &
                                                              GFAC, &
                                                             DIRCT, &
                                                               TR1, &
                                                               RF1, &
                                                                TR, &
                                                                RF


    REAL :: taup, optdp, seczen, sscatp, gfacp, denom, diffuse
    REAL :: upflx, dnflx
    REAL :: expp, expm, en, lam, unew, alpha, gam



    DO k = kts, kte
!       sscat(k)= 0.9
!       gfac(k) = 0.65
!       fscat(k)= 0.4
       sscat(k) = 0.92
       gfac(k)  = 0.55
       fscat(k) = 0.3
       heat(k) = 0.
    END DO
    grnflx=0.

    IF( mu0 > 0. )  THEN
       direct(kts)= 1.
    ELSE
       direct(kts)= 0.
    END IF

    ! Similarity scaling

    DO k = kts, kte
       sscatp= sscat(k)*( 1.-fscat(k) )/(1.-sscat(k)*fscat(k))      !c-mm  Equal to omega_free in KDM
       gfacp = (gfac(k)-fscat(k))/(1.-fscat(k))                     !c-mm  Equal to asymmetry_free in KDM

       ! Rescale optical depths:
       optdp =(1.-sscat(k)*fscat(k))*optd_array(k+1)                !c-mm  Equal to tau_free in KDM
       taup  =(1.-sscat(k)*fscat(k))*(optd_array(k+1)-optd_array(k))
       seczen= 35./sqrt( 1.224e3*mu0*mu0 + 1. )

       ! Form transmission and reflectivities for the individual layers
       lam  =sqrt(3.*(1.-sscatp)*(1.-sscatp*gfacp))
       unew = 1.5*(1.-sscatp*gfacp)/lam
       denom= 1./( 1. - (lam*mu0)**2 )
       alpha= 0.75*mu0*sscatp*(1.+gfacp*(1.-sscatp))*denom
       gam  = 0.5*sscatp*(1.+3.*gfacp*(1.-sscatp)*mu0*mu0)*denom

       ! tr1 and rf1 are independent of zenith angle
       IF( mu0 > 0. ) THEN
          direct(k+1)= exp( -optdp*seczen )
          dirct (k  )= exp(  -taup*seczen )
          expp= exp( lam*taup )
          expm= exp(-lam*taup )
          en  = ( (unew+1.)**2 * expp ) - ( (unew-1.)**2 * expm )
          tr1(k)= 4.*unew/en
          rf1(k)= (unew+1.)*(unew-1.)*(expp-expm)/en
          tr(k) = (alpha+gam)*tr1(k) + &
                  ( (alpha-gam)*rf1(k)-(alpha+gam-1.) )*dirct(k)
          rf(k) = (alpha-gam)*tr1(k)*dirct(k) + (alpha+gam)*rf1(k) -(alpha-gam)
       ELSE
          direct(k+1)= 0.
          dirct (k  )= 0.
          tr1(k)= 1.
          rf1(k)= 0.
          tr (k)= 1.
          rf (k)= 0.
       END IF
    END DO

    ! Now march upward and downward, combining levels:

    ! dirct = exp(-taup/mu0)! transmission of direct beam through level k only
    ! direct= exp(-optd/mu0)! % transmission of direct beam through top
                            !   k-1 levels

    ! -------------First, do the upward marching--------
    ! Combine levels by marching upwards from the surface

    rup (kte+1)= albedo      ! surface reflectivity : direct beam
    rup1(kte+1)= albedo      ! surface reflectivity : diffuse radiation

    DO k= kte, kts, -1
       denom  = 1./( 1.-rf1(k)*rup1(k+1) ) ! factor for multiple refections
       diffuse= tr(k)-dirct(k) ! total transmission - direct
       rup1(k)= rf1(k) + ( tr1(k)*rup1(k+1)*tr1(k) )*denom
       rup(k) = rf(k)+tr1(k)*(diffuse*rup1(k+1)+dirct(k)*rup(k+1))*denom
    END DO

   ! -----------Now, do the downward marching----------

    tdn (kts)= 1.         ! total transmission through the top layer:
                          !  only direct beam is involved
    rdn1(kts)= 0.         ! reflection of top layer to radiation from 
                          !  below: only diffuse radiation
    tdn (kts+1)= tr(kts)
    rdn1(kts+1)= rf1(kts)

    DO k= kts+1, kte
       diffuse  = tdn(k)-direct(k) ! % total transmission
                                   ! - direct through top k-1 levels
       denom    = 1./( 1.-rdn1(k)*rf1(k)) ! % factor for multiple reflections
       tdn(k+1) = direct(k)*tr(k) + tr1(k)* &
                  ( diffuse + direct(k)*rdn1(k)*rf(k) )*denom
       rdn1(k+1)= rf1(k) + ( tr1(k)*rdn1(k)*tr1(k) )*denom
    END DO
      
    IF( mu0 > 0. )  THEN
       DO k= kts, kte+1
          diffuse= tdn(k) - direct(k)
          denom= 1.0 / ( 1.0 - rdn1(k)*rup(k) )
          upflx=( direct(k)*rup(k)+diffuse*rup1(k) ) * denom
          dnflx=direct(k)+(diffuse+direct(k)*rup(k)*rdn1(k))*denom
          netflx(k)= upflx-dnflx
       END DO

       grnflx= -netflx(kte+1)
       sfc_dir=direct(kte+1)/(direct(kte+1)+diffuse)
       sfc_diff=diffuse/(direct(kte+1)+diffuse)   !c-mm  Diffuse is a scalar, but last iteration through the above loop is for kte+1, so this is OK...
       DO k= kts, kte
          heat(k)=  netflx(k+1)-netflx(k)
       END DO

    ELSE

       grnflx= 0.0
       sfc_dir=0.0
       sfc_diff=0.0
       DO k= kts, kte+1
          heat(k)= 0.0
       END DO

    END IF

    RETURN
  END SUBROUTINE SW_AEROSOL_SCATTER

!====================================================================
  SUBROUTINE LW_AEROSOL_HEAT( DELPI, T, TG, DUST, HR_DUST, SFLX_DUST, &
                              G, CP, kts, kte )
!--------------------------------------------------------------------
!
!     Calculates the heating due to aerosols in the longwave
!     Ordering of arrays in the vertical is level 1 is the top
!     and level K_MAX is the surface.  This is inverse from
!     WRF-standard, but reversing of arrays is handled when calling
!     this subroutine.
!
!--------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------

    INTEGER, INTENT(IN) :: kts, kte
    REAL, DIMENSION(kts:kte),   INTENT(IN) :: delpi, t
    REAL, DIMENSION(kts:kte+1), INTENT(IN) :: dust
    REAL, INTENT(IN) :: tg, g, cp
    REAL, DIMENSION(kts:kte), INTENT(OUT) :: hr_dust
    REAL, INTENT(OUT) :: sflx_dust

    INTEGER :: k,kk
    REAL, DIMENSION(kts:kte+1,kts:kte+1) :: dep
    REAL, DIMENSION(kts:kte+1) :: planckfn, optd, dnflux, upflux, netflux
    REAL :: diffac, gcp, boltz, emiss_fact

    ! Look in K.N.Liou's "Radiation and Cloud Processes in the
    ! Atmosphere" for diffac factor (p. 42)
    ! = 1/cosine(mean emergent angle)
    diffac= 1.66

    gcp=G/CP

    boltz = 5.67E-8
    DO k= kts, kte
       planckfn(k)= boltz * t(k)*t(k)*t(k)*t(k)
    END DO
    planckfn(kte+1) = boltz * tg*tg*tg*tg

    ! Should be 1/2.5, but this works better because the LW radiation
    ! doesn't consider dust opacity in the 15 micron band (which is
    ! difficult to implement).  We have increased our opacity to
    ! reflect the extra emissive capability of the 15 micron band that
    ! we are not calculating.

    emiss_fact= 0.7 

    emiss_fact= emiss_fact * diffac

    DO k = kts, kte+1
       IF (k < kte+1) hr_dust(k)=0.
       DO kk = kts, kte+1
          dep(k,kk)=0.
       END DO
    END DO

    CALL DEMISS3( emiss_fact, dust, dep, kts, kte )

    ! Boundary conditions
    upflux(kte+1)= planckfn(kte+1)
    dnflux(kts)= 0.

    DO k= kts, kte
       upflux(k)= planckfn(kte+1)*dep(k,kte+1)  
       DO kk= k, kte
          upflux(k)= upflux(k) + planckfn(kk)*dep(k,kk) 
       END DO

       dnflux(k+1)= 0.
       DO kk= kts, k
          dnflux(k+1)= dnflux(k+1)-planckfn(kk)*dep(k+1,kk)
       END DO
    END DO

    DO k= kts, kte+1
       netflux(k)= upflux(k)-dnflux(k)
    END DO
    sflx_dust= dnflux(kte+1)

    DO k= kts, kte
       hr_dust(k)= (netflux(k+1)-netflux(k))*delpi(k)*gcp
    ENDDO

    RETURN
  END SUBROUTINE LW_AEROSOL_HEAT

!====================================================================
  SUBROUTINE DEMISS3( EFACT, DUST, X, kts, kte )
!--------------------------------------------------------------------
!
!     Calculates the emissivity as a fucntion of solar optical depth as
!     described in Haberle et al., 1982 (Icarus, v.50, p. 322)
!     Emissivity is returned as a 3d variable (i,k,k) as a sort of
!     "transmission" function between different levels.
!     Ordering of vertical indices is the same as the calling subroutine
!     (i.e., LW_AEROSOL_HEAT)
!
!--------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------

    INTEGER, INTENT(IN) :: kts, kte
    REAL, INTENT(IN) :: efact
    REAL, DIMENSION(kts:kte+1), INTENT(IN) :: dust
    REAL, DIMENSION(kts:kte+1,kts:kte+1), INTENT(OUT) :: x

    INTEGER :: j, k, kk
    INTEGER, PARAMETER :: n=5

    REAL, DIMENSION(kts:kte+1,kts:kte+1) :: emiss
    REAL, DIMENSION(kts:kte+1) :: optd
    REAL, DIMENSION(n) :: c

    c = (/  -2.365e-4,  2.082e-1, 6.409e-2, -4.7e-2, 7.217e-3 /)
    ! A non-zero c(1) gives absorption even in a clear atmosphere ...
    !c = (/        0.0,  2.082e-1, 6.409e-2, -4.7e-2, 7.217e-3 /)

    DO k=kts, kte+1
       optd(k)= efact*dust(k)
    END DO

    ! Take advantage of tau(k+1) > tau(k); no need for absolute value

    DO k=kts, kte+1
       DO kk= kts, k
          x(k,kk)= LOG( 1.0 + (optd(k) - optd(kk)) )
          emiss(k,kk)= c(n)*x(k,kk) + c(n-1)
          DO j= n-2, 1, -1
             emiss(k,kk)= emiss(k,kk)*x(k,kk) + c(j)
          END DO
       END DO
    END DO

    ! Invoke symmetry 
    DO k= kts, kte+1
       DO kk= 1, k
          emiss(kk,k)= emiss(k,kk)
       END DO
    END DO

    DO k=kts,kte+1
       DO kk= kts, kte
          x(k,kk)= emiss(k,kk+1) - emiss(k,kk)
       END DO
       x(k,kte+1)= 1. - emiss(k,kte+1)
    END DO

    RETURN
  END SUBROUTINE DEMISS3

!====================================================================
   SUBROUTINE dust_tes_limb_init(limb_data_rec)
!---------------------------------------------------------------------
! Initialize optical depth arrays for TES limb observation derived
! database
!---------------------------------------------------------------------
     IMPLICIT NONE
!---------------------------------------------------------------------

     INTEGER, INTENT(IN) :: limb_data_rec

     INTEGER :: i, j, k, l, m
     INTEGER :: idum, jdum

     SELECT CASE(limb_data_rec)

       CASE(1)
         write(0,*) 'Forcing the model with TES limb dust opacities (climatology)'
         OPEN(1969,FILE='./Data/radiation/TES_limb_tau.txt',STATUS='old')
     
       CASE(2)
         write(0,*) 'Forcing the model with MCS limb dust opacities (experimental)'
         OPEN(1969,FILE='./Data/radiation/MCS_limb_tau.txt',STATUS='old')

       CASE(3)
         write(0,*) 'Forcing the model with TES limb dust opacities (MY24)'
         OPEN(1969,FILE='./Data/radiation/TES_limb_tau_MY24.txt',STATUS='old')
     
       CASE(4)
         write(0,*) 'Forcing the model with TES limb dust opacities (MY25)'
         OPEN(1969,FILE='./Data/radiation/TES_limb_tau_MY25.txt',STATUS='old')
     
       CASE(5)
         write(0,*) 'Forcing the model with TES limb dust opacities (MY26)'
         OPEN(1969,FILE='./Data/radiation/TES_limb_tau_MY26.txt',STATUS='old')
     
       CASE(6)
         write(0,*) 'Forcing the model with TES limb dust opacities (MY27)'
         OPEN(1969,FILE='./Data/radiation/TES_limb_tau_MY27.txt',STATUS='old')

       CASE(99)
         write(0,*) 'Forcing the model with TES limb dust opacities (MY26 with Ls=310-360 from MY24)'
         OPEN(1969,FILE='./Data/radiation/TES_limb_tau_MY2426blend.txt',STATUS='old')

       CASE DEFAULT
         write(0,*) 'You said you wanted to force with limb aerosol data,'
         write(0,*) 'but you did not give me a valid aerosol file flag:'
         write(0,*) 'limb_data_rec=',limb_data_rec,' valid options: 1-6 or 99.'
         STOP
       
     END SELECT
     
     ! Modifications:
     ! 1) Pad the edges in latitude/longitude space so the 2D horizontal
     ! interpolation will *always* have points around it for "reasonable"
     ! interpolations
     ! 2) Pad the edges in altitude space so the 1D vertical interpolation
     ! will *always* have points around it for "reasonable" interpolations

     READ(1969,*) dseason_avg_tes
     READ(1969,*) n_ls_tes
     ALLOCATE(ls_tes(n_ls_tes+2))
     DO i=2,n_ls_tes+1
        READ(1969,*) idum, ls_tes(i)
     END DO
     ls_tes(1)          = ls_tes(2)           - dseason_avg_tes
     ls_tes(n_ls_tes+2) = ls_tes(n_ls_tes+1) + dseason_avg_tes
     n_ls_tes = n_ls_tes + 2

     READ(1969,*) n_times_tes
     ALLOCATE(local_times_tes(n_times_tes))
     DO i=1,n_times_tes
        READ(1969,*) idum, local_times_tes(i)
     END DO
     time_offset_tes = 4.
     dtime_tes = 24./n_times_tes

     READ(1969,*) n_altitudes_tes
     ALLOCATE(altitudes_tes(n_altitudes_tes+2))
     DO i=2,n_altitudes_tes+1
        READ(1969,*) idum, altitudes_tes(i)
     END DO
     ! Altitudes in file are in km; convert to m for consistency with
     ! WRF model units
     altitudes_tes = altitudes_tes * 1.e3
     altitudes_tes(1)                = -1.e5 !  -10 km, below everything
     altitudes_tes(n_altitudes_tes+2)=  1.e6 ! 1000 km, above everything
     n_altitudes_tes = n_altitudes_tes + 2

     READ(1969,*) n_latitudes_tes
     ALLOCATE(latitudes_tes(n_latitudes_tes+2))
     DO i=2,n_latitudes_tes+1
        READ(1969,*) idum, latitudes_tes(i)
     END DO
     dlat_tes = 180./n_latitudes_tes
     latitudes_tes(1) = latitudes_tes(2) - dlat_tes
     latitudes_tes(n_latitudes_tes+2) = &
          latitudes_tes(n_latitudes_tes+1) + dlat_tes
     n_latitudes_tes = n_latitudes_tes + 2

     READ(1969,*) n_longitudes_tes
     ALLOCATE(longitudes_tes(n_longitudes_tes+2))
     DO i=2,n_longitudes_tes+1
        READ(1969,*) idum, longitudes_tes(i)
     END DO
     dlon_tes = 360./n_longitudes_tes
     longitudes_tes(1) = longitudes_tes(2) - dlon_tes
     longitudes_tes(n_longitudes_tes+2) = &
          longitudes_tes(n_longitudes_tes+1) + dlon_tes
     n_longitudes_tes = n_longitudes_tes + 2

     ALLOCATE(tes_avg_tau(n_longitudes_tes,     &
                          n_latitudes_tes,      &
                          n_altitudes_tes,      &
                          n_ls_tes,             &
                          n_times_tes,          &
                          n_aerosol_types    ), & ! 1=dust, 2=ice
              tes_avg_tau_err(n_longitudes_tes, &
                          n_latitudes_tes,      &
                          n_altitudes_tes,      &
                          n_ls_tes,             &
                          n_times_tes,          &
                          n_aerosol_types    ) )

     DO m = 2, n_ls_tes-1
        DO l = 1, n_times_tes
           DO k = 2, n_altitudes_tes-1
              DO j = 2, n_latitudes_tes-1
                 DO i = 2, n_longitudes_tes-1
                    READ(1969,*) idum, idum, jdum, idum, idum, &
                         tes_avg_tau(i,j,k,m,l,1), &
                         tes_avg_tau(i,j,k,m,l,2), &
                         tes_avg_tau_err(i,j,k,m,l,1), &
                         tes_avg_tau_err(i,j,k,m,l,2)
                 END DO
                 tes_avg_tau(1,j,k,m,l,:) = &
                      tes_avg_tau(n_longitudes_tes-1,j,k,m,l,:)
                 tes_avg_tau(n_longitudes_tes,j,k,m,l,:) = &
                      tes_avg_tau(2,j,k,m,l,:)
                 tes_avg_tau_err(1,j,k,m,l,:) = &
                      tes_avg_tau_err(n_longitudes_tes-1,j,k,m,l,:)
                 tes_avg_tau_err(n_longitudes_tes,j,k,m,l,:) = &
                      tes_avg_tau_err(2,j,k,m,l,:)
              END DO
              tes_avg_tau(:,1,k,m,l,:) = tes_avg_tau(:,2,k,m,l,:)
              tes_avg_tau(:,n_latitudes_tes,k,m,l,:) = &
                   tes_avg_tau(:,n_latitudes_tes-1,k,m,l,:)
              tes_avg_tau_err(:,1,k,m,l,:) = tes_avg_tau_err(:,2,k,m,l,:)
              tes_avg_tau_err(:,n_latitudes_tes,k,m,l,:) = &
                   tes_avg_tau_err(:,n_latitudes_tes-1,k,m,l,:)
           END DO
           tes_avg_tau(:,:,1,m,l,:) = tes_avg_tau(:,:,2,m,l,:)
           tes_avg_tau(:,:,n_altitudes_tes  ,m,l,:) = 0.
           tes_avg_tau_err(:,:,1,m,l,:) = tes_avg_tau_err(:,:,2,m,l,:)
           tes_avg_tau_err(:,:,n_altitudes_tes  ,m,l,:) = 0.
        END DO
     END DO
     tes_avg_tau(:,:,:,1,:,:)         = tes_avg_tau(:,:,:,n_ls_tes-1,:,:)
     tes_avg_tau(:,:,:,n_ls_tes,:,:) = tes_avg_tau(:,:,:,2,:,:)
     tes_avg_tau_err(:,:,:,1,:,:)         = &
          tes_avg_tau_err(:,:,:,n_ls_tes-1,:,:)
     tes_avg_tau_err(:,:,:,n_ls_tes,:,:) = &
          tes_avg_tau_err(:,:,:,2,:,:)
     CLOSE(1969)

     RETURN
   END SUBROUTINE dust_tes_limb_init

subroutine get_montabone_ls(curr_my,start_my,ls,mls)

    REAL, INTENT(IN) :: ls
    INTEGER, INTENT(IN) :: curr_my,start_my
    REAL, INTENT(OUT):: mls

! this routine uses the actual ls and a running counter of the mars year to
! create a faux ls that runs from 0->360*n where n is the number of years in
! the chosen selection of the montabone record (i.e. loop_my selected in the 
! namelist). The code will then endlessly loop over this portion of the
! montabone record.

! clock is now "properly" handled in pre_radiation_driver
    mls = real(curr_my-start_my)*360. + ls 

    RETURN
end subroutine get_montabone_ls

subroutine dust_montabone(curr_my, start_my, l_s,           &
                          xlat, xlong, p,                   &
                          optical_depth_multiplier,         &
                          ids,ide, jds,jde, kds,kde,        &
                          ims,ime, jms,jme, kms,kme,        &
                          its,ite, jts,jte, kts,kte,        &
                          band,                             &
                          do_imposed_storm,                 &
                          optdpth_array,                    &
                          tes_cdod_sfc                      )

    IMPLICIT NONE
!---------------------------------------------------------------------
    INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                           ims,ime, jms,jme, kms,kme, &
                                           its,ite, jts,jte, kts,kte

! incoming, in order in which they come in:
    REAL,                                         INTENT(IN   ) :: l_s
    REAL, DIMENSION( ims:ime, jms:jme ),          INTENT(IN   ) :: xlat, xlong
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: p
    REAL,                                         INTENT(IN   ) :: optical_depth_multiplier

    INTEGER, INTENT(IN) :: curr_my, start_my

    CHARACTER(len=3), INTENT(IN) :: band

    LOGICAL :: do_imposed_storm

! outgoing optionals - need to use keyword to trigger
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT), OPTIONAL ::    &
         optdpth_array
    REAL, DIMENSION( ims:ime, jms:jme ),          INTENT(  OUT), OPTIONAL ::    &
         tes_cdod_sfc

    REAL :: ml_s

 !--------------------------------------------------------------------
    ! Local variables
    INTEGER:: i,j,k
    REAL :: pi, gls, zls
    REAL :: taueq, tauS, tauN
    REAL :: topdust, tauref, zp
    REAL :: ratio_fac

    ! source data is IR absorption optical depth
    !  -- assume that the conversion from absorption to extinction is 1.3
    !  -- assume that the conversion from ir to vis is 2.

    IF(band=="TIR") THEN
            ratio_fac = 1.3   ! if request TIR, dust_montabone returns the thermal IR extinction optical depth
    ELSE IF(band=="VIS") THEN
            ratio_fac = 2.6   ! if request VIS, it returns the VIS extinction optical depth
    ELSE
       WRITE ( wrf_err_message , * ) 'dust_montabone was set invalid band request: ',band
       CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
    ENDIF

    CALL get_montabone_ls(curr_my,start_my,l_s,ml_s)

    pi = ACOS(-1.)
    gls = ml_s*pi/180.
    zls = SIN(gls-2.76)

    DO j=jts,jte
    DO i=its,ite
       if(present(optdpth_array)) optdpth_array(i,kte+1,j) = 0.
       topdust = 60.+18.*zls - (32.+18.*zls)*(SIN(xlat(i,j))**4) &
                 - 8.*zls*(SIN(xlat(i,j))**5)

       if(present(tes_cdod_sfc)) tes_cdod_sfc(i,j) = 0.
       
       ! source data is optical depth at 610Pa, but want it at 700Pa for consistency with other MarsWRF code
       tauref = (700./610.)*LINEAR_INTERPOLATE(xlong(i,j),     &
                           xlat(i,j),       &
                           ml_s,            &
                           montabone_longitude, &
                           montabone_latitude,  &
                           montabone_ls,        &
                           montabone_cdod610)

       if(do_imposed_storm) then
         call impose_storm(xlong(i,j),xlat(i,j),ml_s,tauref)
       endif

       if(present(optdpth_array)) then
       DO k=kts,kte
          zp=(700./p(i,k,j))**(70./topdust)
          optdpth_array(i,k,j)= ratio_fac * optical_depth_multiplier * &
                                (tauref/700.) * p(i,k,j) * &
                                MAX( EXP(.007*(1.-MAX(zp,1.))) , 1.e-3 )
       END DO
       endif
       
       if(present(tes_cdod_sfc)) tes_cdod_sfc(i,j) = optical_depth_multiplier * &    ! for particle lifting, want optical
                                (tauref/700.) * p(i,kts,j)                           ! depth at actual surface

    END DO
    END DO

end subroutine dust_montabone

subroutine impose_storm(xlong,xlat,ml_s,tauref)

  implicit none

  real, intent(in   ) :: xlong, xlat, ml_s
  real, intent(inout) :: tauref

  real :: storm_opac

  if(ml_s > minval(storm_ls) .and. ml_s < maxval(storm_ls)) then

   storm_opac = LINEAR_INTERPOLATE(xlong,           &
                                   xlat,            &
                                   ml_s,            &
                                   storm_longitude, &
                                   storm_latitude,  &
                                   storm_ls,        &
                                   storm_cdod610    )

   tauref = max( tauref, storm_opac)

  endif

  return

end subroutine impose_storm

subroutine dust_montabone_init(start_my,loop_my)
  implicit none

  INTEGER, intent(in) :: start_my, loop_my

  INTEGER, parameter :: nlon    = 120, &
                        nlat    = 60, &
                        ntime   = 669, &
                        npadding=2
! temp buffer arrays to read in individual montabone files
  REAL, DIMENSION(nlon,nlat,ntime) :: tmp_cdod610
  REAL, DIMENSION(nlon)            :: tmp_longitude
  REAL, DIMENSION(nlat)            :: tmp_latitude
  REAL, DIMENSION(ntime)           :: tmp_ls

  INTEGER :: nlon_in, nlat_in, ntime_in
  real :: dlong, dlat
  REAL, ALLOCATABLE, DIMENSION(:,:) :: cdod_reverse

  integer, parameter :: n_years_in_files = 12
  character(len=52), dimension(n_years_in_files) :: file_name
  INTEGER, dimension(n_years_in_files) :: mars_years, start_of_rec, end_of_rec

  INTEGER :: climatology_length, first_year, year_count, first_year_m1
  INTEGER :: i, j, k, l, m, iyear, ipad, prior_days
  INTEGER :: idum, jdum
  INTEGER :: ncid, varid
  INTEGER :: istart, iend

  LOGICAL :: read_in_files = .true.

  ! Some years have repeats from prior years. Here we start from the start_of_rec
  ! record in each file and hardwire how many unique records exist in each
  ! file that are not repeated in the next file (since Montabone use 669 
  ! days per year, this number is always going to be either 669 or 668). As 
  ! of the time of writing, MY33 was the last year, so this has 669 by 
  ! default - If you're adding a new year, you will need to check this.
  ! Years 25 and 30 have initial Ls values of near Ls=359
  ! To add extra years, just ammend the mars_year, start_of_rec, end_of_rec
  ! and file_name lists - shouldn't need to futz with anything else in the
  ! rest of this routine.

  ! To obtain the input files: 
  ! http://www-mars.lmd.jussieu.fr/mars/dust_climatology/dataset_v2/krigedcdod_v2.tar.gz
  ! if that's dead, try http://www-mars.lmd.jussieu.fr/mars/dust_climatology

  ! Added preliminary MY34 file, provided by Luca Montabone early via Dropbox

  data mars_years    / 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35/
  data start_of_rec  /  1,  2,  1,  1,  1,  1,  2,  1,  1,  1, 1, 1/
  data end_of_rec    /669,669,669,668,669,669,669,669,668,669,669,669/

  file_name(1)  = "./Data/radiation/montabone/dustscenario_MY24_v2-0.nc"
  file_name(2)  = "./Data/radiation/montabone/dustscenario_MY25_v2-0.nc"
  file_name(3)  = "./Data/radiation/montabone/dustscenario_MY26_v2-0.nc"
  file_name(4)  = "./Data/radiation/montabone/dustscenario_MY27_v2-0.nc"
  file_name(5)  = "./Data/radiation/montabone/dustscenario_MY28_v2-0.nc"
  file_name(6)  = "./Data/radiation/montabone/dustscenario_MY29_v2-0.nc"
  file_name(7)  = "./Data/radiation/montabone/dustscenario_MY30_v2-0.nc"
  file_name(8)  = "./Data/radiation/montabone/dustscenario_MY31_v2-0.nc"
  file_name(9)  = "./Data/radiation/montabone/dustscenario_MY32_v2-0.nc"
  file_name(10) = "./Data/radiation/montabone/dustscenario_MY33_v2-1.nc"
  file_name(11) = "./Data/radiation/montabone/dustscenario_MY34_v2-5.nc"
  file_name(12) = "dustscenario_MY35_v2-5_partial_Ls300.nc"

  ! montabone_rec_len is stored away for use in the montabone reading subroutines
  montabone_rec_len = loop_my

  ! if start_my is a number prior to MY24, it might be a blended year
  if((start_my .ge. 1) .and. (start_my .le. 8)) then
  ! in these cases, we're making a fake average year and so we're setting this up to use the
  ! normal machinery (below) for reading and storing the montabone climatology

          if(loop_my /= 1) then ! don't want to mess with designated "in" variables and this also
                  ! provides a useful check to see if the user really knew what they wanted
               write(wrf_err_message,*) "ERROR: loop_my needs to =1 for average montabone ", &
                  "option using start_my=",start_my," but was loop_my=",loop_my
               CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
          endif
 
          ! pre-blend the composite montabone options elsewhere:
          call get_average_montabone(file_name,n_years_in_files,start_my,  &
               nlon,nlat,ntime,tmp_longitude,tmp_latitude,tmp_ls,tmp_cdod610)

          read_in_files   = .false.
          mars_years(1)   = start_my
          mars_years(2)   = start_my+1
          start_of_rec(1) = 1
          end_of_rec(1)   = 669
  endif

  first_year_m1 = mars_years(1)-1
  climatology_length=0
  do iyear=(start_my-first_year_m1),((start_my-first_year_m1)+(loop_my-1))
    climatology_length = climatology_length+(1+end_of_rec(iyear)-start_of_rec(iyear))
  enddo

  if(start_my .lt. mars_years(1)) then
          write(wrf_err_message,*) "ERROR: trying to start montabone before MY=",mars_years(1), &
                  "but there's nothing in the record prior to this."
          CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
  endif

  if( ((start_my-first_year_m1)+(loop_my-1)) .gt. n_years_in_files) then
          write(wrf_err_message,*) "ERROR: would run off the end of montabone, decrease loop_my ", &
                  "trying to access to MY ",start_my+loop_my-1," but only have to ", &
                  mars_years(n_years_in_files)
          CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
  endif

  !allocate memory for data
  allocate(montabone_cdod610(nlon+npadding*2,nlat+npadding*2,climatology_length+npadding*2))
  allocate(montabone_longitude(nlon+npadding*2))
  allocate(montabone_latitude(nlat+npadding*2))
  allocate(montabone_ls(climatology_length+npadding*2))
  allocate(cdod_reverse(nlon+npadding*2,nlat+npadding*2))
 
  montabone_cdod610=0.

  prior_days=0
  do iyear=(start_my-first_year_m1),((start_my-first_year_m1)+(loop_my-1))
    year_count = iyear - (start_my-mars_years(1))  ! just counting number of yrs used

    if(read_in_files) then

    !open netcdf file
    write(wrf_err_message,*) "MONTABONE CLIM: Loading MY: ",mars_years(iyear),file_name(iyear)
    CALL wrf_message(TRIM(wrf_err_message))
    call check_ncdf(nf90_open(file_name(iyear), nf90_NoWrite, ncid))
    !read the size of the long, lat, time dimensions
    call dimension_size(ncid, "longitude", nlon_in)
    call dimension_size(ncid, "latitude", nlat_in)
    call dimension_size(ncid, "Time", ntime_in)
    if((nlon_in /= nlon) .or. (nlat_in /= nlat) .or. (ntime_in /= ntime)) then
       WRITE ( wrf_err_message , * ) 'montabone_init: input file size unexpected for ',&
               file_name(iyear),' we got ',nlon_in,nlat_in,ntime_in,' but expected ',nlon,nlat,ntime
       CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
    endif
 
    !read in data
    !longitude
    call check_ncdf(nf90_inq_varid(ncid,"longitude",varid),"longitude id")
    call check_ncdf(nf90_get_var(ncid,varid,tmp_longitude),"longitude data")
    !latitude
    call check_ncdf(nf90_inq_varid(ncid,"latitude",varid),"latitude id")
    call check_ncdf(nf90_get_var(ncid,varid,tmp_latitude),"latitude data")
    !ls
    call check_ncdf(nf90_inq_varid(ncid,"Ls",varid),"ls id")
    call check_ncdf(nf90_get_var(ncid,varid,tmp_ls),"ls data")
    !cdod610
    call check_ncdf(nf90_inq_varid(ncid,"cdod610",varid),"cdod610 id")
    call check_ncdf(nf90_get_var(ncid,varid,tmp_cdod610),"cdod610 data")

    endif ! end of the check to see if we're reading in the files rather than using the avg
          ! from get_average_montabone

    !inject the tmp data into appropriate area of the full arrays
    ! the purpose of the pad is to generate:
    !   - lat binning that runs npadding boxes beyond map both north and south (hence 2*npadding)
    !   - lon binning that runs npadding boxes beyond map both east and west
    !   - ls binning that runs npadding boxes before Ls=0 in first year and after 360 in last
    ! Montabone info comes with lat, lon, ls info defining the box centres not the edges

    ! start and stop indices in the time direction in the final arrays
    istart = 1 + npadding + prior_days
    iend   = 1 + npadding + prior_days + (end_of_rec(iyear)-start_of_rec(iyear))
    ! the Ls is going to increase from 0->n*360 where n = loop_my
    montabone_ls(istart:iend) = tmp_ls(start_of_rec(iyear):end_of_rec(iyear)) + (real(year_count)-1.)*360.
    montabone_cdod610((1+npadding):(nlon+npadding),(1+npadding):(nlat+npadding),istart:iend) = &
            tmp_cdod610(:,:,start_of_rec(iyear):end_of_rec(iyear))

    !padding east-west cdod - periodic continuation
    montabone_cdod610(1:npadding,:,istart:iend) = &
            montabone_cdod610(nlon+1:nlon+npadding,:,istart:iend)
    montabone_cdod610(nlon+npadding+1:nlon+2*npadding,:,istart:iend) = &
            montabone_cdod610(npadding+1:2*npadding,:,istart:iend)

    if(iyear == (start_my-first_year_m1)) then  ! only need to fill the lat and lon arrays once
      
      !box steps

      montabone_longitude((1+npadding):(nlon+npadding)) = tmp_longitude(:)
      montabone_latitude((1+npadding):(nlat+npadding))  = tmp_latitude(:)

      dlat  = tmp_latitude(2)  - tmp_latitude(1)

      !padding east-west cdod - periodic continuation
      montabone_longitude(1:npadding) = montabone_longitude(nlon+1:nlon+npadding) - 360.
      montabone_longitude(nlon+npadding+1:nlon+2*npadding) = &
            360. + montabone_longitude(npadding+1:2*npadding)

      !padding north-south
      do i=1,npadding
        montabone_latitude(i) = montabone_latitude(1+npadding) - dlat*real(npadding+1-i)
        montabone_latitude(nlat+npadding+i) = montabone_latitude(nlat+npadding)+dlat*real(i)
      enddo

    endif

    !south pole cdod
    !north pole cdod
    do j=1,npadding
    do i=1, nlon+2*npadding
      montabone_cdod610(i,              j,istart:iend) = &
           montabone_cdod610(nlon+2*npadding-i+1,nlat+npadding,istart:iend)
      montabone_cdod610(i,nlat+npadding+j,istart:iend) = &
          montabone_cdod610(nlon+2*npadding-i+1,   1+npadding,istart:iend)
    end do
    end do

    do i=istart, iend
    !reverse cdod
      cdod_reverse(:,:) = montabone_cdod610(:,npadding*2+nlat:1:-1,i)
      montabone_cdod610(:,:,i) = cdod_reverse(:,:)
    end do

    if(read_in_files) call check_ncdf(nf90_close(ncid))
    prior_days = prior_days + (1 + end_of_rec(iyear) - start_of_rec(iyear))

  enddo  ! end of loop of the years

  !first shall be last; and the last shall be first
  do ipad=1,npadding
    montabone_ls(ipad) = montabone_ls(climatology_length+1-ipad) - real(loop_my)*360.
    montabone_ls(climatology_length+npadding+ipad) = montabone_ls(npadding+ipad) + real(loop_my)*360.
    montabone_cdod610(:,:,ipad) = montabone_cdod610(:,:,climatology_length+1-ipad)
    montabone_cdod610(:,:,climatology_length+npadding+ipad) = montabone_cdod610(:,:,npadding+ipad)
  enddo
  if(montabone_ls(1) .ge. montabone_ls(2)) montabone_ls(1)=montabone_ls(2)-1.
  if(montabone_ls(climatology_length+2*npadding) .le. montabone_ls(climatology_length-1+2*npadding)) &
         montabone_ls(climatology_length+2*npadding) = montabone_ls(climatology_length-1+2*npadding) + 1.

  !reverse lat
  cdod_reverse(1,:)=montabone_latitude(npadding*2+nlat:1:-1)
  montabone_latitude(:) = cdod_reverse(1,:)
  if(allocated(cdod_reverse)) deallocate(cdod_reverse)
 
!convert to radians
  montabone_latitude = montabone_latitude * pi2/360.
  montabone_longitude = montabone_longitude * pi2/360.

end subroutine dust_montabone_init

subroutine prescribed_dust_storm_init(imposed_storm_file)
  implicit none

  INTEGER :: nlon_in, nlat_in, ntime_in, ncat_in
  real :: dlong, dlat
  REAL, ALLOCATABLE, DIMENSION(:) :: lon_in, category_optd, time_in
  INTEGER, ALLOCATABLE, DIMENSION(:) :: category_index
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: mask

  character(len=*), intent(in) :: imposed_storm_file

  INTEGER :: i, j, k, l
  INTEGER :: ncid, varid

  LOGICAL :: index_found

  !open netcdf file
  write(wrf_err_message,*) "READING IMPOSED DUST STORM FILE: ",imposed_storm_file
  CALL wrf_message(TRIM(wrf_err_message))
  call check_ncdf(nf90_open(imposed_storm_file, nf90_NoWrite, ncid))
  !read the size of the long, lat, time dimensions
  call dimension_size(ncid, "x", nlon_in)
  call dimension_size(ncid, "y", nlat_in)
  call dimension_size(ncid, "t", ntime_in)
  call dimension_size(ncid, "category", ncat_in)
  write(wrf_err_message,*) "...file size lat,lon,ls=",nlat_in,nlon_in,ntime_in
  CALL wrf_message(TRIM(wrf_err_message))
 
  allocate(storm_cdod610(nlon_in+2,nlat_in,ntime_in+2))
  allocate(storm_longitude(nlon_in+2))
  allocate(storm_latitude(nlat_in))
  allocate(storm_ls(ntime_in+2))

  allocate(lon_in(nlon_in))
  allocate(time_in(ntime_in))
  allocate(mask(nlon_in,nlat_in,ntime_in))
  allocate(category_optd(ncat_in))
  allocate(category_index(ncat_in))

  !read in data
  !longitude
  call check_ncdf(nf90_inq_varid(ncid,"longitude",varid),"longitude id")
  call check_ncdf(nf90_get_var(ncid,varid,lon_in),"longitude data")
  !latitude
  call check_ncdf(nf90_inq_varid(ncid,"latitude",varid),"latitude id")
  call check_ncdf(nf90_get_var(ncid,varid,storm_latitude),"latitude data")
  !ls
  call check_ncdf(nf90_inq_varid(ncid,"l_s",varid),"ls id")
  call check_ncdf(nf90_get_var(ncid,varid,time_in),"ls data")
  !mask
  call check_ncdf(nf90_inq_varid(ncid,"mask",varid),"mask id")
  call check_ncdf(nf90_get_var(ncid,varid,mask),"mask data")
  !categories of opacty
  call check_ncdf(nf90_inq_varid(ncid,"category_optd",varid),"category_optd id")
  call check_ncdf(nf90_get_var(ncid,varid,category_optd),"category_optd data")
  call check_ncdf(nf90_inq_varid(ncid,"category_index",varid),"category_index id")
  call check_ncdf(nf90_get_var(ncid,varid,category_index),"category_index data")

  storm_ls(2:ntime_in+1) = time_in(:)
  storm_ls(1) = storm_ls(2) - ( storm_ls(3) - storm_ls(2))  ! make a fake time delta_t before first time
  storm_ls(ntime_in+2) = storm_ls(ntime_in+1) + ( storm_ls(ntime_in+1) - storm_ls(ntime_in)) ! and after last

  do i=1,nlon_in
    storm_longitude(i+1) = lon_in(i)
    do j=1,nlat_in
      do k=2,ntime_in+1
        index_found = .false.
        do l=1,ncat_in

          if(mask(i,j,k) == category_index(l) ) then
            storm_cdod610(i+1,j,k) = category_optd(l)
            index_found = .true.
          endif

          if(.not.index_found) then
            write(wrf_err_message,*) "READING IMPOSED DUST STORM FILE: unrecognized index"
            CALL wrf_error_fatal(TRIM(wrf_err_message))
          endif

        enddo
      enddo
    enddo
  enddo

  do k=2,ntime_in+2
    if(storm_ls(k) .lt. storm_ls(k-1)) then
      write(wrf_err_message,*) "...READING IMPOSED DUST STORM FILE: does not currently support time series wrap around ls=360"
      CALL wrf_error_fatal(TRIM(wrf_err_message))
    endif
  enddo

  storm_cdod610(1,:,:) = storm_cdod610(nlon_in+1,:,:)    ! allow the index to wrap around
  storm_cdod610(nlon_in+2,:,:) = storm_cdod610(2,:,:)
  storm_longitude(1) = lon_in(nlon_in+1) - 360.
  storm_longitude(nlon_in+2) = lon_in(2) + 360.

  storm_cdod610(:,:,1) = 0.               ! opacity before start and after end of storm is zero
  storm_cdod610(:,:,ntime_in+2) = 0.

  deallocate(lon_in)
  deallocate(time_in)
  deallocate(mask)
  deallocate(category_optd)
  deallocate(category_index)


end subroutine prescribed_dust_storm_init

subroutine get_average_montabone(file_name,n_years_in_files,start_my,  &
               nlon_in,nlat_in,ntime_in,longitude,latitude,ls, cdod610)

  IMPLICIT NONE

  integer, intent(in) :: n_years_in_files, start_my, nlon_in, nlat_in, ntime_in
  character(len=*), intent(in), dimension(n_years_in_files) :: file_name

  integer, dimension(n_years_in_files) :: use_this_year
  integer :: i, nlon,nlat,ntime, number_of_years, j, k, l
  INTEGER :: ncid, varid

  REAL, DIMENSION(nlon_in,nlat_in,ntime_in) :: tmp_cdod610
  REAL, DIMENSION(nlon_in)                  :: tmp_longitude
  REAL, DIMENSION(nlat_in)                  :: tmp_latitude
  REAL, DIMENSION(ntime_in)                 :: tmp_ls

  REAL, INTENT(OUT), DIMENSION(nlon_in,nlat_in,ntime_in) :: cdod610
  REAL, INTENT(OUT), DIMENSION(nlon_in)                  :: longitude
  REAL, INTENT(OUT), DIMENSION(nlat_in)                  :: latitude
  REAL, INTENT(OUT), DIMENSION(ntime_in)                 :: ls

  cdod610   = 0.
  longitude = 0.
  latitude  = 0.
  ls        = 0.

  SELECT CASE(start_my)

    CASE(1) ! this is the whole 10 years
            use_this_year = (/ 1,1,1,1,1,1,1,1,1,1,1/) 

    CASE(2) ! this is the clear MCS era (MY29-33)
            use_this_year = (/ 0,0,0,0,0,1,1,1,1,1,0/) 

    CASE(3) ! this is all the clear years (MY24, 26-27, 29-33)
            use_this_year = (/ 1,0,1,1,0,1,1,1,1,1,0/) 

    CASE(4) ! this is the clear TES years (MY24, 26)
            use_this_year = (/ 1,0,1,0,0,0,0,0,0,0,0/) 

    CASE(5) ! this is clear TES plus MO (MY24, 26, 27)
            use_this_year = (/ 1,0,1,1,0,0,0,0,0,0,0/) 

    CASE(6) ! this is the YEAR-FROM-HELL (MY25+28) - twice the dust storms, twice the fun!
            ! simple averge
            use_this_year = (/ 0,1,0,0,1,0,0,0,0,0,0/) 

    CASE(7) ! this is HELL x2 (MY25+28) weighted to give full storms in each period
            ! you pay for the whole seat, but you're only... going... to need... THE EDGE!
            ! the 2001 (MY25) storm started near Ls=180 and the 2007 (MY28) nearer Ls=270
            ! this bad boy weights the climatology so neither storm is diluted
            use_this_year = (/ 0,1,0,0,1,0,0,0,0,0,0/)

    CASE(8) ! this is all the clear years (MY24, 26-27, 29-33) but take MIN value
            use_this_year = (/ 1,0,1,1,0,1,1,1,1,1,0/) 


    CASE DEFAULT

            write(wrf_err_message,*) "ERROR: unknown start_my option ",start_my
            CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
   
  END SELECT

  number_of_years = sum(use_this_year)
  write(wrf_err_message,*) "MONTABONE CLIM: Average Case ",start_my," Uses ",number_of_years, &
                           " years in pattern (MY24-34) ",use_this_year
  CALL wrf_message(TRIM(wrf_err_message))

  DO i = 1,n_years_in_files
    if(use_this_year(i) .eq. 1) then

      !open netcdf file
      write(wrf_err_message,*) "MONTABONE CLIM: Loding MY: ",file_name(i)
      CALL wrf_message(TRIM(wrf_err_message))
      call check_ncdf(nf90_open(file_name(i), nf90_NoWrite, ncid))
      !read the size of the long, lat, time dimensions
      call dimension_size(ncid, "longitude", nlon)
      call dimension_size(ncid, "latitude", nlat)
      call dimension_size(ncid, "Time", ntime)
      if((nlon_in /= nlon) .or. (nlat_in /= nlat) .or. (ntime_in /= ntime)) then
         WRITE ( wrf_err_message , * ) 'montabone_init: input file size unexpected for ',&
                 file_name(i),' we got ',nlon,nlat,ntime,' but expected ',nlon_in,nlat_in,ntime_in
         CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
      endif
   
      !read in data
      !longitude
      call check_ncdf(nf90_inq_varid(ncid,"longitude",varid),"longitude id")
      call check_ncdf(nf90_get_var(ncid,varid,tmp_longitude),"longitude data")
      !latitude
      call check_ncdf(nf90_inq_varid(ncid,"latitude",varid),"latitude id")
      call check_ncdf(nf90_get_var(ncid,varid,tmp_latitude),"latitude data")
      !ls
      call check_ncdf(nf90_inq_varid(ncid,"Ls",varid),"ls id")
      call check_ncdf(nf90_get_var(ncid,varid,tmp_ls),"ls data")
      !cdod610
      call check_ncdf(nf90_inq_varid(ncid,"cdod610",varid),"cdod610 id")
      call check_ncdf(nf90_get_var(ncid,varid,tmp_cdod610),"cdod610 data")

      ! deal with ls wrap around:
      if(tmp_ls(1)    .gt. 350.) tmp_ls(1)     = tmp_ls(1)     - 360.
      if(tmp_ls(ntime).lt. 10.)  tmp_ls(ntime) = tmp_ls(ntime) + 360.

      if(start_my .lt. 7) then
          longitude = longitude + tmp_longitude/real(number_of_years)
          latitude  = latitude  + tmp_latitude/real(number_of_years)
          ls        = ls        + tmp_ls/real(number_of_years)
          cdod610   = cdod610   + tmp_cdod610/real(number_of_years)
      else if(start_my .eq. 7) then
          longitude = longitude + tmp_longitude/real(number_of_years)
          latitude  = latitude  + tmp_latitude/real(number_of_years)
          ls        = ls        + tmp_ls/real(number_of_years)
          do j=1,ntime  ! this should be a 1-669 loop
             if((tmp_ls(j) .gt. 175) .and. (tmp_ls(j) .lt. 250)) then
               if(i == 2) then ! this is 2001 storm
                cdod610(:,:,j)   = tmp_cdod610(:,:,j)
               endif
             else if((tmp_ls(j) .gt. 265) .and. (tmp_ls(j) .lt. 325)) then
               if(i == 5) then ! this is 2007 storm
                cdod610(:,:,j)   = tmp_cdod610(:,:,j)
               endif
             else
                ! put the average in all the other times
                cdod610(:,:,j)   = cdod610(:,:,j) + tmp_cdod610(:,:,j)/2.
             endif
          enddo
      else if(start_my .eq. 8) then  ! take minimum value over all years
          longitude = longitude + tmp_longitude/real(number_of_years)
          latitude  = latitude  + tmp_latitude/real(number_of_years)
          ls        = ls        + tmp_ls/real(number_of_years)
          do j=1,ntime
            do k=1,nlon_in
              do l=1,nlat_in
                if(cdod610(k,l,j)==0.) then
                  cdod610(k,l,j)=tmp_cdod610(k,l,j)
                else
                  cdod610(k,l,j)=min(cdod610(k,l,j),tmp_cdod610(k,l,j))
                endif
              enddo
            enddo
          enddo
      endif

      call check_ncdf(nf90_close(ncid))

    endif ! close off the check to see if this is a year we're interested in
  ENDDO

end subroutine get_average_montabone

SUBROUTINE mars_msr_fit(l_s,p8w, p, optdpth_array,         &
                                ids,ide, jds,jde, kds,kde, &
                                ims,ime, jms,jme, kms,kme, &
                                its,ite, jts,jte, kts,kte )
   implicit none

   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde, &
                          ims,ime, jms,jme, kms,kme, &
                          its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   )            ::  l_s
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN ) ::                                        p, &
                                                             p8w

   REAL, INTENT(OUT), DIMENSION( ims:ime, kms:kme, jms:jme ) :: optdpth_array        

! locals
   REAL    :: vscale_fac
   REAL    :: a,b,c,od,ls_shift             ! Used only in the MSR opacity fit option

        !c-mm The MSR project is doing a basic fit to opacity observations from
        !c-mm  Ls=330-193 that consists of three segments (330-30, 30-120, 120-193)
        !c-mm  The purpose of this block of code is to essentially calculate tau
        !c-mm  as a function of L_s using quadratic coefficients that were provided
        !c-mm  to MM on 1-22-20 by Mike Lisano.  This block, then, ignores what
        !c-mm  is in the namelist for optical_depth, and generates a local variable
        !c-mm  'od' which is fed, then, into the dust_distrib_fixed routine.
        !c-mm  We don't care about vertical distribution, only integrated tau.

        !c-mm The Lisano spreadsheet (1-22-20) fits a quadratic to opacity data from Ls=330-193
        !c-mm  but they shift the label by 30 degrees, so it's -30 to 193.  I need to account for
        !c-mm  this shift by converting all Ls values from 330-360 into -30 to 0.

        IF (l_s.GT.325) THEN
           ls_shift=l_s-360.
        ELSE
           ls_shift=l_s
        ENDIF
        SELECT CASE (INT(ls_shift))
            CASE (-35:29)
              !c-mm  Max tau values
!c-mm              a=2.9167e-4
!c-mm              b=-0.0175
!c-mm              c=0.9125
              !c-mm  Min tau values
              a=1.0e-5
              b=-0.0018
              c=0.1779
           CASE (30:44)
              !c-mm  Max tau values
!c-mm              a=0.
!c-mm              b=0.
!c-mm              c=0.65
              !c-mm  Min tau values
              a=1.0e-5
              b=-0.0018
              c=0.1779
           CASE (45:119)
              !c-mm  Max tau values
!c-mm              a=0.
!c-mm              b=0.
!c-mm              c=0.65
              !c-mm  Min tau values
              a=1.0e-5
              b=-0.0018
              c=0.1779
           CASE (120:199)
              !c-mm  Max tau values
!c-mm              a=2.0833e-4
!c-mm              b=-0.05
!c-mm              c=3.65
              !c-mm  Min tau values
              a=1.0e-5
              b=-0.0018
              c=0.1779
           CASE (200:325) !c-mm These values are approximate, but for a season we don't care about
              !c-mm  Max tau values
!c-mm              a=0.
!c-mm              b=0.
!c-mm              c=1.7
              !c-mm  Min tau values
              a=0.
              b=0.
              c=0.6
           CASE DEFAULT
              a=0.
              b=0.
              c=0.
        END SELECT
        od=a*ls_shift*ls_shift+b*ls_shift+c ! Quadratic formula from Lisano email 1-22-20
        od=3.0    !c-mm  For a test case where tau=3
        write(40,*) ls_shift,od
        ! Constant parameter optical depth is now in the namelist
        ! and will be written to the header of the WRFOUT file
        vscale_fac= 0.007
!c-mm        vscale_fac= 0.03
        CALL dust_distrib_fixed(p8w, p, od,                &
                                vscale_fac, optdpth_array, &
                                ids,ide, jds,jde, kds,kde, &
                                ims,ime, jms,jme, kms,kme, &
                                its,ite, jts,jte, kts,kte )

END SUBROUTINE mars_msr_fit

SUBROUTINE get_lower_pbl_ir_heat(lower_pbl_ir_heat,p8w,t,hr_ir, &
                                 ids,ide, jds,jde, kds,kde,     &
                                 ims,ime, jms,jme, kms,kme,     &
                                 its,ite, jts,jte, kts,kte      )

   USE module_model_constants, only: g, r_d, cp

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde, &
                          ims,ime, jms,jme, kms,kme, &
                          its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN ) ::                                        t, &
                                                             p8w, &
                                                           hr_ir

   REAL, DIMENSION(ims:ime, jms:jme ), INTENT(OUT) :: lower_pbl_ir_heat

   REAL :: height, mass

   INTEGER :: i,j,k

   do j=jts,min(jte,jde-1)
   do i=its,min(ite,ide-1)
   height = 0.
   lower_pbl_ir_heat(i,j) = 0.
   do k=kts,kte-1

      height = height + (r_d*t(i,k,j)/g)*log(p8w(i,k,j)/p8w(i,k+1,j))
      if (height < 1500.) then  ! LES shows heating always falls to negligible by 1.5km

         mass = (p8w(i,k,j)-p8w(i,k+1,j))/g  ! kg/m2 in layer
                                             !    W/m2  =     K/s       j/kgK   kg/m2
         lower_pbl_ir_heat(i,j) = lower_pbl_ir_heat(i,j) + (hr_ir(i,k,j)* cp*    mass)

      else
         exit
      endif

   enddo
   enddo
   enddo   

END SUBROUTINE get_lower_pbl_ir_heat

END MODULE module_ra_mars_common
