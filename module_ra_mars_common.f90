!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_mars_common

  USE module_wrf_error
  USE module_nrutils, only: linear_interpolate
  use module_mars24
!c-mm  use netcdf
  use module_model_constants, only: pi2
!c-mm        use module_planet_utilities, only: read_wrf_profile_file, &

  IMPLICIT NONE

  PRIVATE

public :: newton20, dust_distrib_fixed, dust_distribution, oxford_dust, &
     mcd_mgs, mcd_mgsx2, mcd_viking, mcs_dust, init_mcd_mgs, &
     sw_aerosol_scatter, lw_aerosol_heat, demiss3, &
         inject_mcd_mgs, delta_mcd_mgs, da_mcd_mgs, diagnose_tau_2d, &
         mars_msr_fit, chicagodust, get_lower_pbl_ir_heat, cp_mars, &
!c-mm         dust_distrib_truly_fixed, forced_opac_profile, fop_init, &
                  dust_distrib_truly_fixed

  PUBLIC :: dust_tes_limb
  PUBLIC :: dust_tes_limb_init

  PUBLIC :: dust_montabone
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

!c-mm!forced opacity profile
!c-mm  INTEGER, SAVE :: fop_n_heights, fop_n_times
!c-mm  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: fop_pres, fop_opac
!c-mm  REAL, SAVE, ALLOCATABLE, DIMENSION(:)   :: fop_time
!c-mm!end forced opacity profile

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

!c-mm   !====================================================================
!c-mm   SUBROUTINE forced_opac_profile(p_in, optdpth_array, looper,     &
!c-mm                                 ids,ide, jds,jde, kds,kde,        &
!c-mm                                 ims,ime, jms,jme, kms,kme,        &
!c-mm                                 its,ite, jts,jte, kts,kte        )
!c-mm
!c-mm! this routine forces the dust opacity to be that read in from a horizontal
!c-mm! averge profile file, but with time and vertical variation
!c-mm!---------------------------------------------------------------------
!c-mm     IMPLICIT NONE
!c-mm!---------------------------------------------------------------------
!c-mm     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
!c-mm                                            ims,ime, jms,jme, kms,kme, &
!c-mm                                            its,ite, jts,jte, kts,kte
!c-mm
!c-mm     INTEGER, INTENT(IN)  :: looper ! 0 = match time exactly to the second
!c-mm                                    ! 1 = day looper, match local time to second
!c-mm                                    ! 2 = year looper, match ls
!c-mm!
!c-mm     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::    &
!c-mm                                                        optdpth_array
!c-mm
!c-mm     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::    &
!c-mm                                                                 p_in
!c-mm 
!c-mm     INTEGER:: i,j,k, it
!c-mm
!c-mm     ! put code in here to find right it, or instead do bilinear interp
!c-mm
!c-mm     DO j=jts,jte
!c-mm     DO i=its,ite
!c-mm        optdpth_array(i,kts:kte,j)=0.
!c-mm        DO k= kts, kte+1
!c-mm           optdpth_array(i,k,j)=linear_interpolate(x0=p_in(i,k,j), &
!c-mm                   x=fop_pres(:,it),y=fop_opac(:,it),              &
!c-mm                   log_x=.true.,                                   &
!c-mm                   out_of_range_use_nearest_edge=.true.)
!c-mm        END DO
!c-mm     END DO
!c-mm     END DO
!c-mm
!c-mm     RETURN
!c-mm   END SUBROUTINE forced_opac_profile
!c-mm
!c-mm   SUBROUTINE fop_init(filename)
!c-mm
!c-mm   implicit none
!c-mm
!c-mm   character(len=*), intent(in) :: filename
!c-mm
!c-mm   CALL read_wrf_profile_file(filename=filename,       &
!c-mm                              memo="forced opacity",   &
!c-mm                              n_times=fop_n_times,     &
!c-mm                              n_heights=fop_n_heights, &
!c-mm                              pres=fop_pres,           &
!c-mm                              opac=fop_opac,           &
!c-mm                              time=fop_time            )
!c-mm
!c-mm   END SUBROUTINE fop_init

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
