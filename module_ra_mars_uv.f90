!WRF:MEDIATION_LAYER:PHYSICS

!-----------------------------------------------------------------------
! This module contains parameters needed to deal with the UV section   !
! of the Hadley radiation scheme.  The method is derived from a paper  !
! by Gonzalez-Galindo et al. (2005)                                    !
! Species number for purposes of this code:                            !
!            1:   CO2                                                  !
!            2:   O2                                                   !
!            3:   O(3P)                                                !
!            4:   H2O                                                  !
!            5:   H2                                                   !
!            6:   H2O2                                                 !
!            7:   O3                                                   !
!-----------------------------------------------------------------------
   MODULE uv_parameters

      USE module_model_constants
      USE module_wrf_error
      USE module_ra_mars_common, only: cp_mars
      USE module_state_description, ONLY: UV_GONGAL

      IMPLICIT NONE

      REAL(KIND(0.d0)), PARAMETER :: hc           = dble(h_planck)*dble(light_speed) ! 1.986D-25 ! Planck constant * speed of light
      REAL(KIND(0.d0)), PARAMETER :: local_mw_co2 = dble(mw_co2)/1000.d0             ! 0.044     ! [kg/mol]
      REAL,             PARAMETER :: solar_date   = 1993.4     ! Year of solar medium conditions.  Solar min=1996.4, Solar max=1990.6
      REAL(KIND(0.d0)), PARAMETER :: pi           = dble(pi2)/2.d0      ! Pi
      REAL(KIND(0.d0)), PARAMETER :: scale_height = dble(r_d)*dble(t0)/dble(g) ! this gives about 15 km, which may be too big
      REAL(KIND(0.d0)), PARAMETER :: radius       = 1.d0/dble(reradius) ! radius

         INTEGER, PARAMETER ::                                         &
                               n_colden_levs =                    253, &
                                n_uv_species =                      7, &
                              n_uv_intervals =                     36

         REAL(KIND(0.d0)), PARAMETER ::                                &
                                 uv_heat_eff =                   0.22, & ! UV+EUV heating efficiency (from Fox, 1988)
                               mol_constants =   dble(N_avogadro)*     &
                                                 co2_mixing_ratio/     & 
                                                 mw_co2/g           

         REAL(KIND(0.d0)), DIMENSION(n_uv_intervals) ::                &
                                                        uv_bandwidths= &
                           (/ 4.9, 5.3, 6.1, 12.8, 2.0, 6.9, 11.8,     & ! Width of each of the 36 intervals [nm]
                              0.4, 1.9, 4.9, 0.9, 4.0, 6.9, 4.9, 4.9,  &
                              10.2, 5.6, 1.9, 0.9, 0.9, 0.9, 0.9, 3.7, &
                              13.1, 3.9, 6.0, 4.8, 23.9, 8.4, 8.5,     &
                              26.8, 7.5, 20.8, 9.0, 97.6, 462.2 /)

         REAL(KIND(0.d0)), DIMENSION(n_uv_intervals+1) ::              &
                                                          lambda_edge= &
                        (/ 0.1, 5.0, 10.4, 16.5, 29.4, 31.5, 38.5,     &
                            50.4, 51.4, 53.4, 58.4, 59.4, 63.5, 70.5,  &
                            75.5, 80.5, 90.8, 96.5, 98.5, 99.5, 100.5, &
                            101.5, 102.5, 106.3, 119.5, 123.5, 129.6,  &
                            134.5, 158.5, 167.0, 175.6, 202.5, 210.0,  &
                            230.9, 240.0, 337.7, 800.0 /)

         REAL(KIND(0.d0)), DIMENSION(n_uv_intervals) ::                &
                                                           lambda_avg= &
                        (/ 2.55, 7.65, 13.45, 22.9, 30.4, 34.95, 44.4, & ! Median wavelength of each of the 36
                           50.6, 52.35, 55.85, 58.85, 61.4, 66.95,     & !    intervals [nm]
                           72.95, 77.95, 85.6, 93.6, 97.45, 98.95,     &
                           99.95, 100.95, 101.95, 104.35, 112.85,      &
                           121.45, 126.5, 132.0, 146.45, 162.8,        &
                           171.25, 189.0, 206.25, 220.4, 235.4, 288.8, &
                           568.8 /)

         REAL(KIND(0.d0)), DIMENSION(n_uv_intervals) ::                &
                                                         uv_flux_band 

         INTEGER, DIMENSION(n_uv_species) ::                           &
                                                            int_start= &
                       (/ 2, 1, 1, 25, 1, 25, 34 /)                      ! This only works because there are no 'gaps' in the
                                                                         !    absorption bands of these 7 gases in the UV
         INTEGER, DIMENSION(n_uv_species) ::                           &
                                                              int_end= &
                       (/ 32, 34, 16, 31, 15, 35, 36 /)

   END MODULE uv_parameters

!-----------------------------------------------------------------------
! This module contains all the UV-specific routines, following the     !
!    work of Gonzalez-Galindo et al. (2005).  One thing to note is     !
!    that the order of the vertical dimension follows WRF (i.e. k=1 is !
!    the surface-most layer, and k=n_layer is at TOA.                  !
!-----------------------------------------------------------------------
   MODULE module_ra_mars_uv

      USE uv_parameters
      USE module_ra_chapman

      IMPLICIT NONE

      REAL(KIND(0.d0)), DIMENSION(n_colden_levs,n_uv_intervals) ::     &
                                                               colden    ! Column density array from Gonzalez-Galindo et al. (2005)

      REAL(KIND(0.d0)), DIMENSION(n_colden_levs,n_uv_intervals,n_uv_species) ::  &
                                                                j_Xn     ! Photoabsorption coefficients for all gases from Gonzalez-Galindo et al. (2005)

      REAL(KIND(0.d0)), DIMENSION(24,4) ::                             &
                                                          flux_coeffs    ! Flux coefficients c1, p1, c2, p2 from Gonzalez-Galindo et al., (2005) Table 1.
   CONTAINS

     SUBROUTINE mars_uv(uv_physics,coszen,p3d,pf3d,   &
                        pi3d, t3d, dz8w,              &
                        hr_g_uv, rthraten,            &
                        ids, ide, jds, jde, kds, kde, & 
                        ims, ime, jms, jme, kms, kme, &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::          ids, ide, jds, jde, kds, kde, &
                                         ims, ime, jms, jme, kms, kme, &
                                         its, ite, jts, jte, kts, kte, &
                                                           uv_physics

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::   &
                                                                  p3d, &
                                                                 pf3d, &
                                                                 pi3d, &
                                                                  t3d, &
                                                                 dz8w


      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) ::    coszen 

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::   &
                                                              hr_g_uv     ! Heating rate due to gas in UV wavelengths (kinetic)
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) ::   &
                                                             rthraten     !  Total potential temperature heating rate

!c Local variables
      INTEGER, PARAMETER ::          number_of_fake_layers_uv    = 1    !  i.e. number of layers above p_top that need rad tran

      INTEGER            ::                             n_layer_wrf  , &
                                                         n_layer_uv    

      REAL, DIMENSION(ite-its+1, kte-kts+1) ::              hr_g_uv2d    ! 2-D UV heating rate returned from uv_heating

      REAL(KIND(0.d0)), DIMENSION(ite-its+1,kte-kts+1+number_of_fake_layers_uv+1) ::   pf2duv    ! Two-dimensional level (edge) pressure for UV

      REAL(KIND(0.d0)), DIMENSION(ite-its+1,kte-kts+1+number_of_fake_layers_uv)   ::    p2duv, & ! Two-dimensional layer (centre) pressure for UV
                                                                t2duv, & ! Two-dimensional layer temperature for UV
                                                                 dz2d    ! Two-dimensional layer thickness
      REAL(kind(0.d0)), PARAMETER ::        pressure_TOA =      1.d-8    ! top of atmosphere pressure (in Pa) (machine small but non-zero) 10^-8 is about 300km on Mars

      LOGICAL, PARAMETER ::  verbose_debug = .false.
      INTEGER :: i,j,ii,k,kk,n_profile,j_end

      n_layer_wrf     = kte-kts+1
      n_layer_uv      = n_layer_wrf + number_of_fake_layers_uv

      hr_g_uv(its:ite,kts:kte,jts:jte)=0.0  ! this puppy is a wrf domain variable - don't write outside of local tile!

#ifdef mpas
      n_profile=ite-its+1
      j_end=jte
#else
      n_profile=MIN(ite,ide-1)-its+1
      j_end=MIN(jte,jde-1)
#endif

      j_loop: DO j=jts,j_end

         DO i=1,n_profile
            ! Convert from lower bound of 1 (local) to lower bound of
            ! its (WRF distributed memory, or serial with its /= 1)
            ii = i+its-1
            ! Variables on the z-stagger
  
! The UV only needs one fake layer to absorb anything between ptop and the Sun
            pf2duv(i,1) = REAL(pf3d(ii,1,j),KIND(0.d0))           ! surface pressure (pressure at bottom edge of lowest level)
            DO k=1,n_layer_wrf
               pf2duv(i,k+1) = REAL(pf3d(ii,k+1,j),KIND(0.d0))    ! pf2duv is pressure not flipped, since it goes into UV code only
               t2duv (i,k)   = REAL(t3d(ii,k,j),KIND(0.d0))       ! Used in the UV code, so not flipped.  Same orientation as WRF (bottom to top)
               p2duv (i,k)   = REAL(p3d(ii,k,j),KIND(0.d0))       ! Used in the UV code, so not flipped.  Same orientation as WRF
               dz2d  (i,k)   = REAL(dz8w(ii,k,j),KIND(0.d0))      ! Used in UV code, so not flipped
            ENDDO
            IF(number_of_fake_layers_uv .eq. 1) THEN
               k=n_layer_uv
               pf2duv(i,k+1) = pressure_TOA
               t2duv (i,k) = REAL(t3d(ii,k-1,j),KIND(0.d0)) 
               p2duv (i,k) = REAL((pf2duv(i,k+1)+pf2duv(i,k))/2.,KIND(0.d0)) ! intent: centre of fake layer at half mass of layer
               ! not going to define dz2d at the top - Leonard's Law: if there's an attempt to use it, you're already someplace you shouldn't be
             ELSE
               WRITE ( wrf_err_message ,*) 'ERROR: KDM has not been setup to use other than one fake uv layer ',& 
               number_of_fake_layers_uv 
               CALL wrf_error_fatal ( wrf_err_message )
            ENDIF

            SELECT CASE (uv_physics)
             CASE (UV_GONGAL)
! Only do UV if the Sun is up.
               IF (coszen(ii,j) > 0.00001) THEN
                  CALL uv_heating(SNGL(t2duv(i,:)),SNGL(p2duv(i,:)),   &
                                  SNGL(pf2duv(i,:)),hr_g_uv2d(i,:),    &
                                  coszen(ii,j),SNGL(dz2d(i,:)),n_layer_uv,   &
                                  n_layer_wrf)
                  ii=i+its-1
! Here we don't have to do any flipping in the vertical, 
!  since UV and WRF have the same orientation.
                  DO kk=kts,kte
                     hr_g_uv(ii,kk,j)=hr_g_uv2d(i,kk) 

                     rthraten(ii,kk,j)=rthraten(ii,kk,j)+(hr_g_uv(ii,kk,j)/pi3d(ii,kk,j))

                  ENDDO
               ENDIF
IF(verbose_debug) THEN
               if((ii==(its+ite)/2) .and. (j==jts)) then
                  DO kk=kts,kte
                        write(0,*) "uv rthraten: ",+hr_g_uv(ii,kk,j)/   &
                                       pi3d(ii,kk,j)
                  ENDDO
               endif
ENDIF
             CASE DEFAULT
!c-mm                CALL wrf_debug( 300 , "in mars_uv, but not doing any uv as uv_physics option not recognized.")
             END SELECT
         ENDDO

      ENDDO j_loop

     END SUBROUTINE mars_uv

!-----------------------------------------------------------------------
!+ Subroutine to determine the UV heating rate as per method described !
!  in Gonzalez-Galindo et al. (2005).                                  !
!                                                                      !
! Method:                                                              !
!       Loop over the individual 'intervals', summing up contributions !
!       from each gas in each interval.  Right now it's only for CO2   !
!-----------------------------------------------------------------------
      SUBROUTINE uv_heating(t,p,pf,hr_g_uv2d,mubar,dz,n_layer_uv,n_layer_wrf)

      IMPLICIT NONE


      REAL, INTENT(IN   ) ::                                           &
                                                                mubar

      
      INTEGER, INTENT(IN) ::                               n_layer_uv, &
                                                          n_layer_wrf

      REAL, DIMENSION(n_layer_uv), INTENT(IN   ) ::                    &
                                                                    t, & ! Temperature of each layer
                                                                    p, & ! Pressure of each layer
                                                                   dz

      REAL, DIMENSION(n_layer_uv+1), INTENT(IN   ) ::                 pf  ! Pressure at each level (layer edges)

      REAL,                                                            &
          DIMENSION(n_layer_wrf), INTENT(  OUT) ::                     &
                                                            hr_g_uv2d

!c-mm Local variables

      INTEGER                                                 i,j,k,l

      REAL(KIND(0.d0)) ::                                     airmass

      REAL(KIND(0.d0)), DIMENSION(n_layer_uv) ::                       &
                                                           int_d_mass, &
                                                             molperm2, &
                                                             test

      REAL(KIND(0.d0)),                                                &
       DIMENSION(n_layer_uv, n_uv_intervals, n_uv_species) ::          &
                                                              j_array

      REAL(KIND(0.d0)), DIMENSION(24:32) ::                            &
                                                            sigma_195= & ! Average cross-section for intervals 24-32 @ 195 K
                    (/ 2.05864D-17, 5.90557D-20, 3.10270D-19,          &
                       6.70653D-19, 4.55132D-19, 8.87122D-20,          &
                       1.32138D-20, 7.22244D-23, 2.88002D-26 /)

      REAL(KIND(0.d0)), DIMENSION(24:32) ::                            &
                                                            sigma_295= & ! Average cross-section for intervals 24-32 @ 295 K
                    (/ 2.05897D-17, 6.71104D-20, 3.45509D-19,          &
                       7.45711D-19, 4.82752D-19, 1.11594D-19,          &
                       1.98308D-20, 1.38530D-22, 2.14140D-25 /)

      REAL(KIND(0.d0)), DIMENSION(24:32) ::                            &
                                                                alpha

      REAL(KIND(0.d0)), PARAMETER ::                                   &
                        t_min =                                  195., &
                        t_max =                                  295.

      j_array=0.d0                                                        ! Initialize j_array to zero everywhere

      alpha=((sigma_295/sigma_195)-1.)/(295.-195.)

      airmass=SQRT((radius/scale_height*mubar)**2+ 2.*radius/          &  ! Modeling the atmosphere as a simple spherical shell
              scale_height+1)-radius/scale_height*mubar                   !    to improve accuracy near horizon (vs. sec approx.)
                                                                          ! From Schoenberg (1929) (via Wikipedia for 'Air mass')
                                                                          ! Yeah, that's right, we went there.  That's a Wikipedia
                                                                          ! reference. Wanna make something of it? huh?

      DO k=1,n_layer_uv  ! This is going bottom (k=1) to top (k=n_layer)
!c-mm--------------------------
!c-mm Next two statements are to be used if we employ a Chapman function
!c-mm    instead of simple SZA.
!         molperm2(k) = (pf(k)-pf(k+1))*mol_constants*                  &
!                       atm_chapman((1./reradius/1000.+SUM(dz(1:k)))/   &
!                       11000.,SNGL(acos(mubar)*180./pi))
!         int_d_mass(k) = pf(k)*mol_constants*                          &
!                         atm_chapman((1./reradius/1000.+SUM(dz(1:k)))/ &
!                         11000.,SNGL(acos(mubar)*180./pi))
!c-mm--------------------------
!c-mm In conversation with Francisco Gonzalez-Galindo, he indicated that we
!c-mm    needed to divide by cos(SZA), or mubar.  I'm finding that we need
!c-mm    to multiply instead, which makes more intuitive sense.  Results look
!c-mm    better, too, especially around the twilight regions.
!c-mm Here's my argument...incoming flux as W/m2 is multiplied by cos(SZA)
!c-mm    to account for angle, thus it =0 at SZA=90.  When calculating the
!c-mm    heating rate at the end, we are calculating W/m3, or, another way,
!c-mm    W/m2/m.  The W/m2 term portion is obtained from hc/lambda*j*molperm2
!c-mm    Notice that we're just incorporating mubar (=cos(SZA)) into molperm2
!c-mm    but that it's equivalent to multiplying the W/m2 portion by cos(SZA)
!c-mm    at the very end.  In a sense, we're just 'hiding' the angular
!c-mm    dependence inside molperm2.  We could just as easily pull it out
!c-mm    and explicitly multiply the heating rate term by cos(SZA).  The way
!c-mm    it's done is, indeed, correct (mm 10-12-10)
         molperm2(k) = (pf(k)-pf(k+1))*mol_constants*mubar                   ! Convert from pressure to number of molecules per m2. Equivalent to n(z)dz at level z
         int_d_mass(k) = pf(k)*mol_constants*mubar                           ! Number of molecules above and including layer k. Both these are units of m-2
      ENDDO                   

      FORALL (k=1:n_layer_uv,                                          & ! For just CO2 in the appropriate
              j=int_start(1):int_end(1))                               & !    spectral intervals for each point in the GCM grid,
              j_array(k,j,1)=get_j(int_d_mass(k)/1.D4,j,1)               !    j_array holds the partial photoabsorption coeff. @ 195 K. 1.D4 goes from m-2 to cm-2

      FORALL (k=1:n_layer_uv, j=1:24)                                  & ! This does the scaling for solar flux.  Currently, solar date is hardwired
              j_array(k,j,1)=j_array(k,j,1)*scale_flux(j,solar_date)              !    to solar medium conditions.

      DO k=1,n_layer_wrf  ! only need to calculate heat rating for the prognostic wrf layers
!c-mm            DO l=1,n_uv_gases
         DO l=1,1                                                     ! Only for CO2 right now
            DO j=24,32                                                ! Only bands 24-32 have CO2 temperature dependence.
               IF (l==1) THEN
                  j_array(k,j,l)=j_array(k,j,l)*                    &
                                   EXP(-sigma_195(j)*alpha(j)*      &
                                   SUM((MAX(MIN(DBLE(t(k:n_layer_uv)), & ! <---This MAX(MIN(...)) block limits the temperature
                                   t_max),t_min)-t_min)*            & !    dependence to only between 195 and 295 K.  Above/
                                   molperm2(k:n_layer_uv)))*        & !    below this range, value is pegged at 295/195 K
                                   (1+alpha(j)*(MAX(MIN(DBLE(t(k)), & !    respectively.
                                   t_max),t_min)-t_min))                   
               ELSE
                  j_array(k,j,l)=j_array(k,j,l)*                    & ! The purpose of these two j_array blocks is to adjust
                                   EXP(-sigma_195(j)*alpha(j)*      & !    the photoabsorption coefficients based on the varying
                                   SUM((MAX(MIN(DBLE(t(k:n_layer_uv)), & !    CO2 cross-section with temperature.  Slightly
                                   t_max),t_min)-t_min)*            & !    different adjustment for CO2 vs. other gases.
                                   molperm2(k:n_layer_uv)))
               ENDIF
            ENDDO
         ENDDO

!c-mm The t(k)/p(k)/4. term comes about because we want to convert the heating rate [W/m3] into [K/s]
!c-mm    this requires dividing the former by density*cp.  Density is p/RT, and R=cp/4. This
!c-mm    cancels out the cp values, leaving us with t/4p
         hr_g_uv2d(k)=1.d9*uv_heat_eff*hc*molperm2(k)/dz(k)*(t(k)/  &
                      p(k)/4.)*sum(j_array(k,1:n_uv_intervals,1)/   &
                      lambda_avg(1:n_uv_intervals))
      ENDDO

      END SUBROUTINE uv_heating

!-----------------------------------------------------------------------
!+ Function to interpolate the photoabsorption coefficient from data   !
!  tables as described in Gonzalez-Galindo et al. (2005).              !
!                                                                      !
! Method:                                                              !
!       Searches down the corresponding 'interval' column in coln.dat  !
!       and interpolates column number density.  Takes the appropriate !
!       interval and interpolates the proper photoabsorption           !
!       coefficient for each gas and each of the 36 spectral intervals !
!       in the UV (only for those in which the gas is active)          !
!-----------------------------------------------------------------------
      PURE FUNCTION get_j(int_d_mass, j, l)                              ! Pure function has no side effects
                                                                         !   (i.e. doesn't open files or change
      IMPLICIT NONE                                                      !   global variable values)

      INTEGER, INTENT(IN   ) ::                                        &
                                                                 j, l

      REAL(KIND(0.d0)), INTENT(IN   ) ::                               &
                                                           int_d_mass

!c-mm Local variables

      INTEGER                                                       i

      REAL(KIND(0.d0))                                                 &
                                                                get_j, &
                                                                 frac, &
                                                       log_int_d_mass


      log_int_d_mass=DLOG10(int_d_mass)

      DO i=1,n_colden_levs,1
         IF (colden(i,j) > int_d_mass) EXIT                              ! Find the index for point greater than int_d_mass
      ENDDO
      ! If int_d_mass > colden(n_colden_levs,j), then the loop
      ! will finish, and i will end up with a value of n_colden_levs+1
      ! and the equations for frac and get_j below will (or should...)
      ! fail with an illegal array address
      IF (i > n_colden_levs) i=n_colden_levs

      frac=(DLOG10(colden(i,j))-log_int_d_mass)/(DLOG10(colden(i,j))-  &
            DLOG10(colden(i-1,j)))                                       ! Fractional difference

      get_j=(1-frac)*j_Xn(i,j,l)+frac*j_Xn(i-1,j,l)                      ! Interpolate within j_Xn for the appropriate

      END FUNCTION get_j

!-----------------------------------------------------------------------
!+ Function to calculate the scale factor for the solar flux depending !
!  on season                                                           !
!                                                                      !
! Method:                                                              !
!       Table in solar_flux.dat contains 4 parameters (taken from      !
!       Gonzalez-Galindo et al., (2005) for each wavelength bin from   !
!       1-24.  A parameterized function is used to calculate a scaling !
!       value (which should be ~between 0.5-1.5) which is multiplied   !
!       to the value of j_array                                        !
!-----------------------------------------------------------------------
      PURE FUNCTION scale_flux(j,solar_date)

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                                    j

      REAL, INTENT(IN   ) ::                                           &
                                                           solar_date

!c-mm Local variables

      REAL(KIND(0.d0))                                                 &
                                                           scale_flux

      scale_flux=(flux_coeffs(j,1)+flux_coeffs(j,2)*solar_date)*       &
                 (SIN(2*pi/11*(solar_date-1985.-pi)))+(flux_coeffs(j,3)+        &
                 flux_coeffs(j,4)*solar_date)

      END FUNCTION scale_flux

!-----------------------------------------------------------------------
!+ Subroutine calculating the fraction of blackbody radiation          !
!  between two wavelengths, assuming this is for solar wavelengths.    ! 
!                                                                      !
! Method:                                                              !
!        Uses series formulae from Houghton 'Physics of Atmospheres'   !
!-----------------------------------------------------------------------
      SUBROUTINE solar_fraction(wl1, wl2, denom, fst4)

      IMPLICIT NONE

      REAL(KIND(0.d0)), PARAMETER ::    			       &
                                          c  =         1.53989733D-01, &
                                         c2  =            1.43883D+04, &
                                       temp  =             5777.0D+00    ! Equivalent blackbody temperature of the Sun

      REAL(KIND(0.d0)), INTENT(IN   ) ::			       &
                                                                  wl1, &
                                                                  wl2, &
                                                                denom    ! Term to calculate *fraction* of total flux rather
                                                                         !    than flux itself.
      REAL(KIND(0.d0)), INTENT(  OUT) ::                         fst4

!c-mm  Local variables

      REAL(KIND(0.d0)), DIMENSION(2) ::   		               &
                                                       	           wl, &
                                                                    f

      REAL(KIND(0.d0)) ::                                              &
                                                                 facc, &
                                                                   ff, &
                                                                    w, &
                                                                    v, &
                                                                   wv

      INTEGER         					               &
                                                                 i, m


      wl(1)=wl1
      wl(2)=wl2
      facc=0.00001
      DO i=1,2
         f(i)=0.
         IF (wl(i) <= 0.0) CYCLE
         v=c2/(wl(i)*temp)
         DO m=1,1000
            w=FLOAT(m)
            wv=w*v
            ff=(w**(-4))*(EXP(-wv))*(((wv+3.)*wv+6.)*wv+6.)
            f(i)=f(i)+ff
            IF (f(i) == 0.0) CYCLE
            IF (ff/f(i) <= facc) EXIT
         ENDDO
      ENDDO
      fst4=(c*(f(2)-f(1))*5.67D-08*temp*temp*temp*temp)/denom
      END SUBROUTINE

!-----------------------------------------------------------------------
!+ Subroutine to read in the column densities and UV photoabsorption   !
!  cross-sections.                                                     !
!                                                                      !
! Method:                                                              !
!	Straightforward.                                               !
!       For now, hardwire these data filenames.                        !
!-----------------------------------------------------------------------
      SUBROUTINE uvinit

      IMPLICIT NONE

!c-mm Local variables

      INTEGER                                                    i, j

      CHARACTER (LEN=80) ::                                      junk

      OPEN(UNIT=2010,FILE='./Data/radiation/coln_full.dat',STATUS='old')
      READ(2010,'(A80)') junk
      READ(2010,'(A80)') junk
      READ(2010,*) ((colden(i,j), j=1,n_uv_intervals), i=1,n_colden_levs)
      CLOSE(2010)

      OPEN(UNIT=2011,FILE='./Data/radiation/j1_Xn.dat',STATUS='old')
      READ(2011,'(A80)') junk
      READ(2011,'(A80)') junk
      READ(2011,*) ((j_Xn(i,j,1), j=1,n_uv_intervals), i=1,n_colden_levs)
      CLOSE(2011)

      OPEN(UNIT=2012,FILE='./Data/radiation/j2_Xn.dat',STATUS='old')
      READ(2012,'(A80)') junk
      READ(2012,'(A80)') junk
      READ(2012,*) ((j_Xn(i,j,2), j=1,n_uv_intervals), i=1,n_colden_levs)
      CLOSE(2012)

      OPEN(UNIT=2013,FILE='./Data/radiation/j3_Xn.dat',STATUS='old')
      READ(2013,'(A80)') junk
      READ(2013,'(A80)') junk
      READ(2013,*) ((j_Xn(i,j,3), j=1,n_uv_intervals), i=1,n_colden_levs)
      CLOSE(2013)

      OPEN(UNIT=2014,FILE='./Data/radiation/j4_Xn.dat',STATUS='old')
      READ(2014,'(A80)') junk
      READ(2014,'(A80)') junk
      READ(2014,*) ((j_Xn(i,j,4), j=1,n_uv_intervals), i=1,n_colden_levs)
      CLOSE(2014)

      OPEN(UNIT=2015,FILE='./Data/radiation/j5_Xn.dat',STATUS='old')
      READ(2015,'(A80)') junk
      READ(2015,'(A80)') junk
      READ(2015,*) ((j_Xn(i,j,5), j=1,n_uv_intervals), i=1,n_colden_levs)
      CLOSE(2015)

      OPEN(UNIT=2016,FILE='./Data/radiation/j6_Xn.dat',STATUS='old')
      READ(2016,'(A80)') junk
      READ(2016,'(A80)') junk
      READ(2016,*) ((j_Xn(i,j,6), j=1,n_uv_intervals), i=1,n_colden_levs)
      CLOSE(2016)

      OPEN(UNIT=2017,FILE='./Data/radiation/j7_Xn.dat',STATUS='old')
      READ(2017,'(A80)') junk
      READ(2017,'(A80)') junk
      READ(2017,*) ((j_Xn(i,j,7), j=1,n_uv_intervals), i=1,n_colden_levs)
      CLOSE(2017)

      OPEN(UNIT=2018, FILE='./Data/radiation/solar_flux.dat',STATUS='old')
      READ(2018,'(A80)') junk
      READ(2018,'(A80)') junk
      READ(2018,*) ((flux_coeffs(i,j), j=1,4), i=1,24)
      CLOSE(2018)

      END SUBROUTINE 

    END MODULE module_ra_mars_uv
