!WRF:MODEL_LAYER:PHYSICS
!c-mm!XREF:MEDIATION_LAYER:PHYSICS
!

!-----------------------------------------------------------------------
! This module contains most parameters needed to operate the Hadley    !
! radiation scheme.  For each variable type, there are two blocks of   !
! variables.  The first block contains fixed values that should never  !
! change.  The second block contains values that can be modified       !
! according to need.                                                   !
!-----------------------------------------------------------------------
      MODULE module_kdm_parameters
      USE, INTRINSIC :: IEEE_ARITHMETIC  ! to catch NaN and Inf

#ifndef mpas
      USE module_model_constants
      USE module_diag_common, only: test_for_bad, test_for_bad_detailed, test_for_bad_detailed_2d
#endif

      IMPLICIT NONE
!                               -------- hardwired case-independent constants --- DO NOT CHANGE THESE!!!

         INTEGER, PARAMETER ::                                              & ! THESE VARIBLES SHOULD REMAIN FIXED
                             ip_scf_solar_up =                           1, & ! Index of source coefficient for upward solar beam
                           ip_scf_solar_down =                           2, & ! Index of source coefficient for downward solar beam
                                ip_scf_ir_1d =                           1, & ! Index of source coefficient for 1st difference of Planckian
                                ip_scf_ir_2d =                           2, & ! Index of source coefficient for 2nd difference of Planckian
                            npd_source_coeff =                           2, & ! Number of coefficients for two-stream sources
                                    ip_solar =                           1, & ! Index for solar loop = 1 
                                ip_infra_red =                           2, & ! Index for IR loop = 2
                                    ip_clear =                           1, & ! Index for loop with clear-sky conditions = 1
                                     ip_dust =                           2, & ! Index for loop with dusty conditions = 2
                                  kdm_do_vis =                           0, & ! flag name:  if (kdm_mode_flag == kdm_do_vis) => we do vis
                                  kdm_do_ir  =                           1, & ! flag name:  if (kdm_mode_flag == kdm_do_ir)  => we do ir
                                  n_band_ir  =                           7, & !c-mm Number of bands in the IR.  Hardwired in KDM.  Originally in module_driver_constants as max_n_band_ir
                                  n_band_solar =                         7    !c-mm Number of bands in the solar.  Hardwired in KDM.  Originally in module_driver_constants as max_n_band_solar

         REAL(KIND(0.d0)), PARAMETER ::                                     & ! THESE VARIABLES SHOULD REMAIN FIXED
                                 tol_machine =                2.220446D-16, & ! Machine tolerance, also equal to EPSILON(x) (intrinsic) where x is virtually any 
                                     tol_div =         3.2D+01*tol_machine, & !    conceivable number used in this code.  Architecture dependent, but this
                            sqrt_tol_machine =                1.490116D-08, & !    is a pretty good number for most.
                             log_tol_machine =                    -36.0437, &
                          quad_correct_limit =                   6.828D-06    ! Equal to EXP(3.3D-01*log_tol_machine)

!                               -------- case-dependent constants --- CHANGE THESE IF YOU WILL, BUT ONLY IF YOU KNOW WHAT YOU'RE DOING!!

         INTEGER, PARAMETER ::                                              &
                                    npd_band = max(n_band_ir,n_band_solar), & ! Basically equals MAX(n_band_ir,n_band_solar).  Used for setting up arrays.
                                 n_tot_bands =      n_band_ir+n_band_solar    ! Number of bands in k-coefficient array

#if ( WRF_MARS == 1 )
         INTEGER, PARAMETER ::                                              & ! THESE VARIABLES CAN BE MODIFIED
                                 npd_species =                           2, & ! Number of "gases" in atmospere.  One for each variable gas plus one total for fixed gases
                         npd_aerosol_species =                           3    ! Max number of "aerosols" in the atmosphere.  Presently set to three. 1=dust, 2=ice, 3=co2ice
#elif ( WRF_TITAN == 1 )
         INTEGER, PARAMETER ::                                              & ! THESE VARIABLES CAN BE MODIFIED
                                 npd_species =                           2, & ! Number of "gases" in atmospere.  One for each variable gas plus one total for fixed gases
                         npd_aerosol_species =                           1    ! Max number of "aerosols" in the atmosphere.  Presently set to one. 1=tholin
#else
         INTEGER, PARAMETER ::                                              & ! THESE VARIABLES CAN BE MODIFIED
                                 npd_species =                           1, & ! Number of "gases" in atmospere.  One for each variable gas plus one total for fixed gases
                         npd_aerosol_species =                           1    ! Max number of "aerosols" in the atmosphere.  Presently set to three. 1=dust, 2=ice, 3=co2ice
#endif

         REAL(KIND(0.d0)), PARAMETER ::                                     & ! THESE VARIABLES SHOULD REMAIN FIXED
#if ( WRF_TITAN == 1 )
                                   p_spacing =                      0.3333, &
                                 mix_spacing =                         1.0, &
#elif ( WRF_MARS == 1 )
                                   p_spacing =                         0.4, &
                                 mix_spacing =                         1.0, &
                                       zs_lv =                        -5.5, & ! Altitude of the NLTE transition region in terms of ALOG(p), where p is in nbar
                                       zw_lv =                         0.5, & ! Width of the NLTE transition region.  The above value of zs_lv corresponds to
                                                                              !    approximately 82 km for reference atmosphere 1 in Lopez-Valverde.
#else
                                   p_spacing =                      0.3333, &
                                 mix_spacing =                         1.0, &
#endif
                                pressure_TOA =                       1.d-8, & ! top of atmosphere pressure (in Pa) (machine small but non-zero) 10^-8 is about 300km on Mars
#if ( WRF_MARS == 1 )
                               pressure_AWDC =                       5.d-4, & ! pressure at which to stop adding fake layers (awdc - above which dont care)
                                  fake_scale =                        0.85, & ! 0.85 is exp(-dz/H) where (e.g., for Mars 0.85 maps to dz=1.5km for H=10km)
#elif ( WRF_TITAN == 1 )
                               pressure_AWDC =                       2.d-3, & ! pressure at which to stop adding fake layers (awdc - above which dont care)
                                  fake_scale =                        0.5,  & ! 0.5 is exp(-dz/H) where (e.g., for Titan 0.5 maps to dz=15km for H=30km)
#else
                               pressure_AWDC =                       2.d-4, & ! pressure at which to stop adding fake layers (awdc - above which dont care)
                                  fake_scale =                        0.5,  & ! 0.5 is exp(-dz/H) where (e.g., for Titan 0.5 maps to dz=15km for H=30km)
#endif
                                     math_pi =                        pi_d


      INTEGER, PARAMETER :: max_fake_layers_visir                = 50

      LOGICAL            ::     verbose_debug = .false.

!                               -------- case-dependent constants that are set in kdminit

      REAL ::                                       T_ref, &   !  these values are now set in kdminit
                                                    T_iso, &   !  and pulled in from 1d state variables
                                               lapse_rate

      REAL(KIND(0.d0)), DIMENSION(n_band_solar+1) ::    band_wn_solar
      REAL(KIND(0.d0)), DIMENSION(n_band_ir+1)    ::    band_wn_ir

      INTEGER :: n_temps,                  & ! Number of temperatures in k-coefficient array
                 n_press,                  & ! Number of pressures in k-coefficient array
                 n_mix,                    & ! Number of mixing ratios in k-coefficient array
                 n_terms                     ! Number of gaussian quadrature terms in each wl band

      REAL(KIND(0.d0)), ALLOCATABLE, DIMENSION(:) :: kval_temps, &   ! temperature values in kvals array
                                                     kval_press, &   ! pressure values in kvals array
                                                     kval_mix,   &   ! mixing ratio values in kvals array
                                                     weights         ! gaussian quadrature weights
      REAL(KIND(0.d0))                            :: kval_press_min, &   ! min pressure value in kvals array
                                                     kval_press_max, &   ! max pressure value in kvals array
                                                     kval_mix_min,   &   ! min mixing ratio value in kvals array
                                                     kval_mix_max,   &   ! max mixing ratio value in kvals array
                                                     kval_temps_min, &   ! min temperature in kvals array
                                                     kval_temps_max      ! max temperature in kvals array


#if ( WRF_MARS == 1 )
      integer, parameter ::                            ko_debug  =  0, &
                                               ko_multimie_dust  =  1, &
                                              ko_multimie_water  =  2, &
                                                ko_multimie_co2  =  3, &
                                                     ko_dummy04  =  4, &
                                                     ko_dummy05  =  5, &
                                                     ko_dummy06  =  6, &
                                                     ko_dummy07  =  7, &
                                                     ko_dummy08  =  8, &
                                                     ko_dummy09  =  9, &
                                                     ko_dummy10  = 10, &
                                                     ko_dummy11  = 11, &
                                                     ko_dummy12  = 12, &
                                                     ko_dummy13  = 13, &
                                                     ko_dummy14  = 14, &
                                                     ko_dummy15  = 15
#endif

      LOGICAL, PARAMETER :: gas_mix_fatal = .false.

      END MODULE module_kdm_parameters

!-----------------------------------------------------------------------
! This is essentially the Hadley/Unified Model radiation code modified !
!    to use a correlated-k scheme.  The subroutine planetary_kdm is in !
!    many ways a "harness" which sets up the variables needed in the   !
!    radiation code.  Subroutine flux_calc is the nominal start of the !
!    Hadley scheme.  Note that it is called twice per implementation,  !
!    once for the solar, and once for the IR.                          !
! There will be a loop over 'j', or latitude, doing latitude strips    !
!    one at a time.                                                    !
! In this code, n_layer is at the surface, and 0 is above the top      !
!    layer.                                                            !
!-----------------------------------------------------------------------
   MODULE module_ra_kdm

      USE module_wrf_error
      USE module_kdm_parameters
      USE module_planet_utilities, ONLY:  rsat_h2o
#if ( WRF_MARS == 1 )
      USE module_state_description, ONLY : TES_LIMB_AVG
      USE module_ra_valverde
      USE module_state_description, ONLY: VALVERDE_LOOKUP,             &
                                          VALVERDE_FIT
#endif

      IMPLICIT NONE

      ! Make sure the declared header common variables don't interfere
      ! with similarly-named variables in routines that USE this module
      ! by making all variables PRIVATE and then making only the
      ! necessary subroutines PUBLICally available
      PRIVATE
      PUBLIC :: planetary_kdm, kdminit,                  &
#if ( WRF_MARS == 1 )
                aerosol_obs_scaling, build_ko_options,   &
                aerosol_obs_scaling_1particle,           &
                reference_opacity, reference_extinction, &
#endif
                solar_fraction

      REAL(KIND(0.d0)), ALLOCATABLE, DIMENSION(:,:,:,:,:) ::           &
!           DIMENSION(n_temps,n_press,n_terms,n_tot_bands,n_mix) ::     & ! Full array of k-coefficients
                                                           kval_array

      REAL(KIND(0.d0)), DIMENSION(n_band_solar) ::    solar_flux_band, &
                                                 rayleigh_coefficient

#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), DIMENSION(n_band_solar) ::   kCext_dust_solar, &
                                                     omega_dust_solar, &
                                                 asymmetry_dust_solar, &
                                                      kCext_ice_solar, &
                                                      omega_ice_solar, &
                                                  asymmetry_ice_solar, &
                                                   kCext_co2ice_solar, &
                                                   omega_co2ice_solar, &
                                               asymmetry_co2ice_solar
                                         
      REAL(KIND(0.d0)), DIMENSION(n_band_ir) ::                        &
                                                        kCext_dust_IR, &
                                                        omega_dust_IR, &
                                                    asymmetry_dust_IR, &
                                                         kCext_ice_IR, &
                                                         omega_ice_IR, &
                                                     asymmetry_ice_IR, &
                                                      kCext_co2ice_IR, &
                                                      omega_co2ice_IR, &
                                                  asymmetry_co2ice_IR
#endif

#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)) ::                              kCext_dust_ref, &
                                                    kCext_dust_ref_ir, & ! reference at dust IR diagnostic wavelength
                                                        kCext_ice_ref, &
                                                     kCext_co2ice_ref
#endif

      INTEGER ::                           n_fake_layers_needed_visir  
      REAL(KIND(0.d0)), DIMENSION(max_fake_layers_visir) ::            &
                                                    p_fake_lower_edge, &
                                                        p_fake_centre, &
                                                    t_fake_lower_edge, &
                                                        t_fake_centre

      LOGICAL ::                                             l_double

      LOGICAL ::                                           kval_debug = .false.

#if ( WRF_MARS == 1 )
! Chicago water ice cloud layers
   REAL, SAVE :: chi_cloud_prs, chi_ls1, chi_ls2, chi_lat1, chi_lat2, chi_lon1, chi_lon2, chi_ice_tau
   REAL, SAVE :: chi_tval1, chi_tval2
   INTEGER, SAVE :: chi_ncloud, chi_dayornight
 
 !multiple particle size mie tables
     real(kind(0.d0)), save, dimension(:,:,:,:), allocatable ::     aerosol_mie
     real(kind(0.d0)), save, dimension(:,:,:), allocatable ::     aerosol_Qref
     real(kind(0.d0)), save, dimension(:), allocatable :: Qref
     real(kind(0.d0)) :: reffstart, reffstop, reffstep
     integer :: nbins_mie

     integer :: kdm_options
     logical :: is_two_moment_dust, is_two_moment_water_ice
     integer :: n_dust_data, n_ice_data, n_co2ice_data ! Number of rows in the various ice and dust parameter data files
#endif
#if ( WRF_TITAN == 1 )
     integer, parameter :: n_tholin_rad_bins = 20
     real, dimension(n_tholin_rad_bins) :: tholin_rad_bins
     real(kind(0.d0)), save, dimension(n_band_solar,n_tholin_rad_bins) ::  &
                                                     tholin_vis_ext_xsect, &
                                                           tholin_vis_ssa, &
                                                          tholin_vis_asym

     real(kind(0.d0)), save, dimension(n_band_ir,   n_tholin_rad_bins) ::  &
                                                      tholin_ir_ext_xsect, &
                                                            tholin_ir_ssa, &
                                                           tholin_ir_asym
#endif

   CONTAINS
!=======================================================================
      SUBROUTINE planetary_kdm(rthratensw, rthratenlw,                 &
                        gsw, glw,                                      &
                        toa_sw_u, toa_sw_d, toa_lw_u, toa_lw_d,        &
                        flux_sw_u,flux_sw_d,flux_lw_u,flux_lw_d,       &
                        xlat, xlong, albedo,                           &
                        emiss, ha, coszen, angslope, azmslope,         &
                        dhr, diurnavg, rho_phy, t3d, tf3d, tsk, qv3d,  &
                        qc3d, qr3d, qi3d, qs3d, qg3d, p3d, pf3d, pi3d, &
                        dz8w, julday,                                  &
#ifdef mpas
                        p2si,                                          &
#endif
                        hr_g_vis, hr_a_vis, hr_vis,                    &
                        hr_g_ir,  hr_a_ir,  hr_ir,                     &
                        sunfrac, julian, declin, solcon,               &
                        l_s, hrang,                                    &
                        f_qv, f_qc, f_qr, f_qi, f_qs, f_qg,            & ! flags for whether moisture variables are present in WRF
                        i_2stream_solar, i_2stream_ir,                 &
                        nlte_physics,                                  &
                        sw_correct, ra_du_physics,                     &
                        do_fake_layers_in,                             &
                        kdm_mode_flag,                                 &
                        t_sounding, polar, p_scale,                    &
                        verbose_debug_in,nan_detail_check_in,          &
#if ( WRF_MARS == 1 )
                        dust_array, cloud_array,co2_cloud_array,       &
                        dust_array_ir,                                 &
                        do_chicago_ice_layers,                         &
                        f_qst_co2ice,                                  &
                        include_dust, include_water_ice,               &
                        include_water_vap, include_co2_ice,            &
                        rhod_h2o,rhod_co2,h2oice_radius,co2ice_radius, &
                        qst_co2ice,                                    &
                        reff_dust, ndust,                              &
                        reff_ice, nice,                                &
                        cdod_scale,                                    &
#endif
#if ( WRF_TITAN == 1 )
                        mck_haze_numden, mck_haze_radius,              &
                        mck_haze_tau64,                                &
                        include_tholin, mp_physics_tholin,             &
                        q_ch4_3d, f_q_ch4, include_ch4_vap,            &
#endif
                        ra_kdm_gray_k_adhoc_ir,                        & 
                        ra_kdm_gray_k_adhoc_vis,                       &
                        ids, ide, jds, jde, kds, kde,                  & 
                        ims, ime, jms, jme, kms, kme,                  &
                        its, ite, jts, jte, kts, kte,                  & 
                        ra_kdm_taulayer_calculation,                   &
                        ra_kdm_clearsky_calculation,                   &
                        swdown, swdowndir,                             &
                        swdownunscat                                   &
                        )

      USE module_ra_utils
#ifdef mpas
      USE mpas_atmphys_constants,only: pa2bar, p0, rcp
#else
      USE module_model_constants
#endif
#if ( WRF_MARS == 1 )
      USE module_ra_mars_common, only: cp_mars
#endif
      USE module_planet_utilities

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::          ids, ide, jds, jde, kds, kde, &
                                         ims, ime, jms, jme, kms, kme, &
                                         its, ite, jts, jte, kts, kte, &
                                                      i_2stream_solar, & ! Two-stream method for solar wavelengths
                                                         i_2stream_ir, & ! Two-stream method for IR wavelengths
                                                               julday, &
                                                         nlte_physics, &
                                                        ra_du_physics

      LOGICAL, INTENT(IN   ) ::                            sw_correct, &
#if ( WRF_MARS == 1 )
                                                         include_dust, &
                                                    include_water_ice, &
                                                    include_water_vap, &
                                                      include_co2_ice, &
#endif
#if ( WRF_TITAN == 1 )
                                                       include_tholin, &
                                                      include_ch4_vap, &
#endif
                                                    do_fake_layers_in

      LOGICAL, OPTIONAL ::     f_qv,f_qc,f_qr,f_qi,f_qs,f_qg
#if ( WRF_MARS == 1 )
      LOGICAL, OPTIONAL ::                      f_qst_co2ice
#endif
#if ( WRF_TITAN == 1 )
      LOGICAL, OPTIONAL ::                           f_q_ch4
#endif

      LOGICAL, INTENT(IN   ) ::                              diurnavg
      logical, intent(in)    ::           ra_kdm_taulayer_calculation, &
                                          ra_kdm_clearsky_calculation

      INTEGER, INTENT(IN   ) ::                         kdm_mode_flag  ! 0=vis, 1=IR see kdm_do_vis and kdm_do_ir definitions

      integer, intent(in   ) ::                            t_sounding  ! method for temperature sounding, default = 2 (constant lapse rate), also used by fake layers
      logical, intent(in   ) ::                                 polar  ! gcm (.true.) or les (.false.) modeling, also used by fake layers
      real,    intent(in   ) ::                               p_scale  ! pressure scalar, default = 1., also used by fake layers

#if ( WRF_MARS == 1 )
      logical, intent(in   ) ::                 do_chicago_ice_layers  ! use water ice opacity information from presc_aero block of the namelist
#endif
#if ( WRF_TITAN == 1 )
      integer, intent(in   ) ::                     mp_physics_tholin
#endif
                                                 
      REAL(KIND(0D0)), INTENT(IN   ) ::                           dhr

 
      REAL, DIMENSION(n_band_solar), INTENT(IN    ) ::                    &
                                         ra_kdm_gray_k_adhoc_vis

      REAL, DIMENSION(n_band_ir), INTENT(IN    ) ::                    &
                                         ra_kdm_gray_k_adhoc_ir

      REAL, INTENT(IN    ) ::                                  julian, &
                                                               declin, &
                                                               solcon, &
                                                                  l_s
#ifdef mpas
      REAL, INTENT(IN    ) ::                                    p2si
#endif
#if ( WRF_MARS == 1 )
      REAL, INTENT(IN) ::      rhod_h2o,rhod_co2,h2oice_radius,co2ice_radius
#endif
      
      LOGICAL, INTENT(IN) ::  verbose_debug_in, nan_detail_check_in


      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) ::            &
                                                                 xlat, &
                                                                xlong, &
                                                               albedo, & ! Surface albedo
                                                                emiss, & ! Surface emissivity--generally pegged at 1.0 for IR
                                                             angslope, &
                                                             azmslope, &
                                                              sunfrac, &
                                                                  tsk, &
                                                                  ha,  &
                                                               coszen, &
                                                                hrang

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) ::   &
                                                                  p3d, &
                                                                 pf3d, &
                                                                 pi3d, &
                                                                  t3d, &
                                                                 tf3d, &
                                                              rho_phy, &
                                                                 dz8w

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL,           &
                                                    INTENT(IN   ) ::   &
                                                                 qv3d, &
                                                                 qc3d, &
                                                                 qr3d, &
                                                                 qi3d, &
                                                                 qs3d, &
                                                                 qg3d
#if ( WRF_MARS == 1 )
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL,           &
                                                    INTENT(IN   ) ::   &
                                                           qst_co2ice

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) ::   &
                                                           dust_array, &  ! dust_array is definitionally VIS (0.67micron) for MarsWRF
                                                          cloud_array, &
                                                      co2_cloud_array
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::   &
                                                        dust_array_ir     ! dust_array_ir is non-trivial for 2-moment (diagnostic)
#endif
#if ( WRF_TITAN == 1 )
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT(IN   ) ::   &
                                                      mck_haze_numden, &
                                                      mck_haze_radius
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT(INOUT) ::   &
                                                       mck_haze_tau64
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL,           &
                                                    INTENT(IN   ) ::   &
                                                             q_ch4_3d
#endif

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) ::   &
                                                           rthratensw, &
                                                           rthratenlw, &
                                                            flux_sw_u, &
                                                            flux_sw_d, &
                                                            flux_lw_u, &
                                                            flux_lw_d

      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::            &
                                                                  gsw, &
                                                                  glw, &
                                                             toa_sw_d, & ! Top of Atmosphere Short Wave (TOASW) flux (down)for heat balance analysis.
                                                             toa_sw_u, & ! Top of Atmosphere Short Wave (TOASW) flux (up  ) for heat balance analysis.
                                                             toa_lw_d, & ! Top of Atmosphere Long Wave (TOALW) flux (down) for heat balance analysis.
                                                             toa_lw_u    ! Top of Atmosphere Long Wave (TOALW) flux (up) for heat balance analysis.


      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) ::   &
                                                             hr_g_vis, & ! Heating rate due to gas in solar wavelengths
                                                             hr_a_vis, & ! Heating rate due to aerosol in solar wavelengths
                                                             hr_vis,   & ! Total heating rate in solar wavelengths
                                                              hr_g_ir, & ! Heating rate due to gas in IR wavelengths
                                                              hr_a_ir, & ! Heating rate due to aerosol in IR wavelengths
                                                              hr_ir      ! Total heating rate in IR wavelengths

      REAL, DIMENSION ( ims:ime, jms:jme), INTENT(INOUT), OPTIONAL ::  &
                                                               swdown, & ! sw down at the surface (i.e. not net)
                                                            swdowndir, & ! direct beam sw down at the surface
                                                         swdownunscat    ! unscattered beam sw down at the surface

#if ( WRF_MARS == 1 )
      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), intent(in), optional &
                                        ::              ndust, reff_dust

      REAL, DIMENSION(ims:ime,kms:kme,jms:jme), intent(in), optional &
                                        ::              nice, reff_ice

      REAL, DIMENSION(ims:ime,jms:jme), optional, intent(in)  &
                                        ::              cdod_scale
#endif


!c-mm Local variables

!     These variables are going to take the values:
!                   n_layer_wrf = kte-kts+1
!               n_layer_visir_m = n_layer_wrf + max_fake_layers_visir
      INTEGER            ::                              n_layer_wrf  , &
                                                     n_layer_visir_m    ! for memory, not for looping

      INTEGER :: n_layer_visir_l                                        !  looping variable, needs creating from needed fake layers from init


      ! ADT: Enforced consistent dimensioning of arrays to ensure
      ! correct behavior on single and multiple processors

      REAL(KIND(0.d0)), DIMENSION( ite-its+1, kte-kts+1+max_fake_layers_visir  ) ::   & 
                                                                 qv2d, &
                                                                 qc2d, &
                                                                 qr2d, &
                                                                 qi2d, &
                                                                 qs2d, &
#if ( WRF_MARS == 1 )
                                                              ndust2d, &
                                                          reff_dust2d, &
                                                               nice2d, &
                                                           reff_ice2d, &
#endif
#if ( WRF_TITAN == 1 )
                                                             n_tholin, &
                                                             r_tholin, &
#endif
                                                                 qg2d


      REAL(KIND(0.d0)),                                                &
           DIMENSION( ite-its+1, kte-kts+1+max_fake_layers_visir, npd_species ) ::           & ! Array of gas mixing ratios for all gases in atmosphere
                                                        gas_mix_ratio

      REAL, DIMENSION( kte-kts+1 ) ::                     hr_co2_nlte

!      REAL, DIMENSION( kte-kts+2 ) ::                             t1d, &
!                                                                 pf1d

      REAL(KIND(0.d0)),                                                &
           DIMENSION( ite-its+1, 0:kte-kts+1+max_fake_layers_visir, npd_aerosol_species )::  & ! Array of aerosol opacities for all aerosol species in atmosphere
                                                      aerosol_opacity
#if ( WRF_TITAN == 1 )                                               
      REAL, DIMENSION( ite-its+1, kte-kts+1+max_fake_layers_visir ):: tau_haze_64
#endif
      REAL, DIMENSION( ite-its+1, kte-kts+1+max_fake_layers_visir, jte-jts+1 ) ::        &
                                                           hr_tot_vis, & ! Total (gas+dust+ice) heating rate in atmosphere for solar wavelengths
                                                            hr_tot_ir    ! Total (gas+dust+ice) heating rate in atmosphere for IR wavelengths

      REAL, DIMENSION( ite-its+1 ) ::                          xlat1d, & !-------!
                                                              xlong1d, &     !
                                                                gsw1d, &     !
                                                                 ha1d, &     ! One-dimensional subarrays of above 2-D arrays
                                                           alphaslp1d, &     !
                                                           gammaslp1d, &     !
                                                            sunfrac1d    !-------!

      REAL(KIND(0.d0)), DIMENSION(ite-its+1) ::                        &
                                                            solar_toa, & ! Incoming solar flux at top of the atmosphere
                                                                sec_0, & ! Solar zenith angle
                                                                tsk1d, & ! One-dimensional surface temperature
                                                                alb1d, & ! One-dimensional surface albedo
                                                                 em1d, & ! One-dimensional surface emissivity
                                                               cosfac    ! 

      REAL(KIND(0.d0)), DIMENSION(kte-kts+1) ::        xltecorrection

      REAL(KIND(0.d0)), DIMENSION(ite-its+1,kte-kts+1) ::       pnorm    ! Two-dimensional layer pressure in normal model coordinates

      REAL(KIND(0.d0)), DIMENSION(ite-its+1,kte-kts+1+max_fake_layers_visir) ::   p2d, & ! Two-dimensional layer pressure
                                                                  t2d    ! Two-dimensional layer temperature

#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), DIMENSION(kte-kts+1+max_fake_layers_visir) ::            z_lv, &
                                                             sigma_lv
#endif

      REAL(KIND(0.d0)), DIMENSION(ite-its+1,0:kte-kts+1+max_fake_layers_visir) ::pf2d, & ! Two-dimensional level pressure
                                                                 tf2d    ! Two-dimensional level temperature
 
      REAL(KIND(0.d0)), DIMENSION(ite-its+1,0:kte-kts+1+max_fake_layers_visir,2) ::    &
                                                          flux_direct, & ! Direct fluxes for aerosol/no aerosol
                                                            flux_down, & ! Downward fluxes for aerosol/no aerosol
                                                              flux_up, & ! Upward fluxes for aerosol/no aerosol
                                                     flux_unscattered    ! Unscattered downward fluxes for aerosol/no aerosol

      REAL(KIND(0.d0)), DIMENSION(ite-its+1,kte-kts+1+max_fake_layers_visir) ::  d_mass    ! Mass of an atmospheric layer
      REAL(KIND(0.d0)), DIMENSION(ite-its+1,kte-kts+1+max_fake_layers_visir) ::  d_z       ! Vertical distance thickness of an atmospheric layer

      REAL(KIND(0.d0)), DIMENSION(ite-its+1,kte-kts+1+max_fake_layers_visir,npd_band) ::     &
                                                          trans_layer    

      REAL ::                                              juldayfrac, &
                                                             tsunfrac, &
                                                             cosfac0d, &
                                                                  cpm ! heat capacity of actual air parcel inc. moisture

      REAL(KIND(0.d0)) :: ref_extinction, top_n_potl, t_frost, top_n_kin, avg_weight
      
      INTEGER ::                               i, j, k, ii, jj, kk, l, &
                                                            n_profile, &
                                                              n_layer, &
                                                                  pos, &
                                                         n_to_average, &
                                                             n, j_end

      LOGICAL ::                                     l_rayleigh, bool, &
                                                     include_aerosols, &
                                                        verbose_debug, &
                                                     nan_detail_check, &
                                                              do_warn
 
      REAL(KIND(0.d0)) ::                                 p_top_check


      INTEGER, DIMENSION(2) :: locations
      INTEGER :: i_at_min, i_at_max, k_at_min, k_at_max
 
      LOGICAL, PARAMETER ::      from_profile = .false., & ! take T(z) from profile or (if false, then) from fit to top 5 prognostic layers
                                  use_kinetic = .true.     ! use kinetic temps from top 5 layers or (if false, then) use potential temps
      LOGICAL, PARAMETER ::     l_rayleigh_ir = .false.    ! to treat rayleigh scattering set equal to true
      LOGICAL, PARAMETER ::    l_rayleigh_vis = .true.     ! to treat rayleigh scattering set equal to true
      LOGICAL            ::    do_fake_layers

!      write(wrf_err_message,*) "solcon=",solcon
!      call wrf_message(trim(wrf_err_message))

      do_warn = .false.
      nan_detail_check = nan_detail_check_in
      verbose_debug = verbose_debug_in
      if((wrf_debug_level==50) .or. (wrf_debug_level>= 999)) verbose_debug=.true.
      if((wrf_debug_level==60) .or. (wrf_debug_level>= 999)) kval_debug=.true.
      if(verbose_debug) nan_detail_check = .true.

      include_aerosols = .false.  ! false unless we actually request aerosols via include_dust, etc.
#if ( WRF_TITAN == 1 )
      include_aerosols =  include_tholin
#elif ( WRF_MARS == 1 )
      include_aerosols =  (include_dust).OR.(include_water_ice).OR.(include_co2_ice)
#endif

      do_fake_layers = do_fake_layers_in

      n_layer_wrf     = kte-kts+1
!c-mm    write(0,*) 'zz',ra_kdm_gray_k_adhoc_ir
#if ( WRF_MARS == 1 )
      if((.not.include_dust).or.is_two_moment_dust) then
        dust_array(its:ite,kts:kte,jts:jte) = 0.
        dust_array_ir(its:ite,kts:kte,jts:jte) = 0.
      endif
      if((.not.include_water_ice).or.is_two_moment_water_ice) then
        cloud_array(its:ite,kts:kte,jts:jte) = 0.
      endif
      if(.not.include_co2_ice) then
        co2_cloud_array(its:ite,kts:kte,jts:jte) = 0.
      endif
#endif

      IF(do_fake_layers) THEN
!c-yl: prepare fake layers for KDM: set number of fake layers, and corresponding pressure and temperature
! (in wrf, pf3d at k=kde (not kme!) is the model top - kme might be the top, or there might be pad)
        call set_fake_layers(p_top=minval(pf3d(its:ite,kde,jts:jte)), t_sounding=t_sounding, polar=polar, p_scale=p_scale)
      ELSE
        n_fake_layers_needed_visir = 0
      ENDIF

      p_top_check=MINVAL(pf3d(its:ite,kde,jts:jte))

      n_layer_visir_l = n_layer_wrf + n_fake_layers_needed_visir
      n_layer_visir_m = n_layer_wrf + max_fake_layers_visir      ! for memory, not for looping

      IF(do_fake_layers) THEN
      if(n_layer_visir_l .gt. n_layer_visir_m) then
          WRITE ( wrf_err_message ,*) 'ERROR: KDM has not set up enough total (real+fake) layers in memory, ',&
                  n_layer_visir_l,' > ',n_layer_visir_m
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF

      if(n_fake_layers_needed_visir .le. 0) then
          WRITE ( wrf_err_message ,*) 'ERROR: KDM needs at least one ("fake") layer (for vis and ir) above the model ptop, ',&
                  n_fake_layers_needed_visir
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      ENDIF ! -- end of fake layers bit --

      qv2d(:,:)= 0.d0   ! default all local water types to zero before reading any in over them
      qc2d(:,:)= 0.d0   ! implicit looping is threadsafe as these are local and only defined on tile
      qr2d(:,:)= 0.d0
      qi2d(:,:)= 0.d0
      qs2d(:,:)= 0.d0
      qg2d(:,:)= 0.d0

      gas_mix_ratio(:,:,:) = 0.d0 ! best not to just pick up junk lying around in memory

!c-mm      call wrf_debug(500,'---> Entered SUBROUTINE planetary_kdm')

      ! Typically, in a WRF run, on the processors occupying the edge
      ! of a particular axis:
      !   kte = kde = # of Z points including extra stagger
      !   ite = ide = # of X points including extra stagger
      ! Don't calculate radiation at the extra x-stagger point:

#ifdef mpas
      n_profile=ite-its+1
#else
      n_profile=MIN(ite,ide-1)-its+1
#endif

      ! On the off chance the model is being decomposed in the 
      !    Z-direction, kte may not equal kde.  If kte does equal kde,
      !    then we want one less than the number including the extra 
      !    z-stagger.

      IF((kde-1) /= kte) THEN
          WRITE ( wrf_err_message ,*) 'ERROR: KDM has not been setup to vertically tile, kte= ',kte,' kde= ',kde
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF

#if ( WRF_MARS == 1 )
! local arrays defined to have full extent same as tile, so no need to worry about writing outside of tile:
     ndust2d(:,:) = 0.d0
     reff_dust2d(:,:) = 0.d0
     nice2d(:,:) =0.d0
     reff_ice2d(:,:) = 0.d0
#endif
#if ( WRF_TITAN == 1 )
     n_tholin(:,:) = 0.d0
     r_tholin(:,:) = 1.d-3  ! minimum so scattering doesn't go crazy (in microns)
#endif

!c-mm  I've determined the best place to calculate the cloud_array values is here, before we
!c-mm     start looping.  By doing it here, we have all the variables we need within the
!c-mm     module, and don't have to go passing it around.  qi3d and/or qst_co2ice are passed into
!c-mm     this subroutine, which makes it the easiest way to do things.
!c-mm  The data in kCext_ice_ref comes from a source which actually provides Qext (instead of k*Cext), so the
!c-mm     name is a misnomer (although it's correct for the CO2 ice below).  Just need to pass this value
!c-mm     into the subroutine directly.
#if ( WRF_MARS == 1 )
      if((ra_du_physics /= TES_LIMB_AVG) .and. present(f_qi) .and. present(qi3d) .and.(.not.do_chicago_ice_layers)) then
       if(include_water_ice .and. f_qi) then
!c-mm         call wrf_debug(500,'---> doing calc_ice_tau for h2o')
         CALL CALC_ICE_TAU(qi3d, cloud_array, pf3d,         &
                           kCext_ice_ref,                   &
                           H2OICE_RADIUS, RHOD_H2O,         &
                           ids,ide, jds,jde, kds,kde,       &
                           ims,ime, jms,jme, kms,kme,       &
                           its,ite, jts,jte, kts,kte)
       endif
      endif

      if(do_chicago_ice_layers .and. include_water_ice) then
!c-mm          call wrf_debug(500,'---> doing ice_tau_chicago')
          CALL ice_tau_chicago(cloud_array, pf3d,            &
                              l_s, xlat, xlong, hrang,       &
                              ids,ide, jds,jde, kds,kde,     &
                              ims,ime, jms,jme, kms,kme,     &
                              its,ite, jts,jte, kts,kte)
      endif

!c-mm  The 1.e6 multiplier of part_radius is because lmie code (offline) calculates in microns, and CO2ICE_RADIUS is in m.
!c-mm     The value of kCext_ref is essentially in microns.  Since we ultimately want Qext, which is kCext_ref/pi*r^2
!c-mm     we need to convert the units of part_radius to coincide with those of kCext_ref.  
!c-mm   Unlike for water ice above, the data in the CO2 ice data file is truly k*Cext, so we need to calculate Qext from
!c-mm     it and the geometric cross section of the ice particles (pi*r^2)
      if(present(f_qst_co2ice) .and. present(qst_co2ice)) then
       if(include_co2_ice .and. f_qst_co2ice) then
!c-mm        call wrf_debug(500,'---> doing calc_ice_tau for co2')
        CALL CALC_ICE_TAU(qst_co2ice, co2_cloud_array, pf3d, &
                        kCext_co2ice_ref,                    &
                        CO2ICE_RADIUS, RHOD_CO2,             &
                        ids,ide, jds,jde, kds,kde,           &
                        ims,ime, jms,jme, kms,kme,           &
                        its,ite, jts,jte, kts,kte)
       endif
      endif
#endif

#ifdef mpas
      j_end=jte
#else
      j_end=MIN(jte,jde-1)
#endif
!c-mm      call wrf_debug(500,"--- about to start j loop")
      j_loop: DO j=jts,j_end
         ! Convert to lower bound of 1 (local) from lower bound of
         ! jts (WRF distributed memory, or serial with jts /= 1)
         jj = j-jts+1

IF (kval_debug) then
  write(wrf_err_message,*) '--- ---> j = ',j,' out of ',MIN(jte,jde-1)-jts+1
  call wrf_message(trim(wrf_err_message))
ENDIF

      XLTECORRECTION = 1.0

!c-mm  Within the KDM code, the index 'i' is the local loop counter 
!c-mm     (only for those variables that are local to KDM).  
!c-mm     The variable 'ii' is the global loop counter, used for those 
!c-mm     variables that are passed into or out of KDM.
!c-mm  **NOTE**, however, that the nomenclature is opposite for the j 
!c-mm     values. The index 'j' is the global loop counter and 'jj' is 
!c-mm     the local loop counter.  I have no idea why the hell I did 
!c-mm     it this way...

         DO i=1,n_profile
            ! Convert from lower bound of 1 (local) to lower bound of
            ! its (WRF distributed memory, or serial with its /= 1)
            ii = i+its-1
            ! Variables on the z-stagger

! Now setup the vis and ir (main kdm) layer structure including 'fake' layers to treat the atmosphere above the dynamical p_top

            ! The KDM model layer structure and how it relates to the wrf layer structure:
            ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !
            ! For the main vis and ir kdm (i.e. not the LTE correction)
            !
            !  ===== is a layer edge
            !  ----- is a layer centre
            !
            !                                                     KDM layer                   <==>       WRF layer
            !
            !    space aka p=0     k=0                            =========                 
            !     layer centre     k=1                            - - - - -
            !     layer edge       k=1                            =========
            !        ....
            !     layer centre     k=n_fake_layer_needed_visir-1  - - - - -
            !     layer edge       k=n_fake_layer_needed_visir-1  ========= 
            !     layer centre     k=n_fake_layer_needed_visir    - - - - -
            !     layer edge       k=n_fake_layer_needed_visir    =========          0                   =========  k=kte  model top AKA p_top
            !     layer centre     k=n_fake_layer_needed_visir+1  - - - - -          1                   - - - - -  k=kte-1
            !     layer edge       k=n_fake_layer_needed_visir+1  =========          1                   =========  k=kte-1
            !     layer centre     k=n_fake_layer_needed_visir+2  - - - - -          2                   - - - - -  k=kte-2
            !        ....
            !     layer centre     k=n_fake_layer_needed_visir-1  - - - - -                              - - - - -  k=kts+1
            !     layer edge       k=n_fake_layer_needed_visir-1  =========                              =========  k=kts+1
            !     layer centre     k=n_fake_layer_needed_visir    - - - - -                              - - - - -  k=kts  first layer centre
            !     surface          k=n_layer_visir_l              =========     n_layer_wrf              =========  k=kts  model bottom AKA psurf
            !
            !                        (above ascii arrow, these are KDM "real" layers /\ )

! figure out the simple (naive) average potential temperature in top n layers - where n = min( ! 1+nint(real(n_fake_layers_needed_visir)/4) , 5 )
!   -- this formula has no deep meaning, it's arbitrary and just defined to grab a "reasonable" number of layers.  It's max'ed to 5
!      so as not to grab a crazy number of layers below the dynamical model top.

            IF(do_fake_layers) THEN
            n_to_average = min(1+int(real(n_fake_layers_needed_visir)/4),5)
            top_n_potl = 0.d0
            top_n_kin  = 0.d0
            avg_weight = 1.d0/real(n_to_average,kind(0.d0))
            do kk=kde-1,kde-n_to_average,-1  ! this should be a loop from top down to n below top in wrf-space (layer centres)
              top_n_potl =  top_n_potl + REAL( (avg_weight * tf3d(ii,kk,j) * ( p0 / pf3d(ii,kk,j))**rcp) ,KIND(0.d0))
              top_n_kin  =  top_n_kin  + avg_weight * tf3d(ii,kk,j)
            enddo

            if(verbose_debug .and. (ii==(its+ite)/2) .and. (j==jts)) then
              write(0,*) i,j,' number of points in top average=',n_to_average
              write(0,*) '  p_top=',minval(pf3d(its:ite,kde,jts:jte)),' pressure_TOA=',pressure_TOA
              write(0,*) '  top n potl, kin',top_n_potl,top_n_kin
              write(0,*) '  at p_top, this is',real(top_n_potl)*((pf3d(ii,kde,j)/p0)**rcp)
            endif

! setup fake vis+ir layer bottom edges
            do k=1,n_fake_layers_needed_visir
                pf2d(i,k)  = p_fake_lower_edge(1+n_fake_layers_needed_visir - k)
                if(from_profile) then
                   tf2d(i,k) = t_fake_lower_edge(1+n_fake_layers_needed_visir - k)
                 else
                   if(use_kinetic) then
                     tf2d(i,k)  = top_n_kin
                   else 
                     tf2d(i,k)  = top_n_potl*((pf2d(i,k)/real(p0,kind(0.d0)))**real(rcp,kind(0.d0)))
                   endif
                endif
                aerosol_opacity(i,k,:)=0.
#if ( WRF_MARS == 1 )
                t_frost = real(tfrost_co2(real(pf2d(i,k))),kind(0.d0))  !  function from module_planet_utilities
                if(tf2d(i,k) .lt. t_frost) tf2d(i,k) = t_frost
#endif
            enddo
            ENDIF  ! -- fake layer bit --

! Real layer edge info for the vis and ir
            DO k=n_fake_layers_needed_visir+1,n_layer_visir_l   ! this is a loop from the top->down for KDM
               kk=kde-(k-n_fake_layers_needed_visir)
               tf2d(i,k)  = REAL(tf3d(ii,kk,j),KIND(0.d0))
               pf2d(i,k)  = REAL(pf3d(ii,kk,j),KIND(0.d0))               ! pf is pressure at the full levels (p8w in the regular code)
               aerosol_opacity(i,k,:)=0.                                 ! Initialize aerosol_opacity to zero everywhere (i.e. 'clear-sky' case)

#if ( WRF_MARS == 1 )
!c-mm  So, what's the intent here?  If we are including dust, then we want the internal KDM routine
!c-mm     to receive the opacity values (from dust_array, which is passed in from radiation_driver)
!c-mm     and operate accordingly.  If we are *not* including dust, then we want to overwrite the
!c-mm     dust_array values with zeroes at this stage.  Just above, we initialize aerosol_opacity to
!c-mm     zero, so the KDM routine will run with aerosol_opacity set to zero anyway.  But, by setting
!c-mm     dust_array to zero here, it passes it back into radiation_driver.  Then, when we get to the
!c-mm     diagnose_tau routines, the values of dust_array are zero, and it spits out zeros in the wrfout
!c-mm     files.  Otherwise, the KDM code will assume tau=0, but the wrfout output will still reflect
!c-mm     the values of tau that were calculated in the du_phys routines.
               IF (include_dust) THEN
                  aerosol_opacity(i,k,1)=                                 &
                                     REAL(dust_array(ii,kk,j),KIND(0.d0))   ! Take only a 2-D slice of dust_array
               ! can diagnose the IR (9.3micron) optical depth from the VIS (0.67micron) as we know that the optical depth
               ! is DEFINED as VIS for non-two-moment cases:
                  dust_array_ir(ii,kk,j) = dust_array(ii,kk,j) * kCext_dust_ref_ir/kCext_dust_ref
               ELSE
                  dust_array(ii,kk,j)=0.0
                  dust_array_ir(ii,kk,j)=0.0
               ENDIF
!c-mm  Same as above, but for cloud_array.  The values in cloud_array are either determined from the tes_limb
!c-mm     routine (if du_phys=49) or else they're calculated above in this routine (from qi3d).  Either way
!c-mm     if water ice is not being used, we want to zero this out in cloud_array so radiation_driver properly
!c-mm     sees zeroes.
               aerosol_opacity(i,k,2)=0.d0
               IF (include_water_ice) THEN  ! can't check f_qi here since TES_LIMB can also set water ice opac
                  aerosol_opacity(i,k,2)=                                 &
                                     REAL(cloud_array(ii,kk,j),KIND(0.d0))  ! Take only a 2-D slice of cloud_array
               ELSE
                  cloud_array(ii,kk,j)=0.0
               ENDIF
!c-mm  Same as above, but for co2_cloud_array.  The values in co2_cloud_array are always calculated above, in this
!c-mm     routine (from qst_co2ice).  Make sure they get returned as zeroes if we are not considering co2 ice.
               aerosol_opacity(i,k,3)=0.d0
               IF (include_co2_ice .and. present(qst_co2ice) .and. present(f_qst_co2ice)) THEN
                if(f_qst_co2ice) then
                  aerosol_opacity(i,k,3)=                                 &
                                     REAL(co2_cloud_array(ii,kk,j),KIND(0.d0))  ! Take only a 2-D slice of cloud_array
                else
                  co2_cloud_array(ii,kk,j)=0.0
                endif
               ELSE
                  co2_cloud_array(ii,kk,j)=0.0
               ENDIF
#endif
            ENDDO

!c-mm            call wrf_debug(500,"-- just about to do fake layers")

            IF(do_fake_layers) THEN
#if ( WRF_MARS == 1 )
! don't need checks of f_qv, f_qst_co2ice, etc., below, since if everything is not 100% kosher, the respective _array = 0.
              IF (include_dust) THEN     ! fill fake layers with top layer aerosol so that we don't get a big 
                    aerosol_opacity(i,0:n_fake_layers_needed_visir,1)=   & ! "transparency" jump at top of prognostic model
                                REAL(dust_array(ii,kde-1,j),KIND(0.d0))
              ENDIF
              IF (include_water_ice) THEN                    
                    aerosol_opacity(i,0:n_fake_layers_needed_visir,2)=   &
                                REAL(cloud_array(ii,kde-1,j),KIND(0.d0))
              ENDIF
              IF (include_co2_ice) THEN
                    aerosol_opacity(i,0:n_fake_layers_needed_visir,3)=   &
                                REAL(co2_cloud_array(ii,kde-1,j),KIND(0.d0))
              ENDIF
#endif
            ENDIF  ! -- end fake layers bit --

            ! Set Top Of Atmosphere values i.e. at top of a fake layer for the vis and ir
            k=0  ! k=0 is the top of the atmosphere - no layer centre information is stored here, but the edge info (pf and tf) is stored here
            tf2d(i,k)  = tf2d(i,1)
            pf2d(i,k)  = pressure_TOA                     ! set TOA at vey small number (approx 0 Pa)
            aerosol_opacity(i,k,:)=0.                     ! sat TOA opacity to zero for all aerosols

!c-mm            aerosol_opacity(i,:,1)=aerosol_opacity(i,:,1)*0.05  !c-mm  Hardwire an 95% reduction in opacity for testing, 12-30-17
! Varaibles in the middle of a layer
            DO k = 1,n_layer_wrf
               pnorm(i,k) = REAL(p3d(ii,k,j),KIND(0.d0))          ! pnorm, used in the NLTE code, has the same orientation as WRF, opposite of KDM
            ENDDO

            DO k=1,n_layer_visir_l

              kk= kts + (n_layer_visir_l - k)

              if (k.le.n_fake_layers_needed_visir) then

                  p2d(i,k) = p_fake_centre(n_fake_layers_needed_visir - (k-1))
                  if(from_profile) then
                    t2d(i,k) = t_fake_centre(n_fake_layers_needed_visir - (k-1))
                  else
                    if(use_kinetic) then
                      t2d(i,k) = top_n_kin
                    else
                      t2d(i,k) = top_n_potl*((p2d(i,k)/real(p0,kind(0.d0)))**real(rcp,kind(0.d0)))
                    endif
                  endif
#if ( WRF_MARS == 1 )
                  t_frost = real(tfrost_co2(real(p2d(i,k))),kind(0.d0))  !  function from module_planet_utilities
                  if(tf2d(i,k) .lt. t_frost) t2d(i,k) = t_frost
#endif
                  qv2d(i,k) = 0.d0   ! default all local water types to zero before reading any in over them
                  qc2d(i,k) = 0.d0
                  qr2d(i,k) = 0.d0
                  qi2d(i,k) = 0.d0
                  qs2d(i,k) = 0.d0
                  qg2d(i,k) = 0.d0

              else
                 p2d(i,k)     = REAL(p3d(ii,kk,j),KIND(0.d0))
                 t2d(i,k)     = REAL(t3d(ii,kk,j),KIND(0.d0))

                 if(qv3d(ii,kk,j) .gt. 0.1) write(0,*) "qv3d bad at ",ii,kk,k,qv3d(ii,kk,j)

                 qv2d(i,k) = 0.d0   ! default all local water types to zero before reading any in over them
                 qc2d(i,k) = 0.d0
                 qr2d(i,k) = 0.d0
                 qi2d(i,k) = 0.d0
                 qs2d(i,k) = 0.d0
                 qg2d(i,k) = 0.d0

#if ( WRF_MARS == 1 )              
                 IF (present(f_qv) .and. (f_qv) .and. include_water_vap) qv2d(i,k)= REAL(qv3d(ii,kk,j),KIND(0.d0))
                 IF (present(f_qi) .and. (f_qi) .and. include_water_ice .and. (.not.do_chicago_ice_layers)) &
                                                                         qi2d(i,k)= REAL(qi3d(ii,kk,j),KIND(0.d0))
#endif

! Mars KDM is not setup to handle liquid water clouds, rain, snow, nor grouple, if it is in
! the future, the following can be uncommented:
!                 IF (present(f_qc) .and. (f_qc)) qc2d(i,k)= REAL(qc3d(ii,kk,j),KIND(0.d0))
!                 IF (present(f_qr) .and. (f_qr)) qr2d(i,k)= REAL(qr3d(ii,kk,j),KIND(0.d0))
!                 IF (present(f_qs) .and. (f_qs)) qs2d(i,k)= REAL(qs3d(ii,kk,j),KIND(0.d0))
!                 IF (present(f_qg) .and. (f_qg)) qg2d(i,k)= REAL(qg3d(ii,kk,j),KIND(0.d0))

                 if(qv2d(i,k) .gt. 0.1) write(0,*) "qv2d immediately bad at ",i,k,qv2d(i,k)
              endif

#if ( WRF_MARS == 1 )
! mass mixing ratio incoming is mmr (mass), but need nmr (number) to lookup k-values: mmr = nmr * molec_mass_species / molec_mass_air
              IF(present(f_qv)) THEN
               IF((f_qv) .and. include_water_vap) then  ! check if we actually have any water vapour
                  gas_mix_ratio(i,k,1)=qv2d(i,k) * real(mwair_o_mwh2o,KIND(0.d0))     ! Already swapped in above statement.
                  IF(.not.gas_mix_fatal) THEN
                    gas_mix_ratio(i,k,1)=MIN(gas_mix_ratio(i,k,1),kval_mix(n_mix))      ! So that code doesn't choke on max mmr lookup
                  ENDIF
               ENDIF
              ENDIF
              gas_mix_ratio(i,k,2)=1.d0-gas_mix_ratio(i,k,1)
#endif
              d_mass(i,k) = REAL((pf2d(i,k)-pf2d(i,k-1))/g,KIND(0.d0))  ! Calculation of d_mass from p_level
              d_z(i,k)    = (REAL(r_d,KIND(0.d0))*t2d(i,k)/REAL(g,KIND(0.d0))) * &
                            log(pf2d(i,k)/pf2d(i,k-1))

#if ( WRF_TITAN == 1 )
              if (k.gt.n_fake_layers_needed_visir) then
                n_tholin(i,k) = REAL(mck_haze_numden(ii,kk,j),KIND(0.d0)) * d_z(i,k) !-> !scale from density (#/m3) to (#/m2)
                r_tholin(i,k) = REAL(mck_haze_radius(ii,kk,j),KIND(0.d0))
                if(present(f_q_ch4) .and. f_q_ch4 .and. include_ch4_vap) then
! mass mixing ratio incoming is mmr (mass), but need nmr (number) to lookup k-values: mmr = nmr * molec_mass_species / molec_mass_air
                  gas_mix_ratio(i,k,1)=REAL(q_ch4_3d(ii,kk,j),KIND(0.d0)) * REAL( mwair_o_mwch4 , kind(0.d0))
                  IF(.not.gas_mix_fatal) THEN
                    gas_mix_ratio(i,k,1)=MIN(gas_mix_ratio(i,k,1),kval_mix(n_mix)) ! So that code doesn't choke on max mmr lookup
                  ENDIF
                endif
              endif
              gas_mix_ratio(i,k,2)=1.d0-gas_mix_ratio(i,k,1)

#endif
#if ( WRF_MARS == 1 )              
              IF ( is_two_moment_dust ) THEN  ! here we construct *diagnostic* dust_array field from two-moment variables
!c-mm              call wrf_debug(500,"...constructing dust_array for 2mom")
              if(btest(kdm_options, ko_multimie_dust).and.include_dust) then
                 if (k.gt.n_fake_layers_needed_visir) then
!c-mm                   call wrf_debug(500,"...doing dust_array after btest")
                   if(present(ndust)) ndust2d(i,k) = REAL(ndust(ii,kk,j),KIND(0.d0)) * d_mass(i,k) !-> !scale ndust from density to number
                   if(present(cdod_scale)) ndust2d(i,k) = cdod_scale(ii,j) * ndust2d(i,k)
                   if(present(reff_dust)) reff_dust2d(i,k) = REAL(reff_dust(ii,kk,j),KIND(0.d0))
                   !update dust_array with the new optical_depth
!c-mm                   call wrf_debug(500,"...doing ref extinction")
                   call reference_extinction(reff_dust2d(i,k),ref_extinction,aerosol_number=1,reference_band=1) ! generate VIS tau for dust
!c-mm                   call wrf_debug(500,"...doing dust_array after ref extinction")
                   dust_array(ii,kk,j) = ndust2d(i,k)*ref_extinction
!c-mm                   call wrf_debug(500,"...doing dust_array_ir after ref extinction")
                   call reference_extinction(reff_dust2d(i,k),ref_extinction,aerosol_number=1,reference_band=2) ! generate IR tau for dust
                   dust_array_ir(ii,kk,j) = ndust2d(i,k)*ref_extinction
                        !if(j==1) write(0,*) 'dust_array: ',ii,j,k,kk,ndust2d(i,k), ref_extinction, dust_array(ii,kk,j)
                 else
                   if(present(ndust)) ndust2d(i,k) = REAL(ndust(ii,kde-1,j),KIND(0.d0)) * d_mass(i,k) !-> !scale ndust from density to number
                   if(present(cdod_scale)) ndust2d(i,k) = cdod_scale(i,j) * ndust2d(i,k)
                   if(present(reff_dust)) reff_dust2d(i,k) = REAL(reff_dust(ii,kde-1,j),KIND(0.d0))
                !   write(*,*) 'fake: ',ii,j,k,kk,n_fake_layers_needed_visir - (k-1),ndust2d(i,k), ref_extinction, dust_array(ii,kk,j)
                 endif
!c-mm                 write(*,*) ii,j,k,kk,ndust2d(i,k), reff_dust2d(i,k), dust_array(ii,kk,j)
              endif
              ENDIF

              IF ( is_two_moment_water_ice ) THEN  ! here we construct *diagnostic* cloud_array field from two-moment variables
              if(btest(kdm_options, ko_multimie_water).and.include_water_ice) then
                 if (k.gt.n_fake_layers_needed_visir) then
                   !-2m ice-
                   if(present(nice)) nice2d(i,k) = REAL(nice(ii,kk,j),KIND(0.d0)) * d_mass(i,k) !-> !scale ndust from density to number
                   if(present(reff_ice)) reff_ice2d(i,k) = REAL(reff_ice(ii,kk,j),KIND(0.d0))
                   call reference_extinction(reff_ice2d(i,k),ref_extinction,aerosol_number=2,reference_band=3) ! IR 12 microns tau for water ice
                   cloud_array(ii,kk,j) = nice2d(i,k)*ref_extinction
                 else
                   if(present(nice)) nice2d(i,k) = REAL(nice(ii,kde-1,j),KIND(0.d0)) * d_mass(i,k) !-> !scale ndust from density to number
                   if(present(reff_ice)) reff_ice2d(i,k) = REAL(reff_ice(ii,kde-1,j),KIND(0.d0))
                 endif
              endif
              ENDIF  ! end of is_two_moment_water_ice check
#endif

            ENDDO ! end of loop in the vertical 1 => n_layer_visir_l

!c-mm            call wrf_debug(500,"just about to finish up the optical depths")

#if ( WRF_MARS == 1 )              
            IF ( is_two_moment_dust ) THEN
            !finish integrating opacity
            if(btest(kdm_options, ko_multimie_dust).and.include_dust) then 
                do k=kte-1,kts,-1
                    dust_array(ii,k,j) = dust_array(ii,k,j) + dust_array(ii,k+1,j)
                    dust_array_ir(ii,k,j) = dust_array_ir(ii,k,j) + dust_array_ir(ii,k+1,j)
!c-mm                    write(0,*) 'dust_array ',k,dust_array_ir(ii,k,j)
                 end do
            endif
            ENDIF
            IF ( is_two_moment_water_ice ) THEN
            if(btest(kdm_options, ko_multimie_water).and.include_water_ice) then 
                do k=kte-1,kts,-1
                    cloud_array(ii,k,j) = cloud_array(ii,k,j) + cloud_array(ii,k+1,j)
                end do
            endif
            ENDIF
!c-mm            if(btest(kdm_options, ko_multimie_co2ice)) then 
!c-mm                do k=kte-1,kts,-1
!c-mm                    co2_cloud_array(ii,k,j) = co2_cloud_array(ii,k,j) + co2_cloud_array(ii,k+1,j)
!c-mm                end do
!c-mm            endif
#endif

!c-mm            call wrf_debug(500,"now loading up location, surface properties, surface temps, etc.")
            tsk1d(i) = REAL(tsk(ii,j),KIND(0.d0))
            xlat1d(i) = xlat(ii,j)
            xlong1d(i) = xlong(ii,j)
            alb1d(i) = REAL(albedo(ii,j),KIND(0.d0))
            em1d(i) = REAL(emiss(ii,j),KIND(0.d0))
! Albedo and Emissivity are no longer modified here because of surface ice.
! Modifications happen in the microphysics schemes and are applied to 
! the albedo/emissivity fields. 
! ALBBCK and EMBCK retain the 'ice free' surface properties to allow the 
! surface to return to its original state.
! Data Assimilation requires the actual Albedo/Emissivity to perform radiance forward
! models correctly.
!           IF (xlat1d(i) > 0.) alb1d(i) = 0.12
            ha1d(i)       = ha(ii,j)
            alphaslp1d(i) = angslope(ii,j)
            gammaslp1d(i) = azmslope(ii,j)
            juldayfrac=julian-INT(julian)

          if(verbose_debug .and. (ii==(its+ite)/2) .and. (j==jts)) then
            write(0,*) 'column i=',i,' j=',j
          IF(do_fake_layers) THEN
            write(0,*) 'fake layers, needed=',n_fake_layers_needed_visir
            write(0,*) '  p_centre=',p_fake_centre(1:n_fake_layers_needed_visir)
            write(0,*) '  t_centre=',t_fake_centre(1:n_fake_layers_needed_visir)
            write(0,*) '  p_lower =',p_fake_lower_edge(1:n_fake_layers_needed_visir)
            write(0,*) '  t_lower =',t_fake_lower_edge(1:n_fake_layers_needed_visir)
          ENDIF
            write(0,*) ' '
            write(0,*) 'n_layer_visir_l',n_layer_visir_l,' kts,kte=',kts,kte,n_layer_wrf,' needed_visir=',n_fake_layers_needed_visir
            write(0,*) '  ----'
            write(0,*) 'p3d: ',p3d(ii,1:n_layer_wrf,j)
            write(0,*) '  ----'
            write(0,*) 'p2d: ',p2d(i,1:n_layer_visir_l)
            write(0,*) '  ----'
            write(0,*) 't3d: ',t3d(ii,1:n_layer_wrf,j)
            write(0,*) '  ----'
            write(0,*) 't2d: ',t2d(i,1:n_layer_visir_l)
            write(0,*) '  ----'
            write(0,*) 'pf3d:',pf3d(ii,1:n_layer_wrf+1,j)
            write(0,*) '  ----'
            write(0,*) 'pf2d:',pf2d(i,0:n_layer_visir_l)
            write(0,*) '  ----'
            write(0,*) 'tf3d:',tf3d(ii,1:n_layer_wrf+1,j)
            write(0,*) '  ----'
            write(0,*) 'tf2d:',tf2d(i,0:n_layer_visir_l)
            write(0,*) '  ----'
#if ( WRF_MARS == 1 )
            write(0,*) 'vap nmr:',gas_mix_ratio(i,1:n_layer_visir_l,1)
            write(0,*) '  ----'
            write(0,*) 'co2 nmr:',gas_mix_ratio(i,1:n_layer_visir_l,2)
#endif
#if ( WRF_TITAN == 1 )
            write(0,*) 'ch4 nmr:',gas_mix_ratio(i,1:n_layer_visir_l,1)
            write(0,*) '  ----'
            write(0,*) 'n2 nmr:',gas_mix_ratio(i,1:n_layer_visir_l,2)
#endif
            write(0,*) '  ----'
            write(0,*) 'total gas mass: ',d_mass(i,1:n_layer_visir_l)
            write(0,*) '  ----'
          endif

IF(kdm_mode_flag == kdm_do_vis) THEN
!c-mm  call wrf_debug(500,"...doing solar band KDM")

            tsunfrac = sunfrac(ii,j)
            cosfac(i)=1.                                                 ! Initialize cosfac to have no impact on gsw later on unless modified below
            solar_toa(i) = REAL(solcon*tsunfrac,KIND(0.d0))              ! KDM code requires double precision variables, so start converting input values
            sec_0(i) = REAL(1./(coszen(ii,j)+tol_machine),KIND(0.d0))
ENDIF ! end of the local vis check

         ENDDO  ! this is the end of the 'i' loop

         call kdm_debug('--- ---> Done set up',kval_debug)

          !  Check if pressure is out of bounds of the table.
          IF ((ANY(p2d(:,1:n_layer_visir_l) < kval_press_min)).OR.(ANY(p2d(:,1:n_layer_visir_l) > kval_press_max))) THEN
            locations=MINLOC(p2d(:,1:n_layer_visir_l))
            i_at_min = locations(1)+its-1
            k_at_min = locations(2)
            locations=MAXLOC(p2d(:,1:n_layer_visir_l))
            i_at_max = locations(1)+its-1
            k_at_max = locations(2)
            WRITE ( wrf_err_message ,*) 'KDM: table min and max pressure:',kval_press_min,kval_press_max
            call wrf_message ( wrf_err_message )
            WRITE ( wrf_err_message ,*) 'KDM: Pressure out of bounds.',MINVAL(p2d(:,1:n_layer_visir_l)), &
                                                                       MAXVAL(p2d(:,1:n_layer_visir_l)), &
                                                                 " loc min i,k,j: ",i_at_min,k_at_min,j, &
                                                                 " loc max i,k,j: ",i_at_max,k_at_max,j
            if(verbose_debug) then
                write(0,*) "pressures in column with minimum:"
                write(0,*) p2d(i_at_min,1:n_layer_visir_l)
                write(0,*) "pressures in column with maximum:"
                write(0,*) p2d(i_at_max,1:n_layer_visir_l)
            endif
            CALL wrf_error_fatal ( wrf_err_message )
          ENDIF

#if ( WRF_MARS == 1 )
          !  Check if water vapor mixing ratio exeeds top of table.
          if(include_water_vap) then
          IF (ANY(gas_mix_ratio(:,1:n_layer_visir_l,1) > kval_mix_max)) THEN
            locations=MAXLOC(gas_mix_ratio(:,1:n_layer_visir_l,1))
            i_at_max = locations(1)+its-1
            k_at_max = locations(2)
            if(verbose_debug) then
              WRITE ( wrf_err_message ,*) "Water vap column with biggest error: ",gas_mix_ratio(i_at_max,1:n_layer_visir_l,1)
              CALL wrf_message ( TRIM( wrf_err_message ) )
            endif
            if(gas_mix_fatal) then
              WRITE ( wrf_err_message ,*) 'Water vapor mixing ratio greater than ',& 
              kval_mix_max,' max at i,k,j: ',i_at_max,k_at_max,j
              CALL wrf_error_fatal ( wrf_err_message )
            endif
          ENDIF
          IF((.not.gas_mix_fatal) .and. do_warn .and. (ANY(gas_mix_ratio(:,1:n_layer_visir_l,1) .eq. kval_mix_max))) THEN
            WRITE ( wrf_err_message ,*) 'Water vap mmr restricted to ',kval_mix_max,' in KDM, but likely exceeds this value (WARNING!)'
            CALL wrf_message ( TRIM( wrf_err_message ) )
          ENDIF
          ENDIF
#endif
#if ( WRF_TITAN == 1 )
          !  Check if methane vapor mixing ratio exeeds top of table.
          if(include_ch4_vap) then
          IF (ANY(gas_mix_ratio(:,1:n_layer_visir_l,1) > kval_mix_max)) THEN
            locations=MAXLOC(gas_mix_ratio(:,1:n_layer_visir_l,1))
            i_at_max = locations(1)+its-1
            k_at_max = locations(2)
            if(verbose_debug) then
              WRITE ( wrf_err_message ,*) "Methane vap column with biggest error: ",gas_mix_ratio(i_at_max,1:n_layer_visir_l,1)
              CALL wrf_message ( TRIM( wrf_err_message ) )
            endif
            if(gas_mix_fatal) then
              WRITE ( wrf_err_message ,*) 'Methane vapor mixing ratio greater than ktable max at i,k,j, max mix: ',i_at_max,k_at_max,j,kval_mix_max
              CALL wrf_error_fatal ( wrf_err_message )
            endif
          ENDIF
          IF((.not.gas_mix_fatal) .and. do_warn .and. (ANY(gas_mix_ratio(:,1:n_layer_visir_l,1) .eq. kval_mix_max))) THEN
            WRITE ( wrf_err_message ,*) 'Methane vap mmr set to K-table max in KDM, but likely exceeds this value (WARNING!)', &
                                    kval_mix_max
            CALL wrf_message ( TRIM( wrf_err_message ) )
          ENDIF
          ENDIF
#endif
          !  Check if temperature is out of bounds of the table.
          IF ((ANY(t2d(:,1:n_layer_visir_l) < kval_temps_min)).OR.(ANY(t2d(:,1:n_layer_visir_l) > kval_temps_max))) THEN
            locations=MINLOC(t2d(:,1:n_layer_visir_l))
            i_at_min = locations(1)+its-1
            k_at_min = locations(2)
            locations=MAXLOC(t2d(:,1:n_layer_visir_l))
            i_at_max = locations(1)+its-1
            k_at_max = locations(2)
            write(wrf_err_message,*) "T bounds on kval-table in KDM: ",kval_temps_min,kval_temps_max
            CALL wrf_message ( TRIM( wrf_err_message ) )
            if(verbose_debug) then
              write(wrf_err_message,*) "T with min: ",t2d(i_at_min,1:n_layer_visir_l)
              CALL wrf_message ( TRIM( wrf_err_message ) )
              write(wrf_err_message,*) "P with min: ",p2d(i_at_min,1:n_layer_visir_l)
              CALL wrf_message ( TRIM( wrf_err_message ) )
              write(wrf_err_message,*) "T with max: ",t2d(i_at_max,1:n_layer_visir_l)
              CALL wrf_message ( TRIM( wrf_err_message ) )
              write(wrf_err_message,*) "P with max: ",p2d(i_at_max,1:n_layer_visir_l)
              CALL wrf_message ( TRIM( wrf_err_message ) )
            endif
            WRITE ( wrf_err_message ,*) 'Temperature out of bounds for kval table', &
                                                                     MINVAL(t2d(:,1:n_layer_visir_l)), &
                                                                     MAXVAL(t2d(:,1:n_layer_visir_l)), &
                                                               " loc min i,k,j: ",i_at_min,k_at_min,j, &
                                                               " loc max i,k,j: ",i_at_max,k_at_max,j
            CALL wrf_error_fatal ( wrf_err_message )
          ENDIF

IF (kdm_mode_flag == kdm_do_vis) THEN

          call kdm_debug('--- ---> going to flux_calc for solar',kval_debug)

!          write(0,*) "before solar flux_calc call, l_rayleigh=",l_rayleigh_vis

          CALL flux_calc(l_rayleigh_vis, i_2stream_solar,              & ! First call to flux_calc is for the shortwave
                         n_band_solar,                                 &
                         n_profile, n_layer_visir_l,                   & !   and it contains the full argument list
                         p2d(:,1:n_layer_visir_l),                     &
                         t2d(:,1:n_layer_visir_l),                     &     
                         tsk1d,                                        &
                         tf2d(:,0:n_layer_visir_l),                    & 
                         d_mass(:,1:n_layer_visir_l),                  &
                         ip_solar,                                     &
                         gas_mix_ratio(:,1:n_layer_visir_l,:),         &
#if ( WRF_MARS == 1 )
                         include_dust,                                 &
                         include_water_ice,                            &
                         include_co2_ice,                              &
#endif
#if ( WRF_TITAN == 1 )
                         include_tholin,                               &
                         n_tholin(:,1:n_layer_visir_l),                &
                         r_tholin(:,1:n_layer_visir_l),                &
                         tau_haze_64,                                  &
#endif
                         aerosol_opacity(:,0:n_layer_visir_l,:),       &
                         alb1d, em1d,                                  &
                         flux_direct(:,0:n_layer_visir_l,:),           &
                         flux_down(:,0:n_layer_visir_l,:),             &
                         flux_up(:,0:n_layer_visir_l,:),               &
                         flux_unscattered(:,0:n_layer_visir_l,:),      &
                         band_wn_ir, weights,                          &
                         trans_layer(:,1:n_layer_visir_l,:),           &
                         j,                                            & ! j is not used and my guess is passed in just for debugging
                         ra_kdm_gray_k_adhoc_ir,                       & 
                         ra_kdm_gray_k_adhoc_vis,                      &
                         sec_0=sec_0,                                  & ! including (not needed in ir):  sec_0, 
                         solar_toa=solar_toa,                          & !                                solar_toa, 
                         solar_flux_band=solar_flux_band,              & !                                solar_flux_band, 
                         rayleigh_coefficient=rayleigh_coefficient,    & !                                rayleigh_coefficient.
#if ( WRF_MARS == 1 )
                         reff_dust=reff_dust2d(:,1:n_layer_visir_l),   &
                         ndust    =ndust2d(:,1:n_layer_visir_l),       &
                         reff_ice =reff_ice2d(:,1:n_layer_visir_l),    &
                         nice     =nice2d(:,1:n_layer_visir_l),        &
#endif
                         p_top_check=p_top_check,                      &
                         ra_kdm_taulayer_calculation=ra_kdm_taulayer_calculation, &
                         ra_kdm_clearsky_calculation=ra_kdm_clearsky_calculation, &
                         nan_detail_check=nan_detail_check,            &
                         kdm_mode_flag=kdm_mode_flag)

         call kdm_debug('--- ---> Done solar flux_calc',kval_debug)

!c-mm  Within the KDM code, the index 'i' is the local loop counter 
!c-mm     (only for those variables that are local to KDM).  
!c-mm     The variable 'ii' is the global loop counter, used for those 
!c-mm     variables that are passed into or out of KDM.
!c-mm  **NOTE**, however, that the nomenclature is opposite for the j 
!c-mm     values. The index 'j' is the global loop counter and 'jj' is 
!c-mm     the local loop counter.  I have no idea why the hell I did 
!c-mm     it this way...
          DO i=1,n_profile
             ! Convert from lower bound of 1 (local) to lower bound of
             ! its (WRF distributed memory, or serial with its /= 1)
             ii = i+its-1

             IF(include_aerosols) THEN
                if(present(swdown))             swdown(ii,j) = REAL(flux_down(i,n_layer_visir_l,ip_dust))        ! SW down at surface
                if(present(swdowndir))       swdowndir(ii,j) = REAL(flux_direct(i,n_layer_visir_l,ip_dust))      ! component of SW down at surface due to the direct solar beam
                if(present(swdowndir) .and. nan_detail_check)  call test_for_bad(swdowndir(ii,j),"flux_direct",.false.,0.)
                if(present(swdownunscat)) swdownunscat(ii,j) = REAL(flux_unscattered(i,n_layer_visir_l,ip_dust)) ! component of SW down at surface due to the direct solar beam
                gsw(ii,j) = REAL(flux_down(i,n_layer_visir_l,ip_dust)-         & ! Surface values of flux_down/flux_up are at the n_layer_visir_l point (from the viewpoint 
                                 flux_up(i,n_layer_visir_l,ip_dust))             !    of flux_down and flux_up which are returned from KDM code)  
                                                                         ! The value of cosfac found in SWSLOPE accounts for both the solar zenith angle 
                                                                         !    and the surface slope effects for calculating GSW (cf. LMD radiation scheme).
                                                                         !    In KDM scheme, we have already accounted for solar zenith angle in the
                                                                         !    calculations of flux_up and flux_down, so 'cosfac' should only account for
                                                                         !    surface slope alone.  That is what this last COS term (currently 
                                                                         !    commented out) does.
                toa_sw_d(ii,j) = REAL(flux_down(i,0,ip_dust))            ! This is the short wave flux entering at the top of the atmosphere.
                toa_sw_u(ii,j) = REAL(flux_up(i,0,ip_dust))              ! This is the short wave flux leaving at the top of the atmosphere.
             ELSE
                gsw(ii,j) = REAL(flux_down(i,n_layer_visir_l,ip_clear)-        & ! Surface values of flux_down/flux_up are at the n_layer_visir_l point (from the viewpoint 
                                 flux_up(i,n_layer_visir_l,ip_clear))            !    of flux_down and flux_up which are returned from KDM code)   
                toa_sw_d(ii,j) = REAL(flux_down(i,0,ip_clear))              ! This is the short wave flux entering at the top of the atmosphere.
                toa_sw_u(ii,j) = REAL(flux_up(i,0,ip_clear))                ! This is the short wave flux leaving at the top of the atmosphere.
             ENDIF
             flux_sw_u(ii,kte+1,j) = toa_sw_u(ii,j)   !c-mm  These started out with k=0, but I think I needed to invert it.
             flux_sw_d(ii,kte+1,j) = toa_sw_d(ii,j)
#if ( WRF_MARS == 1 )
             IF (sw_correct) THEN
                CALL swcorrection(pnorm(i,:),xltecorrection,kms,kme,   &
                                  kts,kte)
             ENDIF
#endif

             DO k=n_fake_layers_needed_visir+1,n_layer_visir_l     ! k=1 is centre of first fake layer, n_fake_layers_needed_visir+1 first real layer, k runs downward
                                                                   !        in the atmosphere with increasing index value (see ascii art figure, many lines above)
               kk= kts + (n_layer_visir_l - k)                     ! exact same kk defn as when we loaded the t and p arrays, kk runs upward with increasing index value

#if ( WRF_MARS == 1 )
!               cpm = cp_mars(t2d(i,k)) * (1. +((1800./cp_mars(t2d(i,k)))-1.) * max(gas_mix_ratio(i,k,1),0.) )  ! functional form from module_diffusion_em (from Michael Mischna)
               cpm = cp_mars(t2d(i,k))
#elif ( WRF_TITAN == 1 )
               cpm = cp
#else
               call wrf_error_fatal("KDM: model configured for unrecognized planet, cannot set cp (you probably have other issues, too")
#endif

#if ( WRF_TITAN == 1 )
                 mck_haze_tau64(ii,kk,j) = tau_haze_64(i,k)
#endif

               IF(include_aerosols) THEN
                   hr_tot_vis(i,k,jj) = REAL(((flux_down(i,k,ip_dust)- & 
                                           flux_up(i,k,ip_dust))-      & 
                                           (flux_down(i,k-1,ip_dust)-  & 
                                           flux_up(i,k-1,ip_dust)))*   &
                                           g/(pf2d(i,k-1)-pf2d(i,k))/  &
                                           cpm)
                   hr_tot_vis(i,k,jj) = hr_tot_vis(i,k,jj)/xltecorrection(kk)

                   flux_sw_u(ii,kk,j) = flux_up(i,k,ip_dust)     ! flux_sw is a WRF index array (hence kk), flux_up is a KDM index array (hence k)
                   flux_sw_d(ii,kk,j) = flux_down(i,k,ip_dust)
                   
                   if(ra_kdm_clearsky_calculation) then
                   hr_g_vis(ii,kk,j) = REAL(((flux_down(i,k,ip_clear)- & ! WRF works (by default) in single precision, so convert back from double precision
                                       flux_up(i,k,ip_clear))-         & ! This code is still safe if WRF is compiled with double precision as default real
                                       (flux_down(i,k-1,ip_clear)-     & ! Fluxes are KDM oriented, hr needs to be WRF oriented, so flip.
                                       flux_up(i,k-1,ip_clear)))*      &
                                       g/(pf2d(i,k-1)-pf2d(i,k))/      &
                                       cpm)

                   hr_g_vis(ii,kk,j) = hr_g_vis(ii,kk,j)/xltecorrection(kk)
                   hr_a_vis(ii,kk,j) = hr_tot_vis(i,k,jj)-hr_g_vis(ii,kk,j)
                   else
                      hr_g_vis(ii,kk,j) = 0.
                      hr_a_vis(ii,kk,j) = 0.
                   endif
                ELSE
                   hr_tot_vis(i,k,jj) = REAL(((flux_down(i,k,ip_clear)-    & ! WRF works (by default) in single precision, so convert back from double precision
                                        flux_up(i,k,ip_clear))-         & ! This code is still safe if WRF is compiled with double precision as default real
                                        (flux_down(i,k-1,ip_clear)-     & ! Fluxes are KDM oriented, hr needs to be WRF oriented, so flip.
                                        flux_up(i,k-1,ip_clear)))*      &
                                        (DBLE(g)/(pf2d(i,k-1)-pf2d(i,k))) /      &
                                        DBLE(cpm))

                   hr_tot_vis(i,k,jj) = hr_tot_vis(i,k,jj)/              &
                                                      xltecorrection(kk)   ! Vertical orientation of NLTE code is same as WRF, so no need to flip.
                   hr_a_vis(ii,kk,j) = 0.0
                   hr_g_vis(ii,kk,j) = hr_tot_vis(i,k,jj)-hr_a_vis(ii,kk,j)

                   flux_sw_u(ii,kk,j) = flux_up(i,k,ip_clear)
                   flux_sw_d(ii,kk,j) = flux_down(i,k,ip_clear)
                ENDIF
                rthratensw(ii,kk,j)=hr_tot_vis(i,k,jj)/pi3d(ii,kk,j)
                hr_vis(ii,kk,j)    =hr_tot_vis(i,k,jj)
             ENDDO
IF(verbose_debug) THEN
                if((ii==(its+ite)/2) .and. (j==jts)) then
                        write(0,*) "vis down raw flux on full (real + fake grid):"
                        write(0,*) flux_down(i,:,ip_clear)
                        write(0,*) "vis toa values, up and down:",toa_sw_u(ii,j),toa_sw_d(ii,j)
                        write(0,*) "vis flux down:",flux_sw_d(ii,kts:kte-1,j)
                        write(0,*) "vis rthraten: ",rthratensw(ii,kts:kte-1,j)
                        write(0,*) "gsw: ",gsw(ii,j)
                endif
ENDIF
          ENDDO

          call kdm_debug('--- ---> Done visible heating rates',kval_debug)

ENDIF ! this is the end of the vis check -------------------------------------

IF(kdm_mode_flag == kdm_do_ir) THEN

          call kdm_debug('--- ---> going to flux_calc for IR',kval_debug)

          CALL flux_calc(l_rayleigh_ir, i_2stream_ir,                  &
                         n_band_ir, n_profile,                         & ! Second call to flux_calc is for the longwave, and only contains a partial argument list
                         n_layer_visir_l,                              & !    the remainder of the variables are defined as 'optional' in flux_calc and can therefore
                         p2d(:,1:n_layer_visir_l),                     & !    be excluded.  This is always a tricky syntax, so be aware of what's happening....
                         t2d(:,1:n_layer_visir_l),                     &     
                         tsk1d,                                        &
                         tf2d(:,0:n_layer_visir_l),                    & 
                         d_mass(:,1:n_layer_visir_l),                  &
                         ip_infra_red,                                 &
                         gas_mix_ratio(:,1:n_layer_visir_l,:),         &
#if ( WRF_MARS == 1 )
                         include_dust,                                 &
                         include_water_ice,                            &
                         include_co2_ice,                              &
#endif
#if ( WRF_TITAN == 1 )
                         include_tholin,                               &
                         n_tholin(:,1:n_layer_visir_l),                &
                         r_tholin(:,1:n_layer_visir_l),                &
                         tau_haze_64,                                  &
#endif
                         aerosol_opacity(:,0:n_layer_visir_l,:),       &
                         alb1d, em1d,                                  &
                         flux_direct(:,0:n_layer_visir_l,:),           &
                         flux_down(:,0:n_layer_visir_l,:),             &
                         flux_up(:,0:n_layer_visir_l,:),               &
                         flux_unscattered(:,0:n_layer_visir_l,:),      &
                         band_wn_ir, weights,                          &
                         trans_layer(:,1:n_layer_visir_l,:),           &
                         j,                                            &
                         ra_kdm_gray_k_adhoc_ir,                       & 
                         ra_kdm_gray_k_adhoc_vis,                      &
#if ( WRF_MARS == 1 )
                         reff_dust=reff_dust2d(:,1:n_layer_visir_l),   &
                         ndust    =ndust2d(:,1:n_layer_visir_l),       &
                         reff_ice =reff_ice2d(:,1:n_layer_visir_l),    &
                         nice     =nice2d(:,1:n_layer_visir_l),        &
#endif
                         p_top_check=p_top_check,                      &
                         ra_kdm_taulayer_calculation=ra_kdm_taulayer_calculation, &
                         ra_kdm_clearsky_calculation=ra_kdm_clearsky_calculation, &
                         nan_detail_check=nan_detail_check,            &
                         kdm_mode_flag=kdm_mode_flag)

          call kdm_debug('--- ---> Done IR flux_calc',kval_debug)

          DO i=1,n_profile
             ! Convert from lower bound of 1 (local) to lower bound of
             ! its (WRF distributed memory, or serial with its /= 1)
             ii = i+its-1
             IF(include_aerosols) THEN
                glw(ii,j) = REAL(flux_down(i,n_layer_visir_l,ip_dust))      ! Surface values of flux_down/flux_up are at the n_layer_visir_l point (from the viewpoint of 
                                                                            !    flux_down and flux_up which are returned from KDM code)
                toa_lw_u(ii,j) = REAL(flux_up(i,0,ip_dust))                 ! This is the long wave flux leaving at the top of the atmosphere.
                toa_lw_d(ii,j) = REAL(flux_down(i,0,ip_dust))               ! This is the long wave flux entering at the top of the atmosphere.
             ELSE
                glw(ii,j) = REAL(flux_down(i,n_layer_visir_l,ip_clear))           
                toa_lw_u(ii,j) = REAL(flux_up(i,0,ip_clear))                ! This is the long wave flux leaving at the top of the atmosphere.
                toa_lw_d(ii,j) = REAL(flux_down(i,0,ip_clear))              ! This is the long wave flux entering at the top of the atmosphere.

!                write(0,*) "TOA values (u and d): ",toa_lw_u(ii,j),toa_lw_d(ii,j)
             ENDIF

             flux_lw_u(ii,kte+1,j) = toa_lw_u(ii,j)
             flux_lw_d(ii,kte+1,j) = toa_lw_d(ii,j)

             SELECT CASE (nlte_physics)
#if ( WRF_MARS == 1 )
                CASE (VALVERDE_LOOKUP)
                   CALL VALVERDE_CO2_IR_LOOKUP(hr_co2_nlte,            &
                          p3d(ii,:,j),t3d(ii,:,j),kms,kme,kts,kte)       ! We pass in 3-D versions of p,T rather than 2-D because 2-D version is 'flipped' for
                CASE (VALVERDE_FIT)                                      !    KDM use, 3-D version is not.  NLTE code is not 'flipped', so it's just easier to
                   CALL VALVERDE_CO2_IR_FIT(hr_co2_nlte,p3d(ii,:,j),   & !    do it this way.
                          t3d(ii,:,j),kms,kme,kts,kte )
#endif
                CASE DEFAULT
             END SELECT

             DO k=n_fake_layers_needed_visir+1,n_layer_visir_l    ! we don't need to load the fake layers into the wrf variables
               kk= kts + (n_layer_visir_l - k)  ! exact same kk defn as when we loaded the t and p arrays

#if ( WRF_MARS == 1 )
                ! Combine the LTE and NLTE cooling rates to yield a 
                !    contnuous profile. LTE is valid from 0 to about 
                !    60 km, above that NLTE becomes important.
                ! z_lv is the vertical variable in the original paper
                ! The expression in parenthesis has to be in nanobars
                z_lv(kk) = -log(pf3d(ii,kk,j)*pa2bar*1.E9)           ! kk references GCM vertical coord not the (flipped) KDM, pf3d is in GCM coord
                ! there had been a long-standing bug concern above the above line, the code to the end of the do-loop is from CL's 'devmort' branch at commit 5967930
                ! The blending coefficient (a function of height)
                sigma_lv(kk) = 0.5*(1.+tanh((z_lv(kk)-zs_lv)/zw_lv))

!               cpm = cp_mars(t2d(i,k)) * (1. +((1800./cp_mars(t2d(i,k)))-1.) * max(gas_mix_ratio(i,k,1),0.) )  ! functional form from module_diffusion_em (from Michael Mischna)
               cpm = cp_mars(t2d(i,k))
#elif ( WRF_TITAN == 1 )
               cpm = cp
#else
               call wrf_error_fatal("KDM: model configured for unrecognized planet, cannot set cp (you probably have other issues, too")
#endif

                IF(include_aerosols) THEN
                   hr_tot_ir(i,k,jj) = REAL(((flux_down(i,k,ip_dust)-  &  
                                          flux_up(i,k,ip_dust))-       &
                                          (flux_down(i,k-1,ip_dust)-   &
                                          flux_up(i,k-1,ip_dust)))*    &
                                          g/(pf2d(i,k-1)-pf2d(i,k))/   &
                                          cpm)

                   flux_lw_u(ii,kk,j) = flux_up(i,k,ip_dust)
                   flux_lw_d(ii,kk,j) = flux_down(i,k,ip_dust)
                                  
                   if(ra_kdm_clearsky_calculation) then 
                   hr_g_ir(ii,kk,j) = REAL(((flux_down(i,k,ip_clear)-     & ! Fluxes are KDM oriented, hr needs to be WRF oriented, so flip.
                                          flux_up(i,k,ip_clear))-      &
                                          (flux_down(i,k-1,ip_clear)-   &
                                          flux_up(i,k-1,ip_clear)))*g/ &
                                          (pf2d(i,k-1)-pf2d(i,k))/       &
                                          cpm)
                   else
                    hr_g_ir(ii,kk,j) = 0.
                   endif
                   SELECT CASE (nlte_physics)
#if  ( WRF_MARS == 1 )
                      CASE (VALVERDE_LOOKUP)
                         hr_tot_ir(i,k,jj) = hr_tot_ir(i,k,jj)*(1.-sigma_lv(kk))+            &
                                             hr_co2_nlte(kk)*sigma_lv(kk)
                        if(ra_kdm_clearsky_calculation) then 
                         hr_g_ir(ii,kk,j) = hr_g_ir(ii,kk,j)*(1.-sigma_lv(kk))+            &
                                           hr_co2_nlte(kk)*sigma_lv(kk)
                        endif

                      CASE (VALVERDE_FIT)                                      !    KDM use, 3-D version is not.  NLTE code is not 'flipped', so it's just easier to
                         hr_tot_ir(i,k,jj) = hr_tot_ir(i,k,jj)*(1.-sigma_lv(kk))+            &
                                             hr_co2_nlte(kk)*sigma_lv(kk)
                         if(ra_kdm_clearsky_calculation) then                     
                         hr_g_ir(ii,kk,j) = hr_g_ir(ii,kk,j)*(1.-sigma_lv(kk))+            &
                                           hr_co2_nlte(kk)*sigma_lv(kk)
                         endif
#endif
                      CASE DEFAULT
                   END SELECT
                   if(ra_kdm_clearsky_calculation) then 
                   hr_a_ir(ii,kk,j) = hr_tot_ir(i,k,jj)-hr_g_ir(ii,kk,j)
                   else
                    hr_a_ir(ii,kk,j) = 0.
                   endif

                ELSE
                   hr_tot_ir(i,k,jj) = REAL(((flux_down(i,k,ip_clear)-     & ! Fluxes are KDM oriented, hr needs to be WRF oriented, so flip.
                                        flux_up(i,k,ip_clear))-            &
                                        (flux_down(i,k-1,ip_clear)-        &
                                        flux_up(i,k-1,ip_clear)))*         &
                                        (dble(g)/(pf2d(i,k-1)-pf2d(i,k)))/ &
                                        dble(cpm))

!IF(verbose_debug) THEN
!                   if((ii==(its+ite)/2) .and. (j==jts)) then
!                        write(0,*) k,"hr_tot_ir: ",hr_tot_ir(i,k,jj),flux_down(i,k,ip_clear), &
!                          cpm,pf2d(i,k-1)-pf2d(i,k),g,flux_up(i,k,ip_clear)
!                   endif
!ENDIF

                   flux_lw_u(ii,kk,j) = flux_up(i,k,ip_clear)
                   flux_lw_d(ii,kk,j) = flux_down(i,k,ip_clear)

                   SELECT CASE (nlte_physics)
#if  ( WRF_MARS == 1 )
                      CASE (VALVERDE_LOOKUP)
                         hr_tot_ir(i,k,jj) = hr_tot_ir(i,k,jj)*(1.-sigma_lv(kk))+            &
                                             hr_co2_nlte(kk)*sigma_lv(kk)
                      CASE (VALVERDE_FIT)                                      !    KDM use, 3-D version is not.  NLTE code is not 'flipped', so it's just easier to
                         hr_tot_ir(i,k,jj) = hr_tot_ir(i,k,jj)*(1.-sigma_lv(kk))+            &
                                             hr_co2_nlte(kk)*sigma_lv(kk)
#endif
                      CASE DEFAULT
                   END SELECT
                   hr_a_ir(ii,kk,j) = 0.0
                   hr_g_ir(ii,kk,j) = hr_tot_ir(i,k,jj)-hr_a_ir(ii,kk,j)

IF(verbose_debug) THEN
                if((ii==(its+ite)/2) .and. (j==jts)) then
                   write(0,*) k," terms: ",flux_down(i,k,ip_clear),flux_up(i,k,ip_clear),        &
                                           flux_down(i,k-1,ip_clear),flux_up(i,k-1,ip_clear),   &
                                           (g/((pf2d(i,k-1)-pf2d(i,k))*cpm)), &
                                           hr_tot_ir(i,k,jj)
                endif
ENDIF
                ENDIF
                rthratenlw(ii,kk,j)=hr_tot_ir(i,k,jj)/pi3d(ii,kk,j)
                hr_ir(ii,kk,j)     =hr_tot_ir(i,k,jj)
             ENDDO
IF(verbose_debug) THEN
                if((ii==(its+ite)/2) .and. (j==jts)) then
                        write(0,*) "ir down raw flux on full (real + fake grid):"
                        write(0,*) flux_down(i,:,ip_clear)
                        write(0,*) "ir toa values, up and down:",toa_lw_u(ii,j),toa_lw_d(ii,j)
                        write(0,*) "ir flux down:",flux_lw_d(ii,kts:kte-1,j)
                        write(0,*) "ir rthraten: ",rthratenlw(ii,kts:kte-1,j)
                        write(0,*) "glw: ",glw(ii,j)
                endif
ENDIF
         ENDDO

         call kdm_debug('--- ---> Done IR heating rates',kval_debug)

ENDIF ! end of the IR check -------------------------------

      ENDDO j_loop

      if(nan_detail_check) &
      call do_nan_checking(kdm_mode_flag, rthratensw, rthratenlw, gsw, glw          &
#if ( WRF_MARS == 1 )
                          ,dust_array, cloud_array, co2_cloud_array, dust_array_ir  &
#endif
#if ( WRF_TITAN == 1 )
                          ,mck_haze_tau64                                           &
#endif
                          ,ims, ime, jms, jme, kms, kme,                            &
                           its, ite, jts, jte, kts, kte,                            &
                           swdown, swdowndir, swdownunscat                          )

      call kdm_debug('---> Exited SUBROUTINE planetary_kdm',kval_debug)

      END SUBROUTINE planetary_kdm
      
!-----------------------------------------------------------------------
!+ Subroutine to calculate radiative fluxes.                           !
!                                                                      !
! Method:                                                              !
!       Properties independent of the spectral bands are set.          !
!       A loop over bands is then entered. Gray optical properties     !
!       are set. The final fluxes are assigned.                        !
!-----------------------------------------------------------------------

      SUBROUTINE flux_calc(l_rayleigh, i_2stream, n_band, n_profile,   &
                           n_layer, p, t, t_ground, t_level, d_mass,   &
                           isolir, gas_mix_ratio,                      &
#if ( WRF_MARS == 1 )
                           include_dust,                               &
                           include_water_ice,                          &
                           include_co2_ice,                            &
#endif
#if ( WRF_TITAN == 1 )
                           include_tholin,                             &
                           n_tholin, r_tholin,                         &
                           tau_haze_64,                                &
#endif
                           aerosol_opacity, surface_albedo,            &
                           emissivity_ground, flux_direct, flux_down,  &
                           flux_up, flux_unscattered,                  &
                           band_wn_ir, weights, trans_layer,           &
                           j,                                          & 
                           ra_kdm_gray_k_adhoc_ir,                     & 
                           ra_kdm_gray_k_adhoc_vis,                    &
                           sec_0, solar_toa, solar_flux_band,          &
                           rayleigh_coefficient,                       &
#if ( WRF_MARS == 1 )
                           reff_dust, ndust,                           &
                           reff_ice, nice,                             &
#endif
                           p_top_check,                                &
                           ra_kdm_taulayer_calculation,                &
                           ra_kdm_clearsky_calculation,                &
                           nan_detail_check,                           &
                           kdm_mode_flag                               )  

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                               isolir, &
                                                            i_2stream, &
                                                            n_profile, &
                                                              n_layer, &
                                                               n_band, &
                                                                    j

      REAL(KIND(0.d0)), DIMENSION(n_terms), INTENT(IN   ) ::  weights

      REAL(KIND(0.d0)), DIMENSION(n_band+1), INTENT(IN   ) ::          &
                                                           band_wn_ir    ! Has the extra index

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                    d_mass

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer, npd_species),    &
           INTENT(IN   ) ::                                            &
                                                        gas_mix_ratio

      REAL(KIND(0.d0)),                                                &
           DIMENSION(n_profile, 0:n_layer, npd_aerosol_species),       &
           INTENT(INOUT) ::                           aerosol_opacity

      REAL(KIND(0.d0)), DIMENSION(n_band_solar), INTENT(IN   ),        &
           OPTIONAL ::                                solar_flux_band, & 
                                                 rayleigh_coefficient

      REAL, DIMENSION(n_band_solar), INTENT(IN   ) ::      &
                                         ra_kdm_gray_k_adhoc_vis

      REAL, DIMENSION(n_band_ir), INTENT(IN   ) ::         &
                                         ra_kdm_gray_k_adhoc_ir

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN   ),           &
           OPTIONAL ::                                      solar_toa, &
                                                                sec_0

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN   ) ::         &
                                                       surface_albedo, &
                                                    emissivity_ground, &
                                                             t_ground                                 

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                         p, &
                                                                    t

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer),               &
           INTENT(IN   ) ::                                   t_level

      REAL(KIND(0.d0)), INTENT(IN   ) ::                  p_top_check


      LOGICAL, INTENT(IN   ) ::                            l_rayleigh, &
                                                     nan_detail_check

#if ( WRF_MARS == 1 )
      LOGICAL, INTENT(IN   ) ::                          include_dust, &
                                                    include_water_ice, &
                                                      include_co2_ice
#endif
#if ( WRF_TITAN == 1 )
      LOGICAL, INTENT(IN   ) ::                        include_tholin
#endif
      REAL(KIND(0.d0)), DIMENSION(n_profile,0:n_layer,2),              &
           INTENT(  OUT) ::                               flux_direct, &
                                                            flux_down, &
                                                              flux_up, &
                                                     flux_unscattered

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer, npd_band),       &
           INTENT(  OUT)::                                trans_layer

      INTEGER, INTENT(IN), OPTIONAL ::                  kdm_mode_flag

#if ( WRF_MARS == 1 )
      REAL(kind(0.d0)), DIMENSION(n_profile, n_layer), intent(in), optional &
                                        ::               ndust, reff_dust, &
                                                         nice, reff_ice
#endif
#if ( WRF_TITAN == 1 )
      REAL(kind(0.d0)), DIMENSION(n_profile, n_layer), intent(in) ::   &  
                                                             n_tholin, &
                                                             r_tholin
      REAL, DIMENSION(n_profile, n_layer), intent(inout) :: tau_haze_64
#endif
      LOGICAL, INTENT(IN) ::              ra_kdm_taulayer_calculation, &
                                          ra_kdm_clearsky_calculation
!c-mm Local variables:

      INTEGER ::                                               i_band, & ! Index of which band we're doing
                                                                    i, &
                                                                    k

      REAL(KIND(0.d0)), DIMENSION(n_profile,2*n_layer+2,2) ::          &
                                                           flux_total, & ! total flux
                                                      flux_total_band    ! total flux in band


      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                     diff_planck_band, & ! differential thermal source in band
                                                   diff_planck_band_2, & ! 2 x 2nd diff. thermal source in band
                                                     trans_layer_band

      REAL(KIND(0.d0)), dimension(n_profile, n_layer) ::               &
                                                           k_gray_ext, & ! free total gray extinction
                                                          k_gray_scat, & ! free scattering extinction
                                                        asymmetry_gas, & ! free asymmetries
                                                      forward_scatter    ! free forward scattering function
#if ( WRF_TITAN == 1 )
      REAL(KIND(0.d0)), dimension(n_profile, n_layer) ::               &
                                                         omega_tholin, &
                                                     asymmetry_tholin
#endif
#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), dimension(n_profile, n_layer) ::               &
                                                           omega_dust, &
                                                       asymmetry_dust, &
                                                            omega_ice, &
                                                        asymmetry_ice, &
                                                         omega_co2ice, &
                                                     asymmetry_co2ice
    
      REAL(kind(0.d0)), dimension(n_profile, n_layer) ::               &
                                                        reff_dust_index, &
                                                        reff_ice_index
#endif

      REAL(KIND(0.d0)), DIMENSION(n_profile) ::                        &
                                                       albedo_surface, & ! direct surface albedo
                                                  inc_solar_flux_band, & ! incident solar flux in band
                                                  thermal_ground_band    ! ground source function in band

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer,2) ::           &
                                                     flux_direct_band, &    ! direct flux in band
                                                flux_unscattered_band       ! unscattered direct flux in band

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer) ::             &
                                                          planck_flux, & ! Plankian flux in band
                                                   planck_source_band    ! Planck function in band at levels

      REAL(KIND(0.d0)) ::                  rayleigh_coefficient_local    ! Local variable needed to handle optional variable

      ! Local variable for the adhoc gray calculations.
      REAL(KIND(0.d0)) ::                                  adhoc_gray

   call kdm_debug('---> Entered SUBROUTINE flux_calc',kval_debug)

#if ( WRF_MARS == 1 )
   IF ( is_two_moment_dust ) THEN
    if(btest(kdm_options, ko_multimie_dust).and.include_dust) then
      if (present(reff_dust)) then
        call mie_findex(reff_dust, reff_dust_index)
      endif
    endif
   ENDIF

   IF( is_two_moment_water_ice ) THEN
    if(btest(kdm_options, ko_multimie_water).and.include_water_ice) then
      if (present(reff_ice)) then
        call mie_findex(reff_ice, reff_ice_index)
      endif
    endif
   ENDIF

!c-mm    if(btest(kdm_options, ko_multimie_co2ice)) then
!c-mm      if (present(reff_co2ice)) then
!c-mm        call mie_findex(reff_co2ice, reff_co2ice_index)
!c-mm      endif
!c-mm    endif
#endif

      flux_direct=0.                                                  ! Initialize the total fluxes.
      flux_unscattered=0.
      planck_flux=0.                                                  ! Initialize the Planckian fluxes if required.

      IF (isolir == ip_solar) THEN
         albedo_surface=surface_albedo
      ELSE IF (isolir == ip_infra_red) THEN
         albedo_surface=1.-emissivity_ground
      ENDIF
      flux_total=0.
      DO i_band=1, n_band                                                ! Solve transfer equation in each band & increment flux

         IF (PRESENT(rayleigh_coefficient)) THEN
            rayleigh_coefficient_local=rayleigh_coefficient(i_band)
         ELSE
            rayleigh_coefficient_local=0.
         END IF

         IF (isolir == ip_solar) THEN
             adhoc_gray = DBLE(ra_kdm_gray_k_adhoc_vis(i_band))
         ELSE IF (isolir == ip_infra_red) THEN
	     adhoc_gray = DBLE(ra_kdm_gray_k_adhoc_ir(i_band))
         END IF

         CALL gray_extinction(l_rayleigh, rayleigh_coefficient_local,  &
                              k_gray_ext, k_gray_scat, asymmetry_gas,  & ! Calculate gray extinction within the band
                              forward_scatter,                         &
                              n_profile, n_layer,                      &
                              i_band, isolir, aerosol_opacity,         &
#if ( WRF_TITAN == 1 )
                              n_tholin, r_tholin, include_tholin,      &
                              omega_tholin,   asymmetry_tholin,        &
#endif
#if ( WRF_MARS == 1 )
                              omega_dust,   asymmetry_dust,            &
                              omega_ice,    asymmetry_ice,             &
                              omega_co2ice, asymmetry_co2ice,          &
                              reff_dust_index=reff_dust_index,         &
                              reff_ice_index=reff_ice_index,           &
                              ndust=ndust, nice=nice,                  &
                              include_dust=include_dust,               &
                              include_water_ice=include_water_ice,     &
                              include_co2_ice=include_co2_ice,         &
#endif
                              adhoc_gray=adhoc_gray                    )


#if ( WRF_TITAN == 1 )
         IF( (isolir == ip_solar) .and. (i_band == 4)) THEN
           do i=1,n_profile
             tau_haze_64(i,1) = aerosol_opacity(i,1,1)
             do k=2,n_layer
               tau_haze_64(i,k) = tau_haze_64(i,k-1) + aerosol_opacity(i,k,1)
             enddo
           enddo
         ENDIF
#endif

         IF (isolir == ip_solar) THEN
            IF (PRESENT(solar_flux_band).AND.PRESENT(sec_0).AND.       &
                PRESENT(solar_toa)) THEN
               inc_solar_flux_band=solar_toa*solar_flux_band(i_band)/  &
                                   sec_0                                 ! Convert normalized solar fluxes to energy fluxes.
            ELSE
               inc_solar_flux_band=0.                                    ! Optional variables should be passed when isolir=ip_solar
            END IF                                                       ! but define variable for safety anyway.
            flux_direct_band=0.                                          ! Initialize the flux in the band to zero.
            flux_unscattered_band=0.
         ELSE IF (isolir == ip_infra_red) THEN
            CALL cfst4(band_wn_ir(i_band), band_wn_ir(i_band+1), t,    & ! Calculate the change in the thermal source function 
                       t_level, t_ground, planck_source_band,          & !    across each layer for the infra-red part of the 
                       thermal_ground_band, diff_planck_band,          & !    spectrum.
                       diff_planck_band_2, n_profile, n_layer, j)
         ENDIF

         IF (isolir == ip_solar) THEN                                    ! Set the emissivity here.  This IF-THEN is extraneous
            thermal_ground_band=0.                                       !    as thermal_ground_band is only used to set
         ELSE IF (isolir == ip_infra_red) THEN                           !    source_ground later, which is itself set to zero
            thermal_ground_band=emissivity_ground*thermal_ground_band    !    in the solar bands.
         ENDIF
         flux_total_band=0.

        CALL solve_band_one_gas(n_profile, n_layer, p, t, d_mass,      &
                                i_2stream, i_band, gas_mix_ratio,      &
#if ( WRF_TITAN == 1 )
                                include_tholin,                        &
#endif
#if ( WRF_MARS == 1 )
                                include_dust,                          &
                                include_water_ice,                     &
                                include_co2_ice,                       &
#endif
                                aerosol_opacity, isolir,               &
                                inc_solar_flux_band,                   &
                                planck_source_band(:,0),               &
                                planck_source_band(:,n_layer),         &
                                diff_planck_band, diff_planck_band_2,  &
                                albedo_surface, thermal_ground_band,   &
                                k_gray_ext, k_gray_scat,               &
                                asymmetry_gas, forward_scatter,        &
#if ( WRF_TITAN == 1 )
                                omega_tholin, asymmetry_tholin,        &
#endif
#if ( WRF_MARS == 1 )
                                omega_dust, asymmetry_dust,            &
                                omega_ice, asymmetry_ice,              &
                                omega_co2ice, asymmetry_co2ice,        &
#endif
                                flux_direct_band, flux_total_band,     &
                                flux_unscattered_band,                 &
                                weights, trans_layer_band, j,          &
                                nan_detail_check,                      &
                                sec_0=sec_0,  p_top_check=p_top_check, &
                                ra_kdm_taulayer_calculation=ra_kdm_taulayer_calculation, &
                                ra_kdm_clearsky_calculation=ra_kdm_clearsky_calculation, &
                                kdm_mode_flag=kdm_mode_flag )

         trans_layer(:,:,i_band) = trans_layer_band                      ! Stores transmissivities for all bands, but model can see only one at a time
         if(.not.ieee_is_normal(trans_layer(1,1,i_band))) write(0,*) "trans_layer has gone bad"

         call kdm_debug('---> Entering SUBROUTINE augment_total_flux',kval_debug)
         CALL augment_total_flux(isolir, planck_source_band,           & ! Increment the total fluxes in the band
                                 flux_direct,  flux_total,             & 
                                 flux_unscattered,                     &
                                 flux_direct_band, flux_total_band,    &
                                 flux_unscattered_band,                &
                                 planck_flux, n_profile, n_layer)

      ENDDO

      CALL assign_flux(n_layer, n_profile, flux_total, isolir,         & ! Pass the calculated fluxes into the output arrays
                       planck_flux, flux_down, flux_up,                &
#if ( WRF_MARS == 1 )
                       include_dust,include_water_ice,include_co2_ice, &
#endif
                       j)

      call kdm_debug('---> Exited SUBROUTINE flux_calc',kval_debug)

      END SUBROUTINE flux_calc

!-----------------------------------------------------------------------
!+ Subroutine to calculate gray extinctions.                           !
!                                                                      !
! Method:                                                              !
!       For each activated optical process, excluding gaseous          !
!       absorption, increments are calculated for the total and        !
!       scattering extinctions, these increments are summed, and the   !
!       gray total and scattering extinctions are thus calculated.     !
!       At this point, only rayleigh scattering and basic dust + ice   !
!       scattering are treated here.                                   !
!-----------------------------------------------------------------------
      SUBROUTINE gray_extinction(l_rayleigh, rayleigh_coeff,           &
                                 k_gray_ext, k_gray_scat,              &
                                 asymmetry_gas, forward_scatter,       &
                                 n_profile, n_layer, i_band, isolir,   &
                                 aerosol_opacity,                      &
#if ( WRF_TITAN == 1 )
                                 n_tholin, r_tholin, include_tholin,   &
                                 omega_tholin,   asymmetry_tholin,     &
#endif
#if ( WRF_MARS == 1 )
                                 omega_dust,   asymmetry_dust,         &
                                 omega_ice,    asymmetry_ice,          &
                                 omega_co2ice, asymmetry_co2ice,       &
                                 reff_dust_index, ndust,               &
                                 reff_ice_index, nice,                 &
                                 include_dust,include_water_ice,       &
                                 include_co2_ice,                      &
#endif
                                 adhoc_gray                            )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                            n_profile, &
                                                              n_layer, &
                                                               i_band, &
                                                               isolir
       
      REAL(KIND(0.d0)), INTENT(IN   ) ::               rayleigh_coeff

      LOGICAL, INTENT(IN   ) ::                            l_rayleigh


      REAL(KIND(0.d0)), dimension(n_profile, n_layer), &
                                            INTENT(  OUT) ::           &
                                                          k_gray_scat, &
                                                           k_gray_ext, &
                                                        asymmetry_gas, &
                                                      forward_scatter

#if ( WRF_TITAN == 1 )
      LOGICAL, INTENT(IN) ::                           include_tholin
      real(kind(0.d0)), intent(in), dimension(n_profile, n_layer) ::   &
                                                             n_tholin, &
                                                             r_tholin
      REAL(KIND(0.d0)), dimension(n_profile, n_layer), &
                                            INTENT(  OUT) ::           &
                                                         omega_tholin, &
                                                     asymmetry_tholin
#endif
#if ( WRF_MARS == 1 )
      real(kind(0.d0)), intent(in), optional, &
                            dimension(n_profile, n_layer) ::           &
                                              reff_dust_index, ndust,  &
                                              reff_ice_index, nice

      LOGICAL, INTENT(IN), OPTIONAL ::                   include_dust, &
                                                    include_water_ice, &
                                                      include_co2_ice

      REAL(KIND(0.d0)), dimension(n_profile, n_layer), &
                                            INTENT(  OUT) ::           &
                                                           omega_dust, &
                                                       asymmetry_dust, &
                                                            omega_ice, &
                                                        asymmetry_ice, &
                                                         omega_co2ice, &
                                                     asymmetry_co2ice
#endif

      ! adhoc gray == an additional amount of absorption added to the
      ! gray extinction absorption coefficient variable, k_gray_ext.
      REAL(KIND(0.d0)), INTENT(IN   ), OPTIONAL ::         adhoc_gray


      REAL(KIND(0.d0)),                                                &
           DIMENSION(n_profile, 0:n_layer, npd_aerosol_species),       &
           INTENT(OUT   ) ::                           aerosol_opacity
      real(kind(0.d0)), DIMENSION(n_profile, n_layer) :: coefficient
     
      call kdm_debug('---> Entered SUBROUTINE gray_extinction',kval_debug)

      k_gray_ext(:,:)=0.d0                                               ! Initialize the extinction coefficients and the asymmetry 
      k_gray_scat(:,:)=0.d0                                              !    product.    
      asymmetry_gas(:,:)=0.d0                                            ! Forward scattering only required in visible when delta-
      forward_scatter(:,:)=0.d0                                          !    rescaling is performed.

      IF (l_rayleigh) THEN                                               ! Include rayleigh scattering.
         k_gray_scat(:,:)=k_gray_scat(:,:) + rayleigh_coeff              ! "rayleigh_coeff" will be 0. when l_rayleigh is false,
         k_gray_ext(:,:)=k_gray_ext(:,:) + k_gray_scat(:,:)              ! which also handles the use of the variable on this line
      ENDIF                                                              ! With Rayleigh scattering only, k_gray_ext will equal k_gray_scat
      
      ! The following adds a band-dependent gray absorption to simulate
      ! an unknown species. The adhoc_gray variables is in units of
      ! m^2/kg. (Soto)
      if(present(adhoc_gray)) k_gray_ext(:,:)=k_gray_ext(:,:)+adhoc_gray

      WHERE (k_gray_scat > tol_div)                                     ! This IF-THEN block is a bit extraneous right now since
         asymmetry_gas=asymmetry_gas/k_gray_scat                         !    asymmetry of gas is always 0.0.  This code can
      ELSEWHERE                                                               !    potentially be used for other gray processes, so I'm
         asymmetry_gas=0.                                                !    keeping it here for heritage purposes.  But again
      END WHERE                                                          !    asymmetry_gas=0.0 always.
      forward_scatter=asymmetry_gas*asymmetry_gas                        ! Also an extraneous statement

#if ( WRF_TITAN == 1 )
      call get_tholin_tau_omega_asym(r_tholin,n_tholin,i_band,isolir,             &
                                     n_profile,n_layer,npd_aerosol_species,       &
                                     aerosol_opacity,omega_tholin,asymmetry_tholin)
#endif

#if ( WRF_MARS == 1 ) 
      IF(present(include_dust)) THEN
      IF(include_dust) THEN 
      IF (is_two_moment_dust .and. (btest(kdm_options, ko_multimie_dust))) THEN
         !-dust-
         IF (.NOT.present(reff_dust_index)) THEN
            WRITE(wrf_err_message,*) "This kdm_options options requires reff_dust_index, but it wasn't present.", __LINE__
            CALL debug_ko_options(kdm_options)
            CALL wrf_error_fatal(TRIM(wrf_err_message))
         ENDIF
         !interpolate mie scattering table to fix dust properties
         CALL mie_interpolation(reff_dust_index, i_band, isolir, 1,1, coefficient)
         CALL mie_interpolation(reff_dust_index, i_band, isolir, 1,2, omega_dust) !interpolate OMEGA
         CALL mie_interpolation(reff_dust_index, i_band, isolir, 1,3, asymmetry_dust) !interpolate asymmetry
         !write(*,*) aerosol_opacity(:,1:n_layer,1)
         !opacity is per (mean) particle, scale by number of particles.
         aerosol_opacity(:,1:n_layer,1) = coefficient * ndust
      ELSE
         IF (isolir == ip_solar) THEN
            omega_dust     = omega_dust_solar(i_band)
            asymmetry_dust = asymmetry_dust_solar(i_band)
         ELSE IF (isolir == ip_infra_red) THEN
            omega_dust     = omega_dust_ir(i_band)
            asymmetry_dust = asymmetry_dust_ir(i_band)
         ENDIF
      ENDIF
      ELSE !not include dust
        omega_dust = 0
        asymmetry_dust = 0
      ENDIF
      ENDIF

      IF(present(include_water_ice)) THEN
      IF(include_water_ice) THEN 
      IF (is_two_moment_water_ice .and. (btest(kdm_options, ko_multimie_water))) THEN
         !-ice-
         IF (.NOT.present(reff_ice_index)) THEN
            WRITE(wrf_err_message,*) "This kdm_options options requires reff_ice_index, but it wasn't present.", __LINE__
            CALL debug_ko_options(kdm_options)
            CALL wrf_error_fatal(TRIM(wrf_err_message))
         ENDIF
         !interpolate mie scattering table to fix dust properties
         CALL mie_interpolation(reff_ice_index, i_band, isolir, 2,1, coefficient)
         CALL mie_interpolation(reff_ice_index, i_band, isolir, 2,2, omega_ice) !interpolate OMEGA
         CALL mie_interpolation(reff_ice_index, i_band, isolir, 2,3, asymmetry_ice) !interpolate asymmetry
         !write(*,*) aerosol_opacity(:,1:n_layer,1)
         !opacity is per (mean) particle, scale by number of particles.
         aerosol_opacity(:,1:n_layer,2) = coefficient * nice
         !write(*,*) aerosol_opacity(:,1,2)
      ELSE
         IF (isolir == ip_solar) THEN
            omega_ice      = omega_ice_solar(i_band)
            asymmetry_ice  = asymmetry_ice_solar(i_band)
         ELSE IF (isolir == ip_infra_red) THEN
            omega_ice      = omega_ice_ir(i_band)
            asymmetry_ice  = asymmetry_ice_ir(i_band)
         ENDIF
      ENDIF
      ELSE !not include ice
        omega_ice = 0
        asymmetry_ice = 0
      ENDIF
      ENDIF

!c-mm      IF (btest(kdm_options, ko_multimie_co2ice)) THEN
!c-mm         !-co2ice-
!c-mm         IF (.NOT.present(reff_co2ice_index)) THEN
!c-mm            WRITE(0,*) "This kdm_options options requires reff_co2ice_index, but it wasn't present.", __LINE__
!c-mm            CALL debug_ko_options(kdm_options)
!c-mm            STOP
!c-mm         ENDIF
!c-mm         !interpolate mie scattering table to fix dust properties
!c-mm         CALL mie_interpolation(reff_co2ice_index, i_band, isolir, 3,1, coefficient)
!c-mm         CALL mie_interpolation(reff_co2ice_index, i_band, isolir, 3,2, omega_co2ice) !interpolate OMEGA
!c-mm         CALL mie_interpolation(reff_co2ice_index, i_band, isolir, 3,3, asymmetry_co2ice) !interpolate asymmetry
!c-mm         !opacity is per (mean) particle, scale by number of particles.
!c-mm         aerosol_opacity(:,1:n_layer,3) = coefficient * nco2ice
!c-mm      ELSE
      IF(present(include_co2_ice)) THEN
      IF(include_co2_ice) THEN 
         IF (isolir == ip_solar) THEN
            omega_co2ice      = omega_co2ice_solar(i_band)
            asymmetry_co2ice  = asymmetry_co2ice_solar(i_band)
         ELSE IF (isolir == ip_infra_red) THEN
            omega_co2ice      = omega_co2ice_ir(i_band)
            asymmetry_co2ice  = asymmetry_co2ice_ir(i_band)
         ENDIF
      ELSE !not include co2ice
        omega_co2ice = 0
        asymmetry_co2ice = 0
      ENDIF
      ENDIF
#endif

      call kdm_debug('---> Exited SUBROUTINE gray_extinction',kval_debug)

      END SUBROUTINE gray_extinction

!-----------------------------------------------------------------------
!+ Subroutine calculating the fraction of blackbody radiation          !
!  between two wavelengths, and the difference in blackbody radiation  !
!  between adjacent layers.                                            !
!                                                                      !
! Method:                                                              !
!        Uses series formulae from Houghton 'Physics of Atmospheres'   !
!-----------------------------------------------------------------------
      SUBROUTINE cfst4(wl1, wl2, temp, t_level, t_ground, fst4,        &
                       fst4_ground, diff_planck, diff_planck_2,        &
                       n_profile, n_layer, j)

      IMPLICIT NONE

      REAL(KIND(0.d0)), PARAMETER ::                                   &
                                          c  =         1.53989733D-01, &
                                         c2  =            1.43883D+04

      INTEGER, INTENT(IN   ) ::                                        &
                                                            n_profile, &
                                                              n_layer

      REAL(KIND(0.d0)), INTENT(IN   ) ::                               &
                                                                  wl1, &
                                                                  wl2  

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer),               &
           INTENT(IN   ) ::                                   t_level

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                      temp

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN   ) ::         &
                                                             t_ground                                                

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(  OUT) ::         &
                                                          fst4_ground

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(  OUT) ::                               diff_planck, &
                                                        diff_planck_2

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer),               &
           INTENT(  OUT) ::                                      fst4

!c-mm  Local variables

      REAL(KIND(0.d0)), DIMENSION(2) ::                                &
                                                                   wl, &
                                                                    f

      REAL(KIND(0.d0)) ::                                              &
                                                                 facc, &
                                                                   ff, &
                                                                    w, &
                                                                    v, &
                                                                   wv

      INTEGER ::                                                       &
                                                       i, m, ii, k, j

IF (kval_debug) PRINT*, '---> Entered SUBROUTINE cfst4'

      IF (wl1 == 0.) THEN                                                ! Swapping of wl's going from wavenumber to wavelength.
         wl(2)=1000.
      ELSE
           wl(2)=10000./wl1
      ENDIF
      wl(1)=10000./wl2
      facc=0.00001
      DO ii=1,n_profile
         DO i=1,2
            f(i)=0.
            IF (wl(i) <= 0.) CYCLE 
            v=c2/(wl(i)*t_ground(ii))
            DO m=1,1000
               w=FLOAT(m)
               wv=w*v
               ff=(w**(-4))*(EXP(-wv))*(((wv+3.)*wv+6.)*wv+6.)
               f(i)=f(i)+ff
               IF (f(i) == 0.) CYCLE
               IF (ff/f(i) <= facc) EXIT
            ENDDO
         ENDDO
         fst4_ground(ii)=c*(f(2)-f(1))*5.67D-08*t_ground(ii)*          &
                         t_ground(ii)*t_ground(ii)*t_ground(ii)
!         fst4_ground(ii)=planck_function(REAL(t_ground(ii)))            ! Gets fluxes in 15um band only, a la LMD model
         DO k=0,n_layer
            DO i=1,2
               f(i)=0.
               IF (wl(i) <= 0.) CYCLE
               v=c2/(wl(i)*t_level(ii,k))
               DO m=1,1000
                  w=FLOAT(m)
                  wv=w*v
                  ff=(w**(-4))*(EXP(-wv))*(((wv+3.)*wv+6.)*wv+6.)
                  f(i)=f(i)+ff
                  IF (f(i) == 0.) CYCLE
                  IF (ff/f(i) <= facc) EXIT
               ENDDO
            ENDDO
            fst4(ii,k)=c*(f(2)-f(1))*5.67D-08*t_level(ii,k)*           &
                       t_level(ii,k)*t_level(ii,k)*t_level(ii,k)
!            fst4(ii,k)=planck_function(REAL(t_level(ii,k)))             ! Gets fluxes in 15um band only, a la LMD model.
         ENDDO
         DO k=1,n_layer
            diff_planck(ii,k)=fst4(ii,k)-fst4(ii,k-1)
            DO i=1,2
               f(i)=0.
               IF (wl(i) <= 0.) CYCLE
               v=c2/(wl(i)*temp(ii,k))
               DO m=1,1000
                  w=FLOAT(m)
                  wv=w*v
                  ff=(w**(-4))*(EXP(-wv))*(((wv+3.)*wv+6.)*wv+6.)
                  f(i)=f(i)+ff
                  IF (f(i) == 0.) CYCLE
                  IF (ff/f(i) <= facc) EXIT
               ENDDO
            ENDDO
            diff_planck_2(ii,k)=c*(f(2)-f(1))*5.67D-08*temp(ii,k)*     & ! Use diff_planck_2 as temporary storage of planck function
                                temp(ii,k)*temp(ii,k)*temp(ii,k)         !    at layer temperatures...
            diff_planck_2(ii,k)=2.*(fst4(ii,k)+fst4(ii,k-1)-2.*        & ! ...and solve for it in this line.
                                diff_planck_2(ii,k))
            diff_planck_2(ii,k)=0.                                       ! This line turns off diff_planck_2
         ENDDO
      ENDDO

IF (kval_debug) PRINT*, '---> Exited SUBROUTINE cfst4'

      END SUBROUTINE cfst4

!-----------------------------------------------------------------------
!+ Subroutine to calculate the fluxes within the band with one gas.    !
!                                                                      !
! Method:                                                              !
!       Here we obtain the appropriate k-coefficient for each of the   !
!       individual terms and pass it along for "monochromatic"         !
!       calculations.  The results are then summed.                    !
!-----------------------------------------------------------------------
      SUBROUTINE solve_band_one_gas(n_profile, n_layer, p, t, d_mass,  &
                                    i_2stream, i_band, gas_mix_ratio,  &
#if ( WRF_TITAN == 1 )
                                    include_tholin,                    &
#endif
#if ( WRF_MARS == 1 )
                                    include_dust,                      &
                                    include_water_ice,                 &
                                    include_co2_ice,                   &
#endif
                                    aerosol_opacity, isolir,           &
                                    solar_flux,                        &
                                    planck_source_top,                 &
                                    planck_source_bottom,              &
                                    diff_planck_band,                  &
                                    diff_planck_band_2,                &
                                    albedo_surface,                    &
                                    thermal_ground_band, k_gray_ext,   &
                                    k_gray_scat, asymmetry_gas,        &
                                    forward_scatter,                   &
#if ( WRF_TITAN == 1 )
                                    omega_tholin, asymmetry_tholin,    &
#endif
#if ( WRF_MARS == 1 )
                                    omega_dust, asymmetry_dust,        &
                                    omega_ice, asymmetry_ice,          &
                                    omega_co2ice, asymmetry_co2ice,    &
#endif
                                    flux_direct_band,                  &
                                    flux_total_band,                   &
                                    flux_unscattered_band,             &
                                    weights,                           &
                                    trans_layer_band, j,               &
                                    nan_detail_check,                  &
                                    sec_0, p_top_check,                &
                                    ra_kdm_taulayer_calculation,       &
                                    ra_kdm_clearsky_calculation,       &
                                    kdm_mode_flag                      )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                            n_profile, &
                                                              n_layer, &
                                                            i_2stream, &
                                                               i_band, &
                                                               isolir, &
                                                                    j

      REAL(KIND(0.d0)), DIMENSION(n_terms), INTENT(IN   ) ::           &
                                                              weights

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer, npd_species),    &
              INTENT(IN   ) ::                                         &
                                                        gas_mix_ratio

      REAL(KIND(0.d0)),                                                &
           DIMENSION(n_profile, 0:n_layer, npd_aerosol_species),       &
           INTENT(IN   ) ::                           aerosol_opacity

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                         p, &
                                                                    t

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                    d_mass, &
                                                     diff_planck_band, &
                                                   diff_planck_band_2                                       

      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                           k_gray_ext, &
                                                          k_gray_scat, &
                                                        asymmetry_gas, &
                                                      forward_scatter   
#if ( WRF_TITAN == 1 )
      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                          omega_tholin, &
                                                      asymmetry_tholin
#endif
#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                           omega_dust, &
                                                       asymmetry_dust, &
                                                            omega_ice, &
                                                        asymmetry_ice, &
                                                         omega_co2ice, &
                                                     asymmetry_co2ice
#endif

      INTEGER, INTENT(IN), OPTIONAL ::                  kdm_mode_flag

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN   ) ::         &
                                                           solar_flux, &
                                                    planck_source_top, &
                                                 planck_source_bottom, &
                                                       albedo_surface, & 
                                                  thermal_ground_band            

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN), OPTIONAL ::  &
                                                                sec_0

      REAL(KIND(0.d0)), INTENT(IN   ) ::                  p_top_check                                                                

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer,2),             &
           INTENT(  OUT) ::                          flux_direct_band, &
                                                flux_unscattered_band

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2,2),           &
           INTENT(  OUT) ::                           flux_total_band                                      

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(  OUT) ::                          trans_layer_band

      LOGICAL, INTENT(IN) ::              ra_kdm_taulayer_calculation, &
                                          ra_kdm_clearsky_calculation, &
                                                     nan_detail_check

#if ( WRF_TITAN == 1 )
      LOGICAL, INTENT(IN) ::                           include_tholin
#endif
#if ( WRF_MARS == 1 )
      LOGICAL, INTENT(IN) ::                             include_dust, &
                                                    include_water_ice, &
                                                      include_co2_ice
#endif

!c-mm  Local variables:

      INTEGER                                                          &
                                                                  iex, &
                                                              i, k, l
      REAL(KIND(0.d0)), DIMENSION(n_terms, n_profile, n_layer) ::       &
                                                            k_gas_abs

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                  tau_layer_band_term

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer, 2) ::          &
                                                     flux_direct_part, &
                                                flux_unscattered_part

      REAL(KIND(0.d0)), DIMENSION(n_profile) ::                        &
                                                        source_ground, &
                                                      flux_inc_direct, &
                                                        flux_inc_down                 

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2, 2) ::        &
                                                      flux_total_part

      REAL(KIND(0.d0)) ::                                       start, &
                                                               finish, &
                                                        sat_mix_ratio

IF (kval_debug) PRINT*, '---> Entered SUBROUTINE solve_band_one_gas'

      trans_layer_band=0.
      !these need to be zerod because they're always added onto and don't get zerod always.
      flux_total_band=0.
      flux_direct_band=0.
      flux_unscattered_band=0.
      !moving the get_kval here and calculating all gauss terms at once is worth 2 dynamics timesteps.
        DO l=1,n_profile
            DO k=1,n_layer
               k_gas_abs(:,l,k)=get_kvals(p(l,k),t(l,k),                 & ! Assign monochromatic absorption coeff.
                       gas_mix_ratio(l,k,1),i_band,isolir)

!             IF ( (isolir == ip_solar) .and. (i_band /= 4) ) THEN ! mir - temporarily turn off all gas abs except band4
!                k_gas_abs(:,l,k) = 0.
!             endif
               
            ENDDO
         ENDDO

        
      DO iex=1, n_terms                                                  ! Loop over all the terms in the band

         IF (isolir == ip_solar) THEN                                    !                            -----
            source_ground=0.                                             !                              |
            flux_inc_down=solar_flux                                     !                              |
            flux_inc_direct=solar_flux                                   !                              |
         ELSE IF (isolir == ip_infra_red) THEN                           ! Set the appropriate boundary terms for the total upward 
            flux_inc_direct=0.                                           !                     and downward fluxes.
            flux_inc_down=-planck_source_top                             !                              |
            source_ground=thermal_ground_band-(1.-                     & !                              |
                          albedo_surface)*planck_source_bottom           !                              |
         ENDIF                                                           !                            -----

         CALL monochromatic_flux(n_layer, n_profile, d_mass,            &
                                 i_2stream, k_gas_abs(iex,:,:), isolir, &
                                 diff_planck_band, diff_planck_band_2,  &
#if ( WRF_TITAN == 1 )
                                 include_tholin,                        &
#endif
#if ( WRF_MARS == 1 )
                                 include_dust,                          &
                                 include_water_ice,                     &
                                 include_co2_ice,                      &
#endif
                                 flux_inc_direct,                      &
                                 flux_inc_down, albedo_surface,        &
                                 source_ground, k_gray_ext,            &
                                 k_gray_scat, asymmetry_gas,           &
                                 forward_scatter,                      &
#if ( WRF_TITAN == 1 )
                                 omega_tholin, asymmetry_tholin,       &
#endif
#if ( WRF_MARS == 1 )
                                 omega_dust, asymmetry_dust,           &
                                 omega_ice, asymmetry_ice,             &
                                 omega_co2ice, asymmetry_co2ice,       &
#endif
                                 flux_direct_part,                     &
                                 flux_unscattered_part,                &
                                 flux_total_part, gas_mix_ratio,       &
                                 aerosol_opacity, i_band, iex,         &
                                 tau_layer_band_term, j,               &
                                 nan_detail_check,                     &
                                 sec_0=sec_0,                          &
                                 ra_kdm_clearsky_calculation=ra_kdm_clearsky_calculation, &
                                 kdm_mode_flag=kdm_mode_flag )

if (ra_kdm_taulayer_calculation) then
        where(tau_layer_band_term < 300)
         trans_layer_band=trans_layer_band+weights(iex)*               &
                          EXP(-tau_layer_band_term)
        end where
endif
         CALL augment_flux(isolir, weights(iex), flux_direct_band,     & ! Increment the fluxes within the band
                           flux_unscattered_band,                      &
                           flux_total_band, flux_direct_part,          & 
                           flux_unscattered_part,                      &
                           flux_total_part, n_profile, n_layer)          

!         band_extinction = band_extinction + ( weight(iex) * exp(-k_gas_abs))

      ENDDO

IF (kval_debug) PRINT*, '---> Exited SUBROUTINE solve_band_one_gas'

      END SUBROUTINE solve_band_one_gas

!-----------------------------------------------------------------------
!+ Function to interpolate within pressure and mixing ratio to find    !
!  the appropriate k-values from an array.                             !
!                                                                      !
! Method:                                                              !
!        Find where in the press array the given pressure falls, and   !
!        interpolate the true k-value from adjacent values. Likewise   !
!        for the mixing ratio.                                         !
!-----------------------------------------------------------------------
      FUNCTION get_kvals(p, t, gas_mix_ratio, band, isolir)

      IMPLICIT NONE
    
!     input variables:
      REAL(KIND(0.d0)), INTENT(IN   ) ::                               &
                                                                    p, &
                                                                    t, &
                                                        gas_mix_ratio            

      INTEGER, INTENT(IN   ) ::                                        &
                                                                 band, &
                                                               isolir

!     output variables:
      REAL(KIND(0.d0)), DIMENSION(n_terms) ::                          &
                                                            get_kvals

!     local variables
      INTEGER  ::                                              i_term, &
                                                               i_band
      REAL(KIND(0.d0)) ::                                log_p, log_m
      INTEGER,DIMENSION(1)  :: i1d, i2d, i3d
      INTEGER               :: i1, i2, i3
      REAL(KIND=KIND(0.D0)) :: d1, d2, d3


      IF(isolir == ip_solar) THEN
        i_band = n_band_ir + band
      ELSE IF(isolir == ip_infra_red) THEN
        i_band = band
      ELSE
        call wrf_error_fatal("KDM get_kvals, vis ir band setup error")
      ENDIF

      log_p = log(p)
      if(gas_mix_ratio > kval_mix_min) then  ! protect against .le. 0. tracer
       log_m = log(gas_mix_ratio)
      else
       log_m = log(kval_mix_min)
      endif

      ! it is assumed that kval_press and kval_mix have already been converted to log values
      ! it is assumed that kval_press, kval_mix, kval_temps are all monotonic and increasing

      ! find the lower index bracketing the value we want:
      i1d = minloc( (kval_temps-t),     MASK=(kval_temps-t)     .gt. 0.)
      i2d = minloc( (kval_press-log_p), MASK=(kval_press-log_p) .gt. 0.)      
      i3d = minloc( (kval_mix-log_m),   MASK=(kval_mix-log_m)   .gt. 0.)

      ! make sure min value is 1 and max value is dimension-1
      i1 = max ( min ( (i1d(1) - 1 ), (n_temps-1) ), 1)
      i2 = max ( min ( (i2d(1) - 1 ), (n_press-1) ), 1)
      i3 = max ( min ( (i3d(1) - 1 ), (n_mix-1)   ), 1)
    
      ! find the distance from the lower index on scale from 0. -> 1.
      d1 = max (min ( (t - kval_temps(i1))   /(kval_temps(i1+1)-kval_temps(i1)) , 1.d0), 0.d0)
      d2 = max (min ((log_p - kval_press(i2))/(kval_press(i2+1)-kval_press(i2)) , 1.d0), 0.d0)
      d3 = max (min ( (log_m - kval_mix(i3)) /(kval_mix(i3+1)  -kval_mix(i3))   , 1.d0), 0.d0)
      ! value is limited to between 0. and 1. and so effectively uses lowest or highest value for 
      !                                                                any beyond-bounds lookup

      DO i_term = 1,n_terms

        get_kvals(i_term) =                                                       &
         (1.d0-d1)*(1.d0-d2)*(1.d0-d3)*kval_array(i1  ,i2  ,i3  ,i_band,i_term) + &
         (     d1)*(1.d0-d2)*(1.d0-d3)*kval_array(i1+1,i2  ,i3  ,i_band,i_term) + &
         (1.d0-d1)*(     d2)*(1.d0-d3)*kval_array(i1  ,i2+1,i3  ,i_band,i_term) + &
         (     d1)*(     d2)*(1.d0-d3)*kval_array(i1+1,i2+1,i3  ,i_band,i_term) + &
         (1.d0-d1)*(1.d0-d2)*(     d3)*kval_array(i1  ,i2  ,i3+1,i_band,i_term) + &
         (     d1)*(1.d0-d2)*(     d3)*kval_array(i1+1,i2  ,i3+1,i_band,i_term) + &
         (1.d0-d1)*(     d2)*(     d3)*kval_array(i1  ,i2+1,i3+1,i_band,i_term) + &
         (     d1)*(     d2)*(     d3)*kval_array(i1+1,i2+1,i3+1,i_band,i_term)
        
      ENDDO

      END FUNCTION get_kvals

!-----------------------------------------------------------------------
!+ Subroutine to solve for the monochromatic fluxes.                   !
!                                                                      !
! Method:                                                              !
!       The final single scattering properties are calculated and      !
!       rescaled.                                                      !
!-----------------------------------------------------------------------
      SUBROUTINE monochromatic_flux(n_layer, n_profile, d_mass,        &
                                    i_2stream, k_gas_abs, isolir,      &
                                    diff_planck, diff_planck_2,        &
#if ( WRF_TITAN == 1 )
                                    include_tholin,                    &
#endif
#if ( WRF_MARS == 1 )
                                    include_dust,                      &
                                    include_water_ice,                 &
                                    include_co2_ice,                   &
#endif
                                    flux_inc_direct,                   &
                                    flux_inc_down, albedo_surface,     &
                                    source_ground, k_gray_ext,         &
                                    k_gray_scat, asymmetry_gas,        &
                                    forward_scatter,                   &
#if ( WRF_TITAN == 1 )
                                    omega_tholin, asymmetry_tholin,    &
#endif
#if ( WRF_MARS == 1 )
                                    omega_dust, asymmetry_dust,        &
                                    omega_ice, asymmetry_ice,          &
                                    omega_co2ice, asymmetry_co2ice,    &
#endif
                                    flux_direct,                       &
                                    flux_unscattered,                  &
                                    flux_total, gas_mix_ratio,         &
                                    aerosol_opacity, i_band, iex,      &
                                    tau_free, j,                       &
                                    nan_detail_check,                  &
                                    sec_0,                             &
                                    ra_kdm_clearsky_calculation,       &
                                    kdm_mode_flag                      )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                              n_layer, &
                                                            n_profile, &
                                                            i_2stream, &
                                                               isolir, &
                                                               i_band, &
                                                                  iex, &
                                                                    j

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer, npd_species),    &
              INTENT(IN   ) ::                                         &
                                                        gas_mix_ratio

      REAL(KIND(0.d0)),                                                &
           DIMENSION(n_profile, 0:n_layer, npd_aerosol_species),       &
           INTENT(IN   ) ::                           aerosol_opacity

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                    d_mass, &
                                                            k_gas_abs, &
                                                          diff_planck, &
                                                        diff_planck_2                                            

      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                           k_gray_ext, &
                                                          k_gray_scat, &
                                                        asymmetry_gas, &
                                                      forward_scatter
#if ( WRF_TITAN == 1 )
      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                          omega_tholin, &
                                                      asymmetry_tholin
#endif
#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                           omega_dust, &
                                                       asymmetry_dust, &
                                                            omega_ice, &
                                                        asymmetry_ice, &
                                                         omega_co2ice, &
                                                     asymmetry_co2ice
#endif

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN   ) ::         &
                                                      flux_inc_direct, &
                                                        flux_inc_down, &
                                                       albedo_surface, &
                                                        source_ground

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN), OPTIONAL ::  &
                                                                sec_0

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer, 2),            &
           INTENT(  OUT) ::                               flux_direct, &
                                                     flux_unscattered

      REAL(KIND(0.d0)), DIMENSION(n_profile,2*n_layer+2,2),            &
           INTENT(  OUT) ::                                flux_total                                               

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 & ! This is used to determine the layer transmissivity
           INTENT(  OUT) ::                                  tau_free    !    as a diagnostic output.

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer)                  & ! This is used to determine the layer transmissivity
                         ::                              tau_unscaled
   
      LOGICAL, INTENT(IN) ::              ra_kdm_clearsky_calculation, &
                                                     nan_detail_check

      INTEGER, INTENT(IN), OPTIONAL :: kdm_mode_flag

#if ( WRF_TITAN == 1 )
      LOGICAL, INTENT(IN) ::                           include_tholin
#endif
#if ( WRF_MARS == 1 )
      LOGICAL, INTENT(IN) ::                             include_dust, &
                                                    include_water_ice, &
                                                       include_co2_ice
#endif
!c-mm  Local variables:

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               & 
                                                           omega_free, &
                                                       asymmetry_free, &
                                                 forward_scatter_free

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer) ::   flux_dir, & ! Direct fluxes returned from two_stream
                                                             flux_uns    ! Unscattered fluxes

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2) :: flux_tot    ! Total fluxes returned from two_stream

      INTEGER ::                                                  i,k, &
                                                          ip_dust_end, &
                                                          ip_dust_start

      REAL, DIMENSION(n_profile, 0:n_layer) :: aaa

      LOGICAL :: include_aerosol

! If there's no aerosol, we only need to loop once
!   for the clear-sky case.  
! If there is aerosol we loop twice, once "aerosol-free" 
!    and once with the aerosol.  This is necessary to extract
!    the heating due to the aerosol alone.

      !These are needed because we have an optional entry in the array (clear sky)
      !if that isn't calculated these never get filled. fill them with zeros now.
      flux_direct = 0.d0
      flux_unscattered = 0.d0
      flux_total = 0.d0

      flux_dir=0.d0
      flux_uns=0.d0
      flux_tot=0.d0

      include_aerosol = .false.
#if ( WRF_TITAN == 1 )
      include_aerosol = include_tholin
#endif
#if ( WRF_MARS == 1 )
      include_aerosol = (include_dust).OR.(include_water_ice).OR.(include_co2_ice)
#endif
      IF (include_aerosol) THEN
         ip_dust_end = 2
      ELSE
         ip_dust_end = 1
      ENDIF

      if (ra_kdm_clearsky_calculation) then
        ip_dust_start=1
      else 
        ip_dust_start=ip_dust_end
      ENDIF
                                                          
      DO i=ip_dust_start,ip_dust_end
         CALL single_scattering(d_mass, k_gray_ext, k_gray_scat,       &
                                k_gas_abs, asymmetry_gas,              &
#if ( WRF_TITAN == 1 )
                                omega_tholin, asymmetry_tholin,        &
                                include_tholin,                        &
#endif
#if ( WRF_MARS == 1 )
                                omega_dust, asymmetry_dust,            &
                                omega_ice, asymmetry_ice,              &
                                omega_co2ice, asymmetry_co2ice,        &
                                include_dust,                          &
                                include_water_ice,                     &
                                include_co2_ice,                       &
#endif
                                tau_free, omega_free,                  &
                                asymmetry_free, gas_mix_ratio,         &
                                aerosol_opacity, i_band, iex,          & 
                                n_profile, n_layer, isolir, i          )

         tau_unscaled = tau_free
                        
         forward_scatter_free=asymmetry_free*asymmetry_free              ! Delta-rescale asymmetry, tau and omega.
         asymmetry_free=(asymmetry_free-forward_scatter_free)/         & ! *NOTE* This happens after we convolve
                        (1.-forward_scatter_free)                        !    dust and gas values
         tau_free=tau_free*(1.-omega_free*forward_scatter_free)
         omega_free=omega_free*(1.-forward_scatter_free)/              &
                               (1.-omega_free*forward_scatter_free)

         CALL two_stream(n_layer, n_profile, i_2stream, isolir,        & ! Solve the equations using a two-
                         diff_planck, diff_planck_2, flux_inc_down,    & !    stream scheme with no clouds.
                         flux_inc_direct, albedo_surface,              &
                         source_ground, tau_free,                      &
                         tau_unscaled,                                 &
                         omega_free,                                   &
                         asymmetry_free, flux_dir, flux_uns,           &
                         flux_tot, j, nan_detail_check,                &
                         sec_0=sec_0, kdm_mode_flag=kdm_mode_flag)

         flux_direct(:,:,i)=flux_dir
         aaa=real(flux_dir)
         if(nan_detail_check)                                             &
              call test_for_bad_detailed_2d(aaa,"flux_dir in monochro",   &
                                           1,n_profile,0,n_layer,         &
                                           1,n_profile,0,n_layer,         &
                                        .false.,0.,kdm_mode=kdm_mode_flag )
         flux_unscattered(:,:,i)=flux_uns
         flux_total(:,:,i)=flux_tot

      ENDDO

      END SUBROUTINE monochromatic_flux

!-----------------------------------------------------------------------
!+ Subroutine to find the optical depth and single scattering albedo.  !
!                                                                      !
! Method:                                                              !
!       Depending on the treatment of scattering, the optical and      !    
!       single scattering albedo are determined from the extinctions   !
!       supplied.                                                      !
!-----------------------------------------------------------------------
      SUBROUTINE single_scattering(d_mass, k_gray_ext, k_gray_scat,    &
                                   k_gas_abs, asymmetry_gas,           &
#if ( WRF_TITAN == 1 )
                                   omega_tholin, asymmetry_tholin,     &
                                   include_tholin,                     &
#endif
#if ( WRF_MARS == 1 )
                                   omega_dust, asymmetry_dust,         &
                                   omega_ice, asymmetry_ice,           &
                                   omega_co2ice, asymmetry_co2ice,     &
                                   include_dust,                       &
                                   include_water_ice,                  &
                                   include_co2_ice,                    &
#endif
                                   tau,                                &
                                   omega, asymmetry, gas_mix_ratio,    &
                                   aerosol_opacity, i_band, iex,       &
                                   n_profile, n_layer, isolir, i       )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                               i_band, &
                                                                  iex, &
                                                            n_profile, &
                                                              n_layer, &
                                                               isolir, &
                                                                    i

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                    d_mass, &
                                                            k_gas_abs

      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                           k_gray_ext, &
                                                          k_gray_scat, &
                                                        asymmetry_gas

#if ( WRF_TITAN == 1 )
      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                          omega_tholin, &
                                                      asymmetry_tholin
      LOGICAL, INTENT(IN   ) ::                                        &
                                                        include_tholin
#endif
#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), dimension(n_profile, n_layer), INTENT(IN   ) :: &
                                                           omega_dust, &
                                                       asymmetry_dust, &
                                                            omega_ice, &
                                                        asymmetry_ice, &
                                                         omega_co2ice, &
                                                     asymmetry_co2ice

      LOGICAL, INTENT(IN   ) ::                                        &
                                                         include_dust, &
                                                    include_water_ice, &
                                                      include_co2_ice
#endif

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer, npd_species),    &
              INTENT(IN   ) ::                                         &
                                                        gas_mix_ratio

      REAL(KIND(0.d0)),                                                &
           DIMENSION(n_profile, 0:n_layer, npd_aerosol_species),       &
           INTENT(IN   ) ::                           aerosol_opacity

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(  OUT) ::                                       tau, &
                                                                omega, &
                                                            asymmetry
!      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer), intent(in), optional &
!                                        ::               ndust, reff_index

                                            

!c-mm Local variables     

      INTEGER ::                                                    k,ii,kk

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                              tau_gas, &
                                                         tau_gray_gas, & !c-mm  Just scattering (Rayleigh), not absorption.
                                                            omega_gas, &
                                                          tau_aerosol
#if ( WRF_TITAN == 1 )
      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                             tau_tholin
#endif
#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                             tau_dust, &
                                                              tau_ice, &
                                                           tau_co2ice
#endif

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) :: coefficient


!c-mm      tau_gas=(k_gray_ext+(k_gas_abs*gas_mix_ratio(:,:,2)))*d_mass       
!c-mm      omega_gas=k_gray_scat/(k_gray_ext+(k_gas_abs*                    & 
!c-mm                gas_mix_ratio(:,:,2))+tol_machine)                       
      tau_gas=(k_gray_ext+k_gas_abs)*d_mass
      tau_gray_gas=k_gray_scat*d_mass
      omega_gas=k_gray_scat/(k_gray_ext+k_gas_abs                      & 
                +tol_machine)                                            
      

#if ( WRF_TITAN == 1 )
      tau_tholin = aerosol_opacity(:, 1:n_layer, 1)
#endif
     
#if ( WRF_MARS == 1 )
      IF (include_dust) THEN
      !do DUST
         IF (is_two_moment_dust .and. (btest(kdm_options, ko_multimie_dust))) THEN
            !pre-calculated in gray_extinction
            tau_dust = aerosol_opacity(:, 1:n_layer, 1)
!c-mm            write(0,*) 'tau dust ',tau_dust
         ELSE IF(.not.is_two_moment_dust) THEN
            IF (isolir == ip_solar) THEN                                       !c-mm
               coefficient = kCext_dust_solar(i_band)/kCext_dust_ref
            ELSE IF (isolir == ip_infra_red) THEN
               coefficient = kCext_dust_IR(i_band)/kCext_dust_ref
            ENDIF
            DO k=0,n_layer-1                                                
               tau_dust(:,k+1)=(aerosol_opacity(:,k+1,1)-               &
                             aerosol_opacity(:,k,1))*                   &
                             coefficient(:,k+1)
            ENDDO
         ENDIF
      !end dust
      ELSE
        tau_dust = 0.0d0
      ENDIF
      
      IF (include_water_ice) THEN
      !do ICE
         IF (is_two_moment_water_ice .and. (btest(kdm_options, ko_multimie_water))) THEN
            !pre-calculated in gray_extinction
            tau_ice = aerosol_opacity(:, 1:n_layer, 2)
         ELSE IF(.not.is_two_moment_water_ice) THEN
            IF (isolir == ip_solar) THEN                                       !c-mm
               coefficient = kCext_ice_solar(i_band)/kCext_ice_ref
            ELSE IF (isolir == ip_infra_red) THEN
               coefficient = kCext_ice_IR(i_band)/kCext_ice_ref
            ENDIF
            DO k=0,n_layer-1                                                
               tau_ice(:,k+1)=(aerosol_opacity(:,k+1,2)-               &
                              aerosol_opacity(:,k,2))*                  &
                              coefficient(:,k+1)
           ENDDO
         ENDIF 
        !end water ice
      ELSE
        TAU_ICE = 0.d0
      ENDIF

      IF (include_co2_ice) THEN
      !do CO2 ICE
!c-mm         IF (btest(kdm_options, ko_multimie_water)) THEN
!c-mm            !pre-calculated in gray_extinction
!c-mm            tau_co2ice = aerosol_opacity(:, 1:n_layer, 3)
!c-mm         ELSE
            IF (isolir == ip_solar) THEN                                       !c-mm
               coefficient = kCext_co2ice_solar(i_band)/kCext_co2ice_ref
            ELSE IF (isolir == ip_infra_red) THEN
               coefficient = kCext_co2ice_IR(i_band)/kCext_co2ice_ref
            ENDIF
            DO k=0,n_layer-1                                                
               tau_co2ice(:,k+1)=(aerosol_opacity(:,k+1,3)-               &
                              aerosol_opacity(:,k,3))*                  &
                              coefficient(:,k+1)
            ENDDO
!c-mm         ENDIF 
        !end co2 ice
      ELSE
        tau_co2ice = 0.0d0
      ENDIF
#endif

      tau_aerosol = 0.d0
#if ( WRF_MARS == 1 )
      if(include_dust)      tau_aerosol = tau_aerosol + tau_dust
      if(include_water_ice) tau_aerosol = tau_aerosol + tau_ice
      if(include_co2_ice)   tau_aerosol = tau_aerosol + tau_co2ice
#endif
#if ( WRF_TITAN == 1 )
      if(include_tholin)    tau_aerosol = tau_aerosol + tau_tholin
#endif

      IF (i == ip_dust) THEN                                             ! If we have dust and gas in atmosphere
         WHERE (tau_aerosol == 0.0d0)                                    ! Even if number of aerosols > 0, there may be layers without
            tau=tau_gas
            omega=omega_gas                                              !    dust.  In this case, just set omega/asymmetry to
            asymmetry=asymmetry_gas                                      !    gas values to avoid having to divide by zero below
         ELSEWHERE                                                       !    (as tau=0.0 if tau_dust AND tau_gas = 0.0)
            tau = tau_gas 
#if ( WRF_TITAN == 1 )
            tau = tau + tau_tholin
#endif
#if ( WRF_MARS == 1 )
            tau = tau + tau_dust      +    &
                        tau_ice       +    &
                        tau_co2ice
#endif

!c-mm            omega    =(tau_gas /tau)*omega_gas  +                     &
!c-mm                      (tau_dust/tau)*omega_dust +                     &  ! Values of omega and asymmetry are "weighted" by 
!c-mm                      (tau_ice /tau)*omega_ice  +                     &
!c-mm                      (tau_co2ice /tau)*omega_co2ice
            omega = tau_gray_gas                                 ! first add in the gas
#if ( WRF_TITAN == 1 )
            omega = omega + (omega_tholin*tau_tholin)            ! then the aerosols
#endif
#if ( WRF_MARS == 1 )
            omega = omega + (omega_dust*tau_dust)       +  &     ! then the aerosols
                            (omega_ice*tau_ice)         +  &
                            (omega_co2ice*tau_co2ice)
#endif
            omega = omega/tau                                    ! finally normalize by tau

!c-mm            asymmetry=(tau_gas /tau)*asymmetry_gas  +                 &
!c-mm                      (tau_dust/tau)*asymmetry_dust +                 &  !    relative opacities of dust and gas
!c-mm                      (tau_ice /tau)*asymmetry_ice  +                 &
!c-mm                      (tau_co2ice /tau)*asymmetry_co2ice
            asymmetry = tau_gas*omega_gas*asymmetry_gas                                ! first add in the gas
#if ( WRF_TITAN == 1 )
            asymmetry = asymmetry + (tau_tholin*  omega_tholin*  asymmetry_tholin)     ! then the aerosols
#endif
#if ( WRF_MARS == 1 )
            asymmetry = asymmetry + (tau_dust*  omega_dust*  asymmetry_dust) + &       ! then the aerosols
                                    (tau_ice*   omega_ice*   asymmetry_ice)  + &
                                    (tau_co2ice*omega_co2ice*asymmetry_co2ice)
#endif
            asymmetry = asymmetry/(tau*omega)                                          ! finally normalized by tau*omega
         ENDWHERE
                  
      ELSE IF (i == ip_clear) THEN                                       ! If we have gas only
         tau=tau_gas
         omega=omega_gas
         asymmetry=asymmetry_gas
      ENDIF

      END SUBROUTINE single_scattering

!-----------------------------------------------------------------------
!+ Subroutine to solve the two-stream equations in a column.           !
!                                                                      !
! Method:                                                              !
!       The coefficients of the two-stream equations are calculated.   !
!       From these we obtain the transmission and reflection           !
!       coefficients and the source terms. Depending on the solver     !
!       selected, an appropriate set of matrix equations is formulated !
!       and solved to give the fluxes.                                 !
!-----------------------------------------------------------------------
      SUBROUTINE two_stream(n_layer, n_profile, i_2stream, isolir,     &
                            diff_planck, diff_planck_2, flux_inc_down, &
                            flux_inc_direct,                           &
                            albedo_surface, source_ground, tau,        &
                            tau_unscaled, omega,                       &
                            asymmetry, flux_direct, flux_unscattered,  &
                            flux_total, j,                             &
                            nan_detail_check,                          &
                            sec_0, kdm_mode_flag                       )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                              n_layer, &
                                                            n_profile, &
                                                               isolir, &
                                                            i_2stream, &
                                                                    j

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                       tau, &
                                                            asymmetry, &
                                                          diff_planck, &
                                                        diff_planck_2, &
                                                         tau_unscaled

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(INOUT) ::                                     omega

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN   ) ::         &
                                                       albedo_surface, &
                                                        source_ground, &
                                                        flux_inc_down, &
                                                      flux_inc_direct

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN), OPTIONAL ::  &
                                                                sec_0

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2),             &
           INTENT(  OUT) ::                                flux_total                                                 

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer),               &
           INTENT(  OUT) ::                               flux_direct, &
                                                     flux_unscattered
      INTEGER, INTENT(IN), OPTIONAL :: kdm_mode_flag
      LOGICAL, INTENT(IN)           :: nan_detail_check
!c-mm  Local variables:

      INTEGER ::                                                       &
                                                                  i,jj

      REAL(KIND(0.d0)),                                                &
           DIMENSION(n_profile, n_layer, npd_source_coeff) ::          &
                                                         source_coeff

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                                trans, &
                                                              reflect, &
                                                              trans_0, &
                                                               s_down, &
                                                                 s_up, &
                                                       trans_unscaled

      REAL, DIMENSION(n_profile, 0:n_layer)  :: aaa
      REAL :: bbb

                                                               
      CALL two_coeff(n_layer, n_profile, i_2stream, asymmetry, omega,  & ! Calculate the two-stream coefficients.
                     tau, tau_unscaled, isolir, trans, reflect,        &
                     trans_0, trans_unscaled,                          &
                     source_coeff, j, nan_detail_check, sec_0=sec_0    ) 

      aaa = real(trans_0)
      if(nan_detail_check)                                         &
        call test_for_bad_detailed_2d(aaa,"trans_0 in two_stream", &
                                 1,n_profile,0,n_layer,            &
                                 1,n_profile,0,n_layer,            &
                                 .false.,0.,kdm_mode=kdm_mode_flag )

      IF (isolir == ip_solar) THEN
           do jj=1,n_profile
            bbb=flux_inc_direct(jj)
            write(wrf_err_message,*) "flux_inc_direct j=",jj
            if(nan_detail_check) &
              call test_for_bad(bbb,trim(wrf_err_message), &
                        .false.,0.)
           enddo

         flux_direct(:,0)=flux_inc_direct
         flux_unscattered(:,0)=flux_inc_direct     
         DO i=1, n_layer                                                 !                    -----
            flux_direct(:,i)=flux_direct(:,i-1)*trans_0(:,i)             !                      |
            flux_unscattered(:,i)=flux_unscattered(:,i-1)*trans_unscaled(:,i)   !               |
            s_up(:,i)=source_coeff(:,i,ip_scf_solar_up)*               & ! Calculate the appropriate source terms.
                      flux_direct(:,i-1)                                 !                      |    
            s_down(:,i)=source_coeff(:,i,ip_scf_solar_down)*           & !                      |
                        flux_direct(:,i-1)                               !                      |
         ENDDO                                                           !                    -----
      ELSE IF (isolir == ip_infra_red) THEN
        flux_direct = 0
        flux_unscattered = 0
         s_up=source_coeff(:,:,ip_scf_ir_1d)*                          & ! Multiply the source coefficients by the 
              diff_planck+source_coeff(:,:,ip_scf_ir_2d)*diff_planck_2   !    Planckian differences to the order 
         s_down=-source_coeff(:,:,ip_scf_ir_1d)*diff_planck+           & !    required.
                 source_coeff(:,:,ip_scf_ir_2d)*diff_planck_2
      !write(0,*) "2s s_down:",s_down(:,46),diff_planck(:,46),diff_planck_2(:,46), &
      !        source_coeff(:,46,ip_scf_ir_1d),source_coeff(:,46,ip_scf_ir_2d), &
      !        ip_scf_ir_1d,ip_scf_ir_2d

      ENDIF

      aaa = real(flux_direct)
      if(nan_detail_check) &
        call test_for_bad_detailed_2d(aaa,"flux_dir in two_stream",        &
                                 1,n_profile,0,n_layer,         &
                                 1,n_profile,0,n_layer,         &
                                 .false.,0.,kdm_mode=kdm_mode_flag )


      CALL solver_homogen_direct(n_layer, n_profile, trans, reflect,   &
                                 s_down, s_up, albedo_surface,         &
                                 flux_direct(:,n_layer),               &
                                 flux_inc_down, source_ground,         &
                                 flux_total, j, isolir)

      aaa = real(flux_direct)
      if(nan_detail_check) &
              call test_for_bad_detailed_2d(aaa,"flux_dir in two_stream after homo",        &
                                 1,n_profile,0,n_layer,         &
                                 1,n_profile,0,n_layer,         &
                                 .false.,0.,kdm_mode=kdm_mode_flag )

      !write(0,*) "2s k=46, ",flux_total(:,(2*46+2))

      END SUBROUTINE two_stream

!-----------------------------------------------------------------------
!+ Subroutine to calculate coefficients in the two-stream equations.   !
!                                                                      !
! Method:                                                              !
!       The basic two-stream coefficients in the differential          !
!       equations are calculated. These are then used to determine the !
!       transmission and reflection coefficients. Coefficients for     !
!       determining the solar or IR source terms are calculated.       !
!-----------------------------------------------------------------------
      SUBROUTINE two_coeff(n_layer, n_profile, i_2stream, asymmetry,   &
                           omega, tau, tau_unscaled,                   &
                           isolir, trans, reflect,                     &
                           trans_0, trans_unscaled, source_coeff, j,   &
                           nan_detail_check, sec_0                     )

      USE module_state_description, ONLY: eddington, discrete_ord,     &
                                     elsasser, hemi_mean, pifm_sw,     &
                                     pifm_lw

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                              n_layer, &
                                                            n_profile, &
                                                               isolir, &
                                                            i_2stream, &
                                                                    j

      LOGICAL, INTENT(IN   ) ::                      nan_detail_check 

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN), OPTIONAL ::  &
                                                                sec_0

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                 asymmetry, &
                                                                  tau, &
                                                         tau_unscaled

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(INOUT) ::                                     omega

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(  OUT) ::                                   reflect, &
                                                              trans_0

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(  OUT) :: trans, trans_unscaled

      REAL(KIND(0.d0)),                                                &
           DIMENSION(n_profile, n_layer, npd_source_coeff),            &
              INTENT(  OUT) ::                                         &
                                                         source_coeff

!c-mm  Local variables:

      INTEGER ::                                                       &
                                                                  i,k           

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                               lambda, &
                                                                 sums, & ! Call it sums so it doesn't interfere with intrinsic SUM
                                                                 diff, &
                                                             gamma_up, &
                                                           gamma_down

      REAL(KIND(0.d0)) ::                                      target              

      REAL(KIND(0.d0)), PARAMETER ::                                   &
                             elsasser_factor =               1.66D+00

      target=1.-tol_div                                                  ! Perturb the single scattering albedo away from 1 to avoid later 
      DO i=1, n_layer                                                    !    division by 0.
         WHERE (omega(:,i) > target) omega(:,i)=target
      ENDDO
      IF (i_2stream == EDDINGTON) THEN                                   !                    -----
         sums=1.5*(1.-omega*asymmetry)                                   !                      |
         diff=2.*(1.-omega)                                              !                      |
      ELSE IF (i_2stream == DISCRETE_ORD) THEN                           !                      |
         sums=SQRT(3.)*(1.-omega*asymmetry)                              !                      |
         diff=SQRT(3.)*(1.-omega)                                        !                      |
      ELSE IF (i_2stream == ELSASSER) THEN                               ! Calculate the basic two-stream coefficients.
         sums=elsasser_factor-1.5*omega*asymmetry                        !                      |
         diff=elsasser_factor*(1.-omega)                                 !                      |
      ELSE IF (i_2stream == HEMI_MEAN) THEN                              !                      |
         sums=2.*(1.-omega*asymmetry)                                    !                      |
         diff=2.*(1.-omega)                                              !                      |
      ELSE IF (i_2stream == PIFM_SW) THEN                                !                      |
         sums=2.-1.5*omega*asymmetry                                     !                      |
         diff=2.*(1.-omega)                                              !                      |
      ELSE IF (i_2stream == PIFM_LW) THEN                                !                      |
         sums=2.-1.5*omega*asymmetry                                     !                      |
         diff=2.*(1.-omega)                                              !                      |
      ENDIF                                                              !                    -----
      lambda=SQRT(sums*diff)                                             ! Lambda is now calculated.
      IF (isolir == ip_solar) THEN

         CALL solar_coefficient_basic(n_layer, n_profile, omega,       & ! Calculate the basic coefficients for the solar source 
                                      asymmetry, i_2stream,            & !    terms.  Lambda may be perturbed by this routine to 
                                      sums, diff, lambda, gamma_up,    & !    avoid ill-conditioning for a singular zenith angle.
                                      gamma_down, sec_0=sec_0)

      ENDIF
      CALL trans_source_coeff(n_layer, n_profile, isolir, tau, sums,   & ! Determine the transmission and reflection coefficients
                              diff, lambda, gamma_up,                  &
                              gamma_down, trans, reflect, trans_0,     &
                              source_coeff, tau_unscaled,              &
                              trans_unscaled,                          &
                              nan_detail_check,                        &
                              j, sec_0=sec_0)

      END SUBROUTINE two_coeff

!-----------------------------------------------------------------------
!+ Subroutine to calculate the basic coefficients for the solar beam.  !
!                                                                      !
! Method:                                                              !
!       Straightforward.                                               !
!-----------------------------------------------------------------------
      SUBROUTINE solar_coefficient_basic(n_layer, n_profile, omega,    &
                                         asymmetry, i_2stream,         &
                                         sums, diff, lambda, gamma_up, &
                                         gamma_down, sec_0)

      USE module_state_description, ONLY: discrete_ord

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                              n_layer, &
                                                            n_profile, &
                                                            i_2stream         

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                     omega, &
                                                            asymmetry                                            

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN), OPTIONAL ::  &
                                                                sec_0

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(INOUT) ::                                      sums, &
                                                                 diff, &
                                                               lambda

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(  OUT) ::                                  gamma_up, &
                                                           gamma_down

!c-mm  Local variables:

      INTEGER ::                                                    i

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::        ksi_0

      REAL(KIND(0.d0)), DIMENSION(n_profile) ::                        &
                                                               factor, &
                                                           test_array

      DO i=1, n_layer
         IF (PRESENT(sec_0)) THEN
            test_array=abs(lambda(:,i)-sec_0)
         ELSE
            test_array=abs(lambda(:,i))                                  ! This subroutine should only be called if sec_0
         END IF                                                          ! is present, but still need a safety statement
         WHERE (test_array < tol_div)                                    ! If lambda is too close to sec_0, perturb it.
            sums(:,i)=(1.+tol_div)*sums(:,i)
            diff(:,i)=(1.+tol_div)*diff(:,i)
            lambda(:,i)=(1.+tol_div)*lambda(:,i)
         ENDWHERE
      ENDDO

      IF (i_2stream == DISCRETE_ORD) THEN
         DO i=1, n_layer
            ksi_0(:,i)=0.
            IF (PRESENT(sec_0)) THEN
               WHERE (sec_0 /= 1.0D+20) ksi_0(:,i)=SQRT(3.)*           &
                     asymmetry(:,i)/sec_0
            END IF
         ENDDO
      ELSE
         DO i=1, n_layer
            ksi_0(:,i)=0.
            IF (PRESENT(sec_0)) THEN
               WHERE (sec_0 /= 1.0D+20) ksi_0(:,i)=1.5*                &
                     asymmetry(:,i)/sec_0
            END IF
         ENDDO
      ENDIF

      DO i=1, n_layer
         factor=0.
         IF (PRESENT(sec_0)) THEN
            WHERE (sec_0 /= 1.0D+20) factor=0.5*omega(:,i)*sec_0/      & ! Determine the basic solar coefficients for the 
               ((lambda(:,i)-sec_0)*(lambda(:,i)+sec_0))                 !    two-stream equations.
            gamma_up(:,i)=factor*(sums(:,i)-sec_0-ksi_0(:,i)*          &
                          (diff(:,i)-sec_0))
            gamma_down(:,i)=factor*(sums(:,i)+sec_0+ksi_0(:,i)*        &
                            (diff(:,i)+sec_0))
         ELSE
            gamma_up(:,i)   = 0.                                         ! This subroutine should only be called if sec_0
            gamma_down(:,i) = 0.                                         ! is present, but still need a safety statement
         END IF
      ENDDO

      END SUBROUTINE solar_coefficient_basic

!-----------------------------------------------------------------------
!+ Subroutine to calculate transmission and reflection coefficients.   !
!                                                                      !
! Method:                                                              !
!        Refer to Hadley Centre notes of radiation scheme to           !
!        understand what's happening here.                             !
!-----------------------------------------------------------------------
      SUBROUTINE trans_source_coeff(n_layer, n_profile, isolir, tau,   &
                                    sums, diff, lambda,                &
                                    gamma_up, gamma_down, trans,       &
                                    reflect, trans_0, source_coeff,    &
                                    tau_unscaled, trans_unscaled,      &
                                    nan_detail_check,                  &
                                    j, sec_0                           )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                              n_layer, &
                                                            n_profile, &
                                                               isolir, &
                                                                    j

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                       tau, &
                                                                 sums, &
                                                                 diff, &
                                                               lambda, &
                                                             gamma_up, &
                                                           gamma_down, &
                                                         tau_unscaled

      LOGICAL, INTENT(IN) ::                         nan_detail_check

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN), OPTIONAL ::  &
                                                                sec_0

      REAL(KIND(0.d0)),                                                &
           DIMENSION(n_profile, n_layer, npd_source_coeff),            &
           INTENT(  OUT) ::                              source_coeff

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(  OUT) ::                                     trans, &
                                                              reflect, &
                                                              trans_0, &
                                                       trans_unscaled
!c-mm  Local variables:

      INTEGER ::                                                 i,jj

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                                gamma, &
                                                          exponential

      real :: aaa

      exponential = 0.0
      where(lambda*tau < 300)
        exponential=EXP(-lambda*tau)
      end where
      gamma=(sums-lambda)/(sums+lambda)
      trans=(exponential*(1.-gamma**2)/(1.-(exponential*               & ! Determine the diffuse transmission and reflection 
            gamma)**2))                                                  !    coefficients.
      reflect=gamma*(1.-exponential**2)/(1.-(exponential*              &
              gamma)**2) 
      !write(0,*) "in trans source ",trans(:,46),reflect(:,46),tau(:,46),lambda(:,46)

      IF (isolir == ip_solar) THEN
         DO i=1,n_layer
            trans_0(:,i)=0.
            IF (PRESENT(sec_0)) THEN
               WHERE (sec_0 /= 1.0D+20)
                  trans_0(:,i)=EXP(-tau(:,i)*sec_0)
                  trans_unscaled(:,i) = EXP(-tau_unscaled(:,i)*sec_0)
               END WHERE
               if(nan_detail_check) then
                  do jj=1,n_profile
                    aaa = real(tau(jj,i))
                    write(wrf_err_message,*) "tau in trans src coef, j=",jj
                    call test_for_bad(aaa,wrf_err_message,.false.,0.)
                    aaa = real(trans_0(jj,i))
                    write(wrf_err_message,*) "trans_0 in trans src coef, j=",jj
                    call test_for_bad(aaa,wrf_err_message,.false.,0.)
                    aaa = real(sec_0(jj))
                    write(wrf_err_message,*) "sec_0 in trans src coef, j=",jj
                    call test_for_bad(aaa,wrf_err_message,.false.,0.)
                  enddo
               endif
            END IF
            source_coeff(:,i, ip_scf_solar_up)=(gamma_up(:,i)-         & ! Calculate the direct transmission and the 
                                               reflect(:,i)*(1.+       & !    source coefficients for the solar beam: 
                                               gamma_down(:,i)))-      & !    in the solar case these are the 
                                               gamma_up(:,i)*          & !    coefficients which will multiply the 
                                               trans(:,i)*trans_0(:,i)   !    direct flux at the top of the layer to
            source_coeff(:,i, ip_scf_solar_down)=trans_0(:,i)*         & !    give the source terms for the upward    
                                                 (1.+                  & !    diffuse flux and the total downward flux.
                                                 gamma_down(:,i)-      &
                                                 gamma_up(:,i)*        & 
                                                 reflect(:,i))-        &
                                                 (1.+                  &
                                                 gamma_down(:,i))*     &
                                                 trans(:,i)
         END DO
      ELSE IF (isolir == ip_infra_red) THEN
         source_coeff(:,:,ip_scf_ir_1d)=((1.-trans+reflect)+           & ! In the case of infra-red radiation, the first 
                                        (tol_machine/                  & !    source coefficient holds the multiplier for
                                        (sqrt_tol_machine+tau)))/      & !    the first difference of the planck function
                                        ((tol_machine/                 & !    across the layer, and the second that for
                                        (sqrt_tol_machine+tau))+tau*   & !    the second difference.
                                        sums)                            ! A tolerance is added to the numerator and the 
                                                                         !    denomiator to avoid ill-conditioning at 
                                                                         !    small optical depths.
        !write(0,*) "trans source ",source_coeff(:,46,ip_scf_ir_1d),trans(:,46),reflect(:,46),tau(:,46),sums(:,46)
         source_coeff(:,:,ip_scf_ir_2d)=-2.+diff*tau
         WHERE (tau > quad_correct_limit)                              & ! Quadratic correction to source 
               source_coeff(:,:,ip_scf_ir_2d)=-2.*(1.-trans-           & !    function. This correction is very 
                                              reflect+                 & !    ill-conditioned for small optical 
                                              sqrt_tol_machine)/(diff* & !    depths so the asymptotic form is 
                                              tau+sqrt_tol_machine)      !    then used.
         source_coeff(:,:,ip_scf_ir_2d)=-(1.+reflect+trans+            & 
                                      source_coeff(:,:,ip_scf_ir_2d))/ &
                                      (sums*tau+sqrt_tol_machine)
      END IF

      END SUBROUTINE trans_source_coeff

!-----------------------------------------------------------------------
!+ Subroutine to calculate fluxes in a homogeneous column directly.    !
!                                                                      !
! Method:                                                              !
!       Straightforward.                                               !
!-----------------------------------------------------------------------
      SUBROUTINE solver_homogen_direct(n_layer, n_profile, trans,      &
                                       reflect, s_down, s_up,          &
                                       albedo_surface,                 &
                                       flux_direct_ground,             &
                                       flux_inc_down, source_ground,   &
                                       flux_total, j, isolir)

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                              n_layer, &
                                                            n_profile, &
                                                                    j, &
                                                               isolir

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer),                 &
           INTENT(IN   ) ::                                     trans, &
                                                              reflect, &
                                                               s_down, &
                                                                 s_up                                                  

      REAL(KIND(0.d0)), DIMENSION(n_profile), INTENT(IN   ) ::         &
                                                       albedo_surface, &
                                                        source_ground, &
                                                        flux_inc_down, &
                                                   flux_direct_ground

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2),             &
           INTENT(  OUT) ::                                flux_total

!c-mm  Local variables:

      INTEGER ::                                                    i

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer) ::               &
                                                                 beta, &
                                                                gamma, &
                                                                    h

      REAL(KIND(0.d0)), DIMENSION(n_profile, n_layer+1) ::             &
                                                                alpha, &
                                                           s_up_prime

      alpha(:,n_layer+1)=albedo_surface                                  ! Initialization at the 
      s_up_prime(:,n_layer+1)=source_ground                              !    bottom for upward elimination:

      DO i=n_layer, 1, -1                                                ! Eliminating loop:
         beta(:,i)=1./(1.-alpha(:,i+1)*reflect(:,i))
         gamma(:,i)=alpha(:,i+1)*trans(:,i)
         h(:,i)=s_up_prime(:,i+1)+alpha(:,i+1)*s_down(:,i)
         alpha(:,i)=reflect(:,i)+beta(:,i)*gamma(:,i)*trans(:,i)
         s_up_prime(:,i)=s_up(:,i)+beta(:,i)*trans(:,i)*h(:,i)
      ENDDO

      flux_total(:,2)=flux_inc_down                                      ! Initialize for backward 
      flux_total(:,1)=alpha(:,1)*flux_total(:,2)+s_up_prime(:,1)         !    substitution.

      DO i=1, n_layer                                                    ! Backward substitution:
         flux_total(:,2*i+1)=beta(:,i)*(h(:,i)+gamma(:,i)*             & ! Upward flux
                             flux_total(:,2*i))                          
         flux_total(:,2*i+2)=s_down(:,i)+trans(:,i)*flux_total(:,2*i)+ & ! Downward flux
                             reflect(:,i)*flux_total(:,2*i+1)            
      ENDDO

      !i=46
      !write(0,*) "shd k=46, ",flux_total(:,2*i+2),s_down(:,i),trans(:,i),flux_total(:,2*i), &
!                             reflect(:,i),flux_total(:,2*i+1)

      END SUBROUTINE solver_homogen_direct

!-----------------------------------------------------------------------
!+ Subroutine to increment a sum of fluxes.                            !
!                                                                      !
! Method:                                                              !
!       The arrays holding the summed fluxes are incremented by a      !
!       weighted sum of the variables suffixed with _INCR. Arguments   !
!       specify which arrays are to be incremented.                    !
!-----------------------------------------------------------------------
      SUBROUTINE augment_flux(isolir, weight_incr, flux_direct,        &
                              flux_unscattered,                        &
                              flux_total, flux_direct_incr,            &
                              flux_unscattered_incr,                   &
                              flux_total_incr, n_profile, n_layer)

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                               isolir, &
                                                            n_profile, &
                                                              n_layer

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer, 2),            &
           INTENT(IN   ) ::                          flux_direct_incr, &
                                                flux_unscattered_incr

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2, 2),          &
           INTENT(IN   ) ::                           flux_total_incr                                        

      REAL(KIND(0.d0)), INTENT(IN   ) ::                  weight_incr

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer, 2),            &
           INTENT(INOUT) ::                               flux_direct, &
                                                     flux_unscattered

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2, 2),          &
           INTENT(INOUT) ::                                flux_total                                             


      !write(0,*) "b k=45, ",flux_total(:,(2*45+2),1),weight_incr,flux_total_incr(:,(2*45+2),1)
      !write(0,*) "b k=46, ",flux_total(:,(2*46+2),1),weight_incr,flux_total_incr(:,(2*46+2),1)
      !write(0,*) "b k=47, ",flux_total(:,(2*47+2),1),weight_incr,flux_total_incr(:,(2*47+2),1)
      IF (isolir == ip_solar) THEN
         flux_direct=flux_direct+weight_incr*flux_direct_incr            ! Increment the actual fluxes.
         flux_unscattered=flux_unscattered+weight_incr*flux_unscattered_incr            ! Increment the unscattered fluxes.
      ENDIF
      flux_total=flux_total+weight_incr*flux_total_incr
      !write(0,*) "a k=45, ",flux_total(:,(2*45+2),1),weight_incr,flux_total_incr(:,(2*45+2),1)
      !write(0,*) "a k=46, ",flux_total(:,(2*46+2),1),weight_incr,flux_total_incr(:,(2*46+2),1)
      !write(0,*) "a k=47, ",flux_total(:,(2*47+2),1),weight_incr,flux_total_incr(:,(2*47+2),1)

      END SUBROUTINE augment_flux

!-----------------------------------------------------------------------
!+ Subroutine to increment the total flux within a spectral band.      !
!                                                                      !
! Method:                                                              !
!       The total flux is incremented by a multiple of the flux within !
!       a spectral band. This routine is similar to AUGMENT_FLUX, but  !
!       here the Planckian flux must be incremented in the IR.         !
!-----------------------------------------------------------------------
      SUBROUTINE augment_total_flux(isolir, planck_source_band,        &
                                    flux_direct, flux_total,           &
                                    flux_unscattered,                  &
                                    flux_direct_band, flux_total_band, &
                                    flux_unscattered_band,             &
                                    planck_flux, n_profile, n_layer)

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                               isolir, &
                                                            n_profile, &
                                                              n_layer

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer),               &
           INTENT(IN   ) ::                         planck_source_band                                         

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer, 2),            &
           INTENT(IN   ) ::                           flux_direct_band,&
                                                 flux_unscattered_band

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2, 2),          &
           INTENT(IN   ) ::                            flux_total_band                                        

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer, 2),            &
           INTENT(INOUT) ::                                flux_direct,&
                                                      flux_unscattered

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer),               &
           INTENT(INOUT) ::                                planck_flux

      REAL(KIND(0.d0)), DIMENSION(n_profile, 2*n_layer+2, 2),          &
           INTENT(INOUT) ::                                 flux_total

      !write(0,*) "b k=46, ",flux_total(:,(2*46+2),1),flux_total_band(:,(2*46+2),1)
      IF (isolir == ip_solar) THEN
         flux_direct=flux_direct+flux_direct_band                        ! Increment the direct fluxes.
         flux_unscattered=flux_unscattered+flux_unscattered_band         ! Increment the unscattered fluxes
      ENDIF
      flux_total=flux_total+flux_total_band

      IF (isolir == ip_infra_red) planck_flux=planck_flux+             & ! Sum the Planckian fluxes for later addition to the 
                                              planck_source_band         !    differential fluxes. Not necessary for a net-flux 
                                                                         !    scheme
      !write(0,*) "a k=46, ",flux_total(:,(2*46+2),1),flux_total_band(:,(2*46+2),1)

      END SUBROUTINE augment_total_flux

!-----------------------------------------------------------------------
!+ Subroutine to assign fluxes to the final arrays.                    !
!                                                                      !
! Method:                                                              !
!       The array FLUX_TOTAL holds the calculated total fluxes         !
!       (differential fluxes in the IR). Upward and downward fluxes    !
!       are passed to FLUX_UP and FLUX_DOWN and incremented by the     !
!       Planckian flux in the IR.                                      !
!-----------------------------------------------------------------------
      SUBROUTINE assign_flux(n_layer, n_profile, flux_total, isolir,   &
                             planck_flux, flux_down, flux_up,          &
#if ( WRF_MARS == 1 )
                             include_dust,                             &
                             include_water_ice,                        &
                             include_co2_ice,                          &
#endif
                             j)

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                        &
                                                              n_layer, &
                                                            n_profile, &
                                                               isolir, &
                                                                    j

      REAL(KIND(0.d0)), DIMENSION(n_profile,2*n_layer+2,2),            &
           INTENT(IN   ) ::                                flux_total                                          

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer),               &
           INTENT(IN   ) ::                               planck_flux   

      REAL(KIND(0.d0)), DIMENSION(n_profile, 0:n_layer, 2),            &
           INTENT(  OUT) ::                                 flux_down, &
                                                              flux_up

#if ( WRF_MARS == 1 )
      LOGICAL, INTENT(IN   ) ::                          include_dust, &
                                                    include_water_ice, &
                                                      include_co2_ice
#endif
!c-mm  Local variables:

      INTEGER ::                                                 i, k, &
                                                          ip_dust_end

! If there's no aerosol, we only need to loop once
!   for the clear-sky case.  
! If there is aerosol we loop twice, once "aerosol-free" 
!    and once with the aerosol.  This is necessary to extract
!    the heating due to the aerosol alone.

      ip_dust_end = 1
#if ( WRF_MARS == 1 )
      IF ((include_dust).OR.(include_water_ice).OR.(include_co2_ice)) THEN
         ip_dust_end = 2
      ENDIF
#endif

      IF (isolir == ip_solar) THEN
         DO k=0, n_layer
            flux_up(:,k,:)=flux_total(:,2*k+1,:)
            flux_down(:,k,:)=flux_total(:,2*k+2,:)
         ENDDO
      ELSE IF (isolir == ip_infra_red) THEN
         DO k=0, n_layer
!            DO i=ip_clear,ip_dust_end
            DO i=ip_clear,ip_dust
               flux_up(:,k,i)=flux_total(:,2*k+1,i)+planck_flux(:,k)
               flux_down(:,k,i)=flux_total(:,2*k+2,i)+planck_flux(:,k)
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE assign_flux

!c-yl:
!-------------------------------------------------------------------
!---- set_fake_layers: sets pressure and temperature          -----|
!---- of the fake layers above the model top                  -----|
!---- model top is defined as the smallest pressure           -----|
!---- in a horizontal plane minval(pf3d(its:ite,kme,jts:jte)) -----|
!------------------------------------------------------------------|
      SUBROUTINE set_fake_layers(p_top, t_sounding, polar, p_scale)

#ifdef mpas
      USE mpas_atmphys_constants, only : g=>gravity, r_d, p0
#else
      USE module_model_constants
#endif
      USE module_setup_tpz

      IMPLICIT NONE

      ! the following are needed to setup the 'fake' layers above the model ptop
      real, intent(in)                                ::        p_top
      real, intent(in)                                ::      p_scale
      integer, intent(in)                             ::   t_sounding
      logical, intent(in)                             ::        polar

      integer, parameter :: error_val = 1000
#if ( WRF_MARS == 1 )
      integer, parameter                              :: gcm_t_sounding = 2
#elif ( WRF_TITAN == 1 )      
      integer, parameter                              :: gcm_t_sounding = 63
#else
      integer, parameter                              :: gcm_t_sounding = error_val
#endif

      logical,save :: first_pass = .true.

      real :: p_fake_top

      if(gcm_t_sounding == error_val) &
        call wrf_error_fatal("KDM: fake layers not setup for this planet")

      ! get the t_profile for this sim if not a gcm
      if(first_pass) then
      if(.not.polar) then  ! we're not a gcm and can use the les_scm t_sounding
              WRITE( wrf_err_message , * ) 'KDM is creating fake layers based on t_sounding=',&
                      t_sounding,' from les_scm block - make sure this is what you want!'
              CALL wrf_message ( TRIM( wrf_err_message ) )
      else
              WRITE( wrf_err_message , * ) 'KDM is creating fake layers based on t_sounding=2 ',&
                      'as a default as it thinks this is a GCM run (ignoring les_scm block) - make sure this is what you want!'
              CALL wrf_message ( TRIM( wrf_err_message ) )
      endif
      endif


      n_fake_layers_needed_visir = 0
      p_fake_top = REAL(p_top,KIND(0.d0))                 ! set to ptop initially

!      write(0,*) 'kdminit thinks p_top=',p_top

!      write(0,*) (REAL(p_top,KIND(0.d0)) .ge. pressure_AWDC), n_fake_layers_needed_visir == 0
!      write(0,*) (p_fake_top .ge. pressure_AWDC), n_fake_layers_needed_visir .gt. 0
      if((REAL(p_top,KIND(0.d0)) .ge. pressure_AWDC) .or. n_fake_layers_needed_visir == 0)  then  ! awdc = 0.005, in this case the model top is at about 100km or higher
              do while ((p_fake_top .ge. pressure_AWDC) .or. n_fake_layers_needed_visir == 0) ! always do at least one fake layer above the model ptop

                 if (n_fake_layers_needed_visir .eq. max_fake_layers_visir) exit !c-yl: add additional condition in case the surface pressure is too high

                 n_fake_layers_needed_visir = n_fake_layers_needed_visir + 1

                 if(n_fake_layers_needed_visir .gt. 1) then
                   ! put fake layers every 1/2 scale height from ptop to 0.005
                   ! 0.6 is exp(-z/H) where z=5km, H-10km, so it's a drop of 0.5 * scale height
                   p_fake_lower_edge(n_fake_layers_needed_visir) = p_fake_lower_edge(n_fake_layers_needed_visir-1)*fake_scale
                   p_fake_centre(    n_fake_layers_needed_visir) = p_fake_lower_edge(n_fake_layers_needed_visir)  *(1.d0+fake_scale)/2.d0
                 else
                   p_fake_lower_edge(n_fake_layers_needed_visir) = REAL(p_top,KIND(0.d0))
                   p_fake_centre(    n_fake_layers_needed_visir) = REAL(p_top,KIND(0.d0))*(1.d0+fake_scale)/2.d0
                 endif
                 !c-yl: p_fake_centre is always half layer above p_fake_lower_edge
                 p_fake_top = p_fake_centre(n_fake_layers_needed_visir)
                 !c-yl: if p_fake_centre is going out of range of pressure in k-table (minimum 1.e-4), make an average between
                 !c-yl: smallest p_fake_lower_edge and 1.e-4. This hardwired 1.e-4 shall be made to be consistent to press(1)
                 !c-yl: perhaps just define 1.e-4 as the lower limit (parameter) of "press" instead of using it in the calculation of "press"
                 !c-yl: This only avoids the case where fake_scale and p_top are too small
                 if (p_fake_top .lt. 1.e-4) p_fake_centre(n_fake_layers_needed_visir) = (p_fake_lower_edge(n_fake_layers_needed_visir) + 1.e-4)/2.

                 if(.not.polar) then  ! we're not a gcm and can use the les_scm t_sounding
                    t_fake_lower_edge(n_fake_layers_needed_visir) = REAL(P2T(t_sounding,                                            &
                                                                             real(p_fake_lower_edge(n_fake_layers_needed_visir)),   &
                                                                             T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale          &
                                                                             ),KIND(0.d0))
                    t_fake_centre(    n_fake_layers_needed_visir) = REAL(P2T(t_sounding,                                            &
                                                                             real(p_fake_centre(n_fake_layers_needed_visir)),       &
                                                                             T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale          &
                                                                             ),KIND(0.d0))
                 else ! if we are a gcm, then use the t_sounding=2 used to set up GCM's in module_initialize_simple_cyl
                    t_fake_lower_edge(n_fake_layers_needed_visir) = REAL(P2T(gcm_t_sounding,                                        &
                                                                             real(p_fake_lower_edge(n_fake_layers_needed_visir)),   &
                                                                             T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale          &
                                                                              ),KIND(0.d0))
                    t_fake_centre(    n_fake_layers_needed_visir) = REAL(P2T(gcm_t_sounding,                                        &
                                                                             real(p_fake_centre(n_fake_layers_needed_visir)),       &
                                                                             T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale          &
                                                                            ),KIND(0.d0))
                 endif

              enddo
      endif
      if(first_pass) then
        WRITE( wrf_err_message , * ) 'KDM has determined that this run needs ',n_fake_layers_needed_visir,&
             ' layers above p_top (AKA that number of fake layers)'
        CALL wrf_message ( TRIM( wrf_err_message ) )
        first_pass = .false.
      endif

     END SUBROUTINE set_fake_layers
!-----------------------------------------------------------------------
!+ Subroutine to initialize the heating rate tendencies and read in    !
!  k-coefficients and specific heat values.                            !
!                                                                      !
! Method:                                                              !
!       Straightforward.                                               !
!-----------------------------------------------------------------------
      SUBROUTINE kdminit(    restart, k_coeff_file,                    &
                             l_double,                                 &
                             t_sounding,polar,                         &
                             do_fake_layers,                           &
                             rs_file_name,                             &
                             mp_physics,                               &
                             namelist_n_band_ir,                       &
                             namelist_n_band_solar,                    &
#if ( WRF_MARS == 1 )
                             dust_params_file,                         & 
                             l_simpledust,                             &
                             simpledust_omega, simpledust_asymmetry,   &
                             ice_params_file, co2ice_params_file,      &
                             include_dust, include_water_ice,          &
                             include_water_vap, include_co2_ice,       &
                             mp_physics_du, ra_du_physics,             &
                             do_chicago_ice_layers,                    &
                             chi_cloud_prs_in, chi_ncloud_in,          &
                             chi_ls1_in, chi_ls2_in, chi_lat1_in,      &
                             chi_lat2_in, chi_lon1_in, chi_lon2_in,    &
                             chi_tval1_in, chi_tval2_in,               &
                             chi_dayornight_in, chi_ice_tau_in,        &
#endif
#if ( WRF_TITAN == 1 )
                             include_tholin, mp_physics_tholin,        &
                             include_ch4_vap,                          &
#endif
                             reg_t_ref, reg_t_iso, reg_lapse_rate,     &
                             ids, ide, jds, jde, kds, kde,             &
                             ims, ime, jms, jme, kms, kme,             &
                             its, ite, jts, jte, kts, kte)

      USE module_setup_tpz
      USE module_state_description

      IMPLICIT NONE

#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), PARAMETER ::                     ref_wvl_dust=0.67,    &  ! Reference wavelength for vis dust opacity (in um)
                                                      ref_wvl_ir_dust=9.3,     &  ! Reference wavelength for IR dust opacity (in um)
                                                    ref_wvl_water_ice=12.1212, &  ! Reference wavelength for IR water ice opacity (in um)
                                                      ref_wvl_co2_ice=0.67        ! Reference wavelength for IR water ice opacity (in um)
#endif

      INTEGER,  INTENT(IN   ) ::         ids, ide, jds, jde, kds, kde, &
                                         ims, ime, jms, jme, kms, kme, &
                                         its, ite, jts, jte, kts, kte

      INTEGER, INTENT(IN   ) ::                    namelist_n_band_ir, &
                                                namelist_n_band_solar

      LOGICAL, INTENT(IN   ) ::                               restart, &
                                                             l_double, &
                                                       do_fake_layers

#if ( WRF_MARS == 1 )
      LOGICAL, INTENT(IN   ) ::                          l_simpledust, &
                                                         include_dust, &
                                                    include_water_ice, &
                                                    include_water_vap, &
                                                      include_co2_ice
#endif
#if ( WRF_TITAN == 1 )
      LOGICAL, INTENT(IN   ) ::                        include_tholin, &
                                                      include_ch4_vap
      INTEGER, INTENT(IN   ) ::                     mp_physics_tholin
#endif

      REAL, INTENT(IN) ::     reg_t_ref, reg_t_iso, reg_lapse_rate

#if ( WRF_MARS == 1 )
      real, intent(in), optional ::                  simpledust_omega, &
                                                   simpledust_asymmetry

      integer, intent(in) :: mp_physics_du, ra_du_physics
#endif
      integer, intent(in) :: mp_physics

#if ( WRF_MARS == 1 )
      real, intent(in) ::    chi_cloud_prs_in,                         &
                             chi_ls1_in, chi_ls2_in, chi_lat1_in,      &
                             chi_lat2_in, chi_lon1_in, chi_lon2_in,    &
                             chi_tval1_in, chi_tval2_in,               &
                             chi_ice_tau_in
      integer, intent(in) :: chi_ncloud_in, chi_dayornight_in
#endif

      ! the following are needed to setup the 'fake' layers above the model ptop
      integer, intent(in)                             ::   t_sounding
      logical, intent(in)                             ::        polar
      character(len = 100), intent(in)                :: rs_file_name
      

      CHARACTER (LEN=80), INTENT(IN   ) ::               k_coeff_file

#if ( WRF_MARS == 1 )
      CHARACTER (LEN=80), INTENT(IN   ) ::           dust_params_file, &
                                                      ice_params_file, &
                                                   co2ice_params_file

      LOGICAL, INTENT(IN) ::                    do_chicago_ice_layers
#endif

      ! Local variables

#if ( WRF_MARS == 1 )
      REAL(KIND(0.d0)), DIMENSION(:), ALLOCATABLE ::         wvl_in_d, &
                                                            cext_in_d, &
                                                           cscat_in_d, &
                                                               g_in_d, &
                                                            qext_in_d, &
                                                           omega_in_d

      REAL(KIND(0.d0)), DIMENSION(:), ALLOCATABLE::          wvl_in_i, &
                                                            cext_in_i, &
                                                           cscat_in_i, &
                                                               g_in_i, &
                                                            qext_in_i, &
                                                           omega_in_i

      REAL(KIND(0.d0)), DIMENSION(:), ALLOCATABLE::          wvl_in_c, &
                                                            cext_in_c, &
                                                           cscat_in_c, &
                                                               g_in_c, &
                                                            qext_in_c, &
                                                           omega_in_c

      REAL(KIND(0.d0)), ALLOCATABLE, DIMENSION(:,:,:,:,:) ::           &
                                                       old_kval_array
#endif

      REAL, DIMENSION(:), ALLOCATABLE :: wavenumbers

      INTEGER                                                 i, j, k, &
                                                        itf, jtf, ktf, &
                                                        c1,c2,c3,c4,c5

      REAL               :: r_in_t,r_in_p,r_in_term,r_in_band,r_in_mix
      INTEGER            :: n_in_t,n_in_p,n_in_term,n_in_band,n_in_mix,n_in_tot
      CHARACTER(LEN=100) :: target_object, header1, header2

      REAL(KIND(0.d0))                          total_solar_flux, jnk, &
                                                           p_fake_top
      LOGICAL            :: old_kvals_file_style


      IF ( namelist_n_band_solar /= n_band_solar ) THEN
         CALL wrf_error_fatal("max_n_band_solar in frame/module_driver_constants and kdm_n_band_solar in Registry not equal.")
      ENDIF
      IF ( namelist_n_band_ir /= n_band_ir ) THEN
         CALL wrf_error_fatal("max_n_band_ir in frame/module_driver_constants and kdm_n_band_ir in Registry not equal.")
      ENDIF


      T_iso = reg_t_ref
      T_ref = reg_t_iso
      lapse_rate = reg_lapse_rate

#ifdef mpas
      jtf=jte
      ktf=kte
      itf=ite
#else
      jtf=MIN(jte,jde-1)
      ktf=MIN(kte,kde-1)
      itf=MIN(ite,ide-1)
#endif

      IF(do_fake_layers) THEN
         CALL wrf_message(trim("KDM: running code with fake layers above model TOA"))
      ELSE
         CALL wrf_message(trim("KDM: running code WITHOUT fake layers above model TOA"))
      ENDIF

#if ( WRF_MARS == 1 )
      kdm_options = 0
      IF ( (ra_du_physics == active_dust_rad) .and. &
           (mp_physics_du == mars_two_moment_dust)) THEN  
        is_two_moment_dust = .true.
        kdm_options=3
        if (btest(kdm_options, ko_debug)) then
          call debug_ko_options(kdm_options)
        endif
      ELSE  
        is_two_moment_dust = .false.
      ENDIF

      IF(mp_physics==mars_two_moment_ice) THEN
          is_two_moment_water_ice = .true.
          kdm_options = kdm_options+4
      ELSE
          is_two_moment_water_ice = .false.
      ENDIF
#endif

! get the t_profile for this sim if not a gcm
      if(.not.polar) then  ! we're not a gcm and can use the les_scm t_sounding
#if ( WRF_MARS == 1 )
              if(t_sounding .eq. 42) call INIT_RADIO_SCIENCE_PROFILE(rs_file_name)
              if(t_sounding .eq. 43) call INIT_GCM_PROFILE(rs_file_name)
#elif ( WRF_TITAN == 1 )
              if(t_sounding .eq. 63) call INIT_LH_PROFILE()
#endif
      endif

! figure out if the kvals file is old format or new
      OPEN(UNIT=2008, FILE=TRIM(k_coeff_file))
      READ(2008,*) target_object
      CLOSE(2008)
      old_kvals_file_style = .true.
#if ( WRF_MARS == 1 )
      IF(trim(adjustl(target_object)) == "MARS")  old_kvals_file_style = .false.
#elif ( WRF_TITAN == 1 )
      IF(trim(adjustl(target_object)) == "TITAN") old_kvals_file_style = .false.
#endif

      IF ( old_kvals_file_style ) THEN

#if ( WRF_MARS == 1 )
                                     n_temps =                          18 ! Number of temperatures in k-coefficient array
                                     n_press =                          26 ! Number of pressures in k-coefficient array
                                       n_mix =                           8 ! Number of mixing ratios in k-coefficient array
                                     n_terms =                          32 ! Number of gaussian quadrature terms in each wl band
 
                ! desired ordering kval_array ( temps, press, mix, wl_band, term )
                ALLOCATE(old_kval_array(n_temps,n_press,n_terms,n_tot_bands,n_mix))
                ALLOCATE(    kval_array(n_temps,n_press,n_mix,n_tot_bands,n_terms)) ! Full array of k-coefficients
                IF(.NOT.ALLOCATED(weights)) ALLOCATE(weights(n_terms))

                call wrf_message("KDM is assuming you passed it an old format kvals file")
                OPEN(UNIT=2008, FILE=TRIM(k_coeff_file))
                READ(2008,*) old_kval_array
                CLOSE(2008)

                DO c5=1,n_mix
                  DO c4=1,n_tot_bands
                    DO c3=1,n_terms
                      DO c2=1,n_press
                        DO c1=1,n_temps

! desired ordering kval_array ( temps, press, mix, wl_band, term )
                          kval_array(c1,c2,c5,c4,c3) = old_kval_array(c1,c2,c3,c4,c5)

                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO

                DEALLOCATE(old_kval_array)

                               band_wn_solar = (/ 2222.22, 3087.37, 4030.63, 5370.57,           &
                                                  7651.11, 12500.0, 25000.0, 41666.67 /)

                                  band_wn_ir = (/ 10.0, 166.667, 416.667, 625.370, 710.020,     &
                                                  833.333, 1250.0, 2222.22 /)

                IF(.NOT.ALLOCATED(kval_temps)) ALLOCATE(kval_temps(n_temps))
                                  kval_temps = (/ 50.,  70.,  90.,  110., 130., 150., 170., 190., 210., &
                                                  230., 250., 270., 290., 310., 330., 350., 370., 390. /) ! temperature values in kvals array


                IF(.NOT.ALLOCATED(kval_press)) ALLOCATE(kval_press(n_press))
                             DO i=1,n_press
                                  kval_press(i)=10**(-4.+(i-1)*0.4)  ! p_spacing = 0.4
                             ENDDO
!                                  kval_press = (/ 1.e-4, 2.5119e-4, 6.3096e-4, 1.5849e-3, 3.9822e-3,    &
!                                                  1.e-2, 2.5119e-2, 6.3096e-2, 1.5849e-1, 3.9822e-1,    &
!                                                  1.e-0, 2.5119e-0, 6.3096e-0, 1.5849e+1, 3.9822e+1,    &
!                                                  1.e+2, 2.5119e+2, 6.3096e+2, 1.5849e+3, 3.9822e+3,    &
!                                                  1.e+4, 2.5119e+4, 6.3096e+4, 1.5849e+5, 3.9822e+5,    &
!                                                  1.e+6                                                 /) ! pressure values in kvals array

                  IF(.NOT.ALLOCATED(kval_mix)) ALLOCATE(kval_mix(n_mix))
                               kval_mix(1)=1.0D-25
                               DO i=2,n_mix
                                    kval_mix(i)=10**(-7.+(i-2)*1.0)    ! mix_spacing = 1.0
                               ENDDO
!                                    kval_mix = (/ 1.e-25, 1.e-7, 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1 /) ! mixing ratio values
#else
        CALL wrf_error_fatal("...KDM init trying to read old style kvals, but only valid for Mars")
#endif

      ELSE   ! is new style k-table file

        call wrf_message(new_line('a')//"KDM is assuming you passed it a new format kvals file")

        OPEN(UNIT=2008, FILE=TRIM(k_coeff_file))

        READ(2008,*) target_object
          WRITE(wrf_err_message,*) "...reading k-table for ",trim(adjustl(target_object))
          CALL wrf_message(TRIM(wrf_err_message))

        READ(2008,*) header1
        READ(2008,*) n_temps,n_press,n_terms,n_in_band,n_mix,n_in_tot
         
          WRITE(wrf_err_message,*)   " n_temps=",n_temps,new_line('a'), &
                                   "   n_press=",n_press,new_line('a'), &
                                   "   n_terms=",n_terms,new_line('a'), &
                                   "   n_mix  =",n_mix,new_line('a')
          CALL wrf_message(wrf_err_message)

          IF( n_tot_bands /= n_in_band ) THEN
              WRITE(wrf_err_message,*) "KDM radtran: KVAL table error:",new_line('a'), &
                                       " unexpected number of wavelength bands in file",new_line('a'), &
                                       " expected=",n_tot_bands," got=",n_in_band
              CALL wrf_error_fatal(wrf_err_message)
          ENDIF

          ALLOCATE(wavenumbers(n_in_band+1))

          ALLOCATE(    kval_array(n_temps,n_press,n_mix,n_tot_bands,n_terms)) ! Full array of k-coefficients

          IF(.NOT.ALLOCATED(kval_temps)) ALLOCATE(kval_temps(n_temps))
          IF(.NOT.ALLOCATED(kval_press)) ALLOCATE(kval_press(n_press))
          IF(.NOT.ALLOCATED(kval_mix))   ALLOCATE(kval_mix(n_mix))
          IF(.NOT.ALLOCATED(weights))    ALLOCATE(weights(n_terms))

        READ(2008,*) header1
        READ(2008,*) kval_temps
        READ(2008,*) header1
        READ(2008,*) kval_press

        READ(2008,*) header1
        READ(2008,*) weights   ! these are NOT the weights but are the abscissae (cumulative fraction)

        READ(2008,*) header1
        READ(2008,*) weights   ! weights are read in with new-style kvals tables

        READ(2008,*) header1
        READ(2008,*) wavenumbers

        band_wn_ir    = wavenumbers(1:n_band_ir+1)/100.           ! divide by 100 since file has wn
        band_wn_solar = wavenumbers(n_band_ir+1:n_in_band+1)/100. ! in m-1 and we want it in cm-1

        DEALLOCATE(wavenumbers)

        READ(2008,*) header1
        READ(2008,*) kval_mix

        write(wrf_err_message,*) " Temps=",kval_temps
        CALL wrf_message(wrf_err_message)
        write(wrf_err_message,*) " Press=",kval_press
        CALL wrf_message(wrf_err_message)
        write(wrf_err_message,*) " Mix=  ",kval_mix
        CALL wrf_message(wrf_err_message)
        write(wrf_err_message,*) " WN IR= ",band_wn_ir
        CALL wrf_message(wrf_err_message)
        write(wrf_err_message,*) " WN VIS=",band_wn_solar
        CALL wrf_message(wrf_err_message)

        read(2008,*) header2

        do c5=1,n_mix
          do c4=1,n_tot_bands
            do c3=1,n_terms
              do c2=1,n_press
                do c1=1,n_temps
                  ! desired ordering kval_array ( temps, press, mix, wl_band, term )
                  read(2008,*) r_in_t,r_in_p,r_in_term,r_in_band,r_in_mix,kval_array(c1,c2,c5,c4,c3)
                enddo
              enddo
            enddo
          enddo
        enddo
        CLOSE(2008)

      ENDIF ! end of if for which type of k-table file

      ! set the min and max values for the k-table dimensions:
      kval_press_min=kval_press(1)
      kval_press_max=kval_press(n_press)
      kval_mix_min=kval_mix(1)
      kval_mix_max=kval_mix(n_mix)
      kval_temps_min=kval_temps(1)
      kval_temps_max=kval_temps(n_temps)
      ! pre-convert the kval pressure and mixing ratio axes into log values for interpolation in log space:
      kval_press = log(kval_press)
      kval_mix   = log(kval_mix)

      ! check that kval_array is monotonically increasing in temp, pres, and mix axes 
      ! - get_kvals will not work properly otherwise
      do c5=1,n_mix-1
        if(kval_mix(c5) .ge. kval_mix(c5+1)) CALL wrf_error_fatal("kval mix axes not monotonic increasing")
      enddo
      do c2=1,n_press-1
        if(kval_press(c2) .ge. kval_press(c2+1)) CALL wrf_error_fatal("kval press axes not monotonic increasing")
      enddo
      do c1=1,n_temps-1
        if(kval_temps(c1) .ge. kval_temps(c1+1)) CALL wrf_error_fatal("kval temps axes not monotonic increasing")
      enddo

      CALL wrf_message(" ")
      CALL wrf_message("Namelist activated raditive aerosol options (ra_include_...) for this run:")
      CALL wrf_message(" (if these aren't what you want, change it in the namelist)")
      CALL wrf_message(" ")
      CALL wrf_message(" You requested: ")

#if ( WRF_MARS == 1 )
      IF(include_dust) THEN
        CALL wrf_message (" --> radiatively active dust")
      ELSE
        CALL wrf_message (" --> radiatively inactive dust - DUST *NOT* ACTIVE!")
      ENDIF

      IF(include_water_ice) THEN
        CALL wrf_message (" --> radiatively active water ice clouds")
      ELSE
        CALL wrf_message (" --> radiatively inactive water ice clouds - H2O ICE *NOT* ACTIVE!")
      ENDIF

      IF(include_co2_ice) THEN
        CALL wrf_message (" --> radiatively active co2 ice clouds")
      ELSE
        CALL wrf_message (" --> radiatively inactive co2 ice clouds - CO2 ICE *NOT* ACTIVE!")
      ENDIF

      IF(include_water_vap) THEN
        CALL wrf_message (" --> radiatively active water vapour")
      ELSE
        CALL wrf_message (" --> radiatively inactive water vapour - WATER VAPOUR *NOT* ACTIVE!")
      ENDIF
#elif ( WRF_TITAN == 1 )
      IF(include_ch4_vap) THEN
        CALL wrf_message (" --> radiatively active methane vapour")
      ELSE
        CALL wrf_message (" --> radiatively inactive methane vapour - CH4 GAS *NOT* ACTIVE!")
      ENDIF

      IF(include_tholin) THEN
        CALL wrf_message (" --> radiatively active tholin haze")
      ELSE
        CALL wrf_message (" --> radiatively inactive tholin - THOLIN HAZE *NOT* ACTIVE!")
      ENDIF
#endif
      CALL wrf_message(" ")

#if ( WRF_MARS == 1 )
      write(wrf_err_message,*) "KDM init thinks dust 2mom:",is_two_moment_dust," ice 2mom:",is_two_moment_water_ice, &
                    " kdm_options:",kdm_options, ko_multimie_dust, btest(kdm_options, ko_multimie_dust)
      call wrf_message(trim(wrf_err_message))

      IF (include_dust) THEN
         IF (is_two_moment_dust .and. btest(kdm_options, ko_multimie_dust)) THEN
            CALL wrf_message ("...reading in dust aerosol parameters for two moment scheme")
            CALL read_multisize_miedata(trim(dust_params_file),1)
         ELSE 
            CALL wrf_message ("...reading in dust aerosol parameters for non-moment scheme")
            CALL read_singlesize_miedata(trim(dust_params_file),1,n_dust_data, &
                    wvl_in_d,cext_in_d,cscat_in_d,g_in_d,qext_in_d,omega_in_d, &
                    ref_wvl_dust, ref_wvl_ir_dust )
         ENDIF
      ENDIF

      IF (include_water_ice) THEN
         IF (is_two_moment_water_ice .and. btest(kdm_options, ko_multimie_water)) THEN
            CALL wrf_message ("...reading in ice aerosol parameters for two moment scheme")
            CALL read_multisize_miedata(trim(ice_params_file),2)
         ELSE 
            CALL wrf_message ("...reading in water ice aerosol parameters for non-moment scheme")
            CALL read_singlesize_miedata(trim(ice_params_file),2,n_ice_data,   &
                    wvl_in_i,cext_in_i,cscat_in_i,g_in_i,qext_in_i,omega_in_i, &
                    ref_wvl_water_ice )
         ENDIF
      ENDIF

      IF (include_co2_ice) THEN
!         IF (is_two_moment .and. btest(kdm_options, ko_multimie_co2ice)) THEN
!             WRITE(0,*) "...reading in co2ice aerosol parameters with two moment scheme"
!                CALL read_multisize_miedata(trim(ice_params_file),3)
!         ELSE 
            CALL wrf_message ("...reading in CO2 ice aerosol parameters")
            CALL read_singlesize_miedata(trim(co2ice_params_file),3,n_co2ice_data, &
                    wvl_in_c,cext_in_c,cscat_in_c,g_in_c,qext_in_c,omega_in_c,     &
                    ref_wvl_co2_ice )
!         ENDIF
      ENDIF
#endif
#if ( WRF_TITAN == 1 )
      IF(include_tholin) THEN

        IF(mp_physics_tholin == mckay_tholin_haze) THEN
          call read_mckay_tholin()
        ELSE
          write(wrf_err_message,*) "KDM: invalid mp_physics_tholin choice for radiative haze, ",mp_physics_tholin
          call wrf_error_fatal(trim(wrf_err_message))
        ENDIF

      ENDIF
#endif

      IF ( old_kvals_file_style ) THEN
      CALL prep_gaussian_weights(n_terms,l_double,weights)               ! This calculates the weighting terms, which
                                                                         !    depend only on n_terms.  Do this ONCE
      ENDIF

      CALL solar_fraction(band_wn_solar(1),                            & ! Argument #3 ('1.0D0') lets this code calculate the
                          band_wn_solar(n_band_solar+1),1.0D0,         & !    total solar flux the first time through so
                          total_solar_flux)                              !    the subsequent loop can use total_solar_flux
                        
      DO i=1,n_band_solar                                                !    to calculate the fractional flux in each band.
         CALL solar_fraction(band_wn_solar(i),band_wn_solar(i+1),      & ! This calculates the solar fraction in each
                             total_solar_flux,solar_flux_band(i))        !    spectral band.  Do this ONCE.
         CALL rayleigh_coeffs(band_wn_solar(i),band_wn_solar(i+1),     & ! This calculates the rayleigh scattering cross-section
                              rayleigh_coefficient(i))                   !    in each spectral band.  Do this ONCE.
      ENDDO

#if ( WRF_MARS == 1 )
      IF ((include_dust)                                .and. &
          (.not.(btest(kdm_options, ko_multimie_dust))) .and. &
          (.not.is_two_moment_dust)) THEN
! dust solar processing
            CALL process_aerosol(band_wn_solar,n_band_solar,wvl_in_d,omega_in_d, &
                              g_in_d,cext_in_d,6000.D0,omega_dust_solar,         &
                              asymmetry_dust_solar,kCext_dust_solar,             &
                              l_simpledust,n_dust_data,                          &
                              simpledust_omega, simpledust_asymmetry)
! dust thermal IR processing
            CALL process_aerosol(band_wn_ir,n_band_ir,wvl_in_d,omega_in_d,g_in_d, &
                              cext_in_d,215.D0,omega_dust_ir,                     &
                              asymmetry_dust_ir,kCext_dust_ir,                    &
                              .FALSE.,n_dust_data)

         deallocate(wvl_in_d)
         deallocate(cext_in_d)
         deallocate(cscat_in_d)
         deallocate(g_in_d)
         deallocate(qext_in_d)
         deallocate(omega_in_d)
      ELSE
         omega_dust_solar = 0
         asymmetry_dust_solar = 0
         kCext_dust_solar = 0
         omega_dust_ir = 0
         asymmetry_dust_ir = 0
         kCext_dust_ir = 0
      ENDIF
       
      IF ((include_water_ice)                            .and. &
          (.not.(btest(kdm_options, ko_multimie_water))) .and. &
          (.not.is_two_moment_water_ice)) THEN
! ice solar processing
               CALL process_aerosol(band_wn_solar,n_band_solar,wvl_in_i,omega_in_i, &
                                 g_in_i,cext_in_i,6000.D0,omega_ice_solar,          &
                                 asymmetry_ice_solar,kCext_ice_solar,               &
                                 .FALSE.,n_ice_data)
! ice thermal IR processing
               CALL process_aerosol(band_wn_ir,n_band_ir,wvl_in_i,omega_in_i,g_in_i, &
                                 cext_in_i,215.D0,omega_ice_ir,                      &
                                 asymmetry_ice_ir,kCext_ice_ir,                      &
                                 .FALSE.,n_ice_data)
         deallocate(wvl_in_i)
         deallocate(cext_in_i)
         deallocate(cscat_in_i)
         deallocate(g_in_i)
         deallocate(qext_in_i)
         deallocate(omega_in_i)
      ELSE
         omega_ice_solar = 0
         asymmetry_ice_solar = 0
         kCext_ice_solar = 0
         omega_ice_ir = 0
         asymmetry_ice_ir = 0
         kCext_ice_ir = 0
      ENDIF

      IF (include_co2_ice) THEN
! CO2 ice solar processing
               CALL process_aerosol(band_wn_solar,n_band_solar,wvl_in_c,omega_in_c, &
                                 g_in_c,cext_in_c,6000.D0,omega_co2ice_solar,       &
                                 asymmetry_co2ice_solar,kCext_co2ice_solar,         &
                                 .FALSE.,n_co2ice_data)
! CO2 ice thermal IR processing
               CALL process_aerosol(band_wn_ir,n_band_ir,wvl_in_c,omega_in_c,g_in_c, &
                                 cext_in_c,215.D0,omega_co2ice_ir,                   &
                                 asymmetry_co2ice_ir,kCext_co2ice_ir,                &
                                 .FALSE.,n_co2ice_data)
         deallocate(wvl_in_c)
         deallocate(cext_in_c)
         deallocate(cscat_in_c)
         deallocate(g_in_c)
         deallocate(qext_in_c)
         deallocate(omega_in_c)
      ELSE
         omega_co2ice_solar = 0
         asymmetry_co2ice_solar = 0
         kCext_co2ice_solar = 0
         omega_co2ice_ir = 0
         asymmetry_co2ice_ir = 0
         kCext_co2ice_ir = 0
      ENDIF

      if(do_chicago_ice_layers) then
        call init_ice_tau_chicago(chi_cloud_prs_in, chi_ncloud_in,          &
                                  chi_ls1_in, chi_ls2_in, chi_lat1_in,      &
                                  chi_lat2_in, chi_lon1_in, chi_lon2_in,    &
                                  chi_tval1_in, chi_tval2_in,               &
                                  chi_dayornight_in, chi_ice_tau_in         )
      endif
#endif

      END SUBROUTINE kdminit

      SUBROUTINE read_singlesize_miedata(params_file,itype,n_data, &
                    wvl_in,cext_in,cscat_in,g_in,qext_in,omega_in, &
                    ref_wvl, ref_wvl_ir )

            IMPLICIT NONE

            CHARACTER(LEN=*), INTENT(IN)       ::         params_file

            !ins
            INTEGER, INTENT(IN) ::                              itype

            REAL(KIND(0.d0)), INTENT(IN)                ::    ref_wvl
            REAL(KIND(0.d0)), INTENT(IN), OPTIONAL      :: ref_wvl_ir

            !outs
            REAL(KIND(0.d0)), DIMENSION(:), INTENT(OUT),               &
                                      ALLOCATABLE ::           wvl_in, &
                                                              cext_in, &
                                                             cscat_in, &
                                                                 g_in, &
                                                              qext_in, &
                                                             omega_in

            INTEGER, INTENT(OUT) ::                            n_data

            !locals
            real(kind(0.d0)) :: d1, d2, d3, d4, d5, kCext_ref, kCext_ref_ir
            character(len=20) :: c_type
            integer :: file_status, i

            select case (itype)
              case(1)
                c_type = "dust"
              case(2)
                c_type = "water ice"
              case(3)
                c_type = "co2 ice"
              case default
                WRITE (wrf_err_message,*) "Unrecognized aerosol index in read_singlesize_miedata: ",itype
                CALL wrf_error_fatal(trim(wrf_err_message))
            end select
            n_data = 0
            OPEN(UNIT=2009, FILE=TRIM(params_file))   ! open and read to end to count number of lines
            DO
              IF(itype == 1) THEN
                      READ(2009,*,IOSTAT=file_status)  d1, d2, d3, d4, d5
              ELSE
                      READ(2009,*,IOSTAT=file_status)  d1, d2, d3, d4
              ENDIF
              IF ((file_status /= 0) .and. (n_data == 0))  THEN
                WRITE (wrf_err_message,*) "Error reading ",trim(c_type)," parameter file: ", trim(params_file)
                CALL wrf_error_fatal(trim(wrf_err_message))
              ELSE IF (file_status /= 0) THEN  ! file_status -ve for true EOF, but some files have +ve for text junk at the end
                EXIT  !... end of file reached ...
              ELSE
                n_data = n_data + 1
              END IF
            END DO
            CLOSE(2009)

            if(allocated(wvl_in)) then
              CALL wrf_error_fatal("read_singlesize_miedata: wvl_in already allocated for "//trim(c_type))
             else
              allocate(wvl_in  (n_data))
            endif

            if(allocated(cext_in)) then
              CALL wrf_error_fatal("read_singlesize_miedata: cext_in already allocated for "//trim(c_type))
             else
              allocate(cext_in (n_data))
            endif

            if(allocated(cscat_in)) then
              CALL wrf_error_fatal("read_singlesize_miedata: cscat_in already allocated for "//trim(c_type))
             else
              allocate(cscat_in(n_data))
            endif

            if(allocated(g_in)) then
              CALL wrf_error_fatal("read_singlesize_miedata: g_in already allocated for "//trim(c_type))
             else
              allocate(g_in    (n_data))
            endif

            if(allocated(qext_in)) then
              CALL wrf_error_fatal("read_singlesize_miedata: qext_in already allocated for "//trim(c_type))
             else
              allocate(qext_in (n_data))
            endif

            if(allocated(omega_in)) then
              CALL wrf_error_fatal("read_singlesize_miedata: omega_in already allocated for "//trim(c_type))
             else
              allocate(omega_in(n_data))
            endif

            OPEN(UNIT=2009, FILE=TRIM(params_file))
            IF(itype == 1) THEN  ! read the first line
              READ(2009,*) wvl_in(1),cext_in(1),cscat_in(1),g_in(1),qext_in(1)
              omega_in(1)=cscat_in(1)/cext_in(1)
            ELSE
              READ(2009,*) wvl_in(1),cext_in(1),omega_in(1),g_in(1)
            ENDIF

            DO i=2,n_data
               IF(itype == 1) THEN
                 READ(2009,*) wvl_in(i),cext_in(i),cscat_in(i),g_in(i),  &
                              qext_in(i)
                 omega_in(i)=cscat_in(i)/cext_in(i)
               ELSE
                 READ(2009,*) wvl_in(i),cext_in(i),omega_in(i),g_in(i)
               ENDIF
#if ( WRF_MARS == 1 )
               IF ((wvl_in(i) > ref_wvl).AND.(wvl_in(i-1) <= ref_wvl))    & ! We want to pull out the value of kCext at the 
                      kCext_ref=(wvl_in(i)-ref_wvl)/(wvl_in(i)-           & !    reference wavelength (0.67 um).  The approach
                                 wvl_in(i-1))*cext_in(i-1)+(ref_wvl-      & !    used here searches for the appropriate bounds
                                 wvl_in(i-1))/(wvl_in(i)-wvl_in(i-1))*    & !    and interpolates between them to 0.67 um.
                                 cext_in(i)
                      if(itype==1) kCext_dust_ref   = kCext_ref
                      if(itype==2) kCext_ice_ref    = kCext_ref
                      if(itype==3) kCext_co2ice_ref = kCext_ref
               IF(present(ref_wvl_ir) .and. (itype==1)) THEN
                IF ((wvl_in(i) > ref_wvl_ir).AND.(wvl_in(i-1) <= ref_wvl_ir))  &
                   kCext_ref_ir=(wvl_in(i)-ref_wvl_ir)/(wvl_in(i)-             & !    reference wavelength (9.3 um).  The approach
                                 wvl_in(i-1))*cext_in(i-1)+(ref_wvl_ir-        & !    used here searches for the appropriate bounds
                                 wvl_in(i-1))/(wvl_in(i)-wvl_in(i-1))*         & !    and interpolates between them to 9.3 um.
                                 cext_in(i)
                   kCext_dust_ref_ir = kCext_ref_ir
               ENDIF
#endif
            ENDDO
            CLOSE(2009)

      END SUBROUTINE read_singlesize_miedata


#if ( WRF_TITAN == 1 )
      SUBROUTINE read_mckay_tholin()

            IMPLICIT NONE

            CHARACTER(LEN=41) :: params_file="./Data/radiation/mckay_tholin_scatter.txt"

            integer, save :: n_radii, n_wvls
            real, dimension(:),   allocatable :: mck_wvl, mck_rad
            real, dimension(:,:), allocatable :: mck_qext, mck_ssa, mck_asym

            character(len=3) :: cjunk  ! header junk
            integer :: file_status, i, j
            real :: a,b,c,d,e,f

            OPEN(UNIT=2009, FILE=TRIM(params_file))   ! open and read to end to count number of lines

            read(2009,*,IOSTAT=file_status) cjunk
            IF (file_status /= 0 )  THEN
               WRITE (wrf_err_message,*) "Error reading tholin parameter file: ", trim(params_file)
               CALL wrf_error_fatal(trim(wrf_err_message))
            END IF

            read(2009,*) n_wvls, n_radii
            if((n_wvls /= 381) .or. (n_radii /= n_tholin_rad_bins )) then
               WRITE (wrf_err_message,*) "expected 381 wavelengsth and ",n_tholin_rad_bins, &
                                         " radii in table, but got: ", n_wvls, n_radii
               CALL wrf_error_fatal(trim(wrf_err_message))
            endif

            read(2009,*,IOSTAT=file_status) cjunk
            IF (file_status /= 0 )  THEN
               WRITE (wrf_err_message,*) "Error reading tholin parameter file: ", trim(params_file)
               CALL wrf_error_fatal(trim(wrf_err_message))
            END IF

            allocate(mck_wvl(n_wvls))
            allocate(mck_rad(n_radii))

            allocate(mck_qext(n_wvls,n_radii))
            allocate(mck_ssa(n_wvls,n_radii))
            allocate(mck_asym(n_wvls,n_radii))

            do i=1,n_wvls
              do j=1,n_radii

                read(2009,*,IOSTAT=file_status) a,b,c,d,e
                
                mck_rad(j)    = b
                mck_qext(i,j) = c
                mck_ssa(i,j)  = d
                mck_asym(i,j) = e
              enddo

              mck_wvl(i)    = a
            enddo

            CLOSE(2009)

            tholin_rad_bins = mck_rad

            do j=1,n_radii
              call mckay_wavelength_int(band_wn_solar,n_band_solar,n_wvls,         &
                                        dble(mck_wvl), dble(mck_qext(:,j)),        &
                                        dble(mck_ssa(:,j)),dble(mck_asym(:,j)),    &
                                        6000.D0,                                   &
                                        tholin_vis_ext_xsect(:,j),                 & 
                                        tholin_vis_ssa(:,j), tholin_vis_asym(:,j)  )

              call mckay_wavelength_int(band_wn_ir,n_band_ir,n_wvls,               &
                                        dble(mck_wvl), dble(mck_qext(:,j)),        &
                                        dble(mck_ssa(:,j)),dble(mck_asym(:,j)),    &
                                        100.D0,                                    &
                                        tholin_ir_ext_xsect(:,j),                  & 
                                        tholin_ir_ssa(:,j), tholin_ir_asym(:,j)    )
            enddo

            deallocate(mck_wvl)
            deallocate(mck_rad)
            deallocate(mck_qext)
            deallocate(mck_ssa)
            deallocate(mck_asym)

      END SUBROUTINE read_mckay_tholin


      SUBROUTINE mckay_wavelength_int(band_wn,n_band,n_aerosol_data,wvl_in, &
                                      cext_in,omega_in,g_in, planck_temp,   &
                                      kCext_out, omega_out, asymmetry_out   )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                n_band, &
                                                       n_aerosol_data

      REAL(KIND(0.d0)), INTENT(IN   ) ::                   planck_temp

      REAL(KIND(0.d0)), INTENT(IN   ), DIMENSION(n_band+1) ::          &
                                                               band_wn

      REAL(KIND(0.d0)), INTENT(IN   ), DIMENSION(n_aerosol_data) ::    &
                                                               wvl_in, &
                                                             omega_in, &
                                                                 g_in, &
                                                              cext_in

      REAL(KIND(0.d0)), INTENT(  OUT), DIMENSION(n_band) ::            &
                                                            omega_out, &
                                                        asymmetry_out, &
                                                            kCext_out
!c-mm  Local variables

      INTEGER ::                                                  i,k, &
                                                                 kmin, &
                                                                 kmax

      REAL(KIND(0.d0)) ::                                         tmp, &
                                                                 dtmp, &
                                                                 ntmp

      REAL(KIND(0.d0)), DIMENSION(n_band+1) ::                band_wl


      DO i=1,n_band+1
          band_wl(i)=1.d4/band_wn(n_band+2-i)
      ENDDO

      kmin=1
      kmax=1
      DO i=1,n_band

        DO k=kmin,n_aerosol_data
           IF (wvl_in(k) > band_wl(i)) THEN
              kmin=MAX(k-1,1)
              EXIT
           ENDIF
        ENDDO

        DO k=kmin,n_aerosol_data
           IF (wvl_in(k) > band_wl(i+1)) THEN
              kmax=k
              EXIT
           elseif (k == n_aerosol_data) then
              kmax = k
           ENDIF
        ENDDO

        !c-mm  Approximate integrals with trapezoidal integration

        !c-mm  Do Cext first
        dtmp = 0.0D0
        DO k=kmin,kmax-1
           tmp=0.5D0*( planck(wvl_in(k),planck_temp) +                 &
                       planck(wvl_in(k+1),planck_temp) )*              &
                       (wvl_in(k+1)-wvl_in(k))
           dtmp = dtmp+tmp
        ENDDO
        ntmp=0.0D0
        DO k=kmin,kmax-1
           tmp=0.5D0*( cext_in(k)*planck(wvl_in(k),planck_temp) +      &
                       cext_in(k+1)*planck(wvl_in(k+1),planck_temp) )* &
                       (wvl_in(k+1)-wvl_in(k))
           ntmp=ntmp+tmp
        ENDDO
        kCext_out(n_band-i+1)=ntmp/dtmp                                 ! Put the value back in its correct 'wavenumber' index

        !c-mm  Do omega next.  Denominator of omega is the same as
        !c-mm     numerator of Cext.
        dtmp=ntmp
        ntmp=0.D0
        DO k=kmin,kmax-1
           tmp=0.5D0*( omega_in(k)*planck(wvl_in(k),planck_temp)*      &
               cext_in(k)+omega_in(k+1)*                               &
               planck(wvl_in(k+1),planck_temp)*cext_in(k+1) )*         &
               (wvl_in(k+1)-wvl_in(k))
           ntmp=ntmp+tmp
        ENDDO
        omega_out(n_band-i+1)=ntmp/dtmp

        !c-mm  Do asymmetry next.  Denominator of asymmetry is the 
        !c-mm     same as numerator of omega
        dtmp=ntmp
        ntmp=0.D0
        DO k=kmin,kmax-1
           tmp=0.5D0*( g_in(k)*omega_in(k)*                            &
               planck(wvl_in(k),planck_temp)*cext_in(k)+g_in(k+1)*     &
               omega_in(k+1)*planck(wvl_in(k+1),planck_temp)*          &
               cext_in(k+1) )*(wvl_in(k+1)-wvl_in(k))
           ntmp=ntmp+tmp
        ENDDO
        asymmetry_out(n_band-i+1)=ntmp/dtmp

      ENDDO

      END SUBROUTINE mckay_wavelength_int
#endif

 
!-----------------------------------------------------------------------
!+ Function to calculate the saturation vapor pressure at a given      !
!  temperature.                                                        !
!                                                                      !
! Method:                                                              !
!       Straightforward.  Uses Clausius-Clapeyron equation.            !
!-----------------------------------------------------------------------
      FUNCTION sat_vap(t)

#ifdef mpas
      USE mpas_atmphys_constants, only : xlv, r_v 
#else
      USE module_model_constants, ONLY : xlv, r_v
#endif

      IMPLICIT NONE

      REAL(KIND(0.d0)), INTENT(IN   ) ::                            t

! Local variables

      REAL(KIND(0.d0)) ::                                     sat_vap

      sat_vap=611.*EXP(xlv/r_v*(3.661D-03-(1/t)))

      END FUNCTION sat_vap

!-----------------------------------------------------------------------
!+ Subroutine to prepare the domain for the correct quadrature method  !
!     (single or double gaussian), followed by a call to               !
!     get_gaussian_weights                                             !
!  points.                                                             !
!                                                                      !
! Method:                                                              !
!    Registry variable l_double determines whether we are using a      !
!       single or double gaussian quadrature.  Depending on the choice,!
!       we call get_gaussian_weights once or twice, with the correct   !
!       domain bounds (x1,x2) for the particular gaussian.             !
!-----------------------------------------------------------------------
      SUBROUTINE prep_gaussian_weights(num_terms,l_double,weights)

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) :: num_terms

      LOGICAL, INTENT(IN   ) :: l_double

      REAL(KIND(0.d0)), INTENT(   OUT) :: weights(num_terms)

! Local variables

      INTEGER :: i

      REAL(KIND(0.d0)) :: x1, x2

!c-mm      REAL(KIND(0.d0)) :: weight_block(num_terms)
      REAL(KIND(0.d0)), ALLOCATABLE :: weight_block(:)

! Initialize weights and set other initial values
      weights=0.0
      SELECT CASE (l_double)
         CASE (.FALSE.)
            x2=1.0
            x1=0.0
            ALLOCATE(weight_block(num_terms))
            CALL get_gaussian_weights(num_terms,x1,x2,weight_block)
            weights(1:num_terms)=weight_block
            DEALLOCATE(weight_block)
         CASE (.TRUE.)
            DO i=1,2
               SELECT CASE (i)
                  CASE (1)
                     x2=0.95
                     x1=0.0
                  CASE (2)
                     x2=1.0
                     x1=0.95
                  CASE DEFAULT
                   WRITE (wrf_err_message,*) &
                      'Crashed in get_gaussian_weights'
                   CALL wrf_error_fatal ( wrf_err_message )
               END SELECT
               ALLOCATE(weight_block(num_terms/2))
               CALL get_gaussian_weights(num_terms/2,x1,x2,            &
                                         weight_block)
               weights((i-1)*(num_terms/2)+1:i*num_terms/2)=           &
                                                          weight_block
               DEALLOCATE(weight_block)
            ENDDO
         CASE DEFAULT
            WRITE (wrf_err_message,*) &
               'Crashed in get_gaussian_weights'
          CALL wrf_error_fatal ( wrf_err_message )
      END SELECT

      END SUBROUTINE prep_gaussian_weights


!-----------------------------------------------------------------------
!+ Subroutine to get gaussian weights given a number of quadrature     !
!  points.                                                             !
!                                                                      !
! Method:                                                              !
!       Has roots in Numerical Recipes, and was originally coded by    !
!          Mike Smith at GSFC                                          !
!-----------------------------------------------------------------------

      SUBROUTINE get_gaussian_weights(num_terms,x1,x2,weights)

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) :: num_terms

      REAL(KIND(0.d0)), INTENT(IN   ) :: x1, x2

      REAL(KIND(0.d0)), INTENT(   OUT) :: weights(num_terms)

! Local variables
      REAL(KIND(0.d0)), PARAMETER ::                  small = 5.0D-10, &
                                                      pi = 3.14159D+0
      INTEGER  :: m, j, k
      REAL(KIND(0.d0)) :: xm, xl, z, z1, p1, p2, p3, pp, rn
      m=INT((num_terms+1)/2)
      xm=(x2+x1)/2.
      xl=(x2-x1)/2.
      rn=REAL(num_terms)

      DO j=1,m
         z=COS(PI*(REAL(j)-0.25)/(rn+0.5)) 
         z1=1.0D6
         DO WHILE (ABS(z-z1)>=small)
            p1=1.0
            p2=0.0
            DO k=1,num_terms
               p3=p2
               p2=p1
               p1=(REAL(2*k-1)*z*p2-REAL(k-1)*p3)/REAL(k)
            END DO
            pp=rn*(z*p1-p2)/(z*z-1.0)
            z1=z
            z=z1-p1/pp
         END DO
         weights(j)=2.0*xl/((1.0-z*z)*pp*pp)
         weights(num_terms+1-j)=weights(j)
      END DO

      END SUBROUTINE get_gaussian_weights


!-----------------------------------------------------------------------
!+ Subroutine calculating the fraction of blackbody radiation          !
!  between two wavelengths, assuming this is for solar wavelengths.    ! 
!                                                                      !
! Method:                                                              !
!        Uses series formulae from Houghton 'Physics of Atmospheres'   !
!-----------------------------------------------------------------------
      SUBROUTINE solar_fraction(wl1, wl2, denom, fst4)

!c-mm      USE module_ra_mars_wbm                                              ! For call to planck_function

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


      IF (wl1 == 0.0) THEN						 ! Swapping of wl's going from wavenumber to wavelength.
	   wl(2)=1000.
      ELSE
           wl(2)=10000./wl1
      ENDIF
      wl(1)=10000./wl2
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
!+ Subroutine calculating the Rayleigh scattering cross-section for a  !
!    martian atmosphere composed of CO2 and N2.  The output of this    !
!    program is the scattering coefficient in m^2/kg, calculated at    !
!    the center of each band. Everything is intially done in um before !
!    the final conversion into m.                                      !
!                                                                      !
! Method:                                                              !
!        Uses standard equation for Rayleigh scattering, as found in   !
!           Vardavas and Taylor (2007) 'Radiation and Climate' Eq. 6.68!
!        Values for a and b are from Allen's Astrophysical Quantities  !
!           and depol_fac are taken from Young (1980), Appl. Opt. 19   !
!-----------------------------------------------------------------------
      SUBROUTINE rayleigh_coeffs(wn1,wn2,rayleigh_coefficient)

#ifdef mpas
      USE mpas_atmphys_constants, only : co2_mixing_ratio, n2_mixing_ratio
#else
      USE module_model_constants, ONLY: co2_mixing_ratio, n2_mixing_ratio
#endif

      IMPLICIT NONE

      REAL(KIND(0.d0)), INTENT(IN) ::                             wn1, &
                                                                  wn2

      REAL(KIND(0.d0)), INTENT(OUT) ::           rayleigh_coefficient

      REAL(KIND(0.d0)), DIMENSION(2) ::                                &
                                   depol_fac = (/0.078, 0.021/),       & ! First element is for CO2, second
                                           a = (/4.390D-4, 2.906D-4/), & !    is for N2
                                           b = (/6.4D-3, 7.7D-3/),     &
#if ( WRF_MARS == 1 )
                                      volmix = (/co2_mixing_ratio,     &
                                              (1.-co2_mixing_ratio)/), &
#elif ( WRF_TITAN == 1 )
                                      volmix = (/0., 1./),             & ! Need values for CH4 to do correctly
#else
                                      ! Dummy variables for other WRF cases.
                                      volmix = (/0., 0./),             &
#endif
                                      mol_wt = (/44.D-3, 14.D-3/)        ! In kg/mol

      REAL(KIND(0.d0)), PARAMETER ::                     pi = 3.14159, &
                                                        Na = 6.022D23, & ! Avogadro's number in molecules/mol
                                                         Lo = 2.687D7, & ! Loschmidt's number in molecules/um^3
                         constant = (128/3)*pi**5*((1./(2*pi*Lo))**2), & ! constant (=4.58D-13 um^6/molecule^2)
                                                    um2_to_m2 = 1.D12    ! um2_to_m2 converts from um^2 to m^2   

      INTEGER ::                                                 band, &
                                                                    i

      REAL(KIND(0.d0)) ::                                          wl, &
                                                                 wl2i, &
                                                                  sum, &
                                                                aniso, &
                                                                alpha, &
                                                           avg_mol_wt

#if ( WRF_MARS != 1 ) && ( WRF_TITAN != 1 )
      call wrf_error_fatal("Trying to setup Rayleigh damping, but only defined for Mars and Titan")
#endif

      wl=(1.D4/wn1+1.D4/wn2)/2.                                          ! Take average waveLENGTH across band
      wl2i=1./wl/wl
      sum=0.0
      DO i=1,2                                                           ! Loop over both gases
         aniso = (6.D0+3.D0*depol_fac(i))/(6.D0-7.D0*depol_fac(i))       ! Anisotropy factor, sometimes known as 'delta'
         alpha = (a(i)*(1.+b(i)*wl2i))                                   ! A 'quasi' alpha term, since (1/(2*pi*Lo))^2 is folded into 'constant'.
         sum = sum+volmix(i)*aniso*alpha*alpha
      ENDDO
      avg_mol_wt=mol_wt(1)*volmix(1)+mol_wt(2)*volmix(2)
      rayleigh_coefficient=Na/avg_mol_wt*constant*wl2i*wl2i*sum/       &
                           um2_to_m2                                     ! After this step, coefficient should be in m2/kg units

      END SUBROUTINE

!-----------------------------------------------------------------------
!+ This subroutine generates the aerosol optical properties for the    !
!     bands being used in the KDM code.  It takes data from a          !
!     scattering file and numerically integrates over the values       !
!     falling within the spectral band being used. For herritage       !
!     reasons, this code uses the term 'dust' a lot, when the more     !
!     general 'aerosol' (=dust plus ice plus whatever) is now approp.  !
!-----------------------------------------------------------------------
      SUBROUTINE process_aerosol(band_wn,n_band,wvl_in,omega_in,g_in,  &
                              cext_in,planck_temp,                     &
                              omega_dust,                              &
                              asymmetry_dust,kCext_dust,               &
                              l_simpledust,  &
                              n_aerosol_data,                          &
                              simpledust_omega, simpledust_asymmetry   )

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) ::                                n_band, &
                                                       n_aerosol_data

      REAL(KIND(0.d0)), INTENT(IN   ) ::                   planck_temp

      REAL(KIND(0.d0)), INTENT(IN   ), DIMENSION(n_band+1) ::          &
                                                               band_wn

      REAL(KIND(0.d0)), INTENT(IN   ), DIMENSION(n_aerosol_data) ::    &
                                                               wvl_in, &
                                                             omega_in, &
                                                                 g_in, &
                                                              cext_in

! the l_simpledust flag causes the visible optical properties for dust to be set to the same as
! those used in the WBM (Breigleb)
      real, intent(in), optional :: simpledust_omega, simpledust_asymmetry
      real :: simpledust_omega_value, simpledust_asymmetry_value
      LOGICAL, INTENT(IN   ) ::                          l_simpledust
      REAL(KIND(0.d0)), INTENT(  OUT), DIMENSION(n_band) ::            &
                                                           omega_dust, &
                                                       asymmetry_dust, &
                                                           kCext_dust

!c-mm  Local variables

      INTEGER ::                                                  i,k, &
                                                                 kmin, &
                                                                 kmax

      REAL(KIND(0.d0)) ::                                         tmp, &
                                                                 dtmp, &
                                                                 ntmp

      REAL(KIND(0.d0)), DIMENSION(n_band+1) ::                band_wl

  
#if ( WRF_MARS == 1 )
      if(l_simpledust .and. .not. (present(simpledust_omega) .and. present(simpledust_asymmetry))) then
          WRITE ( wrf_err_message ,*) ' KDM process_aerosol : got l_simpledust but missing simpledust values in function call. ',present(simpledust_omega), present(simpledust_asymmetry)
          CALL wrf_error_fatal ( wrf_err_message )
      endif

      if (present(simpledust_omega)) &
           simpledust_omega_value = simpledust_omega
      if (present(simpledust_asymmetry)) &
           simpledust_asymmetry_value = simpledust_asymmetry

IF (l_simpledust) WRITE(0,*) 'Setting up KDM vis to use fixed (ssa=', &
                              simpledust_omega_value,                 &
                              ', asym=',                              &
                              simpledust_asymmetry_value,             &
                              ') dust optical properties'
#endif

      DO i=1,n_band+1
         band_wl(i)=1.d4/band_wn(n_band+2-i)
      ENDDO

      kmin=1
      kmax=1
      DO i=1,n_band

        DO k=kmin,n_aerosol_data
           IF (wvl_in(k) > band_wl(i)) THEN
!c-mm              kmin=k-1
!c-mm  For KDM code, the lowest wavelength is 0.24 um, but for the ice scattering
!c-mm     parameters, the table goes down only to 0.25 um.  That means this loop
!c-mm     would stop with k=1, setting kmin=0, and crashing the code, since it
!c-mm     would be looking for wvl_in(0), which doesn't exist.  So we will just
!c-mm     force it to start at the first row of the table.
!c-mm  My question is, why hadn't this been a problem before?  Did we never use
!c-mm     the ice data before?  8-22-17
              kmin=MAX(k-1,1)
              EXIT
           ENDIF
        ENDDO

        DO k=kmin,n_aerosol_data
           IF (wvl_in(k) > band_wl(i+1)) THEN
              kmax=k
              EXIT
           elseif (k == n_aerosol_data) then
              kmax = k
           ENDIF
        ENDDO

        !c-mm  Approximate integrals with trapezoidal integration

        !c-mm  Do Cext first
        dtmp = 0.0D0
        DO k=kmin,kmax-1
           tmp=0.5D0*( planck(wvl_in(k),planck_temp) +                 &
                       planck(wvl_in(k+1),planck_temp) )*              &
                       (wvl_in(k+1)-wvl_in(k))
           dtmp = dtmp+tmp
        ENDDO
        ntmp=0.0D0
        DO k=kmin,kmax-1
           tmp=0.5D0*( cext_in(k)*planck(wvl_in(k),planck_temp) +      &
                       cext_in(k+1)*planck(wvl_in(k+1),planck_temp) )* &
                       (wvl_in(k+1)-wvl_in(k))
           ntmp=ntmp+tmp
        ENDDO
        kCext_dust(n_band-i+1)=ntmp/dtmp                                 ! Put the value back in its correct 'wavenumber' index

#if ( WRF_MARS == 1 )
IF ( l_simpledust) kCext_dust(n_band-i+1) = kCext_dust_ref     ! eliminates any scaling between bands
!IF ( l_simpledust) kCext_dust(n_band-i+1) = 0.     ! this effectively sets VIS (not IR) dust to zero
#endif

        !c-mm  Do omega next.  Denominator of omega is the same as
        !c-mm     numerator of Cext.
        dtmp=ntmp
        ntmp=0.D0
        DO k=kmin,kmax-1
           tmp=0.5D0*( omega_in(k)*planck(wvl_in(k),planck_temp)*      &
               cext_in(k)+omega_in(k+1)*                               &
               planck(wvl_in(k+1),planck_temp)*cext_in(k+1) )*         &
               (wvl_in(k+1)-wvl_in(k))
           ntmp=ntmp+tmp
        ENDDO
        omega_dust(n_band-i+1)=ntmp/dtmp

#if ( WRF_MARS == 1 )
IF ( l_simpledust) omega_dust(n_band-i+1) = simpledust_omega_value  ! got to love random hardwired values!
#endif

        !c-mm  Do asymmetry next.  Denominator of asymmetry is the 
        !c-mm     same as numerator of omega
        dtmp=ntmp
        ntmp=0.D0
        DO k=kmin,kmax-1
           tmp=0.5D0*( g_in(k)*omega_in(k)*                            &
               planck(wvl_in(k),planck_temp)*cext_in(k)+g_in(k+1)*     &
               omega_in(k+1)*planck(wvl_in(k+1),planck_temp)*          &
               cext_in(k+1) )*(wvl_in(k+1)-wvl_in(k))
           ntmp=ntmp+tmp
        ENDDO
        asymmetry_dust(n_band-i+1)=ntmp/dtmp

#if ( WRF_MARS == 1 )
IF ( l_simpledust) asymmetry_dust(n_band-i+1) = simpledust_asymmetry_value  ! got to love random hardwired values!
#endif

      ENDDO

      END SUBROUTINE process_aerosol
!-----------------------------------------------------------------------
!+ This function calculates the planck function for a given temperature!
!     and wavelength                                                   !
!-----------------------------------------------------------------------
      FUNCTION planck(wvl,T)

      IMPLICIT NONE

      REAL(KIND(0.d0)), PARAMETER ::             h_planck = 6.626D-34, &  ! Planck constant (J/s)
                                                  c_light = 2.998D+08, &  ! Speed of light (m/s)
                                                  k_boltz = 1.381D-23     ! Boltzmann constant (J/K)

      REAL(KIND(0.d0)), INTENT(IN   ) :: wvl, T

      REAL(KIND(0.d0)) ::                planck

      planck = (1.0D30*2.0D0*h_planck*c_light*c_light *wvl**(-5))/     &
               (DEXP((1.0D6*h_planck*c_light)/(k_boltz*T*wvl))-1.0D0)
      END FUNCTION PLANCK

#if ( WRF_MARS == 1 )
    subroutine read_multisize_miedata(mie_filename,aerosol_number)

#ifdef mpas
       use mpas_atmphys_constants, only : pi2
#endif
        character(len=*), intent(in) :: mie_filename
        integer, intent(in) :: aerosol_number
        character(len=256) :: dummy, dummy2
        integer :: nbands_solar_dust, nbands_ir_dust, nref
        real :: veff, reff
        integer :: i,j
        integer, parameter :: expected_nref=5
        real, dimension(expected_nref), parameter :: expected_ref_wn=(/0.670000, 9.300000, 12.121212, 11.862396, 21.598272/)
        real, dimension(expected_nref) :: got_ref_wn
        
        OPEN(UNIT=2011, FILE=TRIM(mie_filename))
        
        !header line
        read(2011,*) dummy
        !number of bins, bandsSol,bandsIR, NR
        read(2011,*) nbins_mie, nbands_solar_dust, nbands_ir_dust, nref
        if(nref /= expected_nref) then
          write(wrf_err_message,*) "unexpected number of reference bands in read_multisize_miedata,", &
           " expected=",expected_nref,", got=",nref
          call wrf_error_fatal(trim(wrf_err_message))
        endif
        !veff label
        read(2011,*) dummy
        !veff data
        read(2011,*) veff
        !reff header
        read(2011,*) dummy
        !reff data
        read(2011,*) reffstart,reffstop,reffstep
        !ir band label
        read(2011,*) dummy
        !ir band data
        read(2011,*) dummy
        !solar band label
        read(2011,*) dummy
        !solar band data
        read(2011,*) dummy
        !qref label
        read(2011,*) dummy
        !qref data
        read(2011,*) got_ref_wn
        do i=1,nref
          if(got_ref_wn(i) /= expected_ref_wn(i)) then
             write(wrf_err_message,*) "index ",i," of Qref reference wl read in=",got_ref_wn(i), &
              "but expected=",expected_ref_wn(i)
            call wrf_error_fatal(trim(wrf_err_message))
          endif
        enddo
        !allocate table
        if(.not. allocated(Qref)) allocate(Qref(nref))
        if (.not. allocated(aerosol_mie)) allocate(aerosol_mie(npd_aerosol_species, nbins_mie,nbands_solar_dust+nbands_ir_dust,3)) !3=qext, w, g
        if (.not. allocated(aerosol_Qref)) allocate(aerosol_Qref(npd_aerosol_species, nbins_mie, nref)) ! just Q->0.67E, 9.0A, 20.E
        !loop to read table data
        do i=1, nbins_mie
            read(2011,*) reff, Qref !header line = reff
            aerosol_Qref(aerosol_number, i,:) = Qref(:) * (0.5 * pi2 * reff * reff * 1e-12)
            do j=1, nbands_solar_dust+nbands_ir_dust
                read(2011,*) aerosol_mie(aerosol_number, i,j,:)
                !convert extinction efficiency to extinction coefficient
                aerosol_mie(aerosol_number, i,j,1) = aerosol_mie(aerosol_number, i,j,1) * (0.5 *  pi2 * reff * reff * 1e-12) !scale by r.r to get extinction coeff.
            end do
        end do
    
        close(2011)
    
    end subroutine read_multisize_miedata

    subroutine mie_findex(reff, indices) 
    !find the indices into the mie array based on reff
    
    real(kind(0.d0)), dimension(:,:), intent(in) :: reff
    real(kind(0.d0)), dimension(:,:),  intent(out) :: indices
    
    indices = 1+(reff-reffstart)/reffstep
    
    where(indices<1)
        indices=1
    elsewhere(indices>=(nbins_mie-1d-5))
        indices=nbins_mie-1d-5
    end where
    !that was too easy
    end subroutine mie_findex

    subroutine mie_interpolation(reff_index, i_band, isolir, aerosol_number, dataindex, data_out)
        real(kind(0.d0)), dimension(:,:),  intent(in) :: reff_index
        integer, intent(in) :: i_band, isolir, dataindex, aerosol_number
        real(kind(0.d0)), dimension(:,:),  intent(out) :: data_out
        
        integer band
        integer :: i, j, l
        real(kind(0.d0)) :: w 
        
        if (isolir == ip_solar) then
            band = i_band
        else
            band = i_band + n_band_solar
        endif

        do i=1, size(data_out,1)
            do j=1, size(data_out,2)
                l = floor(reff_index(i,j))
                w = 1. - (reff_index(i,j) - l)
                data_out(i,j) = w*aerosol_mie(aerosol_number, l,band,dataindex) &
                          + (1-w)*aerosol_mie(aerosol_number, l+1,band,dataindex)
            end do
        end do
        
    end subroutine mie_interpolation

    subroutine reference_extinction(reff, ext,aerosol_number, reference_band)
        !repeating code is bad, but this function is much simpler than the full
        !mie_findex -> mie_interpolation and the reason they are seperate is because
        !the optimization comes from only calling findex once for every 4 interpolations.
        !that optimization doesn't occur here

        !  - reference band 1 = 0.67 micron (red) vis band
        !                   2 = 9.3 micron (silicate) thermal IR band
        !                   3 = 12.12 micron (water ice?) thermal IR
        !                   4 = 11.86 micron (water ice?) thermal IR
        !                   5 = 21.6 micron (co2 ice?) thermal IR

        ! - aerosol_number 1 = dust
        !                  2 = water ice
        !                  3 = co2 ice

        real(kind(0.d0)), intent(in) :: reff
        real(kind(0.d0)), intent(out) :: ext
        integer, intent(in)::aerosol_number
        integer, intent(in)::reference_band
        real (kind(0.d0)) :: index, w
        integer :: l

        index=1+(reff-reffstart)/reffstep
        index=min(nbins_mie-1d-5,max(1d0,index))
        l = floor(index)
        w = 1. - (index-l)
        ext =    w*aerosol_Qref(aerosol_number,l,reference_band) + &
             (1-w)*aerosol_Qref(aerosol_number,l+1,reference_band)
end subroutine reference_extinction


subroutine reference_opacity(tau_3d, tau_2d, naerosol, reff, pf3d,    &
                                ids,ide, jds,jde, kds,kde,         &
                                ims,ime, jms,jme, kms,kme,         &
                                its,ite, jts,jte, kts,kte,         &
                                aerosol_number,                    &
                                reference_band)
#ifdef mpas
       use mpas_atmphys_constants, only : g => gravity
#endif

     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN  ) ::    &
                                                        naerosol, reff,pf3d
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT ) ::    &
                                                        tau_3d
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT ) ::    &
                                                        tau_2d
    integer, intent(in), optional ::reference_band
    integer, intent(in) :: aerosol_number
    integer :: reference_band_use

     real, dimension( ims:ime, kms:kme, jms:jme ) :: temp
     real(kind(0.d0)) :: ref_extinction
     integer :: i,j,k


    reference_band_use = 1
    if(present(reference_band)) reference_band_use = reference_band
    tau_3d(:,kte,:) = 0.0
    do j=jts,jte
      do i=its, ite
        do k=kts, kte
            call reference_extinction(1.0d0*reff(i,k,j),ref_extinction, aerosol_number, reference_band=reference_band_use)
            tau_3d(i,k,j) = naerosol(i,k,j)*ref_extinction*(pf3d(i,k,j)-pf3d(i,k+1,j))/g
!            write(*,*) "b",i, j,k, ndust(i,k,j), reff(i,k,j), ref_extinction, tau_3d(i,k,j),(pf3d(i,k,j)-pf3d(i,k+1,j))/g
        ENDDO
            !finish integrating opacity
        do k=kte-1,kts,-1
!          write(*,*) i,j,k, tau_3d(i,k,j), ndust(i,k,j),ref_extinction,(pf3d(i,k,j)-pf3d(i,k+1,j))/g
          tau_3d(i,k,j) = tau_3d(i,k,j) + tau_3d(i,k+1,j)
        enddo
      enddo
    enddo
    
    tau_2d(its:ite,jts:jte) = tau_3d(its:ite,kts,jts:jte)

   end SUBROUTINE reference_opacity

   subroutine aerosol_obs_scaling(tau_3d, tau_2d, ndust, qdust, reff, pf3d,    &
                                ids,ide, jds,jde, kds,kde,         &
                                ims,ime, jms,jme, kms,kme,         &
                                its,ite, jts,jte, kts,kte,         &
                                tracking_sols,dustincr,fixed_incr, &
                                reference_band, target_cdod,       &
                                coszen,day_only,                   &
                                lifted_reff,lifted_veff,obs_scaling_limit,     &
                                dust_size_distribution,                        &
                                o_lif1,o_lif2,o_sed1,o_sed2,scale_mode, &
                                actual_rad_dt,do_nudge_over_co2ice,co2ice)

!c-mm     use module_mp_morr_two_moment, only: average_mass_of_particle
#ifdef mpas
       use mpas_atmphys_constants, only : g => gravity
#endif

     USE module_state_description, ONLY: MONTA_BRUTE_OVERRIDE,           &
                                         MONTA_LIFT_ONLY,                &
                                         MONTA_LIFT_AND_FALL,            &
                                         MONTA_TEST_NO_ACTION,           &
                                         WANG_MDGM_NUDGE

     IMPLICIT NONE
!---------------------------------------------------------------------
     INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT  ) ::  &
                                                        ndust, qdust
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN  ) ::    &
                                                        reff,pf3d
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT ) ::    &
                                                        tau_3d
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT ) ::    &
                                                        tau_2d
    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) ::    &
                                                        target_cdod, coszen
    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) ::    &
                                        o_lif1,o_lif2,o_sed1,o_sed2
                                integer, intent(in) :: scale_mode

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) ::    &
                                                        co2ice

    real(kind(0.d0)) , intent(in) :: lifted_reff, lifted_veff
    real , intent(in) :: obs_scaling_limit, dustincr
    integer, intent(in), optional ::reference_band
    real, intent(in) :: actual_rad_dt, tracking_sols
    integer, intent(in) :: dust_size_distribution
    logical, intent(in) :: do_nudge_over_co2ice, fixed_incr, day_only

    integer :: reference_band_use
    real  :: obs_scaling_limit_rate, mass_needed, number_needed
    real, dimension( ims:ime, kms:kme, jms:jme ) :: temp
    real(kind(0.d0)) :: scale
    real(kind(0.d0)) :: lifted_reff_extinction
    real(kind(0.d0)), dimension(kms:kme) :: extinction, mass, weight
    real(kind(0.d0)) :: delta_cdod, column_mass, ref_extinction, re2
    real :: avg_mass_particle, tracking_rate, day_fac
    integer :: i,j,k
    logical :: do_inject

    do_inject = .false. ! do not change to .true. !!! this is an auto-switch used in the scale_mode case select

    obs_scaling_limit_rate = obs_scaling_limit/86400.
    reference_band_use = 1
    if(present(reference_band)) reference_band_use = reference_band
    tau_3d(:,kte,:) = 0.0
    tau_2d(:,:) = 0.

    ! scale_mode selection options:
    !
    ! 1 = brute force scaling of q and n and place lifted / settled dust in o_lif/o_sed arrays (lif and sed diagnostic)
    ! 2 = calculate a rate of lifting to be applied and placed in o_lif (lif prognostic), don't do any forced sedimentation (zero sed)
    ! 3 = calculate a rate of lifting to be applied and placed in o_lif (lif prognostic), immediately settle and place in sed (sed diagnostic)
    ! 4 = test the plumbing by passing through this routine but don't modify any injection / fall / dust arrays

    ! calculate model-predicted opacity in observation band
    !  - reference band 1 = 0.67 micron (red) vis band
    !                   2 = 9.3 micron (silicate) thermal IR band
    !                   3 = 12.12 micron (water ice?) thermal IR
    !                   4 = 11.86 micron (water ice?) thermal IR
    !                   5 = 21.6 micron (co2 ice?) thermal IR

    ! - aerosol_number 1 = dust
    !                  2 = water ice
    !                  3 = co2 ice

    do j=jts,jte
      do i=its, ite

        do k=kts, kte   ! here we calculate the model-actual total column optical depth from the 2-moment tracer fields:
            re2 = reff(i,k,j)  ! reff is in microns by default
            call reference_extinction(re2,ref_extinction, &
                    aerosol_number=1, reference_band=reference_band_use)
            mass(k) = (pf3d(i,k,j)-pf3d(i,k+1,j))/g
            extinction(k) = ndust(i,k,j)*ref_extinction*mass(k)
            tau_3d(i,k,j) = extinction(k)
        ENDDO
        !finish integrating opacity

        do k=kte,kts,-1
          tau_2d(i,j) = tau_2d(i,j) + tau_3d(i,k,j)  ! this is optical depth down to actual surface pressure (not 610 or 700Pa)
        enddo
!        write(0,*) i,j," tau in model=",tau_2d(i,j)," target cdod=",target_cdod(i,j)
      enddo
    enddo

    ! calculate the avg mass of a particle in the dust surface distribution - this is the mass defined such that
    !     number density [#/kg] * avg_mass_particle [kg] = mass mixing ratio [kg/kg]
    ! we're only going to need these if we're injecting, but we don't need to keep recalcuating them for
    ! every its->ite and jts->jte, so let's just calculate it now:
    call reference_extinction(lifted_reff*1.d6, lifted_reff_extinction, &  ! reff input needs to be in microns
                    aerosol_number=1, reference_band=reference_band_use)
    avg_mass_particle = average_mass_of_particle(dust_size_distribution,real(lifted_reff),real(lifted_veff))

    ! let's make sure to always initialize to zero (and do it in one place!)
    do j=jts,jte
    do i=its,ite

        o_lif1(i,j)=0.
        o_lif2(i,j)=0.
        o_sed1(i,j)=0.
        o_sed2(i,j)=0.

    enddo
    enddo

    select case (scale_mode)

      case(MONTA_BRUTE_OVERRIDE)  ! just brute force the model dust to agree with observations:

       do j=jts,jte
       do i=its,ite

        scale = 1.
        if (tau_2d(i,j) > 0) then
          scale = target_cdod(i,j)/tau_2d(i,j)
          scale = min(obs_scaling_limit, max(scale, 1./obs_scaling_limit))
          scale = 1.0 + (scale - 1.)*0.1
!          write(0,*) i,j," scale=",scale," model tau=",tau_2d(i,j)," target tau=",target_cdod(i,j)
        endif

!CEN: Note that the lines below attribute an increase in column dust to net lifting (lift>dep) and a decrease to net deposition (lift<dep)
!in the former case, the increase in dust is shown by o_sed=0 and o_lif>0, where o_lif represents the *net* dust that would have been lifted in this tstep
!in the latter case, the decrease in dust is shown by o_lif=0 and o_sed>0, where o_sed represents the *net* dust that would have been deposited in this tstep
!but even this interpretation ignores transport between columns in the atmosphere itself
!Hence for these two reasons, 'o_lif' and 'o_sed' need to be interpreted with great caution when brute force scaling is used!

        if(scale > 1.) then

          do k=kts,kte
            mass(k) = (pf3d(i,k,j)-pf3d(i,k+1,j))/g   ! "mass" is actually in kg/m2
            o_lif1(i,j)=o_lif1(i,j)+mass(k)*(scale-1.)*qdust(i,k,j)/actual_rad_dt
            o_lif2(i,j)=o_lif2(i,j)+mass(k)*(scale-1.)*ndust(i,k,j)/actual_rad_dt
          enddo

        else

          do k=kts,kte
            mass(k) = (pf3d(i,k,j)-pf3d(i,k+1,j))/g   ! "mass" is actually in kg/m2
            o_sed1(i,j)=o_sed1(i,j)+mass(k)*(scale-1.)*qdust(i,k,j)/actual_rad_dt
            o_sed2(i,j)=o_sed2(i,j)+mass(k)*(scale-1.)*ndust(i,k,j)/actual_rad_dt
          enddo

        endif

        ndust(i,kts:kte,j) = ndust(i,kts:kte,j)*scale
        qdust(i,kts:kte,j) = qdust(i,kts:kte,j)*scale

       enddo
       enddo

       ! --------------------

      case(MONTA_LIFT_ONLY)  ! let's set a rate for injection for use in PBL scheme, but don't do forced sedimentation
               !  - this causes dust to be injected if the model tau < observed tau, but it lets
               !    sedimentation work at it's own rate (i.e., it won't act to augment the dust fall-out)
               !    if model tau > observed tau.


       do_inject = .true.  ! do the lifting in the "do_inject" conditional (below)

      case(MONTA_LIFT_AND_FALL) ! let's set a rate for injection for use in PBL scheme, but brute force sedimentation
               !  - this causes dust to be injected if the model tau < observed tau, and forces
               !    excess dust to be taken out by augmented sedimentation.

       do_inject = .true. ! do the lifting in the "do_inject" conditional (below)

       do j=jts,jte
       do i=its,ite

       scale = 1.
       if (tau_2d(i,j) < 0.) then
         scale = target_cdod(i,j)/tau_2d(i,j)  ! we're only going to use "scale" to see if it's > or < 1.
       endif

       if(scale < 1.) then  ! we need to remove some dust

          scale = min(obs_scaling_limit, max(scale, 1./obs_scaling_limit))
          scale = 1.0 + (scale - 1.)*0.1

          do k=kts,kte
            mass(k) = (pf3d(i,k,j)-pf3d(i,k+1,j))/g   ! "mass" is actually in kg/m2
            o_sed1(i,j)=o_sed1(i,j)+mass(k)*(scale-1.)*qdust(i,k,j)/actual_rad_dt
            o_sed2(i,j)=o_sed2(i,j)+mass(k)*(scale-1.)*ndust(i,k,j)/actual_rad_dt
          enddo

          ndust(i,kts:kte,j) = ndust(i,kts:kte,j)*scale
          qdust(i,kts:kte,j) = qdust(i,kts:kte,j)*scale

       endif

       enddo
       enddo

      case(MONTA_TEST_NO_ACTION)  ! nada case - test to flow through scaling but without scaling

       do_inject = .false.  ! do the lifting in the "do_inject" conditional (below)

      case default

       call wrf_error_fatal ( "invalid scale_mode in obs_scaling" )

      end select

      !  -- the injection, below, is automatically triggered for some scale_mode cases, above

      if(do_inject) then  ! doing the injection calc here so as not to replicate code for options 2 and 3, above

       tracking_rate = tracking_sols*86400.  ! this is the time constant for tracking to observations

       do j=jts,jte
       do i=its,ite

        if ((co2ice(i,j).le.0.).or.(do_nudge_over_co2ice)) then

          scale = 1.
          if (tau_2d(i,j) > 0.) then
            scale = target_cdod(i,j)/tau_2d(i,j)  ! we're only going to use "scale" to see if it's > or < 1.
          else if(target_cdod(i,j) > 0.) then
            scale = 1.5  ! value doesn't matter, just set a > 1 so that we lift into clear atmosphere
          endif

          if(scale > 1.) then  ! if dust is needing to be lifted (i.e., model < obs)

            if(fixed_incr) then
     !                         [ kg/m^2/s]               tracking_rate is multiplied by as we divide by it to get o_lif1
              mass_needed   = dustincr * tracking_rate    ! dustincr is in kg/m2/s
              number_needed = mass_needed / avg_mass_particle                         ! in #/m^2/timestep
            else
              number_needed = (target_cdod(i,j) - tau_2d(i,j))/lifted_reff_extinction  ! in #/m^2/timestep
              mass_needed   = number_needed * avg_mass_particle                        ! in kg/m^2/timestep
            endif

            if(day_only) then  ! only do obs nudging the day when sun up past a certain amount
              if(coszen(i,j) > 0.25) then
                day_fac = coszen(i,j)
                number_needed = number_needed * day_fac
                mass_needed   = mass_needed   * day_fac
              else
                mass_needed = 0.
                number_needed = 0.
              endif
            endif

            o_lif1(i,j)=mass_needed/tracking_rate   ! kg/m^2/s
            o_lif2(i,j)=number_needed/tracking_rate ! #/m^2/s

          endif

        endif

       enddo
       enddo

      endif


  END SUBROUTINE aerosol_obs_scaling

  SUBROUTINE aerosol_obs_scaling_1particle(dust3d,        &
                               pplev, ddpp, rhod,         &
                               totcoldust,                &
                               fixed_incr, dustincr,      &
                               scale_mode,                &
                               tracking_sols, dt,         &
                               o_lif1, o_lif2,            &
                               o_sed1, o_sed2,            &
                               ids,ide, jds,jde, kds,kde, &
                               ims,ime, jms,jme, kms,kme, &
                               its,ite, jts,jte, kts,kte  )

   USE module_state_description, ONLY: MONTA_BRUTE_OVERRIDE,           &
                                       MONTA_LIFT_ONLY,                &
                                       MONTA_LIFT_AND_FALL,            &
                                       MONTA_TEST_NO_ACTION

   IMPLICIT NONE

   INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

   REAL, DIMENSION(ims:ime, kts:kte, jms:jme), INTENT(IN) :: &
                                                     dust3d, &
                                                      pplev

   REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) ::          &
                                                 totcoldust

   REAL, INTENT(IN) :: rhod, ddpp, dustincr, tracking_sols, dt

   LOGICAL, INTENT(IN) :: fixed_incr
   INTEGER, INTENT(IN) :: scale_mode

   REAL, DIMENSION(ims:ime, jms:jme), INTENT(OUT) ::         &
                                             o_lif1, o_lif2, &
                                             o_sed1, o_sed2

   REAL, PARAMETER :: dustmin = 1.e-8, &
                          qext = 3.04

   REAL, DIMENSION(ims:ime,jms:jme) :: modeldust, target_tau_od2d
   REAL :: aerosol, qextrhor, dustincr_local
   INTEGER :: i,j,k

   o_lif1(its:ite,jts:jte) = 0.
   o_lif2(its:ite,jts:jte) = 0.
   o_sed1(its:ite,jts:jte) = 0.
   o_sed2(its:ite,jts:jte) = 0.

   !Find total column dust present above surface in model currently
    qextrhor=0.75*qext/(rhod*(0.5*ddpp))
    DO j = jts, jte
      DO i = its,ite
        modeldust(i,j)=0.
        DO k = kts,kte
          aerosol = max(dustmin,qextrhor*dust3d(i,k,j)*(pplev(i,k,j)-pplev(i,k+1,j))/g)
          modeldust(i,j)=modeldust(i,j)+aerosol
        ENDDO
      ENDDO
    ENDDO

    SELECT CASE(scale_mode)

      CASE(MONTA_LIFT_ONLY)

        DO j = jts, jte    ! nudge lifting to get to observations
          DO i = its,ite

             target_tau_od2d(i,j) = totcoldust(i,j)*700./pplev(i,kts,j) ! tau_od2d are defined at 700mb

             if(fixed_incr) then
                     dustincr_local = dustincr   ! dustincr is in kg/m^2/s
             else
                     dustincr_local = max(0.,totcoldust(i,j)-modeldust(i,j)) ! tau difference from obs
                     dustincr_local = dustincr_local / (qextrhor*(pplev(i,kts,j)-pplev(i,kts+1,j))/g) ! convert to delta kg/m2/s
                     dustincr_local = (dt/(tracking_sols*86400.)) * dustincr_local  ! scale to take tracking_sols to get there
             endif

             if (totcoldust(i,j).gt.modeldust(i,j)) then
                o_lif1(i,j)=dustincr_local
             else
                o_lif1(i,j)=0.
             end if

          ENDDO
        ENDDO

      CASE(MONTA_TEST_NO_ACTION)
        o_lif1(its:ite,jts:jte) = 0.

      CASE DEFAULT
       call wrf_error_fatal("invalid scale_mode in obs_scaling for 1 particle")

      END SELECT

  END SUBROUTINE aerosol_obs_scaling_1particle

  function build_ko_options(use_multimie_dust, use_multimie_water, use_debug) result(return_value)
    logical, intent(in), optional :: use_multimie_dust, &
                                     use_multimie_water, &
                                     use_debug

    integer :: return_value

    return_value = 0
    if(present(use_debug)) then
      if(use_debug) return_value = return_value + 2**ko_debug
    endif

    if(present(use_multimie_dust)) then
      if(use_multimie_dust) return_value = return_value + 2**ko_multimie_dust
    endif

    if(present(use_multimie_water)) then
      if(use_multimie_water) return_value = return_value + 2**ko_multimie_water
    endif

      return
  end function build_ko_options
 
  subroutine debug_ko_options(kdm_options_value)
    integer, intent(in) :: kdm_options_value
    integer :: pos
    logical :: bool
    character(len=255), dimension(16) :: strings
    strings(1) = "ko: Debug Flag Enabled"
    strings(2) = "ko: MultiMie Dust enabled"
    strings(3) = "ko: MultiMie Water Ice enabled"
    strings(4) = "dummy"
    strings(5) = "dummy"
    strings(6) = "dummy"
    strings(7) = "dummy"
    strings(8) = "dummy"
    strings(9) = "dummy"
    strings(10) = "dummy"
    strings(11) = "dummy"
    strings(12) = "dummy"
    strings(13) = "dummy"
    strings(14) = "dummy"
    strings(15) = "dummy"
    strings(16) = "dummy"

    do pos=0,15
      bool = btest(kdm_options_value, pos)
      if (bool) then 
        write(0,*) pos, bool, trim(strings(pos+1))
      else
        write(0,*) pos, bool
      endif
    end do

    !if (btest(kdm_options_value,ko_debug)) write(0,*) "ko: Debug Flag Enabled"
    !if (btest(kdm_options_value,ko_multimie_dust)) write(0,*) "ko: MultiMie Dust enabled"
    !if (btest(kdm_options_value,ko_multimie_water)) write(0,*) "ko: MultiMie Water Ice enabled"
    
  end subroutine debug_ko_options

!c-mm----------------------------------------------------------------
!c-mm  This subroutine takes the reference Cext value for CO2 ice
!c-mm     and calcualtes a global map of CO2 cloud opacity, given
!c-mm     the amount of CO2 ice in each model grid/layer.
!c-mm  This is done uniquely for CO2 ice (versus dust or H2O ice)
!c-mm     because we don't have observations of opacity values in
!c-mm     hand and we need to calculate them from the scattering
!c-mm     parameters.
!c-mm----------------------------------------------------------------
   SUBROUTINE calc_ice_tau(qst_ice, cloud_array, pf3d,      &
                             kCext_ref,                     &
                             part_radius, part_density,     &
                             ids,ide, jds,jde, kds,kde,     &
                             ims,ime, jms,jme, kms,kme,     &
                             its,ite, jts,jte, kts,kte)

#ifdef mpas
   USE mpas_atmphys_constants, only : g => gravity
#else
   USE module_model_constants
#endif

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qst_ice
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pf3d


   REAL(KIND(0.d0)), INTENT(IN   ) :: kCext_ref

   REAL, INTENT(IN   ) ::            part_radius, &
                                     part_density

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) :: cloud_array

   !Local variables

   INTEGER :: i,j,k
         
   REAL, PARAMETER :: pi=3.14159265d0

   REAL :: Cext_ref

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) :: ice_mass

   !!! LJS
   INTEGER :: cloudlev, error
   LOGICAL :: addcloud

      ! Get extinction coefficient for water ice (iflag=1) or CO2 ice (iflag=2)
      !if (iflag == 1) Cext_ref = kCext_ref
      !if (iflag == 2) Cext_ref = kCext_ref/pi/(part_radius*1.e6)**2
      Cext_ref = kCext_ref/pi/(part_radius*1.e6)**2 ! new parameter file for ice formatted same as for CO2

      ! Set optical depth to zero
      cloud_array(its:ite,kts:kte+1,jts:min(jte,jde-1)) = 0.0

        ! Do normal ice opacity calculation
        DO j=jts,MIN(jte,jde-1)
          DO i=its,ite
            DO k=kte,kts,-1
              ice_mass(i,k,j) = qst_ice(i,k,j)*REAL((pf3d(i,k,j)-pf3d(i,k+1,j))/g,KIND(0.d0))
              cloud_array(i,k,j) = cloud_array(i,k+1,j) + 0.75*Cext_ref*ice_mass(i,k,j)/part_density/part_radius
            ENDDO
          ENDDO
        ENDDO

   END SUBROUTINE calc_ice_tau

   SUBROUTINE ice_tau_chicago(cloud_array, pf3d,             &
                              l_s, xlat, xlong, hrang,       & !!! LJS
                              ids,ide, jds,jde, kds,kde,     &
                              ims,ime, jms,jme, kms,kme,     &
                              its,ite, jts,jte, kts,kte)

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   REAL, INTENT(IN    )    ::  l_s                                    !!! LJS
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: xlat, xlong, hrang

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pf3d

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) :: cloud_array

   !Local variables

   INTEGER :: i,j,k
   REAL :: opac, loctime
   INTEGER :: cloudlev
   LOGICAL :: addcloud

      ! Set optical depth to zero
      cloud_array(its:ite,kts:kte+1,jts:jte) = 0.0

      ! Add layer of cloud. In each column, sum(cloud_array) = ice_tau value from file !!! LJS
        if ((l_s >= chi_ls1) .and. (l_s <= chi_ls2)) then
          DO j=jts,MIN(jte,jde-1)
            DO i=its,MIN(ite,ide-1)

              loctime = mod((12.+(24.*hrang(i,j)/pi2)),24.)  ! convert hour angle in rads to local time in hours
              if(loctime < 0.) loctime = loctime + 24.

              if ((xlat(i,j) >= chi_lat1) .and. (xlat(i,j) <= chi_lat2) .and.  &
                  (xlong(i,j) >= chi_lon1) .and. (xlong(i,j) <= chi_lon2))      &
              then

                addcloud = .false.
                if (((chi_dayornight == 1) .and. ((loctime > chi_tval1) .and. (loctime < chi_tval2)))) addcloud = .true.
                if (((chi_dayornight == 2) .and. ((loctime < chi_tval1) .or.  (loctime > chi_tval2)))) addcloud = .true.

                ! Find pressure level for base of clouds
                DO k=kte,kts,-1
                  if ((pf3d(i,k,j) >= chi_cloud_prs) .and. (pf3d(i,k+1,j) <= chi_cloud_prs)) then
                    cloudlev = k
                    exit
                  endif
                ENDDO

                ! Add in optical depth for number of specified levels above cloud base
                DO k=cloudlev+chi_ncloud-1,kts,-1
                  opac = 0.0
                  if ((k >= cloudlev) .and. addcloud) then
                    opac = chi_ice_tau*real((pf3d(i,k,j)-pf3d(i,k+1,j)),kind(0.d0))/  &
                                   real((pf3d(i,cloudlev,j)-pf3d(i,cloudlev+chi_ncloud,j)),kind(0.d0)) 
                  endif
                  cloud_array(i,k,j) = cloud_array(i,k+1,j) + opac
                ENDDO

              endif
            ENDDO
          ENDDO
        endif

!        DO j=jts,MIN(jte,jde-1)
!         DO i=its,MIN(ite,ide-1)
!          DO k=kts,kte

!            if(cloud_array(i,k,j) .lt. 0.) write(0,*) "negative chi cloud: ",i,j,k,cloud_array(i,k,j)
!            if(cloud_array(i,k,j) .gt. 5.) write(0,*) "thick chi cloud: ",i,j,k,cloud_array(i,k,j)

!          enddo
!          if (j==5 .and. i==37) then
!                  write(0,*) "--> cloud :",cloud_array(i,kts,j),pf3d(i,kts,j),chi_ice_tau,chi_cloud_prs
!                  write(0,*) "  ",chi_ncloud,chi_ls1,chi_ls2,chi_lat1,chi_lat2,chi_lon1,chi_lon2,chi_tval1
!                  write(0,*) "  ",chi_tval2, chi_dayornight, chi_ice_tau,loctime
!          endif

!         enddo
!        enddo

    END SUBROUTINE ice_tau_chicago

    SUBROUTINE init_ice_tau_chicago(chi_cloud_prs_in, chi_ncloud_in,          &
                                    chi_ls1_in, chi_ls2_in, chi_lat1_in,      &
                                    chi_lat2_in, chi_lon1_in, chi_lon2_in,    &
                                    chi_tval1_in, chi_tval2_in,               &
                                    chi_dayornight_in, chi_ice_tau_in         )

    IMPLICIT NONE

    real, intent(in) ::    chi_cloud_prs_in,                         &
                           chi_ls1_in, chi_ls2_in, chi_lat1_in,      &
                           chi_lat2_in, chi_lon1_in, chi_lon2_in,    &
                           chi_tval1_in, chi_tval2_in,               &
                           chi_ice_tau_in
    integer, intent(in) :: chi_ncloud_in, chi_dayornight_in

    chi_cloud_prs  = chi_cloud_prs_in
    chi_ls1        = chi_ls1_in
    chi_ls2        = chi_ls2_in
    chi_lat1       = chi_lat1_in
    chi_lat2       = chi_lat2_in
    chi_lon1       = chi_lon1_in
    chi_lon2       = chi_lon2_in
    chi_tval1      = chi_tval1_in
    chi_tval2      = chi_tval2_in
    chi_ice_tau    = chi_ice_tau_in

    chi_ncloud     = chi_ncloud_in
    chi_dayornight = chi_dayornight_in
    
    END SUBROUTINE
#endif

    SUBROUTINE kdm_debug(message,do_print)

      implicit none

      character(len=*), intent(in) :: message
      logical,          intent(in) :: do_print

!c-mm      call wrf_debug(500,trim(message))
      IF (do_print) write(0,*) trim(message)

    END SUBROUTINE kdm_debug

#if ( WRF_TITAN == 1 )
      subroutine get_tholin_tau_omega_asym(tholin_radius,tholin_numden,i_band,isolir, &
                                           n_profile,n_layer,npd_aerosol_species,     &
                                         aerosol_opacity,omega_tholin,asymmetry_tholin)

        use module_nrutils, only: linear_interpolate

        implicit none

        INTEGER, INTENT(IN) :: n_profile,n_layer,npd_aerosol_species, i_band, isolir

        REAL(KIND(0.d0)), DIMENSION( n_profile, n_layer ),  INTENT(IN   ) ::     &
                                                      tholin_numden, &  ! in #/m2 (NOT m3)
                                                      tholin_radius

        REAL(KIND(0.d0)), dimension(n_profile, n_layer),             &
                                            INTENT(  OUT) ::         &
                                                       omega_tholin, &
                                                   asymmetry_tholin
        REAL(KIND(0.d0)),                                            &
           DIMENSION(n_profile, 0:n_layer, npd_aerosol_species),     &
           INTENT(OUT   ) ::                        aerosol_opacity

        INTEGER :: i,j

        REAL(KIND(0.d0)) :: cross_section

        IF(npd_aerosol_species /= 1) THEN
          call wrf_error_fatal("Titan KDM: only setup for one radiative aerosol, but called for more.")
        endif

        DO i=1,n_profile
          DO j=1,n_layer

            if(isolir == ip_solar) then

               omega_tholin(i,j) = linear_interpolate(x0=tholin_radius(i,j),             &
                                                       x=dble(tholin_rad_bins),          &
                                                       y=tholin_vis_ssa(i_band,:),       &
                           out_of_range_use_nearest_edge=.true.                          )

           asymmetry_tholin(i,j) = linear_interpolate(x0=tholin_radius(i,j),             &
                                                       x=dble(tholin_rad_bins),          &
                                                       y=tholin_vis_asym(i_band,:),      &
                           out_of_range_use_nearest_edge=.true.                          )

                   cross_section = linear_interpolate(x0=tholin_radius(i,j),             &
                                                       x=dble(tholin_rad_bins),          &
                                                       y=tholin_vis_ext_xsect(i_band,:), &
                           out_of_range_use_nearest_edge=.true.                          )

!cross_section = cross_section *1.d0

            elseif(isolir == ip_infra_red) then

               omega_tholin(i,j) = linear_interpolate(x0=tholin_radius(i,j),             &
                                                       x=dble(tholin_rad_bins),          &
                                                       y=tholin_ir_ssa(i_band,:),        &
                           out_of_range_use_nearest_edge=.true.                          )

           asymmetry_tholin(i,j) = linear_interpolate(x0=tholin_radius(i,j),             &
                                                       x=dble(tholin_rad_bins),          &
                                                       y=tholin_ir_asym(i_band,:),       &
                           out_of_range_use_nearest_edge=.true.                          )

                   cross_section = linear_interpolate(x0=tholin_radius(i,j),             &
                                                       x=dble(tholin_rad_bins),          &
                                                       y=tholin_ir_ext_xsect(i_band,:),  &
                           out_of_range_use_nearest_edge=.true.                          )

            else

              call wrf_error_fatal("KDM Titan: band is neither solar nor ir")
                           
            endif
                          aerosol_opacity(i,j,1) = cross_section * tholin_numden(i,j)
          ENDDO
        ENDDO

      end subroutine get_tholin_tau_omega_asym
#endif

      SUBROUTINE do_nan_checking(kdm_mode_flag, rthratensw, rthratenlw, gsw, glw    &
#if ( WRF_MARS == 1 )
                          ,dust_array, cloud_array, co2_cloud_array, dust_array_ir  &
#endif
#if ( WRF_TITAN == 1 )
                          ,mck_haze_tau64                                           &
#endif
                          ,ims, ime, jms, jme, kms, kme,                            &
                           its, ite, jts, jte, kts, kte,                            &
                           swdown, swdowndir, swdownunscat                          )

      implicit none

      INTEGER, INTENT(IN   ) ::          ims, ime, jms, jme, kms, kme, &
                                         its, ite, jts, jte, kts, kte, &
                                         kdm_mode_flag

#if ( WRF_MARS == 1 )
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) ::   &
                                                           dust_array, &  ! dust_array is definitionally VIS (0.67micron) for MarsWRF
                                                          cloud_array, &
                                                      co2_cloud_array, &
                                                        dust_array_ir
#endif

#if ( WRF_TITAN == 1 )
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT(INOUT) ::   &
                                                       mck_haze_tau64
#endif

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) ::   &
                                                           rthratensw, &
                                                           rthratenlw

      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::            &
                                                                  gsw, &
                                                                  glw
      
      REAL, DIMENSION ( ims:ime, jms:jme), INTENT(INOUT), OPTIONAL ::  &
                                                               swdown, & ! sw down at the surface (i.e. not net)
                                                            swdowndir, & ! direct beam sw down at the surface
                                                         swdownunscat    ! unscattered beam sw down at the surface

      call test_for_bad_detailed(rthratensw,"kdm rthratensw", &  ! check both regardless of kdm_mode
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte, &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

      call test_for_bad_detailed(rthratenlw,"kdm rthratenlw", &
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte, &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

      call test_for_bad_detailed_2d(gsw,"kdm gsw",        &
                                 ims,ime,jms,jme,         &
                                 its,ite,jts,jte,         &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

      call test_for_bad_detailed_2d(glw,"kdm glw",        &
                                 ims,ime,jms,jme,         &
                                 its,ite,jts,jte,         &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

     if(present(swdown)) &
      call test_for_bad_detailed_2d(swdown,"kdm swdown",  &
                                 ims,ime,jms,jme,         &
                                 its,ite,jts,jte,         &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

     if(present(swdowndir)) &
      call test_for_bad_detailed_2d(swdowndir,"kdm swdowndir", &
                                 ims,ime,jms,jme,         &
                                 its,ite,jts,jte,         &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

     if(present(swdownunscat)) &
      call test_for_bad_detailed_2d(swdown,"kdm swdownuscar",  &
                                 ims,ime,jms,jme,         &
                                 its,ite,jts,jte,         &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

#if ( WRF_MARS == 1 )
      call test_for_bad_detailed(dust_array,"kdm dust_array", &
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte, &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

      call test_for_bad_detailed(cloud_array,"kdm cloud_array", &
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte, &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

      call test_for_bad_detailed(co2_cloud_array,"kdm co2_cloud_array", &
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte, &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )

      call test_for_bad_detailed(dust_array_ir,"kdm dust_array_ir", &
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte, &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )
#endif

#if ( WRF_TITAN == 1 )
      call test_for_bad_detailed(mck_haze_tau64,"kdm mck_haze_tau64", &
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte, &
                                 .true.,0.,               &
                                 kdm_mode=kdm_mode_flag   )
#endif

      END SUBROUTINE do_nan_checking

real function average_mass_of_particle(dust_distribution,reff,veff)

        use module_model_constants, ONLY: pi2

        implicit none

! in-n-out:
        integer, intent(in) :: dust_distribution
        real, intent(in)    :: reff, veff

! locals:
        real :: mu_dust, sigma_squared, local_const
        real :: rhod

!c-mm        rhod = DUST_CONS_RHO
        rhod = 2500.   !c-mm  Hardwiring this in...when's it ever going to change??!
        
! average_mass_of_particle is the average particle mass in the distribution. it is thus also the
! conversion function between the number/kg and mass/kg (mmr) moments
!
! q-array = avg_mass_particle * n-array
!
! see Lee et al Icarus 311 (2018) 23–34 Eqn. 9
!
! the code is called out here as this routine is needed outside of mp_morr_two_moment, but it
! replecates code in here - candidate for rolling more of this stuff up into subroutines like this

        select case (dust_distribution)

          case(0) !gamma
            mu_dust = (1.-3*veff)/veff
            sigma_squared = 0
            local_const =  ( (mu_dust + 2) * (mu_dust + 1) ) / ( (mu_dust + 3) * (mu_dust + 3) )
            average_mass_of_particle = ((4/3.) * 0.5*pi2*rhod*reff**3) * local_const

          case(1) !lognormal
            sigma_squared = alog(veff+1)
            local_const = exp(3*sigma_squared)
            average_mass_of_particle = ((4/3.)*0.5*pi2*rhod*reff**3) / local_const
        
          case default
            call wrf_error_fatal('error in average_mass_of_particle: dust_distribution not valid')

        end select
        return

end function average_mass_of_particle

subroutine numberdensity_from_opacity(opacity, reff, pf3d, naerosol, &
                                      ids,ide,jds,jde,kds,kde,       &
                                      ims,ime,jms,jme,kms,kme,       &
                                      its,ite,jts,jte,kts,kte,       &
                                      aerosol_number, reference_band)
     USE module_wrf_error
     implicit none

     INTEGER, INTENT(IN   ) :: ids,ide,jds,jde,kds,kde,              &
                               ims,ime,jms,jme,kms,kme,              &
                               its,ite,jts,jte,kts,kte

     REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::    &
                               opacity, reff, pf3d
     REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(  OUT) ::    &
                               naerosol

     INTEGER, INTENT(IN   ), OPTIONAL :: reference_band
     INTEGER, INTENT(IN   ) :: aerosol_number
     INTEGER :: reference_band_use

     REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: temp
     REAL(KIND(0.d0)) :: ref_extinction
     INTEGER :: i,j,k

     reference_band_use = 1
     IF(PRESENT(reference_band)) reference_band_use = reference_band

     DO j=jts,jte
        DO i=its,ite
           DO k=kts,kte
              IF (btest(kdm_options,ko_multimie_dust)) THEN
                 CALL reference_extinction(1.0d0*reff(i,k,j),        &
                      ref_extinction,aerosol_number,                 &
                      reference_band=reference_band_use)
              ELSE
                 ref_extinction=kCext_dust_ref
              ENDIF
              naerosol(i,k,j)=(opacity(i,k,j)-opacity(i,k+1,j))/(ref_extinction*(pf3d(i,k,j)-pf3d(i,k+1,j))/g)
           ENDDO
        ENDDO
     ENDDO

END subroutine numberdensity_from_opacity

END MODULE module_ra_kdm
