PROGRAM main
      USE module_radiation_driver
      USE module_ra_kdm
      USE module_planet_utilities, ONLY : get_julian
      USE namelist
        
      IMPLICIT NONE
      REAL, PARAMETER :: p2si = 1.027491252! 1.0274907
      REAL, PARAMETER :: degrad = 0.01745329
      REAL, PARAMETER :: rcp = 0.25
      REAL, PARAMETER :: p1000mb = 610.
      REAL, PARAMETER :: r_d = 189.
      REAL, PARAMETER :: grav = 3.72
      INTEGER::  ids,ide, jds,jde, kds,kde,        &
                 ims,ime, jms,jme, kms,kme,        &
                 its,ite, jts,jte, kts,kte          

      REAL :: l_s

      INTEGER :: i,j,k,kk,xx

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rthratensw, rthratenlw
      LOGICAL :: restart, l_double, l_simpledust, include_dust, include_water_ice, include_water_vap, include_co2_ice, polar
      REAL :: simpledust_omega, simpledust_asymmetry
      CHARACTER(LEN=80) :: k_coeff_file, dust_params_file, ice_params_file, co2ice_params_file

      REAL :: p_scale
      INTEGER :: kdm_options_incoming
      INTEGER :: t_sounding
      INTEGER :: dust_opacity_option, ice_opacity_option           
      CHARACTER(LEN=100)       :: rs_file_name

      REAL :: fys_factor, sunbody
      REAL :: obliquity, eccentricity, semimajoraxis, equinox_fraction, zero_date
      REAL :: l_s_input, time_of_day           
      REAL(KIND(0.d0)) :: ref_ext

      INTEGER :: p_qv, p_qc, p_qr, &
                 p_qi, p_qs, p_qg, &
                 p_first_scalar,   &
                 i_2stream_solar,  & ! Two-stream method for solar wavelengths
                 i_2stream_ir,     & ! Two-stream method for IR wavelengths
                 julday,           &
                 nlte_physics,     &
                 uv_physics,       &
                 ra_du_physics

      LOGICAL ::                            sw_correct
      LOGICAL ::                              diurnavg
      LOGICAL ::           ra_kdm_taulayer_calculation, &
                           ra_kdm_clearsky_calculation

      INTEGER ::                         kdm_mode_flag  ! 0=vis, 1=IR see kdm_do_vis and kdm_do_ir definitions
      INTEGER ::                   dust_reference_band, &
                                  cloud_reference_band
                                                 
      REAL(KIND(0D0)) ::                           dhr

      REAL ::                                   julian, &
                                                declin, &
                                                solcon, &
                                                   gmt
                                                                 
      REAL, DIMENSION( :, :), ALLOCATABLE :: &
            xlat,          &
            xlong,         &
            albedo,        & ! Surface albedo
            emiss,         & ! Surface emissivity--generally pegged at 1.0 for IR
            angslope,      &
            azmslope,      &
            co2ice,        &
            h2oice,        &
            sunfrac,       &
            tsk,           &
            psfc,          &
            ha,            &
            coszen,        &
            hrang,         &
            opt_depth,     &
            opt_depth_ice, &
            rnorm,         &
            tau

      REAL, DIMENSION( :, :, :), ALLOCATABLE :: &
           p3d,            &
           pf3d,           &
           pi3d,           &
           t3d,            &
           tf3d,           &
           qv3d,           &
           qc3d,           &
           qr3d,           &
           qi3d,           &
           qs3d,           &
           qg3d,           &
           rho_phy,        &
           dz8w,           &
           qst_co2ice

      REAL, DIMENSION( :, :, :), ALLOCATABLE :: &
           dust_array,      &
           cloud_array,     &
           co2_cloud_array, &
           dust_array_ir,   &
           opt,             &
           optdpth_array,   &
           dtau,            &
           dtau_ice
           
      REAL, DIMENSION(:), ALLOCATABLE :: &
           ra_kdm_gray_k_adhoc_ir, &
           ra_kdm_gray_k_adhoc_vis
           
      REAL, DIMENSION( :, :), ALLOCATABLE :: &
           gsw,   &
           glw,   &
           toasw, & ! Top of Atmosphere Short Wave (TOASW) flux for heat balance analysis.
           toalw    ! Top of Atmosphere Long Wave (TOALW) flux for heat balance analysis.
      
      REAL, DIMENSION( :, :, :), ALLOCATABLE :: &
           hr_g_vis, & ! Heating rate due to gas in solar wavelengths
           hr_a_vis, & ! Heating rate due to aerosol in solar wavelengths
           hr_g_ir,  & ! Heating rate due to gas in IR wavelengths
           hr_a_ir,  & ! Heating rate due to aerosol in IR wavelengths
           hr_g_uv     ! Heating rate due to gas in UV wavelengths

      REAL, DIMENSION( :, :), ALLOCATABLE :: &
           toa_sw_d, &
           toa_sw_u, &
           toa_lw_d, &
           toa_lw_u
           
      REAL, DIMENSION( :, :, :), ALLOCATABLE :: &
           flux_sw_d, &
           flux_sw_u, &
           flux_lw_d, &
           flux_lw_u, &
           hr_vis, &
           hr_ir

      REAL, DIMENSION ( :, :), ALLOCATABLE :: &
           swdown,       &  ! sw down at the surface (i.e. not net)
           swdowndir,    &  ! direct beam sw down at the surface
           swdownunscat

      REAL, DIMENSION ( :, :, :), ALLOCATABLE :: &
           ndust, reff_dust, nice, reff_ice

      REAL, DIMENSION(:, :), ALLOCATABLE :: cdod_scale

      INTEGER :: julian_int

      REAL ::     chi_ls1_in, &
                  chi_ls2_in, &
                 chi_lat1_in, &
                 chi_lat2_in, &
                 chi_lon1_in, &
                 chi_lon2_in, &
                chi_tval1_in, &
                chi_tval2_in, &
            chi_cloud_prs_in, &
              chi_ice_tau_in

      INTEGER :: chi_dayornight_in, &
                     chi_ncloud_in

      REAL ::     reg_t_ref, &
                  reg_t_iso, &
             reg_lapse_rate

      INTEGER ::    namelist_n_band_ir, &
                 namelist_n_band_solar

      LOGICAL ::        do_fake_layers, &
                     do_fake_layers_in, &
                 do_chicago_ice_layers, &
                      verbose_debug_in, &
                   nan_detail_check_in

      INTEGER ::    mp_physics, &
                 mp_physics_du

      LOGICAL :: f_qv, f_qi, f_qc, f_qr, f_qs, f_qg, f_qst_co2ice

      REAL ::     rhod_h2o, &
                  rhod_co2, &
             h2oice_radius, &
             co2ice_radius

      CALL namelist_init()
!c-mm Here, we define the dimensions of the model in 3-D (i,j,k)
!c-mm  For a reason beyond the scope of this comment, MarsWRF
!c-mm sets the 's'tart index to 1, and the 'e'nd index to 2, for
!c-mm a 1-D model.  Because of its parallel construction, there
!c-mm are other ways in which dimensions are calculated,
!c-mm (both the 'd' and 'm' variables).  Not going to get into it
!c-mm but this is done correctly.
      
      its = 1
      ite = 2
      jts = 1
      jte = 2
      kts = 1
      kte = n_layers_nml

      ids = 1
      ide = ite+1
      jds = 1
      jde = jte+1     
      kds = 1
      kde = kte+1

      ims = 1
      ime = ide
      jms = 1
      jme = jde
      kms = 1
      kme = kde

!c-mm VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
!c-mm This block contains parameters that are needed to get
!c-mm  the KDM code to operate, but are not needed for this
!c-mm  1-D version of the code.  I'm lumping these all here
!c-mm  but they DO NOT need to be touched.

      namelist_n_band_ir = 7.
      namelist_n_band_solar = 7.
      i_2stream_solar = 2
      i_2stream_ir = 12
      verbose_debug_in = .false.
      nan_detail_check_in = .false.
      restart = .FALSE.
      l_double = .TRUE.
      l_simpledust = .TRUE.
      simpledust_omega = 0.88
      simpledust_asymmetry = 0.55
      k_coeff_file      = "./co2_h2o_v4.dat"
      dust_params_file  = "./dust_gamma0.25v7.dat"
      ice_params_file   = "./water_gamma0.25v7.dat"
      co2ice_params_file = "./co2ice_10um_params.dat"
      rs_file_name = "" 
      kdm_options_incoming = 0 
      t_sounding = 2 
      polar = .TRUE.
      p_scale = 1
      ! co2ice = 0                !c-mm  XXXXXX  Keep here or move up?
      ! h2oice = 0
      fys_factor = 1
      reg_t_ref = 120.
      reg_t_iso = 120.
      reg_lapse_rate = 1.
      do_fake_layers = .false.
      f_qv = .false.
      f_qi = .false.
      f_qc = .false.
      f_qr = .false.
      f_qs = .false.
      f_qg = .false.
      f_qst_co2ice = .false.
      dhr = 0.
      do_fake_layers_in = do_fake_layers
      rhod_h2o = 1000.
      rhod_co2 = 970.
      h2oice_radius = 1.e-6
      co2ice_radius = 5.e-6
      mp_physics = mp_physics_nml
      mp_physics_du = mp_physics_du_nml
      ra_du_physics = ra_du_physics_nml
      p_first_scalar = 0 ! because we don't properly clear memory
      p_qv = 1
      p_qc = -1
      p_qr = -1
      p_qi = -1
      p_qs = -1
      p_qg = -1
      ra_kdm_taulayer_calculation = .FALSE.
      ra_kdm_clearsky_calculation = .FALSE.
      nlte_physics = 41
      sw_correct = .FALSE.
      diurnavg= .false.
      chi_ls1_in = 0.
      chi_ls2_in = 360.
      chi_lat1_in = -90.
      chi_lat2_in = 90.
      chi_lon1_in = -180.
      chi_lon2_in = 180.
      chi_tval1_in = 0.
      chi_tval2_in = 24.
      chi_dayornight_in = 1    !c-mm  Value of 1 indicates 'yes' for cloud if >tval1 and <tval2
      chi_cloud_prs_in = 0.1
      chi_ncloud_in = 1
      chi_ice_tau_in = 0.2
      do_chicago_ice_layers = .false.
      julday = 1
      ALLOCATE(xlong(ims:ime,jms:jme))
      xlong = 0.0
      ALLOCATE(sunfrac(ims:ime,jms:jme))
      sunfrac = 0.0
      ALLOCATE(ha(ims:ime,jms:jme))
      ha = 0.0
      ALLOCATE(coszen(ims:ime,jms:jme))
      coszen = 0.0
      ALLOCATE(hrang(ims:ime,jms:jme))
      hrang = 0.0
      ALLOCATE(qc3d(ims:ime, kms:kme, jms:jme))
      qc3d = 0
      ALLOCATE(qr3d(ims:ime, kms:kme, jms:jme))
      qr3d = 0
      ALLOCATE(qs3d(ims:ime, kms:kme, jms:jme))
      qs3d = 0
      ALLOCATE(qg3d(ims:ime, kms:kme, jms:jme))
      qg3d = 0
      ALLOCATE(qv3d(ims:ime, kms:kme, jms:jme))
      qv3d = 0
      ALLOCATE(qi3d(ims:ime, kms:kme, jms:jme))
      qi3d = 0
      ALLOCATE(angslope(ims:ime,jms:jme))
      angslope = 0.0
      ALLOCATE(azmslope(ims:ime,jms:jme))
      azmslope = 00
      ALLOCATE(ra_kdm_gray_k_adhoc_ir(namelist_n_band_ir))
      ra_kdm_gray_k_adhoc_ir=0.0
      ALLOCATE(ra_kdm_gray_k_adhoc_vis(namelist_n_band_solar))
      ra_kdm_gray_k_adhoc_vis=0.0
         
!c-mm ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!c-mm VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
!c-mm  In this block are all the variables we will need to run
!c-mm  the KDM code.  They are being properly allocated here.
      
      ALLOCATE(dust_array(ims:ime,kms:kme,jms:jme))
      dust_array = 0.
      ALLOCATE(dust_array_ir(ims:ime,kms:kme,jms:jme))
      dust_array_ir = 0.
      ALLOCATE(rthratensw(ims:ime,kms:kme,jms:jme))
      rthratensw = 0
      ALLOCATE(rthratenlw(ims:ime,kms:kme,jms:jme))
      rthratenlw = 0
      ALLOCATE(xlat(ims:ime,jms:jme))
      xlat = xlat_nml
      ALLOCATE(albedo(ims:ime,jms:jme))
      albedo = albedo_nml
      ALLOCATE(emiss(ims:ime,jms:jme))
      emiss = emiss_nml
      ALLOCATE(co2ice(ims:ime,jms:jme))
      co2ice = 0
      ALLOCATE(h2oice(ims:ime,jms:jme))
      h2oice = 0
      ALLOCATE(tsk(ims:ime,jms:jme))
      tsk = tsk_nml
      ALLOCATE(psfc(ims:ime,jms:jme))
      psfc=psfc_nml
      ALLOCATE(p3d(ims:ime, kms:kme, jms:jme))
      p3d = 0
      ALLOCATE(pf3d(ims:ime, kms:kme+1, jms:jme))
      pf3d = 0
      ALLOCATE(pi3d(ims:ime, kms:kme, jms:jme))
      pi3d = 0 
      ALLOCATE(t3d(ims:ime, kms:kme, jms:jme))
      t3d = 0
      ALLOCATE(tf3d(ims:ime, kms:kme+1, jms:jme))
      tf3d = 0
      ALLOCATE(rho_phy(ims:ime, kms:kme, jms:jme))
      rho_phy = 0
      ALLOCATE(dz8w(ims:ime, kms:kme, jms:jme))
      dz8w = 0
      ALLOCATE(qst_co2ice(ims:ime, kms:kme, jms:jme))
      qst_co2ice = 0
      ALLOCATE(cloud_array(ims:ime, kms:kme, jms:jme)) !water
      cloud_array = 0
      ALLOCATE(co2_cloud_array(ims:ime, kms:kme, jms:jme)) !co2
      co2_cloud_array = 0
      ALLOCATE(gsw(ims:ime,jms:jme))
      gsw = 0
      ALLOCATE(glw(ims:ime,jms:jme))
      glw = 0
      ALLOCATE(toasw(ims:ime,jms:jme))
      toasw = 0
      ALLOCATE(toalw(ims:ime,jms:jme))
      toalw = 0 
      ALLOCATE(hr_a_vis(ims:ime, kms:kme, jms:jme))
      hr_a_vis = 0
      ALLOCATE(hr_g_vis(ims:ime, kms:kme, jms:jme))
      hr_g_vis = 0
      ALLOCATE(hr_g_ir(ims:ime, kms:kme, jms:jme)) 
      hr_g_ir  = 0
      ALLOCATE(hr_a_ir(ims:ime, kms:kme, jms:jme)) 
      hr_a_ir = 0
      ALLOCATE(hr_g_uv(ims:ime, kms:kme, jms:jme)) 
      hr_g_uv = 0
      ALLOCATE(swdown(ims:ime, jms:jme))
      swdown = 0
      ALLOCATE(swdowndir(ims:ime, jms:jme))
      swdowndir = 0
      ALLOCATE(swdownunscat(ims:ime, jms:jme))
      swdownunscat = 0
      ALLOCATE(opt_depth(ims:ime, jms:jme))
      opt_depth = opt_depth_nml
      ALLOCATE(opt(ims:ime, kms:kme+1, jms:jme))
      opt = 0
      ALLOCATE(optdpth_array(ims:ime, kms:kme+1, jms:jme))
      optdpth_array = 0
      ALLOCATE(dtau(ims:ime, kms:kme, jms:jme))
      dtau = 0
      ALLOCATE(dtau_ice(ims:ime, kms:kme, jms:jme))
      dtau_ice = 0
      ALLOCATE(rnorm(ims:ime, jms:jme))
      rnorm = 0
      ALLOCATE(tau(ims:ime, jms:jme))
      tau = 0
      ALLOCATE(ndust(ims:ime, kms:kme, jms:jme))
      ndust = 0
      ALLOCATE(reff_dust(ims:ime, kms:kme, jms:jme))
      reff_dust = 0
      ALLOCATE(nice(ims:ime, kms:kme, jms:jme))
      nice = 0
      ALLOCATE(reff_ice(ims:ime, kms:kme, jms:jme))
      reff_ice = 0
      ALLOCATE(cdod_scale(ims:ime, jms:jme))
      cdod_scale = cdod_scale_nml
      ALLOCATE(toa_sw_d(ims:ime, jms:jme)) 
      toa_sw_d  = 0
      ALLOCATE(toa_sw_u(ims:ime, jms:jme)) 
      toa_sw_u = 0
      ALLOCATE(toa_lw_d(ims:ime, jms:jme)) 
      toa_lw_d = 0
      ALLOCATE(toa_lw_u(ims:ime, jms:jme)) 
      toa_lw_u = 0
      ALLOCATE(flux_sw_d(ims:ime, kms:kme, jms:jme)) 
      flux_sw_d  = 0
      ALLOCATE(flux_sw_u(ims:ime, kms:kme, jms:jme)) 
      flux_sw_u = 0
      ALLOCATE(flux_lw_d(ims:ime, kms:kme, jms:jme)) 
      flux_lw_d = 0
      ALLOCATE(flux_lw_u(ims:ime, kms:kme, jms:jme)) 
      flux_lw_u = 0
      ALLOCATE(hr_vis(ims:ime, kms:kme, jms:jme)) 
      hr_vis = 0
      ALLOCATE(hr_ir(ims:ime, kms:kme, jms:jme)) 
      hr_ir = 0
!c-mm ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
!c-mm  These are just dimensionless values--I don't need to copy them to new
!c-mm  variables, but I'm doing it just for consistency in nomenclature.

      include_dust = include_dust_nml
      include_water_ice = include_water_ice_nml
      include_water_vap = include_water_vap_nml
      include_co2_ice = include_co2_ice_nml
      dust_opacity_option = dust_opacity_option_nml
      ice_opacity_option = ice_opacity_option_nml
      
!c-mm  Here is where we set up the vertical grid using the values input from
!c-mm  the namelist.  Note that the grid indexes from the surface going up
!c-mm  and that the bottom layer is comprised of the surface values (psfc/tsk)
      
      pf3d(:,kts,:)=psfc
      tf3d(:,kts,:)=tsk
      DO k=kte,kts,-1
         pf3d(:,k+1,:)=pf3d_nml(k)
         tf3d(:,k+1,:)=tf3d_nml(k)
         t3d(:,k,:)=t3d_nml(k)
         p3d(:,k,:)=p3d_nml(k)
         rho_phy(:,k,:)=p3d(:,k,:)/r_d/t3d(:,k,:)
         pi3d(:,k,:)=(p3d(:,k,:)/p1000mb)**rcp
         dz8w(:,k,:)=(r_d*t3d(:,k,:)/grav)*log(pf3d(:,k+1,:)/pf3d(:,k,:))
         ndust(:,k,:)=ndust_nml(k)
         reff_dust(:,k,:)=reff_dust_nml(k)
         nice(:,k,:)=nice_nml(k)
         reff_ice(:,k,:)=reff_ice_nml(k)
      ENDDO
write(0,*) 'before init'
!c-mm  The kdminit subroutine initializes the KDM code
      CALL kdminit(restart,k_coeff_file,l_double,t_sounding,polar,do_fake_layers,           &
                   rs_file_name,mp_physics,namelist_n_band_ir,namelist_n_band_solar,        &
                   dust_params_file,l_simpledust,simpledust_omega,simpledust_asymmetry,     &
                   ice_params_file,co2ice_params_file,include_dust,include_water_ice,       &
                   include_water_vap,include_co2_ice,mp_physics_du,ra_du_physics,           &
                   do_chicago_ice_layers,chi_cloud_prs_in,chi_ncloud_in,chi_ls1_in,         &
                   chi_ls2_in,chi_lat1_in,chi_lat2_in,chi_lon1_in,chi_lon2_in,chi_tval1_in, &
                   chi_tval2_in,chi_dayornight_in,chi_ice_tau_in,reg_t_ref,reg_t_iso,       &
                   reg_lapse_rate,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,          &
                   its,ite,jts,jte,kts,kte                                                  &
                   )
write(0,*) 'after init'
!c-mm  To facilitate use of this code with more intuitive inputs, I prefer to use
!c-mm  L_s rather than sol-of-year (as carried by 'julian').  Therefore, the user
!c-mm  will input 'l_s_input' and the function get_julian will convert that to
!c-mm  the corresponding value of julian (sol-of-year).  This is needed by radconst.
!c-mm  Also, the user can input time of day which will be converted directly
!c-mm  into hour angle (ha or ha1), which is needed to calculate sza.
!c-mm  The resulting output of radconst will be, among other things, coszen, which
!c-mm  is needed by the KDM code.
      obliquity = obliquity_nml
      eccentricity = eccentricity_nml
      semimajoraxis = semimajoraxis_nml
      equinox_fraction = equinox_fraction_nml
      zero_date = zero_date_nml
      l_s_input = l_s_input_nml
      time_of_day = time_of_day_nml
      julian=get_julian(l_s_input,eccentricity,equinox_fraction,zero_date)

!c-mm  Right here is where you might want to put in some loop statements to generate
!c-mm  output for a range of Ls values, times of day, latitudes, etc.  Enddo
!c-mm  would go after the second call to planetary_kdm, near the bottom.  These would
!c-mm  overwrite what is in the namelist.

!c-mm    DO l_s_input=1,360
!c-mm       DO time_of_day=0,24

!c-mm  VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
!c-mm  Radconst spits out a number of parameters about the radiation field (declination, solar constant
!c-mm  zenith angle, etc.)  These are used by the KDM code to inform it of where the Sun is in the sky.
write(0,*) 'before radconst'      
      CALL radconst(0.0,  declin, solcon, julian,                         &
                    time_of_day, degrad, 0.0,                             &
                    obliquity, eccentricity, semimajoraxis,               &
                    equinox_fraction, zero_date, fys_factor,              &
                    xlat, ha, p2si, coszen, sunfrac, l_s, sunbody,        &
                    ids, ide, jds, jde, kds, kde,                         &
                    ims, ime, jms, jme, kms, kme,                         &
                    its, ite, jts, jte, kts, kte                          &
                    )
write(0,*) 'after radconst'                    
!c-mm      write(0,*) 'coszen',l_s,get_julian(l_s_input,eccentricity,equinox_fraction,zero_date),coszen(1,1),time_of_day
!c-mm                    ENDDO
!c-mm ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!c-mm VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
!c-mm This block of code calculates delta(tau) across each layer (as dtau).
!c-mm  This will allow the user to input a column opacity (as opt_depth)
!c-mm  and then let the code back out what ndust/nice should be, rather
!c-mm  than having the user have to put ndust/nice in the namelist, which
!c-mm  is awkward.

!c-mm There are two options, plus an implicit third option.  Option #1
!c-mm  allows you to input a fixed column opacity, using a Conrath
!c-mm  profile (hardwired with nu=0.01).  The code will then
!c-mm  determine the breakdown of opacities in each layer (dtau) and then
!c-mm  calculate the number of particles in that layer to produce
!c-mm  that dtau.  Option #2 allows you to prescribe dtau for each
!c-mm  layer directly.  This is perhaps more useful for water ice
!c-mm  which is often restricted to individual layers rather than
!c-mm  spread through the column.  The code will then convert dtau
!c-mm  into a number of particles.  The third option lets you input
!c-mm  ndust/nice directly.  It will default to this if option 1 or
!c-mm  2 is not selected.
!c-mm       Option #1:  Assigned total column opacity
!c-mm       Option #2:  Assigned layer opacities
!c-mm       Option #3:  Assigned layer number density
        
      IF (include_dust) THEN
         IF (dust_opacity_option.EQ.1) THEN
            opt(:,kte+1,:)=0.0
            DO k=kte,kts,-1
               tau=exp(0.01*(1.0-(pf3d(:,kts,:)/p3d(:,k,:))**(70./55.)))  !c-mm  Using Conrath profile
               opt(:,k,:)=opt(:,k+1,:)+tau*(pf3d(:,k,:)-pf3d(:,k+1,:))
            ENDDO
            rnorm=opt_depth/opt(:,kts,:)
            DO k=kts,kte+1
               optdpth_array(:,k,:)=opt(:,k,:)*rnorm
            ENDDO
            DO k=kts,kte
               dtau(:,k,:)=optdpth_array(:,k,:)-optdpth_array(:,k+1,:)   !c-mm  Calculates layer opacity
            ENDDO
         ENDIF
         IF (dust_opacity_option.EQ.2) THEN
            DO k=kts,kte
               dtau(:,k,:) = dtau_dust_nml(k)   !c-mm  Reads in layer opacities from namelist
            ENDDO
         ENDIF
         IF ((dust_opacity_option.EQ.1).OR.(dust_opacity_option.EQ.2)) THEN
            DO k=kts,kte
               !c-mm  Tau is defined in a reference band, so we need to scale accordingly.
               CALL reference_extinction(1.0d0*reff_dust(1,k,1),ref_ext,aerosol_number=1,reference_band=1)
               ndust(:,k,:)=dtau(:,k,:)/(ref_ext*(pf3d(:,k,:)-pf3d(:,k+1,:))/grav)
            ENDDO
         ENDIF
      ENDIF

      !c-mm  Now do the ice.  Same as dust
      IF (include_water_ice) THEN
         IF (ice_opacity_option.EQ.1) THEN
            opt(:,kte+1,:)=0.0
            DO k=kte,kts,-1
               tau=exp(0.01*(1.0-(pf3d(:,kts,:)/p3d(:,k,:))**(70./55.)))
               opt(:,k,:)=opt(:,k+1,:)+tau*(pf3d(:,k,:)-pf3d(:,k+1,:))
            ENDDO
            rnorm=opt_depth_ice/opt(:,kts,:)
            DO k=kts,kte+1
               optdpth_array(:,k,:)=opt(:,k,:)*rnorm
            ENDDO
            DO k=kts,kte
               dtau_ice(:,k,:)=optdpth_array(:,k,:)-optdpth_array(:,k+1,:)
            ENDDO
         ENDIF
         IF (ice_opacity_option.EQ.2) THEN
            DO k=kts,kte
               dtau_ice(:,k,:) = dtau_ice_nml(k)
            ENDDO
         ENDIF
         IF ((ice_opacity_option.EQ.1).OR.(ice_opacity_option.EQ.2)) THEN
            DO k=kts,kte
               !c-mm  Different reference band and aerosol number for ice.
               CALL reference_extinction(1.0d0*reff_ice(1,k,1),ref_ext,aerosol_number=2,reference_band=3)
               nice(:,k,:)=dtau_ice(:,k,:)/(ref_ext*(pf3d(:,k,:)-pf3d(:,k+1,:))/grav)
            ENDDO
         ENDIF
      ENDIF
      
      !c-mm  Print out the aerosol profile
      write(0,*) '             Aerosol Profiles'
      write(0,*) '        Layer    Dust           Ice'
      DO k=kts,kte
         write(0,*) kte+1-k,ndust(1,k,1),nice(1,k,1)
      ENDDO
      write(0,*)

!c-mm ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
write(0,*) 'before sw'
      kdm_mode_flag = 0   !c-mm  Start with the shortwave
        
      CALL planetary_kdm(rthratensw,rthratenlw,gsw,glw,toa_sw_u,toa_sw_d,                         &
                         toa_lw_u,toa_lw_d,flux_sw_u,flux_sw_d,flux_lw_u,flux_lw_d,xlat,          &
                         xlong,albedo,emiss,ha,coszen,angslope,azmslope,dhr,diurnavg,             &
                         rho_phy,t3d,tf3d,tsk,qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,p3d,pf3d,             &
                         pi3d,dz8w,julday,hr_g_vis,hr_a_vis,hr_vis,hr_g_ir,hr_a_ir,hr_ir,         &
                         sunfrac,julian,declin,solcon,l_s,hrang,f_qv,                             &
                         f_qc,f_qr,f_qi,f_qs,f_qg,i_2stream_solar,i_2stream_ir,                   &
                         nlte_physics,sw_correct,ra_du_physics,do_fake_layers_in,                 &
                         kdm_mode_flag,t_sounding,polar,p_scale,verbose_debug_in,                 &
                         nan_detail_check_in,dust_array,cloud_array,co2_cloud_array,              &
                         dust_array_ir,do_chicago_ice_layers,f_qst_co2ice,include_dust,           &
                         include_water_ice,include_water_vap,include_co2_ice,rhod_h2o,            &
                         rhod_co2,h2oice_radius,co2ice_radius,qst_co2ice,reff_dust,ndust,         &
                         reff_ice,nice,cdod_scale,ra_kdm_gray_k_adhoc_ir,ra_kdm_gray_k_adhoc_vis, &
                         ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte, &
                         ra_kdm_taulayer_calculation,ra_kdm_clearsky_calculation,                 &
                         swdown,swdowndir,swdownunscat                                            &
                         )
write(0,*) 'after sw'        
      write(0,*) 'TOA SW Down [W/m2]  TOA SW Up [W/m2]'
      write(0,*) toa_sw_d(1,1),toa_sw_u(1,1)
      write(0,*)
      write(0,'(a73)') 'Pressure [Pa]  Flux SW Down [W/m2]  Flux SW Up [W/m2]  Integrated Opacity'
      write(0,'(2x,f4.0,13x,f5.1,16x,f5.1,3x,a31,i3)') pf3d(1,kte+1,1),flux_sw_d(1,kte+1,1),flux_sw_u(1,kte+1,1),'------------------------- Level',kte
      DO k=kte,kts,-1
         write(0,'(55x,g12.5,a9,i3)') dust_array(1,k,1),'    LAYER',k
         write(0,'(2x,f4.0,13x,f5.1,16x,f5.1,3x,a31,i3)') pf3d(1,k,1),flux_sw_d(1,k,1),flux_sw_u(1,k,1),'------------------------- Level',k-1
      ENDDO
      write(0,*)
      write(0,*) 'GSW [W/m2]'
      write(0,*) gsw(1,1)
      write(0,*)
      write(0,*)
write(0,*) 'before lw'           
      kdm_mode_flag=1      !c-mm  Now do the IR

      CALL planetary_kdm(rthratensw,rthratenlw,gsw,glw,toa_sw_u,toa_sw_d,                         &
                         toa_lw_u,toa_lw_d,flux_sw_u,flux_sw_d,flux_lw_u,flux_lw_d,xlat,          &
                         xlong,albedo,emiss,ha,coszen,angslope,azmslope,dhr,diurnavg,             &
                         rho_phy,t3d,tf3d,tsk,qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,p3d,pf3d,             &
                         pi3d,dz8w,julday,hr_g_vis,hr_a_vis,hr_vis,hr_g_ir,hr_a_ir,hr_ir,         &
                         sunfrac,julian,declin,solcon,l_s,hrang,f_qv,                             &
                         f_qc,f_qr,f_qi,f_qs,f_qg,i_2stream_solar,i_2stream_ir,                   &
                         nlte_physics,sw_correct,ra_du_physics,do_fake_layers_in,                 &
                         kdm_mode_flag,t_sounding,polar,p_scale,verbose_debug_in,                 &
                         nan_detail_check_in,dust_array,cloud_array,co2_cloud_array,              &
                         dust_array_ir,do_chicago_ice_layers,f_qst_co2ice,include_dust,           &
                         include_water_ice,include_water_vap,include_co2_ice,rhod_h2o,            &
                         rhod_co2,h2oice_radius,co2ice_radius,qst_co2ice,reff_dust,ndust,         &
                         reff_ice,nice,cdod_scale,ra_kdm_gray_k_adhoc_ir,ra_kdm_gray_k_adhoc_vis, &
                         ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte, &
                         ra_kdm_taulayer_calculation,ra_kdm_clearsky_calculation,                 &
                         swdown,swdowndir,swdownunscat                                             &
                         )
write(0,*) 'after lw'
      write(0,*) 'TOA LW Down [W/m2]  TOA LW Up [W/m2]'
      write(0,*) toa_lw_d(1,1),toa_lw_u(1,1)
      write(0,*)
      write(0,'(a73)') 'Pressure [Pa]  Flux LW Down [W/m2]  Flux LW Up [W/m2]  Integrated Opacity'
      write(0,'(2x,f4.0,13x,f5.1,16x,f5.1,3x,a31,i3)') pf3d(1,kte+1,1),flux_lw_d(1,kte+1,1),flux_lw_u(1,kte+1,1),'------------------------- Level',kte
      DO k=kte,kts,-1
         write(0,'(55x,g12.5,a9,i3)') dust_array_ir(1,k,1),'    LAYER',k
         write(0,'(2x,f4.0,13x,f5.1,16x,f5.1,3x,a31,i3)') pf3d(1,k,1),flux_lw_d(1,k,1),flux_lw_u(1,k,1),'------------------------- Level',k-1
      ENDDO
      write(0,*)
      write(0,*) 'GLW [W/m2]'
      write(0,*) glw(1,1)

!c-mm      ENDDO     !c-mm  End of the loops over time/season/location/etc.
!c-mm   ENDDO
      

END PROGRAM main
