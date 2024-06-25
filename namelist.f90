module namelist
  IMPLICIT NONE

  INTEGER ::                            n_layers_nml
  
  LOGICAL, DIMENSION(7) ::           use_band_ir_nml, &
                                    use_band_vis_nml

  LOGICAL ::               do_chicago_ice_layers_nml

  REAL*8 ::                               albedo_nml, &
                                            xlat_nml, &
                                       l_s_input_nml, &
                                     time_of_day_nml, &
                                           emiss_nml, &
                                              ha_nml, &
                                          coszen_nml, &
                                             tsk_nml, &
                                            psfc_nml, &
                                       obliquity_nml, &
                                    eccentricity_nml, &
                                       zero_date_nml, &
                                   semimajoraxis_nml, &
                                equinox_fraction_nml

  REAL*8 ::                     chi_cloud_prs_in_nml, &
                                   chi_ncloud_in_nml, &
                                  chi_ice_tau_in_nml

  REAL*8, DIMENSION(:), ALLOCATABLE ::       t3d_nml, &
                                            tf3d_nml, &
                                             p3d_nml, &
                                            pf3d_nml, &
                                            pi3d_nml, &
                                            dz8w_nml, &
                                       dtau_dust_nml, &
                                           ndust_nml, &
                                       reff_dust_nml, &
                                        dtau_ice_nml, &
                                            nice_nml, &
                                        reff_ice_nml

  LOGICAL ::                        include_dust_nml, &
                               include_water_ice_nml, &
                               include_water_vap_nml, &
                                 include_co2_ice_nml

  REAL*8 ::                        h2oice_radius_nml, &
                                   co2ice_radius_nml, &
                                      cdod_scale_nml, &
                                       opt_depth_nml, &
                                   opt_depth_ice_nml

  INTEGER ::                          mp_physics_nml, &
                                   mp_physics_du_nml, &
                                   ra_du_physics_nml, &
                             dust_opacity_option_nml, &
                              ice_opacity_option_nml

  PUBLIC namelist_init, get_lun

CONTAINS

  SUBROUTINE namelist_init()

    character*256 :: filename
    integer :: unit
    logical :: exs
    integer :: i, &
         count
    namelist /layer_data_nml/ n_layers_nml

    namelist /band_nml/ use_band_ir_nml, use_band_vis_nml

    namelist /data_nml/ albedo_nml,xlat_nml,l_s_input_nml,time_of_day_nml, &
         emiss_nml,ha_nml,coszen_nml, &
         t3d_nml,tf3d_nml,tsk_nml,psfc_nml,p3d_nml, &
         pf3d_nml,pi3d_nml,dz8w_nml,   &
         include_dust_nml,include_water_ice_nml,             &
         include_water_vap_nml,include_co2_ice_nml,  &
         opt_depth_nml,opt_depth_ice_nml,dust_opacity_option_nml, &
         ice_opacity_option_nml,dtau_dust_nml,ndust_nml,reff_dust_nml, &
         dtau_ice_nml,nice_nml,reff_ice_nml,cdod_scale_nml


    namelist /orbital_nml/ obliquity_nml,eccentricity_nml,zero_date_nml, &
         semimajoraxis_nml,equinox_fraction_nml

    namelist /physics_params_nml/ mp_physics_nml,mp_physics_du_nml,ra_du_physics_nml

    !c-mm Read in layer_data_nml values
    filename="input.nml"
    INQUIRE(file=TRIM(filename), exist=exs)
    IF (exs) THEN
       CALL get_lun(unit)
       OPEN(unit,file=TRIM(filename))
       READ(unit, nml=layer_data_nml)
       CLOSE(unit)
    ENDIF

    ALLOCATE(p3d_nml(n_layers_nml))
    ALLOCATE(pf3d_nml(n_layers_nml+1))
    ALLOCATE(t3d_nml(n_layers_nml))
    ALLOCATE(tf3d_nml(n_layers_nml))
    ALLOCATE(dtau_dust_nml(n_layers_nml))
    ALLOCATE(dtau_ice_nml(n_layers_nml))
    ALLOCATE(ndust_nml(n_layers_nml))
    ALLOCATE(reff_dust_nml(n_layers_nml))
    ALLOCATE(nice_nml(n_layers_nml))
    ALLOCATE(reff_ice_nml(n_layers_nml))
    
    !c-mm Read in band_nml values
    filename="input.nml"
    INQUIRE(file=TRIM(filename), exist=exs)
    IF (exs) THEN
       CALL get_lun(unit)
       OPEN(unit,file=TRIM(filename))
       READ(unit, nml=band_nml)
       CLOSE(unit)
    ENDIF

    !c-mm Read in data_nml values
    filename="input.nml"
    INQUIRE(file=TRIM(filename), exist=exs)
    IF (exs) THEN
       CALL get_lun(unit)
       OPEN(unit,file=TRIM(filename))
       READ(unit, nml=data_nml)
       CLOSE(unit)
    ENDIF

    !c-mm Read in orbital_nml values
    filename="input.nml"
    INQUIRE(file=TRIM(filename), exist=exs)
    IF (exs) THEN
       CALL get_lun(unit)
       OPEN(unit,file=TRIM(filename))
       READ(unit, nml=orbital_nml)
       CLOSE(unit)
    ENDIF

    !c-mm Read in physics_params_nml values
    filename="input.nml"
    INQUIRE(file=TRIM(filename), exist=exs)
    IF (exs) THEN
       CALL get_lun(unit)
       OPEN(unit,file=TRIM(filename))
       READ(unit, nml=physics_params_nml)
       CLOSE(unit)
    ENDIF

  END SUBROUTINE namelist_init

  SUBROUTINE get_lun(unit)
    IMPLICIT NONE

    INTEGER :: lun
    INTEGER, INTENT(OUT) :: unit
    LOGICAL :: exs, opn
    DO lun=10,99
       INQUIRE(unit=lun, exist=exs, opened=opn)
       IF (exs.AND..NOT.opn) THEN
          unit=lun
          RETURN
       ENDIF
    ENDDO
    write(*,*) "There are no free Fortran logical units available"
    RETURN
    END SUBROUTINE get_lun

  END MODULE namelist
  
