!IDEAL:MODULE_LAYER:INITIALIZATION

!-----------------------------------------------------------------------------
!
! This module contains 4 function:
! Z2P : takes a height as input and returns a pressure
! Z2T : takes a height as input and returns a temperature
! P2Z : takes a pressure as input and returns a height
! P2T : takes a pressure as input and returns a temperature
!
! All four funtions need the "t_sounding" variable, which describes how
! T is related to P or Z ... either there is an analytic function (isothermal,
! constant lapse rate, etc.) or there is some hardwired input values.
! Each subroutine breaks down into a "CASE" statement based on the value
! of "t_sounding"
!
! Also, the functions require a few extra parameters, like T_ref, T_iso,
! lapse_rate, g, r_d, p0, and p_scale, since theses variables are necessary
! for the math in each of the funcations.
!
! All four functions work off the same principle: Calculate the return value
! by analytic function when possible, otherwise integrate the hydrostatic
! equation to get the third of the missing trio of variables (P,T,Z).  If
! using the hydrostatic equation, generate a vector of hardwired values, and
! then interpolate to find the result for the input value.  Sometimes, based
! on what is hardwired, it is easier to use the other functions within one
! of the functions, e.g., if T(p) is hardwired, z -> T is done as z -> P and
! then P -> T.
!
! If any new hardwired profiles are put in, copy and pasting from the examples
! should be sufficient (with the exception of the eta->pressure level kludge)
! and having the new hardwired profile entered in the common area.
!
!-----------------------------------------------------------------------------

MODULE module_setup_tpz

   USE module_wrf_error          ! frame
   USE module_init_utilities     ! dyn_em
   USE module_nrutils            ! share
   USE module_read_soundings     ! share
   ! note we use p0 as a pass-in invariable in the subroutines of this
   ! module, so a namespace collision seems very likely if we were to
   ! just use the normal name for p0 here, so it's renamed
#ifdef mpas
   USE mpas_constants , ONLY : p0_in_front => p0
#else
   USE module_model_constants , ONLY : p0_in_front => p0
#endif

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: Z2P, P2Z, Z2T, P2T 
   PUBLIC :: INIT_RADIO_SCIENCE_PROFILE, INIT_GCM_PROFILE,  &
             INIT_CARMA_PROFILE, INIT_LH_PROFILE
#ifdef USE_CARMA
   PUBLIC :: get_t_ref_carma
#endif

   ! Give a unique name to each hardwired profile here

   ! radio science is defined here, but has to be read in
   REAL, DIMENSION(1000) :: t_hardwire_rs, z_hardwire_rs, p_hardwire_rs
   REAL, DIMENSION(201) :: t_hardwire_lh, z_hardwire_lh, p_hardwire_lh ! Lellouch and Hunten for Titan
   INTEGER :: n_profile
   ! MER Opportunity (Meridiani) T(z) hardwire
   REAL, DIMENSION(  7), PARAMETER :: t_hardwire_mer = &
        (/ 200., 214., 220., 223.,  223.,  222.,   221. /)
   REAL, DIMENSION(  7), PARAMETER :: z_hardwire_mer = &
        (/   0., 100., 200., 500., 1300., 2000., 20000. /) - 1375.

   ! Venus profile for Lee 2007 case
   REAL, DIMENSION( 31), PARAMETER :: t_hardwire_venus = &
        (/ 725.387,  714.349,  697.988,  678.074,  655.313,  630.378, &
           603.783,  576.020,  547.050,  517.772,  488.579,  459.782, &
           431.657,  404.357,  378.163,  353.205,  329.580,  307.405, &
           286.918,  268.068,  250.948,  235.318,  221.176,  208.874, &
           198.228,  189.304,  181.047,  173.733,  169.569,  167.574, &
           166.711 /)

   REAL, DIMENSION( 31), PARAMETER :: p_hardwire_venus = &
        (/ 8.80463009e+06, 7.94351027e+06, 6.95819049e+06, 5.91468071e+06, &
           4.89693094e+06, 3.96014114e+06, 3.13467132e+06, 2.43156147e+06, &
           1.84920160e+06, 1.37931170e+06, 1.00905778e+06, 7.23236843e+05, &
           5.07174890e+05, 3.47416924e+05, 2.31979950e+05, 1.50638167e+05, &
           9.48700794e+04, 5.77917874e+04, 3.39752926e+04, 1.92414958e+04, &
           1.04899977e+04, 5.50819880e+03, 2.79511939e+03, 1.38050470e+03, &
           6.69459855e+02, 3.24459930e+02, 1.61044965e+02, 8.39949822e+01, &
           4.77124901e+01, 3.10144937e+01, 1.39369974e+01 /)

   ! Titan, T(p) hardwire
   REAL, DIMENSION( 55), PARAMETER :: t_hardwire_titan = &
        (/ 176.0000,   175.8000,   175.6762,   175.5546, &
           175.4352,   175.3181,   175.2516,   175.1952, &
           175.1400,   175.0717,   174.9656,   174.8618, &
           174.7204,   174.5217,   174.2547,   173.8438, &
           172.9802,   172.2556,   171.3730,   170.1814, &
           168.1232,   165.5046,   162.5942,   160.0195, &
           156.0198,   152.6676,   149.6165,   145.3064, &
           140.3683,   135.1095,   127.5253,   118.4804, &
           104.6954,   85.91945,   76.20904,   72.70622, &
           71.58616,   71.22580,   71.10616,   71.41412, &
           71.85484,   73.07359,   74.47054,   76.04823, &
           77.83385,   79.59901,   81.48008,   83.26595, &
           85.22786,   86.98520,   88.86466,   90.35686, &
           91.65218,   92.90943,   93.90000/)
   ! Normally this shouldn't be necessary... p will be p and not eta
   ! 0.000503915 is ptop given ztop=500000. and the other Titan consts.
   REAL, DIMENSION( 55), PARAMETER :: p_hardwire_titan = &
        (/ 8.1944445E-06,  2.1736112E-05, 2.7912094E-05,  3.5699370E-05, &
           4.5540364E-05,  5.7946905E-05, 7.3761097E-05,  9.3533381E-05, &
           1.1799204E-04,  1.4831639E-04, 1.8623908E-04,  2.3315514E-04, &
           2.9094610E-04,  3.6190671E-04, 4.5187571E-04,  5.6924188E-04, &
           7.0812576E-04,  8.7710674E-04, 1.0837355E-03,  1.3325177E-03, &
           1.6378304E-03,  2.0037096E-03, 2.4518929E-03,  2.9939252E-03, &
           3.6594970E-03,  4.4580889E-03, 5.4337773E-03,  6.6158851E-03, &
           8.0593424E-03,  9.7946636E-03, 1.2013729E-02,  1.4802283E-02, &
           1.8443536E-02,  2.3743227E-02, 3.1609237E-02,  4.2539738E-02, &
           5.7112142E-02,  7.5864673E-02, 9.9377923E-02,  0.1293146, &
           0.1653590,      0.2074817,     0.2570405,      0.3127847, &
           0.3742619,      0.4404446,     0.5104085,      0.5821730, &
           0.6543787,      0.7230682,     0.7927206,      0.8550193, &
           0.9097435,      0.9586097,     1.000000/) &
           * (p0_in_front-0.000503915)+0.000503915

   REAL, DIMENSION( 68), PARAMETER :: t_hardwire_jupiter = &
        (/1.22010e+02, 1.22010e+02, 1.22010e+02, 1.19398e+02, 1.14732e+02, &
          1.11588e+02, 1.10111e+02, 1.10350e+02, 1.12253e+02, 1.15680e+02, &
          1.20436e+02, 1.26310e+02, 1.33088e+02, 1.40588e+02, 1.48646e+02, &
          1.57136e+02, 1.65952e+02, 1.74884e+02, 1.83834e+02, 1.92789e+02, &
          2.01755e+02, 2.10726e+02, 2.19711e+02, 2.28708e+02, 2.37702e+02, &
          2.46706e+02, 2.55720e+02, 2.64732e+02, 2.73759e+02, 2.82789e+02, &
          2.91826e+02, 3.00862e+02, 3.09900e+02, 3.18946e+02, 3.27993e+02, &
          3.37046e+02, 3.46099e+02, 3.55158e+02, 3.64217e+02, 3.73283e+02, &
          3.82352e+02, 3.91417e+02, 4.00487e+02, 4.09557e+02, 4.18635e+02, &
          4.27718e+02, 4.36800e+02, 4.45879e+02, 4.54962e+02, 4.64047e+02, &
          4.73129e+02, 4.82219e+02, 4.91313e+02, 5.00407e+02, 5.09500e+02, &
          5.18587e+02, 5.27680e+02, 5.36774e+02, 5.45868e+02, 5.54972e+02, &
          5.64071e+02, 5.73175e+02, 5.82282e+02, 5.91389e+02, 6.00495e+02, &
          6.09597e+02, 6.18706e+02, 6.27807e+02 /)
   REAL, DIMENSION( 68), PARAMETER :: p_hardwire_jupiter = &
        (/1.16355e+03, 2.06397e+03, 2.96439e+03, 3.86480e+03, 5.09360e+03, &
          6.76497e+03, 9.01851e+03, 1.20148e+04, 1.59283e+04, 2.09404e+04, &
          2.72321e+04, 3.49861e+04, 4.43739e+04, 5.55721e+04, 6.87435e+04, &
          8.40693e+04, 1.01719e+05, 1.21883e+05, 1.44773e+05, 1.70577e+05, &
          1.99524e+05, 2.31818e+05, 2.67718e+05, 3.07453e+05, 3.51191e+05, &
          3.99227e+05, 4.51814e+05, 5.09133e+05, 5.71532e+05, 6.39206e+05, &
          7.12444e+05, 7.91436e+05, 8.76483e+05, 9.67905e+05, 1.06591e+06, &
          1.17083e+06, 1.28288e+06, 1.40243e+06, 1.52968e+06, 1.66503e+06, &
          1.80871e+06, 1.96093e+06, 2.12211e+06, 2.29248e+06, 2.47251e+06, &
          2.66245e+06, 2.86252e+06, 3.07294e+06, 3.29424e+06, 3.52667e+06, &
          3.77044e+06, 4.02615e+06, 4.29406e+06, 4.57442e+06, 4.86747e+06, &
          5.17343e+06, 5.49298e+06, 5.82638e+06, 6.17388e+06, 6.53626e+06, &
          6.91325e+06, 7.30566e+06, 7.71376e+06, 8.13781e+06, 8.57808e+06, &
          9.03481e+06, 9.50892e+06, 1.0e+07 /)

CONTAINS

!-----------------------------------------------------------------------------

  REAL FUNCTION Z2P(t_sounding, z, T_ref, T_iso, lapse_rate, g, r_d, &
                    p0, p_scale, topo,do_write) &
  RESULT (p)

    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: t_sounding
    REAL, INTENT(IN) :: z, T_ref, T_iso, lapse_rate, g, R_d, p0, p_scale
    REAL, INTENT(IN), OPTIONAL :: topo
    LOGICAL, INTENT(IN), OPTIONAL :: do_write
  
    ! Local variables
    REAL :: power, z_iso, p_iso, z0

    if(present(do_write)) then
       if(do_write) then
          WRITE( wrf_err_message , * ) "incoming to z2p t_sounding, z, ", &
                 " T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale"
          CALL wrf_message ( TRIM( wrf_err_message ) )
          WRITE( wrf_err_message , * ) t_sounding, z, T_ref, T_iso,  &
                  lapse_rate, g, r_d, p0, p_scale 
          CALL wrf_message ( TRIM( wrf_err_message ) )
       endif
    endif

    z0 = 0.
    IF (PRESENT(topo)) z0 = topo

    SELECT CASE (t_sounding)
    CASE (1) ! Isothermal
       p = p_scale * p0 * EXP(-(g*z)/(r_d*T_ref))
    CASE (2) ! Constant lapse rate
       z_iso = (T_ref-T_iso)/lapse_rate
       power = g/(R_d*lapse_rate)
       IF (z < z_iso) THEN
          p = p_scale * p0 * (1.-(lapse_rate/t_ref)*z)**power
       ELSE
          p_iso = p_scale * p0*(1.-(lapse_rate/t_ref)*z_iso)**power
          p = p_iso * EXP(-(g*(z-z_iso))/(r_d*T_ref))
       END IF
    CASE (3) ! isothermal at specified temperature
       p = p_scale * p0 * EXP(-(g*z)/(r_d*220.))
    CASE (21) ! Venus profile
       p = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_venus, &
                           t_hardwire_in = t_hardwire_venus, &
                           x             = z,           &
                           p0            = p0,          &
                           g             = g,           &
                           R_d           = R_d,         &
                           Z2P           = .TRUE.       )
    CASE (41) ! MER Opportunity (Meridiani) T(z) hardwire
       p = HARDWIRE_T_of_Z(z_hardwire_in = z_hardwire_mer+z0, &
                           t_hardwire_in = t_hardwire_mer,    &
                           x             = z,                 &
                           p0            = p0,                &
                           g             = g,                 &
                           R_d           = R_d,               &
                           p_scale       = p_scale,           &
                           Z2P           = .TRUE.             )
    CASE (42) ! Radio Science T(z) hardwire
       p = HARDWIRE_T_of_Z(z_hardwire_in = z_hardwire_rs(1:n_profile)+z0, &
                           t_hardwire_in = t_hardwire_rs(1:n_profile),    &
                           x             = z,                             &
                           p0            = p0,                            &
                           g             = g,                             &
                           R_d           = R_d,                           &
                           p_scale       = p_scale,                       &
                           Z2P           = .TRUE.                         )
    CASE (43) ! Mars GCM profile
       p = HARDWIRE_T_of_P(p_hardwire_in = z_hardwire_rs(1:n_profile), & ! it says 'z' but it's really'p'
                           t_hardwire_in = t_hardwire_rs(1:n_profile), &
                           x             = z,                          &
                           p0            = p0*p_scale,                 &
                           g             = g,                          &
                           R_d           = R_d,                        &
                           Z2P           = .TRUE.                      )
    CASE (51) ! Jupiter, T(p) hardwire
       IF (p0 /= 1.e7) THEN
          WRITE( wrf_err_message , * ) 'p0 set inappropriately for Jupiter. Should be 1.e7, is ',p0
          CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
       END IF
       p = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_jupiter, &
                           t_hardwire_in = t_hardwire_jupiter, &
                           x             = z,                  &
                           p0            = p0,                 &
                           g             = g,                  &
                           R_d           = R_d,                &
                           Z2P           = .TRUE.              )
    CASE (61) ! Titan, T(p) hardwire
       p = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_titan,  &
                           t_hardwire_in = t_hardwire_titan,  &
                           x             = z,                 &
                           p0            = p0,                &
                           g             = g,                 &
                           R_d           = R_d,               &
                           Z2P           = .TRUE.             )
    CASE (62) ! Titan, read profile using init_carma_profile
       p = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_rs(1:n_profile),  &
                           t_hardwire_in = t_hardwire_rs(1:n_profile),  &
                           x             = z,                           &
                           p0            = p0,                          &
                           g             = g,                           &
                           R_d           = R_d,                         &
                           Z2P           = .TRUE.                       )
    CASE (63) ! Titan, read profile using init_lh_profile
       p = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_lh(1:201),        &      
                           t_hardwire_in = t_hardwire_lh(1:201),        &
                           x             = z,                           &
                           p0            = p0,                          &
                           g             = g,                           &
                           R_d           = R_d,                         &
                           Z2P           = .TRUE.                       )
    CASE DEFAULT
       p = 0.
    END SELECT

  END FUNCTION Z2P

!-----------------------------------------------------------------------------

  REAL FUNCTION P2Z(t_sounding, p, T_ref, T_iso, lapse_rate, g, r_d, &
                    p0, p_scale, topo) &
  RESULT (z)
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: t_sounding
    REAL, INTENT(IN) :: p, T_ref, T_iso, lapse_rate, g, R_d, p0, p_scale
    REAL, INTENT(IN), OPTIONAL :: topo

    ! Local variables
    REAL :: power, rpower, z_iso, p_iso, z0

    z0 = 0.
    IF (PRESENT(topo)) z0 = topo

    SELECT CASE (t_sounding)
    CASE (1) ! Isothermal
       z = -r_d*T_ref*ALOG(p/(p_scale*p0))/g
    CASE (2) ! Constant lapse rate
       z_iso = (T_ref-T_iso)/lapse_rate
       power = g/(R_d*lapse_rate)
       rpower= 1./power
       p_iso = p_scale * p0*(1.-(lapse_rate/t_ref)*z_iso)**power
       IF (p > p_iso) THEN
          z = (1. - (p/(p_scale*p0))**rpower)*T_ref/lapse_rate
       ELSE
          z = -r_d*T_ref*ALOG(p/p_iso)/g + z_iso
       END IF
    CASE (3) ! isothermal at specified temperature
       z = -r_d*220.*ALOG(p/(p_scale*p0))/g
    CASE (21) ! Venus profile
       z = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_venus,  &
                           t_hardwire_in = t_hardwire_venus,  &
                           x             = p,                 &
                           p0            = p0,                &
                           g             = g,                 &
                           R_d           = R_d,               &
                           P2Z           = .TRUE.             )
    CASE (41) ! MER Opportunity (Meridiani) T(z) hardwire
       z = HARDWIRE_T_of_Z(z_hardwire_in = z_hardwire_mer+z0, &
                           t_hardwire_in = t_hardwire_mer,    &
                           x             = p,                 &
                           p0            = p0,                &
                           g             = g,                 &
                           R_d           = R_d,               &
                           p_scale       = p_scale,           &
                           p2z           = .TRUE.             )
    CASE (42) ! Radio Science T(z) hardwire
       z = HARDWIRE_T_of_Z(z_hardwire_in = z_hardwire_rs(1:n_profile)+z0, &
                           t_hardwire_in = t_hardwire_rs(1:n_profile),    &
                           x             = p,                 &
                           p0            = p0,                &
                           g             = g,                 &
                           R_d           = R_d,               &
                           p_scale       = p_scale,           &
                           p2z           = .TRUE.             )
    CASE (43) ! Mars GCM profile
       z = HARDWIRE_T_of_P(p_hardwire_in = z_hardwire_rs(1:n_profile),  & ! The 'z' is really 'p'
                           t_hardwire_in = t_hardwire_rs(1:n_profile),  &
                           x             = p,                           &
                           p0            = p0*p_scale,                  &
                           g             = g,                           &
                           R_d           = R_d,                         &
                           P2Z           = .TRUE.                       )
    CASE (51) ! Jupiter, T(p) hardwire
       IF (p0 /= 1.e7) THEN
          WRITE( wrf_err_message , * ) 'p0 set inappropriately for Jupiter. Should be 1.e7, is ',p0
          CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
       END IF
       z = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_jupiter,  &
                           t_hardwire_in = t_hardwire_jupiter,  &
                           x             = p,                   &
                           p0            = p0,                  &
                           g             = g,                   &
                           R_d           = R_d,                 &
                           p2z           = .TRUE.               )
    CASE (61) ! Titan, T(p) hardwire
       z = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_titan,    &
                           t_hardwire_in = t_hardwire_titan,    &
                           x             = p,                   &
                           p0            = p0,                  &
                           g             = g,                   &
                           R_d           = R_d,                 &
                           p2z           = .TRUE.               )
    CASE (62) ! Titan, read profile using init_carma_profile
       z = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_rs(1:n_profile),    &
                           t_hardwire_in = t_hardwire_rs(1:n_profile),    &
                           x             = p,                             &
                           p0            = p0,                            &
                           g             = g,                             &
                           R_d           = R_d,                           &
                           p2z           = .TRUE.                         )
    CASE (63) ! Titan, read profile using init_lh_profile
       z = HARDWIRE_T_of_P(p_hardwire_in = p_hardwire_lh(1:201),          &
                           t_hardwire_in = t_hardwire_lh(1:201),          &      
                           x             = p,                             &
                           p0            = p0,                            &
                           g             = g,                             &
                           R_d           = R_d,                           &
                           p2z           = .TRUE.                         )

    CASE DEFAULT
       z = 0.
    END SELECT

  END FUNCTION P2Z

!-----------------------------------------------------------------------------

  REAL FUNCTION Z2T(t_sounding, z, T_ref, T_iso, lapse_rate, g, r_d, &
                    p0, p_scale, topo) &
  RESULT (t)

    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: t_sounding
    REAL, INTENT(IN) :: z, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale
    REAL, INTENT(IN), OPTIONAL :: topo

    ! Local variables
    REAL :: z_iso, p, z0

    z0 = 0.
    IF (PRESENT(topo)) z0 = topo

    SELECT CASE (t_sounding)
    CASE (1) ! Isothermal
       t = T_ref
    CASE (2) ! Constant lapse rate
       z_iso = (T_ref-T_iso)/lapse_rate
       IF (z < z_iso) THEN
          t = T_ref-lapse_rate*z
       ELSE
          t = T_iso
       END IF
    CASE (3) ! Isothermal
       t = 220.
    CASE (21) ! Venus profile
       p = Z2P(t_sounding, z, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
       t = P2T(t_sounding, p, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
    CASE (41) ! MER Opportunity (Meridiani) T(z) hardwire
       t = interp_0(t_hardwire_mer, z_hardwire_mer+z0, z, 7)
    CASE (42) ! Radio Science T(z) hardwire
       t = interp_0(t_hardwire_rs(1:n_profile), &
                    z_hardwire_rs(1:n_profile)+z0, z, n_profile)
    CASE (43,51,61,62,63) ! Mars GCM, jupiter, and Titan profiles
       p = Z2P(t_sounding, z, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
       t = P2T(t_sounding, p, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
    CASE DEFAULT
       t = 0.
    END SELECT

  END FUNCTION Z2T

!-----------------------------------------------------------------------------

  REAL FUNCTION P2T(t_sounding, p, T_ref, T_iso, lapse_rate, g, r_d, &
                    p0, p_scale, topo, debug_print) &
      RESULT (t)

    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: t_sounding
    REAL, INTENT(IN) :: p, T_ref, T_iso, lapse_rate, g, R_d, p0, p_scale
    REAL, INTENT(IN), OPTIONAL :: topo
    LOGICAL, INTENT(IN), OPTIONAL :: debug_print
 
    ! Local variables

    REAL :: z, z_iso

    SELECT CASE (t_sounding)
    CASE (1) ! isothermal
       T = T_ref
    CASE (2) ! constant lapse rate
       z = P2Z(t_sounding, p, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
       t = Z2T(t_sounding, z, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
    CASE (21) ! Hardwired Venus profile
       t = interp_0(t_hardwire_venus, ALOG(p_hardwire_venus), ALOG(p), 31)
    CASE (41) ! MER Opportunity (Meridiani) T(z) hardwire
       z = P2Z(t_sounding, p, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
       t = Z2T(t_sounding, z, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
    CASE (42) ! Radio Science T(z) hardwire
       z = P2Z(t_sounding, p, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
       t = Z2T(t_sounding, z, T_ref, T_iso, lapse_rate, g, r_d, p0, p_scale, &
               topo=topo)
    CASE (43) ! Hardwired Mars GCM profile
       t = interp_0(t_hardwire_rs(1:n_profile), ALOG(z_hardwire_rs(1:n_profile)), ALOG(p), n_profile)
    CASE (51) ! Jupiter, T(p) hardwire
       IF (p0 /= 1.e7) THEN
          WRITE( wrf_err_message , * ) 'p0 set inappropriately for Jupiter. Should be 1.e7, is ',p0
          CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
       END IF
       ! Interpolate with respect to ln(p) (more linear)
       t = interp_0(t_hardwire_jupiter, ALOG(p_hardwire_jupiter), ALOG(p), 68)
    CASE (61) ! Hardwired Titan profile
       ! Interpolate with respect to ln(p) (more linear)
       t = interp_0(t_hardwire_titan, ALOG(p_hardwire_titan), ALOG(p), 55)
    CASE (62) ! Hardwired CARMA Titan profile
       t = interp_0(t_hardwire_rs(1:n_profile), ALOG(p_hardwire_rs(1:n_profile)), ALOG(p), n_profile)
    CASE (63) ! Hardwired CARMA Titan profile
       if(present(debug_print)) then
         if(debug_print) then
           write(0,*) "t_hardwire: ",t_hardwire_lh(1:201)
           write(0,*) "log of p_hardwire:",ALOG(p_hardwire_lh(1:201))
           write(0,*) "log of target pressure: ",ALOG(p)
         endif
       endif
       t = interp_0(t_hardwire_lh(1:201), ALOG(p_hardwire_lh(1:201)), ALOG(p), 201)
    CASE DEFAULT
       t = 0.
    END SELECT

  END FUNCTION P2T

!-----------------------------------------------------------------------------
! There can be no T -> z or T -> P mapping since these are non-unique
! (especially in the isothermal case!)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  REAL FUNCTION HARDWIRE_T_of_P(p_hardwire_in, t_hardwire_in, x, p0, g, R_d, &
       z2p, p2z) RESULT (y)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(IN) :: p_hardwire_in, t_hardwire_in
    REAL,               INTENT(IN) :: x, g, R_d, p0
    LOGICAL, OPTIONAL,  INTENT(IN) :: z2p, p2z

    REAL, DIMENSION(SIZE(p_hardwire_in)) :: p_hardwire
    REAL, DIMENSION(SIZE(t_hardwire_in)) :: t_hardwire
    REAL, DIMENSION(SIZE(p_hardwire_in)) :: lp_hardwire
    REAL, DIMENSION(SIZE(p_hardwire_in)+1) :: lp_f, pf, zf
    INTEGER :: nlevels, k, k0

    if ((.not.PRESENT(z2p)) .and. (.not.PRESENT(p2z))) then
      write(wrf_err_message,*) "hardwire_t_of_p function illegally used without purpose (no z2p or p2z)."
      call wrf_error_fatal(TRIM(wrf_err_message))
    else if((PRESENT(z2p)) .and. (PRESENT(p2z))) then
      write(wrf_err_message,*) "hardwire_t_of_p function illegally overloaded (z2p and p2z present)."
      call wrf_error_fatal(TRIM(wrf_err_message))
    endif

    y = 0.

    nlevels = SIZE(p_hardwire_in)
    IF (nlevels /= SIZE(t_hardwire_in)) THEN
       WRITE(wrf_err_message,*) 'module_setup_tpz:hardwire_t_of_p: ', &
            'length of p_hardwire and t_hardwire vectors are not the same:', &
            nlevels, SIZE(t_hardwire_in)
       CALL wrf_error_fatal(TRIM(wrf_err_message))
    END IF
    p_hardwire = p_hardwire_in
    t_hardwire = t_hardwire_in

    ! Make sure P monotonically decreases with Z
    ! We'll send sort2 -p so that it sorts p in decreasing order
    p_hardwire=-p_hardwire
    CALL SORT2(p_hardwire,t_hardwire)
    p_hardwire=-p_hardwire
    ! Find the pressure on the "borders" of a vertical grid that has the
    ! give T(p) at the center.  We will use this to integrate the
    ! hydrostatic equation
    lp_hardwire(:)     = ALOG(p0/p_hardwire(:))
    lp_f(2:nlevels) = 0.5*(lp_hardwire(1:nlevels-1)+lp_hardwire(2:nlevels))
    lp_f(1)         = 2*lp_hardwire(1)-lp_f(2)
    lp_f(nlevels+1) = 2*lp_hardwire(nlevels)-lp_f(nlevels)
    pf(:)           = p0*EXP(-lp_f(:))
    IF (ANY(pf == p0)) THEN
       k0 = 1
       DO k=2,nlevels+1
          IF (p0 == pf(k)) THEN
             k0 = k
             EXIT
          END IF
       END DO
       zf(k0) = 0. ! By definition of p0
       ! Integrate upwards
       DO k = k0, nlevels
          zf(k+1) = zf(k)+R_d*t_hardwire(k)*ALOG(pf(k)/pf(k+1))/g
       END DO
       ! Integrate downwards
       DO k = k0-1, 1, -1
          zf(k) = zf(k+1)-R_d*t_hardwire(k)*ALOG(pf(k)/pf(k+1))/g
       END DO
    ELSE
       k0 = LOCATE(pf,p0)
       ! p0 will be between pf(k0) and pf(k0+1)
       ! Integrate upwards
       ! First from p0 to pf(k0+1)
       zf(k0+1) = R_d*t_hardwire(k0)*ALOG(p0/pf(k0+1))/g
       DO k = k0+1, nlevels
          zf(k+1) = zf(k)+R_d*t_hardwire(k)*ALOG(pf(k)/pf(k+1))/g
       END DO
       ! Integrate downwards
       ! First from p0 to pf(k0)
       zf(k0) = -R_d*t_hardwire(k0)*ALOG(pf(k0)/p0)/g
       DO k = k0-1, 1, -1
          zf(k) = zf(k+1)-R_d*t_hardwire(k)*ALOG(pf(k)/pf(k+1))/g
       END DO
    END IF

    IF (PRESENT(z2p)) THEN
       IF (z2p) THEN
          ! Interp with respect to ln(p) (more linear)
          ! x - > z
          ! y - > p
          y = EXP(interp_0(ALOG(pf), zf, x, nlevels+1))
       END IF
    END IF

    IF (PRESENT(p2z)) THEN
       IF (p2z) THEN
          ! Interp with respect to ln(p) (more linear)
          ! x -> p
          ! y -> z
          y = interp_0(zf, ALOG(pf), ALOG(x), nlevels+1)
       END IF
    END IF

    ! Clean up
  END FUNCTION HARDWIRE_T_of_P

!-----------------------------------------------------------------------------
  REAL FUNCTION HARDWIRE_T_of_Z(z_hardwire_in, t_hardwire_in, x, p0, g, R_d, &
       p_scale, z2p, p2z) RESULT (y)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(IN) :: z_hardwire_in, t_hardwire_in
    REAL,               INTENT(IN) :: x, g, R_d, p0, p_scale
    LOGICAL, OPTIONAL,  INTENT(IN) :: z2p, p2z

    REAL, DIMENSION(SIZE(z_hardwire_in)) :: z_hardwire
    REAL, DIMENSION(SIZE(t_hardwire_in)) :: t_hardwire
    REAL, DIMENSION(SIZE(z_hardwire_in)+1) :: pf, zf
    INTEGER :: nlevels, k, k0

    if ((.not.PRESENT(z2p)) .and. (.not.PRESENT(p2z))) then
      write(wrf_err_message,*) "hardwire_t_of_z function illegally used without purpose (no z2p or p2z)."
      call wrf_error_fatal(TRIM(wrf_err_message))
    else if((PRESENT(z2p)) .and. (PRESENT(p2z))) then
      write(wrf_err_message,*) "hardwire_t_of_z function illegally overloaded (z2p and p2z present)."
      call wrf_error_fatal(TRIM(wrf_err_message))
    endif

    y = 0.

    nlevels = SIZE(z_hardwire_in)
    IF (nlevels /= SIZE(t_hardwire_in)) THEN
       WRITE(wrf_err_message,*) 'module_setup_tpz:hardwire_t_of_z: ', &
            'length of z_hardwire and t_hardwire vectors are not the same:', &
            nlevels, SIZE(t_hardwire_in)
       CALL wrf_error_fatal(TRIM(wrf_err_message))
    END IF
    z_hardwire = z_hardwire_in
    t_hardwire = t_hardwire_in

    ! Make sure z is in increasing order
    CALL SORT2(z_hardwire,t_hardwire)

    zf(2:nlevels) = 0.5*(z_hardwire(1:nlevels-1)+z_hardwire(2:nlevels))
    zf(1)         = 2*z_hardwire(1)-zf(2)
    zf(nlevels+1) = 2*z_hardwire(nlevels)-zf(nlevels)

    pf(1) = p_scale*p0*EXP(-g*zf(1)/(R_d*t_hardwire(1)))
    DO k = 1, nlevels-1
          pf(k+1) = pf(k)*EXP(-g*(zf(k+1)-zf(k))/(R_d*t_hardwire(k)))
    END DO
      
    IF (PRESENT(z2p)) THEN
       IF (z2p) THEN
          ! Interp with respect to ln(p) (more linear)
          ! x - > z
          ! y - > p
          y = EXP(interp_0(ALOG(pf), zf, x, nlevels+1))
       END IF
    END IF

    IF (PRESENT(p2z)) THEN
       IF (p2z) THEN
          ! Interp with respect to ln(p) (more linear)
          ! x -> p
          ! y -> z
          y = interp_0(zf, ALOG(pf), ALOG(x), nlevels+1)
       END IF
    END IF

    ! Clean up

  END FUNCTION HARDWIRE_T_of_Z

  SUBROUTINE INIT_RADIO_SCIENCE_PROFILE(in_file)

    implicit none

    CHARACTER(LEN=100) :: in_file 

    CALL READ_SOUNDINGS( in_file, 3, 42,z_profile_out=z_hardwire_rs, &
                                        t_profile_out=t_hardwire_rs, &
                                        n_profile = n_profile)

    if(n_profile .le. 0) then
          CALL wrf_error_fatal('Something went wrong reading the profile, no points were read.')
    endif

    RETURN

  END SUBROUTINE INIT_RADIO_SCIENCE_PROFILE

  SUBROUTINE INIT_GCM_PROFILE(in_file)

    ! so as not to define a new variable, we reuse z_profile, but these are
    ! really pressure profiles T(p) (not the T(z) used in RS)

    implicit none

    CHARACTER(LEN=100) :: in_file

    WRITE( wrf_err_message , * ) "setting up input profile from prior GCM output: ",in_file
    CALL wrf_message ( TRIM( wrf_err_message ) )
 

    CALL READ_SOUNDINGS( in_file, 3, 43,z_profile_out=z_hardwire_rs, &
                                        t_profile_out=t_hardwire_rs, &
                                        n_profile = n_profile)

    if(n_profile .le. 0) then
          CALL wrf_error_fatal('Something went wrong reading the profile, no points were read.')
    endif

    RETURN

  END SUBROUTINE INIT_GCM_PROFILE

  SUBROUTINE INIT_LH_PROFILE

    ! so as not to define a new variable, we reuse z_profile, but these are
    ! really pressure profiles T(p) (not the T(z) used in RS)

    implicit none

    CHARACTER(LEN=41) :: in_file = "./Data/titan/CARMA/Titan/LH_rec_model.txt"
    CHARACTER(len=3) :: cjunk
    REAL :: a,b
    integer ::  luni, k

    WRITE( wrf_err_message , * ) "setting up input profile from Lellouch and Hunten 1987: ",in_file
    CALL wrf_message ( TRIM( wrf_err_message ) )

    luni = 10
    open(unit=luni,file=in_file,status='old')
    do k=1,3
      read(luni,*) cjunk
    enddo
    do k=1,201
      read(luni,*) z_hardwire_lh(k),p_hardwire_lh(k),a, t_hardwire_lh(k),b
    enddo
    close(luni)

    z_hardwire_lh(:) = z_hardwire_lh(:) * 1000.  ! convert km->m
    p_hardwire_lh(:) = p_hardwire_lh(:) /10.     ! convert dyn/cm2 to Pa

!    do k=1,201
!      write(0,*) k, z_hardwire_lh(k),p_hardwire_lh(k),t_hardwire_lh(k)
!    enddo

    RETURN

  END SUBROUTINE INIT_LH_PROFILE


  SUBROUTINE INIT_CARMA_PROFILE(carma_profile  &
#ifdef USE_CARMA
                        ,ingest_carma_outfile  &
#endif 
                  )

    implicit none

    integer, intent(in) :: carma_profile
#ifdef USE_CARMA
    character(len=*), intent(in) :: ingest_carma_outfile
#endif 

    WRITE( wrf_err_message , * ) "setting up input profile using CARMA code: ",carma_profile
    CALL wrf_message ( TRIM( wrf_err_message ) )
 

    CALL READ_CARMA_PROFILES( load_atmos= carma_profile, &
                                       p=p_hardwire_rs,  &
                                       t=t_hardwire_rs,  &
                                       n=n_profile       &
#ifdef USE_CARMA
             ,ingest_carma_outfile=ingest_carma_outfile  &
#endif 
                                       )

    if(n_profile .le. 0) then
          CALL wrf_error_fatal('Something went wrong reading the profile, no points were read.')
    endif

    RETURN

  END SUBROUTINE INIT_CARMA_PROFILE

#ifdef USE_CARMA
  SUBROUTINE get_t_ref_carma(ingest_carma_outfile,         &
                             titan_carma_load_atmos,       &
                             t_ref_carma_c,                &
                             t_ref_carma_e,                &
                             p_ref_carma_c,                &
                             p_ref_carma_e,                &
                             znw, znu, p_top, p_surf,      &
                             kms, kme, kts, kte)

   implicit none

   ! incoming variables
   character(len=*), intent(in) :: ingest_carma_outfile
   integer, intent(in) ::          titan_carma_load_atmos
   integer, intent(in) :: kms, kme, kts, kte
   real, intent(in) :: p_surf, p_top
   real, dimension(kms:kme) :: znw, znu

   ! outgoing variables
   real, intent(out), dimension(kms:kme) :: t_ref_carma_e  ! layer edge temps
   real, intent(out), dimension(kms:kme) :: t_ref_carma_c  ! layer centre temps
   real, intent(out), dimension(kms:kme) :: p_ref_carma_e  ! layer edge pressures
   real, intent(out), dimension(kms:kme) :: p_ref_carma_c  ! layer centre pressures

   ! locals
   real :: p_surf_lower, p_top_higher
   integer :: k

   ! add some margin - this is reference profile, not the true WRF (or CARMA) profile:
   p_surf_lower = p_surf*1.2   ! 20% higher pres than nominal surface
   p_top_higher = p_top*0.8    ! 80% of the pres at the nominal wrf top

   call  INIT_CARMA_PROFILE(titan_carma_load_atmos, ingest_carma_outfile)
   do k=kts,kte-1
     p_ref_carma_e(k) = (znw(k)*(p_surf_lower - p_top_higher)) + p_top_higher
     t_ref_carma_e(k) = P2T(t_sounding = 62,               &
                            p          = p_ref_carma_e(k), & 
                            T_ref      = 0.,               &  ! zeros set as unsed for case 62
                            T_iso      = 0.,               &
                            lapse_rate = 0.,               &
                            g          = 0.,               &
                            r_d        = 0.,               &
                            p0         = 0.,               &
                            p_scale    = 0.                )
            
     p_ref_carma_c(k) = (znw(k)*(p_surf_lower - p_top_higher)) + p_top_higher
     t_ref_carma_c(k) = P2T(62,                            &
                            p          = p_ref_carma_c(k), & 
                            T_ref      = 0.,               &  ! zeros set as unsed for case 62
                            T_iso      = 0.,               &
                            lapse_rate = 0.,               &
                            g          = 0.,               &
                            r_d        = 0.,               &
                            p0         = 0.,               &
                            p_scale    = 0.                )
   enddo
   p_ref_carma_e(kte) = (znw(kte)*(p_surf_lower - p_top_higher)) + p_top_higher
   t_ref_carma_e(kte) = P2T(t_sounding = 62,               &
                            p          = p_ref_carma_e(kte), & 
                            T_ref      = 0.,               &  ! zeros set as unsed for case 62
                            T_iso      = 0.,               &
                            lapse_rate = 0.,               &
                            g          = 0.,               &
                            r_d        = 0.,               &
                            p0         = 0.,               &
                            p_scale    = 0.                )
   return

  END SUBROUTINE get_t_ref_carma
#endif

END MODULE module_setup_tpz
