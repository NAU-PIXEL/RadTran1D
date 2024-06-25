MODULE module_planet_utilities

      USE module_model_constants
      USE module_wrf_error

        ! Parameter variables
  REAL, PARAMETER :: AA0    = 6.107799961
  REAL, PARAMETER :: AA1    = 4.436518521e-01
  REAL, PARAMETER :: AA2    = 1.428945805e-02
  REAL, PARAMETER :: AA3    = 2.650648471e-04
  REAL, PARAMETER :: AA4    = 3.031240396e-06
  REAL, PARAMETER :: AA5    = 2.034080948e-08
  REAL, PARAMETER :: AA6    = 6.136820929e-11
  REAL, PARAMETER :: CC1    = 9.09718
  REAL, PARAMETER :: CC2    = 3.56654
  REAL, PARAMETER :: CC3    = 0.876793
  REAL, PARAMETER :: EIS    = 6.1071
  REAL, PARAMETER :: TC     = 273.16
  REAL, PARAMETER :: LN10   = LOG(10.)

  REAL, PARAMETER :: pi = pi_s

CONTAINS

      
!----------------------------------------------------------------
      REAL FUNCTION get_julian(ls,eccentricity,  &
      equinox_fraction,zero_date) RESULT (julian)
!----------------------------------------------------------------
!
! Calculate Julian day for a given Ls.
! The definition used for this subroutine (and similarly for
! the calculation in phys/module_radation_driver.F:radconst is
! midnight at the beginning of the first day of the year is
! equal to a "Julian day" of 0.0....
!
! Input: Solar longitude (degrees)
!
! Output: Julian day (sols, fractional)
!
!----------------------------------------------------------------
    USE module_model_constants
    IMPLICIT NONE
!----------------------------------------------------------------

!   Input/Ouptut variables
      REAL, INTENT(IN) :: ls, & ! Day of the year, sols
                eccentricity, &
            equinox_fraction, &
                   zero_date

!   Parameter variables
    REAL(KIND (0d0)), PARAMETER :: SMALL_VALUE = 1.D-6
    INTEGER, PARAMETER :: PLANET_YEAR = 669
      
!   Local Variables
    REAL(KIND (0d0)) :: pi
    REAL(KIND (0d0)) :: deleqn, date_dbl
    REAL(KIND (0d0)) :: er, qq, e, cd0, ep, em
    REAL(KIND (0d0)) :: eq, w, als, ajulian, dp_date

    PI=ACOS(-1.d0)

    deleqn = equinox_fraction * REAL(PLANET_YEAR) 

    er = SQRT( (1.d0+eccentricity)/(1.d0-eccentricity) )

    !  qq is the mean anomaly
    qq = 2.d0 * (pi * deleqn / REAL(PLANET_YEAR))

    !  determine true anomaly at equinox:  eq
    !  Iteration for eq
    e = 1.d0
    cd0 = 1.d0
    DO WHILE (cd0 > SMALL_VALUE)
       ep = e - (e-eccentricity*SIN(e)-qq)/(1.d0-eccentricity*COS(e))
       cd0 = ABS(e-ep)
       e = ep
    END DO
    eq = 2.d0 * ATAN( er * TAN(0.5d0*e) )

    w = eq + real(ls,kind(0.d0))*pi/180.

    ! e is the eccentric anomaly
    e = 2.0*atan((tan(w*0.5))/er)

    em = e - eccentricity*sin(e)

    dp_date = em * planet_year / (2.0 * pi)

    ajulian = dp_date + zero_date

    if(ajulian .lt. 0)           ajulian = ajulian + planet_year
    if(ajulian .gt. planet_year) ajulian = ajulian - planet_year
    
    julian = REAL(ajulian)
   
  END FUNCTION get_julian

!----------------------------------------------------------------
  REAL FUNCTION emissivity_lw(emiss, n_ices, surface_ice, &
#if ( WRF_MARS == 1 )
                       CO2ICE_THRESHOLD, &
                       H2OICE_THRESHOLD, &
#endif
#if ( WRF_PLUTO )
                       N2ICE_THRESHOLD, &
                       CH4ICE_THRESHOLD, &
#endif 
                       latitude, temperature, emiss_ice)
!----------------------------------------------------------------
!
! Calculate the long-wave infrared surface emissivity.
!
! Wrapper function for calculations based on surface ices
! and latitudes appropriate to each planet
!
! Input: emiss (unitless)
!          2D field of (typically bedrock) surface emissivity
!        suface ice (kg/m^2)
!          2D field of surface ice column density
!        latitude (degrees N)
!          latitude of grid point
!        temperature (K)
!          surface temperature of grid point
!
! Output: emissivity (unitless)
!
!----------------------------------------------------------------
    USE module_wrf_error
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input variables
    INTEGER,                 INTENT(IN) :: n_ices
    REAL,                    INTENT(IN) :: emiss
    REAL,                    INTENT(IN) :: latitude
    REAL,                    INTENT(IN) :: temperature
    REAL, DIMENSION(n_ices), INTENT(IN) :: surface_ice
    REAL, DIMENSION(n_ices), INTENT(IN), OPTIONAL :: emiss_ice
#if ( WRF_MARS == 1 )
    REAL,                    INTENT(IN) :: CO2ICE_THRESHOLD, &
                                           H2OICE_THRESHOLD
#endif
#if ( WRF_PLUTO == 1 ) 
    REAL,                    INTENT(IN) :: N2ICE_THRESHOLD, &
                                           CH4ICE_THRESHOLD
#endif

! Local variables
    ! ice_threshold (kg/m^2)
    !  critical value of surface ice column density
    !  above which the emissivity of the ice (as opposed
    !  to say underlying bedrock) should be used
    REAL, DIMENSION(n_ices) :: ice_threshold
    LOGICAL :: use_ice_emissivities

#ifdef WRF_MARS
    ! Infrared emissivity of CO2 ice
    !REAL, PARAMETER :: emiss_mars_co2ice= 0.8
    ! Better fits to Viking lander pressure curves
    REAL, PARAMETER :: emiss_mars_co2ice_north = 0.485
    REAL, PARAMETER :: emiss_mars_co2ice_south = 0.785

    ! Infrared emissivity of H2O ice
    REAL, PARAMETER :: emiss_mars_h2oice= 1.

    ice_threshold = (/ CO2ICE_THRESHOLD, &
                       H2OICE_THRESHOLD /)
#endif
#ifdef WRF_PLUTO
    IF (n_ices < 2) THEN
       WRITE(wrf_err_message,*) 'Not enough ices and/or ice emissivities '//&
            'specified!'
       CALL wrf_error_fatal(TRIM(wrf_err_message))
    END IF
    ice_threshold = (/ N2ICE_THRESHOLD, &
                       CH4ICE_THRESHOLD /)
    use_ice_emissivities = .FALSE.
    IF (PRESENT(emiss_ice)) use_ice_emissivities = .TRUE.
#endif

    ! Default value
    emissivity_lw = emiss

#if defined ( WRF_MARS )
    ! Check on CO2 ice
    IF ( surface_ice(1) > ice_threshold(1) ) THEN
       IF (latitude > 0.) THEN
          emissivity_lw = emiss_mars_co2ice_north
       ELSE
          emissivity_lw = emiss_mars_co2ice_south
       END IF
    ! Modify longwave surface emissivity for H2O ice here, if necessary
    !ELSE IF (surface_ice(2) > ice_threshold(2)) THEN
    !   emissivity_lw=emiss_h2oice
    END IF
#elif defined ( WRF_PLUTO )
#if 0
    IF ( surface_ice(1) > ice_threshold(1) ) THEN
       ! Stansberry, J. A., and R. V. Yelle, "Emissivity and the Fate of
       ! Pluto's Atmosphere", Icarus, 141, pp. 299-306, 1999
       IF (temperature < 35.6) THEN
          emissivity_lw = 0.30 ! (N2(s), T<35.6K).
       ELSE
          emissivity_lw = 0.75 ! (N2(s), T>35.6K)
       END IF
    END IF
#else
    IF (use_ice_emissivities) THEN
       IF ( surface_ice(2) > ice_threshold(2) ) THEN
          emissivity_lw = emiss_ice(2)
       ELSE IF ( surface_ice(1) > ice_threshold(1) ) THEN
          emissivity_lw = emiss_ice(1)
       ELSE
          emissivity_lw = emiss
       END IF
    END IF
#endif
#elif defined ( WRF_TITAN )
#elif defined ( WRF_VENUS )
#elif defined ( WRF_EARTH )
#endif

  END FUNCTION emissivity_lw

!----------------------------------------------------------------
   REAL FUNCTION esat(temperature)
!----------------------------------------------------------------
!
! Calculate saturation vapour pressure of water using a formulas
! that are valid over ice at low temperatures and over ice and
! liquid at higher pressures
!
! Input: temperature (K)
!
! Output: Saturation vapor pressure (Pa)

!
!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
   REAL, INTENT(IN   ) ::  temperature

    esat = PFROST_H2O(temperature,.FALSE.)

   END FUNCTION esat

!----------------------------------------------------------------
  REAL FUNCTION rsat_h2o(temperature,pressure)
!----------------------------------------------------------------
!
! Calculate saturation vapor pressure of water using a formulas
! that are valid over ice at low temperatures and over ice and
! liquid at higher pressures
!
! Input: temperature (K), pressure (Pa)
!
! Output: Saturation mass mixing ratio (kg/kg)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL, INTENT(IN   ) ::  temperature, pressure

!  Local variables
    REAL :: esat

    esat = PFROST_H2O(temperature,.FALSE.)
    ! With the following definition of rsat, the assumption is that the
    ! moisture variables (Q_vapor, Q_ice, etc.) refer to mass mixing
    ! ratio per mass of the _TOTAL_ atmosphere, not to merely the dry
    ! atmosphere.  This issue arises in planetary atmospheres because,
    ! for Mars in particular, the saturation vapor pressure of water
    ! is nearly equal to the total atmospheric pressure at some temperatures,
    ! and thus using the defintion of mass mixing ration per mass of dry
    ! air would give a value of infinity ( something/(x-x) = 1/0 = infinity)
    ! which would be useless (as well as crash the model).  Hopefully
    ! the only place this difference in definition should matter is here
    ! (with the saturation values), in the surface pressure fluxes (we
    ! account for dry air mass changes, e.g., CO2, but not for pressure
    ! changes due to H2O vapor pressure) and in adding the water vapor
    ! amounts to the total atmospheric pressure.
    rsat_h2o = MW_H2O * esat / (MW_AIR * pressure)
    ! WRONG!: Limit saturation between 0 and 1
    !rsat_h2o = MAX(MIN(rsat_h2o,1.),0.)
    ! Limit saturation to positive values
    rsat_h2o = MAX(rsat_h2o,0.)

  END FUNCTION rsat_h2o

!----------------------------------------------------------------
  REAL FUNCTION rsat_ch4(temperature,pressure)
!----------------------------------------------------------------
!
! Calculate saturation mixing ratio of methane
!
! Input: temperature (K), pressure (Pa)
!
! Output: Saturation mass mixing ratio (kg/kg)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL, INTENT(IN   ) ::  temperature, pressure

!  Local variables
    REAL :: esat

    ! See full comments on methodology in function RSAT_H2O
    esat = PFROST_CH4(temperature,.FALSE.)
    rsat_ch4 = MW_CH4 * esat / (MW_AIR * pressure)
    ! WRONG!: Limit saturation between 0 and 1
    !rsat_ch4 = MAX(MIN(rsat_ch4,1.),0.)
    ! Limit saturation to positive values
    rsat_ch4 = MAX(rsat_ch4,0.)

  END FUNCTION rsat_ch4

!----------------------------------------------------------------
  REAL FUNCTION molecular_diffusivity(temperature, pressure) &
       RESULT (diffusivity)
!----------------------------------------------------------------
!
! Calculate the molecular self-diffusivity of a gas at a given
! temperature and pressure.
!
! Input: Temperature (K), Pressure (Pa), Planet identifier ID
!
! Output: diffusivity (m^2/s)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!   Input/Ouptut variables
    REAL,    INTENT(IN) :: temperature  ! K
    REAL,    INTENT(IN) :: pressure     ! Pa

!   Parameter variables
    ! Also prevent roundoff and over/underflow errors
    REAL, PARAMETER :: c1 = k_boltzmann/(molecular_diameter*molecular_diameter)
    REAL, PARAMETER :: c2 = k_boltzmann/(pi*mw_air*m_amu)

!   Local Variables
    REAL :: term

    term = SQRT(c2 * temperature)
    diffusivity = term * c1 * (3. * temperature) / (8. * pressure)

  END FUNCTION molecular_diffusivity

!----------------------------------------------------------------
  REAL FUNCTION tfrost(pressure)
!----------------------------------------------------------------
!
! Calculate the deposition ("frost" = "meteorological deposition",
! <http://en.wikipedia.org/wiki/Sublimation_(chemistry)>)
! temperature of a gas as a function of pressure.
!
! Wrapper function for the gas appropriate to each planet
! Input: pressure (Pa)
!
! Output: frost temperature (K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL, INTENT(IN   ) :: pressure

#if defined WRF_MARS
    tfrost = tfrost_co2(pressure)
#elif defined WRF_PLUTO
    tfrost = tfrost_n2(pressure)
#elif defined WRF_TITAN
    ! Return flag value
    tfrost = 0.
#elif defined WRF_VENUS
    ! Return flag value
    tfrost = 0.
#elif defined WRF_EARTH
    ! Return flag value
    tfrost = 0.
#endif

  END FUNCTION tfrost

!----------------------------------------------------------------
  REAL FUNCTION tfrost_h2o(pressure)
!----------------------------------------------------------------
!
! Calculate the deposition ("frost" = "meteorological deposition",
! <http://en.wikipedia.org/wiki/Sublimation_(chemistry)>)
! temperature of H2O as a function of pressure.
!
! Input: pressure (Pa)
!
! Output: H2O frost temperature (K)
!
!----------------------------------------------------------------
    USE module_nrutils, ONLY : zroots, zbrak, rtsafe
    USE module_wrf_error
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL, INTENT(IN   ) :: pressure

!  Local variables
    REAL(KIND(0d0)) :: p
    COMPLEX(KIND((0d0,0d0))), DIMENSION(7) :: coeffs
    COMPLEX(KIND((0d0,0d0))), DIMENSION(6) :: roots
    REAL(KIND(0d0)) :: x1, x2, xl, xh
    INTEGER :: n, nb
    REAL(KIND(0d0)), DIMENSION(:), POINTER :: xb1, xb2

    !  esat is in mbar -- need to convert to Pa
    p = pressure * 1.d-2

    IF ( p > aa0 ) THEN
       coeffs(1) = CMPLX(REAL(AA0,KIND(0d0))-p,KIND=KIND(0d0))
       coeffs(2) = CMPLX(AA1,KIND=KIND(0d0))
       coeffs(3) = CMPLX(AA2,KIND=KIND(0d0))
       coeffs(4) = CMPLX(AA3,KIND=KIND(0d0))
       coeffs(5) = CMPLX(AA4,KIND=KIND(0d0))
       coeffs(6) = CMPLX(AA5,KIND=KIND(0d0))
       coeffs(7) = CMPLX(AA6,KIND=KIND(0d0))
       CALL ZROOTS(coeffs, roots, .TRUE.)
       ! The correct root is the larger of the two real roots (i.e., no
       ! imaginary part)
       tfrost_h2o = TC + MAXVAL(REAL(roots,KIND(0.)), &
                                MASK=(ABS(AIMAG(roots)) < 1.d-6))
    ELSE
       ! Find a 1 K bracketing range to search for root to invert equation
       x1 = 50.d0
       x2 = 274.d0
       n  = 224
       CALL ZBRAK(PFROST_H2O_ROOT,x1,x2,n,p,xb1,xb2,nb)
       IF (nb == 0) THEN
          ! No roots found in our range!  Something massively wrong!
          WRITE(wrf_err_message,*) 'No roots found in the range:',x1, &
               ' to ',x2
          CALL wrf_error_fatal(TRIM(wrf_err_message))
       ELSE IF (nb > 1) THEN
          ! We probably have a root that sits exactly on the edge of our
          ! attempted bracketing range.  Compensate by extending range by
          ! 1 K on either side
          xl = xb1(1) - 1
          xh = xb2(1) + 1
          WRITE(wrf_err_message,*) 'Found more than one root with zbrak:',nb
!c-mm          CALL WRF_debug(100,TRIM(wrf_err_message))
       ELSE
          xl = xb1(1)
          xh = xb2(1)
       END IF
       DEALLOCATE(xb1,xb2)
       NULLIFY(xb1,xb2)
       tfrost_h2o = REAL(RTSAFE(PFROST_H2O_ROOT_DERIV,xl,xh,1.d-6,p), &
                         KIND=KIND(0.))
    END IF

  END FUNCTION tfrost_h2o

!----------------------------------------------------------------
  REAL FUNCTION tfrost_co2(pressure)
!----------------------------------------------------------------
!
! Calculate the deposition ("frost" = "meteorological deposition",
! <http://en.wikipedia.org/wiki/Sublimation_(chemistry)>)
! temperature of CO2 as a function of pressure for a range that is
! typical of Martian pressures.
!
! Input: pressure (Pa)
!
! Output: CO2 frost temperature (K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL, INTENT(IN   ) :: pressure

!  Parameter variables
   ! for equation of the form y = A/(B-ln(C*x))

    REAL, PARAMETER :: C = 0.9532  ! This represents the fraction of Mars
                                   ! atmospheric gas that is CO2 (as opposed
                                   ! to Ar, N2, etc.)  If necessary, this
                                   ! could become an input variable as well...

    ! James et al., 1992 formula
    !REAL, PARAMETER :: A = 3182.48
    !REAL, PARAMETER :: B = 23.3493

    ! MCS team formula (Armin Kleinbohl fit to CRC Handbook data)
    ! This forumula used base 10 logarithm, which can be converted to a 
    ! natural logarithm by noting that log10 x = ln x / ln 10.
    ! Thus y = A2/(B2-log10(C2*C*x)) is equivalent to y = A/(B-ln(C*x))
    ! where:
    ! A = A2 ln 10
    ! B = B2 ln 10 - ln C2
    REAL, PARAMETER :: A = 3148.4278 ! (1367.3448 * ln 10)
    REAL, PARAMETER :: B = 27.707244 ! (9.9082 * ln 10 - ln 0.00750061674)

    tfrost_co2 = A / ( B - LOG( C * pressure ) )

    ! Old formula, diverges from others at low pressures
    !tfrost_co2 = 149.2 + 6.48 * LOG( 1.35e-3 * pressure )

  END FUNCTION tfrost_co2

!----------------------------------------------------------------
  REAL FUNCTION tfrost_n2(pressure)
!----------------------------------------------------------------
!
! Calculate the deposition ("frost" = "meteorological deposition",
! <http://en.wikipedia.org/wiki/Sublimation_(chemistry)>)
! temperature of N2 as a function of pressure for a range that is
! typical of the outer planets and satellites
!
! Input: pressure (Pa)
!
! Output: N2 frost temperature (K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL, INTENT(IN   ) :: pressure

!  Parameter variables
    ! For equation of the form log_10(P) = A + B/T + C/T^2
    ! Invert to solve for T(P):
    ! T = 2*C / (-B - SQRT(B^2 - 4*C*(A - LOG10(P))))

    REAL, PARAMETER :: A =  7.98698
    REAL, PARAMETER :: B = -161.597 ! K
    REAL, PARAMETER :: C = -5363.80 ! K^2

    IF (pressure <= 0.) THEN
       tfrost_n2 = 0. ! in an asymptotic sense...
    ELSE
       tfrost_n2 = 2.d0*C / (-B - SQRT(B*B - 4.d0*C*(A - LOG10(pressure))))
    END IF

  END FUNCTION tfrost_n2

!----------------------------------------------------------------
  REAL FUNCTION tfrost_ch4(pressure)
!----------------------------------------------------------------
!
! Calculate the deposition ("frost" = "meteorological deposition",
! <http://en.wikipedia.org/wiki/Sublimation_(chemistry)>)
! temperature of CH4 as a function of pressure for a range that is
! typical of the outer planets and satellites
!
! Input: pressure (Pa)
!
! Output: CH4 frost temperature (K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL, INTENT(IN   ) :: pressure

!  Parameter variables
    ! For equation of the form P = A e^(-B/(T+C))
    !pfrost_ch4 = A * EXP(-B/(temperature+C))
    ! Invert to solve for T(P):
    ! T = -C + B/(ln(A) - ln(P))
    REAL, PARAMETER :: A = 133.3224*EXP(6.61184*LN10) ! Pa
    REAL, PARAMETER :: B = 389.93*LN10 ! K
    REAL, PARAMETER :: C = 266.-273.15 ! K

    IF (pressure <= 0.) THEN
       tfrost_ch4 = 0. ! in an asymptotic sense...
    ELSE
       tfrost_ch4 = -C + B/(LOG(A)-LOG(pressure))
    END IF

  END FUNCTION tfrost_ch4

!----------------------------------------------------------------
  REAL FUNCTION pfrost(temperature, dpdt)
!----------------------------------------------------------------
!
! Calculate the equilibrium vapor pressure of a gas as a function
! of temperature, unless the "dpdt" optional parameter is set to
! .TRUE., in which case the function returns the value of dP/dT
! (the Clausius-Clapeyron relation) of the equilibrium vapor
! pressure equation at the input temperature
!
! Wrapper function for the gas appropriate to each planet
! Input: temperature (K)
!
! Output: pressure (Pa)
!         OR
!         pressure/temperature (Pa/K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL,    INTENT(IN)           :: temperature
    LOGICAL, INTENT(IN), OPTIONAL :: dpdt

    LOGICAL :: calculate_dpdt
    calculate_dpdt = .FALSE.
    IF (PRESENT(dpdt)) calculate_dpdt = dpdt

#if defined WRF_MARS
    pfrost = pfrost_co2(temperature, calculate_dpdt)
#elif defined WRF_PLUTO
    pfrost = pfrost_n2( temperature, calculate_dpdt)
#elif defined WRF_TITAN
    ! Return flag value
    pfrost = 0.
#elif defined WRF_VENUS
    ! Return flag value
    pfrost = 0.
#elif defined WRF_EARTH
    ! Return flag value
    pfrost = 0.
#endif

  END FUNCTION pfrost

!----------------------------------------------------------------
  REAL FUNCTION pfrost_h2o(temperature,calculate_dpdt)
!----------------------------------------------------------------
!
! Calculate the equilibrium vapor pressure of H2O gas as a function
! of temperature, unless the "calculate_dpdt" optional parameter is
! set to .TRUE., in which case the function returns the value of
! dP/dT (the Clausius-Clapeyron relation) of the equilibrium vapor
! pressure equation at the input temperature
!
! Input: temperature (K)
!
! Output: pressure (Pa)
!         OR
!         pressure/temperature (Pa/K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL,    INTENT(IN) :: temperature
    LOGICAL, INTENT(IN) :: calculate_dpdt

!  Local variables
    REAL :: t1, esat, rhs

    IF ( temperature > TC ) THEN
       t1 = temperature - TC
       esat =  AA0 + t1 * (AA1 + t1 * (AA2 + t1 * (AA3 + t1 * &
              (AA4 + t1 * (AA5 + t1 *  AA6)))))
       esat = MAX(esat,0.)
       IF (calculate_dpdt) THEN
          esat = AA1 + t1 * (2.*AA2 + t1 * (3.*AA3 + t1 * &
                 (4.*AA4 + t1 * (5.*AA5 + t1 * 6.*AA6))))
       END IF
    ELSE
       rhs =  -CC1 * (TC / temperature - 1.)  &
             - CC2 * LOG10(TC / temperature)  &
             + CC3 * (1. - temperature / TC) + LOG10(EIS)
       esat = 10.**rhs
       esat = MAX(esat,0.)
       IF (calculate_dpdt) THEN
          esat = esat * &
                 (CC1*TC*LN10/(temperature*temperature) + &
                  CC2/temperature - CC3*LN10/TC)
       END IF
    END IF
    !  esat is in mbar -- need to convert to Pa
    pfrost_h2o = 100.*esat

  END FUNCTION pfrost_h2o

!----------------------------------------------------------------
  REAL FUNCTION pfrost_co2(temperature,calculate_dpdt)
!----------------------------------------------------------------
!
! Calculate the equilibrium vapor pressure of CO2 gas as a function
! of temperature for a range that is typical for Mars, unless the
! "calculate_dpdt" optional parameter is set to .TRUE., in which
! case the function returns the value of dP/dT (the
! Clausius-Clapeyron relation) of the equilibrium vapor pressure
! equation at the input temperature
!
! Input: temperature (K)
!
! Output: pressure (Pa)
!         OR
!         pressure/temperature (Pa/K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL,    INTENT(IN) :: temperature
    LOGICAL, INTENT(IN) :: calculate_dpdt

!  Parameter variables
   ! for equation of the form y = A/(B-ln(C*x))

    REAL, PARAMETER :: C = 0.9532  ! This represents the fraction of Mars
                                   ! atmospheric gas that is CO2 (as opposed
                                   ! to Ar, N2, etc.)  If necessary, this
                                   ! could become an input variable as well...

    ! James et al., 1992 formula
    !REAL, PARAMETER :: A = 3182.48
    !REAL, PARAMETER :: B = 23.3493

    ! MCS team formula (Armin Kleinbohl fit to CRC Handbook data)
    ! This forumula used base 10 logarithm, which can be converted to a 
    ! natural logarithm by noting that log10 x = ln x / ln 10.
    ! Thus y = A2/(B2-log10(C2*C*x)) is equivalent to y = A/(B-ln(C*x))
    ! where:
    ! A = A2 ln 10
    ! B = B2 ln 10 - ln C2
    REAL, PARAMETER :: A = 3148.4278 ! (1367.3448 * ln 10)
    REAL, PARAMETER :: B = 27.707244 ! (9.9082 * ln 10 - ln 0.00750061674)

    pfrost_co2 = EXP(B - A/temperature) / C
    IF (calculate_dpdt) THEN
       pfrost_co2 = pfrost_co2 * A/temperature/temperature
    ENDIF

    ! Old formula, diverges from others at low pressures
    !pfrost_co2 = EXP( (temperature - 149.2) / 6.48 ) / 1.35e-3

  END FUNCTION pfrost_co2

!----------------------------------------------------------------
  REAL FUNCTION pfrost_n2(temperature,calculate_dpdt)
!----------------------------------------------------------------
!
! Calculate the equilibrium vapor pressure of N2 gas as a function
! of temperature for a range that is typical of the outer planets
! and satellites, unless the "calculate_dpdt" optional parameter
! is set to .TRUE., in which case the function returns the value
! of dP/dT (the Clausius-Clapeyron relation) of the equilibrium
! vapor pressure equation at the input temperature
!
! Input: temperature (K)
!
! Output: pressure (Pa)
!         OR
!         pressure/temperature (Pa/K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL,    INTENT(IN) :: temperature
    LOGICAL, INTENT(IN) :: calculate_dpdt

!  Parameter variables
    ! log_10(P) = A + B/T + C/T^2
    ! Rewrite in natural logarithm:
    ! ln P = A ln 10 + B ln 10 / T + C ln 10 / T^2
    ! ln P = A' + B'/T + C'/T^2
    
    !REAL, PARAMETER :: A =  7.98698
    !REAL, PARAMETER :: B = -161.597 ! K
    !REAL, PARAMETER :: C = -5363.80 ! K^2
    REAL, PARAMETER :: A =     18.3907
    REAL, PARAMETER :: B =   -372.090731 ! K
    REAL, PARAMETER :: C = -12350.598    ! K^2

    IF (temperature == 0.) THEN
       pfrost_n2 = 0.
    ELSE
       pfrost_n2 = EXP(A  + B/temperature + (C/temperature)/temperature)
    END IF
    IF (calculate_dpdt) THEN
       pfrost_n2 = pfrost_n2 * (-B/(temperature*temperature) &
                                - 2.*C/(temperature**3))
    ENDIF

  END FUNCTION pfrost_n2

!----------------------------------------------------------------
  REAL FUNCTION pfrost_ch4(temperature,calculate_dpdt)
!----------------------------------------------------------------
!
! Calculate the equilibrium vapor pressure of CH4 gas as a function
! of temperature for a range that is typical of the outer planets
! and satellites, unless the "calculate_dpdt" optional parameter
! is set to .TRUE., in which case the function returns the value
! of dP/dT (the Clausius-Clapeyron relation) of the equilibrium
! vapor pressure equation at the input temperature
!
! Input: temperature (K)
!
! Output: pressure (Pa)
!         OR
!         pressure/temperature (Pa/K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL,    INTENT(IN) :: temperature
    LOGICAL, INTENT(IN) :: calculate_dpdt

!  Parameter variables
    ! P = A e^(-B/(T+C))
    REAL, PARAMETER :: A =  133.3224*EXP(6.61184*LN10) ! Pa
    REAL, PARAMETER :: B = 389.93*LN10 ! K
    REAL, PARAMETER :: C = 266.-273.15 ! K

    IF (temperature == 0.) THEN
       pfrost_ch4 = 0.
    ELSE
       pfrost_ch4 = A * EXP(-B/(temperature+C))
    END IF
    IF (calculate_dpdt) THEN
       ! P = A e^(-B/(T+C))
       ! dP/dT = P * b / (T+C)^2
       pfrost_ch4 = pfrost_ch4 * B / ((temperature+C)*(temperature+C))
    ENDIF

  END FUNCTION pfrost_ch4

!----------------------------------------------------------------
  REAL(KIND(0d0)) FUNCTION pfrost_h2o_root(temperature,pressure)
!----------------------------------------------------------------
!
! Function to help invert the P_eq(T) equation for H2O when T < 273.16
! and a non-linear (i.e., non-polynomial) equation.
!
! Input: temperature (K)
!        pressure (Pa)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL(KIND(0d0)), INTENT(IN) :: temperature, pressure

    pfrost_h2o_root = (CC1*TC*LN10)/temperature - CC2*LOG(temperature) + &
                      (CC3*LN10/TC)*temperature +                        &
                      (-CC1*LN10 + CC2*LOG(TC) - CC3*LN10 - LOG(EIS) +   &
                       LOG(pressure))

  END FUNCTION pfrost_h2o_root

!----------------------------------------------------------------
  SUBROUTINE pfrost_h2o_root_deriv(temperature,pressure,f,df)
!----------------------------------------------------------------
!
! Function to help invert the P_eq(T) equation for H2O when T < 273.16
! and a non-linear (i.e., non-polynomial) equation.
!
! Input: temperature (K)
!        pressure (Pa)
! Output: f: Value of root-finding function with given inputs (ln[Pa])
!         df: Derivative of root-finding function with repsect to
!             temperature at the given input values (ln[Pa] K^-1)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL(KIND(0d0)), INTENT(IN   ) :: temperature, pressure
    REAL(KIND(0d0)), INTENT(  OUT) :: f, df

    f  = PFROST_H2O_ROOT(temperature, pressure)
    df = -(CC1*TC*LN10)/(temperature*temperature) - (CC2/temperature) + &
         CC3*LN10/TC
  END SUBROUTINE pfrost_h2o_root_deriv

!----------------------------------------------------------------
  SUBROUTINE pfrost_cooling_root(temperature,b,c,f,df)
!----------------------------------------------------------------
!
! Function to help solve for the energy balance when surface
! temperature must remain in equilibrium balance with surface
! pressure when cooling (due to radiation, etc.).  Equation
! to be solved is non-linear.
!
! Input:  temperature (K)
!         b, c: parameters used in defining equation
! Output: f:   Value of root-finding function with given inputs (K)
!         df:  Derivative of root-finding function with repsect to
!              temperature at the given input values (unitless, K/K)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL(KIND(0d0)), INTENT(IN   ) :: temperature, b, c
    REAL(KIND(0d0)), INTENT(  OUT) :: f, df

    f  = temperature + b*PFROST(REAL(temperature)) + c
    df = 1.d0 + b*PFROST(REAL(temperature),dpdt=.TRUE.) 

  END SUBROUTINE pfrost_cooling_root

!----------------------------------------------------------------
  REAL FUNCTION VISCOSITY(temperature)
!----------------------------------------------------------------
!
! Calculate the dynamic viscosity of air at a given temperature.
! Wrapper function for the gas appropriate to each planet.
!
! Input: temperature (K)
!
! Output: dynamic viscosity (Pa s == kg/m/s)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!  Input/Ouptut variables
    REAL, INTENT(IN) :: temperature

#if defined WRF_MARS
    viscosity = VISCOSITY_CO2(temperature)
#elif defined WRF_PLUTO
    viscosity = VISCOSITY_N2( temperature)
#elif defined WRF_TITAN
    ! Return flag value
    viscosity = 0.
#elif defined WRF_VENUS
    ! Return flag value
    viscosity = 0.
#elif defined WRF_EARTH
    ! Return flag value
    viscosity = 0.
#endif

  END FUNCTION VISCOSITY

!----------------------------------------------------------------
  REAL FUNCTION VISCOSITY_CO2(temperature)
!----------------------------------------------------------------
!
! Calculate the dynamic viscosity of CO2 gas
! Reference: Fenghour, A., W. A. Wakeham, and V. Vesovic, 
! "The Viscosity of Carbon Dioxide", J. Phys. Chem. Ref. Data, 27 (1),
! pp. 31--44, 1998
!
! Input:  temperature (K)
!
! Output: viscosity (Pa s == kg/m/s)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------
!  Input variables
    REAL, INTENT(IN) :: temperature

!  Local variables
    REAL, PARAMETER  :: c0 =  1.006970e-6 ! Pa s K^-(1/2)
    REAL, PARAMETER  :: a0 =  0.235156    ! unitless
    REAL, PARAMETER  :: a1 = -0.491266    ! unitless
    REAL, PARAMETER  :: a2 =  5.211155e-2 ! unitless
    REAL, PARAMETER  :: a3 =  5.347906e-2 ! unitless
    REAL, PARAMETER  :: a4 = -1.537102e-2 ! unitless
    REAL, PARAMETER  :: eovk = 251.196    ! K
    REAL             :: gstarn, lntstar

    lntstar = LOG(temperature/eovk)
    gstarn = EXP(a0 + lntstar*(a1 + lntstar*(a2 + lntstar*(a3 + lntstar*a4))))
    viscosity_co2 = c0*SQRT(temperature)/gstarn

    ! Alternate "Sutherland's formula":
    ! viscosity = K * T^(3/2) / (T + C)
    ! where K and C are constants characteristic of the gas

  END FUNCTION VISCOSITY_CO2

!----------------------------------------------------------------
  REAL FUNCTION VISCOSITY_N2(temperature)
!----------------------------------------------------------------
!
! Calculate the dynamic viscosity of N2 gas
!
! Input:  temperature (K)
!
! Output: viscosity (Pa s == kg/m/s)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------
!  Input variables
    REAL, INTENT(IN) :: temperature

!  Local variables
    REAL, PARAMETER  :: lambda0 = 1.406732195d-6 ! Pa s K^-2
    REAL, PARAMETER  :: C0      = 300.55         ! K

    viscosity_n2 = lambda0*temperature*SQRT(temperature)/(temperature+C0)

  END FUNCTION VISCOSITY_N2

!----------------------------------------------------------------
  REAL FUNCTION THERMAL_VELOCITY(temperature)
!----------------------------------------------------------------
!
! Calculate the mean thermal velocity of an air molecule
!
! Input:  temperature (K)
!
! Output: molecular speed (m/s)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------
!  Input variables
    REAL, INTENT(IN) :: temperature

!  Parameter variables
    REAL, PARAMETER :: term = 8.*k_boltzmann/(pi*mw_air*m_amu)

    !thermal_velocity = SQRT(8.*k_boltzmann*temperature/(pi*mw_air*m_amu))
    thermal_velocity = SQRT(term*temperature)

  END FUNCTION THERMAL_VELOCITY

!----------------------------------------------------------------
  REAL FUNCTION MEAN_FREE_PATH(temperature, pressure)
!----------------------------------------------------------------
!
! Calculate the mean free path of an air molecule at a given 
! temperature and pressure 
!
! Input:  temperature (K)
!         pressure    (Pa)
!
! Output: mean free path (m)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------
!  Input variables
    REAL, INTENT(IN) :: temperature
    REAL, INTENT(IN) :: pressure

!  Local variables
    REAL :: mu, v, rho

    mu  = VISCOSITY(temperature)        ! Pa s = kg/m/s
    v   = THERMAL_VELOCITY(temperature) ! m/s
    rho = pressure/(R_d*temperature)    ! kg/m^3
    mean_free_path = 2.*mu/(rho*v)

  END FUNCTION MEAN_FREE_PATH

!----------------------------------------------------------------
  REAL FUNCTION CUNNINGHAM_CORRECTION_FACTOR(Kn)
!----------------------------------------------------------------
!
! Calculate G, the Cunningham correction factor, to allow the
! "Stokes flow" equation to be applicable in the "slip flow" regime
!
! Input:  Knudsen number (unitless)
!
! Output: correction factor (unitless)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------
!  Input variables
    REAL, INTENT(IN) :: Kn

!  Local variables
    REAL, PARAMETER :: A_CCF = 1.257
    REAL, PARAMETER :: B_CCF = 0.400
    REAL, PARAMETER :: C_CCF = 1.100
    cunningham_correction_factor = 1. + Kn*(A_CCF+B_CCF*EXP(-C_CCF/Kn))

  END FUNCTION CUNNINGHAM_CORRECTION_FACTOR

!----------------------------------------------------------------
  REAL FUNCTION PARTICLE_FALL_VELOCITY(radius, temperature, pressure, rho_p)
!----------------------------------------------------------------
!
! Calculate the fall velocity of a particle of a given radius (e.g.,
! aerosol particle) through the air at a given temperature and pressure
!
! Input:  particle radius  (m)
!         temperature      (K)
!         pressure         (Pa)
!         particle density (kg/m^3)
!
! Output: fall velocity   (m/s)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------
!  Input variables
    REAL, INTENT(IN) :: radius
    REAL, INTENT(IN) :: temperature
    REAL, INTENT(IN) :: pressure
    REAL, INTENT(IN) :: rho_p

!  Local variables
    REAL :: rho_a, mu, Kn, Gfac, v

    rho_a         = pressure/(R_d*temperature)                    ! kg/m^3
    mu            = VISCOSITY(temperature)                        ! kg/m/s
    Kn            = MEAN_FREE_PATH(temperature,pressure)/radius   ! unitless
    Gfac          = CUNNINGHAM_CORRECTION_FACTOR(Kn)              ! unitless
    particle_fall_velocity = 2.*radius*radius*g*(rho_p-rho_a)*Gfac/(9.*mu) !m/s
    ! No upward velocities! (for now...)
    IF (particle_fall_velocity < 0.) particle_fall_velocity = 0.

  END FUNCTION PARTICLE_FALL_VELOCITY

!----------------------------------------------------------------
  REAL FUNCTION REYNOLDS_NUMBER_FALL_VELOCITY(radius, temperature, pressure, &
                                              rho_p)
!----------------------------------------------------------------
!
! Calculate the Reynolds number using the fall velocity of a
! particle of a given radius (e.g., aerosol particle) through
! the air at a given temperature and pressure
!
! Input:  particle radius (m)
!         temperature     (K)
!         pressure        (Pa)
!         particle density (kg/m^3)
!
! Output: Reynolds number (unitless)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------
!  Input variables
    REAL, INTENT(IN) :: radius
    REAL, INTENT(IN) :: temperature
    REAL, INTENT(IN) :: pressure
    REAL, INTENT(IN) :: rho_p

!  Local variables
    REAL :: v, mu, rho_a

    v     = PARTICLE_FALL_VELOCITY(radius, temperature, pressure, rho_p) ! m/s
    mu    = VISCOSITY(temperature)                                    ! kg/m/s
    rho_a = pressure/(R_d*temperature)                                ! kg/m^3

    reynolds_number_fall_velocity = 2.*radius*v*rho_a/mu            ! unitless

  END FUNCTION REYNOLDS_NUMBER_FALL_VELOCITY

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

      USE module_wrf_error

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) :: num_terms

      LOGICAL, INTENT(IN   ) :: l_double

      REAL(KIND(0.d0)), INTENT(   OUT) :: weights(2,num_terms)

! Local variables

      INTEGER :: i

      REAL(KIND(0.d0)) :: x1, x2

!c-mm      REAL(KIND(0.d0)) :: weight_block(num_terms)
      REAL(KIND(0.d0)), ALLOCATABLE :: weight_block(:,:)

! Initialize weights and set other initial values
      weights=0.0
      SELECT CASE (l_double)
         CASE (.FALSE.)
            x2=1.0
            x1=0.0
            ALLOCATE(weight_block(2,num_terms))
            CALL get_gaussian_weights(num_terms,x1,x2,weight_block)
            weights(:,1:num_terms)=weight_block
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
               END SELECT
               ALLOCATE(weight_block(2,num_terms/2))
               CALL get_gaussian_weights(num_terms/2,x1,x2,            &
                                         weight_block)
               weights(:,(i-1)*(num_terms/2)+1:i*num_terms/2)=           &
                                                          weight_block
               DEALLOCATE(weight_block)
            ENDDO
         CASE DEFAULT
            WRITE (wrf_err_message,*) &
               'Crashed in get_gaussian_weights'
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

      REAL(KIND(0.d0)), INTENT(   OUT) :: weights(2,num_terms)

! Local variables
      REAL(KIND(0.d0)), PARAMETER ::                  small = 5.0D-10, &
                                                      pi = 3.14159D+0
      INTEGER  :: m, j, k
      REAL(KIND(0.d0)) :: xm, xl, z, z1, p1, p2, p3, pp, rn
      m=INT((num_terms+1)/2)
      xm=(x2+x1)/2.
      xl=(x2-x1)/2.
      rn=REAL(num_terms)

      DO j=1,num_terms
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
         !write(*,*) j, num_terms, xl, z, z*xl+xm
         weights(1,num_terms+1-j)=2.0*xl/((1.0-z*z)*pp*pp)
         weights(2,num_terms+1-j)=z*xl+xm
         !weights(num_terms+1-j)=weights(j)

      END DO

      END SUBROUTINE get_gaussian_weights

END MODULE module_planet_utilities
