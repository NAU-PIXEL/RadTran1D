!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_utils

CONTAINS

!====================================================================
  SUBROUTINE SWSLOPE( COSFAC, MU0, JULDAYFRAC, DECLIN, XLAT, XLONG, &
                      ALPHASLP, GAMMASLP )
!--------------------------------------------------------------------
!   The subroutine calculates the "cosine factor" or reduction in
!   incident solar radiation due to the presence of a slope.
!   N.B. This subroutine only accounts for the difference in vector
!   to the sun and the normal vector of the surface element.
!   It does **NOT** account for shadowing by slopes.  That is handled
!   in a different set of WRF subroutines.
!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------

    REAL , INTENT(IN   ) :: mu0, xlat, xlong, alphaslp, gammaslp
    REAL , INTENT(IN   ) :: juldayfrac, declin
    REAL , INTENT(  OUT) :: cosfac

    REAL :: zsin, jdayfrac, hra, rxlat
    REAL :: betaazm, pi

    pi = ACOS(-1.0)

    ! ZSIN is the SIN of the solar zenith angle, always positive.
    ! MU0 is the COS of the solar zenith angle, so use an identity.
    ! The formula from Pielke (2002) for COS Z represents the instantaneous
    ! value of MU0.  The subroutine ZENITH appropriately handles
    ! sunrises and sunsets and averaging over the radiation timeperiod
    ! to get the "effective" MU0.  So it is best to use the above MU0
    ! to calculate SIN Z.

    ZSIN = SQRT(1.-MU0*MU0)

    ! Check for roundoff error

    IF (MU0*MU0 > 1.) ZSIN=0.

    ! Need to recalculate fractional Julian day to account for longitude
    ! For HRA, 0 is local noon (JDAYFRAC is 0.5 at noon)

    JDAYFRAC=JULDAYFRAC+XLONG/360.
    IF (JDAYFRAC.LT.0.) JDAYFRAC=JDAYFRAC+1.
    JDAYFRAC=JDAYFRAC-INT(JDAYFRAC)
    HRA = JDAYFRAC-0.5
    IF (HRA < 0.) HRA = HRA + 1.
    HRA = HRA * 2 * PI

    RXLAT=PI*XLAT/180.

    ! BETAAZM is the azimuth of the sun.
    ! Can not use the following formula directly in case of roundoff error.
    ! See fix below
    ! BETAAZM = ASIN(COS(DECLIN)*SIN(HRA)/ZSIN)
    BETAAZM=COS(DECLIN)*SIN(HRA)/ZSIN
    ! Need to check for roundoff error since we will take an ARC SIN
    IF (ABS(BETAAZM) > 1.) BETAAZM=SIGN(1.,BETAAZM)
    BETAAZM=ASIN(BETAAZM)
    IF (RXLAT < DECLIN) BETAAZM = PI - BETAAZM

    ! Calculate the factor as per Pielke (2002), p. 405
    ! Reverts to MU0(I) if ALPHASLP is 0 (no slope)

    COSFAC = COS(ALPHASLP)*MU0+SIN(ALPHASLP)*ZSIN*COS(BETAAZM-GAMMASLP)

    ! Correct for unreal values
    IF (COSFAC < 0.) COSFAC = 0.
    IF (COSFAC > 1.) COSFAC = 1.

    RETURN

  END SUBROUTINE SWSLOPE

!====================================================================
  SUBROUTINE ZENITH( FJD, DHR, DECLIN, XLAT, XLONG, HA, COSZ, FRAC, HRANG )
!------------------------------------------------------------------
!
!   zenith computes effective mean cosine of zenith angle and daylight
!   fraction from latitude and parameters given by subroutine solar
!   by means of the following equation -
!
!      cosz = sin(xlat)*sin(dlt) + cos(xlat)*cos(dlt)*sin(arg)/arg
!
!                                  ^ missing a *cos(hloc)*
!
!     definition of variables appearing in this subroutine
!
!     fjd   = day fraction beyond integral julian day(begins at noon ut)
!     dlt   = declination of sun
!     xlat  = latitude
!     ha    = hour angle of sun at sunset
!     dhr   = half the integration (or averaging) period in radians
!     nlng  = number of gridpoints in a row
!     cosz  = mean cosine of zenith angle for all longitudes
!     frac  = daylight fraction at all longitudes
!     cvpr  = difference in longitude between gridpoints
!     gha   = hour angle of sun at greenwich (west of meridian is plus)
!     arg   = half of integration period in radians
!     xlng  = longitude of gridpoint
!     hloc  = hour angle of sun plus arg at longitude xlng
!     delsh = truncated averaging period resulting from late sunrise or
!               early sunset (dele and delw are similarly truncated)
!
!   input arguments to circular functions are in radians
!
!------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------

    REAL, INTENT(IN) ::  FJD, XLAT, XLONG, HA, DECLIN
    REAL(KIND(0d0)), INTENT(IN) :: DHR
    REAL, INTENT(OUT) :: COSZ, FRAC, HRANG

    ! Local variables
    LOGICAL :: RISE, SET
    INTEGER :: I
    REAL(KIND(0d0)) :: PI, TPI, DLT, GHA, ARG, SINFAC, RXLAT
    REAL(KIND(0d0)) :: SS, CC, CONS, HLOC, HLPAR, ARMHL
    REAL(KIND(0d0)) :: DELE, DELW, DELSH

    PI  = 4.0*ATAN(1.d0)
    TPI = 2.0*PI
    DLT = DECLIN

    ! compute basic constants used by this subroutine

    ! Our Julian days begin at midnight so we need to add pi to be
    ! consistent with the above definition of gha
    GHA    = FJD * TPI + PI
    ARG    = DHR
    ! This term apparently represents an averaging over the "radiation"
    ! timeperiod (dhr), allowing us to get an "average" mu0 for use in
    ! the radiation schemes.
    SINFAC = SIN(ARG) / ARG

    RXLAT=PI*XLAT/180.
    SS     = SIN(RXLAT) * SIN(DLT)
    CC     = COS(RXLAT) * COS(DLT)
    IF (HA > 0.0) THEN
       CONS = SS + CC * SIN(HA) / HA
    END IF

    ! compute cosz and frac at all longitudes

    HLOC = GHA + PI*XLONG/180.
    RISE = .FALSE.
    SET  = .FALSE.

    ! reset hloc to within plus and minus pi
    HLOC = MOD(HLOC,TPI)
    IF (HLOC .GT. PI) HLOC = HLOC - TPI

    ! determine if sun rises or sets during averaging period

    HLPAR = HLOC + ARG
    ARMHL = ARG - HLOC
    IF (HLPAR > HA) SET  = .TRUE.
    IF (ARMHL > HA) RISE = .TRUE.
    IF (RISE .AND. SET) THEN
       IF (HA <= 0.0) THEN
          COSZ = 0.0
          FRAC = 0.0
       ELSE

          ! averaging period covers the entire duration of daylight

          COSZ = CONS
          FRAC = HA / ARG
       END IF
    ELSE
       IF (HLPAR > PI) THEN
          DELE = 0.5d0 * MAX(HLPAR + REAL(HA,KIND(0D0)) - TPI, 0.d0)
          DELW = 0.5d0 * MAX(REAL(HA,KIND(0D0)) + ARMHL, 0.d0)
       ELSE IF (ARMHL > PI) THEN
          DELE = 0.5d0 * MAX(REAL(HA,KIND(0D0)) + HLPAR, 0.d0)
          DELW = 0.5d0 * MAX(ARMHL + REAL(HA,KIND(0D0)) - TPI, 0.d0)
       END IF
       IF ( (HLPAR > PI) .OR. (ARMHL > PI) ) THEN

          ! averaging interrupted by both sunset and sunrise - done in 2 parts

          FRAC = (DELE + DELW) / ARG
          IF (FRAC == 0.0) THEN
             COSZ = 0.0
          ELSE
             COSZ = SS + CC*(COS(HA - DELE)*SIN(DELE) &
                       + COS(HA - DELW) * SIN(DELW))  &
                         / (DELE + DELW)
          END IF
       END IF
    END IF
    IF ( (.NOT.(RISE.AND.SET)) .AND. (HLPAR <= PI) .AND. &
         (ARMHL.LE.PI) ) THEN
       IF (SET)  DELSH = 0.5 * (HA + ARMHL)
       IF (RISE) DELSH = 0.5 * (HA + HLPAR)

       ! either a sunrise or sunset occurs during the averaging period

       IF (RISE .OR. SET) THEN
          IF (DELSH <= 0.0) THEN
             COSZ = 0.0
             FRAC = 0.0
          ELSE
             COSZ = SS + CC * COS(HA -DELSH) * SIN(DELSH) / DELSH
             FRAC = DELSH / ARG
          END IF
       ELSE

          ! normal case where sun is up during entire averaging period

          COSZ = SS + CC * COS(HLOC) * SINFAC
          FRAC = 1.0
       END IF
    END IF

    ! limit cosz and fraction to within range of zero to one

    COSZ = MIN( 1.0, COSZ )
    COSZ = MAX( 0.0, COSZ )
    FRAC = MIN( 1.0, FRAC )
	HRANG = HLOC
	
    RETURN

  END SUBROUTINE ZENITH

END MODULE module_ra_utils
