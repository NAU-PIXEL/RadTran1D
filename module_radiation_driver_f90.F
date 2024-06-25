MODULE module_radiation_driver
CONTAINS

SUBROUTINE radconst(secs_into_sol,DECLIN,SOLCON,JULIAN,           &
                       DEGRAD,DPD,                                   &
                       obliquity,eccentricity,semimajoraxis,         &
                       equinox_fraction,zero_date,                   &
                       fys_factor,                                   &
                       xlat,ha,P2SI,sza0,sunfrac,l_s,sunbody,        &
                       ids,ide, jds,jde, kds,kde,                    &
                       ims,ime, jms,jme, kms,kme,                    &
                       its,ite, jts,jte, kts,kte                     )
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   REAL, INTENT(IN   )      ::       DEGRAD,DPD,secs_into_sol,JULIAN!, gmt
!   INTEGER, intent(in)      ::       julday
   REAL, INTENT(OUT  )      ::       DECLIN,SOLCON
   REAL                     ::       OBECL,SINOB,SXLONG,ARG,  &
                                     DECDEG,DJUL,RJUL,ECCFAC

   INTEGER,  INTENT(IN   )  ::       ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte
   REAL, INTENT(  OUT)      ::       l_s,sunbody
   REAL, INTENT(IN   )      ::       P2SI, obliquity, eccentricity
   REAL, INTENT(IN   )      ::       semimajoraxis,equinox_fraction,zero_date
   REAL, INTENT(IN   )      ::       fys_factor
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: XLAT
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(  OUT) :: HA, sza0, sunfrac

!
! !DESCRIPTION:
! Compute terms used in radiation physics 
! This subroutine computes radius vector, declination and right ascension of
! sun, equation of time, and hour angle of sun at sunset, the mean
! cosine of the suns zenith angle and daylight fraction for n
! specified latitudes given the julian day and fraction.
!
! definition of the arguments used by this subroutine
!
! julday =         julian day
! radius =         radius vector (distance to sun in a. u.)
! declin =         declination of sun
! alp    =         right ascension of sun
! slag   =         apparent sun lag angle
! sza0   =         mean cosine zenith angle(do not divide by r**2).
! sunfrac=         fraction of daylight
! als    =         solar longitude
!
!     note - all angles are expressed in radians
!
!EOP

   REAL(KIND (0d0)), PARAMETER :: planet_year = 669.
   REAL(KIND (0d0)), PARAMETER :: solar_constant = 1361.5
   
   ! Local Variables
   REAL(KIND (0d0)) :: PI
   REAL(KIND (0d0)) :: DELEQN, DP_DATE

   REAL(KIND (0d0)) :: EM, W, ALS, RADIUS, SINDEC, PH, HA1, SS, CC
   REAL(KIND (0d0)) :: AR, AC, CD0, AP, EPS
   INTEGER I, J

   REAL(KIND (0d0)) :: EP, E, EQ, QQ, EN, ER, CR, SMALL_VALUE

   !REAL(KIND (0d0)) :: alp, slag, tini, tst



   PI=ACOS(-1.d0)
   !small_value = 1.d-7      ! caused problems on certain single precision
                             ! machines when it got stuck at: 4.7683716E-07
   SMALL_VALUE = 1.D-6
   DELEQN = EQUINOX_FRACTION * REAL(PLANET_YEAR)



   ! for short wave radiation

   DECLIN=0.
   SOLCON=0.
   l_s = 0.

   !-----OBECL : OBLIQUITY = 23.5 DEGREE.
        
   OBECL=OBLIQUITY*DEGRAD
   SINOB=SIN(OBECL)
        
   !-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:
        
   ! DP_DATE = DAYS SINCE LAST PERIHELION PASSAGE
   DP_DATE = JULIAN - ZERO_DATE
   DO WHILE (DP_DATE < 0.)
      DP_DATE=DP_DATE+REAL(PLANET_YEAR)
   END DO
   DO WHILE (DP_DATE > REAL(PLANET_YEAR))
      DP_DATE=DP_DATE-REAL(PLANET_YEAR)
   END DO

   ER = SQRT( (1.0+ECCENTRICITY)/(1.0-ECCENTRICITY) )
 
   ! qq is the mean anomaly
   QQ = 2.0 * (PI * DELEQN / REAL(PLANET_YEAR))

   ! Determine true anomaly at equinox: eq
   ! Iteration for eq
   E = 1.0
   CD0 = 1.
   DO WHILE (CD0 > SMALL_VALUE)
      EP = E - (E-ECCENTRICITY*SIN(E)-QQ) / (1.0-ECCENTRICITY*COS(E))
      CD0 = ABS(E-EP)
      E = EP
   END DO
   EQ = 2. * ATAN( ER * TAN(0.5*E) )

   ! Determine true anomaly at current date:  w
   ! Iteration for w
   EM = 2.0 * PI * DP_DATE / REAL(PLANET_YEAR)
   E = 1.0
   CD0 = 1.
   DO WHILE (CD0 > SMALL_VALUE)
      EP = E-(E-ECCENTRICITY*SIN(E)-EM) / (1.0-ECCENTRICITY*COS(E))
      CD0 = ABS(E-EP)
      E = EP
   END DO
   W = 2.0 * ATAN( ER * TAN(0.5*E) )
 
   ! Radius vector ( astronomical units:  AU )
   ALS= (W - EQ)/DEGRAD      !Aerocentric Longitude
   IF (ALS.LT.0.) ALS=ALS+360.
   IF (ALS.GE.360.) ALS=ALS-360.
   l_s = REAL(ALS)
   IF (l_s.LT.0.) l_s=l_s+360.
   IF (l_s.GE.360.) l_s=l_s-360.
   RADIUS = SEMIMAJORAXIS * (1 - ECCENTRICITY * COS(E))
   sunbody= radius

   ! Declination: declin
   ! sindec = sine of angle of inclination of planet's orbit
   SINDEC = SINOB * SIN(w-eq)
   DECLIN = ASIN(SINDEC)
   DECDEG=DECLIN/DEGRAD
   ! Right Ascension: alp
   ! Right ascension is only used for sun lag
   !tini = 0.d0
   !IF (obliquity > 0.) tini = 1.d0/TAN(obecl)
   !alp = ASIN(TAN(declin)*tini)
   !tst = COS(w-eq)
   !IF (tst < 0.d0) alp = pi  - alp
   !IF (alp < 0.d0) alp = alp + 2.d0*pi

   DO j = jts,jte
      DO i = its,ite
         PH = PI*XLAT(I,J)/180.

         IF (DECLIN.NE.0.0) THEN
            AP = ABS(PH)
            ! EPS is absolute angular distance from pole (radians)
            EPS = ABS(AP - 0.5 * PI)

            IF(EPS.LE.SMALL_VALUE) THEN
               ! i.e., if we are at a pole point
               ! Hour angle of sunset at the pole is either
               ! 0 or pi, depending on declination
               HA1 = 0.5 * PI * ABS(AP / PH + ABS(DECLIN) / DECLIN)
               SS = SIN(PH) * SIN(DECLIN)
               CC = 0.0
            ELSE
               ! we are not a pole point
               SS = SIN(PH) * SIN(DECLIN)
               CC = COS(PH) * COS(DECLIN)
               AR = - SS / CC
               AC = ABS(AR)

               IF((AC - 1. + SMALL_VALUE) .LT. 0.0) THEN
                  ! ABS(-SS/CC) was < 1. (can do ACOS)
                  HA1 = ACOS(AR)
               ELSE IF((AC - 1. + SMALL_VALUE) .EQ. 0.0) THEN
                  ! ABS(-SS/CC) is 0., result is +/- PI/2
                  HA1 = (AC - AR) * 0.5 * PI
               ELSE
                  ! ABS(-SS/CC) was > 1. (can't do ACOS)
                  IF(AR .LT. 0.) THEN
                     HA1 = PI
                  ELSE
                     HA1 = 0.0
                  END IF
               END IF
            END IF
         ELSE
            HA1 = 0.5 * PI
            SS = 0.0
            CC= COS(PH)
         END IF

         HA(I,J)=REAL(HA1)
         sunfrac(i,j) = REAL(HA1/PI)
         ! Integrating from -h to +h yields the average zenith angle
         ! sza is zenith angle averaged over daylight hours
         ! SIN(ha1)/ha1 is {integral of cos(theta) dtheta} / 
         !                 {integral of dtheta} ... ????
         if (ha1 == 0.d0) then
            sza0(i,j)=0.
         else
            sza0(i,j)=MAX(REAL(ss+cc*SIN(ha1)/ha1),0.)
         end if
      END DO
   END DO

   SOLCON=(solar_constant/(radius*radius))*fys_factor

   END SUBROUTINE radconst

 END MODULE module_radiation_driver
