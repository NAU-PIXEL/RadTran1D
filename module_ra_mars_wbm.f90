!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_mars_wbm

  USE module_ra_utils
  USE module_wrf_error
  USE module_ra_mars_common
  USE module_model_constants
  USE module_ra_valverde
  USE module_state_description, ONLY : VALVERDE_LOOKUP, &
                                       VALVERDE_FIT

CONTAINS

!==================================================================
  SUBROUTINE WBMVIS(RTHRATENSW,RTHRATEN,GSW,XLAT,                 &
                    XLONG,ALBEDO,HA,                              &
                    ANGSLOPE,AZMSLOPE,                            &
                    DUST_ARRAY,CLOUD_ARRAY,DHR,diurnavg,          &
                    rho_phy,T3D,QV3D,QC3D,QR3D,                   &
                    QI3D,QS3D,QG3D,P3D,PF3D,pi3D,dz8w,            &
                    R,CP,G,JULIAN,                                &
                    hr_vis,hr_a_vis,hr_g_vis,toasw,               &
                    sza0,sunfrac,                                 &
                    DECLIN,SOLCON,                                &
                    scm_is_global,                                &
                    ra_include_dust,                              &
                    du_physics,                                   &
                    sw_correct,                                   &
                    sfc_dir,sfc_diff,                             &
                    ids,ide, jds,jde, kds,kde,                    & 
                    ims,ime, jms,jme, kms,kme,                    &
                    its,ite, jts,jte, kts,kte                     ) 

!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------
    INTEGER,    INTENT(IN   ) ::        ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte
    INTEGER,    INTENT(IN   ) ::        du_physics
    LOGICAL,    INTENT(IN   ) ::        diurnavg
    LOGICAL,    INTENT(IN   ) ::        sw_correct
    LOGICAL,    INTENT(IN   ) ::        ra_include_dust, scm_is_global

    REAL, INTENT(IN    )      ::        JULIAN,DECLIN,SOLCON

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
          INTENT(IN    ) ::                                   P3D, &
                                                             PF3D, &
                                                             pi3D, &
                                                              T3D, &
                                                             QV3D, &
                                                             QC3D, &
                                                             QR3D, &
                                                             QI3D, &
                                                             QS3D, &
                                                             QG3D, &
                                                          rho_phy, &
                                                             dz8w, &
                                                       DUST_ARRAY, &
                                                      CLOUD_ARRAY

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
          INTENT(INOUT)  ::                              RTHRATEN
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
          INTENT(  OUT)  ::                            RTHRATENSW

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
          INTENT(  OUT)  ::                      hr_vis, hr_a_vis, &
                                                         hr_g_vis

    REAL, DIMENSION( ims:ime, jms:jme ),                           &
          INTENT(IN   )  ::                                  XLAT, &
                                                            XLONG, &
                                                           ALBEDO, &
                                                         ANGSLOPE, &
                                                         AZMSLOPE, &
                                                               HA, &
                                                             sza0, &
                                                          sunfrac

    REAL, DIMENSION( ims:ime, jms:jme ),                           &
          INTENT(INOUT)  ::                                   GSW

    REAL, DIMENSION( ims:ime, jms:jme ),                           &
          INTENT(  OUT)  ::                                 TOASW, &
                                                          sfc_dir, &
                                                         sfc_diff
          

    REAL, INTENT(IN   )   ::                          R,CP,G

    REAL(KIND(0d0)), INTENT(IN   )   ::                       DHR

! LOCAL VARS
 
    REAL, DIMENSION( ims:ime, kms:kme ) ::                         &
                                                           TTEN2D, &
                                                           RHO02D, &
                                                              P2D, &
                                                             PF2D, &
                                                           DUST2D, &
                                                          CLOUD2D, &
                                                             DZ2D, &
                                                              T2D, &
                                                             QV2D, &
                                                             QC2D, &
                                                             QR2D, &
                                                             QI2D, &
                                                             QS2D, &
                                                             QG2D, &
                                                         hr_vis2d, &
                                                       hr_a_vis2d, &
                                                       hr_g_vis2d

    REAL, DIMENSION( ims:ime ) ::                          XLAT1D, &
                                                          XLONG1D, &
                                                            ALB1D, &
                                                            GSW1D, &
                                                          TOASW1D, &
                                                             HA1D, &
                                                       ALPHASLP1D, &
                                                       GAMMASLP1D, &
                                                           SZA01D, &
                                                        SUNFRAC1D, &
                                                        sfc_dir1d, &
                                                       sfc_diff1d
    REAL :: JULDAYFRAC

    INTEGER :: I,J,K

!------------------------------------------------------------------
    j_loop: DO J=jts,MIN(jte,jde-1)

       DO i=its,MIN(ite,ide-1)
          DO k=kts,kte
             TTEN2D(i,k)=0.
             T2D(i,k)   = T3D(i,k,J)
             P2D(i,k)   = P3D(i,k,J)
             ! PF is pressure at the full levels (p8w in the regular code)
             PF2D(i,k)  = PF3D(i,k,J)
             RHO02D(i,k)= rho_phy(i,k,j)
             DZ2D(i,k)  = dz8w(i,k,J)
             DUST2D(i,k)= DUST_ARRAY(i,k,J)
             CLOUD2D(i,k)= CLOUD_ARRAY(i,k,J)
             !IF (P_QV >= P_FIRST_SCALAR) QV2D(i,k)=QV3D(i,k,J)
             !IF (P_QC >= P_FIRST_SCALAR) QC2D(i,k)=QC3D(i,k,J)
             !IF (P_QR >= P_FIRST_SCALAR) QR2D(i,k)=QR3D(i,k,J)
             !IF (P_QI >= P_FIRST_SCALAR) QI2D(i,k)=QI3D(i,k,J)
             !IF (P_QS >= P_FIRST_SCALAR) QS2D(i,k)=QS3D(i,k,J)
             !IF (P_QG >= P_FIRST_SCALAR) QG2D(i,k)=QG3D(i,k,J)
          END DO
          PF2D(i,kte+1)  = PF3D(i,kte+1,J)
          DUST2D(i,kte+1)= DUST_ARRAY(i,kte+1,j)
          XLAT1D(i)     = XLAT(i,J)
          XLONG1D(i)    = XLONG(i,J)
          ALB1D(i)      = ALBEDO(i,J)
          HA1D(i)       = HA(i,J)
          ALPHASLP1D(i) = ANGSLOPE(i,J)
          GAMMASLP1D(i) = AZMSLOPE(i,J)
          sza01D(i)     = sza0(i,j)
          sunfrac1D(i)  = sunfrac(i,j)
       END DO

       ! juldayfrac = fraction of a julian day (from midnight)

       JULDAYFRAC=JULIAN-INT(JULIAN)

       CALL WBMVIS2D(TTEN2D,GSW1D,XLAT1D,XLONG1D,ALB1D,           &
                     T2D,QV2D,QC2D,QR2D,QI2D,QS2D,QG2D,P2D,PF2D,  &
                     hr_vis2d,hr_a_vis2d, hr_g_vis2d, toasw1d,    &
                     diurnavg, sza01d, sunfrac1d,                 &
                     DUST2D,CLOUD2D,                              &
                     HA1D, ALPHASLP1D, GAMMASLP1D,                &
                     DHR, RHO02D, DZ2D,                           &
                     R,CP,G,DECLIN,SOLCON,JULDAYFRAC,             &
                     scm_is_global,                               &
                     ra_include_dust,                             &
                     du_physics,                                  &
                     sw_correct,                                  &
                     sfc_dir1d,sfc_diff1d,                        &
                     ims,ime,kms,kme,its,ite,kts,kte,ide,jde      )

       DO i=its,MIN(ite,ide-1)
          GSW(i,J) = GSW1D(i)
          TOASW(i,J) = TOASW1D(i)
          sfc_dir(i,J) = sfc_dir1d(i)
          sfc_diff(i,J) = sfc_diff1d(i)
          DO k=kts,kte
             RTHRATENSW(i,k,J) = TTEN2D(i,k)/pi3D(i,k,J)
             RTHRATEN(i,k,J)   = RTHRATEN(i,k,J) + RTHRATENSW(i,k,J)
             hr_vis(i,k,j)     = hr_vis2d(i,k)
             hr_a_vis(i,k,j)   = hr_a_vis2d(i,k)
             hr_g_vis(i,k,j)   = hr_g_vis2d(i,k)
          END DO
       END DO

    ENDDO j_loop

  END SUBROUTINE WBMVIS

!==================================================================
  SUBROUTINE WBMVIS2D(TTEN,GSW,XLAT,XLONG,ALBEDO,                  &
                      T,QV,QC,QR,QI,QS,QG,P,PF,                    &
                      hr_vis,hr_a_vis,hr_g_vis,toasw,              &
                      diurnavg, sza0, sunfrac,                     &
                      DUST_ARRAY,CLOUD_ARRAY,                      &
                      HA, ALPHASLP, GAMMASLP,                      &
                      DHR, RHO0, DZ,                               &
                      R,CP,G,DECLIN,SOLCON,JULDAYFRAC,             &
                      scm_is_global,                               &
                      ra_include_dust,                             &
                      du_physics,                                  &
                      sw_correct,                                  &
                      sfc_dir,sfc_diff,                            &
                      ims,ime,kms,kme,its,ite,kts,kte,ide,jde      )
!------------------------------------------------------------------
!
!     TO CALCULATE SHORT-WAVE ABSORPTION AND SCATTERING IN CLEAR
!     AIR AND REFLECTION AND ABSORPTION IN DUST LAYERS
!
!------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------

    INTEGER, INTENT(IN) :: ims,ime, kms,kme, its,ite, kts,kte, ide, jde
    INTEGER, INTENT(IN) :: du_physics

    REAL, DIMENSION( ims:ime, kms:kme ), INTENT(IN   ) ::          &
                                                             RHO0, &
                                                                T, &
                                                                P, &
                                                               PF, &
                                                               DZ, &
                                                               QV, &
                                                               QC, &
                                                               QR, &
                                                               QI, &
                                                               QS, &
                                                               QG, &
                                                       DUST_ARRAY, &
                                                      CLOUD_ARRAY

    REAL, DIMENSION( ims:ime, kms:kme ), INTENT(INOUT) ::    TTEN

    REAL, DIMENSION( ims:ime, kms:kme ), INTENT(  OUT) ::  hr_vis, &
                                                         hr_a_vis, &
                                                         hr_g_vis

    REAL, DIMENSION( ims:ime ), INTENT(IN   ) ::             xlat, &
                                                            xlong, &
                                                           albedo, &
                                                               ha, &
                                                         alphaslp, &
                                                         gammaslp, &
                                                             sza0, &
                                                          sunfrac
!                ALPHASLP = slope angle, GAMMASLP = slope azimuth

    REAL, DIMENSION( ims:ime ), INTENT(  OUT) ::              GSW, &
                                                            TOASW, &
                                                          sfc_dir, &
                                                         sfc_diff

    REAL, INTENT(IN  )   ::                         R,CP,G,DECLIN, &
                                                SOLCON,JULDAYFRAC

    REAL(KIND (0d0)), INTENT(IN) ::                           DHR

    LOGICAL, INTENT(IN   ) ::                            diurnavg
    LOGICAL, INTENT(IN   ) ::                          sw_correct
    LOGICAL, INTENT(IN   ) ::                     ra_include_dust
    LOGICAL, INTENT(IN   ) ::                       scm_is_global

!
! LOCAL VARS
!
    REAL, DIMENSION( kms:kme ) ::                            DLPI, &
                                                              HSW, &
                                                         HSW_DUST, &
                                                     dust_array1d

    REAL(KIND(0.d0)), DIMENSION( kms:kme ) ::      xltecorrection

    REAL(KIND(0.d0)), DIMENSION( ims:ime, kms:kme ) ::     p_nlte    ! A quasi-kludge to handle the 'KIND' of p

    REAL :: JDAYFRAC, RADIN, TRANS, MU0, COSFAC, RXLAT, TSUNFRAC, HRANG
    REAL :: PI, dTdt_0, p0ref, pNLTE, r0
    REAL :: radius, mubar

    INTEGER :: i, k, kk
!------------------------------------------------------------------

    PI     = ACOS(-1.)
    dTdt_0 = 1.3/(86400.*P2SI) ! K/s (SI seconds...)
    p0ref  = 700.              ! Pa
    pNLTE  = 0.0075            ! Pa
    r0     = 1.52              ! AU

    ! Instantaneous distance from the Sun to Mars, in AU
    !
    ! Rather than redo the integration of Kepler's equation, take
    ! advantage of the fact that we have as in input to this subroutine
    ! the local solar constant which is just the universal solar constant
    ! (1367.6 W/m^2 at 1 AU) divided by the radius squared.
    radius = SQRT(solar_constant/SOLCON) ! solar_constant now from module_model_consts

    i_loop: DO i = its, ite

       GSW(i)           = 0.
       DO k=kts,kte
          tten(i,k)     = 0.
          hr_vis(i,k)   = 0.
          hr_a_vis(i,k) = 0.
          hr_g_vis(i,k) = 0.
          HSW(k)        = 0.
          hsw_dust(k)   = 0.
       END DO

       ! mu0 = COS(zenith angle of the sun)
       ! tsunfrac = fraction of the timestep that the sun is up 

       ! for a global energy balance model, the product of the sunfrac and the cosine of the zenith angle
       ! should be 0.25, so as to capture the 1:4 ratio between the area of intercepted solar irradiation
       ! and the planetary surface area available for thermal emission
     IF(scm_is_global) THEN

         IF((ide == 2) .AND. (jde == 2)) THEN
           tsunfrac = 0.5
           mu0 = 0.5
!         ELSE
!           CALL wrf_error_fatal('You have scm_is_global set TRUE, but this is not an SCM!')
         ENDIF

     ELSE

       IF (diurnavg) THEN
          mu0 = sza0(i)
          tsunfrac = sunfrac(i)
       ELSE
          CALL ZENITH(JULDAYFRAC, DHR, DECLIN, XLAT(i), XLONG(i), HA(i), &
                      MU0, TSUNFRAC, HRANG)
       END IF

     END IF

       RADIN = SOLCON*TSUNFRAC

       TOASW(i) = MAX((RADIN*MU0), 0.)

       ! Sample values for testing code
       !  MU0=0.
       !  TSUNFRAC = HA / PI

       ! Don't do any of the shortwave code if the sun isn't up! (mu0 < 0)

       sun_up: IF( MU0 > 0.00001) THEN ! roundoff error

          mubar = SQRT((1224.*mu0*mu0+1.)/1225.)

          DO k=kts,kte
             DLPI(k)=PF(i,k)-PF(i,k+1)
             DLPI(k)=1./DLPI(k)
          END DO

          ! HSW is in K/s
          DO k = kts, kte
             HSW(k) = dTdt_0*((r0/radius)**2) * SQRT(mubar*p0ref/p(i,k)) / &
                      (1. + pNLTE/p(i,k))
          END DO

          ! Now adjust for slopes for surface

          ! If not considering slopes, use the following line and don't
          ! call the subroutine SWSLOPE
          !cosfac = mu0
          CALL SWSLOPE(cosfac, mu0, juldayfrac, declin, xlat(i), xlong(i), &
                       alphaslp(i), gammaslp(i))

          IF(ra_include_dust) THEN
          !dust_array1d(kts:kte+1) = dust_array(i,kte+1:kts:-1)
            DO k=kts,kte+1   ! note scatter requires z array to run in op dirn
               kk=kte+1-k+kts
               dust_array1d(kk) = dust_array(i,k)
            ENDDO

            CALL SW_AEROSOL_SCATTER(MU0,albedo(i),dust_array1d,HSW_DUST,TRANS, &
                                  sfc_dir(i),sfc_diff(i),kts,kte)
          ELSE

            HSW_DUST(kts:kte) = 0.

          ENDIF

          ! If dust is turned off for some reason, This definition of "trans"
          ! will make GSW be calculated properly, as F_solar*(1-albedo)
          IF ((du_physics == 0) .or. (.not.ra_include_dust)) trans=1.-albedo(i)

          GSW(i)=RADIN*COSFAC*TRANS

          IF (sw_correct) THEN
             p_nlte(i,:)=REAL(p(i,:),KIND(0.d0))
             CALL SWCORRECTION(p_nlte(i,:),xltecorrection,kms,kme,kts,kte)
             hsw(kts:kte)=hsw(kts:kte)/xltecorrection(kts:kte)
          ENDIF

          DO K=kts, kte
             ! Scatter has vertical levels in MM5 order,
             ! so HSW_DUST needs to be reversed
             kk=kte-k+kts
             hr_a_vis(i,k) = (G/CP)*RADIN*MU0*HSW_DUST(kk)*DLPI(K)
             hr_g_vis(i,k) = HSW(k)
             hr_vis(i,k)   = hr_a_vis(i,k) + hr_g_vis(i,k)
             TTEN(I,K) = hsw(k)+hr_a_vis(i,k)
          END DO

          !! Scatter has vertical in MM5 order,
          !! so HSW_DUST needs to be reversed
          !hr_vis(i,kts:kte)=HSW(kts:kte)
          !hr_a_vis(i,kts:kte)=(G/CP)*RADIN*MU0*HSW_DUST(kte:kts:-1)*DLPI(kts:kte)
          !TTEN(I,kts:kte)=HSW(kts:kte) + hr_a_vis(i,kts:kte)

       ENDIF sun_up ! mu0 > 0, sun is up 

    END DO i_loop

    RETURN

  END SUBROUTINE WBMVIS2D

!==================================================================
  SUBROUTINE WBMIR(RTHRATENLW,RTHRATEN,                           &
                   GLW,XLAT,XLONG,TSK,DUST_ARRAY,CLOUD_ARRAY,     &
                   rho_phy,T3D,QV3D,QC3D,QR3D,QI3D,QS3D,QG3D,     &
                   P3D,PF3D,pi3D,dz8w,R,CP,G,                     &
                   hr_ir,hr_a_ir,hr_g_ir,FNM,FNP,                 &
                   ra_include_dust,                               &
                   isothermal_top,nlte_physics,                   &
                   ids,ide, jds,jde, kds,kde,                     & 
                   ims,ime, jms,jme, kms,kme,                     &
                   its,ite, jts,jte, kts,kte                      ) 

!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------
    INTEGER,    INTENT(IN   ) ::        ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte
    LOGICAL,    INTENT(IN   ) ::        isothermal_top
    LOGICAL,    INTENT(IN   ) ::        ra_include_dust
    INTEGER,    INTENT(IN   ) ::        nlte_physics

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
          INTENT(IN    ) ::                                   P3D, &
                                                             PF3D, &
                                                             pi3D, &
                                                              T3D, &
                                                             QV3D, &
                                                             QC3D, &
                                                             QR3D, &
                                                             QI3D, &
                                                             QS3D, &
                                                             QG3D, &
                                                          rho_phy, &
                                                             dz8w, &
                                                       DUST_ARRAY, &
                                                      CLOUD_ARRAY

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
          INTENT(INOUT)  ::                              RTHRATEN

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
          INTENT(  OUT)  ::                            RTHRATENLW

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
          INTENT(  OUT)  ::                        hr_ir, hr_a_ir, &
                                                          hr_g_ir

    REAL, DIMENSION( ims:ime, jms:jme ),                           &
          INTENT(IN   )  ::                                  XLAT, &
                                                            XLONG, &
                                                              TSK

    REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: fnm, fnp

    REAL, DIMENSION( ims:ime, jms:jme ),                           &
          INTENT(INOUT)  ::                                   GLW

    REAL, INTENT(IN   )   ::                               R,CP,G

!
! LOCAL VARS
!

    REAL, DIMENSION( ims:ime, kms:kme ) ::                         &
                                                           TTEN2D, &
                                                           RHO02D, &
                                                              P2D, &
                                                             PF2D, &
                                                             DZ2D, &
                                                              T2D, &
                                                             QV2D, &
                                                             QC2D, &
                                                             QR2D, &
                                                             QI2D, &
                                                             QS2D, &
                                                             QG2D, &
                                                           DUST2D, &
                                                          hr_ir2d, &
                                                        hr_a_ir2d, &
                                                        hr_g_ir2d

    REAL, DIMENSION( ims:ime ) ::   TSK1D, GLW1D, XLAT1D, XLONG1D

    INTEGER :: I,J,K

!------------------------------------------------------------------

    j_loop: DO J=jts,jte

       DO i=its,ite
          DO k=kts,kte
             T2D(i,k)    = T3D(i,k,J)
             P2D(i,k)    = P3D(i,k,J)
             RHO02D(i,k) = rho_phy(i,k,J)
             DZ2D(i,k)   = dz8w(i,k,J)
             PF2D(i,k)   = PF3D(i,k,J)
             DUST2D(i,k) = DUST_ARRAY(i,k,J)
             !IF (P_QV >= P_FIRST_SCALAR) QV2D(i,k)=QV3D(i,k,J)
             !IF (P_QC >= P_FIRST_SCALAR) QC2D(i,k)=QC3D(i,k,J)
             !IF (P_QR >= P_FIRST_SCALAR) QR2D(i,k)=QR3D(i,k,J)
             !IF (P_QI >= P_FIRST_SCALAR) QI2D(i,k)=QI3D(i,k,J)
             !IF (P_QS >= P_FIRST_SCALAR) QS2D(i,k)=QS3D(i,k,J)
             !IF (P_QG >= P_FIRST_SCALAR) QG2D(i,k)=QG3D(i,k,J)
          END DO
          PF2D(i,kte+1)   = PF3D(i,kte+1,J)
          DUST2D(i,kte+1) = DUST_ARRAY(i,kte+1,J)
          XLAT1D(i)  = XLAT(i,J)
          XLONG1D(i) = XLONG(i,J)
          TSK1D(i)   = TSK(i,J)
       END DO

       CALL WBMIR2D(TTEN2D,GLW1D,TSK1D,XLAT1D,XLONG1D,          &
                    T2D,QV2D,QC2D,QR2D,QI2D,QS2D,QG2D,P2D,PF2D, &
                    DUST2D,CLOUD_ARRAY,fnm,fnp,hr_ir2d,         &
                    hr_a_ir2d, hr_g_ir2d,                       &
                    RHO02D,DZ2D,R,CP,G,isothermal_top,          &
                    nlte_physics, ra_include_dust,              &
                    ims,ime,kms,kme,its,ite,kts,kte             )

       DO i=its,ite
          DO k=kts,kte
             hr_g_ir(i,k,j)  = hr_ir2d(i,k)
             hr_a_ir(i,k,j)  = hr_a_ir2d(i,k)
             RTHRATENLW(i,k,J) = TTEN2D(i,k)/pi3D(i,k,J)
             RTHRATEN(i,k,J) = RTHRATEN(i,k,J) + RTHRATENLW(i,k,J)
          END DO
          GLW(i,j) = GLW1D(i)
       END DO

    ENDDO j_loop                         

  END SUBROUTINE WBMIR

!==================================================================
  SUBROUTINE WBMIR2D(TTEN,GLW,TG,XLAT,XLONG,                      &
                     T,QV,QC,QR,QI,QS,QG,P,PF,                    &
                     DUST_ARRAY,CLOUD_ARRAY,fnm,fnp,hr_ir,        &
                     hr_a_ir,hr_g_ir,                             &
                     RHO0,DZ,R,CP,G,isothermal_top,nlte_physics,  &
                     ra_include_dust,                             &
                     ims,ime,kms,kme,its,ite,kts,kte              )
!------------------------------------------------------------------
!
!     TO CALCULATE LONG WAVE RADIATION
!
!------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------
    INTEGER, INTENT(IN) :: ims,ime, kms,kme, its,ite, kts,kte
    LOGICAL, INTENT(IN) :: ra_include_dust
    LOGICAL, INTENT(IN) :: isothermal_top
    INTEGER, INTENT(IN) :: nlte_physics

    REAL, DIMENSION( ims:ime, kms:kme ), INTENT(IN   ) ::           &
                                                              RHO0, &
                                                                 T, &
                                                                 P, &
                                                                PF, &
                                                                DZ, &
                                                                QV, &
                                                                QC, &
                                                                QR, &
                                                                QI, &
                                                                QS, &
                                                                QG, &
                                                        DUST_ARRAY, &
                                                       CLOUD_ARRAY

    REAL, DIMENSION( ims:ime, kms:kme ), INTENT(INOUT) ::     TTEN

    REAL, DIMENSION( ims:ime, kms:kme ), INTENT(  OUT) ::    hr_ir, &
                                                           hr_a_ir, &
                                                           hr_g_ir

    REAL, DIMENSION( ims:ime ), INTENT(IN   ) :: XLAT, XLONG, TG

    REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: fnm, fnp

    REAL, DIMENSION( ims:ime ), INTENT(  OUT) :: GLW

    REAL, INTENT(IN  )   ::                                 R,CP,G

!
!   Local Variables
!
    REAL, DIMENSION( kms:kme ) ::                          HR_DUST, &
                                                            HR_CO2, &
                                                       HR_CO2_NLTE, &
                                                               T1D, &
                                                               P1D, &
                                                              PF1D, &
                                                              DLPI, &
                                                             DLPIT, &
                                                            DUST1D

    REAL, DIMENSION( kms:kme ) :: bsource, b_diff
    REAL, DIMENSION( kms:kme ) :: td1, tu1, fup, fdn, net_flux
    REAL, DIMENSION( kms:kme, kms:kme ) :: dep, td
    ! vertical variables in the NLTE scheme 
    REAL, DIMENSION( kms:kme ) :: z_lv, sigma_lv
    REAL :: bsource_surf, bsurf_diff, sflx_co2, sflx_dust

    REAL :: GOVCP, PI

    ! Altitude of the NLTE transition region in terms of -ln(p),
    ! where p is in nbars
    REAL, PARAMETER ::  zs_lv = -5.5
    ! Width of the NLTE transition region. The above value of zs_lv
    ! corresponds to approximately 82 km for reference atmosphere 1 in L.-V.
    REAL, PARAMETER ::  zw_lv =  0.5

    INTEGER :: i, k, kk
 
    GOVCP = G/CP
    hr_co2_nlte=0.0
    PI = ACOS(-1.0)
    
    i_loop: DO i = its, ite

       DO K=kts, kte
          TTEN(i,k) = 0.
          DLPI(K)   = PF(i,K)-PF(i,K+1)
          DLPI(K)   = 1./DLPI(K)
          p1d(k)    = p(i,k)
          pf1d(k)   = pf(i,k)
          t1d(k)    = t(i,k)
       END DO
       pf1d(kte+1)  = pf(i,kte+1)

       !TTEN(i,kts:kte) = 0.
       !DLPI(kts:kte)=PF(i,kts:kte)-PF(i,kts+1:kte+1)
       !DLPI(kts:kte)=1./DLPI(kts:kte)

       !t1d(kts:kte) = t(i,kts:kte)
       !The following made no difference...Preserved for posterity
       !DO k = kts+1, kte
       !   t1d(k) = fnm(k)*t(i,k) + fnp(k)*t(i,k-1)
       !END DO
       !t1d(kts)   = (1.+fnm(kts+1))*t(i,kts+1) + fnm(kts+1)*t(i,kts)
       !t1d(kte+1) = (1.+fnp(kte)  )*t(i,kte)   - fnp(kte)  *t(i,kte-1)

       DO k = kts, kte!+1
          !bsource(k) = PLANCK_FUNCTION(t1d(k))
          bsource(k+1) = PLANCK_FUNCTION(t1d(k))
       END DO
       bsource(kts)    = bsource(kts+1)
       bsource_surf    = PLANCK_FUNCTION(tg(i))

       DO k=kts,kte
          b_diff(k) = bsource(k+1) - bsource(k)
       END DO
       bsurf_diff   = bsource_surf - bsource(kts)

       CALL CO2_TRANS_HOURDIN(p1d, pf1d, t1d,    &
                              tu1, td1, td, g,   &
                              kts, kte, kms, kme )

       ! flux boundary conditions:
       fup(kts)   = bsource_surf
       fdn(kte+1) = 0.0
       ! The following should have worked, but didn't
       ! It produced extreme heating at the model top ...
       !IF (isothermal_top) fdn(kte+1) = -bsource(kte+1)

       DO k = kts+1, kte
         fup(k) = bsurf_diff    *tu1(k) + bsource(k)
         fdn(k) = bsource(kte+1)*td1(k) - bsource(k)
       END DO

       DO k = kts+1, kte
          DO kk = kts+2, k
             fup(k) = fup(k) - b_diff(kk-1)*td(kk-1,k)
          END DO
          fup(k)    = fup(k) - b_diff(k)*0.125*(td(k-1,k)+3.0)

          DO kk= k+2, kte+1
            fdn(k)  = fdn(k) - b_diff(kk-1)*td(kk-1,k)
          END DO
          fdn(k)    = fdn(k) - b_diff(k)*0.125*(td(k+1,k)+3.0)
       END DO

       ! The following is not a perfect solution to the above problem,
       ! But at least it works, and keeps things close to what the rate
       ! seems like it should be in this situation
       IF (isothermal_top) fdn(kte+1) = fdn(kte)

       ! Upward flux at the top boundary:
       fup(kte+1) =  bsurf_diff*tu1(kte+1) + bsource(kte+1)
       DO k = kts+1, kte
          fup(kte+1) = fup(kte+1) - b_diff(k)*td(k,kte+1)
       END DO

       ! Downward flux at the bottom boundary
       fdn(kts) = bsource(kte+1)*td1(kts) - bsource(kts+1)
       DO k = kts+1, kte
          fdn(kts) = fdn(kts) - b_diff(k)*td(k,kts)
       END DO

       DO k=kts,kte+1
          net_flux(k) = fup(k) + fdn(k)
       END DO

       DO k = kts, kte
          hr_co2(k)=govcp*(net_flux(k)-net_flux(k+1))*dlpi(k)
       ENDDO

       sflx_co2= -fdn(kts)

       ! ---------------------------------------------------------------------
       ! Non Local Thermodynamic Equilibrium (NLTE) adjustment to LW flux

       SELECT CASE(nlte_physics)
          ! calculate NLTE CO2 cooling rates between 40 and 140 km
          ! using look up tables
          CASE (VALVERDE_LOOKUP)
             CALL VALVERDE_CO2_IR_LOOKUP(hr_co2_nlte,pf1d,t1d,kms,kme,kts,kte)
          ! calculate NLTE CO2 cooling rates between 40 and 140 km
          ! using analytical fit to table data
          CASE(VALVERDE_FIT)
             CALL VALVERDE_CO2_IR_FIT(hr_co2_nlte,pf1d,t1d,kms,kme,kts,kte)
          CASE DEFAULT
             !WRITE(0,*) 'NLTE_PHYSICS option not found'
             DO k = kts, kte
                hr_co2_nlte(k) = 0.
             END DO
       END SELECT

       ! Combine the LTE and NLTE cooling rates to yield a contnuous profile.
       ! LTE is valid from 0 to about 60 km, above that NLTE becomes important.
       IF (nlte_physics /= 0) THEN
       DO k=kts,kte
          ! z_lv is the vertical variable in the original paper
          ! The expression in parenthesis has to be in nanobars
          z_lv(k) = -log(pf1d(k)*pa2bar*1E9)
          ! The blending coefficient (a function of height)
          sigma_lv(k) = 0.5*(1.+tanh((z_lv(k)-zs_lv)/zw_lv))
          ! Only blend heating rates if there is a non-zero NLTE heating rate
          IF ( ABS(hr_co2_nlte(k)) > 0.) THEN
             hr_co2(k) = hr_co2(k)*(1.-sigma_lv(k)) + &
                         hr_co2_nlte(k)*sigma_lv(k)
          END IF
       END DO   
       ENDIF
       ! --------------------------------------------------------------------

       ! Dust goes here
       IF(ra_include_dust) THEN

       ! LW_AEROSOL_HEAT has vertical indices in MM5 order
       DO k=kts,kte
          kk=kte-k+kts
          t1d(kk)   = t(i,k)
          dlpit(kk) = dlpi(k)
       END DO
       DO k=kts,kte+1
          kk=kte+1-k+kts
          dust1d(kk) = dust_array(i,k)
       END DO
       CALL LW_AEROSOL_HEAT(dlpit, t1d, TG(i), dust1d, HR_DUST, SFLX_DUST, &
                            G, CP, kts, kte)
       ELSE

         HR_DUST(kts:kte)=0.
         SFLX_DUST=0.

       ENDIF

       DO k=kts,kte
          kk=kte-k+kts
          TTEN(i,k)    = HR_CO2(k) + HR_DUST(kk)
          hr_g_ir(i,k) = HR_CO2(k)
          hr_a_ir(i,k) = HR_DUST(kk)
          hr_ir(i,k)   = HR_DUST(kk) + HR_CO2(k)
       END DO
       GLW(i)         = SFLX_CO2 + SFLX_DUST

    END DO i_loop

    RETURN
  END SUBROUTINE WBMIR2D

!====================================================================
  REAL FUNCTION PLANCK_FUNCTION(t) RESULT (planck)
!--------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------
    REAL, INTENT(IN) :: t

    REAL, PARAMETER :: h = 6.62607e-34 ! J-sec
    REAL, PARAMETER :: k = 1.38065e-23 ! J/K
    REAL, PARAMETER :: c = 2.99792e+08 ! m/s

    REAL :: nu, nu3, dnu, c1, arg
    REAL :: pi

    ! result will have units of W/m^2

    pi = ACOS(-1.)

    nu  = 64500.        ! m^-1
    nu3 = nu*nu*nu      ! m^-3
    c1  = 2.*pi*h*c*c   ! J m^2 s^-1
    arg = h*c*nu/k      ! K
    dnu = 86500.-50000. ! m^-1

    planck = c1*nu3*dnu/(EXP(arg/t)-1.)

  END FUNCTION PLANCK_FUNCTION

!====================================================================
  SUBROUTINE co2_trans_hourdin( pp, pph, t, tu1, td1, td, g, &
                                kts, kte, kms, kme)
!--------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------
!   pph = pressure at full sigma levels (i.e., borders of vertical levels)
!   pp  = pressure at half sigma levels (i.e., centers of vertical levels)

    INTEGER, INTENT(IN) :: kts, kte, kms, kme
    REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: pp, pph, t
    REAL, DIMENSION( kms:kme ), INTENT(  OUT) :: tu1, td1
    REAL, DIMENSION( kms:kme, kms:kme ), INTENT(OUT) :: td
    REAL, INTENT(IN) :: g

    REAL, PARAMETER :: tzero   = 200.    ! K
    REAL, PARAMETER :: diffrac = 1.66    ! Unitless
    REAL, PARAMETER :: p0ref   = 1.013e5 ! Pa
    REAL, PARAMETER :: bwt1    = 135./365. ! Fractional part of band
    REAL, PARAMETER :: bwt2    =  70./365. ! Fractional part of band
    REAL, PARAMETER :: bwt3    = 160./365. ! Fractional part of band
    ! c1 and c2 have units!!!!  See below for proper use/definition
    REAL, PARAMETER :: c1c = 0.005
    REAL, PARAMETER :: c2c = 0.01
    REAL, PARAMETER :: c1w = 0.015
    REAL, PARAMETER :: c2w = 0.1
    REAL, PARAMETER, DIMENSION(4) ::                                 &
         acent = (/ 0.694E-1,  0.328E-3, 0.544E-2,  0.596E-5 /),     &
         awing = (/ 0.275E-1, -0.705E-3, 0.356E-1, -0.326E-4 /)
    REAL, PARAMETER, DIMENSION(5) ::                                 &
         pcoef   = (/ 2.882E-5, 1.708E-2, -3.397E-2, 1.454E-2, 5.438E-1 /), &
         pdcoef  = (/ 2.924E-5, 1.732E-2, -3.442E-2, 1.475E-2, 5.511E-1 /), &
         pcoefw  = (/ 2.893E-2, 1.906, 3.841, 1.895, 6.004 /), &
         pdcoefw = (/ 2.832E-2, 1.893, 3.792, 1.881, 5.977 /)

    INTEGER :: k, kk
    REAL, DIMENSION( kms:kme ) :: delp, ubar, upbar, ubarw, upbarw
    REAL :: factor_1, factor_2, exf1, exf2, ueq, ueqw
    REAL :: udiff, updiff, updiffw, udiffw, tdw, tdc
    REAL :: deltat

    factor_1= diffrac/g
    factor_2= factor_1/p0ref

    DO k=kts,kte
       delp(k) = pph(k) - pph(k+1)
    END DO
      
    ! Integrate downwards to obtain path lengths
    ubar  (kte+1) = 0.
    upbar (kte+1) = 0.
    ubarw (kte+1) = 0.
    upbarw(kte+1) = 0.

    DO k= kte, kts, -1
       deltat = t(k) - tzero

       ! Calculate path lengths in the band center:
       exf1 = EXP( deltat*( acent(1) + acent(2)*deltat ) )
       exf2 = EXP( deltat*( acent(3) + acent(4)*deltat ) )

       ubar(k)  = ubar(k+1)  + factor_1*exf1 * delp(k)
       upbar(k) = upbar(k+1) + factor_2*exf2 * delp(k)*pp(k)

       ! Repeat for the wings:                 
       exf1 = EXP( deltat*( awing(1) + awing(2)*deltat ) )
       exf2 = EXP( deltat*( awing(3) + awing(4)*deltat ) )

       ubarw(k)  = ubarw(k+1)  + factor_1*exf1 * delp(k)
       upbarw(k) = upbarw(k+1) + factor_2*exf2 * delp(k)*pp(k)

    END DO

    td1(kte+1) = 1.
    tu1(kts)   = 1.

    ! Downward paths:
    ! use top temperature for pade coeffs
    DO k= kts, kte
       ! equivalent absorber amount:
       ! The c1, c2 constants are valid when u[p]bar[w] are in g/cm^2 and 
       ! our u[p]bar[w] are in kg/m^2. Conversion factor: 1 g/cm^2 = 10 kg/m^2
       ueq  = SQRT( 0.1* upbar(k) ) + c1c*( (0.1*ubar (k))**c2c )
       ueqw = SQRT( 0.1*upbarw(k) ) + c1w*( (0.1*ubarw(k))**c2w )

       ! Pade approximants:  
       tdc = PADE_EVAL(pcoef,  ueq )
       tdw = PADE_EVAL(pcoefw, ueqw)

       ! Finally, combine with Planck function weighting:
       td1(k)= tdc*bwt2 + tdw*(bwt1+bwt3)
    END DO

    ! Upward paths:
    ! Use surface temperature for pade coeffs
    DO k = kts+1, kte+1
       updiff = upbar(kts) - upbar(k)
       udiff  = ubar(kts)  - ubar(k)
       ! The c1, c2 constants are valid when u[p]bar[w] are in g/cm^2 and 
       ! our u[p]bar[w] are in kg/m^2. Conversion factor: 1 g/cm^2 = 10 kg/m^2
       ueq    = SQRT( 0.1*updiff ) + c1c*((0.1*udiff)**c2c)
       ! equivalent absorber amount:
       updiff = upbarw(kts) - upbarw(k)
       udiff  = ubarw(kts)  - ubarw(k)
       ! The c1, c2 constants are valid when u[p]bar[w] are in g/cm^2 and 
       ! our u[p]bar[w] are in kg/m^2. Conversion factor: 1 g/cm^2 = 10 kg/m^2
       ueqw   = SQRT( 0.1*updiff ) + c1w*((0.1*udiff)**c2w)

       ! Pade approximants:
       tdc = PADE_EVAL( pcoef , ueq  )
       tdw = PADE_EVAL( pcoefw, ueqw )

       ! Finally, combine with Planck function weighting:
       tu1(k)= tdc*bwt2 + tdw*( bwt1+bwt3 )
    END DO

    ! ----------------  Differential Paths  ----------------:

    DO k = kts, kte+1
       td(k,k)= 1.0
    END DO

    DO k = kts, kte+1
       DO kk = k+1, kte+1
          ! Center:
          udiff  = ubar(k)  -  ubar(kk)
          updiff = upbar(k) - upbar(kk)
          ! The c1, c2 constants are valid when u[p]bar[w] are in g/cm^2 and 
          ! our u[p]bar[w] are in kg/m^2. Conversion factor: 1 g/cm^2=10 kg/m^2
          ueq  = SQRT( 0.1*updiff ) + c1c*((0.1*udiff)**c2c)
          ! Contribution from the Wings:  
          udiffw  =  ubarw(k) -  ubarw(kk)
          updiffw = upbarw(k) - upbarw(kk) 
          ! The c1, c2 constants are valid when u[p]bar[w] are in g/cm^2 and 
          ! our u[p]bar[w] are in kg/m^2. Conversion factor: 1 g/cm^2=10 kg/m^2
          ueqw = SQRT( 0.1*updiffw )+ c1w*((0.1*udiffw)**c2w)

          ! Pade approximants:  
          tdc = PADE_EVAL( pdcoef,  ueq )
          tdw = PADE_EVAL( pdcoefw, ueqw)

          ! Planck function weighting:  
          td(k,kk)= tdc*bwt2 + tdw*(bwt1+bwt3)
          ! td is symmetric
          td(kk,k)= td(k,kk)
       END DO
    END DO

  END SUBROUTINE CO2_TRANS_HOURDIN

!====================================================================
  REAL FUNCTION pade_eval( pc, ueq ) RESULT (td)
!--------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------
    REAL, DIMENSION(5), INTENT(IN) :: pc
    REAL, INTENT(IN) :: ueq

    REAL :: top, denom

    ! Pade approximants:  
    top   = pc(1) + ueq*( pc(2) + ueq*(pc(3)      ) )
    denom = pc(1) + ueq*( pc(4) + ueq*(pc(5) + ueq) )
    td    = top / denom

  END FUNCTION pade_eval

!==================================================================
  SUBROUTINE wbmlwinit(RTHRATEN, RTHRATENLW, restart,&
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte )
!--------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------
    LOGICAL , INTENT(IN)           :: restart
    INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,  &
                                      ims, ime, jms, jme, kms, kme,  &
                                      its, ite, jts, jte, kts, kte

    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: &
                                                           RTHRATEN, &
                                                         RTHRATENLW
    INTEGER :: i, j, k, itf, jtf, ktf

    jtf=MIN(jte,jde-1)
    ktf=MIN(kte,kde-1)
    itf=MIN(ite,ide-1)

    write(wrf_err_message,*) '--- ---> Using Mars Wide Band rad tran for LW'
    call wrf_message(trim(wrf_err_message))

    IF(.not.restart)THEN
       DO j=jts,jtf
       DO k=kts,ktf
       DO i=its,itf
          RTHRATEN(i,k,j)=0.
          RTHRATENLW(i,k,j)=0.
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE wbmlwinit

!==================================================================
  SUBROUTINE wbmswinit(RTHRATEN, RTHRATENSW, restart, &
                       ids, ide, jds, jde, kds, kde,  &
                       ims, ime, jms, jme, kms, kme,  &
                       its, ite, jts, jte, kts, kte  )
!--------------------------------------------------------------------
    IMPLICIT NONE
!--------------------------------------------------------------------
    LOGICAL , INTENT(IN)           :: restart
    INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,  &
                                      ims, ime, jms, jme, kms, kme,  &
                                      its, ite, jts, jte, kts, kte

    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: &
                                                           RTHRATEN, &
                                                         RTHRATENSW
    INTEGER :: i, j, k, itf, jtf, ktf

    jtf=MIN(jte,jde-1)
    ktf=MIN(kte,kde-1)
    itf=MIN(ite,ide-1)

    write(wrf_err_message,*) '--- ---> Using Mars Wide Band rad tran for SW'
    call wrf_message(trim(wrf_err_message))

    IF(.not.restart)THEN
       DO j=jts,jtf
       DO k=kts,ktf
       DO i=its,itf
          RTHRATEN(i,k,j)=0.
          RTHRATENSW(i,k,j)=0.
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE wbmswinit

!==================================================================

END MODULE module_ra_mars_wbm
