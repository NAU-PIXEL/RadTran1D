!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_valverde

  USE module_array_utilities

! This is a WRF adaption of the Miguel A. Lopez Valverde of Instituto de 
! Astrofisica de Andalucia (CSIC), Granada, Spain code delivered for the
! "Martian Environment Models" project (ESA contract 11369/95/nl/jg CCN2).  
! For details see the contract Final technical report and Lopez-Valverde 
! simplified 2 level model (Lopez-Valverde, Lopez-Puertas, "A fast computation
! of radiative heating rates under non-LTE in a CO2 atmosphere", IRS 2000,
! Current Problems in Atmospheric Radiation,Eds. Smith and Timofeev, 2001)
!
! There are 3 ways to use this scheme: 
! VALVERDE_CO2_IR_LOOKUP, VALVERDE_CO2_IR_FIT, or VALVERDE_CO2_IR_TEST.
! VALVERDE_CO2_IR_LOOKUP returns CO2 heating rates in the Non-Local Thermal
! Equilibrium (NLTE) regime using tables with mixing ratios and escape function
! data provided with the code.  This approach increases computation time by a
! factor of 8.
! Alternatively, the heating rates can be calculated with VALVERDE_CO2_IR_FIT,
! which fits the data in the tables with analytical functions. There are small
! differences between the tabulated and fitted data, but they should be
! unimportant given the approximate nature of the scheme to start with. 
! Finally, VALVERDE_CO2_IR_TEST uses tabulated data and temeparture and pressure
! profiles from the reference atmosphere marte_ref1.dat for comparison of 
! cooling rates with published results (Fig. 7 of the Final technical report).
! This call returns zero coling rates but prints the test cooling rates to
! wrf.log file.
!
! The heating rate is calculated in terms of ergs/cm3/s in the original code.

  IMPLICIT NONE

    ! Ensure that variables declared here don't interfere with names of
    ! variables that USE this module by making all variables PRIVATE
    PRIVATE
    ! But make sure the subroutines names *are* PUBLICally available...
    PUBLIC :: valverde_co2_ir_lookup, valverde_co2_ir_fit, valverdeinit
    PUBLIC :: swcorrection


  INTEGER, PARAMETER :: np1 = 110        ! # data points in Table#1
  INTEGER, PARAMETER :: np2 = 114        ! # data points in Table#2
  INTEGER, PARAMETER :: nl=151           ! number of atmospheric layers in the reference atmosphere file

  REAL, DIMENSION (nl) :: co2vmr,o3pvmr,n2vmr,covmr, & ! vol.mixing ratios from the reference file
                                                 tk, & ! temperature in K from reference atmosphere file
                                                pmb, & ! pressure in milibars from reference file
                                                 nt    ! number density from the reference file

  REAL, DIMENSION (nl) ::  x,xx          ! dummy variable

  REAL, DIMENSION (np1) :: pnb1, &       ! Pressure in table #1
                            ef1          ! Esc.funct. #1, tabulated  

  REAL, DIMENSION (np2) :: pnb2, &       ! Pressure in table #2
                            ef2          ! Esc.funct. #2, tabulated

CONTAINS

!==================================================================
  SUBROUTINE VALVERDE_CO2_IR_LOOKUP (hr_co2_nlte,p,t,kms,kme,kts,kte)
!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------
    INTEGER, INTENT(IN) :: kts, kte, kms, kme
    REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: p, t
    REAL, DIMENSION( kms:kme ), INTENT(  OUT) :: hr_co2_nlte
 
    ! LOCAL VARIABLES
 
    REAL, DIMENSION( kms:kme ) ::   co2, &  ! number density of CO2 (cm-3)
                                    o3p, &  ! "  of atomic oxygen
                                   n2co, &  ! "  of N2 + CO
                                     nd, &  ! number density calculated from P and T
                                   p_mb     ! pressure converted to mbars

    REAL, DIMENSION (1:nl) ::  o2vmr

    ! The mixing ratios of various constituents should be calculated internal
    ! to the program. For now using the reference profiles from the L.-V.2001
    ! paper - reference atmosphere #1.

    !
    ! Vectors and indexes for the tabulation of escape functions
    !
    REAL, DIMENSION (kms:kme):: escf1, &      ! Esc.funct.#1, interpolated
                                escf2, &      ! Esc.funct.#2, interpolated
                                  pxx         ! Pressure converted to log nbars
    !
    ! Local Constants
    !
    REAL ::     nu1, nu2                    ! frec. of energy levels
    REAL ::     imr1, imr2                  ! isotopic abundances
    REAL ::     hplanck, gamma, vlight      ! physical constants
    REAL ::     ee
    REAL ::     rfvt                        ! Collisional rate
    REAL ::     rfvto3p                     !     "
    REAL ::     rfvv                        !     "
    REAL ::     pa2bar, bar2pa              ! Pascals to bars conversion factors

    !
    ! Local variables for the main loop
    !
    REAL ::      n1, n2, co2t       ! ground populations
    REAL ::      l1, p1, p12        ! prod & losses
    REAL ::      l2, p2, p21
    REAL ::      tt                 ! dummy variable
    REAL ::      c1, c2             ! molecular constants
    REAL ::      ae1, ae2           ! einstein spontaneous emission 
    REAL ::      a1, a2, a12, a21    
    REAL ::      pl1, pl2
    REAL ::      el1, el2
    REAL ::      hr1, hr2           ! heating rate due to each band

    !
    ! Indexes
    !
    INTEGER :: 	i

    !
    ! rate coefficients
    !
    REAL ::        k19xca, k19xcb
    REAL ::        k19cap1, k19cap2
    REAL ::        k19cbp1, k19cbp2 
    REAL ::        d19c, d19cp1, d19cp2 
    REAL ::        k20xc, k20cp1, k20cp2
    REAL ::        k21xc, k21cp2

    !
    ! Data
    !
    REAL :: kBolts

    REAL, PARAMETER :: cpco2 = 709.          ! pure CO2 heat capacity at 179 K in J/kg/K
    REAL, PARAMETER :: Rco2 = 189.           ! pure CO2 gas constant

    !*************************************************************************
    !       PROGRAM  STARTS
    !*************************************************************************
    nu1 = 667.38
    nu2 = 662.3734
    hplanck = 6.6261E-27
    gamma = 1.191e-5
    vlight = 3.e10
    ee = 1.438769
    kBolts = 1.380662e-23

    pa2bar=1./101325. ! pascals to bars conversion factor
    bar2pa=101325.    ! bars to pascals conversion factor

    imr1 = 0.987
    imr2 = 0.00408 + 0.0112 
    rfvt = 0.1
    rfvto3p = 1.0
    rfvv = 0.1

    !
    ! Read volume mixing ratios for different species from Initial Reference
    ! Atmosphere and interpolate to the desired grid.
    !
    ! They are required the following parameters:
    ! NAME.....................DESCRIPTION........................DEFINITION
    !  nl            ..........   number of layers                   integer
    !  tk(nl)        ..........   kinetic temperature [K]            real
    !  co2(kms:kme)  ..........   number density of CO2  [cm-3]      real
    !  o3p(kms:kme)  ..........   number density of O    [cm-3]      real
    !  n2co(kms:kme) ..........   numer density of N2 + CO [cm-3]    real
    !  nt(nl)        ..........   total number density [cm-3]        real
    !  pnb(nl)       ..........   pressure [nb]                      real

    DO i=kts, kte
       nd(i)=p(i)/t(i)/kBolts*1.e-6        ! calculate number desnity [cm-3]
                                           ! from actual P and T
       p_mb(i)=p(i)*pa2bar*1.e3            ! convert pressure to mbars
    END DO

    ! interpolating vmrs from pmb to p levels here and calculating the number
    ! desnities
      
    CALL INTERPOL1(nl, pmb, co2vmr, n2vmr, covmr, o3pvmr, kms, kme, kts, kte, &
                   p_mb, nd, co2, n2co, o3p, o2vmr)
    !
    ! Read escape functions for the strong and weak bands from tabulated values.
    ! The data are stored in the two following files:
    !           ddts1_1b.dat          ( strong CO2 emission band)
    !           ddts2_1b.dat          ( weak CO2 emission band )
    ! Each table contains two values, pressure and escape function, at
    ! a number of arbitrary levels (currently 110 and 114). These can be
    ! defined at head of code as follows:
    !      parameter (np1 = 110)            !  # levels for esc.funct # 1
    !      parameter (np2 = 114)            !  # levels for esc.funct # 2
    !
    ! NAME.....................DESCRIPTION........................DEFINITION
    !  np1       ...........   # of levels for the strong band    integer
    !  np2       ...........      "                weak   band    integer
    !  pnb1(np1) ...........   alog( press [nb] )                 real
    !  pnb2(np2) ...........   alog( press [nb] )                 real
    !  ef1(np1)  ...........   esc.func for strong band           real
    !  ef2(np2)  ...........   esc.func for weak band             real
    !

    !
    ! Interpolate escape functions to the desired grid ( from pnb(nl)to
    ! p(kms:kmn) )
    ! They will be required the escape functions at the nl levels
    !
    !  NAME.....................DESCRIPTION........................DEFINITION
    !   escf1(nl) ...........   esc.func for strong band           real
    !   escf2(nl) ...........   esc.func for weak band             real
    !
    DO i=kts,kte 
       pxx(i) = LOG( p_mb(i) * 1.e6 )     ! ln p in [nb]
    END DO

    CALL INTERPOL2( escf2,pxx,kts, kte, kts, kte, ef2,pnb2,np2)   
    CALL INTERPOL2( escf1,pxx,kts, kte, kts, kte, ef1,pnb1,np1)
    !------------------------------------------------------------------------
  
    !
    ! MAIN LOOP , for each altitude:
    !

    main_loop: DO i=kte,kts, -1        
       ! molecular populations  
       n1 = co2(i) * imr1 
       n2 = co2(i) * imr2  
       co2t = n1 + n2 

       ! intermediate collisional rates
	
       tt = t(i)*t(i)

       k19xca = 4.2e-12 * exp( -2988.0/t(i) + 303930.0/tt )
       k19xcb = 2.1e-12 * exp( -2659.0/t(i) + 223052.0/tt )

       IF (t(i) <= 175.0) THEN 
          k19xca = 3.3e-15
          k19xcb = 7.6e-16
       END IF

       k19xca = k19xca * rfvt
       k19xcb = k19xcb * rfvt
       k19cap1 = k19xca * 2.0 * exp( -ee*nu1/t(i) ) 
       k19cap2 = k19xca * 2.0 * exp( -ee*nu2/t(i) ) 
       k19cbp1 = k19xcb * 2.0 * exp( -ee*nu1/t(i) ) 
       k19cbp2 = k19xcb * 2.0 * exp( -ee*nu2/t(i) ) 

       d19c = k19xca*co2t + k19xcb*n2co(i)
       d19cp1 = k19cap1*co2t + k19cbp1*n2co(i) 
       d19cp2 = k19cap2*co2t + k19cbp2*n2co(i) 

       k20xc = 3.e-12 * rfvto3p
       k20cp1 = k20xc * 2.0 * exp( -ee/t(i) * nu1 ) 
       k20cp2 = k20xc * 2.0 * exp( -ee/t(i) * nu2 ) 

       k21xc = 2.49e-11 * 0.5 * rfvv
       k21cp2 = k21xc * exp( - ee/t(i) * (nu2-nu1) )

       l1 = d19c + k20xc*o3p(i) + k21cp2*n2
       p1 = ( d19cp1 + k20cp1*o3p(i) ) * n1
       p12 = k21xc*n1

       l2 = d19c + k20xc*o3p(i) + k21xc*n1 
       p2 = ( d19cp2 + k20cp2*o3p(i) ) * n2 
       p21 = k21cp2*n2

       ! radiative rates
       ae1 = 1.3546d00 * 1.66 / 4.0 * escf1(i)
       ae2 = ( 1.3452d00 + 1.1878d00 ) * 1.66 / 4.0 * escf2(i)
       l1 = l1 + ae1
       l2 = l2 + ae2

       ! solving the system
       c1 = gamma*nu1**3. * 0.5
       c2 = gamma*nu2**3. * 0.5
       a1 = c1 * p1 / (n1*l1)
       a2 = c2 * p2 / (n2*l2)
       a12 = (nu1/nu2)**3. * n2/n1 * p12/l1
       a21 = (nu2/nu1)**3. * n1/n2 * p21/l2
       el2 = (a2 + a21 * a1 ) / ( 1.0 - a21 * a12 )
       el1 = a1 + a12 * el2
       pl1 = el1 * n1 / c1 
       pl2 = el2 * n2 / c2 
 
       !  heating rate
       hr1 = - hplanck*vlight * nu1 * ae1 * pl1
       hr2 = - hplanck*vlight * nu2 * ae2 * pl2
       hr_co2_nlte(i) = hr1 + hr2
       ! Convert to K/s from ergs/cm3/s
        
       hr_co2_nlte(i)=hr_co2_nlte(i)*0.1/(p(i)/Rco2/t(i))/cpco2

    END DO main_loop

  END SUBROUTINE VALVERDE_CO2_IR_LOOKUP

!==================================================================
  SUBROUTINE VALVERDE_CO2_IR_FIT (hr_co2_nlte,p,t,kms,kme,kts,kte)
!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------
    INTEGER, INTENT(IN) :: kts, kte, kms, kme
    REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: p, t
    REAL, DIMENSION( kms:kme ), INTENT(  OUT) :: hr_co2_nlte
 
    ! LOCAL VARIABLES
    REAL, DIMENSION( kms:kme ) ::  co2, & ! number density of CO2 (cm-3)
                                   o3p, & ! "  of atomic oxygen
                                  n2co, & ! "  of N2 + CO
                                    nd, & ! number density calculated from P and T
                                  p_mb    ! pressure converted to mbars

    ! The mixing ratios of various constituents should be calculated internal
    ! to the program. For now using the reference profiles from the L.-V.2001
    ! paper - reference atmosphere #1.

    REAL, DIMENSION (kms:kme) :: co2vmrk,o3pvmrk,n2vmrk,covmrk,o2vmrk, & ! vol.mixing ratios from the reference file
                                                                lnpmb    ! log pressure in milibars 

    REAL ::  x     ! dummy variable

    !
    ! Vectors and indexes for the tabulation of escape functions
    !
    REAL, DIMENSION (kms:kme):: escf1, & ! Esc. funct. #1, interpolated
                                escf2, & ! Esc. funct. #2, interpolated
                                pxx      ! Pressure vector converted to log nbars
    !
    ! Local Constants
    !
    REAL ::     nu1, nu2                    ! frec. of energy levels
    REAL ::     imr1, imr2                  ! isotopic abundances
    REAL ::     hplanck, gamma, vlight      ! physical constants
    REAL ::     ee
    REAL ::     rfvt                        ! Collisional rate
    REAL ::     rfvto3p                     !     "
    REAL ::     rfvv                        !     "
    REAL ::     pa2bar, bar2pa              ! Pascals to bars conversion factors

    !
    ! Local variables for the main loop
    !
    REAL ::      n1, n2, co2t       ! ground populations
    REAL ::      l1, p1, p12        ! prod & losses
    REAL ::      l2, p2, p21
    REAL ::      tt                 ! dummy variable
    REAL ::      c1, c2             ! molecular constants
    REAL ::      ae1, ae2           ! einstein spontaneous emission 
    REAL ::      a1, a2, a12, a21    
    REAL ::      pl1, pl2
    REAL ::      el1, el2
    REAL ::      hr1, hr2           ! heating rate due to each band

    !
    ! Indexes
    !
    INTEGER :: 	i

    !
    ! rate coefficients
    !
    REAL ::        k19xca, k19xcb
    REAL ::        k19cap1, k19cap2
    REAL ::        k19cbp1, k19cbp2 
    REAL ::        d19c, d19cp1, d19cp2 
    REAL ::        k20xc, k20cp1, k20cp2
    REAL ::        k21xc, k21cp2

    !
    ! Data
    !
    REAL :: kBolts

    REAL, PARAMETER :: cpco2 = 709.          ! pure CO2 heat capacity at 179 K in J/kg/K
    REAL, PARAMETER :: Rco2 = 189.           ! pure CO2 gas constant

    REAL, PARAMETER :: pxx1_1=-2.00685, &  ! limit parameters for using in the fit formulas
                       pxx1_2=-0.128602, & ! to calculate escape functions
                       pxx1_3=5.90866, &
                       pxx2_1=-2.92685, &
                       pxx2_2=6.93423, &
                       a1_0=9.88608, &
                       a1_1=9.12538, &
                       a1_2=0.900324,&
                       a1_3=0.09369, &
                       a1_4=-0.9027,&
                       a2_0=0.0635917,&
                       a2_1=2.66337,&
                       a2_2=0.996555       ! parameters of the fit functions


    !*************************************************************************
    !       PROGRAM  STARTS
    !*************************************************************************
    nu1 = 667.38
    nu2 = 662.3734
    hplanck = 6.6261E-27
    gamma = 1.191e-5
    vlight = 3.e10
    ee = 1.438769
    kBolts = 1.380662e-23

    pa2bar=1./101325. ! pascals to bars conversion factor
    bar2pa=101325.    ! bars to pascals conversion factor

    imr1 = 0.987
    imr2 = 0.00408 + 0.0112 
    rfvt = 0.1
    rfvto3p = 1.0
    rfvv = 0.1

    ! Instead of reading the vmr data from the file and interpolating to the
    ! needed pressure levels we could use an analytical fit function to speed
    ! up calculations.

    DO i=kts, kte                 
       nd(i)    = p(i)/t(i)/kBolts*1.e-6   ! number desnity [cm-3]
       p_mb(i)  = p(i)*pa2bar*1.e3         ! convert GCM pressure to mbars
       lnpmb(i) = LOG(p_mb(i))             ! form a natural log of p_mb for interpolation formulas

       ! CO2 vmr calculations from fit to tabulated data
       IF (p_mb(i)>=6.5958) THEN 
          co2vmrk(i) = 0.9561
       ELSE IF (p_mb(i)<=3.9511e-08) THEN
          co2vmrk(i) = 8.327e-01   
       ELSE 
          co2vmrk(i) = -4.66E-06*lnpmb(i)**4  &
                       -7.01E-05*lnpmb(i)**3  &
                       -2.52E-04*lnpmb(i)**2  &
                       +5.12E-05*lnpmb(i)     &
                       +9.57E-01
       END IF

       ! N2 vmr calculations from fit to tabulated data
       IF (p_mb(i)>=6.5958) THEN
          n2vmrk(i) = 2.705e-02
       ELSE IF (p_mb(i)<=3.9511e-08) THEN
          n2vmrk(i) = 1.214e-01   
       ELSE 
          n2vmrk(i) = 2.86E-06*lnpmb(i)**4  &
                     +3.02E-05*lnpmb(i)**3  &
                     -1.02E-05*lnpmb(i)**2  &
                     -3.46E-04*lnpmb(i)     &
                     +2.72E-02
       END IF

       ! CO vmr calculations from fit to tabulated data
       IF (p_mb(i)>=6.5958) THEN
          covmrk(i) = 9.617e-04 
       ELSE IF (p_mb(i)<=3.9511e-08) THEN
          covmrk(i) = 7.278e-03    
       ELSE 
          covmrk(i) = 1.60E-07*lnpmb(i)**4  &
                     +2.48E-06*lnpmb(i)**3  &
                     +2.15E-05*lnpmb(i)**2  &
                     +6.27E-05*lnpmb(i)     &
                     +7.68E-04
       END IF

       ! O3P vmr calculations from fit to tabulated data
       IF (p_mb(i)>=6.5958) THEN
          o3pvmrk(i) = 1.583e-10
       ELSE IF (p_mb(i)<=3.9511e-08) THEN
          o3pvmrk(i) = 1.988e-02 
       ELSE 
          o3pvmrk(i) = EXP(-0.0059*lnpmb(i)**3 &
                           -0.192 *lnpmb(i)**2 &
                           -2.3282*lnpmb(i)    &
                           -16.977)
       END IF

       co2(i)  = nd(i)*co2vmrk(i)
       n2co(i) = nd(i)*(n2vmrk(i) + covmrk(i))
       o3p(i)  = nd(i)*o3pvmrk(i)         

    END DO

    ! Instead of opening the files with tabulated escape functions data and
    ! interpolating to the needed pressure leveles we'll use an analytical
    ! function that fits the tabulated data (AAP 9/12/05)

    DO i=kts,kte 
       pxx(i) = LOG( p(i)*pa2bar * 1.e9 )     ! convert pressure to  ln p [nb]

       IF (pxx(i)<=pxx1_1) THEN
          escf1(i)=1.            
       ELSE IF ((pxx(i)<=pxx1_2).AND.(pxx(i)>pxx1_1)) THEN
          escf1(i)=1./(a1_0*a1_1**(pxx(i))+a1_2)
       ELSE IF ((pxx(i)<=pxx1_3).AND.(pxx(i)>pxx1_2)) THEN
          escf1(i)=a1_3*exp(a1_4*pxx(i))
       ELSE IF (pxx(i)>pxx1_3) THEN
          escf1(i)=0.000458
       END IF

       IF (pxx(i)<=pxx2_1) THEN
          escf2(i)=1.            
       ELSE IF ((pxx(i)<=pxx2_2).AND.(pxx(i)>pxx2_1)) THEN
          escf2(i)=1./(a2_0*a2_1**(pxx(i))+a2_2)
       ELSE IF (pxx(i)>pxx2_2) THEN
          escf2(i)=0.00199
       END IF
    END DO

    !----------------------------------------------------------------------
  
    !
    ! MAIN LOOP , for each altitude:
    !
    main_loop: DO i=kte,kts, -1        

       ! molecular populations  
       n1 = co2(i) * imr1 
       n2 = co2(i) * imr2  
       co2t = n1 + n2 

       ! intermediate collisional rates
       tt = t(i)*t(i)

       k19xca = 4.2e-12 * exp( -2988.0/t(i) + 303930.0/tt )
       k19xcb = 2.1e-12 * exp( -2659.0/t(i) + 223052.0/tt )

       IF (t(i) <= 175.0) THEN 
          k19xca = 3.3e-15
          k19xcb = 7.6e-16
       END IF

       k19xca = k19xca * rfvt
       k19xcb = k19xcb * rfvt
       k19cap1 = k19xca * 2.0 * exp( -ee*nu1/t(i) ) 
       k19cap2 = k19xca * 2.0 * exp( -ee*nu2/t(i) ) 
       k19cbp1 = k19xcb * 2.0 * exp( -ee*nu1/t(i) ) 
       k19cbp2 = k19xcb * 2.0 * exp( -ee*nu2/t(i) ) 

       d19c = k19xca*co2t + k19xcb*n2co(i)
       d19cp1 = k19cap1*co2t + k19cbp1*n2co(i) 
       d19cp2 = k19cap2*co2t + k19cbp2*n2co(i) 

       k20xc = 3.e-12 * rfvto3p
       k20cp1 = k20xc * 2.0 * exp( -ee/t(i) * nu1 ) 
       k20cp2 = k20xc * 2.0 * exp( -ee/t(i) * nu2 ) 

       k21xc = 2.49e-11 * 0.5 * rfvv
       k21cp2 = k21xc * exp( - ee/t(i) * (nu2-nu1) )

       l1 = d19c + k20xc*o3p(i) + k21cp2*n2
       p1 = ( d19cp1 + k20cp1*o3p(i) ) * n1
       p12 = k21xc*n1

       l2 = d19c + k20xc*o3p(i) + k21xc*n1 
       p2 = ( d19cp2 + k20cp2*o3p(i) ) * n2 
       p21 = k21cp2*n2

       ! radiative rates
       ae1 = 1.3546d00 * 1.66 / 4.0 * escf1(i)
       ae2 = ( 1.3452d00 + 1.1878d00 ) * 1.66 / 4.0 * escf2(i)
       l1 = l1 + ae1
       l2 = l2 + ae2 

       ! solving the system
       c1 = gamma*nu1**3. * 0.5
       c2 = gamma*nu2**3. * 0.5
       a1 = c1 * p1 / (n1*l1)
       a2 = c2 * p2 / (n2*l2)
       a12 = (nu1/nu2)**3. * n2/n1 * p12/l1
       a21 = (nu2/nu1)**3. * n1/n2 * p21/l2
       el2 = (a2 + a21 * a1 ) / ( 1.0 - a21 * a12 )
       el1 = a1 + a12 * el2
       pl1 = el1 * n1 / c1 
       pl2 = el2 * n2 / c2 
 
       !  heating rate
       hr1 = - hplanck*vlight * nu1 * ae1 * pl1
       hr2 = - hplanck*vlight * nu2 * ae2 * pl2
       hr_co2_nlte(i) = hr1 + hr2

       ! Convert to K/s from ergs/cm3/s
        
       hr_co2_nlte(i)=hr_co2_nlte(i)*0.1/(p(i)/Rco2/t(i))/cpco2

    END DO main_loop

  END SUBROUTINE VALVERDE_CO2_IR_FIT

!==================================================================
  SUBROUTINE INTERPOL1 (nl, pmb, co2vmr, n2vmr, covmr, o3pvmr, kms, kme, &
                        kts, kte, p, nd, co2, n2co, o3p, o2vmr)
!------------------------------------------------------------------
! interpolate vmrs from pmb levels in the refernce atmosphere table to 
! the p levels in the gcm and calculate number densities of different
! species at this levels.
! Ideally, we would want to have the number densities supplied by the
! GCM, but for now we have to use a reference atmospheric profile
! from L.-V. 2001 paper.
!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nl, kms,kme, kts, kte              ! dimensions of the 2 domains
    REAL, DIMENSION(1:nl), INTENT(IN) :: pmb, co2vmr, n2vmr, covmr, o3pvmr ! orginial data in the table
    REAL, DIMENSION(kms:kme), INTENT(IN) :: p, nd ! pressure from GCM converted to mbars and the number density in cm-3
    REAL, DIMENSION(kms:kme), INTENT(OUT) :: co2, n2co, o3p, o2vmr !calculated interpolated number densities

    ! LOCAL VARIABLES
    INTEGER :: j, &    ! index such that x1(j)<x2(j)<x1(j+1)
               i       ! cycle variable   
      
    DO i=kts, kte
       j=LOCATE(pmb,p(i))
       IF (j>=nl) THEN
          co2(i)=nd(i)*co2vmr(nl)
          n2co(i) = nd(i)*( n2vmr(nl) + covmr(nl) )
          o3p(i) = nd(i) * o3pvmr(nl)
          o2vmr(i) = covmr(nl)                 ! not sure what this on is for - it is not used anywhere
       ELSE IF (j<=1) THEN
          co2(i)=nd(i)*co2vmr(1)
          n2co(i) = nd(i)*( n2vmr(1) + covmr(1) )
          o3p(i) = nd(i) * o3pvmr(1)
          o2vmr(i) = covmr(1)
       ELSE 
          co2(i)=nd(i)*(co2vmr(j)+(p(i)-pmb(j))*(co2vmr(j+1)-co2vmr(j))/(pmb(j+1)-pmb(j)))
          n2co(i)=nd(i)*(n2vmr(j)+(p(i)-pmb(j))*(n2vmr(j+1)-n2vmr(j))/(pmb(j+1)-pmb(j)))
          n2co(i)=n2co(i)+nd(i)*(covmr(j)+(p(i)-pmb(j))*(covmr(j+1)-covmr(j))/(pmb(j+1)-pmb(j)))
          o3p(i)=nd(i)*(o3pvmr(j)+(p(i)-pmb(j))*(o3pvmr(j+1)-o3pvmr(j))/(pmb(j+1)-pmb(j)))
          o2vmr(i)=covmr(i)             
       END IF
    END DO

  END SUBROUTINE INTERPOL1

!==================================================================
  SUBROUTINE INTERPOL2 (y2,x2,n2s,n2e,n2ts,n2te,y1,x1,n1)
!------------------------------------------------------------------
! Interpolate escape functions from the provided tables to the 
! p levels in the GCM
!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n1, n2s, n2e, n2ts, n2te      ! dimensions of the two domains
    REAL, DIMENSION(1:n1  ), INTENT(IN   ) :: y1, x1     ! y and x of the orginial data
    REAL, DIMENSION( n2s:n2e ), INTENT(  OUT) :: y2, x2  ! y and x of the data we are interpolating to

    ! LOCAL VARIABLES
    INTEGER :: J, &    ! index such that x1(j)<x2(j)<x1(j+1)
               i       ! cycle variable   
      
    DO i=n2ts, n2te
       j=LOCATE(x1,x2(i))
       IF (j>=n1) THEN
          y2(i)=y1(n1)
       ELSE IF (j<=1) THEN
          y2(i)=y1(1)
       ELSE 
          y2(i)=y1(j)+(x2(i)-x1(j))*(y1(j+1)-y1(j))/(x1(j+1)-x1(j))
       END IF
    END DO

  END SUBROUTINE INTERPOL2

!-------------------------------------------------------------------

 INTEGER FUNCTION FIND_J(x,xx,n)
!-------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL, DIMENSION (1:n) :: xx
      REAL, INTENT (IN) :: x

! LOCAL VARIABLES

      INTEGER :: j, jl, jm, ju

!-------------------------------------------------------------------

      jl=0
      ju=n+1
      
      DO WHILE((ju-jl)>1)
         jm=(ju+jl)/2
         IF ((xx(n)>=xx(1)).EQV.(x>=xx(jm))) THEN
            jl=jm
         ELSE
            ju=jm
         END IF
      END DO
      IF (x==xx(1)) THEN
         j=1
      ELSE IF (x==xx(n)) THEN
         j=n-1
      ELSE 
         j=jl
      END IF
      FIND_J=j
 END FUNCTION FIND_J
!-------------------------------------------------------------------

!==================================================================
  SUBROUTINE SWCORRECTION (p,xltecorrection,kms,kme,kts,kte)
!------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------

!C     XLTEFACTOR - correction factor for over-prediction of heating rates
!C                  by the LTE code.  Occurs at pressures where NLTE processes
!C                  become important
!C
!C     XLTEPRESSURE - pressure regime at which each factor is to be applied
!C
!C     XLTECORRECTION - interpolated correction applied to the sw heating rate

    INTEGER, PARAMETER :: l_nlte = 51

    REAL(KIND(0.d0)), DIMENSION( kms:kme ), INTENT(IN   )  :: p

    INTEGER :: malt,n,kms,kme,kts,kte,k

!C  LTE/NLTE heating rate conversion factors and associated pressures
      
    REAL(KIND(0.d0)) :: xltefactor(l_nlte),xltepressure(l_nlte)
    REAL(KIND(0.d0)), DIMENSION (kms:kme), INTENT(  OUT) :: xltecorrection

!------------------------------------------------------------------

      DATA  xltefactor / 1.007,  1.007,  1.007,  1.007,  1.008,&
                         1.009,  1.013,  1.017,  1.021,  1.027,&
                         1.033,  1.040,  1.047,  1.054,  1.061,&
                         1.069,  1.078,  1.087,  1.099,  1.112,&
                         1.129,  1.149,  1.174,  1.209,  1.263,&
                         1.344,  1.459,  1.608,  1.796,  2.033,&
                         2.343,  2.777,  3.404,  4.366,  5.856,&
                         8.161, 11.908, 17.791, 27.030, 40.897,&
                        60.798, 82.660,111.855,146.365,185.922,&
                       237.435,284.487,346.883,422.071,513.978,&
                       635.594 /

      DATA xltepressure / 1.122E+1, 8.822, 6.912, 5.397, 4.202,            &
                          3.261, 2.523, 1.946, 1.497, 1.149,               &
                          8.793E-1, 6.713E-1, 5.115E-1, 3.890E-1, 2.953E-1,&
                          2.239E-1, 1.696E-1, 1.284E-1, 9.719E-2, 7.355E-2,&
                          5.566E-2, 4.212E-2, 3.187E-2, 2.412E-2, 1.825E-2,&
                          1.381E-2, 1.045E-2, 7.908E-3, 5.984E-3, 4.530E-3,&
                          3.430E-3, 2.600E-3, 1.973E-3, 1.500E-3, 1.145E-3,&
                          8.761E-4, 6.731E-4, 5.192E-4, 4.024E-4, 3.141E-4,&
                          2.478E-4, 1.981E-4, 1.601E-4, 1.306E-4, 1.076E-4,&
                          8.939E-5, 7.462E-5, 6.242E-5, 5.234E-5, 4.397E-5,&
                          3.702E-5 /

!------------------------------------------------------------------

          DO k=kts, kte
             IF (p(k) > xltepressure(1)) THEN
                xltecorrection(k)= 1.0
             ELSE
                DO n=1,l_nlte-1
                   IF (p(k) > xltepressure(n)) THEN
                      EXIT
                   ELSE
                      malt = n
                   ENDIF
                ENDDO
                xltecorrection(k) = xltefactor(malt)+     &
                   (xltefactor(malt+1)-xltefactor(malt))* &
                   DLOG(p(k)/xltepressure(malt))/         &
                   DLOG(xltepressure(malt+1)/             &
                   xltepressure(malt))
              ENDIF
          END DO

 END SUBROUTINE SWCORRECTION
!-------------------------------------------------------------------

 SUBROUTINE VALVERDEINIT

   IMPLICIT NONE

   INTEGER :: i

   ! File unit numbers changed to prevent possible accidental equivalence
   ! with various WRF I/O standard file unit numbers.  Hopefully numbers
   ! in the 2000's should be safe...
   OPEN (UNIT=2020, FILE='./Data/radiation/marte_ref1.dat', STATUS='OLD')
   READ (2020,*)
   DO i=1,nl
      READ (2020,*) x(i), &
                  tk(i), &   ! temperature [K]
                  pmb(i), &  ! pressure [mb] 
                  nt(i),  &  ! number density [cm-3]
                  co2vmr(i), n2vmr(i), covmr(i), o3pvmr(i), xx(i) ! vol. mixing ratios
    END DO
    CLOSE (UNIT=2020)

    OPEN(UNIT=2022, FILE='./Data/radiation/ddts2_1b.dat', STATUS='OLD')
    READ (2022,*)
    DO i=1,np2
       READ (2022,*) pnb2(i),ef2(i)
    END DO
    CLOSE (UNIT=2022)

    OPEN(UNIT=2021, FILE='./Data/radiation/ddts1_1b.dat', STATUS='OLD')
    READ (2021,*)
    DO i=1,np1
       READ (2021,*) pnb1(i), ef1(i)
    END DO
    CLOSE (UNIT=2021)

  END SUBROUTINE VALVERDEINIT

END MODULE module_ra_valverde
