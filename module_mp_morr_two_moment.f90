!WRF:MODEL_LAYER:PHYSICS
!

! THIS MODULE CONTAINS THE TWO-MOMENT MICROPHYSICS CODE DESCRIBED BY
!     MORRISON ET AL. (2009, MWR)

! CHANGES FOR V3.2, RELATIVE TO MOST RECENT (BUG-FIX) CODE FOR V3.1

! 1) ADDED ACCELERATED MELTING OF GRAUPEL/SNOW DUE TO COLLISION WITH RAIN, FOLLOWING LIN ET AL. (1983)
! 2) INCREASED MINIMUM LAMBDA FOR RAIN, AND ADDED RAIN DROP BREAKUP FOLLOWING MODIFIED VERSION
!     OF VERLINDE AND COTTON (1993)
! 3) CHANGE MINIMUM ALLOWED MIXING RATIOS IN DRY CONDITIONS (RH < 90%), THIS IMPROVES RADAR REFLECTIIVITY
!     IN LOW REFLECTIVITY REGIONS
! 4) BUG FIX TO MAXIMUM ALLOWED PARTICLE FALLSPEEDS AS A FUNCTION OF AIR DENSITY
! 5) BUG FIX TO CALCULATION OF LIQUID WATER SATURATION VAPOR PRESSURE (CHANGE IS VERY MINOR)
! 6) INCLUDE WRF CONSTANTS PER SUGGESTION OF JIMY

! bug fix, 5/12/10
! 7) bug fix for saturation vapor pressure in low pressure, to avoid division by zero
! 8) include 'EP2' WRF constant for saturation mixing ratio calculation, instead of hardwire constant

! CHANGES FOR V3.3
! 1) MODIFICATION FOR COUPLING WITH WRF-CHEM (PREDICTED DROPLET NUMBER CONCENTRATION) AS AN OPTION
! 2) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
! 3) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
! 4) CLEAN UP OF COMMENTS IN THE CODE
    
! additional minor bug fixes and small changes, 5/30/2011
! minor revisions by A. Ackerman April 2011:
! 1) replaced kinematic with dynamic viscosity 
! 2) replaced scaling by air density for cloud droplet sedimentation
!    with viscosity-dependent Stokes expression
! 3) use Ikawa and Saito (1991) air-density scaling for cloud ice
! 4) corrected typo in 2nd digit of ventilation constant F2R

! additional fixes:
! 5) TEMPERATURE FOR ACCELERATED MELTING DUE TO COLLIIONS OF SNOW AND GRAUPEL
!    WITH RAIN SHOULD USE CELSIUS, NOT KELVIN (BUG REPORTED BY K. VAN WEVERBERG)
! 6) NPRACS IS NOT SUBTRACTED FROM SNOW NUMBER CONCENTRATION, SINCE
!    DECREASE IN SNOW NUMBER IS ALREADY ACCOUNTED FOR BY NSMLTS 
! 7) fix for switch for running w/o graupel/hail (cloud ice and snow only)

! hm bug fix 3/16/12

! 1) very minor change to limits on autoconversion source of rain number when cloud water is depleted

! WRFV3.5
! hm/A. Ackerman bug fix 11/08/12

! 1) for accelerated melting from collisions, should use rain mass collected by snow, not snow mass 
!    collected by rain
! 2) minor changes to some comments
! 3) reduction of maximum-allowed ice concentration from 10 cm-3 to 0.3
!    cm-3. This was done to address the problem of excessive and persistent
!    anvil cirrus produced by the scheme.

! CHANGES FOR WRFV3.5.1
! 1) added output for snow+cloud ice and graupel time step and accumulated
!    surface precipitation
! 2) bug fix to option w/o graupel/hail (IGRAUP = 1), include PRACI, PGSACW,
!    and PGRACS as sources for snow instead of graupel/hail, bug reported by
!    Hailong Wang (PNNL)
! 3) very minor fix to immersion freezing rate formulation (negligible impact)
! 4) clarifications to code comments
! 5) minor change to shedding of rain, remove limit so that the number of 
!    collected drops can smaller than number of shed drops
! 6) change of specific heat of liquid water from 4218 to 4187 J/kg/K

! CHANGES FOR WRFV3.6.1
! 1) minor bug fix to melting of snow and graupel, an extra factor of air density (RHO) was removed
!    from the calculation of PSMLT and PGMLT
! 2) redundant initialization of PSMLT (non answer-changing)

! CHANGES FOR WRFV3.8.1
! 1) changes and cleanup of code comments
! 2) correction to universal gas constant (very small change)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL/HAIL.

#define KVALUE 20
#define JVALUE 11
#define IVALUE 20
#ifdef MARS_MORR_DEBUG
#define CALL_DEBUG(a,b) call debugval(a,b)
#endif
MODULE MODULE_MP_MORR_TWO_MOMENT
   USE     module_wrf_error
!      USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm  ! GT
!      USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep  ! GT

! USE WRF PHYSICS CONSTANTS
  use module_model_constants, ONLY: CP, G, R => r_d, RV => r_v, EP_2
#ifdef WRF_MARS
  use module_model_constants, only: co2_molmass, co2_molrad, N_avogadro, kboltz => bk, &
                                    mass_h2o, mass_co2, h2o_molrad, PI_S,SQRTPI
  use module_planet_utilities, only : esat, tfrost_co2, prep_gaussian_weights
#endif

!  USE module_state_description

   IMPLICIT NONE

#if ( WRF_MARS == 1 )
   REAL, PARAMETER :: PI = PI_S
#else
   REAL, PARAMETER :: PI = 3.1415926535897932384626434
#endif
   REAL, PARAMETER :: xxx = 0.9189385332046727417803297

   PUBLIC  ::  MP_MORR_TWO_MOMENT
#ifdef WRF_MARS
   PUBLIC  :: constrain_mars_waterice, constrain_mars_dust, average_mass_of_particle
#endif
   PUBLIC  ::  POLYSVP

   PRIVATE :: GAMMA, DERF1
   PRIVATE :: PI, xxx
   PRIVATE :: MORR_TWO_MOMENT_MICRO

#ifdef WRF_MARS
! set this to true to get error messages if we're getting anything other than
! a vapour<->ice<-cloud system (it shouldn't rain on Mars - at least not today!)
   LOGICAL, PARAMETER :: debug = .true.

! should probably have a namelist flag that comes through that allows this code to
! only treat dust for if we're running the two-moment dust but with water off
! - can't do the inverse of this since there is no homogeneous nucleation (i.e. 
! homogeneous nuc is so hard in the Mars environment that it probably does not
! happen)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SWITCHES FOR MICROPHYSICS SCHEME
! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
! IACT = 3, ACTIVATION CALCULATED IN MODULE_MIXACTIVATE

     INTEGER, PRIVATE ::  IACT

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   
! !!!NOTE: PREDICTED DROPLET CONCENTRATION NOT AVAILABLE IN THIS VERSION
! CONTACT HUGH MORRISON (morrison@ucar.edu) FOR FURTHER INFORMATION

     INTEGER, PRIVATE ::  INUM

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (CM-3)
     REAL, PRIVATE ::      NDCNST

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

     INTEGER, PRIVATE ::  ILIQ

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS
#ifdef WRF_MARS
!      = 2, USE HETEROGENEOUS NUC ON DUST FOR MARS INC. CZICZO [2013] CONTACT PARAM VAR WITH T
!      = 3, pseudo homogeneous ice nucleation following simple water scheme
!      = 4, pseudo heterogeneous ice nucleation with growth
!      = 5, pseudo heterogeneous ice nucleation requiring dust presence, with growth
!      = 6, heterogeneous ice nucleation, omg.
    logical, private :: mars_heterogeneous, mars_growth, mars_dust_scavenging, &
         nucleation_monolayer, nucleation_surface_diffusion
    real, private :: nucleation_m
#endif
     INTEGER, PRIVATE ::  INUC

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
!             NON-EQULIBRIUM SUPERSATURATION, 
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION, BASED ON THE 
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0) IN NON-WRF-CHEM VERSION OF CODE

     INTEGER, PRIVATE ::  IBASE

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0) IN NON-WRF-CHEM VERSION OF CODE

     INTEGER, PRIVATE ::  ISUB      

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

     INTEGER, PRIVATE ::  IGRAUP

! HM ADDED NEW OPTION FOR HAIL
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING GICE IS HAIL

     INTEGER, PRIVATE ::  IHAIL

! CLOUD MICROPHYSICS CONSTANTS

     REAL, PRIVATE ::      AI,AC,AS,AR,AG ! 'A' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
     REAL, PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
!     REAL, PRIVATE ::      R           ! GAS CONSTANT FOR AIR
!     REAL, PRIVATE ::      RV          ! GAS CONSTANT FOR WATER VAPOR
!     REAL, PRIVATE ::      CP          ! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
     REAL, PRIVATE ::      RHOSU       ! STANDARD AIR DENSITY AT 850 MB
     REAL, PRIVATE ::      RHOW        ! DENSITY OF LIQUID WATER
     REAL, PRIVATE ::      RHOI        ! BULK DENSITY OF CLOUD ICE
     REAL, PRIVATE ::      RHOSN       ! BULK DENSITY OF SNOW
     REAL, PRIVATE ::      RHOG        ! BULK DENSITY OF GRAUPEL
     REAL, PRIVATE ::      AIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL, PRIVATE ::      BIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL, PRIVATE ::      ECR         ! COLLECTION EFFICIENCY BETWEEN DROPLETS/RAIN AND SNOW/RAIN
     REAL, PRIVATE ::      DCS         ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
     REAL, PRIVATE ::      MI0         ! INITIAL SIZE OF NUCLEATED CRYSTAL
     REAL, PRIVATE ::      MG0         ! MASS OF EMBRYO GRAUPEL
     REAL, PRIVATE ::      F1S         ! VENTILATION PARAMETER FOR SNOW
     REAL, PRIVATE ::      F2S         ! VENTILATION PARAMETER FOR SNOW
     REAL, PRIVATE ::      F1R         ! VENTILATION PARAMETER FOR RAIN
     REAL, PRIVATE ::      F2R         ! VENTILATION PARAMETER FOR RAIN
!     REAL, PRIVATE ::      G           ! GRAVITATIONAL ACCELERATION
     REAL, PRIVATE ::      QSMALL      ! SMALLEST ALLOWED HYDROMETEOR MIXING RATIO
     REAL, PRIVATE ::      CI,DI,CS,DS,CG,DG ! SIZE DISTRIBUTION PARAMETERS FOR CLOUD ICE, SNOW, GRAUPEL
     REAL, PRIVATE ::      EII         ! COLLECTION EFFICIENCY, ICE-ICE COLLISIONS
     REAL, PRIVATE ::      ECI         ! COLLECTION EFFICIENCY, ICE-DROPLET COLLISIONS
     REAL, PRIVATE ::      RIN     ! RADIUS OF CONTACT NUCLEI (M)
! hm, add for V3.2
     REAL, PRIVATE ::      CPW     ! SPECIFIC HEAT OF LIQUID WATER

! CCN SPECTRA FOR IACT = 1

     REAL, PRIVATE ::      C1     ! 'C' IN NCCN = CS^K (CM-3)
     REAL, PRIVATE ::      K1     ! 'K' IN NCCN = CS^K

! AEROSOL PARAMETERS FOR IACT = 2

     REAL, PRIVATE ::      MW      ! MOLECULAR WEIGHT WATER (KG/MOL)
     REAL, PRIVATE ::      OSM     ! OSMOTIC COEFFICIENT
     REAL, PRIVATE ::      VI      ! NUMBER OF ION DISSOCIATED IN SOLUTION
     REAL, PRIVATE ::      EPSM    ! AEROSOL SOLUBLE FRACTION
     REAL, PRIVATE ::      RHOA    ! AEROSOL BULK DENSITY (KG/M3)
     REAL, PRIVATE ::      MAP     ! MOLECULAR WEIGHT AEROSOL (KG/MOL)
     REAL, PRIVATE ::      MA      ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)
     REAL, PRIVATE ::      RR      ! UNIVERSAL GAS CONSTANT
     REAL, PRIVATE ::      BACT    ! ACTIVATION PARAMETER
     REAL, PRIVATE ::      RM1     ! GEOMETRIC MEAN RADIUS, MODE 1 (M)
     REAL, PRIVATE ::      RM2     ! GEOMETRIC MEAN RADIUS, MODE 2 (M)
     REAL, PRIVATE ::      NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
     REAL, PRIVATE ::      NANEW2  ! TOTAL AEROSOL CONCENTRATION, MODE 2 (M^-3)
     REAL, PRIVATE ::      SIG1    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 1
     REAL, PRIVATE ::      SIG2    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 2
     REAL, PRIVATE ::      F11     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL, PRIVATE ::      F12     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL, PRIVATE ::      F21     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     REAL, PRIVATE ::      F22     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2     
     REAL, PRIVATE ::      MMULT   ! MASS OF SPLINTERED ICE PARTICLE
     REAL, PRIVATE ::      LAMMAXI,LAMMINI,LAMMAXR,LAMMINR,LAMMAXS,LAMMINS,LAMMAXG,LAMMING
#ifdef WRF_MARS
	 REAL, PRIVATE ::   LAMMAXD, LAMMIND, LAMDEFD, LAMD, DUST_CONS0, DUST_CONS1, DUST_CONS2, DUST_CONS3, dust_cons3_a, DUST_MU, &
                        DUST_CONS_A_SEDIM, DUST_CONS_B_SEDIM, DUST_CONS_E_SEDIM, DUST_PIRHOD, DUST_CONS_RHO, &
                        DUST_CONS_INV_G1, DUST_CONS_INV_G4, DUST_CONS_G3,&
                        nfac_lammaxd1, nfac_lammind1, nfac_lammaxd4, nfac_lammind4, &
                        nfac_lamdefd1, nfac_lamdefd4
  real(kind(0.d0)) :: dust_mu_r8
  real, private :: ice_consmin_q_n, ice_consmax_q_n, lammind_q_n, lammaxd_q_n, lamdefd_q_n
!    real, private :: lamdefi,nfac_lamdefd1, nfac_lamdefd4
    
    integer, private :: DISTRIBUTION_MARS_DUST
    real, private :: sigma_squared
    real, private :: dust_reff_max, dust_reff_min, dust_reff_default
    real, private :: ln_dust_reff_max, ln_dust_reff_min, ln_dust_reff_default
    
    !ice constants
    real, private :: ice_mu, ice_cons12, ice_cons11, ice_cons10, ice_cons2
    real, private :: ice_cons11min, ice_cons11max, ice_cons10min, ice_cons10max
    real, private :: ice_cons_nuc_m0,ice_cons_nuc_m2,ice_cons_nuc_m3
    real, private :: ice_cons_nuc_m0p5,ice_cons_nuc_m2p5,ice_cons_nuc_m3p5
    real, private :: homogeneous_ice_number_cons1
    !nucleation
    real, private :: ice_cons_prefactor
    real, private :: water_molec_volume, water_molec_radius
    real, dimension(:,:), allocatable :: hc_data
    integer :: hc_ncolumns, hc_nrows
    real :: hc_tmin, hc_deltaT
    real :: photolysis_scrit
#endif
! CONSTANTS TO IMPROVE EFFICIENCY

     REAL, PRIVATE :: CONS1,CONS2,CONS3,CONS4,CONS5,CONS6,CONS7,CONS8,CONS9,CONS10
     REAL, PRIVATE :: CONS11,CONS12,CONS13,CONS14,CONS15,CONS16,CONS17,CONS18,CONS19,CONS20
     REAL, PRIVATE :: CONS21,CONS22,CONS23,CONS24,CONS25,CONS26,CONS27,CONS28,CONS29,CONS30
     REAL, PRIVATE :: CONS31,CONS32,CONS33,CONS34,CONS35,CONS36,CONS37,CONS38,CONS39,CONS40
     REAL, PRIVATE :: CONS41

     logical :: iprint, jprint


!bin scheme
  integer, parameter :: nbins=6
  real(kind(0.d0)), dimension(nbins) :: number, mass
  real(kind(0.d0)), dimension(nbins) :: diam, vol, bin_width
  real(kind(0.d0)), dimension(2,nbins) :: weights_pos
  real(kind(0.d0)), dimension(nbins+1) :: diam_edges


#ifdef MARS_MORR_DEBUG
interface debugval  
  module procedure debugvaldouble, debugvalsingle, debugvalint,debugvallogical
end interface
#endif
CONTAINS
#ifdef MARS_MORR_DEBUG
subroutine debugvallogical(name, value)
      implicit none
        character(len=*), intent(in) :: name
        logical, intent(in) :: value

        if (iprint) write(*,*) name, " = ", value
end subroutine debugvallogical
    subroutine debugvalint(name, value)
      implicit none
        character(len=*), intent(in) :: name
        integer, intent(in) :: value

        if (iprint) write(*,*) name, " = ", value
    end subroutine debugvalint
    subroutine debugvaldouble(name, value)
      implicit none
        character(len=*), intent(in) :: name
        real(kind(0.d0)), intent(in) :: value

        if (iprint) write(*,*) name, " = ", value
    end subroutine debugvaldouble
    subroutine debugvalsingle(name, value)
      implicit none
        character(len=*), intent(in) :: name
        real, intent(in) :: value

        if (iprint) write(*,*) name, " = ", value
    end subroutine debugvalsingle
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MORR_TWO_MOMENT_INIT(                    &
#ifdef WRF_MARS
                                           rho_dust       &
                                          ,rhod_h2oi      &
                                          ,reff           &
                                          ,veff           &
                                          ,a_dust_sedim   &
                                          ,b_dust_sedim   &
                                          ,e_dust_sedim   &
                                          ,dust_distribution &
                                          ,homogeneous_ice_reff &
                                          ,mars_ice_nucleation &
                                          ,mars_nucleation_surface_diffusion &
                                          ,mars_nucleation_monolayer &
                                          ,mars_nucleation_m &
                                          ,nucleation_filename &
#endif
)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS 
! NEEDED BY THE MICROPHYSICS SCHEME.
! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE
#ifdef WRF_MARS
    real, intent(in) :: rho_dust,        &
                        rhod_h2oi,       &
                        reff,            &
                        veff,            &
                        a_dust_sedim,    &
                        b_dust_sedim,    &
                        e_dust_sedim
    real, intent(in) :: homogeneous_ice_reff
    integer, intent(in) :: dust_distribution, mars_ice_nucleation
    logical, intent(in) :: mars_nucleation_monolayer, mars_nucleation_surface_diffusion
    real, intent(in) :: mars_nucleation_m
    character(len=*), intent(in) :: nucleation_filename
    real :: volratio, diamratio
    real :: bin_min_diam, bin_max_diam, bin_edge_min_diam, bin_edge_max_diam
    real(kind(0.d0)), dimension(nbins+1) :: gloc
#endif
      integer n,i
      real(kind(0.d0)) :: pans, qans, fmin, fmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THE FOLLOWING PARAMETERS ARE USER-DEFINED SWITCHES AND NEED TO BE
! SET PRIOR TO CODE COMPILATION

! INUM IS AUTOMATICALLY SET TO 0 FOR WRF-CHEM BELOW,
! ALLOWING PREDICTION OF DROPLET CONCENTRATION
! THUS, THIS PARAMETER SHOULD NOT BE CHANGED HERE
! AND SHOULD BE LEFT TO 1

      INUM = 1

! SET CONSTANT DROPLET CONCENTRATION (UNITS OF CM-3)
! IF NO COUPLING WITH WRF-CHEM

      NDCNST = 250.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE, THE FOLLOWING OPTIONS RELATED TO DROPLET ACTIVATION 
! (IACT, IBASE, ISUB) ARE NOT AVAILABLE IN CURRENT VERSION
! FOR WRF-CHEM, DROPLET ACTIVATION IS PERFORMED 
! IN 'MIX_ACTIVATE', NOT IN MICROPHYSICS SCHEME


! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA

      IACT = 2

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
!             NON-EQULIBRIUM SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, 
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, BASED ON THE 
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)
 
      IBASE = 2

! INCLUDE SUB-GRID VERTICAL VELOCITY (standard deviation of w) IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! currently, sub-grid w is constant of 0.5 m/s (not coupled with PBL/turbulence scheme)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      ISUB = 0      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

      ILIQ = 0

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)
#ifdef WRF_MARS
!      = 2, USE MARS HETERO NUC ON DUST AND CONTACT PARAM TEMPERATURE DEPENDENCE
      
      nucleation_monolayer = .false.
      nucleation_surface_diffusion=.true.
      nucleation_m = 0.954
      
      nucleation_monolayer = mars_nucleation_monolayer
      nucleation_surface_diffusion = mars_nucleation_surface_diffusion

      nucleation_surface_diffusion = mars_nucleation_surface_diffusion
      nucleation_m = mars_nucleation_m

      INUC = mars_ice_nucleation
      mars_heterogeneous = .false.
      mars_growth=.false.
      mars_dust_scavenging = .true.
      if(inuc.ge.4) then 
        mars_growth=.true.
      endif
      if (inuc.eq.6 .or. inuc.eq.7) then
          mars_heterogeneous = .true.
          if(inuc .eq. 7) mars_dust_scavenging = .false.
          !call load_heterogeneous_coefficients(nucleation_filename)
      endif
#else
      INUC = 0
#endif

! SWITCH FOR GRAUPEL/HAIL NO GRAUPEL/HAIL
! IGRAUP = 0, INCLUDE GRAUPEL/HAIL
! IGRAUP = 1, NO GRAUPEL/HAIL

#ifdef WRF_MARS
    IGRAUP=1
#else
      IGRAUP = 0

#endif
! HM ADDED 11/7/07
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL
! NOTE ---> RECOMMEND IHAIL = 1 FOR CONTINENTAL DEEP CONVECTION

      IHAIL = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET PHYSICAL CONSTANTS

! FALLSPEED PARAMETERS (V=AD^B)
         AI = 700.
         AC = 3.E7
         AS = 11.72
         AR = 841.99667
         BI = 1.
         BC = 2.
         BS = 0.41
         BR = 0.8
         IF (IHAIL.EQ.0) THEN
	 AG = 19.3
	 BG = 0.37
         ELSE ! (MATSUN AND HUGGINS 1980)
         AG = 114.5 
         BG = 0.5
         END IF

! CONSTANTS AND PARAMETERS
!         R = 287.15
!         RV = 461.5
!         CP = 1005.
         RHOSU = 85000./(287.15*273.15)
         RHOW = 997.

#ifdef WRF_MARS
         RHOI = rhod_h2oi
#else
         RHOI = 500.
#endif

         RHOSN = 100.
         IF (IHAIL.EQ.0) THEN
	 RHOG = 400.
         ELSE
         RHOG = 900.
         END IF
         AIMM = 0.66
         BIMM = 100.
         ECR = 1.
#ifdef WRF_MARS
! Need to think about what cutoff size for ice vs. clouds. This must
! surely be something we estimate from the stokes-cunnigham to make it
! general - which means we need some cutoff fall speed above which we
! say it is snow
         DCS = 125.E-6 ! cutoff or transition diameter (ice-snow) in meters
         MI0 = 4./3.*PI*RHOI*(10.E-6)**3  ! initial nucleated ice particle mass
#else
         DCS = 125.E-6
         MI0 = 4./3.*PI*RHOI*(10.E-6)**3
#endif
	 MG0 = 1.6E-10
         F1S = 0.86
         F2S = 0.28
         F1R = 0.78
!         F2R = 0.32
! fix 053011
         F2R = 0.308
!         G = 9.806
         QSMALL = 1.E-14
         EII = 0.1
         ECI = 0.7
! HM, ADD FOR V3.2
! hm, 7/23/13
!         CPW = 4218.
         CPW = 4187.

! SIZE DISTRIBUTION PARAMETERS

         CI = RHOI*PI/6.
         DI = 3.
         CS = RHOSN*PI/6.
         DS = 3.
         CG = RHOG*PI/6.
         DG = 3.

! RADIUS OF CONTACT NUCLEI
         RIN = 0.1E-6

         MMULT = 4./3.*PI*RHOI*(5.E-6)**3

! SIZE LIMITS FOR LAMBDA
#ifdef WRF_MARS
         LAMMAXI = 1./5.E-7
         LAMMINI = 1./(2.*DCS+100.E-6)
         LAMMAXR = 1./1.E-6
!         LAMMINR = 1./500.E-6
         LAMMINR = 1./2800.E-6
#else
         LAMMAXI = 1./1.E-6
         LAMMINI = 1./(2.*DCS+100.E-6)
         LAMMAXR = 1./20.E-6
!         LAMMINR = 1./500.E-6
         LAMMINR = 1./2800.E-6

#endif
         LAMMAXS = 1./10.E-6
         LAMMINS = 1./2000.E-6
         LAMMAXG = 1./20.E-6
         LAMMING = 1./2000.E-6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! note: these parameters only used by the non-wrf-chem version of the 
!       scheme with predicted droplet number

! CCN SPECTRA FOR IACT = 1

! MARITIME
! MODIFIED FROM RASMUSSEN ET AL. 2002
! NCCN = C*S^K, NCCN IS IN CM-3, S IS SUPERSATURATION RATIO IN %

              K1 = 0.4
              C1 = 120. 

! CONTINENTAL

!              K1 = 0.5
!              C1 = 1000. 

! AEROSOL ACTIVATION PARAMETERS FOR IACT = 2
! PARAMETERS CURRENTLY SET FOR AMMONIUM SULFATE

         MW = 0.018
         OSM = 1.
         VI = 3.
         EPSM = 0.7
         RHOA = 1777.
         MAP = 0.132
         MA = 0.0284
! hm fix 6/23/16
!         RR = 8.3187
         RR = 8.3145
         BACT = VI*OSM*EPSM*MW*RHOA/(MAP*RHOW)

! AEROSOL SIZE DISTRIBUTION PARAMETERS CURRENTLY SET FOR MPACE 
! (see morrison et al. 2007, JGR)
! MODE 1

         RM1 = 0.052E-6
         SIG1 = 2.04
         NANEW1 = 72.2E6
         F11 = 0.5*EXP(2.5*(LOG(SIG1))**2)
         F21 = 1.+0.25*LOG(SIG1)

! MODE 2

         RM2 = 1.3E-6
         SIG2 = 2.5
         NANEW2 = 1.8E6
         F12 = 0.5*EXP(2.5*(LOG(SIG2))**2)
         F22 = 1.+0.25*LOG(SIG2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CONSTANTS FOR EFFICIENCY

         CONS1=GAMMA(1.+DS)*CS
         CONS2=GAMMA(1.+DG)*CG
         CONS3=GAMMA(4.+BS)/6.
         CONS4=GAMMA(4.+BR)/6.
         CONS5=GAMMA(1.+BS)
         CONS6=GAMMA(1.+BR)
         CONS7=GAMMA(4.+BG)/6.
         CONS8=GAMMA(1.+BG)
         CONS9=GAMMA(5./2.+BR/2.)
         CONS10=GAMMA(5./2.+BS/2.)
         CONS11=GAMMA(5./2.+BG/2.)
         CONS12=GAMMA(1.+DI)*CI
         CONS13=GAMMA(BS+3.)*PI/4.*ECI
         CONS14=GAMMA(BG+3.)*PI/4.*ECI
         CONS15=-1108.*EII*PI**((1.-BS)/3.)*RHOSN**((-2.-BS)/3.)/(4.*720.)
         CONS16=GAMMA(BI+3.)*PI/4.*ECI
         CONS17=4.*2.*3.*RHOSU*PI*ECI*ECI*GAMMA(2.*BS+2.)/(8.*(RHOG-RHOSN))
         CONS18=RHOSN*RHOSN
         CONS19=RHOW*RHOW
         CONS20=20.*PI*PI*RHOW*BIMM
         CONS21=4./(DCS*RHOI)
         CONS22=PI*RHOI*DCS**3/6.
         CONS23=PI/4.*EII*GAMMA(BS+3.)
         CONS24=PI/4.*ECR*GAMMA(BR+3.)
         CONS25=PI*PI/24.*RHOW*ECR*GAMMA(BR+6.)
         CONS26=PI/6.*RHOW
         CONS27=GAMMA(1.+BI)
         CONS28=GAMMA(4.+BI)/6.
         CONS29=4./3.*PI*RHOW*(25.E-6)**3
         CONS30=4./3.*PI*RHOW
         CONS31=PI*PI*ECR*RHOSN
         CONS32=PI/2.*ECR
         CONS33=PI*PI*ECR*RHOG
         CONS34=5./2.+BR/2.
         CONS35=5./2.+BS/2.
         CONS36=5./2.+BG/2.
         CONS37=4.*PI*1.38E-23/(6.*PI*RIN)
         CONS38=PI*PI/3.*RHOW
         CONS39=PI*PI/36.*RHOW*BIMM
         CONS40=PI/6.*BIMM
         CONS41=PI*PI*ECR*RHOW

#ifdef WRF_MARS
         DUST_PIRHOD = pi * rho_dust
         DUST_CONS_RHO = rho_dust
!compiler issues with dust_cons3
         dust_cons3_a = N_avogadro * co2_molrad * co2_molrad
         DUST_CONS3 = 1./ (4. * sqrt(2.) * pi) * co2_molmass / (dust_cons3_a)
         dust_reff_max = 25.e-6
         dust_reff_min = 0.1e-6
         dust_reff_default = reff
         
if(dust_distribution == 0) then !GAMMA function

         DUST_MU = (1-3*veff)/veff !mu_dust
         dust_mu_r8 = dust_mu
         DUST_CONS0 = GAMMA(DUST_MU + 4.) * DUST_PIRHOD / (6.)
         DUST_CONS1 = DUST_CONS0 / GAMMA(DUST_MU+1.)
         DUST_CONS2 = GAMMA(DUST_MU + 4.) / (2*GAMMA(DUST_MU + 3.))

         DUST_CONS_INV_G1 = 1./ GAMMA(DUST_MU + 1.)
         DUST_CONS_INV_G4 = 1./ GAMMA(DUST_MU + 4.)
         DUST_CONS_G3 = GAMMA(DUST_MU + 3.)
         
         DUST_CONS_A_SEDIM = a_dust_sedim
         DUST_CONS_B_SEDIM = b_dust_sedim
         DUST_CONS_E_SEDIM = e_dust_sedim

        LAMMIND = DUST_CONS2/dust_reff_max !25 micron
        LAMMAXD = DUST_CONS2/dust_reff_min !0.1 micron
        LAMDEFD = DUST_CONS2/dust_reff_default !2.0 micron

        !nfac_LAMMIND1 = GAMMA(DUST_MU+1.) * exp(-(DUST_MU+1.)*log(LAMMIND))
        !nfac_LAMMAXD1 = GAMMA(DUST_MU+1.) * exp(-(DUST_MU+1.)*log(LAMMAXD))

        !nfac_LAMMIND4 = exp((DUST_MU+4.)*log(LAMMIND)) / dust_cons0
        !nfac_LAMMAXD4 = exp((DUST_MU+4.)*log(LAMMAXD)) / dust_cons0

        !nfac_LAMDEFD1 = GAMMA(DUST_MU+1.) * exp(-(DUST_MU+1.)*log(LAMDEFD))
        !nfac_LAMDEFD4 = exp((DUST_MU+4.)*log(LAMDEFD)) / dust_cons0
        
        LAMMIND_Q_N = exp((3.)*log(LAMMIND)) / dust_cons1
        LAMMAXD_Q_N = exp((3.)*log(LAMMAXD)) / dust_cons1
        LAMDEFD_Q_N = exp((3.)*log(LAMDEFD)) / dust_cons1


! Bins Gauss Legendre
    gloc=0.
    call prep_gaussian_weights(nbins,.false.,weights_pos) !these are the fractions. yay.
    bin_width=weights_pos(1,:)
    gloc(:nbins) = weights_pos(2,:)
!    write(*,*) "loc", gloc
!    write(*,*) "bw", bin_width
    call gamma_inc_inv_calc(dust_mu_r8,gloc(:nbins),diam,nbins,"diam")
    
    gloc(1)=0.0
    do i=2, nbins+1
      gloc(i) = gloc(i-1) + bin_width(i-1)
    end do
    fmin=1e-9
    fmax=1-fmin
    gloc(1) = fmin
    gloc(nbins+1)=fmax
    !write(*,*) "edges", gloc
    call gamma_inc_inv_calc(dust_mu_r8,gloc(:),diam_edges,nbins+1,"diam_edges")

    call calc_weights(dust_mu_r8,diam_edges, number, nbins) !fraction of total number between edges
    call calc_weights(dust_mu_r8+3,diam_edges, mass, nbins) !fraction of total mass between edges

!done bin calculation
!debug bins
!write(*,*) bin_width
!write(*,*) diam
!write(*,*) diam_edges
!    do i=1, nbins
!        write(*,*) "bin ", i,": ",diam_edges(i), diam(i), diam_edges(i+1), bin_width(i), number(i),mass(i)
!    end do



else if (dust_distribution == 1) then !lognormal
        sigma_squared = alog(veff+1)
        dust_cons1 = (DUST_PIRHOD / 6. )
        !exp(-sigma**2 / 2)
        dust_cons2 = exp(-sigma_squared/2)
        !Cunningham factors, B and E is not used in lognormal
        DUST_CONS_A_SEDIM = a_dust_sedim
        ln_dust_reff_max = alog(dust_reff_max)
        ln_dust_reff_min = alog(dust_reff_min)
        ln_dust_reff_default = alog(dust_reff_default)
endif
        DISTRIBUTION_MARS_DUST=dust_distribution

!ice constants
    ice_mu = dust_mu
    ice_cons2 = GAMMA(1.+DI+ice_mu) / (2.*GAMMA(ice_mu+DI))  !replaces dust_cons2 for ice on Mars
    LAMMINI = ICE_CONS2/100e-6 !100 micron
    LAMMAXI = ICE_CONS2/.1e-6 !.1 micron
!    LAMDEFI = ICE_CONS2/dust_reff_default !2.0 micron
        

    ! DI is the shape parameter for ice that scale size to mass
    ! alternatively, it is describes how density changes with diameter*(1/3), or makes
    ! ice particles behave less spherical if DI != 3.
    ice_cons10 = GAMMA(1+ice_mu)
    !ice_cons10min = GAMMA(1+ice_mu) / exp((1+ice_mu)*log(lammini)) ! gamma(mu+1) / lam**(mu+1)
    !ice_cons10max = GAMMA(1+ice_mu) / exp((1+ice_mu)*log(lammaxi)) ! gamma(mu+1) / lam**(mu+1)
    ice_cons11 = 1./(ci * GAMMA(1+ice_mu+DI))
    !ice_cons11min = ice_cons11 * exp((1+ice_mu+DI)*log(lammini))
    !ice_cons11max = ice_cons11 * exp((1+ice_mu+DI)*log(lammaxi))
    
    ice_cons12 = GAMMA(1.+DI+ice_mu) * ci / GAMMA(1+ice_mu)  !replaces cons12 for Mars

    ice_consmin_Q_N = exp(DI*log(lammini)) / ice_cons12
    ice_consmax_Q_N = exp(DI*log(lammaxi)) / ice_cons12

    homogeneous_ice_number_cons1 = (1./8) * (1./ci) * (1./homogeneous_ice_reff**3)  &
                                * GAMMA(ice_mu+1) * (ice_mu+DI)**3 / gamma(ice_mu+DI+1) !that's horrible
                                !For DI ==3:
                                !(6/8) * (1./pi*rho) * (1./reff**3) * (mu+3)*(mu+3)/((mu+2)*(mu+1))
                                !The '6' is already folded in 'ci'
    
    
    water_molec_volume = 3.22e-29 !m^3
    water_molec_radius = ((3./(4*pi)) * water_molec_volume)**(1./3)
    
    ice_cons_prefactor = sqrt(3./(8*pi*pi*water_molec_radius*rhoi)) * (1./kboltz)
    ice_cons_nuc_m0 = 1. !ice_cons10
    ice_cons_nuc_m2 = 1. ! GAMMA(3 + ice_mu)
    ice_cons_nuc_m3 = 1. ! GAMMA(DI + 1 + ice_mu)

    ice_cons_nuc_m0p5 = GAMMA(1.5 + ice_mu) / GAMMA(1. + ice_mu)  ! GAMMA(1.5+ ice_mu)
    ice_cons_nuc_m2p5 = GAMMA(3.5 + ice_mu) / GAMMA(3 + ice_mu)
    ice_cons_nuc_m3p5 = GAMMA(DI + 1.5 + ice_mu) / GAMMA(DI + 1. + ice_mu)
    
    photolysis_scrit = 10.
#endif

END SUBROUTINE MORR_TWO_MOMENT_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE MORR_TWO_MOMENT_MICRO) 
! WHICH OPERATES ON 1D VERTICAL COLUMNS.
! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
! BACK TO DRIVER MODEL USING THIS INTERFACE.
! MICROPHYSICS TENDENCIES ARE ADDED TO VARIABLES HERE BEFORE BEING PASSED BACK TO DRIVER MODEL.

! THIS CODE WAS WRITTEN BY HUGH MORRISON (NCAR) AND SLAVA TATARSKII (GEORGIA TECH).

! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MP_MORR_TWO_MOMENT(ITIMESTEP,                       &
                TH, QV, QC, QR, QI, QS, QG, NI, NS, NR, NG, &
                RHO, PII, P, DT_IN, DZ, HT, W,          &
                RAINNC, RAINNCV, SR,                    &
		SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV,    & ! hm added 7/13/13
                qrcuten, qscuten, qicuten, mu           & ! hm added
               ,F_QNDROP, qndrop                        & ! hm added, wrf-chem 
               ,IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
               ,IMS,IME, JMS,JME, KMS,KME               & ! memory dims
               ,ITS,ITE, JTS,JTE, KTS,KTE               & ! tile   dims            )
!jdf		   ,C2PREC3D,CSED3D,ISED3D,SSED3D,GSED3D,RSED3D & ! HM ADD, WRF-CHEM
		   ,QLSINK,PRECR,PRECI,PRECS,PRECG &        ! HM ADD, WRF-CHEM
#ifdef WRF_MARS
               ,qdust,ndust                              &
               ,dustnc, dustncv                          &
               ,qice_core                                &
               ,corenc, corencv                          &
               ,reff_dust,effi                           &
               ,do_water                                 &
               ,DUSTQ_MICRO, DUSTQ_SED                   &
               ,DUSTN_MICRO, DUSTN_SED                   &
               ,COREQ_MICRO, COREQ_SED                   &
               ,ICEQ_MICRO, ICEQ_SED                     &
               ,ICEN_MICRO, ICEN_SED                     &
               ,dsed1, dsed2                             &
#endif
                                            )
 
! QV - water vapor mixing ratio (kg/kg)
! QC - cloud water mixing ratio (kg/kg)
! QR - rain water mixing ratio (kg/kg)
! QI - cloud ice mixing ratio (kg/kg)
! QS - snow mixing ratio (kg/kg)
! QG - graupel mixing ratio (KG/KG)
! NI - cloud ice number concentration (1/kg)
! NS - Snow Number concentration (1/kg)
! NR - Rain Number concentration (1/kg)
! NG - Graupel number concentration (1/kg)
#ifdef WRF_MARS
! QDUST - dust mixing ratio (kg/kg)
! NDUST - dust Number concentration (1/kg)
! PRECD - SEDIMENTATION FLUXES (KG/M^2/S) FOR DUST
! QICE_CORE - dust core mixing ratio (kg/kg)
#endif
! NOTE: RHO AND HT NOT USED BY THIS SCHEME AND DO NOT NEED TO BE PASSED INTO SCHEME!!!!
! P - AIR PRESSURE (PA)
! W - VERTICAL AIR VELOCITY (M/S)
! TH - POTENTIAL TEMPERATURE (K)
! PII - exner function - used to convert potential temp to temp
! DZ - difference in height over interface (m)
! DT_IN - model time step (sec)
! ITIMESTEP - time step counter
! RAINNC - accumulated grid-scale precipitation (mm)
! RAINNCV - one time step grid scale precipitation (mm/time step)
! SNOWNC - accumulated grid-scale snow plus cloud ice (mm)
! SNOWNCV - one time step grid scale snow plus cloud ice (mm/time step)
! GRAUPELNC - accumulated grid-scale graupel (mm)
! GRAUPELNCV - one time step grid scale graupel (mm/time step)
! SR - one time step mass ratio of snow to total precip
! qrcuten, rain tendency from parameterized cumulus convection
! qscuten, snow tendency from parameterized cumulus convection
! qicuten, cloud ice tendency from parameterized cumulus convection

! variables below currently not in use, not coupled to PBL or radiation codes
! TKE - turbulence kinetic energy (m^2 s-2), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! NCTEND - droplet concentration tendency from pbl (kg-1 s-1)
! NCTEND - CLOUD ICE concentration tendency from pbl (kg-1 s-1)
! KZH - heat eddy diffusion coefficient from YSU scheme (M^2 S-1), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! HM, ADDED FOR WRF-CHEM COUPLING
! QLSINK - TENDENCY OF CLOUD WATER TO RAIN, SNOW, GRAUPEL (KG/KG/S)
! CSED,ISED,SSED,GSED,RSED - SEDIMENTATION FLUXES (KG/M^2/S) FOR CLOUD WATER, ICE, SNOW, GRAUPEL, RAIN
! PRECI,PRECS,PRECG,PRECR - SEDIMENTATION FLUXES (KG/M^2/S) FOR ICE, SNOW, GRAUPEL, RAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! reflectivity currently not included!!!!
! REFL_10CM - CALCULATED RADAR REFLECTIVITY AT 10 CM (DBZ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! EFFC - DROPLET EFFECTIVE RADIUS (MICRON)
! EFFR - RAIN EFFECTIVE RADIUS (MICRON)
! EFFS - SNOW EFFECTIVE RADIUS (MICRON)
! EFFI - CLOUD ICE EFFECTIVE RADIUS (MICRON)

! ADDITIONAL OUTPUT FROM MICRO - SEDIMENTATION TENDENCIES, NEEDED FOR LIQUID-ICE STATIC ENERGY

! QGSTEN - GRAUPEL SEDIMENTATION TEND (KG/KG/S)
! QRSTEN - RAIN SEDIMENTATION TEND (KG/KG/S)
! QISTEN - CLOUD ICE SEDIMENTATION TEND (KG/KG/S)
! QNISTEN - SNOW SEDIMENTATION TEND (KG/KG/S)
! QCSTEN - CLOUD WATER SEDIMENTATION TEND (KG/KG/S)

! WVAR - STANDARD DEVIATION OF SUB-GRID VERTICAL VELOCITY (M/S)

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde , &
                                       ims, ime, jms, jme, kms, kme , &
                                       its, ite, jts, jte, kts, kte
! Temporary changed from INOUT to IN
#ifdef WRF_PLANET
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          TH

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT), optional:: &
                          qv, qc, qr, qi, qs, qg, ni, ns, nr, NG   

#else
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          qv, qc, qr, qi, qs, qg, ni, ns, nr, TH, NG   

#endif

#ifdef WRF_MARS
   logical, intent(in),optional :: do_water
   logical :: use_water

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          qdust, ndust, qice_core

!   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), optional,INTENT(INOUT):: PRECD
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                                                              DUSTNC, DUSTNCV, &
                                                              CORENC, CORENCV

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), intent(out)::                     &
                    reff_dust, effi  !effective radius for dust

   real, dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: &
                DUSTQ_MICRO, DUSTQ_SED                  &
               ,DUSTN_MICRO, DUSTN_SED                  &
               ,COREQ_MICRO, COREQ_SED                  &
               ,ICEQ_MICRO, ICEQ_SED                    &
               ,ICEN_MICRO, ICEN_SED                     

   real, dimension(ims:ime, jms:jme), intent(inout), optional :: &
                                                        dsed1, dsed2
#endif
!jdf                      qndrop ! hm added, wrf-chem
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), optional,INTENT(INOUT):: qndrop
!jdf  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT):: CSED3D, &
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), optional,INTENT(INOUT):: QLSINK, &
                          PRECI,PRECS,PRECG,PRECR ! HM, WRF-CHEM
!, effcs, effis

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
                          pii, p, dz, rho, w !, tke, nctend, nitend,kzh
   REAL, INTENT(IN):: dt_in
   INTEGER, INTENT(IN):: ITIMESTEP

   REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                          RAINNC, RAINNCV, SR, &
! hm added 7/13/13
                          SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV

!   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &  ! GT
!                          refl_10cm

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::       ht

   ! LOCAL VARIABLES
#ifdef WRF_MARS
   REAL, DIMENSION(its:ite, kts:kte, jts:jte)::                     &
                      effs, effr, EFFG
#else
   REAL, DIMENSION(its:ite, kts:kte, jts:jte)::                     &
                      effi, effs, effr, EFFG
#endif

   REAL, DIMENSION(its:ite, kts:kte, jts:jte)::                     &
                      T, WVAR, EFFC

   REAL, DIMENSION(kts:kte) ::                                                                & 
                            QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D,                      &
                            NI_TEND1D, NS_TEND1D, NR_TEND1D,                                  &
                            QC1D, QI1D, QR1D,NI1D, NS1D, NR1D, QS1D,                          &
#ifdef WRF_MARS
                            QDUST1D, NDUST1D,                              & !concentrations
                            QCORE1D,                                       &
                            QDUST_TEND1D, NDUST_TEND1D, QDUSTSTEN,         & !tendencies
                            QCORE_TEND1D, QCORESTEN,                       &
                            EFFD1D,                                        & !effective radius
                            NDUSTSTEN, &
#endif
                            T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, W1D, WVAR1D,         &
                            EFFC1D, EFFI1D, EFFS1D, EFFR1D,DZ1D,   &
   ! HM ADD GRAUPEL
                            QG_TEND1D, NG_TEND1D, QG1D, NG1D, EFFG1D, &

! ADD SEDIMENTATION TENDENCIES (UNITS OF KG/KG/S)
                            QGSTEN,QRSTEN, QISTEN, QNISTEN, QCSTEN, NISTEN, &
! ADD CUMULUS TENDENCIES
                            QRCU1D, QSCU1D, QICU1D

! add cumulus tendencies

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
      qrcuten, qscuten, qicuten
   REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN):: &
      mu

  LOGICAL, INTENT(IN), OPTIONAL ::                F_QNDROP  ! wrf-chem
  LOGICAL :: flag_qndrop  ! wrf-chem
  integer :: iinum ! wrf-chem

! wrf-chem
   REAL, DIMENSION(kts:kte) :: nc1d, nc_tend1d,C2PREC,CSED,ISED,SSED,GSED,RSED    
! HM add reflectivity      
!                            dbz
                          
   REAL PRECPRT1D, SNOWRT1D, SNOWPRT1D, GRPLPRT1D ! hm added 7/13/13
#ifdef WRF_MARS
    REAL PRECDUST,PRECNDUST, PRECCORE, PRECDUST_ERR, PRECQVRT_ERR, PRECRT_PHOTOLYSIS
#endif
   INTEGER I,K,J

   REAL DT
#ifdef WRF_MARS
    use_water=.true.
    if(present(do_water)) use_water=do_water
#endif
!   LOGICAL:: dBZ_tstep ! GT
! below for wrf-chem
   flag_qndrop = .false.
   IF ( PRESENT ( f_qndrop ) ) flag_qndrop = f_qndrop
!!!!!!!!!!!!!!!!!!!!!!

   ! Initialize tendencies (all set to 0) and transfer
   ! array to local variables
   DT = DT_IN   

   DO I=ITS,ITE
   DO J=JTS,JTE
   DO K=KTS,KTE
       T(I,K,J)        = TH(i,k,j)*PII(i,k,j)

! NOTE: WVAR NOT CURRENTLY USED IN CODE !!!!!!!!!!
! currently assign wvar to 0.5 m/s (not coupled with PBL scheme)

       WVAR(I,K,J)     = 0.5

! currently mixing of number concentrations also is neglected (not coupled with PBL schemes)

   END DO
   END DO
   END DO

   do i=its,ite      ! i loop (east-west)
   do j=jts,jte      ! j loop (north-south)
   !
   ! Transfer 3D arrays into 1D for microphysical calculations
   !

! hm , initialize 1d tendency arrays to zero

      do k=kts,kte   ! k loop (vertical)

          QC_TEND1D(k)  = 0.
          QI_TEND1D(k)  = 0.
          QNI_TEND1D(k) = 0.
          QR_TEND1D(k)  = 0.
          NI_TEND1D(k)  = 0.
          NS_TEND1D(k)  = 0.
          NR_TEND1D(k)  = 0.
          T_TEND1D(k)   = 0.
          QV_TEND1D(k)  = 0.
          nc_tend1d(k) = 0. ! wrf-chem
#ifdef WRF_PLANET
#ifdef WRF_MARS
          if (present(qi)) then
            QI1D(k)       = QI(i,k,j)
          else
            QI1D(k) = 0.0
          endif
          if (present(ni)) then
            NI1D(k)       = NI(i,k,j)
          else
            NI1D(k) = 0.0
          endif

          qc1d(k) = 0.0
          qs1d(k) = 0.0
          qr1d(k) = 0.0
          ns1d(k) = 0.0
          nr1d(k) = 0.0
          qg1d(k) = 0.0
          ng1d(k) = 0.0
          QG_TEND1D(K)  = 0.
          NG_TEND1D(K)  = 0.
          EFFC1D(k)      = 0.0
          EFFI1D(k)      = 0.0
          EFFS1D(k)      = 0.0
          EFFR1D(k)      = 0.0
          EFFG1D(K)      = 0.0

#else
          CALL wrf_error_fatal('Two moment not set up for your chosen planet')    
#endif
#else

          QC1D(k)       = QC(i,k,j)
          QI1D(k)       = QI(i,k,j)
          QS1D(k)       = QS(i,k,j)
          QR1D(k)       = QR(i,k,j)

          NI1D(k)       = NI(i,k,j)

          NS1D(k)       = NS(i,k,j)
          NR1D(k)       = NR(i,k,j)
! HM ADD GRAUPEL
          QG1D(K)       = QG(I,K,j)
          NG1D(K)       = NG(I,K,j)
          QG_TEND1D(K)  = 0.
          NG_TEND1D(K)  = 0.
#endif

#ifdef WRF_MARS
          QDUST1D(k)    = QDUST(i,k,j)
          NDUST1D(k)    = NDUST(i,k,j)
          QCORE1D(K) = qice_core(i,k,j)
          QDUST_TEND1D(k) = 0.
          NDUST_TEND1D(k) = 0.
          QCORE_TEND1D(k) = 0.
          QDUSTSTEN(K) = 0.
          NDUSTSTEN(K) = 0.
          QCORESTEN(K) = 0.
#endif

          T1D(k)        = T(i,k,j)
#ifdef WRF_PLANET
#ifdef WRF_MARS
    if(present(qv)) then
          QV1D(k)       = QV(i,k,j)
    else
          QV1D(k)       = 0.0
    endif
#else
          CALL wrf_error_fatal('Two moment not set up for your chosen planet')    
#endif
#else
          QV1D(k)       = QV(i,k,j)
#endif
          P1D(k)        = P(i,k,j)
          DZ1D(k)       = DZ(i,k,j)
          W1D(k)        = W(i,k,j)
          WVAR1D(k)     = WVAR(i,k,j)
! add cumulus tendencies, decouple from mu
          qrcu1d(k)     = qrcuten(i,k,j)/mu(i,j)
          qscu1d(k)     = qscuten(i,k,j)/mu(i,j)
          qicu1d(k)     = qicuten(i,k,j)/mu(i,j)
      end do  !jdf added this
! below for wrf-chem
   IF (flag_qndrop .AND. PRESENT( qndrop )) THEN
      iact = 3
      DO k = kts, kte
         nc1d(k)=qndrop(i,k,j)
         iinum=0
      ENDDO
   ELSE
      DO k = kts, kte
         nc1d(k)=0. ! temporary placeholder, set to constant in microphysics subroutine
         iinum=1
      ENDDO
   ENDIF
#ifdef MARS_MORR_DEBUG
   jprint=.false.
   if(i.eq.IVALUE .and. j.eq.JVALUE) jprint =.true.
   if(jprint) then
      iprint=.true.
      CALL_DEBUG("start",0)
      CALL_DEBUG("KVALUE",KVALUE)
      CALL_DEBUG("qv1d_start",qv1d(KVALUE))
      CALL_DEBUG("qi1d_start",qi1d(KVALUE))
      CALL_DEBUG("ni1d_start",ni1d(KVALUE))
      CALL_DEBUG("QCORE1D_start",QCORE1D(KVALUE))
      CALL_DEBUG("QDUST1D_start",QDUST1D(KVALUE))
      CALL_DEBUG("ndust1d_start",ndust1d(KVALUE))
      CALL_DEBUG("effd1d", effd1d(KVALUE))
      CALL_DEBUG("effi1d", effi1d(KVALUE))
      iprint=.false.
    endif
#endif
      call MORR_TWO_MOMENT_MICRO(QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D,            &
       NI_TEND1D, NS_TEND1D, NR_TEND1D,                                                  &
       QC1D, QI1D, QS1D, QR1D,NI1D, NS1D, NR1D,                                          &
#ifdef WRF_MARS
       QDUST1D, NDUST1D, QCORE1D,                                                        &
       QDUST_TEND1D, NDUST_TEND1D, QDUSTSTEN, EFFD1D,                                    &
       QCORE_TEND1D, QCORESTEN,                                                          &
       PRECDUST, PRECNDUST, PRECCORE, NDUSTSTEN,  PRECDUST_ERR,PRECQVRT_ERR, PRECRT_PHOTOLYSIS,     &
#endif
       T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, DZ1D, W1D, WVAR1D,                            &
       PRECPRT1D,SNOWRT1D,                                                               &
       SNOWPRT1D,GRPLPRT1D,                 & ! hm added 7/13/13
       EFFC1D,EFFI1D,EFFS1D,EFFR1D,DT,                                                   &
                                            IMS,IME, JMS,JME, KMS,KME,                   &
                                            ITS,ITE, JTS,JTE, KTS,KTE,                   & ! HM ADD GRAUPEL
                                    QG_TEND1D,NG_TEND1D,QG1D,NG1D,EFFG1D, &
                                    qrcu1d, qscu1d, qicu1d, &
! ADD SEDIMENTATION TENDENCIES
                                  QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN, NISTEN, &
                                  nc1d, nc_tend1d, iinum, C2PREC,CSED,ISED,SSED,GSED,RSED & !wrf-chem
                                  ,i,j &
                       )
!write(*,*) __LINE__, rainncv(i,j)
   !
   ! Transfer 1D arrays back into 3D arrays
   !
!    write(*,*) sum(dz1d*(QI(i,:,j)+qv(i,:,j)) * p1d/(R*t1d)), sum(dz1d*(QI1d+qv1d) * p1d/(R*t1d))
#ifdef MARS_MORR_DEBUG    
    if(jprint) then
      iprint=.true.
      CALL_DEBUG("end",0)
      CALL_DEBUG("KVALUE",KVALUE)
      CALL_DEBUG("qv1d_end",qv1d(KVALUE))
      CALL_DEBUG("qi1d_end",qi1d(KVALUE))
      CALL_DEBUG("ni1d_end",ni1d(KVALUE))
      CALL_DEBUG("QCORE1D_end",QCORE1D(KVALUE))
      CALL_DEBUG("QDUST1D_end",QDUST1D(KVALUE))
      CALL_DEBUG("ndust1d_end",ndust1d(KVALUE))    
      CALL_DEBUG("effd1d", effd1d(KVALUE))
      CALL_DEBUG("effi1d", effi1d(KVALUE))
      CALL_DEBUG("QI_TEND", QI_TEND1D(KVALUE))
      CALL_DEBUG("net QI_TEND", QI_TEND1D(KVALUE)-QISTEN(KVALUE))
      CALL_DEBUG("NI_TEND", NI_TEND1D(KVALUE))
      CALL_DEBUG("net NI_TEND", NI_TEND1D(KVALUE)-NISTEN(KVALUE))
      CALL_DEBUG("QDUST_TEND", QDUST_TEND1D(KVALUE))
      CALL_DEBUG("NDUST_TEND", NDUST_TEND1D(KVALUE))
      CALL_DEBUG("QDUSTSTEN", QDUSTSTEN(KVALUE))
      CALL_DEBUG("QCORE_TEND", QCORE_TEND1D(KVALUE))
      CALL_DEBUG("QCORESTEN", QCORESTEN(KVALUE))
      CALL_DEBUG("QISTEN", QISTEN(KVALUE))
      CALL_DEBUG("NISTEN", NISTEN(KVALUE))
      CALL_DEBUG("QV_TEND1D", QV_TEND1D(KVALUE))
       
      iprint=.false.
    endif
#endif
      do k=kts,kte
!        write(*,*) __LINE__, k,QI(i,k,j)+qv(i,k,j), qi1d(k)+qv1d(k), ( qi1d(k)+qv1d(k)) - (QI(i,k,j)+qv(i,k,j))
!        if (( qi1d(k)+qv1d(k)) - (QI(i,k,j)+qv(i,k,j)) .ne. 0) then
!        write(*,*) "-------------------------"
!        endif
! hm, add tendencies to update global variables 
! HM, TENDENCIES FOR Q AND N NOW ADDED IN M2005MICRO, SO WE
! ONLY NEED TO TRANSFER 1D VARIABLES BACK TO 3D
#ifdef WRF_PLANET
#ifdef WRF_MARS
          if(present(qi)) then 
            QI(i,k,j)        = QI1D(k)
          endif
          if(present(ni)) then 
              NI(i,k,j)        = NI1D(k)
          endif
#else
          CALL wrf_error_fatal('Two moment not set up for your chosen planet')  
#endif
#else            
          QC(i,k,j)        = QC1D(k)
          QI(i,k,j)        = QI1D(k)
          QS(i,k,j)        = QS1D(k)
          QR(i,k,j)        = QR1D(k)
          NI(i,k,j)        = NI1D(k)
          NS(i,k,j)        = NS1D(k)          
          NR(i,k,j)        = NR1D(k)
	  QG(I,K,j)        = QG1D(K)
          NG(I,K,j)        = NG1D(K)
#endif

#ifdef WRF_MARS
          QDUST(i,k,j) = QDUST1D(k)
          NDUST(i,k,j) = NDUST1D(k)
          qice_core(i,k,j) = QCORE1D(K)
	      reff_dust(i,k,j) = EFFD1D(k)
#endif
          T(i,k,j)         = T1D(k)
          TH(I,K,J)        = T(i,k,j)/PII(i,k,j) ! CONVERT TEMP BACK TO POTENTIAL TEMP
#ifdef WRF_PLANET
#ifdef WRF_MARS
    if(present(qv)) then
          QV(i,k,j)       = QV1D(k)
    endif
#else
          CALL wrf_error_fatal('Two moment not set up for your chosen planet')    
#endif
#else
          QV(i,k,j)        = QV1D(k)
#endif
!          EFFC(i,k,j)      = EFFC1D(k)
          EFFI(i,k,j)      = EFFI1D(k)
!          EFFS(i,k,j)      = EFFS1D(k)
!          EFFR(i,k,j)      = EFFR1D(k)
!	  EFFG(I,K,j)      = EFFG1D(K)

! wrf-chem
          IF (flag_qndrop .AND. PRESENT( qndrop )) THEN
             qndrop(i,k,j) = nc1d(k)
!jdf         CSED3D(I,K,J) = CSED(K)
          END IF
#ifdef WRF_PLANET
          IF ( PRESENT( QLSINK ) .and. present(QC)) THEN
#else
          IF ( PRESENT( QLSINK ) ) THEN
#endif
             if(qc(i,k,j)>1.e-10) then
                QLSINK(I,K,J)  = C2PREC(K)/QC(I,K,J)
             else
                QLSINK(I,K,J)  = 0.0
             endif
          END IF
          IF ( PRESENT( PRECR ) ) PRECR(I,K,J) = RSED(K)
          IF ( PRESENT( PRECI ) ) PRECI(I,K,J) = ISED(K)
          IF ( PRESENT( PRECS ) ) PRECS(I,K,J) = SSED(K)
          IF ( PRESENT( PRECG ) ) PRECG(I,K,J) = GSED(K)

! EFFECTIVE RADIUS FOR RADIATION CODE (currently not coupled)
! HM, ADD LIMIT TO PREVENT BLOWING UP OPTICAL PROPERTIES, 8/18/07
!          EFFCS(I,K,J)     = MIN(EFFC(I,K,J),50.)
!          EFFCS(I,K,J)     = MAX(EFFCS(I,K,J),1.)
!          EFFIS(I,K,J)     = MIN(EFFI(I,K,J),130.)
!          EFFIS(I,K,J)     = MAX(EFFIS(I,K,J),13.)

      end do

! hm modified so that m2005 precip variables correctly match wrf precip variables
      RAINNC(i,j) = RAINNC(I,J)+PRECPRT1D
      RAINNCV(i,j) = PRECPRT1D
! hm, added 7/13/13
      SNOWNC(i,j) = SNOWNC(I,J)+SNOWPRT1D
      SNOWNCV(i,j) = SNOWPRT1D
      GRAUPELNC(i,j) = GRAUPELNC(I,J)+GRPLPRT1D
      GRAUPELNCV(i,j) = GRPLPRT1D
      SR(i,j) = SNOWRT1D/(PRECPRT1D+1.E-12)

#ifdef WRF_MARS
      DUSTNC(i,j) = DUSTNC(I,J)+PRECDUST
      DUSTNCV(i,j) = PRECDUST
      CORENC(i,j) = CORENC(I,J)+PRECCORE
      CORENCV(i,j) = PRECCORE
!      write(*,*) i,j,dustnc(i,j), dustncv(i,j), corenc(i,j), corencv(i,j)

    !diagnostics
      DUSTQ_MICRO(i,kts:kte,j) = DUSTQ_MICRO(i,kts:kte,j) + QDUST_TEND1D(kts:kte)
      DUSTQ_SED(i,kts:kte,j)   = DUSTQ_SED(i,kts:kte,j)   + QDUSTSTEN(kts:kte)
      DUSTN_MICRO(i,kts:kte,j) = DUSTN_MICRO(i,kts:kte,j) + NDUST_TEND1D(kts:kte)
      DUSTN_SED(i,kts:kte,j)   = DUSTN_SED(i,kts:kte,j)   + NDUSTSTEN(kts:kte)
      COREQ_MICRO(i,kts:kte,j) = COREQ_MICRO(i,kts:kte,j) + QCORE_TEND1D(kts:kte)
      COREQ_SED(i,kts:kte,j)   = COREQ_SED(i,kts:kte,j)   + QCORESTEN(kts:kte)
      ICEQ_MICRO(i,kts:kte,j)  = ICEQ_MICRO(i,kts:kte,j)  + QI_TEND1D(kts:kte)
      ICEQ_SED(i,kts:kte,j)    = ICEQ_SED(i,kts:kte,j)    + QISTEN(kts:kte)
      ICEN_MICRO(i,kts:kte,j)  = ICEN_MICRO(i,kts:kte,j)  + NI_TEND1D(kts:kte)
      ICEN_SED(i,kts:kte,j)    = ICEN_SED(i,kts:kte,j)    + NISTEN(kts:kte)
      dsed1(i,j) = max(precdust/dt,0.)
      dsed2(i,j) = max(precndust/dt,0.)
      !csed1(i,j) = csed1(i,j) + preccore
       
#endif

   end do
   end do   

END SUBROUTINE MP_MORR_TWO_MOMENT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MORR_TWO_MOMENT_MICRO(QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,         &
       NI3DTEN,NS3DTEN,NR3DTEN,QC3D,QI3D,QNI3D,QR3D,NI3D,NS3D,NR3D,              &
#ifdef WRF_MARS
       qdust3d,ndust3d, qcore3d,                                                 &
!tendency                                                                 
       QDUST3DTEN, NDUST3DTEN, QDUSTSTEN, EFFD,                                  &
       QCORE3DTEN, QCORESTEN,                                                    &
       PRECDUSTRT, PRECNDUSTRT, PRECCORERT,  NDUSTSTEN, PRECDUSTRT_ERR, PRECQVRT_ERR,         &
       PRECRT_PHOTOLYSIS, &
#endif
       T3DTEN,QV3DTEN,T3D,QV3D,PRES,DZQ,W3D,WVAR,PRECRT,SNOWRT,            &
       SNOWPRT,GRPLPRT,                & ! hm added 7/13/13
       EFFC,EFFI,EFFS,EFFR,DT,                                                   &
                                            IMS,IME, JMS,JME, KMS,KME,           &
                                            ITS,ITE, JTS,JTE, KTS,KTE,           & ! ADD GRAUPEL
                        QG3DTEN,NG3DTEN,QG3D,NG3D,EFFG,qrcu1d,qscu1d, qicu1d,    &
                        QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN,NISTEN, &
                        nc3d,nc3dten,iinum, & ! wrf-chem
				c2prec,CSED,ISED,SSED,GSED,RSED  &  ! hm added, wrf-chem
        ,thisi, thisj &
                        )

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THIS PROGRAM IS THE MAIN TWO-MOMENT MICROPHYSICS SUBROUTINE DESCRIBED BY
! MORRISON ET AL. 2005 JAS AND MORRISON ET AL. 2009 MWR

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL/HAIL.

! CODE STRUCTURE: MAIN SUBROUTINE IS 'MORR_TWO_MOMENT'. ALSO INCLUDED IN THIS FILE IS
! 'FUNCTION POLYSVP', 'FUNCTION DERF1', AND
! 'FUNCTION GAMMA'.

! NOTE: THIS SUBROUTINE USES 1D ARRAY IN VERTICAL (COLUMN), EVEN THOUGH VARIABLES ARE CALLED '3D'......

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! DECLARATIONS

      IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THESE VARIABLES BELOW MUST BE LINKED WITH THE MAIN MODEL.
! DEFINE ARRAY SIZES

! INPUT NUMBER OF GRID CELLS
  integer, intent(in) :: thisi, thisj
! INPUT/OUTPUT PARAMETERS                                 ! DESCRIPTION (UNITS)
      INTEGER, INTENT( IN)  :: IMS,IME, JMS,JME, KMS,KME,          &
                               ITS,ITE, JTS,JTE, KTS,KTE

      REAL, DIMENSION(KTS:KTE) ::  QC3DTEN            ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QI3DTEN            ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QNI3DTEN           ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QR3DTEN            ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  NI3DTEN            ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  NS3DTEN            ! SNOW NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  NR3DTEN            ! RAIN NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QC3D               ! CLOUD WATER MIXING RATIO (KG/KG)
      REAL, DIMENSION(KTS:KTE) ::  QI3D               ! CLOUD ICE MIXING RATIO (KG/KG)
      REAL, DIMENSION(KTS:KTE) ::  QNI3D              ! SNOW MIXING RATIO (KG/KG)
      REAL, DIMENSION(KTS:KTE) ::  QR3D               ! RAIN MIXING RATIO (KG/KG)
      REAL, DIMENSION(KTS:KTE) ::  NI3D               ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KTS:KTE) ::  NS3D               ! SNOW NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KTS:KTE) ::  NR3D               ! RAIN NUMBER CONCENTRATION (1/KG)
#ifdef WRF_MARS
      REAL, DIMENSION(KTS:KTE) ::  QDUST3D            ! DUST MIXING RATIO (KG/KG)
      REAL, DIMENSION(KTS:KTE) ::  NDUST3D            ! DUST NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KTS:KTE) ::  QDUST3DTEN         ! DUST MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  NDUST3DTEN         ! DUST NUMBER CONCENTRATION TENDENCY (1/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QDUSTSTEN          ! DUST SED TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QCORE3D            ! DUST CORE MIXING RATIO (KG/KG)
      REAL, DIMENSION(KTS:KTE) ::  QCORE3DTEN         ! DUST CORE TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QCORESTEN          ! DUST CORE SED TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  NDUSTSTEN             ! DUST NUMBER SED TENDENCY (/KG/S)
#endif
      REAL, DIMENSION(KTS:KTE) ::  T3DTEN             ! TEMPERATURE TENDENCY (K/S)
      REAL, DIMENSION(KTS:KTE) ::  QV3DTEN            ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  T3D                ! TEMPERATURE (K)
      REAL, DIMENSION(KTS:KTE) ::  QV3D               ! WATER VAPOR MIXING RATIO (KG/KG)
      REAL, DIMENSION(KTS:KTE) ::  PRES               ! ATMOSPHERIC PRESSURE (PA)
      REAL, DIMENSION(KTS:KTE) ::  DZQ                ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
      REAL, DIMENSION(KTS:KTE) ::  W3D                ! GRID-SCALE VERTICAL VELOCITY (M/S)
      REAL, DIMENSION(KTS:KTE) ::  WVAR               ! SUB-GRID VERTICAL VELOCITY (M/S)
! below for wrf-chem
      REAL, DIMENSION(KTS:KTE) ::  nc3d
      REAL, DIMENSION(KTS:KTE) ::  nc3dten
      integer, intent(in) :: iinum

! HM ADDED GRAUPEL VARIABLES
      REAL, DIMENSION(KTS:KTE) ::  QG3DTEN            ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  NG3DTEN            ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QG3D            ! GRAUPEL MIX RATIO (KG/KG)
      REAL, DIMENSION(KTS:KTE) ::  NG3D            ! GRAUPEL NUMBER CONC (1/KG)

! HM, ADD 1/16/07, SEDIMENTATION TENDENCIES FOR MIXING RATIO

      REAL, DIMENSION(KTS:KTE) ::  QGSTEN            ! GRAUPEL SED TEND (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QRSTEN            ! RAIN SED TEND (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QISTEN            ! CLOUD ICE SED TEND (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QNISTEN           ! SNOW SED TEND (KG/KG/S)
      REAL, DIMENSION(KTS:KTE) ::  QCSTEN            ! CLOUD WAT SED TEND (KG/KG/S)      
      REAL, DIMENSION(KTS:KTE) ::  NISTEN             ! ICE NUMBER SED TENDENCY (/KG/S)

! hm add cumulus tendencies for precip
        REAL, DIMENSION(KTS:KTE) ::   qrcu1d
        REAL, DIMENSION(KTS:KTE) ::   qscu1d
        REAL, DIMENSION(KTS:KTE) ::   qicu1d

! OUTPUT VARIABLES

        REAL, intent(out) :: PRECRT                ! TOTAL PRECIP PER TIME STEP (mm)
        REAL, intent(out) ::  SNOWRT                ! SNOW PER TIME STEP (mm)
! hm added 7/13/13
        REAL, intent(out) ::  SNOWPRT      ! TOTAL CLOUD ICE PLUS SNOW PER TIME STEP (mm)
	REAL, intent(out) ::  GRPLPRT	  ! TOTAL GRAUPEL PER TIME STEP (mm)

        REAL, DIMENSION(KTS:KTE) ::   EFFC            ! DROPLET EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KTS:KTE) ::   EFFI            ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KTS:KTE) ::   EFFS            ! SNOW EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KTS:KTE) ::   EFFR            ! RAIN EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KTS:KTE) ::   EFFG            ! GRAUPEL EFFECTIVE RADIUS (MICRON)
#ifdef WRF_MARS
        REAL, DIMENSION(KTS:KTE) ::   EFFD            ! DUST EFFECTIVE RADIUS (MICRON)
        REAL, intent(out) :: PRECDUSTRT, PRECNDUSTRT, PRECCORERT,PRECDUSTRT_ERR, PRECQVRT_ERR,PRECRT_PHOTOLYSIS

#endif

! MODEL INPUT PARAMETERS (FORMERLY IN COMMON BLOCKS)

        REAL DT         ! MODEL TIME STEP (SEC)
#ifdef WRF_MARS
    real :: core_mass_rate, core_number_rate, ice_mass_rate
    logical :: naked_core_removal
#endif

!.....................................................................................................
! LOCAL VARIABLES: ALL PARAMETERS BELOW ARE LOCAL TO SCHEME AND DON'T NEED TO COMMUNICATE WITH THE
! REST OF THE MODEL.

! SIZE PARAMETER VARIABLES

     REAL, DIMENSION(KTS:KTE) :: LAMC          ! SLOPE PARAMETER FOR DROPLETS (M-1)
     REAL, DIMENSION(KTS:KTE) :: LAMI          ! SLOPE PARAMETER FOR CLOUD ICE (M-1)
     REAL, DIMENSION(KTS:KTE) :: LAMS          ! SLOPE PARAMETER FOR SNOW (M-1)
     REAL, DIMENSION(KTS:KTE) :: LAMR          ! SLOPE PARAMETER FOR RAIN (M-1)
     REAL, DIMENSION(KTS:KTE) :: LAMG          ! SLOPE PARAMETER FOR GRAUPEL (M-1)
     REAL, DIMENSION(KTS:KTE) :: CDIST1        ! PSD PARAMETER FOR DROPLETS
     REAL, DIMENSION(KTS:KTE) :: N0I           ! INTERCEPT PARAMETER FOR CLOUD ICE (KG-1 M-1)
     REAL, DIMENSION(KTS:KTE) :: N0S           ! INTERCEPT PARAMETER FOR SNOW (KG-1 M-1)
     REAL, DIMENSION(KTS:KTE) :: N0RR          ! INTERCEPT PARAMETER FOR RAIN (KG-1 M-1)
     REAL, DIMENSION(KTS:KTE) :: N0G           ! INTERCEPT PARAMETER FOR GRAUPEL (KG-1 M-1)
     REAL, DIMENSION(KTS:KTE) :: PGAM          ! SPECTRAL SHAPE PARAMETER FOR DROPLETS

! MICROPHYSICAL PROCESSES

     REAL, DIMENSION(KTS:KTE) ::  NSUBC     ! LOSS OF NC DURING EVAP
     REAL, DIMENSION(KTS:KTE) ::  NSUBI     ! LOSS OF NI DURING SUB.
     REAL, DIMENSION(KTS:KTE) ::  NSUBS     ! LOSS OF NS DURING SUB.
     REAL, DIMENSION(KTS:KTE) ::  NSUBR     ! LOSS OF NR DURING EVAP
     REAL, DIMENSION(KTS:KTE) ::  PRD       ! DEP CLOUD ICE
     REAL, DIMENSION(KTS:KTE) ::  PRE       ! EVAP OF RAIN
     REAL, DIMENSION(KTS:KTE) ::  PRDS      ! DEP SNOW
     REAL, DIMENSION(KTS:KTE) ::  NNUCCC    ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
     REAL, DIMENSION(KTS:KTE) ::  MNUCCC    ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
     REAL, DIMENSION(KTS:KTE) ::  PRA       ! ACCRETION DROPLETS BY RAIN
     REAL, DIMENSION(KTS:KTE) ::  PRC       ! AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KTS:KTE) ::  PCC       ! COND/EVAP DROPLETS
     REAL, DIMENSION(KTS:KTE) ::  NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL, DIMENSION(KTS:KTE) ::  MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL, DIMENSION(KTS:KTE) ::  MNUCCR    ! CHANGE Q DUE TO CONTACT FREEZ RAIN
     REAL, DIMENSION(KTS:KTE) ::  NNUCCR    ! CHANGE N DUE TO CONTACT FREEZ RAIN
     REAL, DIMENSION(KTS:KTE) ::  NPRA      ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
     REAL, DIMENSION(KTS:KTE) ::  NRAGG     ! SELF-COLLECTION/BREAKUP OF RAIN
     REAL, DIMENSION(KTS:KTE) ::  NSAGG     ! SELF-COLLECTION OF SNOW
     REAL, DIMENSION(KTS:KTE) ::  NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KTS:KTE) ::  NPRC1      ! CHANGE NR AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KTS:KTE) ::  PRAI      ! CHANGE Q ACCRETION CLOUD ICE BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  PRCI      ! CHANGE Q AUTOCONVERSIN CLOUD ICE TO SNOW
     REAL, DIMENSION(KTS:KTE) ::  PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
     REAL, DIMENSION(KTS:KTE) ::  NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
     REAL, DIMENSION(KTS:KTE) ::  NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  NPRAI     ! CHANGE N ACCRETION CLOUD ICE
     REAL, DIMENSION(KTS:KTE) ::  NMULTS    ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  NMULTR    ! ICE MULT DUE TO RIMING RAIN BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  QMULTS    ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
     REAL, DIMENSION(KTS:KTE) ::  QMULTR    ! CHANGE Q DUE TO ICE RAIN/SNOW
     REAL, DIMENSION(KTS:KTE) ::  PRACS     ! CHANGE Q RAIN-SNOW COLLECTION
     REAL, DIMENSION(KTS:KTE) ::  NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
     REAL, DIMENSION(KTS:KTE) ::  PCCN      ! CHANGE Q DROPLET ACTIVATION
     REAL, DIMENSION(KTS:KTE) ::  PSMLT     ! CHANGE Q MELTING SNOW TO RAIN
     REAL, DIMENSION(KTS:KTE) ::  EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
     REAL, DIMENSION(KTS:KTE) ::  NSMLTS    ! CHANGE N MELTING SNOW
     REAL, DIMENSION(KTS:KTE) ::  NSMLTR    ! CHANGE N MELTING SNOW TO RAIN
! HM ADDED 12/13/06
     REAL, DIMENSION(KTS:KTE) ::  PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
     REAL, DIMENSION(KTS:KTE) ::  NIACR     ! CHANGE N, ICE-RAIN COLLECTION
     REAL, DIMENSION(KTS:KTE) ::  PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
     REAL, DIMENSION(KTS:KTE) ::  PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KTS:KTE) ::  NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KTS:KTE) ::  PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KTS:KTE) ::  EPRD      ! SUBLIMATION CLOUD ICE
     REAL, DIMENSION(KTS:KTE) ::  EPRDS     ! SUBLIMATION SNOW
! HM ADDED GRAUPEL PROCESSES
     REAL, DIMENSION(KTS:KTE) ::  PRACG    ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  PSACWG    ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  PGSACW    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  PGRACS    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  PRDG    ! DEP OF GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  EPRDG    ! SUB OF GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
     REAL, DIMENSION(KTS:KTE) ::  PGMLT    ! CHANGE Q MELTING OF GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  NSCNG    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  NGRACS    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL, DIMENSION(KTS:KTE) ::  NGMLTG    ! CHANGE N MELTING GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  NGMLTR    ! CHANGE N MELTING GRAUPEL TO RAIN
     REAL, DIMENSION(KTS:KTE) ::  NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  PSACR    ! CONVERSION DUE TO COLL OF SNOW BY RAIN
     REAL, DIMENSION(KTS:KTE) ::  NMULTG    ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  NMULTRG    ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  QMULTG    ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
     REAL, DIMENSION(KTS:KTE) ::  QMULTRG    ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL

! TIME-VARYING ATMOSPHERIC PARAMETERS

     REAL, DIMENSION(KTS:KTE) ::   KAP   ! THERMAL CONDUCTIVITY OF AIR
     REAL, DIMENSION(KTS:KTE) ::   EVS   ! SATURATION VAPOR PRESSURE
     REAL, DIMENSION(KTS:KTE) ::   EIS   ! ICE SATURATION VAPOR PRESSURE
     REAL, DIMENSION(KTS:KTE) ::   QVS   ! SATURATION MIXING RATIO
     REAL, DIMENSION(KTS:KTE) ::   QVI   ! ICE SATURATION MIXING RATIO
     REAL, DIMENSION(KTS:KTE) ::   QVQVS ! SAUTRATION RATIO
     REAL, DIMENSION(KTS:KTE) ::   QVQVSI! ICE SATURAION RATIO
     REAL, DIMENSION(KTS:KTE) ::   DV    ! DIFFUSIVITY OF WATER VAPOR IN AIR
     REAL, DIMENSION(KTS:KTE) ::   XXLS  ! LATENT HEAT OF SUBLIMATION
     REAL, DIMENSION(KTS:KTE) ::   XXLV  ! LATENT HEAT OF VAPORIZATION
     REAL, DIMENSION(KTS:KTE) ::   CPM   ! SPECIFIC HEAT AT CONST PRESSURE FOR MOIST AIR
     REAL, DIMENSION(KTS:KTE) ::   MU    ! VISCOCITY OF AIR
     REAL, DIMENSION(KTS:KTE) ::   SC    ! SCHMIDT NUMBER
     REAL, DIMENSION(KTS:KTE) ::   XLF   ! LATENT HEAT OF FREEZING
     REAL, DIMENSION(KTS:KTE) ::   RHO   ! AIR DENSITY
     REAL, DIMENSION(KTS:KTE) ::   AB    ! CORRECTION TO CONDENSATION RATE DUE TO LATENT HEATING
     REAL, DIMENSION(KTS:KTE) ::   ABI    ! CORRECTION TO DEPOSITION RATE DUE TO LATENT HEATING

! TIME-VARYING MICROPHYSICS PARAMETERS

     REAL, DIMENSION(KTS:KTE) ::   DAP    ! DIFFUSIVITY OF AEROSOL
     REAL    NACNT                    ! NUMBER OF CONTACT IN
     REAL    FMULT                    ! TEMP.-DEP. PARAMETER FOR RIME-SPLINTERING
     REAL    COFFI                    ! ICE AUTOCONVERSION PARAMETER

! FALL SPEED WORKING VARIABLES (DEFINED IN CODE)

      REAL, DIMENSION(KTS:KTE) ::    DUMI,DUMR,DUMFNI,DUMG,DUMFNG
      REAL UNI, UMI,UMR
      REAL, DIMENSION(KTS:KTE) ::    FR, FI, FNI,FG,FNG
      REAL RGVM
      REAL, DIMENSION(KTS:KTE) ::   FALOUTR,FALOUTI,FALOUTNI
      REAL FALTNDR,FALTNDI,FALTNDNI,RHO2
      REAL, DIMENSION(KTS:KTE) ::   DUMQS,DUMFNS
      REAL UMS,UNS
      REAL, DIMENSION(KTS:KTE) ::   FS,FNS, FALOUTS,FALOUTNS,FALOUTG,FALOUTNG
      REAL FALTNDS,FALTNDNS,UNR,FALTNDG,FALTNDNG
      REAL, DIMENSION(KTS:KTE) ::    DUMC,DUMFNC
      REAL UNC,UMC,UNG,UMG
      REAL, DIMENSION(KTS:KTE) ::   FC,FALOUTC,FALOUTNC
      REAL FALTNDC,FALTNDNC
      REAL, DIMENSION(KTS:KTE) ::   FNC,DUMFNR,FALOUTNR
      REAL FALTNDNR
      REAL, DIMENSION(KTS:KTE) ::   FNR

#ifdef WRF_MARS
    REAL, DIMENSION(KTS:KTE) :: N0D
    REAL, DIMENSION(KTS:KTE) :: LAMD
    REAL, DIMENSION(KTS:KTE) :: DUMD, DUMFND
    REAL :: UMD, UND
    REAL, DIMENSION(KTS:KTE) :: FND, FD
    REAL :: FALTNDD, FALTNDND
    REAL, DIMENSION(KTS:KTE) :: FALOUTD, FALOUTND
    REAL, DIMENSION(KTS:KTE) :: FALOUTNDND
    
    !dust cores
    REAL :: UMCORE
    REAL, DIMENSION(KTS:KTE) :: FCORE, FALOUTCORE, FALOUTNDCORE
    REAL, DIMENSION(KTS:KTE) :: MNUCCORE
    REAL :: FALTNDCORE
    REAL, DIMENSION(KTS:KTE) :: DUMCORE
    
    REAL, DIMENSION(KTS:KTE) :: QSUBCORE, QSUBDUST, NSUBDUST
    
    REAL, DIMENSION(KTS:KTE) :: ice_rho, qicetot3d
    real :: n0i_temp, ni3d_temp, scrit
    real :: scrit_nuc, scrit_cond
    real :: dum_nuc, dum_cond
#endif
! FALL-SPEED PARAMETER 'A' WITH AIR DENSITY CORRECTION

      REAL, DIMENSION(KTS:KTE) ::    AIN,ARN,ASN,ACN,AGN
#ifdef WRF_MARS
      REAL, dimension(KTS:KTE) :: ADN, freepath
#endif      

! EXTERNAL FUNCTION CALL RETURN VARIABLES

!      REAL GAMMA,      ! EULER GAMMA FUNCTION
!      REAL POLYSVP,    ! SAT. PRESSURE FUNCTION
!      REAL DERF1        ! ERROR FUNCTION

! DUMMY VARIABLES

     REAL DUM,DUM1,DUM2,DUMT,DUMQV,DUMQSS,DUMQSI,DUMS

! PROGNOSTIC SUPERSATURATION

     REAL DQSDT    ! CHANGE OF SAT. MIX. RAT. WITH TEMPERATURE
     REAL DQSIDT   ! CHANGE IN ICE SAT. MIXING RAT. WITH T
     REAL EPSI     ! 1/PHASE REL. TIME (SEE M2005), ICE
     REAL EPSS     ! 1/PHASE REL. TIME (SEE M2005), SNOW
     REAL EPSR     ! 1/PHASE REL. TIME (SEE M2005), RAIN
     REAL EPSG     ! 1/PHASE REL. TIME (SEE M2005), GRAUPEL

! NEW DROPLET ACTIVATION VARIABLES
     REAL TAUC     ! PHASE REL. TIME (SEE M2005), DROPLETS
     REAL TAUR     ! PHASE REL. TIME (SEE M2005), RAIN
     REAL TAUI     ! PHASE REL. TIME (SEE M2005), CLOUD ICE
     REAL TAUS     ! PHASE REL. TIME (SEE M2005), SNOW
     REAL TAUG     ! PHASE REL. TIME (SEE M2005), GRAUPEL
     REAL DUMACT,DUM3

! COUNTING/INDEX VARIABLES

     INTEGER K,NSTEP,N ! ,I

! LTRUE IS ONLY USED TO SPEED UP THE CODE !!
! LTRUE, SWITCH = 0, NO HYDROMETEORS IN COLUMN, 
!               = 1, HYDROMETEORS IN COLUMN

      INTEGER LTRUE

! DROPLET ACTIVATION/FREEZING AEROSOL


     REAL    CT      ! DROPLET ACTIVATION PARAMETER
     REAL    TEMP1   ! DUMMY TEMPERATURE
     REAL    SAT1    ! DUMMY SATURATION
     REAL    SIGVL   ! SURFACE TENSION LIQ/VAPOR
     REAL    KEL     ! KELVIN PARAMETER
     REAL    KC2     ! TOTAL ICE NUCLEATION RATE

       REAL CRY,KRY   ! AEROSOL ACTIVATION PARAMETERS

! MORE WORKING/DUMMY VARIABLES

     REAL DUMQI,DUMNI,DC0,DS0,DG0
     REAL DUMQC,DUMQR,RATIO,SUM_DEP,FUDGEF

! EFFECTIVE VERTICAL VELOCITY  (M/S)
     REAL WEF

! WORKING PARAMETERS FOR ICE NUCLEATION

      REAL ANUC,BNUC

! WORKING PARAMETERS FOR AEROSOL ACTIVATION

        REAL AACT,GAMM,GG,PSI,ETA1,ETA2,SM1,SM2,SMAX,UU1,UU2,ALPHA

! DUMMY SIZE DISTRIBUTION PARAMETERS

        REAL DLAMS,DLAMR,DLAMI,DLAMC,DLAMG,LAMMAX,LAMMIN
#ifdef WRF_MARS
		REAL DLAMD, invDLAMD, invDLAMI, tf, invtf, se_const, se_m, se_n
		real :: nuc_Fd, nuc_Fh, nuc_ka
    real :: ratio_dust, ratio_water, dum_dust, sum_dep_dust
#endif
        INTEGER IDROP

! FOR WRF-CHEM
	REAL, DIMENSION(KTS:KTE)::C2PREC,CSED,ISED,SSED,GSED,RSED 

	logical :: skip_micro
! comment lines for wrf-chem since these are intent(in) in that case
!       REAL, DIMENSION(KTS:KTE) ::  NC3DTEN            ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S)
!       REAL, DIMENSION(KTS:KTE) ::  NC3D               ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)
#ifdef WRF_MARS
    real :: water_before, water_after, frac
    real :: dust_before, dust_after, core_before, core_after
    real, DIMENSION(KTS:KTE) :: old_qv3d, old_qi3d, old_dust, old_core

    old_qv3d = qv3d
    old_qi3d = QI3D
    old_dust = qdust3d
    old_core = qcore3d
    qicetot3d = 0.
#endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! INITIALIZE PRECIP AND SNOW RATES
      PRECRT = 0.
      SNOWRT = 0.
! hm added 7/13/13
      SNOWPRT = 0.
      GRPLPRT = 0.
#ifdef WRF_MARS
      PRECDUSTRT = 0
      PRECNDUSTRT = 0
      PRECCORERT = 0
      PRECDUSTRT_ERR = 0
      PRECQVRT_ERR = 0
      PRECRT_PHOTOLYSIS = 0.
#endif

! SET LTRUE INITIALLY TO 0

         LTRUE = 0

! ATMOSPHERIC PARAMETERS THAT VARY IN TIME AND HEIGHT
         DO K = KTS,KTE
          iprint=.false.
          if(jprint.and.k.eq.KVALUE) iprint=.true.
! NC3DTEN LOCAL ARRAY INITIALIZED
               NC3DTEN(K) = 0.
! INITIALIZE VARIABLES FOR WRF-CHEM OUTPUT TO ZERO

		C2PREC(K)=0.
		CSED(K)=0.
		ISED(K)=0.
		SSED(K)=0.
		GSED(K)=0.
		RSED(K)=0.
#ifdef WRF_MARS
! 		DSED(k)=0.
#endif

#ifdef WRF_MARS
    XXLV(K) = 0.0
    !XXLS(K) = 0.0
!    XXLS(K) = 3.15E6-2370.*T3D(K)+0.3337E6
!updated from  R. R. Rogers; M. K. Yau (1989)
!equivalent to (2614869.81 +1895.2*temp - 4.*temp**2)
    XXLS(K) = (2834.1 - 0.29*(T3D(K)-273.15) - 0.004*(T3D(K)-273.15)*(T3D(K)-273.15))*1e3 !J/kg
    


#else
! LATENT HEAT OF VAPORATION

            XXLV(K) = 3.1484E6-2370.*T3D(K)

! LATENT HEAT OF SUBLIMATION

            XXLS(K) = 3.15E6-2370.*T3D(K)+0.3337E6

#endif
#ifdef WRF_MARS
    !cl copying from Michael's code elsewhere
           CPM(K) = CP* (1.0 + ((1800./CP)-1.) * QV3D(K))     !c-mm J/(K KG)   Generalized for all planets/values of CPD
#else
      
            CPM(K) = CP*(1.+0.887*QV3D(K))
#endif
! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
#ifdef WRF_MARS
            EVS(K) = esat(T3D(K)) ! PA
            EIS(K) = EVS(K)       ! PA
    !total fraction for mars, not dry fraction to 
    !prevent extreme values where water vapor pressure equals or exceeds pressure
            QVS(K) = EP_2*EVS(K)/(PRES(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K))
#else
            EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
            EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))
#endif

            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

! AIR DENSITY

            RHO(K) = PRES(K)/(R*T3D(K))
#ifdef MARS_MORR_DEBUG
            if (iprint) then
                CALL_DEBUG("PRES",pres(k))
                CALL_DEBUG("temp",t3d(k))
                CALL_DEBUG("qvi",qvi(k))
                CALL_DEBUG("qv3d",qv3d(k))
                CALL_DEBUG("qvqvsi",qvqvsi(k))
                CALL_DEBUG("evs",evs(k))
            endif
#endif
! ADD NUMBER CONCENTRATION DUE TO CUMULUS TENDENCY
! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM RAIN IS 10^7 M^-4
! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM SNOW IS 2 X 10^7 M^-4
! FOR DETRAINED CLOUD ICE, ASSUME MEAN VOLUME DIAM OF 80 MICRON

#ifndef WRF_MARS
            IF (QRCU1D(K).GE.1.E-10) THEN
            DUM=1.8e5*(QRCU1D(K)*DT/(PI*RHOW*RHO(K)**3))**0.25
            NR3D(K)=NR3D(K)+DUM
            END IF
            IF (QSCU1D(K).GE.1.E-10) THEN
            DUM=3.e5*(QSCU1D(K)*DT/(CONS1*RHO(K)**3))**(1./(DS+1.))
            NS3D(K)=NS3D(K)+DUM
            END IF
            IF (QICU1D(K).GE.1.E-10) THEN
            DUM=QICU1D(K)*DT/(CI*(80.E-6)**DI)
            NI3D(K)=NI3D(K)+DUM
            END IF

#endif
! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! hm modify 7/0/09 change limit to 1.e-8

#ifdef WRF_MARS
!negative fixes?
               if(qi3d(k) .lt. 0 .or. qcore3d(k) .lt. 0) then
                 if(-qi3d(k) .gt. qv3d(k)) then
                   !to much negative ice to fix
                   PRECQVRT_ERR = (qv3d(k) + qi3d(k))*rho(k)*dzq(k) !the mass of the error
                   qv3d(k) = 0.
                   qi3d(k) = 0.
                 else !negative but can be filled, or positive
                   qv3d(k) = qv3d(k) + qi3d(k) !remove water to fill hole
                   qi3d(k) = 0. !zero ice
                 endif
                 !fix core and number
                 if(-qcore3d(k) .gt. qdust3d(k)) then
                   PRECDUSTRT_ERR = PRECDUSTRT_ERR + (qdust3d(k) + qcore3d(k))*rho(k)*dzq(k) !mass
                   qcore3d(k) = 0.0
                   qdust3d(k) = 0.0
                   ndust3d(k) = 0.0
                   ni3d(k) = 0.0
                 else
                   qdust3d(k) = qdust3d(k) + qcore3d(k)
                   ndust3d(k) = ndust3d(k) + ni3d(k) !
                   ni3d(k) = 0.
                   qcore3d(k) = 0.
                 endif
               endif
!low value fixes
               IF (QI3D(K).LT.qsmall) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
!                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF
               
               !coreless ice crystals
               if (qcore3d(k).lt.qsmall .and. qi3d(k).gt.qsmall) then
                PRECRT_PHOTOLYSIS = PRECRT_PHOTOLYSIS + (qi3d(k))*rho(k)*dzq(k) !mass
                qi3d(k) = 0.0
                ni3d(k) = 0.0
               endif

               !too much water for the current model 0.1% of atmosphere and 10 time supersaturation.
               if (qv3d(k).gt.1e-4 .and. qvqvsi(k).gt.photolysis_scrit) then
                PRECRT_PHOTOLYSIS = PRECRT_PHOTOLYSIS + (qv3d(k))*rho(k)*dzq(k) !mass
                qv3d(k) = 0.
               endif

               qicetot3d(k) = qcore3d(k) + qi3d(k)
               IF (qicetot3d(k).LT.qsmall) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  QDUST3D(K) = QDUST3d(K)+ QCORE3D(k)
!                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  QCORE3D(K) = 0.
                  QI3D(K)=0.
               END IF

#else
!water sub-saturation low value fixes.
             IF (QVQVS(K).LT.0.9) THEN
               IF (QR3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QR3D(K)
                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  QR3D(K)=0.
               END IF
               IF (QC3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QC3D(K)
                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  QC3D(K)=0.
               END IF
             END IF

             IF (QVQVSI(K).LT.0.9) THEN
               IF (QI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF
               IF (QNI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QNI3D(K)
                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  QNI3D(K)=0.
               END IF
               IF (QG3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QG3D(K)
                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  QG3D(K)=0.
               END IF
             END IF
#endif
            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

#ifdef MARS_MORR_DEBUG
            if (iprint) then
                CALL_DEBUG("LOW_qvi",qvi(k))
                CALL_DEBUG("LOW_qv3d",qv3d(k))
                CALL_DEBUG("LOW_qvqvsi",qvqvsi(k))
                CALL_DEBUG("LOW_evs",evs(k))
            endif
#endif
! HEAT OF FUSION

#ifdef WRF_MARS
    XLF(K) = 0.0
#else
            XLF(K) = XXLS(K)-XXLV(K)
#endif

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(K).LT.QSMALL) THEN
         QC3D(K) = 0.
         NC3D(K) = 0.
         EFFC(K) = 0.
       END IF
       IF (QR3D(K).LT.QSMALL) THEN
         QR3D(K) = 0.
         NR3D(K) = 0.
         EFFR(K) = 0.
       END IF
#ifndef WRF_MARS
!we can't remove ice here, because of condensation onto cores later
       IF (QI3D(K).LT.QSMALL) THEN
         QI3D(K) = 0.
         NI3D(K) = 0.
         EFFI(K) = 0.
       END IF
#endif
       IF (QNI3D(K).LT.QSMALL) THEN
         QNI3D(K) = 0.
         NS3D(K) = 0.
         EFFS(K) = 0.
       END IF
       IF (QG3D(K).LT.QSMALL) THEN
         QG3D(K) = 0.
         NG3D(K) = 0.
         EFFG(K) = 0.
       END IF
#ifdef WRF_MARS    
       IF (qdust3d(K).LT.QSMALL) THEN
         qdust3d(K) = 0.
         ndust3d(K) = 0.
         EFFD(K) = 0.
       END IF
      qicetot3d(k) = qcore3d(k) + qi3d(k)
       if(qicetot3d(k).lt.qsmall) then 
         qcore3d(k) = 0.
         QI3D(K) = 0.
         NI3D(K) = 0.
         EFFI(K) = 0.
       ELSE IF(QCORE3D(K).lt.QSMALL) THEN
           QCORE3D(K) = 0.0
       ELSE IF(QI3D(K).lt.QSMALL) THEN
           QI3D(K)=0.0
       ENDIF

#endif


! INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO

      QRSTEN(K) = 0.
      QISTEN(K) = 0.
      QNISTEN(K) = 0.
      QCSTEN(K) = 0.
      QGSTEN(K) = 0.
      NISTEN(K) = 0.0
    
#ifdef WRF_MARS    
      QDUSTSTEN(K) = 0.0
      QCORESTEN(K) = 0.0
      NDUSTSTEN(K) = 0.0
#endif

!..................................................................
! MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT

#if defined WRF_MARS
!co2, Sutherland's formula (http://books.google.com/books?id=oRx6U4T8zcIC&pg=PA46)
			MU(K) = 1.572e-6*T3D(K)**1.5/(T3D(K)+240.)
!alternate from mp_sedim*
!			MU(K) = 1e-5
!alternate from mp_mars_common
!			MU(K) = sqrt(m*k*temp)/4*pi*r*r ~=2e-5
!alternate from OPUS-V, Lide et al. (1995)
!			MU(K) = (807 * (T/Tc)**0.618 - 357*exp(-0.449*T/Tc) + 340 *exp(-4.056*T/Tc) + 18) / 176*(Tc/m**3*Pc)**(1/6.)
#elif defined WRF_TITAN
!nitrogen
			MU(K) = 1.407e-6**T3D(K)**1.5/(T3D(K)+111.)
#else
!Earth
! fix 053011
            MU(K) = 1.496E-6*T3D(K)**1.5/(T3D(K)+120.)
#endif
! FALL SPEED WITH DENSITY CORRECTION (HEYMSFIELD AND BENSSEMER 2006)

            DUM = (RHOSU/RHO(K))**0.54

! fix 053011
!            AIN(K) = DUM*AI
! AA revision 4/1/11: Ikawa and Saito 1991 air-density correction 
#ifndef WRF_MARS
            AIN(K) = (RHOSU/RHO(K))**0.35*AI
#endif
            ARN(K) = DUM*AR
            ASN(K) = DUM*AS
!            ACN(K) = DUM*AC
! AA revision 4/1/11: temperature-dependent Stokes fall speed
            ACN(K) = G*RHOW/(18.*MU(K))
! HM ADD GRAUPEL 8/28/06
            AGN(K) = DUM*AG
#ifdef WRF_MARS
!uncorrected Stokes, the Cunningham factor is particle dependent and included later
			ADN(K) = G * DUST_CONS_RHO / (18.*MU(K)) !dust
			AIN(K) = G / (18.*MU(K)) !water ice
#endif
!hm 4/7/09 bug fix, initialize lami to prevent later division by zero
            LAMI(K)=0.

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, AND IF SUBSATURATED, THEN SKIP MICROPHYSICS
! FOR THIS LEVEL

            IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                 .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL &
#ifdef WRF_MARS
                 .and.qdust3d(k).lt.QSMALL &
                 .and.qcore3d(k).lt.QSMALL &
                 .and.qicetot3d(k).lt.QSMALL &
#endif
                 ) THEN
#ifdef WRF_MARS
                !nothing to nucleate onto, doesn't matter how much water we have
                 IF (mars_heterogeneous) goto 200 !
                !no saturation, skip the microphysics
                 IF (QVQVSI(K).LT.0.999) goto 200
#else
                 IF (T3D(K).LT.273.15.AND.QVQVSI(K).LT.0.999) goto 200
                 IF (T3D(K).GE.273.15.AND.QVQVS(K).LT.0.999) goto 200
#endif
            END IF

! THERMAL CONDUCTIVITY FOR AIR

! fix 053011
            KAP(K) = 1.414E3*MU(K)

! DIFFUSIVITY OF WATER VAPOR [m2/s]
#ifdef WRF_MARS
! This is the Washburn et al 2003 International Critical Tables (ICT) form
! via Hudson et al., 2008 in discussion above their Eqn 9 for water
! vapour molecular diffusion through a background CO2 gas - fit somewhat 
! dicey as low temp and low pressure data for the fit are not available
! (yes those look like Earth constants, yes this is what should be used
! for Mars, no those should not be changed to Mars-like values - they're
! just fit values directly from the table - leave them alone!)

! The generalized "best" way of doing this for planets (where bulk background
! gas is different (from Mars to Titan, to whatever)) is likely to take 
! the Holman [1997] generalized fit for binary mixture diffusion (again see 
! page 4 of Hudson et al. 2008 (JGR subsurface vapor diffusion paper). For
! now, this is the ICT fit:

            !DV(K) = 0.1387E-4*((T3D(K)/273.15)**2)*(101300./PRES(K))
            !The above fit is horrible in reality. blah blah Knudsen

            dv(k) = ((12./32.)/pi) * &
                        (1./(co2_molrad + h2o_molrad)**2) * &
                        (kboltz*T3D(K)/PRES(K)) * &
                        SQRT(1.+mass_co2/mass_h2o) * &
                        SQRT(8*kboltz*T3D(K)/(pi*mass_co2)) 
            
            !we'll defer the rest until later, becuase we have to correct based on the particle radius
            !in comparison to the mean free path. and there goes my two moment scheme. 




#else
            DV(K) = 8.794E-5*T3D(K)**1.81/PRES(K)
#endif

! SCHMIT NUMBER

! fix 053011
            SC(K) = MU(K)/(RHO(K)*DV(K))

! PSYCHOMETIC CORRECTIONS

! RATE OF CHANGE SAT. MIX. RATIO WITH TEMPERATURE

            DUM = (RV*T3D(K)**2)

            DQSDT = XXLV(K)*QVS(K)/DUM
            DQSIDT =  XXLS(K)*QVI(K)/DUM

            ABI(K) = 1.+DQSIDT*XXLS(K)/CPM(K)
            AB(K) = 1.+DQSDT*XXLV(K)/CPM(K)

! 
!.....................................................................
!.....................................................................
! CASE FOR TEMPERATURE ABOVE FREEZING

#ifndef WRF_MARS
!nothing in the 'above freezing' case applies to present day Mars 
!because the partial pressures are too low.
!search for ! TEMPERATURE < 273.15 to get to the end of this ifndef
            IF (T3D(K).GE.273.15) THEN

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         IF (iinum.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(K)=NDCNST*1.E6/RHO(K)
         END IF

! GET SIZE DISTRIBUTION PARAMETERS

! MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
       IF (QNI3D(K).LT.1.E-6) THEN
          QR3D(K)=QR3D(K)+QNI3D(K)
          NR3D(K)=NR3D(K)+NS3D(K)
          T3D(K)=T3D(K)-QNI3D(K)*XLF(K)/CPM(K)
          QNI3D(K) = 0.
          NS3D(K) = 0.
       END IF
       IF (QG3D(K).LT.1.E-6) THEN
          QR3D(K)=QR3D(K)+QG3D(K)
          NR3D(K)=NR3D(K)+NG3D(K)
          T3D(K)=T3D(K)-QG3D(K)*XLF(K)/CPM(K)
          QG3D(K) = 0.
          NG3D(K) = 0.
       END IF

       IF (QC3D(K).LT.QSMALL.AND.QNI3D(K).LT.1.E-8.AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.1.E-8 &
#ifdef WRF_MARS
			.and.qdust3d(k).lt.QSMALL &
#endif       
		) GOTO 300

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))
#ifdef WRF_MARS
	  ndust3d(K) = MAX(0.,ndust3d(K))
#endif
!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         DUM = PRES(K)/(287.15*T3D(K))
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      END IF
      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF
      END IF

#ifdef WRF_MARS
!DUST
    call constrain_mars_distribution(qdust3d(K), ndust3d(K), n0d(K), lamd(k))
#endif

!.....................................................................
! ZERO OUT PROCESS RATES

            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PSMLT(K) = 0.
            NSMLTS(K) = 0.
            NSMLTR(K) = 0.
            EVPMS(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            NSUBC(K) = 0.
            NSUBR(K) = 0.
            PRACG(K) = 0.
            NPRACG(K) = 0.
            PSMLT(K) = 0.
            PGMLT(K) = 0.
            EVPMG(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NGMLTG(K) = 0.
            NGMLTR(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         IF (QC3D(K).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1(K) = PRC(K)/CONS29
        NPRC(K) = PRC(K)/(QC3D(k)/NC3D(K))

! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

         END IF

!.......................................................................
! HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING
! FORMULA FROM IKAWA AND SAITO (1991)

         IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! hm fix, 2/12/13
! for above freezing conditions to get accelerated melting of snow,
! we need collection of rain by snow (following Lin et al. 1983)
!            PRACS(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
!                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
!                 N0RR(K)*N0S(K)/LAMS(K)**3*                    &
!                  (5./(LAMS(K)**3*LAMR(K))+                    &
!                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
!                  0.5/(LAMS(K)*LAMR(K)**3)))

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &				 
                  0.5/(LAMR(k)*LAMS(k)**3)))

! fix 053011, npracs no longer subtracted from snow
!            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
!                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
!                (1./(LAMR(K)**3*LAMS(K))+                      &
!                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
!                 1./(LAMR(K)*LAMS(K)**3))

         END IF

! ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
! ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
! ASSUME SHED DROPS ARE 1 MM IN SIZE

         IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! PRACG IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+				   &
				  0.5/(LAMR(k)*LAMG(k)**3)))

! ASSUME 1 MM DROPS ARE SHED, GET NUMBER SHED PER SEC

            DUM = PRACG(K)/5.2E-7

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! hm 7/15/13, remove limit so that the number of collected drops can smaller than 
! number of shed drops
!            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
            NPRACG(K)=NPRACG(K)-DUM

	    END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

           DUM=(QC3D(K)*QR3D(K))
           PRA(K) = 67.*(DUM)**1.15
           NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

         END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

         IF (QR3D(K).GE.1.E-8) THEN
! include breakup add 10/09/09
            dum1=300.e-6
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
!            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
            NRAGG(K) = -5.78*dum*NR3D(K)*QR3D(K)*RHO(K)
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP OF RAIN (RUTLEDGE AND HOBBS 1983)

      IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      ELSE
      EPSR = 0.
      END IF

! NO CONDENSATION ONTO RAIN, ONLY EVAP ALLOWED

           IF (QV3D(K).LT.QVS(K)) THEN
              PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K),0.)
           ELSE
              PRE(K) = 0.
           END IF

!.......................................................................
! MELTING OF SNOW

! SNOW MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

          IF (QNI3D(K).GE.1.E-8) THEN

! fix 053011
! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
!             DUM = -CPW/XLF(K)*T3D(K)*PRACS(K)
             DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACS(K)

! hm fix 1/20/15
!             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(273.15-T3D(K))/       &
!                    XLF(K)*RHO(K)*(F1S/(LAMS(K)*LAMS(K))+        &
!                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
!                    SC(K)**(1./3.)*CONS10/                   &
!                   (LAMS(K)**CONS35))+DUM
             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(273.15-T3D(K))/       &
                    XLF(K)*(F1S/(LAMS(K)*LAMS(K))+        &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
                   (LAMS(K)**CONS35))+DUM

! IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES

      IF (QVQVS(K).LT.1.) THEN
        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
! hm fix 8/4/08
        EVPMS(K) = (QV3D(K)-QVS(K))*EPSS/AB(K)    
        EVPMS(K) = MAX(EVPMS(K),PSMLT(K))
        PSMLT(K) = PSMLT(K)-EVPMS(K)
      END IF
      END IF

!.......................................................................
! MELTING OF GRAUPEL

! GRAUPEL MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

          IF (QG3D(K).GE.1.E-8) THEN

! fix 053011
! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
!             DUM = -CPW/XLF(K)*T3D(K)*PRACG(K)
             DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACG(K)

! hm fix 1/20/15
!             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(273.15-T3D(K))/ 		 &
!                    XLF(K)*RHO(K)*(F1S/(LAMG(K)*LAMG(K))+                &
!                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
!                    SC(K)**(1./3.)*CONS11/                   &
!                   (LAMG(K)**CONS36))+DUM
             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(273.15-T3D(K))/ 		 &
                    XLF(K)*(F1S/(LAMG(K)*LAMG(K))+                &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
                   (LAMG(K)**CONS36))+DUM

! IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES

      IF (QVQVS(K).LT.1.) THEN
        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))
! hm fix 8/4/08
        EVPMG(K) = (QV3D(K)-QVS(K))*EPSG/AB(K)
        EVPMG(K) = MAX(EVPMG(K),PGMLT(K))
        PGMLT(K) = PGMLT(K)-EVPMG(K)
      END IF
      END IF

! HM, V3.2
! RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
! TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
! ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

      PRACG(K) = 0.
      PRACS(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! FOR CLOUD ICE, ONLY PROCESSES OPERATING AT T > 273.15 IS
! MELTING, WHICH IS ALREADY CONSERVED DURING PROCESS
! CALCULATION

! CONSERVATION OF QC

      DUM = (PRC(K)+PRA(K))*DT

      IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN

        RATIO = QC3D(K)/DUM

        PRC(K) = PRC(K)*RATIO
        PRA(K) = PRA(K)*RATIO

        END IF

! CONSERVATION OF SNOW

        DUM = (-PSMLT(K)-EVPMS(K)+PRACS(K))*DT

        IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

! NO SOURCE TERMS FOR SNOW AT T > FREEZING
        RATIO = QNI3D(K)/DUM

        PSMLT(K) = PSMLT(K)*RATIO
        EVPMS(K) = EVPMS(K)*RATIO
        PRACS(K) = PRACS(K)*RATIO

        END IF

! CONSERVATION OF GRAUPEL

        DUM = (-PGMLT(K)-EVPMG(K)+PRACG(K))*DT

        IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

! NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
        RATIO = QG3D(K)/DUM

        PGMLT(K) = PGMLT(K)*RATIO
        EVPMG(K) = EVPMG(K)*RATIO
        PRACG(K) = PRACG(K)*RATIO

        END IF

! CONSERVATION OF QR
! HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

        DUM = (-PRACS(K)-PRACG(K)-PRE(K)-PRA(K)-PRC(K)+PSMLT(K)+PGMLT(K))*DT

        IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN

        RATIO = (QR3D(K)/DT+PRACS(K)+PRACG(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K))/ &
                        (-PRE(K))
        PRE(K) = PRE(K)*RATIO
        
        END IF

!....................................

      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-EVPMS(K)-EVPMG(K))

      T3DTEN(K) = T3DTEN(K)+(PRE(K)*XXLV(K)+(EVPMS(K)+EVPMG(K))*XXLS(K)+&
                    (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*XLF(K))/CPM(K)

      QC3DTEN(K) = QC3DTEN(K)+(-PRA(K)-PRC(K))
      QR3DTEN(K) = QR3DTEN(K)+(PRE(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K))
      QNI3DTEN(K) = QNI3DTEN(K)+(PSMLT(K)+EVPMS(K)-PRACS(K))
      QG3DTEN(K) = QG3DTEN(K)+(PGMLT(K)+EVPMG(K)-PRACG(K))
! fix 053011
!      NS3DTEN(K) = NS3DTEN(K)-NPRACS(K)
! HM, bug fix 5/12/08, npracg is subtracted from nr not ng
!      NG3DTEN(K) = NG3DTEN(K)
      NC3DTEN(K) = NC3DTEN(K)+ (-NPRA(K)-NPRC(K))
      NR3DTEN(K) = NR3DTEN(K)+ (NPRC1(K)+NRAGG(K)-NPRACG(K))

! HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC

	C2PREC(K) = PRA(K)+PRC(K)
      IF (PRE(K).LT.0.) THEN
         DUM = PRE(K)*DT/QR3D(K)
           DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
      END IF

        IF (EVPMS(K)+PSMLT(K).LT.0.) THEN
         DUM = (EVPMS(K)+PSMLT(K))*DT/QNI3D(K)
           DUM = MAX(-1.,DUM)
         NSMLTS(K) = DUM*NS3D(K)/DT
        END IF
        IF (PSMLT(K).LT.0.) THEN
          DUM = PSMLT(K)*DT/QNI3D(K)
          DUM = MAX(-1.0,DUM)
          NSMLTR(K) = DUM*NS3D(K)/DT
        END IF
        IF (EVPMG(K)+PGMLT(K).LT.0.) THEN
         DUM = (EVPMG(K)+PGMLT(K))*DT/QG3D(K)
           DUM = MAX(-1.,DUM)
         NGMLTG(K) = DUM*NG3D(K)/DT
        END IF
        IF (PGMLT(K).LT.0.) THEN
          DUM = PGMLT(K)*DT/QG3D(K)
          DUM = MAX(-1.0,DUM)
          NGMLTR(K) = DUM*NG3D(K)/DT
        END IF

         NS3DTEN(K) = NS3DTEN(K)+(NSMLTS(K))
         NG3DTEN(K) = NG3DTEN(K)+(NGMLTG(K))
         NR3DTEN(K) = NR3DTEN(K)+(NSUBR(K)-NSMLTR(K)-NGMLTR(K))

 300  CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)+DT*T3DTEN(K)
      DUMQV = QV3D(K)+DT*QV3DTEN(K)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(k),POLYSVP(DUMT,0))
      DUMQSS = EP_2*dum/(PRES(K)-dum)
      DUMQC = QC3D(K)+DT*QC3DTEN(K)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

      DUMS = DUMQV-DUMQSS
      PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
      IF (PCC(K)*DT+DUMQC.LT.0.) THEN
           PCC(K) = -DUMQC/DT
      END IF

      QV3DTEN(K) = QV3DTEN(K)-PCC(K)
      T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
      QC3DTEN(K) = QC3DTEN(K)+PCC(K)

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS
! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
! DROPLET CONCENTRATION IS SPECIFIED !!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC(K).LT.0.) THEN
!        DUM = PCC(K)*DT/QC3D(K)
!           DUM = MAX(-1.,DUM)
!        NSUBC(K) = DUM*NC3D(K)/DT
!     END IF

! UPDATE TENDENCIES

!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)

!.....................................................................
!.....................................................................
         ELSE  ! TEMPERATURE < 273.15
#endif
#ifdef WRF_MARS
!This is where all present day Mars stuff happens.
#endif

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

#ifndef WRF_MARS
         IF (iinum.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(K)=NDCNST*1.E6/RHO(K)
         END IF
#endif

! CALCULATE SIZE DISTRIBUTION PARAMETERS
! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))
#ifdef WRF_MARS
	  ndust3d(K) = MAX(0.,ndust3d(K))
#endif

!......................................................................
! CLOUD ICE
#ifdef WRF_MARS

    !allow for non-zero mu values
    if(.not.mars_heterogeneous) then !homo
        ice_rho(k) = rhoi
        qicetot3d(k) = qi3d(k)
        if(qi3d(k).ge.qsmall) then 
            call constrain_mars_waterice_homogeneous(qi3d(k), ni3d(k), n0i(k), lami(k))
        endif
    else
        qicetot3d(k) = qcore3d(k) + qi3d(k)
        if(qicetot3d(K).GE.QSMALL) then
            call mean_ice_density(qcore3d(k), qi3d(k), ice_rho(k))
#ifdef MARS_MORR_DEBUG
            CALL_DEBUG("call constrain",__LINE__)
#endif
            call constrain_mars_waterice_heterogeneous(qi3d(k), ni3d(k), qcore3d(k), ice_rho(k), n0i(k), lami(k),thisi, thisj, k)
        endif
    endif
#else
      IF (QI3D(K).GE.QSMALL) THEN
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)
         N0I(K) = NI3D(K)*LAMI(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(K).LT.LAMMINI) THEN

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      ELSE IF (LAMI(K).GT.LAMMAXI) THEN
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      END IF
      END IF
#endif


!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rain on Mars", k
#else
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
#endif
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have cloud water on Mars", k
#else
         DUM = PRES(K)/(287.15*T3D(K))
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

! TO CALCULATE DROPLET FREEZING

        CDIST1(K) = NC3D(K)/GAMMA(PGAM(K)+1.)
#endif
      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should have snow on Mars", k
#else
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      END IF
#endif
      END IF
!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should have Graupel on Mars", k
#else
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF
#endif
      END IF

#ifdef WRF_MARS
!DUST
    call constrain_mars_distribution(qdust3d(K), ndust3d(K), n0d(K), lamd(k))
#endif

!.....................................................................
! ZERO OUT PROCESS RATES

            MNUCCC(K) = 0.
            NNUCCC(K) = 0.
            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            NSAGG(K) = 0.
            PSACWS(K) = 0.
            NPSACWS(K) = 0.
            PSACWI(K) = 0.
            NPSACWI(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NMULTS(K) = 0.
            QMULTS(K) = 0.
            NMULTR(K) = 0.
            QMULTR(K) = 0.
            NMULTG(K) = 0.
            QMULTG(K) = 0.
            NMULTRG(K) = 0.
            QMULTRG(K) = 0.
            MNUCCR(K) = 0.
            NNUCCR(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PRCI(K) = 0.
            NPRCI(K) = 0.
            PRAI(K) = 0.
            NPRAI(K) = 0.
            NNUCCD(K) = 0.
            MNUCCD(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            PRD(K) = 0.
            PRDS(K) = 0.
            EPRD(K) = 0.
            EPRDS(K) = 0.
            NSUBC(K) = 0.
            NSUBI(K) = 0.
            NSUBS(K) = 0.
            NSUBR(K) = 0.
            PIACR(K) = 0.
            NIACR(K) = 0.
            PRACI(K) = 0.
            PIACRS(K) = 0.
            NIACRS(K) = 0.
            PRACIS(K) = 0.
! HM: ADD GRAUPEL PROCESSES
            PRACG(K) = 0.
            PSACR(K) = 0.
	    PSACWG(K) = 0.
	    PGSACW(K) = 0.
            PGRACS(K) = 0.
	    PRDG(K) = 0.
	    EPRDG(K) = 0.
	    NPRACG(K) = 0.
	    NPSACWG(K) = 0.
	    NSCNG(K) = 0.
 	    NGRACS(K) = 0.
	    NSUBG(K) = 0.

#ifdef WRF_MARS
    MNUCCORE(K) = 0.
    NSUBDUST(K) = 0.
    QSUBDUST(K) = 0.
    QSUBCORE(K) = 0.
#endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES
! ACCRETION/AUTOCONVERSION/FREEZING/MELTING/COAG.
!.......................................................................
! FREEZING OF CLOUD DROPLETS
! ONLY ALLOWED BELOW -4 C
        IF (QC3D(K).GE.QSMALL .AND. T3D(K).LT.269.15) THEN

#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : shouldn't have freezing cloud on Mars"
#else
! NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
! FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3

! MEYERS CURVE

           NACNT = EXP(-2.80+0.262*(273.15-T3D(K)))*1000.

! COOPER CURVE
!        NACNT =  5.*EXP(0.304*(273.15-T3D(K)))

! FLECTHER
!     NACNT = 0.01*EXP(0.6*(273.15-T3D(K)))

! CONTACT FREEZING

! MEAN FREE PATH

            DUM = 7.37*T3D(K)/(288.*10.*PRES(K))/100.

! EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
! BASED ON BROWNIAN DIFFUSION

            DAP(K) = CONS37*T3D(K)*(1.+DUM/RIN)/MU(K)
 
           MNUCCC(K) = CONS38*DAP(K)*NACNT*EXP(LOG(CDIST1(K))+   &
                   LOG(GAMMA(PGAM(K)+5.))-4.*LOG(LAMC(K)))
           NNUCCC(K) = 2.*PI*DAP(K)*NACNT*CDIST1(K)*           &
                    GAMMA(PGAM(K)+2.)/                         &
                    LAMC(K)

! IMMERSION FREEZING (BIGG 1953)

!           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
!                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
!                   EXP(AIMM*(273.15-T3D(K)))

!           NNUCCC(K) = NNUCCC(K)+                                  &
!            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
!                *EXP(AIMM*(273.15-T3D(K)))

! hm 7/15/13 fix for consistency w/ original formula
           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
                   (EXP(AIMM*(273.15-T3D(K)))-1.)

           NNUCCC(K) = NNUCCC(K)+                                  &
            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
                *(EXP(AIMM*(273.15-T3D(K)))-1.)

! PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
! MIXING RATIO, SINCE STRICT CONSERVATION NOT CHECKED FOR NUMBER CONC

           NNUCCC(K) = MIN(NNUCCC(K),NC3D(K)/DT)

#endif
        END IF

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         IF (QC3D(K).GE.1.E-6) THEN

#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : shouldn't have cloud on Mars"
#else
! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1(K) = PRC(K)/CONS29
        NPRC(K) = PRC(K)/(QC3D(K)/NC3D(K))

! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

#endif
         END IF

!.......................................................................
! SELF-COLLECTION OF DROPLET NOT INCLUDED IN KK2000 SCHEME

! SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
! THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW

         IF (QNI3D(K).GE.1.E-8) THEN
#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : should we have snow on Mars"
#else
             NSAGG(K) = CONS15*ASN(K)*RHO(K)**            &
            ((2.+BS)/3.)*QNI3D(K)**((2.+BS)/3.)*                  &
            (NS3D(K)*RHO(K))**((4.-BS)/3.)/                       &
            (RHO(K))
#endif
         END IF

!.......................................................................
! ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
! HERE USE CONTINUOUS COLLECTION EQUATION WITH
! SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

! SNOW

         IF (QNI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : should we have snow on Mars?"
#else
           PSACWS(K) = CONS13*ASN(K)*QC3D(K)*RHO(K)*               &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.)
           NPSACWS(K) = CONS13*ASN(K)*NC3D(K)*RHO(K)*              &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.)

#endif
         END IF

!............................................................................
! COLLECTION OF CLOUD WATER BY GRAUPEL

         IF (QG3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : shouldn't have graupel on Mars"
#else
           PSACWG(K) = CONS14*AGN(K)*QC3D(K)*RHO(K)*               &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.)
           NPSACWG(K) = CONS14*AGN(K)*NC3D(K)*RHO(K)*              &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.)
#endif
	    END IF

!.......................................................................
! HM, ADD 12/13/06
! CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
! BEFORE RIMING CAN OCCUR
! ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
! TO HALLET-MOSSOP SPLINTERING

         IF (QI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : shouldn't have rimed snow on Mars"
#else
! PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
! FROM THOMPSON ET AL. 2004, MWR

            IF (1./LAMI(K).GE.100.E-6) THEN

           PSACWI(K) = CONS16*AIN(K)*QC3D(K)*RHO(K)*               &
                  N0I(K)/                        &
                  LAMI(K)**(BI+3.)
           NPSACWI(K) = CONS16*AIN(K)*NC3D(K)*RHO(K)*              &
                  N0I(K)/                        &
                  LAMI(K)**(BI+3.)
           END IF
#endif
         END IF

!.......................................................................
! ACCRETION OF RAIN WATER BY SNOW
! FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998

         IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN
#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : shouldn't have rain on Mars"
#else

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &				 
                  0.5/(LAMR(k)*LAMS(k)**3)))

            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
                (1./(LAMR(K)**3*LAMS(K))+                      &
                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
                 1./(LAMR(K)*LAMS(K)**3))

! MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACS(K) = MIN(PRACS(K),QR3D(K)/DT)

! COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
! ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG

! HM MODIFY FOR WRFV3.1
!            IF (IHAIL.EQ.0) THEN
            IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
            PSACR(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
                 N0RR(K)*N0S(K)/LAMS(K)**3*                               &
                  (5./(LAMS(K)**3*LAMR(K))+                    &
                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
                  0.5/(LAMS(K)*LAMR(K)**3)))            
            END IF
!            END IF

#endif
         END IF

!.......................................................................

! COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990, 
! USED BY REISNER ET AL 1998
         IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : shouldn't have rain on Mars"
#else
            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+				   &
				  0.5/(LAMR(k)*LAMG(k)**3)))

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACG(K) = MIN(PRACG(K),QR3D(K)/DT)

#endif
	    END IF

!.......................................................................
! RIME-SPLINTERING - SNOW
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW AND DROPLET MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS
! THESE THRESHOLDS CORRESPOND WITH GRAUPEL THRESHOLDS IN RH 1984

!v1.4
         IF (QNI3D(K).GE.0.1E-3) THEN
         IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : shouldn't have rain or cloud on Mars"
#else
         IF (PSACWS(K).GT.0..OR.PRACS(K).GT.0.) THEN
            IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

               IF (T3D(K).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                  FMULT = (270.16-T3D(K))/2.
               ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                  FMULT = (T3D(K)-265.16)/3.
               ELSE IF (T3D(K).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW

               IF (PSACWS(K).GT.0.) THEN
                  NMULTS(K) = 35.E4*PSACWS(K)*FMULT*1000.
                  QMULTS(K) = NMULTS(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                  QMULTS(K) = MIN(QMULTS(K),PSACWS(K))
                  PSACWS(K) = PSACWS(K)-QMULTS(K)

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACS(K).GT.0.) THEN
                   NMULTR(K) = 35.E4*PRACS(K)*FMULT*1000.
                   QMULTR(K) = NMULTR(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                   QMULTR(K) = MIN(QMULTR(K),PRACS(K))

                   PRACS(K) = PRACS(K)-QMULTR(K)

               END IF

            END IF
         END IF
#endif
         END IF
         END IF

!.......................................................................
! RIME-SPLINTERING - GRAUPEL 
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS

!         IF (IHAIL.EQ.0) THEN
! v1.4
         IF (QG3D(K).GE.0.1E-3) THEN
#ifdef WRF_MARS
        write(*,*) "micro: ",__LINE__," : shouldn't have rime graupel on Mars"
#else
         IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
         IF (PSACWG(K).GT.0..OR.PRACG(K).GT.0.) THEN
            IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

               IF (T3D(K).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                  FMULT = (270.16-T3D(K))/2.
               ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                  FMULT = (T3D(K)-265.16)/3.
               ELSE IF (T3D(K).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL

               IF (PSACWG(K).GT.0.) THEN
                  NMULTG(K) = 35.E4*PSACWG(K)*FMULT*1000.
                  QMULTG(K) = NMULTG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                  QMULTG(K) = MIN(QMULTG(K),PSACWG(K))
                  PSACWG(K) = PSACWG(K)-QMULTG(K)

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACG(K).GT.0.) THEN
                   NMULTRG(K) = 35.E4*PRACG(K)*FMULT*1000.
                   QMULTRG(K) = NMULTRG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                   QMULTRG(K) = MIN(QMULTRG(K),PRACG(K))
                   PRACG(K) = PRACG(K)-QMULTRG(K)

               END IF
               END IF
               END IF
            END IF
#endif
            END IF
!         END IF

!........................................................................
! CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL/HAIL

!           IF (IHAIL.EQ.0) THEN
	   IF (PSACWS(K).GT.0.) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rimed cloud water on Mars", k
#else
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(K).GE.0.1E-3.AND.QC3D(K).GE.0.5E-3) THEN

! PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
	     PGSACW(K) = MIN(PSACWS(K),CONS17*DT*N0S(K)*QC3D(K)*QC3D(K)* &
                          ASN(K)*ASN(K)/ &
                           (RHO(K)*LAMS(K)**(2.*BS+2.))) 

! MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
	     DUM = MAX(RHOSN/(RHOG-RHOSN)*PGSACW(K),0.) 

! NUMBER CONCENTRAITON OF EMBRYO GRAUPEL FROM RIMING OF SNOW
	     NSCNG(K) = DUM/MG0*RHO(K)
! LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
             NSCNG(K) = MIN(NSCNG(K),NS3D(K)/DT)

! PORTION OF RIMING LEFT FOR SNOW
             PSACWS(K) = PSACWS(K) - PGSACW(K)
             END IF
#endif
	   END IF

! CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL

	   IF (PRACS(K).GT.0.) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rimed rainwater on Mars", k
#else
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
! PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
	      DUM = CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3 &    
                   /(CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3+ &  
                   CONS19*(4./LAMR(K))**3*(4./LAMR(K))**3)
              DUM=MIN(DUM,1.)
              DUM=MAX(DUM,0.)
	      PGRACS(K) = (1.-DUM)*PRACS(K)
            NGRACS(K) = (1.-DUM)*NPRACS(K)
! LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
            NGRACS(K) = MIN(NGRACS(K),NR3D(K)/DT)
            NGRACS(K) = MIN(NGRACS(K),NS3D(K)/DT)

! AMOUNT LEFT FOR SNOW PRODUCTION
            PRACS(K) = PRACS(K) - PGRACS(K)
            NPRACS(K) = NPRACS(K) - NGRACS(K)
! CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
            PSACR(K)=PSACR(K)*(1.-DUM)
            END IF
#endif
	   END IF
!           END IF

!.......................................................................
! FREEZING OF RAIN DROPS
! FREEZING ALLOWED BELOW -4 C

         IF (T3D(K).LT.269.15.AND.QR3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rain to freeze on Mars", k
#else

! IMMERSION FREEZING (BIGG 1953)
!            MNUCCR(K) = CONS20*NR3D(K)*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3 &
!                 /LAMR(K)**3

!            NNUCCR(K) = PI*NR3D(K)*BIMM*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3

! hm fix 7/15/13 for consistency w/ original formula
            MNUCCR(K) = CONS20*NR3D(K)*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3 &
                 /LAMR(K)**3

            NNUCCR(K) = PI*NR3D(K)*BIMM*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3

! PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
            NNUCCR(K) = MIN(NNUCCR(K),NR3D(K)/DT)

#endif
         END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have cloud water on Mars", k
#else
! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

           DUM=(QC3D(K)*QR3D(K))
           PRA(K) = 67.*(DUM)**1.15
           NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

#endif
         END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

         IF (QR3D(K).GE.1.E-8) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rain on Mars", k
#else

! include breakup add 10/09/09
            dum1=300.e-6
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
!            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
            NRAGG(K) = -5.78*dum*NR3D(K)*QR3D(K)*RHO(K)
#endif
         END IF

!.......................................................................
! AUTOCONVERSION OF CLOUD ICE TO SNOW
! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION
#ifdef MARS_SNOW
         IF (QI3D(K).GE.1.E-8 .AND.QVQVSI(K).GE.1.) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should we have snow on Mars?", k
#else

!           COFFI = 2./LAMI(K)
!           IF (COFFI.GE.DCS) THEN
              NPRCI(K) = CONS21*(QV3D(K)-QVI(K))*RHO(K)                         &
                *N0I(K)*EXP(-LAMI(K)*DCS)*DV(K)/ABI(K)
              PRCI(K) = CONS22*NPRCI(K)
              NPRCI(K) = MIN(NPRCI(K),NI3D(K)/DT)

!           END IF
#endif
         END IF

!.......................................................................
! ACCRETION OF CLOUD ICE BY SNOW
! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
! AND DS >> DI FOR CONTINUOUS COLLECTION

         IF (QNI3D(K).GE.1.E-8 .AND. QI3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should we have snow on Mars", k
#else
            PRAI(K) = CONS23*ASN(K)*QI3D(K)*RHO(K)*N0S(K)/     &
                     LAMS(K)**(BS+3.)
            NPRAI(K) = CONS23*ASN(K)*NI3D(K)*                                       &
                  RHO(K)*N0S(K)/                                 &
                  LAMS(K)**(BS+3.)
            NPRAI(K)=MIN(NPRAI(K),NI3D(K)/DT)
#endif
         END IF

!.......................................................................
! HM, ADD 12/13/06, COLLISION OF RAIN AND ICE TO PRODUCE SNOW OR GRAUPEL
! FOLLOWS REISNER ET AL. 1998
! ASSUMED FALLSPEED AND SIZE OF ICE CRYSTAL << THAN FOR RAIN

         IF (QR3D(K).GE.1.E-8.AND.QI3D(K).GE.1.E-8.AND.T3D(K).LE.273.15) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rain on Mars", k
#else

! ALLOW GRAUPEL FORMATION FROM RAIN-ICE COLLISIONS ONLY IF RAIN MIXING RATIO > 0.1 G/KG,
! OTHERWISE ADD TO SNOW

            IF (QR3D(K).GE.0.1E-3) THEN
            NIACR(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACR(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACI(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACR(K)=MIN(NIACR(K),NR3D(K)/DT)
            NIACR(K)=MIN(NIACR(K),NI3D(K)/DT)
            ELSE 
            NIACRS(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACRS(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACIS(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACRS(K)=MIN(NIACRS(K),NR3D(K)/DT)
            NIACRS(K)=MIN(NIACRS(K),NI3D(K)/DT)
            END IF
#endif
         END IF
#endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL
#ifdef WRF_MARS
!homogeneous scrit is 1.0, hetero is more complex
      scrit_nuc  = 1.0
      scrit_cond = 1.0
#endif      
         IF (INUC.EQ.0) THEN

! add threshold according to Greg Thomspon

         if ((QVQVS(K).GE.0.999.and.T3D(K).le.265.15).or. &
              QVQVSI(K).ge.1.08) then

! hm, modify dec. 5, 2006, replace with cooper curve
      kc2 = 0.005*exp(0.304*(273.15-T3D(K)))*1000. ! convert from L-1 to m-3
! limit to 500 L-1
      kc2 = min(kc2,500.e3)
      kc2=MAX(kc2/rho(k),0.)  ! convert to kg-1

          IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
             NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD(K) = NNUCCD(K)*MI0
          END IF

          END IF

          ELSE IF (INUC.EQ.1) THEN

          IF (T3D(K).LT.273.15.AND.QVQVSI(K).GT.1.) THEN

             KC2 = 0.16*1000./RHO(K)  ! CONVERT FROM L-1 TO KG-1
          IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
             NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD(K) = NNUCCD(K)*MI0
          END IF
          END IF

#ifdef WRF_MARS
          ELSE IF (INUC .EQ. 2) THEN
! Mars nuc on dust here - direct heterogeneous nucleation on dust is not
! treated in the Earth scheme since Waicek and Peter [2009] found that almost
! all dust has been processed through liquid water and hence not viable CCN

! implement (likely) equations 1-7 from Cziczo et al [2013] and modified form of Eqn 7 from
! Trainer et al 2009 (modified as in m seems lower than most observational constraints
! tabulated in Cziczo et al - should probably just generate a new fit for m(T)   

! should result in:
           NNUCCD(K) = 0.0! nucleation rate from Cziczo eqn 1 * 4*pi*rnuc (radius used in MI0 calc)
           MNUCCD(K) = 0.0!NNUCCD(K)*MI0

! need to add removal term for dust such that number density decrease in dust = number density inc. in ice
! mass loss for dust = nucleation rate for ice * volume of dust particle, even though we don't add this mass
! to the formed ice particle.

! actually, size of nucleated particle used in MI0 should be the dust size used for nuc and the size 
! should be that part of the dust spectrum that has enough dust to nucleate the required number of ice particles
! might need to iterate to solve for this...
      
     ELSE IF (INUC .eq. 3) THEN
        !Reproduction of the mars simple water scheme where ice is instantaneously created 
!ice mass 
        qicetot3d(k) = qcore3d(k) + qi3d(k)
        call homogeneous_ice_physics(qv3d(k), qi3d(k), ni3d(k), qvi(k), dt,  MNUCCD(K), NNUCCD(K), PRD(K))
     ELSE IF (INUC .eq. 4) THEN
        qicetot3d(k) = qcore3d(k) + qi3d(k)
        call homogeneous_ice_physics(qv3d(k), qi3d(k), ni3d(k), qvi(k), dt,  MNUCCD(K), NNUCCD(K), PRD(K))
    ELSE IF (INUC .eq. 5) THEN
        qicetot3d(k) = qcore3d(k) + qi3d(k)
        !Reproduction of the mars simple water scheme where ice is instantaneously created 
        !but with the requirement that dust is present
        if(qdust3d(k).gt.qsmall) then 
            call homogeneous_ice_physics(qv3d(k), qi3d(k), ni3d(k), qvi(k), dt,  MNUCCD(K), NNUCCD(K), PRD(K))
        endif
    else if (mars_heterogeneous) then !inuc .eq. 6 .or. inuc .eq. 7) then
        !heterogeneous ice nucleation
        qicetot3d(k) = qcore3d(k) + qi3d(k)
        scrit = 1.0!2.75 * exp(-(T3D(K)-150.)/20.) + 1.25
        scrit_nuc = scrit
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("scrit", scrit)
#endif
        if(qdust3d(k).gt.qsmall .and. qvqvsi(k) .gt. scrit) then
            !needs to know the vapor pressure and temperature
#ifdef MARS_MORR_DEBUG
            CALL_DEBUG("qv3d", qv3d(k))
            CALL_DEBUG("t3d", t3d(k))
            CALL_DEBUG("qdust3d", qdust3d(k))
            CALL_DEBUG("qi3d", qi3d(k))
            CALL_DEBUG("ni3d", ni3d(k))
            CALL_DEBUG("qcore3d", qcore3d(k))
            CALL_DEBUG("k",k)
            CALL_DEBUG("qsmall", qsmall)
            CALL_DEBUG("scrit",scrit)
            CALL_DEBUG("qvqvsi",qvqvsi(k))
            CALL_DEBUG("qdust3d",qdust3d(k))
            CALL_DEBUG("lamd", lamd(k))
            if(lamd(k)>0) CALL_DEBUG("reff", 1e6*(dust_mu+3)/(2*lamd(k)))
            CALL_DEBUG("pres",pres(k))
            CALL_DEBUG("ep_2",ep_2)
            CALL_DEBUG("1./ep_2", 1./ep_2)
#endif
            !call heterogeneous_ice_nucleation(pres(k)*qv3d(k)/ep_2, t3d(k), lamd(k), n0d(k), qdust3d(k), ndust3d(k), &
            !    core_number_rate, core_mass_rate, ice_mass_rate)
             !write(*,*) "nucleation:", k, ndust3d(k)
             call heterogeneous_ice_nucleation_binned(pres(k)*qv3d(k)/ep_2, t3d(k), qvqvsi(k),lamd(k), n0d(k), qdust3d(k), ndust3d(k), &
              number(:), diam(:), mass(:),dt,&
                core_number_rate, core_mass_rate, ice_mass_rate)

  
            !now we have the possible rates, we need to assign them to variables for later use
            MNUCCORE(K) = core_mass_rate
            NNUCCD(K) = core_number_rate !core number is stored in the ice.
            MNUCCD(K) = ice_mass_rate
#ifdef MARS_MORR_DEBUG
            CALL_DEBUG("hetero nuc",0.)
            CALL_DEBUG("k",k)
            CALL_DEBUG("MNUCCD",MNUCCD(k))
            CALL_DEBUG("NNUCCD",NNUCCD(k))
            CALL_DEBUG("MNUCCORE",MNUCCORE(k))
            CALL_DEBUG("qdust3d",qdust3d(k))
            CALL_DEBUG("ndust3d",ndust3d(k))
#endif
        endif

#endif
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 101      CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR

! NO VENTILATION FOR CLOUD ICE
#ifdef WRF_MARS
! choice of ventilation comes from fact that ice particles do not fall
! very fast, which comes back to definition of ice vs. snow being one of
! fall speeds
#endif
#ifndef WRF_MARS
        IF (QI3D(K).GE.QSMALL) THEN

         EPSI = 2.*PI*N0I(K)*RHO(K)*DV(K)/(LAMI(K)*LAMI(K))

      ELSE
         EPSI = 0.
      END IF
#endif
      IF (QNI3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should we have snow on Mars?", k
#else

        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
#endif
      ELSE
      EPSS = 0.
      END IF

      IF (QG3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have graupel on Mars", k
#else

        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))


#endif
      ELSE
      EPSG = 0.
      END IF

      IF (QR3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rain on Mars", k
#else

        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
#endif
      ELSE
      EPSR = 0.
      END IF
#ifdef WRF_MARS
!condensation onto ice

    if (mars_growth) then !inuc.ge.4) then 
    prd(k) = 0.0
    if(qicetot3d(k) .ge. qsmall) then !if there is CORE mass, we can condense onto it!
        !condensation changes mass but not number
        !2.pi.rho_i*n0 / (Fd+Fh) * gamma(mu+1)/lambda^(mu+1) * ((s-1)/lambda - M_w.sigma/(Rd.T.rho_i))
        !The denominators Fd and Fh are calculated first
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("condensation",0)
#endif

      call condense_binned(T3D(k),rho(k),qv3d(k),QVI(k), number(:),diam(:),mass(:),dv(k),XXLS(k),eis(k),lami(k), ni3d(k),dt, prd(k))
      !write(*,*) "condense: ", k,prd(k), qi3d(k), qv3d(k), ni3d(k), ndust3d(k), qdust3d(k)
        !value of qvqvsi that would shut down prd.
        !scrit_cond = (4./RV) * (sigvl / (T3D(K)*RHOI))*LAMI(K)/(ice_mu+1) + 1
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("scrit_nuc", scrit_nuc)  
        CALL_DEBUG("scrit_cond", scrit_cond)
        CALL_DEBUG("QVQVSI",QVQVSI(k))
        CALL_DEBUG("QV3D",QV3D(k))
        CALL_DEBUG("QVI",QVI(k))
        CALL_DEBUG("eis",eis(k))
        CALL_DEBUG("ep_2",ep_2)
        CALL_DEBUG("pres",pres(k))
        CALL_DEBUG("rv",rv)
        CALL_DEBUG("rhoi",rhoi)
        CALL_DEBUG("t3d",t3d(k))
        CALL_DEBUG("dv",dv(k))
        CALL_DEBUG("eis",eis(k))
        CALL_DEBUG("sigvl",sigvl)
        CALL_DEBUG("sigvl", sigvl)
        CALL_DEBUG("nuc_ka",nuc_ka)
        CALL_DEBUG("nuc_Fh",nuc_Fh)
        CALL_DEBUG("XXLS",  XXLS(k))
        CALL_DEBUG("rhoi",  rhoi)
        CALL_DEBUG("RV",    RV)
        CALL_DEBUG("EPSI",  EPSI)
        CALL_DEBUG("ni3d",  ni3d(k))
        CALL_DEBUG("nuc_Fd",nuc_Fd)
        CALL_DEBUG("nuc_Fh",nuc_Fh)
        CALL_DEBUG("ice_mu",ice_mu)
        CALL_DEBUG("qvqvsi",qvqvsi(k))
        CALL_DEBUG("lami",  lami(k))
        CALL_DEBUG("epsi", epsi)
        CALL_DEBUG("term1", (ice_mu+1)*(QVQVSI(K) - 1)/LAMI(K))
        CALL_DEBUG("term2",(4./RV) * (sigvl / (T3D(K)*RHOI)))
        CALL_DEBUG("prd",prd(k))
        CALL_DEBUG("k",k)
        CALL_DEBUG("qicetot3d",qicetot3d(k))
        CALL_DEBUG("qi3d",qi3d(k))
        CALL_DEBUG("qcore3d",qcore3d(k))
#endif
    endif
    endif
#else
! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
! DUM IS FRACTION OF D*N(D) < DCS

! LOGIC BELOW FOLLOWS THAT OF HARRINGTON ET AL. 1995 (JAS)
              IF (QI3D(K).GE.QSMALL) THEN              
              DUM=(1.-EXP(-LAMI(K)*DCS)*(1.+LAMI(K)*DCS))
              PRD(K) = EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM
              ELSE
              DUM=0.
              END IF
! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
              IF (QNI3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should we have snow on Mars?", k
#else
              PRDS(K) = EPSS*(QV3D(K)-QVI(K))/ABI(K)+ &
                EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
! OTHERWISE ADD TO CLOUD ICE
#endif
              ELSE
              PRD(K) = PRD(K)+EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
              END IF
#ifndef WRF_MARS
! VAPOR DPEOSITION ON GRAUPEL
              PRDG(K) = EPSG*(QV3D(K)-QVI(K))/ABI(K)

! NO CONDENSATION ONTO RAIN, ONLY EVAP

           IF (QV3D(K).LT.QVS(K)) THEN
              PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K),0.)
           ELSE
              PRE(K) = 0.
           END IF

#endif
#endif
! MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
! FORMULA FROM REISNER 2 SCHEME
#ifndef WRF_MARS

           DUM = (QV3D(K)-QVI(K))/DT

           FUDGEF = 0.9999
           SUM_DEP = PRD(K)+PRDS(K)+MNUCCD(K)+PRDG(K)

           IF( (DUM.GT.0. .AND. SUM_DEP.GT.DUM*FUDGEF) .OR.                      &
               (DUM.LT.0. .AND. SUM_DEP.LT.DUM*FUDGEF) ) THEN
               MNUCCD(K) = FUDGEF*MNUCCD(K)*DUM/SUM_DEP
               PRD(K) = FUDGEF*PRD(K)*DUM/SUM_DEP
               PRDS(K) = FUDGEF*PRDS(K)*DUM/SUM_DEP
               PRDG(K) = FUDGEF*PRDG(K)*DUM/SUM_DEP
           ENDIF
#else
#ifdef MARS_MORR_DEBUG
      if(iprint) THEN
        CALL_DEBUG("A_MNUCCORE",MNUCCORE(K))
        CALL_DEBUG("A_MNUCCD",MNUCCD(K))
        CALL_DEBUG("A_NNUCCD",NNUCCD(K))
        !CALL_DEBUG("A_FUDGEF",FUDGEF)
        CALL_DEBUG("A_PRD",PRD(K))
        !CALL_DEBUG("A_PRDS",PRDS(K))
        !CALL_DEBUG("A_PRDG",PRDG(K))
        CALL_DEBUG("A_DUST",qdust3d(K)/dt)
      endif
#endif
    !super-saturation
    !for nucleation events we are rate limited by two factors.
    !The first rate limit is the ice consumption (qv-qvi*s_nuc)/dt /(PRD+MNUCCD)
    !The second rate limit is the dust consumption (qdust/dt)/(MNUCCORE)
    !The lower rate limit is applied to MNUUCD, and the remainder of the water is available to condensation only.
    !homogeneous nucleation doesn't consume dust and doesn't trigger nucleation, so isn't affected by any of this.
           FUDGEF = 0.9999
           SUM_DEP = PRD(K)+PRDS(K)+MNUCCD(K)+PRDG(K)
           DUM_nuc = max((QV3D(K)-QVI(K)*scrit_nuc)/DT,0.)

           dum_dust = max(qdust3d(k)/DT,0.)
           sum_dep_dust = MNUCCORE(k)
           ratio_water= 1.
           ratio_dust = 1.
#ifdef MARS_MORR_DEBUG
      if(iprint) THEN
           CALL_DEBUG("DUM_nuc", DUM_nuc)
      endif
#endif
           !if we can reduce supersaturation to the nucleation threshold, do it, and scale the nucleation rates appropriately.
           if( (dum_nuc.gt.0 .and. sum_dep.gt.dum_nuc*FUDGEF) .or. &
               (dum_dust.gt.0 .and. sum_dep_dust.gt.dum_dust*FUDGEF) ) then

              if(sum_dep .gt. 0) ratio_water = dum_nuc*FUDGEF/sum_dep
              if(sum_dep_dust .gt. 0) ratio_dust = dum_dust*FUDGEF/sum_dep_dust
              if(ratio_dust.le.ratio_water) then
                !strip the dust
                MNUCCORE(K) = qdust3d(k)/dt !ratio_dust * MNUCCORE = (qdust3d/dt)/mnuccore * mnuccore
                NNUCCD(K) = ndust3d(k)/dt !have to strip the number too.
                MNUCCD(K) = MNUCCD(K)*ratio_dust
              else
                MNUCCORE(K) = MNUCCORE(K) * ratio_water
                NNUCCD(K) = NNUCCD(k) * ratio_water
                MNUCCD(K) = MNUCCD(K)*ratio_water
              endif

             endif
           !now we are no longer dum_nuc.gt.0 at the end of the timestep, so we could be above or below condensation level
           SUM_DEP = PRD(K)+PRDS(K)+PRDG(K)

           !The constraint > 0 for nucleation limits above is somewhat redundant because being above the nucleation level
           !implies being above the condensation level which implies a position Dum_nuc at all times. At no time is there
           !a 'negative nucleation' event that would remove all ice from a dust core. 
           !but constraining dum_cond below would be (was...) bad because if qv<qvi(cond) then evaporation occurs as a 
           !'negative condensation' to remove any or all ice from the dust.
           !if we do have negative condensation, then nucleation (MNUCCD) must be zero.
           dum_cond = (QV3D(K)-QVI(K)*scrit_cond)/DT - MNUCCD(K)
#ifdef MARS_MORR_DEBUG
      if(iprint) THEN
           CALL_DEBUG("DUM_cond", DUM_cond)
           CALL_DEBUG("DUM_cond2", (QV3D(K)-QVI(K)*scrit_cond)/DT - MNUCCD(K))

      endif
#endif
           if( (dum_cond.gt.0 .and. sum_dep.gt.dum_cond*FUDGEF) .or.    &
               (dum_cond.lt.0 .and. sum_dep.lt.dum_cond*FUDGEF) ) then
                PRD(K) = FUDGEF*PRD(K)*dum_cond/sum_dep
                PRDS(K) = FUDGEF*PRDS(K)*dum_cond/sum_dep !should be zero anyway
                PRDG(K) = FUDGEF*PRDG(K)*dum_cond/sum_dep !should be zero anyway
           endif
           !with an assumption here that dum_cond and sum_dep never have a different sign because that's non-physical.
#endif
#ifdef MARS_MORR_DEBUG
      if(iprint) THEN
        CALL_DEBUG("MNUCCORE",MNUCCORE(K))
        CALL_DEBUG("MNUCCD",MNUCCD(K))
        CALL_DEBUG("NNUCCD",NNUCCD(K))
        !CALL_DEBUG("FUDGEF",FUDGEF)
        CALL_DEBUG("DUM",DUM)
        CALL_DEBUG("SUM_DEP",SUM_DEP)
        CALL_DEBUG("PRD",PRD(K))
        !CALL_DEBUG("PRDS",PRDS(K))
        !CALL_DEBUG("PRDG",PRDG(K))
        CALL_DEBUG("DUST",qdust3d(K)/dt)
      endif
#endif

! IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES

           IF (PRD(K).LT.0.) THEN
              EPRD(K)=PRD(K)
              PRD(K)=0.
           END IF
           IF (PRDS(K).LT.0.) THEN
              EPRDS(K)=PRDS(K)
              PRDS(K)=0.
           END IF
           IF (PRDG(K).LT.0.) THEN
              EPRDG(K)=PRDG(K)
              PRDG(K)=0.
           END IF

!.......................................................................
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! CONSERVATION OF WATER
! THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
! ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.

! IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
! THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

! NOTE: CONSERVATION CHECK NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
! BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
! FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
! TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
! STEP

!****SENSITIVITY - NO ICE

      IF (ILIQ.EQ.1) THEN
      MNUCCC(K)=0.
      NNUCCC(K)=0.
      MNUCCR(K)=0.
      NNUCCR(K)=0.
      MNUCCD(K)=0.
      NNUCCD(K)=0.
      END IF

! ****SENSITIVITY - NO GRAUPEL
      IF (IGRAUP.EQ.1) THEN
            PRACG(K) = 0.
            PSACR(K) = 0.
	    PSACWG(K) = 0.
	    PRDG(K) = 0.
	    EPRDG(K) = 0.
            EVPMG(K) = 0.
            PGMLT(K) = 0.
	    NPRACG(K) = 0.
	    NPSACWG(K) = 0.
	    NSCNG(K) = 0.
 	    NGRACS(K) = 0.
	    NSUBG(K) = 0.
	    NGMLTG(K) = 0.
            NGMLTR(K) = 0.
! fix 053011
            PIACRS(K)=PIACRS(K)+PIACR(K)
            PIACR(K) = 0.
! fix 070713
	    PRACIS(K)=PRACIS(K)+PRACI(K)
	    PRACI(K) = 0.
	    PSACWS(K)=PSACWS(K)+PGSACW(K)
	    PGSACW(K) = 0.
	    PRACS(K)=PRACS(K)+PGRACS(K)
	    PGRACS(K) = 0.
       END IF

! CONSERVATION OF QC

      DUM = (PRC(K)+PRA(K)+MNUCCC(K)+PSACWS(K)+PSACWI(K)+QMULTS(K)+PSACWG(K)+PGSACW(K)+QMULTG(K))*DT

      IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have cloud water on Mars?", k
#else
        RATIO = QC3D(K)/DUM

        PRC(K) = PRC(K)*RATIO
        PRA(K) = PRA(K)*RATIO
        MNUCCC(K) = MNUCCC(K)*RATIO
        PSACWS(K) = PSACWS(K)*RATIO
        PSACWI(K) = PSACWI(K)*RATIO
        QMULTS(K) = QMULTS(K)*RATIO
        QMULTG(K) = QMULTG(K)*RATIO
        PSACWG(K) = PSACWG(K)*RATIO
	PGSACW(K) = PGSACW(K)*RATIO
#endif
        END IF
 
#ifdef WRF_MARS
!conservation of dust - This should be redundant because we're a good boy above in rate limiting things.
    DUM = MNUCCORE(K) * DT
    if(DUM .gt. qdust3d(k)) then
        RATIO = qdust3d(k)/ DUM
        MNUCCORE(K) = qdust3d(k)/DT ! MNUCCORE(K)*RATIO by definition
        NNUCCD(K) = ndust3d(k)/DT
        MNUCCD(K) = MNUCCD(K) * RATIO
    endif
#endif

! CONSERVATION OF QI

      DUM = (-PRD(K)-MNUCCC(K)+PRCI(K)+PRAI(K)-QMULTS(K)-QMULTG(K)-QMULTR(K)-QMULTRG(K) &
                -MNUCCD(K)+PRACI(K)+PRACIS(K)-EPRD(K)-PSACWI(K))*DT
!    write(*,*) k, prd(k), eprd(k), mnuccd(k), dum, qi3d(k)
#ifdef WRF_MARS
     IF (DUM.GT.QI3D(K).AND.(PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K)).GT.0) THEN
#else
      IF (DUM.GT.QI3D(K).AND.QI3D(K).GE.QSMALL) THEN
#endif
        RATIO = (QI3D(K)/DT+PRD(K)+MNUCCC(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+ &
                     MNUCCD(K)+PSACWI(K))/ &
                      (PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K))

        PRCI(K) = PRCI(K)*RATIO
        PRAI(K) = PRAI(K)*RATIO
        PRACI(K) = PRACI(K)*RATIO
        PRACIS(K) = PRACIS(K)*RATIO
#ifdef WRF_MARS
        if (ratio .eq. 0) then
            naked_core_removal=.true.
        else
            naked_core_removal=.false.
        endif
#endif
#ifdef MARS_MORR_DEBUG
      if(iprint) THEN
        CALL_DEBUG("EPRDRATIO",EPRD(K))
        CALL_DEBUG("RATIO",RATIO)
      endif
#endif
        EPRD(K) = EPRD(K)*RATIO

        END IF

! CONSERVATION OF QR

      DUM=((PRACS(K)-PRE(K))+(QMULTR(K)+QMULTRG(K)-PRC(K))+(MNUCCR(K)-PRA(K))+ &
             PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))*DT

      IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should we have rain on Mars?", k
#else

        RATIO = (QR3D(K)/DT+PRC(K)+PRA(K))/ &
             (-PRE(K)+QMULTR(K)+QMULTRG(K)+PRACS(K)+MNUCCR(K)+PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))

        PRE(K) = PRE(K)*RATIO
        PRACS(K) = PRACS(K)*RATIO
        QMULTR(K) = QMULTR(K)*RATIO
        QMULTRG(K) = QMULTRG(K)*RATIO
        MNUCCR(K) = MNUCCR(K)*RATIO
        PIACR(K) = PIACR(K)*RATIO
        PIACRS(K) = PIACRS(K)*RATIO
        PGRACS(K) = PGRACS(K)*RATIO
        PRACG(K) = PRACG(K)*RATIO

#endif
        END IF
#ifndef WRF_MARS
! CONSERVATION OF QNI
! CONSERVATION FOR GRAUPEL SCHEME

        IF (IGRAUP.EQ.0) THEN

      DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K))*DT

      IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

        RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K))/(-EPRDS(K)+PSACR(K))

       EPRDS(K) = EPRDS(K)*RATIO
       PSACR(K) = PSACR(K)*RATIO

       END IF

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
       ELSE IF (IGRAUP.EQ.1) THEN

      DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K)-MNUCCR(K))*DT

      IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

       RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))/(-EPRDS(K)+PSACR(K))

       EPRDS(K) = EPRDS(K)*RATIO
       PSACR(K) = PSACR(K)*RATIO

       END IF

       END IF

! CONSERVATION OF QG

      DUM = (-PSACWG(K)-PRACG(K)-PGSACW(K)-PGRACS(K)-PRDG(K)-MNUCCR(K)-EPRDG(K)-PIACR(K)-PRACI(K)-PSACR(K))*DT

      IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should we have graupel on Mars?", k
#else
        RATIO = (QG3D(K)/DT+PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PRDG(K)+MNUCCR(K)+PSACR(K)+&
                  PIACR(K)+PRACI(K))/(-EPRDG(K))

       EPRDG(K) = EPRDG(K)*RATIO
#endif
      END IF
#endif
! TENDENCIES

      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-PRD(K)-PRDS(K)-MNUCCD(K)-EPRD(K)-EPRDS(K)-PRDG(K)-EPRDG(K))

#ifdef WRF_MARS
      T3DTEN(K) = 0.0 !NO HEATING FROM MICROPHYSICS
#else
! BUG FIX HM, 3/1/11, INCLUDE PIACR AND PIACRS
      T3DTEN(K) = T3DTEN(K)+(PRE(K)                                 &
               *XXLV(K)+(PRD(K)+PRDS(K)+                            &
                MNUCCD(K)+EPRD(K)+EPRDS(K)+PRDG(K)+EPRDG(K))*XXLS(K)+         &
               (PSACWS(K)+PSACWI(K)+MNUCCC(K)+MNUCCR(K)+                      &
                QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+PRACS(K) &
                +PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PIACR(K)+PIACRS(K))*XLF(K))/CPM(K)

#endif

#ifdef WRF_MARS
!dust tendency
if (mars_heterogeneous) then !inuc.eq.6) then 
    if (mars_dust_scavenging) then
      QDUST3DTEN(K) = -MNUCCORE(K)
      NDUST3DTEN(K) = -NNUCCD(K)
    endif
!core
    QCORE3dTEN(k) = MNUCCORE(K)
endif
#endif
  !    write(*,*) "ice: ", k, QI3DTEN(K), PRD(K) ,EPRD(K), PSACWI(K), MNUCCC(K), PRCI(K), &
!                  PRAI(K), QMULTS(K), QMULTG(K), QMULTR(K), QMULTRG(K), MNUCCD(K), PRACI(K), PRACIS(K)

      QC3DTEN(K) = QC3DTEN(K)+                                      &
                 (-PRA(K)-PRC(K)-MNUCCC(K)+PCC(K)-                  &
                  PSACWS(K)-PSACWI(K)-QMULTS(K)-QMULTG(K)-PSACWG(K)-PGSACW(K))
      QI3DTEN(K) = QI3DTEN(K)+                                      &
         (PRD(K)+EPRD(K)+PSACWI(K)+MNUCCC(K)-PRCI(K)-                                 &
                  PRAI(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+MNUCCD(K)-PRACI(K)-PRACIS(K))
      QR3DTEN(K) = QR3DTEN(K)+                                      &
                 (PRE(K)+PRA(K)+PRC(K)-PRACS(K)-MNUCCR(K)-QMULTR(K)-QMULTRG(K) &
             -PIACR(K)-PIACRS(K)-PRACG(K)-PGRACS(K))

      IF (IGRAUP.EQ.0) THEN

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K))
      QG3DTEN(K) = QG3DTEN(K)+(PRACG(K)+PSACWG(K)+PGSACW(K)+PGRACS(K)+ &
                    PRDG(K)+EPRDG(K)+MNUCCR(K)+PIACR(K)+PRACI(K)+PSACR(K))
      NG3DTEN(K) = NG3DTEN(K)+(NSCNG(K)+NGRACS(K)+NNUCCR(K)+NIACR(K))

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
      ELSE IF (IGRAUP.EQ.1) THEN

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K)+NNUCCR(K))

      END IF

      NC3DTEN(K) = NC3DTEN(K)+(-NNUCCC(K)-NPSACWS(K)                &
            -NPRA(K)-NPRC(K)-NPSACWI(K)-NPSACWG(K))

      NI3DTEN(K) = NI3DTEN(K)+                                      &
       (NNUCCC(K)-NPRCI(K)-NPRAI(K)+NMULTS(K)+NMULTG(K)+NMULTR(K)+NMULTRG(K)+ &
               NNUCCD(K)-NIACR(K)-NIACRS(K))

      NR3DTEN(K) = NR3DTEN(K)+(NPRC1(K)-NPRACS(K)-NNUCCR(K)      &
                   +NRAGG(K)-NIACR(K)-NIACRS(K)-NPRACG(K)-NGRACS(K))

! HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC

	C2PREC(K) = PRA(K)+PRC(K)+PSACWS(K)+QMULTS(K)+QMULTG(K)+PSACWG(K)+ &
       PGSACW(K)+MNUCCC(K)+PSACWI(K)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION
#ifdef WRF_MARS
!not sure what goes here, everything should be covered in the ice condensation above.
#else
      DUMT = T3D(K)+DT*T3DTEN(K)
      DUMQV = QV3D(K)+DT*QV3DTEN(K)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(k),POLYSVP(DUMT,0))
      DUMQSS = EP_2*dum/(PRES(K)-dum)
      DUMQC = QC3D(K)+DT*QC3DTEN(K)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

      DUMS = DUMQV-DUMQSS
      PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
      IF (PCC(K)*DT+DUMQC.LT.0.) THEN
           PCC(K) = -DUMQC/DT
      END IF

      QV3DTEN(K) = QV3DTEN(K)-PCC(K)
      T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
      QC3DTEN(K) = QC3DTEN(K)+PCC(K)
#endif

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS
! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
! DROPLET CONCENTRATION IS SPECIFIED !!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC(K).LT.0.) THEN
!        DUM = PCC(K)*DT/QC3D(K)
!           DUM = MAX(-1.,DUM)
!        NSUBC(K) = DUM*NC3D(K)/DT
!     END IF
#ifdef WRF_MARS
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("EPRD",  eprd(k))
        CALL_DEBUG("naked_core_removal",  naked_core_removal)
#endif
      IF (EPRD(K).LT.0. .or. naked_core_removal) THEN
        if(qi3d(k).gt.0) then
         DUM = EPRD(K)*DT/QI3D(K)
            DUM = MAX(-1.,DUM)
        else
            DUM = -1
        endif
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("DUM",  DUM)
        CALL_DEBUG("NI3D", NI3D(k))
        CALL_DEBUG("NSUBI",  DUM*NI3D(k)/DT)                
#endif
         NSUBI(K) = DUM*NI3D(K)/DT
         if(mars_heterogeneous) then ! inuc.eq.6) then
             QSUBCORE(K) = DUM*qcore3d(k)/dt !core depletion
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("QSUBCORE",  QSUBCORE(K))
        CALL_DEBUG("QSUBDUST",  QSUBDUST(K))
        CALL_DEBUG("NSUBDUST",  NSUBDUST(K))
#endif
             if (mars_dust_scavenging) then
              QSUBDUST(K) = -QSUBCORE(K)      !dust growth
              NSUBDUST(K) = -NSUBI(K)         !dust growth
             endif
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("QSUBCORE",  QSUBCORE(K))
        CALL_DEBUG("QSUBDUST",  QSUBDUST(K))
        CALL_DEBUG("NSUBDUST",  NSUBDUST(K))
#endif
         endif

      END IF
#else
      IF (EPRD(K).LT.0.) THEN
         DUM = EPRD(K)*DT/QI3D(K)
            DUM = MAX(-1.,DUM)
         NSUBI(K) = DUM*NI3D(K)/DT
      ENDIF
#endif
      IF (EPRDS(K).LT.0.) THEN
         DUM = EPRDS(K)*DT/QNI3D(K)
           DUM = MAX(-1.,DUM)
         NSUBS(K) = DUM*NS3D(K)/DT
      END IF
      IF (PRE(K).LT.0.) THEN
         DUM = PRE(K)*DT/QR3D(K)
           DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
      END IF
      IF (EPRDG(K).LT.0.) THEN
         DUM = EPRDG(K)*DT/QG3D(K)
           DUM = MAX(-1.,DUM)
         NSUBG(K) = DUM*NG3D(K)/DT
      END IF

!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

! UPDATE TENDENCIES
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("nsubi",  nsubi(k))
        CALL_DEBUG("ni3dten",  ni3dten(k))
#endif
!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)
         NI3DTEN(K) = NI3DTEN(K)+NSUBI(K)
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("nsubi",  nsubi(k))
        CALL_DEBUG("ni3dten",  ni3dten(k))
#endif
         NS3DTEN(K) = NS3DTEN(K)+NSUBS(K)
         NG3DTEN(K) = NG3DTEN(K)+NSUBG(K)
         NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)
#ifdef WRF_MARS
         NDUST3DTEN(K) = NDUST3DTEN(K) + NSUBDUST(K)
         QDUST3DTEN(K) = QDUST3DTEN(K) + QSUBDUST(K)
         QCORE3DTEN(K) = QCORE3DTEN(K) + QSUBCORE(K)
#endif

#ifndef WRF_MARS
         END IF !!!!!! TEMPERATURE
#else
!no T>273 on Mars because of partial pressures.
#endif

! SWITCH LTRUE TO 1, SINCE HYDROMETEORS ARE PRESENT
         LTRUE = 1

 200     CONTINUE

        END DO

! INITIALIZE PRECIP AND SNOW RATES
      PRECRT = 0.
      SNOWRT = 0.
! hm added 7/13/13
      SNOWPRT = 0.
      GRPLPRT = 0.
#ifdef WRF_MARS
      PRECDUSTRT = 0
      PRECNDUSTRT = 0
      PRECCORERT = 0
#endif

! IF THERE ARE NO HYDROMETEORS, THEN SKIP TO END OF SUBROUTINE

        IF (LTRUE.EQ.0) GOTO 400

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!.......................................................................
! CALCULATE SEDIMENATION
! THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
! FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
! STABILITY, I.E. COURANT# < 1

!.......................................................................
#ifdef WRF_MARS
!    water_before = sum(rho*dzq*(qv3d+qv3dten*dt+qi3d+)) !no rain
    
    
#endif
      NSTEP = 1

      DO K = KTE,KTS,-1
  iprint=.false.
          if(jprint.and.k.eq.KVALUE) iprint=.true.
#ifdef MARS_MORR_DEBUG
            CALL_DEBUG("dumfni",0)
            CALL_DEBUG("ni3d",ni3d(k))
            CALL_DEBUG("ni3dten",ni3dten(k))
            CALL_DEBUG("dt",dt)
#endif
        DUMI(K) = QI3D(K)+QI3DTEN(K)*DT
        DUMQS(K) = QNI3D(K)+QNI3DTEN(K)*DT
        DUMR(K) = QR3D(K)+QR3DTEN(K)*DT
        DUMFNI(K) = NI3D(K)+NI3DTEN(K)*DT
        DUMFNS(K) = NS3D(K)+NS3DTEN(K)*DT
        DUMFNR(K) = NR3D(K)+NR3DTEN(K)*DT
        DUMC(K) = QC3D(K)+QC3DTEN(K)*DT
        DUMFNC(K) = NC3D(K)+NC3DTEN(K)*DT
	DUMG(K) = QG3D(K)+QG3DTEN(K)*DT
	DUMFNG(K) = NG3D(K)+NG3DTEN(K)*DT

#ifdef WRF_MARS
        !dust
        DUMD(K) = QDUST3D(K) + QDUST3DTEN(K)*DT
        DUMFND(K) = NDUST3D(K) + NDUST3DTEN(K)*DT
        DUMCORE(K) = QCORE3D(K) + QCORE3DTEN(K)*DT
       ! ni3d_temp = ni3d(k) + ni3dten(k)*dt
#endif

! SWITCH FOR CONSTANT DROPLET NUMBER
        IF (iinum.EQ.1) THEN
        DUMFNC(K) = NC3D(K)
        END IF

! GET DUMMY LAMDA FOR SEDIMENTATION CALCULATIONS

! MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
      DUMFNI(K) = MAX(0.,DUMFNI(K))
      DUMFNS(K) = MAX(0.,DUMFNS(K))
      DUMFNC(K) = MAX(0.,DUMFNC(K))
      DUMFNR(K) = MAX(0.,DUMFNR(K))
      DUMFNG(K) = MAX(0.,DUMFNG(K))
#ifdef WRF_MARS
      DUMFND(K) = MAX(0.,DUMFND(K))
#endif

!......................................................................
! CLOUD ICE
#ifdef WRF_MARS
    ni3d_temp = DUMFNI(k)!ni3d(k)
    if(.not.mars_heterogeneous) then !inuc.ne.6) then !homo
        ice_rho(k) = rhoi
        qicetot3d(k) = dumi(k)
        if(dumi(k).ge.qsmall) then 
            call constrain_mars_waterice_homogeneous(dumi(k), ni3d_temp, n0i_temp, dlami)
        endif
    else
        qicetot3d(k) = dumcore(k) + dumi(k)
        if(qicetot3d(K).GE.QSMALL) then
            call mean_ice_density(dumcore(k), dumi(k), ice_rho(k))
#ifdef MARS_MORR_DEBUG
            CALL_DEBUG("call constrain",__LINE__)
#endif
            call constrain_mars_waterice_heterogeneous(dumi(k), ni3d_temp, dumcore(k), ice_rho(k), n0i_temp, dlami,thisi, thisj, k)
        endif
    endif            
!        if(inuc.ne.6) then
!            ice_rho(k) = rhoi
!            qicetot3d(k) = dumi(k)
!            call constrain_mars_waterice_homogeneous(qi3d(k), ni3d_temp, n0i_temp, dlami)
!        else
!            call mean_ice_density(dumcore(k), dumi(k), ice_rho(k))
!            qicetot3d(k) = dumcore(k) + dumi(k)
!            call constrain_mars_waterice_heterogeneous(dumi(k), ni3d_temp, dumcore(k), ice_rho(k), n0i_temp, dlami)
!        endif
!old method        DLAMI = (ICE_CONS12*DUMFNI(K)/DUMI(K))**(1./DI)
!old method        DLAMI=MAX(DLAMI,LAMMINI)
!old method        DLAMI=MIN(DLAMI,LAMMAXI)
#else
      IF (DUMI(K).GE.QSMALL) THEN
        DLAMI = (CONS12*DUMFNI(K)/DUMI(K))**(1./DI)
        DLAMI=MAX(DLAMI,LAMMINI)
        DLAMI=MIN(DLAMI,LAMMAXI)
      END IF

#endif
!......................................................................
! RAIN

      IF (DUMR(K).GE.QSMALL) THEN
        DLAMR = (PI*RHOW*DUMFNR(K)/DUMR(K))**(1./3.)
        DLAMR=MAX(DLAMR,LAMMINR)
        DLAMR=MIN(DLAMR,LAMMAXR)
      END IF
!......................................................................
! CLOUD DROPLETS

      IF (DUMC(K).GE.QSMALL) THEN
         DUM = PRES(K)/(287.15*T3D(K))
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

        DLAMC = (CONS26*DUMFNC(K)*GAMMA(PGAM(K)+4.)/(DUMC(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
        LAMMIN = (PGAM(K)+1.)/60.E-6
        LAMMAX = (PGAM(K)+1.)/1.E-6
        DLAMC=MAX(DLAMC,LAMMIN)
        DLAMC=MIN(DLAMC,LAMMAX)
      END IF
!......................................................................
! SNOW

      IF (DUMQS(K).GE.QSMALL) THEN
        DLAMS = (CONS1*DUMFNS(K)/ DUMQS(K))**(1./DS)
        DLAMS=MAX(DLAMS,LAMMINS)
        DLAMS=MIN(DLAMS,LAMMAXS)
      END IF
!......................................................................
! GRAUPEL

      IF (DUMG(K).GE.QSMALL) THEN
        DLAMG = (CONS2*DUMFNG(K)/ DUMG(K))**(1./DG)
        DLAMG=MAX(DLAMG,LAMMING)
        DLAMG=MIN(DLAMG,LAMMAXG)
      END IF

#ifdef WRF_MARS
!DUST
    if(DUMD(K).ge.QSMALL) THEN
        DLAMD = (DUST_CONS1 * DUMFND(k) / DUMD(k)) ** (1./3)
        DLAMD = MAX(DLAMD,LAMMIND)
        DLAMD = MIN(DLAMD,LAMMAXD)
    endif
#endif
!......................................................................
! CALCULATE NUMBER-WEIGHTED AND MASS-WEIGHTED TERMINAL FALL SPEEDS

! CLOUD WATER

      IF (DUMC(K).GE.QSMALL) THEN
      UNC =  ACN(K)*GAMMA(1.+BC+PGAM(K))/ (DLAMC**BC*GAMMA(PGAM(K)+1.))
      UMC = ACN(K)*GAMMA(4.+BC+PGAM(K))/  (DLAMC**BC*GAMMA(PGAM(K)+4.))
      ELSE
      UMC = 0.
      UNC = 0.
      END IF

#ifndef WRF_MARS
      IF (DUMI(K).GE.QSMALL) THEN
      UNI =  AIN(K)*CONS27/DLAMI**BI
      UMI = AIN(K)*CONS28/(DLAMI**BI)
      ELSE
      UMI = 0.
      UNI = 0.
      END IF

#endif
      IF (DUMR(K).GE.QSMALL) THEN
      UNR = ARN(K)*CONS6/DLAMR**BR
      UMR = ARN(K)*CONS4/(DLAMR**BR)
      ELSE
      UMR = 0.
      UNR = 0.
      END IF

      IF (DUMQS(K).GE.QSMALL) THEN
      UMS = ASN(K)*CONS3/(DLAMS**BS)
      UNS = ASN(K)*CONS5/DLAMS**BS
      ELSE
      UMS = 0.
      UNS = 0.
      END IF

      IF (DUMG(K).GE.QSMALL) THEN
      UMG = AGN(K)*CONS7/(DLAMG**BG)
      UNG = AGN(K)*CONS8/DLAMG**BG
      ELSE
      UMG = 0.
      UNG = 0.
      END IF
#ifdef WRF_MARS
    !mean free path = 1./ (sqrt(2) * pi * n * (2r)**2)
    !               = 1/(4sqrt(2)*pi*d**2) * molmass / rho/ Na
    !               = dust_cons3 / RHO
    !fall speeds
!CLOUD ICE
    if(qicetot3d(K).ge.QSMALL) THEN
!        freepath(k) = 
        tf = 2.*dust_cons3 / rho(K)
        invDLAMI = 1./DLAMI
        invtf=1./tf
        se_const = log((1 + DUST_CONS_E_SEDIM*invDLAMI*invtf )) !first calculate the fraction
        se_m = exp ( -(5+ICE_MU) * se_const ) !then exponentiate it
        se_n = exp ( -(2+ICE_MU) * se_const ) !then exponentiate it
        UMI = AIN(K) * ice_rho(k) * ((ICE_MU+4)*invDLAMI   ) * ((ICE_MU+5)*invDLAMI + DUST_CONS_A_SEDIM*tf + DUST_CONS_B_SEDIM*tf * se_m )
        UNI = AIN(K) * ice_rho(k) * ((ICE_MU+1)*invDLAMI   ) * ((ICE_MU+2)*invDLAMI + DUST_CONS_A_SEDIM*tf + DUST_CONS_B_SEDIM*tf * se_n )
!        write(*,*) "vel: ", ain(k), ice_rho(k), ((ICE_MU+4)*invDLAMI   ) * ((ICE_MU+5)*invDLAMI + DUST_CONS_A_SEDIM*tf + DUST_CONS_B_SEDIM*tf * se_m ),  ((ICE_MU+1)*invDLAMI   ) * ((ICE_MU+2)*invDLAMI + DUST_CONS_A_SEDIM*tf + DUST_CONS_B_SEDIM*tf * se_n )

        !if core mass, set speeds to ice speeds
        if(qcore3d(k).gt.qsmall) then
            UMCORE = UMI
        else
            UMCORE = 0.
        endif
        !but if no ice mass, zero out ice speeds
        if(qi3d(k).lt.qsmall) then
            UNI = 0.0
            UMI = 0.0
        endif
    else
        !would be here if no core or ice.
        UMI = 0.
        UNI = 0.
        UMCORE = 0.0
    endif
    
    
!write(*,*) "spped", k, dumi(k), umi, uni
!DUST
    if(DUMD(K).ge.QSMALL) THEN
    if(DISTRIBUTION_MARS_DUST == 0) then !GAMMA
    
    !integating over particle size:
    !assuming C = 1 + Kn * (a_dust + b_dust * exp(-e_dust / Kn))
    !Kn = D / (2*freepath)
    !fun
        freepath(k) = dust_cons3 / rho(K)
        tf = 2.*freepath(K)
        invtf = 1./tf
        invDLAMD = 1./DLAMD
        se_const = log((1 + DUST_CONS_E_SEDIM*invDLAMD*invtf )) !first calculate the fraction
        se_m = exp ( -(5+DUST_MU) * se_const ) !then exponentiate it
        se_n = exp ( -(2+DUST_MU) * se_const ) !then exponentiate it
        
        UMD = ADN(K) * ((DUST_MU+4)*invDLAMD   ) * ((DUST_MU+5)*invDLAMD + DUST_CONS_A_SEDIM*tf + DUST_CONS_B_SEDIM*tf * se_m )
        UND = ADN(K) * ((DUST_MU+1)*invDLAMD   ) * ((DUST_MU+2)*invDLAMD + DUST_CONS_A_SEDIM*tf + DUST_CONS_B_SEDIM*tf * se_n )
        
    else if(DISTRIBUTION_MARS_DUST == 1) then !lognormal
    
        tf = 2.*freepath(K)
        invtf = 1./tf
        
        se_m = exp(lamd(K) + 4*sigma_squared)
        se_n = exp(lamd(K) + sigma_squared)
        
        umd = adn(k) * se_m * (se_m + tf*DUST_CONS_A_SEDIM*dust_cons2)
        und = adn(k) * se_n * (se_n + tf*DUST_CONS_A_SEDIM*dust_cons2)
    
    endif
    ELSE
        UMD = 0.
        UND = 0.
    endif
#endif
! SET REALISTIC LIMITS ON FALLSPEED

! bug fix, 10/08/09
        dum=(rhosu/rho(k))**0.54
        UMS=MIN(UMS,1.2*dum)
        UNS=MIN(UNS,1.2*dum)
! fix 053011
! fix for correction by AA 4/6/11
#ifdef WRF_MARS
    !what values are appropriate for Mars?
    umi=min(umi,10.)
    uni=min(uni,10.)
    umcore=min(umcore,10.)
#else
        UMI=MIN(UMI,1.2*(rhosu/rho(k))**0.35)
        UNI=MIN(UNI,1.2*(rhosu/rho(k))**0.35)
#endif
        UMR=MIN(UMR,9.1*dum)
        UNR=MIN(UNR,9.1*dum)
        UMG=MIN(UMG,20.*dum)
        UNG=MIN(UNG,20.*dum)
#ifdef WRF_MARS
    umd=min(umd,10.)
    und=min(und,10.)
!        UMD=min(UMD,(ADN*2/DLAMD**2) *10. What's a reasonable limit?)
!        UND=min(UMD,(ADN*20/DLAMD**2) *10. These values are 10 times the uncorrected Stokes velocity)
#endif
      FR(K) = UMR
      FI(K) = UMI
      FNI(K) = UNI
      FS(K) = UMS
      FNS(K) = UNS
      FNR(K) = UNR
      FC(K) = UMC
      FNC(K) = UNC
      FG(K) = UMG
      FNG(K) = UNG
#ifdef WRF_MARS
      FD(K) = UMD
      FND(K) = UND
      FCORE(K) = UMCORE
#endif

! V3.3 MODIFY FALLSPEED BELOW LEVEL OF PRECIP

	IF (K.LE.KTE-1) THEN
        IF (FR(K).LT.1.E-10) THEN
	FR(K)=FR(K+1)
	END IF
        IF (FI(K).LT.1.E-10) THEN
	FI(K)=FI(K+1)
	END IF
        IF (FNI(K).LT.1.E-10) THEN
	FNI(K)=FNI(K+1)
	END IF
        IF (FS(K).LT.1.E-10) THEN
	FS(K)=FS(K+1)
	END IF
        IF (FNS(K).LT.1.E-10) THEN
	FNS(K)=FNS(K+1)
	END IF
        IF (FNR(K).LT.1.E-10) THEN
	FNR(K)=FNR(K+1)
	END IF
        IF (FC(K).LT.1.E-10) THEN
	FC(K)=FC(K+1)
	END IF
        IF (FNC(K).LT.1.E-10) THEN
	FNC(K)=FNC(K+1)
	END IF
        IF (FG(K).LT.1.E-10) THEN
	FG(K)=FG(K+1)
	END IF
        IF (FNG(K).LT.1.E-10) THEN
	FNG(K)=FNG(K+1)
	END IF
#ifdef WRF_MARS
        IF (FD(K).LT.1.E-10) THEN
    FD(K)=FD(K+1)
    END IF
        IF (FND(K).LT.1.E-10) THEN
    FND(K)=FND(K+1)
    END IF
        IF (FCORE(K).LT.1.E-10) THEN
    FCORE(K)=FCORE(K+1)
    END IF
#endif
	END IF ! K LE KTE-1

! CALCULATE NUMBER OF SPLIT TIME STEPS

      RGVM = MAX(FR(K),FI(K),FS(K),FC(K),FNI(K),FNR(K),FNS(K),FNC(K),FG(K),FNG(K))
#ifdef WRF_MARS
!include dust in CFL calculation
      RGVM = MAX(RGVM, FD(K),FND(K),FCORE(K))
#endif
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEP = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEP)

! MULTIPLY VARIABLES BY RHO
      DUMR(k) = DUMR(k)*RHO(K)
      DUMI(k) = DUMI(k)*RHO(K)
      DUMFNI(k) = DUMFNI(K)*RHO(K)
      DUMQS(k) = DUMQS(K)*RHO(K)
      DUMFNS(k) = DUMFNS(K)*RHO(K)
      DUMFNR(k) = DUMFNR(K)*RHO(K)
      DUMC(k) = DUMC(K)*RHO(K)
      DUMFNC(k) = DUMFNC(K)*RHO(K)
      DUMG(k) = DUMG(K)*RHO(K)
      DUMFNG(k) = DUMFNG(K)*RHO(K)
#ifdef WRF_MARS
       DUMD(k) = DUMD(k)*RHO(k)
      DUMFND(k) = DUMFND(k)*RHO(k)
       DUMCORE(k) = DUMCORE(k)*RHO(k)
#endif
      END DO

      DO N = 1,NSTEP

      DO K = KTS,KTE
      FALOUTR(K) = FR(K)*DUMR(K)
      FALOUTI(K) = FI(K)*DUMI(K)
      FALOUTNI(K) = FNI(K)*DUMFNI(K)
      FALOUTS(K) = FS(K)*DUMQS(K)
      FALOUTNS(K) = FNS(K)*DUMFNS(K)
      FALOUTNR(K) = FNR(K)*DUMFNR(K)
      FALOUTC(K) = FC(K)*DUMC(K)
      FALOUTNC(K) = FNC(K)*DUMFNC(K)
      FALOUTG(K) = FG(K)*DUMG(K)
      FALOUTNG(K) = FNG(K)*DUMFNG(K)
#ifdef WRF_MARS
      FALOUTD(K) = FD(K)*DUMD(K)
      FALOUTND(K) = FND(K)*DUMFND(K)
      FALOUTCORE(K) = FCORE(K)*DUMCORE(K)
#endif
      END DO

! TOP OF MODEL

      K = KTE
      FALTNDR = FALOUTR(K)/DZQ(k)
      FALTNDI = FALOUTI(K)/DZQ(k)
      FALTNDNI = FALOUTNI(K)/DZQ(k)
      FALTNDS = FALOUTS(K)/DZQ(k)
      FALTNDNS = FALOUTNS(K)/DZQ(k)
      FALTNDNR = FALOUTNR(K)/DZQ(k)
      FALTNDC = FALOUTC(K)/DZQ(k)
      FALTNDNC = FALOUTNC(K)/DZQ(k)
      FALTNDG = FALOUTG(K)/DZQ(k)
      FALTNDNG = FALOUTNG(K)/DZQ(k)
#ifdef WRF_MARS
      FALTNDD = FALOUTD(K)/DZQ(K)
      FALTNDND = FALOUTND(K)/DZQ(K)
      FALTNDCORE = FALOUTCORE(K)/DZQ(K)
#endif

! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

      QRSTEN(K) = QRSTEN(K)-FALTNDR/NSTEP/RHO(k)
      QISTEN(K) = QISTEN(K)-FALTNDI/NSTEP/RHO(k)
      NISTEN(K) = NISTEN(K)-FALTNDNI/NSTEP/RHO(k)
      QNISTEN(K) = QNISTEN(K)-FALTNDS/NSTEP/RHO(k)
      NS3DTEN(K) = NS3DTEN(K)-FALTNDNS/NSTEP/RHO(k)
      NR3DTEN(K) = NR3DTEN(K)-FALTNDNR/NSTEP/RHO(k)
      QCSTEN(K) = QCSTEN(K)-FALTNDC/NSTEP/RHO(k)
      NC3DTEN(K) = NC3DTEN(K)-FALTNDNC/NSTEP/RHO(k)
      QGSTEN(K) = QGSTEN(K)-FALTNDG/NSTEP/RHO(k)
      NG3DTEN(K) = NG3DTEN(K)-FALTNDNG/NSTEP/RHO(k)
#ifdef WRF_MARS
      QDUSTSTEN(K) = QDUSTSTEN(K)-FALTNDD/NSTEP/RHO(k)
      NDUSTSTEN(K) = NDUSTSTEN(K)-FALTNDND/NSTEP/RHO(k)
      QCORESTEN(K) = QCORESTEN(K)-FALTNDCORE/NSTEP/RHO(k)
#endif

      DUMR(K) = DUMR(K)-FALTNDR*DT/NSTEP
      DUMI(K) = DUMI(K)-FALTNDI*DT/NSTEP
      DUMFNI(K) = DUMFNI(K)-FALTNDNI*DT/NSTEP
      DUMQS(K) = DUMQS(K)-FALTNDS*DT/NSTEP
      DUMFNS(K) = DUMFNS(K)-FALTNDNS*DT/NSTEP
      DUMFNR(K) = DUMFNR(K)-FALTNDNR*DT/NSTEP
      DUMC(K) = DUMC(K)-FALTNDC*DT/NSTEP
      DUMFNC(K) = DUMFNC(K)-FALTNDNC*DT/NSTEP
      DUMG(K) = DUMG(K)-FALTNDG*DT/NSTEP
      DUMFNG(K) = DUMFNG(K)-FALTNDNG*DT/NSTEP
#ifdef WRF_MARS
      DUMD(K) = DUMD(K)-FALTNDD*DT/NSTEP
      DUMFND(K) = DUMFND(K)-FALTNDND*DT/NSTEP
      DUMCORE(K) = DUMCORE(K) - FALTNDCORE*DT/NSTEP
#endif

      DO K = KTE-1,KTS,-1
      FALTNDR = (FALOUTR(K+1)-FALOUTR(K))/DZQ(K)
      FALTNDI = (FALOUTI(K+1)-FALOUTI(K))/DZQ(K)
      FALTNDNI = (FALOUTNI(K+1)-FALOUTNI(K))/DZQ(K)
      FALTNDS = (FALOUTS(K+1)-FALOUTS(K))/DZQ(K)
      FALTNDNS = (FALOUTNS(K+1)-FALOUTNS(K))/DZQ(K)
      FALTNDNR = (FALOUTNR(K+1)-FALOUTNR(K))/DZQ(K)
      FALTNDC = (FALOUTC(K+1)-FALOUTC(K))/DZQ(K)
      FALTNDNC = (FALOUTNC(K+1)-FALOUTNC(K))/DZQ(K)
      FALTNDG = (FALOUTG(K+1)-FALOUTG(K))/DZQ(K)
      FALTNDNG = (FALOUTNG(K+1)-FALOUTNG(K))/DZQ(K)
#ifdef WRF_MARS
      FALTNDD = (FALOUTD(K+1)-FALOUTD(K))/DZQ(K)
      FALTNDND = (FALOUTND(K+1)-FALOUTND(K))/DZQ(K)
      FALTNDCORE = (FALOUTCORE(K+1)-FALOUTCORE(K))/DZQ(K)
#endif
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

      QRSTEN(K) = QRSTEN(K)+FALTNDR/NSTEP/RHO(k)
      QISTEN(K) = QISTEN(K)+FALTNDI/NSTEP/RHO(k)
      NISTEN(K) = NISTEN(K)+FALTNDNI/NSTEP/RHO(k)
      QNISTEN(K) = QNISTEN(K)+FALTNDS/NSTEP/RHO(k)
      NS3DTEN(K) = NS3DTEN(K)+FALTNDNS/NSTEP/RHO(k)
      NR3DTEN(K) = NR3DTEN(K)+FALTNDNR/NSTEP/RHO(k)
      QCSTEN(K) = QCSTEN(K)+FALTNDC/NSTEP/RHO(k)
      NC3DTEN(K) = NC3DTEN(K)+FALTNDNC/NSTEP/RHO(k)
      QGSTEN(K) = QGSTEN(K)+FALTNDG/NSTEP/RHO(k)
      NG3DTEN(K) = NG3DTEN(K)+FALTNDNG/NSTEP/RHO(k)
#ifdef WRF_MARS
      QDUSTSTEN(K) = QDUSTSTEN(K)+FALTNDD/NSTEP/RHO(k)
      NDUSTSTEN(K) = NDUSTSTEN(K)+FALTNDND/NSTEP/RHO(k)
      QCORESTEN(K) = QCORESTEN(K)+FALTNDCORE/NSTEP/RHO(k)
#endif

      DUMR(K) = DUMR(K)+FALTNDR*DT/NSTEP
      DUMI(K) = DUMI(K)+FALTNDI*DT/NSTEP
      DUMFNI(K) = DUMFNI(K)+FALTNDNI*DT/NSTEP
      DUMQS(K) = DUMQS(K)+FALTNDS*DT/NSTEP
      DUMFNS(K) = DUMFNS(K)+FALTNDNS*DT/NSTEP
      DUMFNR(K) = DUMFNR(K)+FALTNDNR*DT/NSTEP
      DUMC(K) = DUMC(K)+FALTNDC*DT/NSTEP
      DUMFNC(K) = DUMFNC(K)+FALTNDNC*DT/NSTEP
      DUMG(K) = DUMG(K)+FALTNDG*DT/NSTEP
      DUMFNG(K) = DUMFNG(K)+FALTNDNG*DT/NSTEP
#ifdef WRF_MARS
      DUMD(K) = DUMD(K)+FALTNDD*DT/NSTEP
      DUMFND(K) = DUMFND(K)+FALTNDND*DT/NSTEP
      DUMCORE(K) = DUMCORE(K)+FALTNDCORE*DT/NSTEP
#endif

! FOR WRF-CHEM, NEED PRECIP RATES (UNITS OF KG/M^2/S)
	  CSED(K)=CSED(K)+FALOUTC(K)/NSTEP
	  ISED(K)=ISED(K)+FALOUTI(K)/NSTEP
	  SSED(K)=SSED(K)+FALOUTS(K)/NSTEP
	  GSED(K)=GSED(K)+FALOUTG(K)/NSTEP
	  RSED(K)=RSED(K)+FALOUTR(K)/NSTEP
      END DO

! GET PRECIPITATION AND SNOWFALL ACCUMULATION DURING THE TIME STEP
! FACTOR OF 1000 CONVERTS FROM M TO MM, BUT DIVISION BY DENSITY
! OF LIQUID WATER CANCELS THIS FACTOR OF 1000

        PRECRT = PRECRT+(FALOUTR(KTS)+FALOUTC(KTS)+FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))  &
                     *DT/NSTEP
        SNOWRT = SNOWRT+(FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))*DT/NSTEP
! hm added 7/13/13
        SNOWPRT = SNOWPRT+(FALOUTI(KTS)+FALOUTS(KTS))*DT/NSTEP
        GRPLPRT = GRPLPRT+(FALOUTG(KTS))*DT/NSTEP

#ifdef WRF_MARS
      PRECDUSTRT = PRECDUSTRT + FALOUTD(KTS)*DT/NSTEP
      PRECNDUSTRT = PRECNDUSTRT + FALOUTND(KTS)*DT/NSTEP
      PRECCORERT = PRECCORERT + FALOUTCORE(KTS)*DT/NSTEP
#endif

      END DO

        DO K=KTS,KTE
  iprint=.false.
          if(jprint.and.k.eq.KVALUE) iprint=.true.
#ifdef MARS_MORR_DEBUG
            CALL_DEBUG("sediment",0)
            CALL_DEBUG("s ni3d",ni3d(k))
            CALL_DEBUG("s ni3dten",ni3dten(k))
            CALL_DEBUG("QISTEN",QISTEN(k))
            CALL_DEBUG("NISTEN",NISTEN(k))
            CALL_DEBUG("QCORESTEN",QCORESTEN(k))
            CALL_DEBUG("QDUSTSTEN",QDUSTSTEN(k))
            CALL_DEBUG("NDUSTSTEN",NDUSTSTEN(k))
            
            CALL_DEBUG("QI3DTEN",QI3DTEN(k))
            CALL_DEBUG("NI3DTEN",NI3DTEN(k)) 
            CALL_DEBUG("QCORE3DTEN",QCORE3DTEN(k)) 
            CALL_DEBUG("QDUST3DTEN",QDUST3DTEN(k))
            CALL_DEBUG("NDUST3DTEN",NDUST3DTEN(k))

#endif
! ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES
        QR3DTEN(K)=QR3DTEN(K)+QRSTEN(K)
        QI3DTEN(K)=QI3DTEN(K)+QISTEN(K)
        NI3DTEN(K)=NI3DTEN(K)+NISTEN(K)
        QC3DTEN(K)=QC3DTEN(K)+QCSTEN(K)
        QG3DTEN(K)=QG3DTEN(K)+QGSTEN(K)
        QNI3DTEN(K)=QNI3DTEN(K)+QNISTEN(K)
#ifdef WRF_MARS
        QDUST3DTEN(K)=QDUST3DTEN(K)+QDUSTSTEN(K)
        QCORE3DTEN(K)=QCORE3DTEN(K)+QCORESTEN(K)
        NDUST3DTEN(K)=NDUST3DTEN(K)+NDUSTSTEN(K)
#endif

! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs
#ifndef WRF_MARS
!hm 4/7/09 bug fix
!        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15) THEN
        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15.AND.LAMI(K).GE.1.E-10) THEN
        IF (1./LAMI(K).GE.2.*DCS) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should we have snow on Mars?", k
#else
           QNI3DTEN(K) = QNI3DTEN(K)+QI3D(K)/DT+ QI3DTEN(K)
           NS3DTEN(K) = NS3DTEN(K)+NI3D(K)/DT+   NI3DTEN(K)
           QI3DTEN(K) = -QI3D(K)/DT
           NI3DTEN(K) = -NI3D(K)/DT
#endif
        END IF
        END IF
#endif
! hm add tendencies here, then call sizeparameter
! to ensure consisitency between mixing ratio and number concentration

          QC3D(k)        = QC3D(k)+QC3DTEN(k)*DT
          QI3D(k)        = QI3D(k)+QI3DTEN(k)*DT
          QNI3D(k)        = QNI3D(k)+QNI3DTEN(k)*DT
          QR3D(k)        = QR3D(k)+QR3DTEN(k)*DT
          NC3D(k)        = NC3D(k)+NC3DTEN(k)*DT
          NI3D(k)        = NI3D(k)+NI3DTEN(k)*DT
          NS3D(k)        = NS3D(k)+NS3DTEN(k)*DT
          NR3D(k)        = NR3D(k)+NR3DTEN(k)*DT

          IF (IGRAUP.EQ.0) THEN
          QG3D(k)        = QG3D(k)+QG3DTEN(k)*DT
          NG3D(k)        = NG3D(k)+NG3DTEN(k)*DT
          END IF
#ifdef WRF_MARS
          qdust3d(k) = qdust3d(k)+QDUST3DTEN(K)*DT
          ndust3d(k) = ndust3d(k)+NDUST3DTEN(K)*DT
          qcore3d(k) = qcore3d(k)+QCORE3DTEN(K)*DT
          qicetot3d(k) = qcore3d(k) + qi3d(k)
#endif

! ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
          T3D(K)         = T3D(K)+T3DTEN(k)*DT
          QV3D(K)        = QV3D(K)+QV3DTEN(k)*DT

! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
#ifdef WRF_MARS
            EVS(K) = esat(T3D(K)) ! PA
            EIS(K) = EVS(K)       ! PA
    !total fraction for mars, not dry fraction to 
    !prevent extreme values where water vapor pressure equals or exceeds pressure
            QVS(K) = EP_2*EVS(K)/(PRES(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K))
#else
            EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
            EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))
#endif
            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! hm 7/9/09 change limit to 1.e-8
#ifdef WRF_MARS
!not sure if the limits are correct here
#else
             IF (QVQVS(K).LT.0.9) THEN
               IF (QR3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QR3D(K)
                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  QR3D(K)=0.
               END IF
               IF (QC3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QC3D(K)
                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  QC3D(K)=0.
               END IF
             END IF

             IF (QVQVSI(K).LT.0.9) THEN
               IF (QI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF
               IF (QNI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QNI3D(K)
                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  QNI3D(K)=0.
               END IF
               IF (QG3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QG3D(K)
                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  QG3D(K)=0.
               END IF
             END IF
#endif
!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(K).LT.QSMALL) THEN
         QC3D(K) = 0.
         NC3D(K) = 0.
         EFFC(K) = 0.
       END IF
       IF (QR3D(K).LT.QSMALL) THEN
         QR3D(K) = 0.
         NR3D(K) = 0.
         EFFR(K) = 0.
       END IF
!               write(*,*) __LINE__, k, qi3d(k), QSMALL
#ifdef WRF_MARS
!We do want to zero out the ice, but we can't do it here because it kills the condensation
!onto a dust core, so we check later for heterogeneous nucleation
    if (.not.mars_heterogeneous) then !inuc.ne.6) then 
#endif
       IF (QI3D(K).LT.QSMALL) THEN
         QI3D(K) = 0.
         NI3D(K) = 0.
         EFFI(K) = 0.
       END IF
#ifdef WRF_MARS
    endif
#endif
       IF (QNI3D(K).LT.QSMALL) THEN
         QNI3D(K) = 0.
         NS3D(K) = 0.
         EFFS(K) = 0.
       END IF
       IF (QG3D(K).LT.QSMALL) THEN
         QG3D(K) = 0.
         NG3D(K) = 0.
         EFFG(K) = 0.
       END IF
#ifdef WRF_MARS
       IF (QDUST3D(K).LT.QSMALL) THEN
           QDUST3D(K) = 0.
           NDUST3D(K) = 0.
           EFFD(K) = 0.
       ENDIF
       IF (QCORE3D(K).LT.QSMALL .and. QI3D(K).LT.QSMALL) THEN
           QCORE3D(K) = 0.
           QI3D(K) = 0.
           NI3D(K) = 0.
           EFFI(K) = 0.
       ELSE IF(QCORE3D(K).lt.QSMALL) THEN
           QCORE3D(K) = 0.0
       ELSE IF(QI3D(K).lt.QSMALL) THEN
           QI3D(K)=0.0
       ENDIF

#endif
!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, THEN SKIP CALCULATIONS

            IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                 .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL &
#ifdef WRF_MARS
                  .and.qdust3d(k).lt.qsmall &
                  .and.qcore3d(k).lt.qsmall &
#endif
                 ) GOTO 500

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE INSTANTANEOUS PROCESSES

! ADD MELTING OF CLOUD ICE TO FORM RAIN
#ifndef WRF_MARS
        IF (QI3D(K).GE.QSMALL.AND.T3D(K).GE.273.15) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rain on Mars", k
#else
           QR3D(K) = QR3D(K)+QI3D(K)
           T3D(K) = T3D(K)-QI3D(K)*XLF(K)/CPM(K)
           QI3D(K) = 0.
           NR3D(K) = NR3D(K)+NI3D(K)
           NI3D(K) = 0.
#endif
        END IF
#endif
! ****SENSITIVITY - NO ICE
        IF (ILIQ.EQ.1) GOTO 778

#ifdef WRF_MARS
!not true
!        WRITE(*,*) "micro: ", __LINE__, "shouldn't ever run without ice on Mars", k
#else
! HOMOGENEOUS FREEZING OF CLOUD WATER

        IF (T3D(K).LE.233.15.AND.QC3D(K).GE.QSMALL) THEN
           QI3D(K)=QI3D(K)+QC3D(K)
           T3D(K)=T3D(K)+QC3D(K)*XLF(K)/CPM(K)
           QC3D(K)=0.
           NI3D(K)=NI3D(K)+NC3D(K)
           NC3D(K)=0.
        END IF

! HOMOGENEOUS FREEZING OF RAIN

        IF (IGRAUP.EQ.0) THEN

        IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
           QG3D(K) = QG3D(K)+QR3D(K)
           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
           QR3D(K) = 0.
           NG3D(K) = NG3D(K)+ NR3D(K)
           NR3D(K) = 0.
        END IF

        ELSE IF (IGRAUP.EQ.1) THEN

        IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
           QNI3D(K) = QNI3D(K)+QR3D(K)
           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
           QR3D(K) = 0.
           NS3D(K) = NS3D(K)+NR3D(K)
           NR3D(K) = 0.
        END IF

        END IF

#endif
 778    CONTINUE

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

#ifdef WRF_MARS
      NDUST3D(K) = MAX(0.,NDUST3D(K))
#endif
!......................................................................
! CLOUD ICE
#ifdef WRF_MARS
    if(.not.mars_heterogeneous) then !inuc.ne.6) then !homo
        ice_rho(k) = rhoi
        qicetot3d(k) = qi3d(k)
        if(qi3d(k).ge.qsmall) then 
            call constrain_mars_waterice_homogeneous(qi3d(k), ni3d(k), n0i(k), lami(k))
        endif
    else
        qicetot3d(k) = qcore3d(k) + qi3d(k)
        if(qicetot3d(K).GE.QSMALL) then
            call mean_ice_density(qcore3d(k), qi3d(k), ice_rho(k))
#ifdef MARS_MORR_DEBUG
            CALL_DEBUG("call constrain",__LINE__)
#endif
            call constrain_mars_waterice_heterogeneous(qi3d(k), ni3d(k), qcore3d(k), ice_rho(k), n0i(k), lami(k),thisi, thisj, k)
        endif
    endif
#else
      IF (QI3D(K).GE.QSMALL) THEN
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(K).LT.LAMMINI) THEN

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      ELSE IF (LAMI(K).GT.LAMMAXI) THEN
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      END IF
      END IF
#endif

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rain on Mars", k
#else
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
#endif
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have cloud on Mars", k
#else
         DUM = PRES(K)/(287.15*T3D(K))
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

#endif
      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "should we have snow on Mars", k
#else
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1
      NS3D(K) = N0S(K)/LAMS(K)
      END IF

#endif
      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have graupel on Mars", k
#else
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF

#endif
      END IF

#ifdef WRF_MARS
!DUST
    call constrain_mars_distribution(qdust3d(K), ndust3d(K), n0d(K), lamd(k))
#endif
 500  CONTINUE

! CALCULATE EFFECTIVE RADIUS
        
#ifdef WRF_MARS
    if(.not.mars_heterogeneous) then !inuc.ne.6) then !homo
        IF (QI3D(K).GE.QSMALL) THEN
      !  write(*,*) k, qsmall, qi3d(k), effi(k), lami(k), (4*pi*1000.*ni3d(k)/qi3d(k))**(1./3), 2e6/(4*pi*1000.*ni3d(k)/qi3d(k))**(1./3)
            EFFI(K) = ( ice_cons2 / LAMI(K) ) * 1.E6
        ELSE
            EFFI(K) = (ice_cons2 / LAMMINI ) * 1.e6
        endif
    else
        qicetot3d(k) = qcore3d(k) + qi3d(k)
        !TODO This means dust and ice have the same density
        IF (qicetot3d(K).GE.QSMALL) THEN
            EFFI(K) = ( ice_cons2 / LAMI(K) ) * 1.E6
        ELSE
            EFFI(K) = (ice_cons2 / LAMMINI ) * 1.e6
        endif
        
    endif            
#else
      IF (QI3D(K).GE.QSMALL) THEN
         EFFI(K) = 3./LAMI(K)/2.*1.E6
      ELSE
         EFFI(K) = 25.
      END IF
#endif

      IF (QNI3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have snow on Mars", k
#else
         EFFS(K) = 3./LAMS(K)/2.*1.E6
#endif
      ELSE
         EFFS(K) = 25.
      END IF

      IF (QR3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have rain on Mars", k
#else
         EFFR(K) = 3./LAMR(K)/2.*1.E6
#endif
      ELSE
         EFFR(K) = 25.
      END IF

      IF (QC3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have cloud water on Mars", k
#else
      EFFC(K) = GAMMA(PGAM(K)+4.)/                        &
             GAMMA(PGAM(K)+3.)/LAMC(K)/2.*1.E6
#endif
      ELSE
      EFFC(K) = 25.
      END IF

      IF (QG3D(K).GE.QSMALL) THEN
#ifdef WRF_MARS
        WRITE(*,*) "micro: ", __LINE__, "shouldn't have graupel on Mars", k
#else
         EFFG(K) = 3./LAMG(K)/2.*1.E6
#endif
      ELSE
         EFFG(K) = 25.
      END IF
#ifdef WRF_MARS
    if(DISTRIBUTION_MARS_DUST == 0) then !GAMMA
      IF (QDUST3D(K).ge.qsmall) THEN
         EFFD(K) = DUST_CONS2 / LAMD(K)
      ELSE
         EFFD(K) = DUST_CONS2 / LAMMAXD !TODO, I don't know what this value should be.
      END IF
    else if(DISTRIBUTION_MARS_DUST == 1) then !lognormal
      IF (QDUST3D(K).ge.qsmall) THEN
         EFFD(K) = 0.5 * exp(lamd(K) + 2.5 * sigma_squared)
      ELSE
         EFFD(K) = 25e-6 !TODO, I don't know what this value should be. NOTE that we multiple by 1e6 later
      END IF
    endif
    effd(k) = effd(k)*1e6 !convert metre to micron
#endif
! HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
! TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
! OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
!          NI3D(K) = MIN(NI3D(K),10.E6/RHO(K))
! HM, 12/28/12, LOWER MAXIMUM ICE CONCENTRATION TO ADDRESS PROBLEM
! OF EXCESSIVE AND PERSISTENT ANVIL
! NOTE: THIS MAY CHANGE/REDUCE SENSITIVITY TO AEROSOL/CCN CONCENTRATION
#ifdef WRF_MARS
!not sure what limit          NI3D(K) = MIN(NI3D(K),0.3E6/RHO(K))
#else
          NI3D(K) = MIN(NI3D(K),0.3E6/RHO(K))
#endif
! ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
          IF (iinum.EQ.0.AND.IACT.EQ.2) THEN
          NC3D(K) = MIN(NC3D(K),(NANEW1+NANEW2)/RHO(K))
          END IF
! SWITCH FOR CONSTANT DROPLET NUMBER
          IF (iinum.EQ.1) THEN 
! CHANGE NDCNST FROM CM-3 TO KG-1
             NC3D(K) = NDCNST*1.E6/RHO(K)
          END IF

      END DO !!! K LOOP

 400         CONTINUE

    !mass conservation
#ifdef WRF_MARS
    !fill in the negatives in water vapor and capture them in the rain
!    where(qv3d < 0)
!      qv3d = 0.
!    end where

    water_before = 0
    water_after = 0
    dust_before = 0
    dust_after = 0
    core_before = 0
    core_after = 0
    do k=kte, kts,-1
        iprint=.false.
          if(jprint.and.k.eq.KVALUE) iprint=.true.
        water_before = water_before + rho(k)*dzq(k) * (old_qv3d(k) + old_qi3d(k))
        water_after = water_after + rho(k)*dzq(k) * (qv3d(k) + qi3d(k))
        dust_before = dust_before + rho(k)*dzq(k) * (old_dust(k))
        dust_after = dust_after + rho(k)*dzq(k) * (qdust3d(k))
        core_before = core_before + rho(k)*dzq(k) * (old_core(k))
        core_after = core_after + rho(k)*dzq(k) * (qcore3d(k))
    end do
    
    snowrt = snowrt + PRECRT_PHOTOLYSIS
    precrt = precrt + PRECRT_PHOTOLYSIS
    
    water_after = water_after + snowrt
    dust_after = dust_after + PRECDUSTRT
    core_after = core_after + PRECCORERT
#ifdef MARS_MORR_DEBUG
    CALL_DEBUG("water_before",water_before)
    CALL_DEBUG("water_after",water_after)
    CALL_DEBUG("water_delta",water_after-water_before)
    CALL_DEBUG("snowrt",snowrt)
    CALL_DEBUG("snowrt_err",PRECQVRT_ERR)
    CALL_DEBUG("tottenW",sum(rho*dzq*(QI3DTEN+QV3DTEN))*dt)
    CALL_DEBUG("dust_before",dust_before+core_before)
    CALL_DEBUG("dust_after",dust_after+core_after)
    CALL_DEBUG("dust_delta",(dust_after+core_after)-(dust_before+core_before))
    CALL_DEBUG("precd",precdustrt)
    CALL_DEBUG("precderr",precdustrt_err)
    CALL_DEBUG("tottenD",sum(rho*dzq*(QDUST3DTEN+QCORE3DTEN))*dt)
#endif
    
    if(water_after .ne. water_before) then 
        snowrt = snowrt + (water_before-water_after)
        precrt = precrt + (water_before-water_after)

!alternate method
!    frac=0.
!    if (water_after .gt. water_before) then 
!        !too much water, remove some
!        !take as much water from the lowest layer need to remove half the deficit
!        if(qv3d(kts) .gt. 0) then 
!            frac = max(min(0.5, rho(k)*dzq(k)*qv3d(kts)/(water_after-water_before)),0.)
!        endif
!    else if (water_after .lt. water_before) then
!        !not enough water, add some back in
!        !put in enough water to equalize 0.5 of the waters, but don't oversaturate
!        frac = max(min(0.5,rho(kts)*dzq(kts)*(qvi(kts)-qv3d(kts))/(water_before-water_after)),0.)
!    endif
!        qv3d(kts) = qv3d(kts) + frac*(water_before-water_after)/(rho(kts)*dzq(kts))
!        precrt = precrt + (1-frac)*(water_before-water_after)
!        snowrt = snowrt + (1-frac)*(water_before-water_after)    
    endif
    if(dust_before+core_before .ne. dust_after+core_after) then 
        precdustrt = precdustrt + ((dust_before+core_before) - (dust_after+core_after))
    endif
#endif

! ALL DONE !!!!!!!!!!!
      RETURN
      END SUBROUTINE MORR_TWO_MOMENT_MICRO

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL FUNCTION POLYSVP (T,TYPE)

!-------------------------------------------

!  COMPUTE SATURATION VAPOR PRESSURE

!  POLYSVP RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

      IMPLICIT NONE

      REAL DUM
      REAL T
      INTEGER TYPE
! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i 
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
	6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/	

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8 

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

! ICE

      IF (TYPE.EQ.1) THEN

!         POLYSVP = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
!          LOG10(273.16/T)+0.876793*(1.-T/273.16)+						&
!          LOG10(6.1071))*100.


      dt = max(-80.,t-273.16)
      polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt))))))) 
      polysvp = polysvp*100.

      END IF

! LIQUID

      IF (TYPE.EQ.0) THEN

       dt = max(-80.,t-273.16)
       polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp = polysvp*100.

!         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
!             5.02808*LOG10(373.16/T)-									&
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+				&
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+				&
!             LOG10(1013.246))*100.

         END IF


      END FUNCTION POLYSVP

!------------------------------------------------------------------------------

      REAL FUNCTION GAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      implicit none
      INTEGER I,N
      LOGICAL PARITY
      REAL                                                          &
          CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      REAL, DIMENSION(7) :: C
      REAL, DIMENSION(8) :: P
      REAL, DIMENSION(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,				   &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,	   &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        END DO
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO I=1,N
            RES=RES*Y
            Y=Y+ONE
          END DO
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
            SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+xxx
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END FUNCTION GAMMA


      REAL FUNCTION DERF1(X)
      IMPLICIT NONE
      REAL X
      REAL, DIMENSION(0 : 64) :: A, B
      REAL W,T,Y
      INTEGER K,I
      DATA A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,	                          &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,							  &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0, 							  &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,							  &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      DATA (B(I), I = 0, 12) /                                  &
         -0.00000000029734388465E0,  0.00000000269776334046E0, 	&
         -0.00000000640788827665E0, -0.00000001667820132100E0,  &
         -0.00000021854388148686E0,  0.00000266246030457984E0, 	&
          0.00001612722157047886E0, -0.00025616361025506629E0, 	&
          0.00015380842432375365E0,  0.00815533022524927908E0, 	&
         -0.01402283663896319337E0, -0.19746892495383021487E0,  &
          0.71511720328842845913E0 /
      DATA (B(I), I = 13, 25) /                                 &
         -0.00000000001951073787E0, -0.00000000032302692214E0,  &
          0.00000000522461866919E0,  0.00000000342940918551E0, 	&
         -0.00000035772874310272E0,  0.00000019999935792654E0, 	&
          0.00002687044575042908E0, -0.00011843240273775776E0, 	&
         -0.00080991728956032271E0,  0.00661062970502241174E0, 	&
          0.00909530922354827295E0, -0.20160072778491013140E0, 	&
          0.51169696718727644908E0 /
      DATA (B(I), I = 26, 38) /                                 &
         0.00000000003147682272E0, -0.00000000048465972408E0,   &
         0.00000000063675740242E0,  0.00000003377623323271E0, 	&
        -0.00000015451139637086E0, -0.00000203340624738438E0, 	&
         0.00001947204525295057E0,  0.00002854147231653228E0, 	&
        -0.00101565063152200272E0,  0.00271187003520095655E0, 	&
         0.02328095035422810727E0, -0.16725021123116877197E0, 	&
         0.32490054966649436974E0 /
      DATA (B(I), I = 39, 51) /                                 &
         0.00000000002319363370E0, -0.00000000006303206648E0,   &
        -0.00000000264888267434E0,  0.00000002050708040581E0, 	&
         0.00000011371857327578E0, -0.00000211211337219663E0, 	&
         0.00000368797328322935E0,  0.00009823686253424796E0, 	&
        -0.00065860243990455368E0, -0.00075285814895230877E0, 	&
         0.02585434424202960464E0, -0.11637092784486193258E0, 	&
         0.18267336775296612024E0 /
      DATA (B(I), I = 52, 64) /                                 &
        -0.00000000000367789363E0,  0.00000000020876046746E0, 	&
        -0.00000000193319027226E0, -0.00000000435953392472E0, 	&
         0.00000018006992266137E0, -0.00000078441223763969E0, 	&
        -0.00000675407647949153E0,  0.00008428418334440096E0, 	&
        -0.00017604388937031815E0, -0.00239729611435071610E0, 	&
         0.02064129023876022970E0, -0.06905562880005864105E0,   &
         0.09084526782065478489E0 /
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T + 	  &
              A(K + 11)) * T + A(K + 12)) * W
      ELSE IF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + 	  &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + 	  &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + 	  &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF (X .LT. 0) Y = -Y
      DERF1 = Y
      END FUNCTION DERF1

!+---+-----------------------------------------------------------------+
#ifdef WRF_MARS

subroutine constrain_mars_waterice(qi3d, ni3d, qcore3d, reff, &
                          ids,ide, jds,jde, kds,kde,               &     
                          ims,ime, jms,jme, kms,kme,               &
                          its,ite, jts,jte, kts,kte                &
                          )
    implicit none
    integer, intent(in) :: ids,ide, jds,jde, kds,kde,               &     
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte            

    real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: qi3d
    real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: ni3d
    real, dimension(ims:ime, kms:kme, jms:jme), intent(in),optional :: qcore3d
    real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: reff

    integer :: i,j,k
    real :: n0i, lami, ice_rho, qtot
    
    if(.not.mars_heterogeneous) then
      !homogeneous
      do i=its,ite
        do k=kts,kte
          do j=jts,jte
            if(qi3d(i,k,j) .ge. qsmall) then
              call constrain_mars_waterice_homogeneous(qi3d(i,k,j), ni3d(i,k,j), n0i, lami)
              reff(i,k,j) = ( ice_cons2 / lami ) * 1.E6
            else
              ni3d(i,k,j) = 0.0
              reff(i,k,j) = ( ice_cons2 / lammini ) * 1.E6
            endif
          enddo
        enddo
      enddo
    else 
      !heterogeneous
      do i=its,ite
        do k=kts,kte
          do j=jts,jte
            qtot = qi3d(i,k,j) + qcore3d(i,k,j)
            if(qtot .ge. qsmall) then
              call mean_ice_density(qcore3d(i,k,j), qi3d(i,k,j), ice_rho)
              call constrain_mars_waterice_heterogeneous(qi3d(i,k,j), ni3d(i,k,j), qcore3d(i,k,j), ice_rho, n0i, lami,i, j, k)
              reff(i,k,j) = ( ice_cons2 / lami ) * 1.E6
            else
              ni3d(i,k,j) = 0.0
              reff(i,k,j) = ( ice_cons2 / lammini ) * 1.E6
            endif
          enddo
        enddo
      enddo
    endif
end subroutine constrain_mars_waterice

subroutine constrain_mars_dust(qdust3d, ndust3d,reff, &
                          ids,ide, jds,jde, kds,kde,               &     
                          ims,ime, jms,jme, kms,kme,               &
                          its,ite, jts,jte, kts,kte                &
                          )
    implicit none
    integer, intent(in) :: ids,ide, jds,jde, kds,kde,               &     
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte            

    real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: qdust3d
    real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: ndust3d
    real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: reff

    integer :: i,j,k
    real :: n0d, lamd, ice_rho, qtot
    
    do i=its,ite
      do k=kts,kte
        do j=jts,jte
          if(qdust3d(i,k,j) .ge. qsmall .and. ndust3d(i,k,j) .gt. 0.) then
            call constrain_mars_distribution(qdust3d(i,k,j), ndust3d(i,k,j), n0d, lamd)
            reff(i,k,j) = ( dust_cons2 / lamd ) * 1.E6
          else
            ndust3d(i,k,j) = 0.0
            reff(i,k,j) = ( dust_cons2 / lammind ) * 1.E6
          endif
        enddo
      enddo
    enddo
    
end subroutine constrain_mars_dust


subroutine constrain_mars_waterice_homogeneous(qi3d, ni3d, n0i, lami)
    real, intent(in) :: qi3d
    real, intent(inout) :: ni3d
    real, intent(out) :: n0i, lami
    real :: reff
    n0i=0.
    if (qi3d .ge. qsmall) then
        lami = (ICE_CONS12 * ni3d/qi3d)**(1./DI)
        
        if(lami .lt. lammini) then
            lami = lammini
            !n0i = qi3d * ice_cons11min !lami**(1+DI+mu) * qi3d*ICE_CONS11
            !ni3d = n0i*ICE_CONS10min ! *exp(-(1+ice_mu)*log(lami))
            ni3d = qi3d * ice_consmin_Q_N
        else if(lami .gt. lammaxi) then
            lami = lammaxi
            !n0i = qi3d * ice_cons11max !lami**(1+DI+mu) * qi3d*ICE_CONS11
            !ni3d = n0i*ICE_CONS10max !*exp(-(1+ice_mu)*log(lami))
            ni3d = qi3d * ice_consmax_Q_N
        else
            !n0i = ni3d*exp((1+ice_mu) * log(lami)) / ICE_CONS10
        endif
    else
      lami=lammini
    ENDIF
    
!    write(*,*) "J: ", qi3d, lami, n0i, ni3d, ice_cons12, ice_cons11min, ice_cons10min, ice_cons11max, ice_cons10max, ice_cons10

end subroutine constrain_mars_waterice_homogeneous

subroutine constrain_mars_waterice_heterogeneous(qi3d, ni3d, qcore3d, meanrho, n0i, lami,i,j,k)
    real, intent(in) :: qi3d, qcore3d, meanrho
    real, intent(inout) :: ni3d
    real, intent(out) ::  n0i, lami
    real :: reff
    integer, intent(IN) :: i,j,k
    real ::qnew3d, qicetot3d
    n0i=0.
    qicetot3d = 0.0
    if (qcore3d.ge.qsmall) qicetot3d=qicetot3d+qcore3d !annoying negatives from advection
    if (qi3d.ge.qsmall) qicetot3d=qicetot3d+qi3d !annoying negatives from advection
    qnew3d = qicetot3d * rhoi / meanrho
    !qnew3d is used to fold the meanrho value into the equation and replace rhoi in the constants defined above.
!    write(*,*) "preconstrain: ", qtot3d,qnew3d, qi3d, ni3d, qcore3d, n0i, meanrho, lami
!    write(*,*) qi3d, ni3d, qcore3d, meanrho, n0i, lami
#ifdef MARS_MORR_DEBUG
    CALL_DEBUG("constrain+++",0)
#endif
    if (qicetot3d .ge. qsmall) then
        lami = (ICE_CONS12 * ni3d/qnew3d)**(1./DI)
#ifdef MARS_MORR_DEBUG
        CALL_DEBUG("lami",lami)
#endif
        if(lami .lt. lammini) then
#ifdef MARS_MORR_DEBUG
    CALL_DEBUG("low lami",0)
    CALL_DEBUG("qicetot3d",qicetot3d)
    CALL_DEBUG("qcore3d",qcore3d)
    CALL_DEBUG("qi3d",qi3d)
    CALL_DEBUG("qnew3d",qnew3d)
    CALL_DEBUG("ni3d",ni3d)
    !CALL_DEBUG("rhoi",rhoi)
    CALL_DEBUG("meanrho",meanrho)
#endif
            lami = lammini
            !n0i = qnew3d * ice_cons11min !lami**(1+DI+mu) * qi3d*ICE_CONS11
            !ni3d = n0i*ICE_CONS10min ! *exp(-(1+ice_mu)*log(lami))
            ni3d = qnew3d * ice_consmin_Q_N
#ifdef MARS_MORR_DEBUG
    CALL_DEBUG("n0i",n0i)
    CALL_DEBUG("ni3d",ni3d)
    CALL_DEBUG("lami",lami)
    if(lami>0) CALL_DEBUG("reff",1e6*(dust_mu + 3)/(2*lami))
    CALL_DEBUG("constrain---",0)
#endif
        else if(lami .gt. lammaxi) then
#ifdef MARS_MORR_DEBUG
    CALL_DEBUG("high lami",0)
    CALL_DEBUG("qicetot3d",qicetot3d)
    CALL_DEBUG("qcore3d",qcore3d)
    CALL_DEBUG("qi3d",qi3d)
    CALL_DEBUG("qnew3d",qnew3d)
    CALL_DEBUG("ni3d",ni3d)
    !CALL_DEBUG("rhoi",rhoi)
    CALL_DEBUG("meanrho",meanrho)
#endif
            lami = lammaxi
            !n0i = qnew3d * ice_cons11max !lami**(1+DI+mu) * qi3d*ICE_CONS11
            !ni3d = n0i*ICE_CONS10max !*exp(-(1+ice_mu)*log(lami))
            ni3d = qnew3d * ice_consmax_Q_N
#ifdef MARS_MORR_DEBUG
    CALL_DEBUG("n0i",n0i)
    CALL_DEBUG("ni3d",ni3d)
    CALL_DEBUG("lami",lami)
    if(lami>0) CALL_DEBUG("reff",1e6*(dust_mu + 3)/(2*lami))
    CALL_DEBUG("constrain---",0)
#endif
        else
            !n0i = ni3d*exp((1+ice_mu) * log(lami)) / ICE_CONS10
        endif
    else
      lami=lammini
    ENDIF
    
end subroutine constrain_mars_waterice_heterogeneous

 subroutine constrain_mars_distribution(qdust3d, ndust3d, n0d, lamd)
    real, intent(in) :: qdust3d
    real, intent(inout) :: ndust3d, n0d, lamd
    real :: reff
    if (qdust3d.ge.qsmall) then
        
        if (DISTRIBUTION_MARS_DUST == 0) then !GAMMA
        lamd = (DUST_CONS1 * ndust3d/qdust3d)**(1/3.)
        
        if (lamd .lt. lammind) then
        !invert equations to get the N that corresponds to this lambda
            !if there are too many particles for the mass
            !lower the number of particles to make slightly larger particles
            lamd = LAMMIND
            !n0d = nfac_lammind4 * qdust3d
            !ndust3d = n0d * nfac_LAMMIND1
            ndust3d = qdust3d * LAMMIND_Q_N
        else if (lamd .gt. lammaxd) then
        !if there 'arent enough particles' for the mass, split the particles into particles
        !with the injection radius
            lamd = LAMDEFD
            !n0d = nfac_lamdefd4 * qdust3d
            !ndust3d = n0d *nfac_LAMDEFD1
            ndust3d = qdust3d * LAMDEFD_Q_N
        else
            !n0d = DUST_CONS_INV_G1 * exp((DUST_MU+1.)*log(lamd))*ndust3d
        endif
        
        else if(DISTRIBUTION_MARS_DUST == 1) then !lognormal

        !lamd is 'mu' here, something related to the effective radius
        !N = N0
        !q = N0 * (4.pi.rho/3) * exp(3mu + 9sigma**2 / 2)
        !mu = (1./3) * [ log(q/n * 3/(4*pi*rho)) - 9 * sigma_squared / 2]
        lamd = (1./3) * ( alog(qdust3d) - alog(ndust3d*dust_cons1) - 9.*sigma_squared / 2. )
        n0d = ndust3d
        reff = alog(0.5) + lamd + (5/2.) * sigma_squared
        
        if(reff .gt. ln_dust_reff_max) then
            reff = ln_dust_reff_default
            lamd = reff - alog(0.5) - (5/2.) * sigma_squared
            !n = q / (pi.rho/6) * exp(-(3*mu - 9./2 sigma**2))
            n0d = (qdust3d / dust_cons1) * exp(- (3*lamd + (9./2.)*sigma_squared))
            ndust3d = n0d
        else if (reff .lt. ln_dust_reff_min) then
            reff = ln_dust_reff_min
            lamd = reff - alog(0.5) - (5/2.) * sigma_squared
            n0d = (qdust3d / dust_cons1) * exp(- (3*lamd + (9./2.)*sigma_squared))
            ndust3d = n0d
        endif
        
        endif
    else if(DISTRIBUTION_MARS_DUST == 0 ) then !gamma
        lamd = lammaxd
    endif
    
 end subroutine constrain_mars_distribution
 
subroutine homogeneous_ice_physics(qv3d, qi3d, ni3d, qvi, dt,  MNUCCD, NNUCCD, PRD)
    implicit none
    real, intent(in) :: qv3d, qi3d, qvi, dt, ni3d
    real, intent(out):: MNUCCD, NNUCCD, PRD

    if (qv3d > qvi) then
    !condensation
        MNUCCD = (qv3d-qvi)/dt! complete in one timestep
        NNUCCD = MNUCCD * homogeneous_ice_number_cons1 !instant. creation of particles
    else if (qv3d + qi3d > qvi) THEN
    !sublimation to saturation
        PRD = (qv3d - qvi) / dt !negative rate to saturate atmosphere
!        NNUCCD = 0.0 !no destruction
    else if (qv3d + qi3d <= qvi) then
    !complete removal of cloud, not enough to saturate
        PRD = -qi3d / dt
!        NNUCCD = -ni3d
    endif
    
end subroutine homogeneous_ice_physics


real function fmx(m,x)
  implicit None
  real(kind(0.d0)),intent(in) :: x,m
  real(kind(0.d0)) :: p1, p2, p3, phi,y
  real :: ux
  if (x> 1e4) then 
    ux=1e4
  else
    ux=x
  endif

  if (x>5e3) then
    fmx = (2+m)*(1-m)**2/4.
  else
    !resolution goes bad above x=~10000
    phi = sqrt(1 - 2*m*x + x*x)
    y=(x-m)/phi
    p1 = ((1-m*x)/phi)**3
    p2 = x**3*(2 - 3*y + y**3)
    p3 = 3 * m * x**2 * (y-1)

    fmx = 0.5*( 1 + p1 + p2 + p3 )
  endif
  return
end function fmx

real(kind(0.d0)) function sigma(T)
  implicit none
  real, intent(in) :: T
  
  sigma=  1e-3*(141. - 0.15*T)
  return
end function sigma

subroutine condense_binned(T,rho,qv,qvi,number,diam,mass,dv_simple,latent,eis, lami,nice,dt,prd)
    implicit none  
    real, intent(in) :: T, rho, qv,qvi,dv_simple, eis, latent, nice,lami, dt
    real, intent(out) :: prd
    real(kind(0.d0)), intent(in), dimension(nbins) :: number, diam, mass
    real :: epsi, factor1, factor2, tf, kn, nuc_fd, nuc_fh, nuc_Fd_simple, nuc_ka, seq, sigvl, d_prd
    integer :: bin
    real(kind(0.d0)) :: diam_bin, number_bin
    !ka = (5.69 + 0.017T)*1e-5
    sigvl = sigma(T)
    !entirely the wrong units -> nuc_ka = (5.69 + 0.017*T3D(K)) * 1.e-5!
    nuc_ka = 4.184e-4 * (5.69-0.07*273.15+0.18*T) ! W/m/K (the 4184 converts from cal->J)
    prd = 0.
    !need to prep
    nuc_Fd_simple = (RV * RHOI * T) / (dv_simple * eis)
    nuc_Fh = ((latent * RHOI) / (T * nuc_ka)) * ((latent)/(RV * T) - 1)
    !write(*,*) "c1: ", nuc_Fd_simple, nuc_Fh
    do bin =1, nbins
      diam_bin = diam(bin) / lami
      number_bin = number(bin) * nice

      tf = 2.*dust_cons3 / rho            
      kn = tf/diam_bin
      factor1 = (1.33+0.71/kn)/(1+1./kn)
      factor2 = (1+factor1*kn) !crazy Jacobson scaling factor for molecular and viscous diffusion.
      !I think the 1.33 and 0.71 should be adjusted to the Knudsen used in sedim, but I don't 
      !have a clear derivation of the equation.
        
      !Fd = (R_universal . rho_i . T_infty) / (nu . M_w . esat(T_infty) )
      nuc_Fd = nuc_Fd_simple * factor2
     
      EPSI = (2*PI*number_bin*RHOI)/(nuc_Fd+nuc_Fh)
     
     !explicit scheme doesn't deplete water vapor during ice formation to slow down
     !here we assume radius change isn't significant over a timestep
      seq = exp((4*sigvl)/(RV*T*RHOI*diam_bin)) 
      !d_prd/dt = epsi * diam_bin * (qv/qvi - seq)
      !d_prd/dt = epsi * diam_bin * ((qv-dqi)/qvi - seq)
      !d_prd/(dt*epsi*diam_bin/qvi) = qv-dqi-seq*qvi
      !d_prd*((qvi/(dt*epsi*diam_bin)) + 1) = qv-seq*qvi
      !d_prd = (qv-seq*qvi)/ (1+(qvi/(dt*epsi*diam_bin)))
      d_prd = (qv-seq*qvi)/ (dt+(qvi/(epsi*diam_bin))) !extra 1./dt converts back to mean rate
      
      PRD = PRD + d_prd
      !write(*,*) "C: ", bin, diam_bin, d_prd, prd, qv
    end do
    !loop over bins done. Sad!

end subroutine condense_binned

subroutine heterogeneous_ice_nucleation_binned(P_water, T, qvqvsi, lamd, n0d, qdust3d, ndust3d, number, diam, mass, dt,&
                              core_number_rate, core_mass_rate, ice_mass_rate)
    implicit none
    real, intent(out) :: core_number_rate, core_mass_rate, ice_mass_rate
    real              :: d_core_number_rate, d_core_mass_rate, d_ice_mass_rate
    real, intent(in) :: P_water, T, lamd, n0d, ndust3d, qdust3d, qvqvsi, dt
    real(kind(0.d0)) :: Aterm, Bterm, Aprime
    real(kind(0.d0)) :: sigma, m, ag, dfkb, fm, sqrt_lam, gs
    real(kind(0.d0)), intent(in),dimension(nbins) :: number, diam, mass
    real(kind(0.d0)) :: number_bin, diam_bin,mass_bin,ct,cr, c1s, dfkbT, c1s_ss, c1s_ml
    real(kind(0.d0)) :: jrkb,jump_rate, g1, g2, Js, vol_h2o, nh2o, rad_h2o, zled
    real(kind(0.d0)) :: normm, normk, mk_mant, mk_exp, mk_exp_sqrt, mk_mant_sqrt, inv_mk_exp_sqrt
 
    integer :: bin
    sigma = 1e-3 * (141. - 0.15 * T)
    !vol_h2o = 18e-3 / 920.
    ag = (2 * sigma) / ( rv * T * RHOI * log(qvqvsi))
    vol_h2o= (mass_h2o / RHOI)
    rad_h2o = ((vol_h2o * 3)/(4*pi))**(1./3)
    gs = 4*pi*ag**3 / (3*vol_h2o)

    normm=mass_h2o*1d20
    normk=kboltz*1d20
    mk_mant = normk*normm !ick
    mk_exp = 1d-40
    mk_exp_sqrt = 1d-20
    inv_mk_exp_sqrt = 1d20
    mk_mant_sqrt = sqrt(mk_mant)

    dfkbT =  (4 * pi * ag * ag * sigma ) / (3.*kboltz*T)
    jump_rate = 1e13 !s-1 
!    jrkb =  1.38054e-10 ! jump_rate * boltzmann  
    g1 = 3e-20
    g2 = g1/10.
    if (nucleation_m > 1) then 
       m=0.954
    else if(nucleation_m > 0) then 
       m = nucleation_m
    else if (nucleation_m > -1) then 
       !trainer 
       m = max((-nucleation_m) - 6005 * exp(-0.065*T),0.6) !0.6 = m at 150K
    endif

    !m=0.954 !oops
!    !These factors don't depend on radius, calculate once
#ifdef MARS_MORR_DEBUG
!  CALL_DEBUG("P_water",P_water)
  CALL_DEBUG("kboltz",kboltz)
  CALL_DEBUG("ag",ag)
  CALL_DEBUG("T",T)
  CALL_DEBUG("dfkbT",dfkbT)
  CALL_DEBUG("G1",G1)
  CALL_DEBUG("G2",G2)
!  CALL_DEBUG("mass_h2o",mass_h2o)
  CALL_DEBUG("jump_rate",jump_rate)
#endif
    nh2o = max(P_water,0.) / (kboltz*T)
!    write(*,'(7(e18.12,","))') P_water, T, qvqvsi, lamd, n0d, qdust3d, ndust3d
    Aterm = 1.

    !ML
    c1s_ml = 2*rad_h2o/vol_h2o
    !write(*,*) "c1s_ml", c1s_ml
    !SS
    !write(*,*) P_water, pi , mass_h2o, kboltz, T, jump_rate, G1, G1/(kboltz*T)
    !write(*,*) P_water, jump_rate, mk_mat_sqrt, pi,T, exp(G1/(kboltz*T)), inv_mk_exp_sqrt
    !write(*,*) "c1s_ss", c1s_ss
    !SD
    !Bterm = exp((2*G1-G2)/(kboltz*T))
    !ct = sqrt(dfkbT/(3*pi))  * (1./gs) * (nh2o)**2 * (kboltz*T) * ag**2 * (1./(mass_h2o * jump_rate)) * Bterm
    !DVD
    !Bterm = 1.0
    !ct = sqrt(dfkbT/(3*pi))  * (1./gs) * (max(P_water,0) / sqrt(2*pi*mass_h2o*kboltz*T)) * pi*ag**2 * Bterm * c1s
    

    zled =  sqrt(dfkbT/(3*pi))  * (1./gs) * ag**2
    
    if(nucleation_surface_diffusion) then
      if (nucleation_monolayer) then 
        !SD+ML
        c1s_ml = 2*rad_h2o/vol_h2o
        Bterm = exp(-G2/(kboltz*T))
        ct = zled * (jump_rate) * Bterm * 2*pi*c1s_ml**2
        !write(*,*) "sd+ml", ct
      else
        !SD+SS 
        !c1s_ss = max(P_water,0.) * (1./jump_rate) * (1./(mk_mant_sqrt * sqrt(2*pi*T))) * exp(G1/(kboltz*T)) * inv_mk_exp_sqrt
        Bterm = exp((2*G1-G2)/(kboltz*T))
        !ct = sqrt(dfkbT/(3*pi))  * (1./gs) *2*pi* ag**2 * jump_rate * Bterm * c1s_ss**2
        ct = zled * (max(P_water,0.)**2) * (1./(mk_mant*T* jump_rate)) * Bterm * 1./mk_exp
        !write(*,*) "sd+ss", ct
      endif
    else 
      if (nucleation_monolayer) then 
        !DVD+ML
        c1s_ml = 2*rad_h2o/vol_h2o
        Bterm = 1.0
        ct = zled * max(P_water,0.) * (1./(mk_mant_sqrt*sqrt(2*pi*T))) * pi * c1s_ml * inv_mk_exp_sqrt 
        !write(*,*) "dvd+ml", ct
      else
        !DVD+SS
        Bterm = 1.0
        !ct = sqrt(dfkbT/(3*pi))  * (1./gs) * (max(P_water,0) / sqrt(2*pi*mass_h2o*kboltz*T)) * pi*ag**2 * Bterm * c1s_ss
        ct = zled * (max(P_water,0.)**2 / (2*pi*mk_mant*T)) * pi * (1./jump_rate) * exp(G1/(kboltz*T)) * 1./mk_exp
        !write(*,*) "dvd+ss", ct
      endif
    endif

      !old
      !Bterm = exp((2*G1-G2)/(kboltz*T))
      !ct = sqrt(dfkbT/(3*pi))  * (1./gs) * (nh2o)**2 * (kboltz*T) * ag**2 * (1./(mass_h2o * jump_rate)) * Bterm
      !write(*,*) "old", ct

!in an effort to move on. I'm keeping this for now. There are accuracy problems because of the large range in numbers involved.

    core_number_rate = 0.
    core_mass_rate = 0.
    ice_mass_rate = 0.
    
    do bin=1, nbins
      number_bin = number(bin) * ndust3d
      diam_bin = diam(bin) / lamd
      mass_bin = mass(bin) * qdust3d

      fm = fmx(m,diam_bin/2./ag)
      if (fm*dfkbT > 500.) then
        cr = 0.
      else
        cr = 1./sqrt(fm) * exp(-fm*dfkbT)
      endif

      Js = ct * cr
      
      !The 'traditional' interpretation of Js is that Js*CCN area = nucleation rate = core_number_rate.
      !this is wrong.
      !Aterm = (4*pi*ag**3)/(3)!volume of germ

      !explicit method overshoots in a single bin and overflows into other bins
      !rewrite as implicit to limit growth to physical limits
      !I think this will be the same as explicit cutoffs but a bit cleaner.

      !N_i+1 - N_i = - dt * Js' * N_i+1
      !N_i+1 = N_i / (1+dt*Js')
      !dN/dt = (N_i+1 - N_i)/dt = (-) N_i * Js' / (1+dt*Js')
      !J * pi r^2 is the total nucleation rate per particle / second
      !A is the area covered by a nucleation event
      !A/(pi r^2) is the fraction of surface area of particle covered per event
      !J * pi * r^2 * A / (pi * r^2) = JA is the fraction of particle covered per second = core_rate
      !core rate is J (nuc rate per second per area) * area covered (m2) by germ
      d_core_number_rate = dt * Js * (4/3.)*pi*ag**3/(2*h2o_molrad) !* pi * diam_bin**2 * dt
      d_core_number_rate = d_core_number_rate / (1+d_core_number_rate)  / dt !extra 1./dt to get rate.
      d_core_mass_rate = d_core_number_rate !same fractional rate of mass loss

      core_number_rate = core_number_rate + d_core_number_rate * number_bin
      core_mass_rate = core_mass_rate + d_core_mass_rate * mass_bin

      !no water vapor consumption anymore
      !technically, we consume a monomer layer * total surface area * density.
      d_ice_mass_rate = 0.0
      ice_mass_rate = ice_mass_rate + d_ice_mass_rate
#ifdef MARS_MORR_DEBUG
    !CALL_DEBUG("pi",pi)
    CALL_DEBUG("bin", bin)
    CALL_DEBUG("diam_bin", diam_bin)
    CALL_DEBUG("number_bin", number_bin)
    CALL_DEBUG("mass_bin", mass_bin)
    CALL_DEBUG("core_number_rate", core_number_rate)
    CALL_DEBUG("d_core_number_rate", d_core_number_rate)
    CALL_DEBUG("core_mass_rate", core_mass_rate)
    CALL_DEBUG("d_core_mass_rate", d_core_mass_rate)
    CALL_DEBUG("ice_mass_rate", ice_mass_rate)
    CALL_DEBUG("d_ice_mass_rate", d_ice_mass_rate)
#endif
    end do
    !end loop over bins, Sad!
    !write(*,*) "fin: ", core_number_rate, core_mass_rate, ndust3d, qdust3d
end subroutine heterogeneous_ice_nucleation_binned

subroutine mean_ice_density(qcore, qice, meanrho)
    real, intent(in) :: qcore, qice
    real, intent(out) :: meanrho
    
    meanrho = (qcore * DUST_CONS_RHO + qice * rhoi) / (qcore + qice)
    meanrho = min(DUST_CONS_RHO, max(rhoi, meanrho))
    
    return
end subroutine mean_ice_density


subroutine gamma_inc_inv_calc(mu,invals, outvals,nvals, name)
      implicit none
      real(kind(0.d0)), intent(in) :: mu
      integer, intent(in) :: nvals
      real(kind(0.d0)), dimension(nvals), intent(in) :: invals
      real(kind(0.d0)), dimension(nvals), intent(out) :: outvals
      character(len=*) :: name
      real(kind(0.d0)) :: pans, qans
      integer :: v,flag

      do v=1,nvals
        call gamma_inc_inv ( mu+1, outvals(v), 0., invals(v), 1-invals(v), flag )
        if (flag<0) then
          write(*,*) name,": Error in inverse incomplete gamma function: ", flag 
          write(*,*) invals(v), 1-invals(v), outvals(v)
          stop
        endif
        
      end do 
end subroutine gamma_inc_inv_calc

subroutine calc_weights(mu,edges, outvals,nvals)
      implicit none
      real(kind(0.d0)), intent(in) :: mu
      integer, intent(in) :: nvals
      real(kind(0.d0)), dimension(nvals+1), intent(in) :: edges
      real(kind(0.d0)), dimension(nvals), intent(out) :: outvals
      real(kind(0.d0)), dimension(nvals+1) :: temp,qtemp
      integer :: v,flag

      do v=1,nvals+1
        call gamma_inc ( mu+1, edges(v), temp(v),  qtemp(v), 0)
        if (temp(v) > 1.99) then
          write(*,*) ": Error in incomplete gamma function: ", v
          stop
        endif
      end do
      do v=1, nvals
        outvals(v) = (temp(v+1) - temp(v)) 
      end do
end subroutine calc_weights

real function average_mass_of_particle(dust_distribution,reff,veff)

        use module_model_constants, ONLY: pi2

        implicit none

! in-n-out:
        integer, intent(in) :: dust_distribution
        real, intent(in)    :: reff, veff

! locals:
        real :: mu_dust, sigma_squared, local_const
        real :: rhod

        rhod = DUST_CONS_RHO

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

#endif
 
END MODULE module_mp_morr_two_moment
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
