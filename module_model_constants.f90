!WRF:MODEL_LAYER:CONSTANTS
!

 MODULE module_model_constants

   !  1. Following are constants for use in defining real number bounds.

   !  A really small number.

#ifdef WRF_PLANET
   REAL    , PARAMETER :: epsilon_wrf     = 1.E-15
   ! The following definition overrides an intrinsic function, and is
   ! super dangerous.  Try to avoid its use if at all possible!
#else
   REAL    , PARAMETER :: epsilon         = 1.E-15
#endif

#if ( WRF_PLANET == 1 )
! some unit conversions
   REAL, PARAMETER :: g_to_kg   = 1.e3,  &   ! g -> kg conversion factor
                      atm_to_cb = 1.01325e2,  &  ! atm -> cb conversion factor
                      mb_to_cb  = 1.e-1,  &   ! mb -> cb conversion factor
                      cb_to_mb  = 1.e1   ! cb -> mb conversion factor
#endif

   !  2. Following is information related to the physical constants.

   !  These are the physical constants used within the model.

   ! Molar gas constant
   REAL    , PARAMETER :: bk           = 1.380649e-23          ! Boltzmann constant (J/K)          [16]
   REAL    , PARAMETER :: k_boltzmann  = bk
   REAL    , PARAMETER :: N_avogadro   = 6.02214076e23         ! Avogadro's number                 [16]
   REAL    , PARAMETER :: r_univ_gas   = bk*N_avogadro         ! "Universal" gas const [J/(mol.K)] [16]
   REAL    , PARAMETER :: R_universal  = (bk*1000.)*N_avogadro ! J kmol^-1 K^-1                    [16]

   !  References for constants used in the model:
   ! [1] Lodders, K., and B. Fegley, Jr. (Eds.) (1998), "The Planetary
   !     Scientist's Companion," Oxford University Press, New York
   ! [2] Cox, A. (Ed.) (2000), "Allen's Astrophysical Quantities," 4th ed.,
   !     Springer-Verlag, New York
   ! [3] Allison, M. (1997), Accurate analytic representations of solar time
   !     and season on Mars with applications to the Pathfinder/Surveyor
   !     missions, Geophys. Res. Lett., 24, 1967-1970.
   ! [4] Hourdin et al. (1995) via Allison (personal communication).
   !     Perihelion is given as Ls=278
   ! [5] Tokano (1999), for porous icy regolith
   !     Cf. Tokano (1999)'s 1.5e3 value for rock-ice mixture, which matches
   !     Earth values
   ! [6] Surkov et al. (1976), Density of surface rock on Venus from data
   !     obtained by the Venera 10 automatic interplanetary station,
   !     Russian:    Kosmicheskie Issledovaniia, 14, Sept.-Oct. 1976, 697-703.
   !     Translated: Cosmic Research, 14 (5), Mar. 1977, 612-618
   ! [7] http://solarsystem.nasa.gov/planets/profile.cfm?Object=Sat_Titan
   !     and
   !     http://solarsystem.nasa.gov/planets/profile.cfm?Object=Sat_Titan&Display=Facts&System=Metric
   ! [8] Buie et al. (2006), "Orbits and Photometry of Pluto's Satellites: 
   !     Charon, S/2005 P1, and S/2005 P2", The Astronomical Journal, 132 (1),
   !     290
   ! [9] Lellouch et al. (2009), "Pluto's lower atmosphere structure and 
   !     methane abundance from high-resolution spectroscopy and stellar
   !      occultations", Astronomy and Astrophysics, 495 (3) L17-L21
   ! [10] Lide, D. (Ed.) (1993), "CRC Handbook of Chemistry and Physics",
   !      73rd ed., CRC Press, Boca Raton
   ! [11] JPL Horizons software/web page
   ! [12] Dewar, J. (1904), "Physical Constants at Low Temperatures. (1)--The
   !      Densities of Solid Oxygen, Nitrogen, Hydrogen, etc.", Proceedings
   !      of the Royal Society of London, 73, pp. 251-261
   ! [13] Brown, G. N., Jr., and W. T. Ziegler (1979), "Vapor pressure and
   !      heats of sublimation of liquids and solids of interest in cryogenics
   !      below 1-atm pressure", Adv. Cryo. Eng. {Advances in Cryogenic
   !      Engineering}, 25, pp. 662-670
   ! [14] http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
   ! [15] Young, L. A. (2013), "Pluto's Seasons: New Predictions for New
   !      Horizons", Ap. J. Lett., 766 (L22)
   !      Hansen, C. J. and D. A. Paige (2013), "Pluto's Climate Modeled with
   !      New Observational Constraints," The Pluto System on the Eve of
   !      Exploration by New Horizons: Perspectives and Predictions, abstract
   ! [16] The NIST Reference on Constants, Units, and Uncertainty. NIST. 20 May 2019.
   !      2018 CODATA https://physics.nist.gov/cuu/Constants/index.html
   ! [17] Ismail, Ahmad Fauzi; Khulbe, Kailash; Matsuura, Takeshi, Gas Separation Membranes:
   !      Polymeric and Inorganic, Springer, 2015 ISBN 3319010956. [via Wikipedia]
   ! [18] Allison, M., & McEwen, M. (2000). A post-Pathfinder evaluation of areocentric solar 
   !      coordinates with improved timing recipes for Mars seasonal/diurnal climate studies. 
   !      Planetary and Space Science, 48(2-3), 215-235.
   ! [19] Li and Wang, Physical Review E, v. 68, 061206, 2003
   ! [20] Kopp, G.; Lean, J. L. (2011). "A new, lower value of total solar irradiance:
   !      Evidence and climate significance". Geophysical Research Letters. [via Wikipedia]
   ! [21] Houghton 'Physics of Atmospheres'
   ! [22] E.W. Lemmon, M.O. McLinden and D.G. Friend, "Thermophysical 
   !      Properties of Fluid Systems" in NIST Chemistry WebBook, NIST Standard 
   !      Reference Database Number 69, Eds. P.J. Linstrom and W.G. Mallard, 
   !      National Institute of Standards and Technology, Gaithersburg MD, 20899, 
   !      http://webbook.nist.gov, (retrieved November 17, 2014). 
   !      Specifically, the fluid properties data from 
   !      http://webbook.nist.gov/cgi/cbook.cgi?Name=nitrogen&Units=SI
   
   ! Universal constants (needed now rather than with others below)
   

   ! Molecular weights and molecular radii of gases:
   ! H2:
   REAL    , PARAMETER :: mw_h2      = 2.016          ! [1]   g/mol
   REAL    , PARAMETER :: h2_molrad  = (2.89e-10)/2.  ! [17]  m
   !  Ar:
   REAL    , PARAMETER :: mw_ar      = 39.948         ! [1]
   REAL    , PARAMETER :: ar_molrad  = (3.4e-10)/2.   ! [17]
   !  CH4:
   REAL    , PARAMETER :: mw_ch4     = 16.043         ! [1]
   REAL    , PARAMETER :: ch4_molrad = (3.8e-10)/2.   ! [17]
   !  CO:
   REAL    , PARAMETER :: mw_co      = 28.0102        ! [1]
   REAL    , PARAMETER :: co_molrad  = (3.76e-10)/2.  ! [17]
   !  CO2:
   REAL    , PARAMETER :: mw_co2     = 44.0096        ! [1]
   REAL    , PARAMETER :: co2_molrad = (3.3e-10)/2.   ! [17]
   !  H2O:
   REAL    , PARAMETER :: mw_h2o     = 18.0153        ! [1]
   REAL    , PARAMETER :: h2o_molrad = (2.65e-10)/2.  ! [17]
   !  N2:
   REAL    , PARAMETER :: mw_n2      = 28.0135        ! [1]
   REAL    , PARAMETER :: n2_molrad  = (3.64e-10)/2.  ! [17]
   !  O2:
   REAL    , PARAMETER :: mw_o2      = 31.999         ! [1]
   REAL    , PARAMETER :: o2_molrad  = (3.46e-10)/2.  ! [17]
   !  NO:
   REAL    , PARAMETER :: mw_no      = 30.006         ! [1]
   REAL    , PARAMETER :: no_molrad  = (3.17e-10)/2.  ! [17]
   !  He:
   REAL    , PARAMETER :: mw_he      = 4.002602       ! [1]
   REAL    , PARAMETER :: he_molrad  = (2.6e-10)/2.   ! [17]
   !  Ne:
   REAL    , PARAMETER :: mw_ne      = 20.179         ! [1]
   REAL    , PARAMETER :: ne_molrad  = (2.75e-10)/2.  ! [17]

! Stellar information:
   REAL    , PARAMETER :: solar_constant = 1361.5     ! W/m2 [20]
   REAL    , PARAMETER :: solar_emiss_t  = 5777.      ! K    [21]  solar black body emission temperature

! JM NOTE -- can we name this grav instead?

#ifdef WRF_PLANET
#if defined WRF_MARS
   !  Mars
   !  Ref. [1], Sec. 2.3
   REAL    , PARAMETER :: g = 3.727 ! acceleration due to gravity (m {s}^-2)
#elif defined WRF_PLUTO
   !  Pluto
   !  Ref. [8] for radius and mass, and g = GM/r^2
   REAL    , PARAMETER :: g = 0.658 ! acceleration due to gravity (m {s}^-2)
#elif defined WRF_TITAN
   !  Titan
   !  Ref. [1], Sec. 2.3
   REAL    , PARAMETER :: g = 1.354 ! acceleration due to gravity (m {s}^-2)
#elif defined WRF_VENUS
   !  Venus
   !  Ref. [1], Sec. 2.3
   REAL    , PARAMETER :: g = 8.870 ! acceleration due to gravity (m {s}^-2)
#elif defined WRF_JUPITER
   !  Jupiter
   REAL    , PARAMETER :: g = 22.88 ! acceleration due to gravity (m {s}^-2)
#elif defined WRF_TRITON
   !  Triton
   !  Wikipedia
   REAL    , PARAMETER :: g = 0.779 ! acceleration due to gravity (m {s}^-2)
#elif defined GEN_PLAN
   REAL    , PARAMETER :: g = 3.727 ! grav in MKS set this to whatever you want - GEN_PLAN setting
#elif defined WRF_EARTH
   ! Earth
   REAL    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)
#endif
#else
   REAL    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)
#endif

#if ( NMM_CORE == 1 )
   REAL    , PARAMETER :: r_d          = 287.04
   REAL    , PARAMETER :: cp           = 1004.6
#else
#ifdef WRF_PLANET
#if defined WRF_MARS
   !  Mars
   !  Molecular weight of atmosphere = 43.34; Ref. [1], Sec. 2.3
   !  Ref [2] for molecular diameter (CO2)
   REAL    , PARAMETER :: mw_air       = 43.34                ! kg/kmol
   REAL    , PARAMETER :: wtmol_air    = mw_air/1000.         ! kg/mol
   REAL    , PARAMETER :: r_d          = R_universal/mw_air   ! J/K/kg
   REAL    , PARAMETER :: cp           = 4.*r_d ! 8/2, linear triatomic
   REAL    , PARAMETER :: molecular_diameter = co2_molrad*2.  ! m
#elif defined WRF_PLUTO
   !  Pluto
   !  Molecular weight of atmosphere = 27.953553
   !  Ref [9], 99.5% N2 (28.0134), 0.5% CH4 (16.044)
   !  Ref [10] for cp (assuming heat capacities combine linearly)
   !  Ref [2] for molecular diameter (N2)
   REAL    , PARAMETER :: mw_air       = 27.953553            ! kg/kmol
   REAL    , PARAMETER :: wtmol_air    = mw_air/1000.         ! kg/mol
   REAL    , PARAMETER :: r_d          = R_universal/mw_air   ! J/K/kg
   REAL    , PARAMETER :: cp           = 1044.5735            ! J/K/kg
   REAL    , PARAMETER :: molecular_diameter = n2_molrad*2.   ! m
#elif defined WRF_TITAN
   !  Titan
   ! FROM http://pds-atmospheres.nmsu.edu/education_and_outreach/encyclopedia/gas_constant.htm:
   !  Ref [1] for mean molecular weight (Sec 9.2, Table 9.5)
   !  Ref [2] for molecular diameter (N2)
   REAL    , PARAMETER :: mw_air       = 28.6                 ! kg/kmol
   REAL    , PARAMETER :: wtmol_air    = mw_air/1000.         ! kg/mol
   REAL    , PARAMETER :: r_d          = 290.                 ! J/K/kg
   REAL    , PARAMETER :: cp           = 1044.0               ! J/K/kg
   REAL    , PARAMETER :: molecular_diameter = n2_molrad*2.   ! m
#elif defined WRF_TRITON
   !  Triton
   !  Molecular weight of atmosphere = 28.012203
   !  Ref [1], 99.99% N2 (28.0134), 0.01% CH4 (16.044)
   !  Ref [10] for cp (assuming heat capacities combine linearly)
   !  Ref [2] for molecular diameter (N2)
   REAL    , PARAMETER :: mw_air       = 28.012203            ! kg/kmol
   REAL    , PARAMETER :: wtmol_air    = mw_air/1000.         ! kg/mol
   REAL    , PARAMETER :: r_d          = R_universal/mw_air   ! J/K/kg
   REAL    , PARAMETER :: cp           = 1040.1151            ! J/K/kg
   REAL    , PARAMETER :: molecular_diameter = n2_molrad*2.   ! m
#elif defined WRF_VENUS
   !  Venus
   !  Molecular weight of atmosphere = 43.45; Ref. [1], Sec. 2.3
   !  Ref [2] for molecular diameter (CO2)
   REAL    , PARAMETER :: mw_air       = 43.45                ! kg/kmol
   REAL    , PARAMETER :: wtmol_air    = mw_air/1000.         ! kg/mol
   REAL    , PARAMETER :: r_d          = r_univ_gas/wtmol_air ! J/K/kg
   REAL    , PARAMETER :: cp           = 1116.8               ! J/K/kg
   REAL    , PARAMETER :: molecular_diameter = co2_molrad*2.  ! m
#elif defined WRF_JUPITER
   ! Jupiter
   REAL    , PARAMETER :: wtmol_air    = r_univ_gas/r_d       ! kg/mol
   REAL    , PARAMETER :: mw_air       = wtmol_air*1000.      ! kg/kmol
   REAL    , PARAMETER :: r_d          = 3770.                ! J/K/kg
   REAL    , PARAMETER :: cp           = 13000.               ! J/K/kg
   REAL    , PARAMETER :: molecular_diameter = h2_molrad*2.   ! m
#elif defined WRF_EARTH
   !  Ref [14] for mean molecular weight
   REAL    , PARAMETER :: mw_air       = 28.97                ! kg/kmol
   REAL    , PARAMETER :: wtmol_air    = mw_air/1000.         ! kg/mol
   REAL    , PARAMETER :: r_d          = 287.                 ! J/K/kg
   REAL    , PARAMETER :: cp           = 7.*r_d/2.            ! J/K/kg
   !  Ref [2] for molecular diameter (78% N2 [2.2] + 21% O2 [2.4] + 1% Ar [3.2])
   REAL    , PARAMETER :: molecular_diameter = 0.78*(n2_molrad*2.) + &
                                               0.21*(o2_molrad*2.) + &
                                               0.01*(ar_molrad*2.)
#endif
#else
#ifdef GEN_PLAN
   ! stand-in values - edit for what you want: (this is Mars)
   REAL    , PARAMETER :: mw_air       = 43.34                ! kg/kmol
   REAL    , PARAMETER :: wtmol_air    = mw_air/1000.         ! kg/mol
   REAL    , PARAMETER :: r_d          = R_universal/mw_air   ! J/K/kg
   REAL    , PARAMETER :: cp           = 4.*r_d ! 8/2, linear triatomic
   REAL    , PARAMETER :: molecular_diameter = co2_molrad*2.  ! m
#else
   REAL    , PARAMETER :: r_d          = 287.
   REAL    , PARAMETER :: cp           = 7.*r_d/2.
   !  Ref [14] for mean molecular weight
   REAL    , PARAMETER :: mw_air       = 28.97                ! kg/kmol
   REAL    , PARAMETER :: wtmol_air    = mw_air/1000.         ! kg/mol
   !  Ref [2] for molecular diameter (78% N2 [2.2] + 21% O2 [2.4] + 1% Ar [3.2])
   REAL    , PARAMETER :: molecular_diameter = 0.78*(n2_molrad*2.) + &
                                               0.21*(o2_molrad*2.) + &
                                               0.01*(ar_molrad*2.)
#endif
#endif
#endif

#if ( WRF_PLANET == 1 )
   REAL    , PARAMETER :: mwair_o_mwh2o = mw_air/mw_h2o
   REAL    , PARAMETER :: mwair_o_mwco2 = mw_air/mw_co2
   REAL    , PARAMETER :: mwair_o_mwch4 = mw_air/mw_ch4
#endif

   REAL    , PARAMETER :: r_v          = 461.6   ! J/K/kg, water
   REAL    , PARAMETER :: cv           = cp-r_d
   REAL    , PARAMETER :: cpv          = 4.*r_v
   REAL    , PARAMETER :: cvv          = cpv-r_v
   REAL    , PARAMETER :: cvpm         = -cv/cp
   REAL    , PARAMETER :: cliq         = 4190.
   REAL    , PARAMETER :: cice         = 2106.
   REAL    , PARAMETER :: psat         = 610.78
   REAL    , PARAMETER :: rcv          = r_d/cv
   REAL    , PARAMETER :: rcp          = r_d/cp
   REAL    , PARAMETER :: rovg         = r_d/g
   REAL    , PARAMETER :: c2           = cp * rcv
   real    , parameter :: mwdry        = 28.966 ! molecular weight of dry air (g/mole)

#ifdef WRF_PLANET
#if defined WRF_MARS
   !  Mars
   REAL    , PARAMETER :: p1000mb      = 610.    ! Pa
   REAL    , PARAMETER :: t0           = 300.    ! K
#elif defined WRF_PLUTO
   !  Pluto
   !  Note p1000mb is just a reference pressure ... no implication that
   !  this is the "standard" surface pressure
   !  Same for t0 ... just a reference temperature ... no implication that
   !  this is a surface or constant upper atmosphere temperature
   REAL    , PARAMETER :: p1000mb      = 2.      ! Pa 
   REAL    , PARAMETER :: t0           = 100.    ! K
#elif defined WRF_TITAN
   !  Titan
   REAL    , PARAMETER :: p1000mb      = 144000.
   REAL    , PARAMETER :: t0           = 120.
#elif defined WRF_TRITON
   !  Triton
   !  Note p1000mb is just a reference pressure ... no implication that
   !  this is the "standard" surface pressure
   !  Same for t0 ... just a reference temperature ... no implication that
   !  this is a surface or constant upper atmosphere temperature
   REAL    , PARAMETER :: p1000mb      = 2.      ! Pa 
   REAL    , PARAMETER :: t0           = 100.    ! K
#elif defined WRF_VENUS
   !  Venus
   REAL    , PARAMETER :: p1000mb      = 9.56e6
   REAL    , PARAMETER :: t0           = 800.
#elif defined WRF_JUPITER
   !  Jupiter
   REAL    , PARAMETER :: p1000mb      = 1.e7
   REAL    , PARAMETER :: t0           = 165.
#elif defined WRF_EARTH
   ! Earth
   REAL    , PARAMETER :: p1000mb      = 100000.
   REAL    , PARAMETER :: t0           = 300.
#endif
#else
   REAL    , PARAMETER :: p1000mb      = 100000.
   REAL    , PARAMETER :: t0           = 300.
#endif
   REAL    , PARAMETER :: p0           = p1000mb
   REAL    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL    , PARAMETER :: cvovcp       = 1./cpovcv
   REAL    , PARAMETER :: rvovrd       = r_v/r_d

#ifdef WRF_PLANET
   !  All planetary radii from Ref. [1], Sec. 2.3
#if defined WRF_MARS
   !  Mars
   REAL    , PARAMETER :: reradius     = 1./3389.92e03 ! 1/m
#elif defined WRF_PLUTO
   !  Pluto
   !  Ref. [8]
   REAL    , PARAMETER :: reradius     = 1./1153.e03   ! 1/m
#elif defined WRF_TITAN
   !  Titan
   REAL    , PARAMETER :: reradius     = 1./2575.0e03
#elif defined WRF_TRITON
   !  Triton
   !  Wikipedia
   REAL    , PARAMETER :: reradius     = 1./1353.4e03  ! 1/m
#elif defined WRF_VENUS
   !  Venus
   REAL    , PARAMETER :: reradius     = 1./6051.84e03
#elif defined WRF_JUPITER
   !  Jupiter
   REAL    , PARAMETER :: reradius     = 1./71492.0e03
#elif defined WRF_EARTH
   ! Earth
   REAL    , PARAMETER :: reradius     = 1./6370.0e03 
#endif
#else
#ifdef GEN_PLAN
   REAL    , PARAMETER :: reradius     = 1./3389.92e03 
#else
   REAL    , PARAMETER :: reradius     = 1./6370.0e03 
#endif
#endif

   REAL    , PARAMETER :: asselin      = .025
!   REAL    , PARAMETER :: asselin      = .0
   REAL    , PARAMETER :: cb           = 25.

#ifdef WRF_PLANET
   ! The following quantities appear to be latent heats of water
   ! XLS: latent heat of sublimation/deposition   (S <-> G)
   ! XLV: latent heat of evaporation/condensation (L <-> G)
   ! XLF: latent heat of freezing/melting         (S <-> L)
   ! The XLV0/1 and XLS0/1 may be constants in a temperature dependent
   ! form of the latent heats.  But see also our definiton of
   ! "h2o_lheat" below (latent heat of sublimation of water) which
   ! may be a duplication of these quantities (or at least a
   ! referenced, explained, and commented duplication).
   ! Units on all: J/kg
#endif
   REAL    , PARAMETER :: XLV0         = 3.15E6
   REAL    , PARAMETER :: XLV1         = 2370.
   REAL    , PARAMETER :: XLS0         = 2.905E6
   REAL    , PARAMETER :: XLS1         = 259.532
#ifdef WRF_PLANET
#ifdef WRF_TITAN
! these refer to methane for Titan
   REAL    , PARAMETER :: XLV_ch4      = 4.9E5
   REAL    , PARAMETER :: XLF_ch4      = 5.868E4
   REAL    , PARAMETER :: XLS_ch4      = XLV_ch4+XLF_ch4
   REAL    , PARAMETER :: RHOMETH      = 422.62 
!these refer to water:
   REAL    , PARAMETER :: XLV          = 2.5E6   ! latent heat vapourization
   REAL    , PARAMETER :: XLF          = 3.34E5  ! latent heat fusion
   REAL    , PARAMETER :: XLS          = XLV+XLF
#else
   REAL    , PARAMETER :: XLS          = 2.85E6
   REAL    , PARAMETER :: XLV          = 2.5E6
   REAL    , PARAMETER :: XLF          = 3.50E5
#endif
#endif
   REAL    , PARAMETER :: rhowater     = 1000.
   REAL    , PARAMETER :: rhosnow      = 100.
#ifdef WRF_PLANET
#if defined WRF_MARS
   !  Mars
   REAL    , PARAMETER :: rhoair0      = 0.0117      ! kg/m^3, P0/(R*273 K)
#elif defined WRF_PLUTO
   !  Pluto
   REAL    , PARAMETER :: rhoair0      = 1.2274e-5   ! kg/m^3, P0/(R*273 K)
#elif defined WRF_TITAN
   !  Titan
   REAL    , PARAMETER :: rhoair0      = 1.89
#elif defined WRF_TRITON
   !  Triton
   REAL    , PARAMETER :: rhoair0      = 1.2341e-5   ! kg/m^3, P0/(R*273 K)
#elif defined WRF_VENUS
   !  Venus
   REAL    , PARAMETER :: rhoair0      = 49960.55
#elif defined WRF_JUPITER
   !  Jupiter
   REAL    , PARAMETER :: rhoair0      = 16.0759
#elif defined WRF_EARTH
   ! Earth
   REAL    , PARAMETER :: rhoair0      = 1.28
#endif
#else
#ifdef GEN_PLAN
   REAL    , PARAMETER :: rhoair0      = 1.909
#else
   REAL    , PARAMETER :: rhoair0      = 1.28
#endif
#endif

!
   REAL    , PARAMETER :: n_ccn0       = 1.0E8
!
   REAL(KIND(0d0)), PARAMETER ::   pi_d = 3.141592653589793238462643383279502884197d0
   REAL           , PARAMETER ::   pi_s = REAL(pi_d)
   REAL           , PARAMETER ::piconst = pi_s
   REAL           , PARAMETER ::    pi2 = REAL(2.d0 * pi_d)
   REAL           , PARAMETER :: DEGRAD = REAL(pi_d/180.d0)
   REAL           , PARAMETER :: SQRTPI = REAL(DSQRT(pi_d))

   REAL(KIND(0.d0)), PARAMETER :: planetary_surface_area = 4.d0*pi_d*((1.d0/(dble(reradius)))**2.d0)

! vast block of dangerously double-defined DPD moved to end of code - dangerous as we were
! effectively defining the length of every planet year separately here and down below as
! planet_year

   REAL    , PARAMETER ::  SVP1=0.6112
   REAL    , PARAMETER ::  SVP2=17.67
   REAL    , PARAMETER ::  SVP3=29.65
   REAL    , PARAMETER ::  SVPT0=273.15
   REAL    , PARAMETER ::  EP_1=R_v/R_d-1.
   REAL    , PARAMETER ::  EP_2=R_d/R_v
   REAL    , PARAMETER ::  KARMAN=0.4
#ifdef WRF_PLANET
   REAL    , PARAMETER ::  SVPT0_METH=90.68   ! methane triple point temp
   REAL    , PARAMETER ::  PSAT_METH=10600.   ! methane triple point vapour pressure
   !  EOMEG uses the length of the sidereal (not solar) day
   !  EOMEG = 2*{pi}/sidereal day

!   saturation vapor pressure equation coefficients:
         ! For pure methane liquid:
   real, parameter :: vplA_CH4 = 3.901408,  &
                      vplB_CH4 = 437.54809,  &
                      vplC_CH4 = 1598.8512,  &
                      vplD_CH4 = 154567.02
      ! For methane/nitrogen mixture:
   real, parameter :: vp_N2_a0 = 0.76903175,  &
                      vp_N2_a1 = -8.1721842e-06,  &
                      vp_N2_a2 = 1.9970681e-08,  &
                      vp_N2_a3 = 1.7292750e-11

#if defined WRF_MARS
   !  Mars
   !  Sidereal day is 88642.6632 s; Ref [2], Sec. 12.2
   REAL    , PARAMETER ::  EOMEG=7.08822E-5     ! radians/s
#elif defined WRF_PLUTO
   !  Pluto
   !  Sidereal day is 551856.7 s; Ref [8]
   REAL    , PARAMETER ::  EOMEG=1.1385538E-5   ! radians/s
#elif defined WRF_TITAN
   !  Titan
   !  Sidereal day is 1377684.346752 s; Ref [2], Sec. 12.2
   REAL    , PARAMETER ::  EOMEG=4.5606857E-6
#elif defined WRF_TRITON
   !  Triton
   !  Sidereal day is 5 d 21 h 2 m 53 s = 507773 s ; Wikipedia
   REAL    , PARAMETER ::  EOMEG=1.2374E-5      ! radians/s
#elif defined WRF_VENUS
   !  Venus
   !  Sidereal day is 20996927.136 s; Ref [2], Sec. 12.2
   REAL    , PARAMETER ::  EOMEG=2.992430876e-7
#elif defined WRF_JUPITER
   !  Jupiter
   ! based on rotation period of 3.573e4 s
   REAL    , PARAMETER ::  EOMEG=1.7585181381e-4
#elif defined WRF_EARTH
   ! Earth
   REAL    , PARAMETER ::  EOMEG=7.2921E-5
#endif
#else
   REAL    , PARAMETER ::  EOMEG=7.2921E-5
#endif
   REAL    , PARAMETER ::  STBOLT=5.670374419E-8 ! [16]

#if ( WRF_PLANET == 1 )
   REAL    , PARAMETER ::  m_diffusivity = 1.18e-6        ! [m^2/s] [22] momentum diffusivity (kinematic viscosity
   REAL    , PARAMETER ::  t_diffusivity = 1.5e-6         ! [m^2/s] [22] thermal diffusivity
   REAL    , PARAMETER ::  ch4_v_diffusivity = 2.05e-6    ! [m^2/s] [22] methane vapor diffusivity
#endif

   REAL    , PARAMETER ::  prandtl = 1./3.0
                                         ! constants for w-damping option
   REAL    , PARAMETER ::  w_alpha = 0.3 ! strength m/s/s
   REAL    , PARAMETER ::  w_beta  = 1.0 ! activation cfl number

       REAL , PARAMETER ::  pq0=379.90516
       REAL , PARAMETER ::  epsq2=0.2
       REAL , PARAMETER ::  a2=17.2693882
       REAL , PARAMETER ::  a3=273.16
       REAL , PARAMETER ::  a4=35.86
       REAL , PARAMETER ::  epsq=1.e-12
       REAL , PARAMETER ::  p608=rvovrd-1.
!#if ( NMM_CORE == 1 )
       REAL , PARAMETER ::  climit=1.e-20
       REAL , PARAMETER ::  cm1=2937.4
       REAL , PARAMETER ::  cm2=4.9283
       REAL , PARAMETER ::  cm3=23.5518
!       REAL , PARAMETER ::  defc=8.0
!       REAL , PARAMETER ::  defm=32.0
       REAL , PARAMETER ::  defc=0.0
       REAL , PARAMETER ::  defm=99999.0
       REAL , PARAMETER ::  epsfc=1./1.05
       REAL , PARAMETER ::  epswet=0.0
       REAL , PARAMETER ::  fcdif=1./3.
#ifdef HWRF
       REAL , PARAMETER ::  fcm=0.0
#else
       REAL , PARAMETER ::  fcm=0.00003
#endif
       REAL , PARAMETER ::  gma=-r_d*(1.-rcp)*0.5
       REAL , PARAMETER ::  p400=40000.0
       REAL , PARAMETER ::  phitp=15000.0
       REAL , PARAMETER ::  plbtm=105000.0
       REAL , PARAMETER ::  plomd=64200.0
       REAL , PARAMETER ::  pmdhi=35000.0
       REAL , PARAMETER ::  q2ini=0.50
       REAL , PARAMETER ::  rfcp=0.25/cp
       REAL , PARAMETER ::  rhcrit_land=0.75
       REAL , PARAMETER ::  rhcrit_sea=0.80
       REAL , PARAMETER ::  rlag=14.8125
       REAL , PARAMETER ::  rlx=0.90
       REAL , PARAMETER ::  scq2=50.0
       REAL , PARAMETER ::  slopht=0.001
       REAL , PARAMETER ::  tlc=2.*0.703972477
       REAL , PARAMETER ::  wa=0.15
       REAL , PARAMETER ::  wght=0.35
       REAL , PARAMETER ::  wpc=0.075
       REAL , PARAMETER ::  z0land=0.10
#ifdef HWRF 
       REAL , PARAMETER ::  z0max=0.01
#else
       REAL , PARAMETER ::  z0max=0.008
#endif
       REAL , PARAMETER ::  z0sea=0.001
!#endif

   !  1. Universal constants
   !     (May only be needed for one planet though)

   !  Mass of unit atomic weight
   REAL    , PARAMETER :: atomicmassu  = 1.66053906660e-27 ! kg [16]
   REAL    , PARAMETER :: m_amu = atomicmassu

   REAL    , PARAMETER :: e_charge     = 1.602176634e-19   ! Charge of an electron [C] [16]
   REAL    , PARAMETER :: grav_const   = 6.6743e-11        ! Gravitation constant [N.m^2/kg^2] [16]

   !  Latent heat of sublimation, CO2
   !  Ref: <http://en.wikipedia.org/wiki/Latent_heat>
   REAL    , PARAMETER :: co2_lheat    = 5.713e5 ! J/kg

   !  Latent heat of sublimation, H2O
   !  Ref: <http://en.wikipedia.org/wiki/Ice_1h> 50911 J/mol
   !  Molecular weight of water: 18.0153 kg/kmol
   REAL    , PARAMETER :: h2o_lheat    = 2.826e6 ! J/kg
   !  Ref: <http://en.wikipedia.org/wiki/Latent_heat>
   !  but only for range -40 to 0 C
   !REAL    , PARAMETER :: h2o_lheat    = 2.8341e6 ! J/kg

   !  Latent heat of sublimation, N2
   !  Ref: [13]
   REAL    , PARAMETER :: n2_lheat     = 2.6e5 ! J/kg
   !  Value to use when comparing with older Triton models
   !REAL    , PARAMETER :: n2_lheat     = 2.5e5   ! J/kg
   ! This is universal, true for any planet
   REAL    , PARAMETER :: bar2pa       = 101325.
   REAL    , PARAMETER :: pa2bar       = 1./bar2pa

   !  Latent heat of sublimation, CH4
   !  Ref: [10] ... adding fusion and vaporization to approximate sublimation
   !  L_vapor = 8.19 kJ/mol, L_fusion = 0.94 kJ/mol, mol. wt. = 16.043 g/mol
   REAL    , PARAMETER :: ch4_lheat    = 5.691e5 ! J/kg

#ifdef WRF_PLANET
   !  Additional constants necessary for planetWRF
   !
!c-mm This following block of constants are used by the regolith code
   INTEGER , PARAMETER :: ndstp        = 14           !c-mm Number of sub-timesteps for regolith code
!   REAL    , PARAMETER :: grid_factor  = 1.80         !c-mm Grid-spacing coefficient
   REAL    , PARAMETER :: pore_size    = 5.e-6        !c-mm Pore size [m]
   REAL    , PARAMETER :: soil_surface_area = 1.7e4   !c-mm  Soil surface area (m2/kg)
   REAL    , PARAMETER :: tortuosity   = 5.0          !c-mm Soil tortuosity
   REAL    , PARAMETER :: heat_flow    = 0.030        !c-mm Areothermal heat flow
   REAL    , PARAMETER :: Aw           = 3.56e12      !c-mm Vapor pressure constant (Pa)
   REAL    , PARAMETER :: Bw           = 6141.7       !c-mm Vapor pressure constant (K)
!c-mm Replace these next two values which I got from the original adsorbate code with values pulled from Zent et al. (1993).
!c-mm **Note** There is a slight change in the ss_partition code to account for the change in sign of the eps term.
!c-mm   REAL    , PARAMETER :: eps          = 2573.9       !c-mm Temperature term in adsorbate isotherm
!c-mm   REAL    , PARAMETER :: Ko           = 1.57e-8      !c-mm Frequency term in adsorbate isotherm (1/Pa)
   REAL    , PARAMETER :: eps_ads      = -2679.8       !c-mm Temperature term in adsorbate isotherm
   REAL    , PARAMETER :: Ko_ads       = 2.043e-8      !c-mm Frequency term in adsorbate isotherm (1/Pa)
   REAL    , PARAMETER :: mass_co2     = 7.306e-26    !c-mm Mass of one CO2 molecule (kg)
   REAL    , PARAMETER :: mass_h2o     = 2.988e-26    !c-mm Mass of one H2O molecule (kg)
   REAL    , PARAMETER :: nu_ads       = 0.50         !c-mm Exponent in adsorbate isotherm
   REAL    , PARAMETER :: sbc          = STBOLT          ! Stefan-Boltzmann constant
   REAL    , PARAMETER :: diam_h2o     = h2o_molrad*2.     ! Collision diameter of H2O (m)
   REAL    , PARAMETER :: h_planck     = 6.62607015e-34    ! Planck's constant (m2kg/s)  [16]
   REAL    , PARAMETER :: light_speed  = 2.99792458e8      ! Speed of light in vacuum (m/s) [16]
   REAL    , PARAMETER :: magnetic_const = 1.e-7*pi2*2. ! Magnetic constant [N/A^2] 
   REAL    , PARAMETER :: electric_const = 1./(light_speed*light_speed*magnetic_const) ! Electric constant [F/m] (fn.2)
   REAL    , PARAMETER :: coulomb_const  = 1./(pi2*2.*electric_const)             ! Coulomb's (force) constant [m/F?]

   !  2. Planet-specific constants
   !     A. Clocks
   !        P2SI:             Ratio of the number of seconds in a planetary sol
   !                          (solar day, i.e., time between consecutive 
   !                           occurences of the sun directly overhead) divided
   !                          by the number of seconds in a terrestrial solar
   !                          day.
   !                          A useful defining equation is:
   !                           1/solar day = 1/sidereal day - 1/sidereal year
   !                          where negative values are used for retrograde
   !                          motion
   !        PLANET_YEAR:      Number of sols (solar days) in one sidereal
   !                          revolution.
   !                          Must be an integer, so round to nearest integer.
   !                          "Planetary leap days" must be handled specially
   !                          if desired (right now they are ignored)
   !     B: Orbital parameters
   !        SEMIMAJOR_AXIS:   Length of semi-major axis of elliptical planetary
   !                          orbit, in AU
   !        ECCENTRICITY:     Eccentricity of the elliptical planetary orbit,
   !                          unitless
   !        ZERO_DATE:        Time of perihelion passage, in number of sols
   !                          past occurence of instant of solar longitude = 0.
   !                          (Typically, we will assume solar longitude = 0
   !                           occurs at the instant of northern vernal
   !                           [spring] equinox.)
   !        EQUINOX_FRACTION: Fraction into the year (from *perihelion*) of the
   !                          occurence of northern spring equinox.  This
   !                          number is degenerate with the ZERO_DATE if the
   !                          reference point of solar longitude = 0 is taken
   !                          to be northern vernal equinox.  If that is true,
   !                          this value is derivable from ZERO_DATE (or vice
   !                          versa) and can be calculated by:
   !                          EQUINOX_FRACTION = 
   !                                    (PLANET_YEAR - ZERO_DATE)/PLANET_YEAR .
   !        OBLIQUITY:        Tilt of planetary rotation axis relative to
   !                          plane of planetary orbit, in degrees.
   !     C: Surface parameters
   !        Generally, these are only needed for a more developed surface
   !        and subsurface physical parameterization.  At present, this is
   !        primarily Mars.  However, if defined for one planet, all planets
   !        need to have a value defined, even if it's not used, for successful
   !        compilation
   !        RHO_GROUND:       Density of soil, kg m^-3
   !        Mars specific (NOW IN NAMELIST!!!!!!!!!!!!!!):
   !        CO2ICE_THRESHOLD: The amount of CO2 ice above which the albedo of
   !                          the surface should be set to the albedo of CO2
   !                          ice.  Empirically determined from test runs and
   !                          matches to observed quantities.  Units: kg m^-2
   !        H2OICE_THRESHOLD: The amount of H2O ice above which the albedo of
   !                          the surface should be set to the albedo of H2O
   !                          ice.  Empirically determined from test runs and
   !                          matches to observed quantities.  Units: kg m^-2
   !        N2ICE_THRESHOLD:  The amount of N2 ice above which the albedo and
   !                          emissivity of the surface should be set to the
   !                          albedo of N2 ice.  Units: kg m^-2
   !        CH4ICE_THRESHOLD: The amount of N2 ice above which the albedo and
   !                          emissivity of the surface should be set to the
   !                          albedo of N2 ice.  Units: kg m^-2
#if defined WRF_MARS
   !  Mars

   REAL    , PARAMETER :: P2SI         = 1.027491252 !  [18]
#if defined MARS24_TIMING
   INTEGER , PARAMETER :: PLANET_YEAR  = 100000
#else
   !  One Martian sidereal year = 686.929711 d = 668.550424 sols; Ref. [2]
   INTEGER , PARAMETER :: PLANET_YEAR  = 669
#endif
   REAL ,    PARAMETER :: SEMIMAJOR_AXIS   = 1.52366231         ! Ref. [2]
   REAL ,    PARAMETER :: ECCENTRICITY     = 0.09341233         ! Ref. [2]
   REAL ,    PARAMETER :: ZERO_DATE        = 488.7045           ! derived
   REAL ,    PARAMETER :: EQUINOX_FRACTION = ((PLANET_YEAR)-ZERO_DATE)/PLANET_YEAR
   REAL ,    PARAMETER :: OBLIQUITY        = 25.19              ! Ref. [2]

   !  Surface / regolith property constants
   REAL , PARAMETER :: CO2_MIXING_RATIO = 0.953
   REAL , PARAMETER :: N2_MIXING_RATIO = 0.026
#endif

#if ( WRF_PLANET == 1 )
!    Using values of A, B, and E reported in Table 1 of Li and Wang:
!    A = 1.25 (avg), B = 0.43 (avg), E = 0.95 (avg) for Cunningham correction:
   REAL, PARAMETER  :: a_part_sedim = 1.25  ! [19]
   REAL, PARAMETER  :: b_part_sedim = 0.43  ! [19]
   REAL, PARAMETER  :: e_part_sedim = 0.95  ! [19]
   REAL, PARAMETER  :: co2_molmass = mw_co2/1000.
#endif

#ifdef WRF_PLUTO
   !  Pluto

   !  One Plutonian solar day = 6d 9h 16m 57.78 s = 551817.78 s; Ref. [8]
   !  for length of sidereal day and year, and defining equation
   REAL ,    PARAMETER :: P2SI             = 6.3867799
   !  One Plutonian sidereal year = 90553.017 d = 14178.196 sols; Ref. [11]
   INTEGER , PARAMETER :: PLANET_YEAR      = 14178
   REAL ,    PARAMETER :: SEMIMAJOR_AXIS   = 40.09855777865089  ! Ref. [11]
   REAL ,    PARAMETER :: ECCENTRICITY     = 0.2566687011667356 ! Ref. [11]
   REAL ,    PARAMETER :: ZERO_DATE        = 99.4131          ! derived
   REAL ,    PARAMETER :: EQUINOX_FRACTION = 1.-(ZERO_DATE/REAL(PLANET_YEAR))
   REAL ,    PARAMETER :: OBLIQUITY        = 119.6              ! Refs. [8,11]
   !REAL ,    PARAMETER :: RHO_GROUND       = 2060.              ! planet mean
   !                                                             ! density
   !                                                            ! N2 ice:
   REAL ,    PARAMETER :: RHO_GROUND       = 1026.              ! Ref. [12]
   REAL ,    PARAMETER :: CO2_MIXING_RATIO = 0.                 ! Dummy, unused
#endif

#ifdef WRF_TITAN
   !  Titan

   !  One Titan solar day = 15.969131 d; Ref. [2], Sec. 12.2
   !     = 1./((1./15.94542068 d)-(1./10739.585 d))
   REAL ,    PARAMETER :: P2SI             = 15.969131
   !  One Titan sidereal year = 672.52 sols ; Ref. [7]
   INTEGER , PARAMETER :: PLANET_YEAR      = 673
   REAL ,    PARAMETER :: SEMIMAJOR_AXIS   = 9.5719             ! Ref. [1] (Saturn)
   REAL ,    PARAMETER :: ZERO_DATE        = 531.2              ! Calculated using ref. [4]
   REAL ,    PARAMETER :: EQUINOX_FRACTION=((PLANET_YEAR)-ZERO_DATE)/PLANET_YEAR
   REAL ,    PARAMETER :: OBLIQUITY        = 26.7               ! Ref. [1?]
   REAL ,    PARAMETER :: ECCENTRICITY     = 0.05415060         ! Ref. [2] (Saturn)
   REAL ,    PARAMETER :: RHO_GROUND       = 800.0              ! Ref. [5]
   REAL ,    PARAMETER :: CO2_MIXING_RATIO = 0.                 ! Dummy, unused
   REAL ,    PARAMETER :: N2_MIXING_RATIO  = 0.97               ! 
#endif

#ifdef WRF_TRITON
   !  Triton

   !  One Tritonian solar day = 5.8762802 d = 507710.61 s; Dick's notes
   REAL ,    PARAMETER :: P2SI             = 5.8762802
   !  One Neptunian sidereal year = 60190.03 d ; Wikipedia
   !  But Dick's notes show a calculation of 10243 sols
   INTEGER , PARAMETER :: PLANET_YEAR      = 10243
   REAL ,    PARAMETER :: SEMIMAJOR_AXIS   = 30.07090           ! Dick's notes
   REAL ,    PARAMETER :: ECCENTRICITY     = 0.00867797         ! Dick's notes
   REAL ,    PARAMETER :: EQUINOX_FRACTION = 0.023766395        ! Dick's notes
   REAL ,    PARAMETER :: ZERO_DATE        = PLANET_YEAR*(1.-EQUINOX_FRACTION) ! derived
   REAL ,    PARAMETER :: OBLIQUITY        = 129.608            ! Refs. [8,11]
   !REAL ,    PARAMETER :: RHO_GROUND       = 2059.              ! Dick's notes
   ! Value used in older Triton models
   REAL ,    PARAMETER :: RHO_GROUND       =  500.              ! Spencer and Moore, 1992
   REAL ,    PARAMETER :: CO2_MIXING_RATIO = 0.                 ! Dummy, unused
#endif

#ifdef WRF_JUPITER
   !  Jupiter

   !  One jupiter sidereal year = 10475.8 sols
   INTEGER , PARAMETER :: PLANET_YEAR      = 10476
   REAL ,    PARAMETER :: SEMIMAJOR_AXIS   = 5.204267
   REAL ,    PARAMETER :: ZERO_DATE        = 0.                 ! INCORRECT NEEDS TO BE CALCUATED
   REAL ,    PARAMETER :: EQUINOX_FRACTION=((PLANET_YEAR)-ZERO_DATE)/PLANET_YEAR
   REAL ,    PARAMETER :: OBLIQUITY        = 3.13
   REAL ,    PARAMETER :: ECCENTRICITY     = 0.048775
   REAL ,    PARAMETER :: RHO_GROUND       = 800.0              ! Ref. [5]
   REAL ,    PARAMETER :: CO2_MIXING_RATIO = 0.                 ! Dummy, unused
#endif

#ifdef WRF_VENUS
   !  Venus

   !  One Venus solar day = 116.70784 d = 10083557.36655 s; Ref. [2], Sec. 12.2
   !     = 1./((1./224.54163805 d) - (1./-243.01999 d))
   REAL , PARAMETER :: P2SI = 116.70784
   !  One Venus sidereal year = 224.69543 d = 1.9252814 sols; Ref. [2]
   INTEGER , PARAMETER :: PLANET_YEAR      = 2
   REAL ,    PARAMETER :: SEMIMAJOR_AXIS   = 0.72333199         ! Ref. [2]
   !  Date of Venus Vernal equinox = Feb 24 2000 19:30 UT
   !    <http://www.eso.org/public/outreach/eduoff/vt-2004/Background/Infol2/EIS-D3.html>
   !  Date of Venus Perihelion = Jul 13 2000
   !    <http://www.shodor.org/MASTER/galaxsee/curriculum/solar_system_lesson.html>
   !  Assuming Ls=0 at Vernal equinox, then:
   !    ZERO_DATE = (Jul 13, 2000 - Feb 24, 2000)/116.70784 d
   REAL ,    PARAMETER :: ZERO_DATE        = 1.1969
   REAL ,    PARAMETER :: EQUINOX_FRACTION=((PLANET_YEAR)-ZERO_DATE)/PLANET_YEAR
   REAL ,    PARAMETER :: OBLIQUITY        = 177.4              ! Ref. [1]
   REAL ,    PARAMETER :: ECCENTRICITY     = 0.00677323         ! Ref. [2]
   REAL ,    PARAMETER :: RHO_GROUND       = 2700.              ! Ref. [6]
   REAL ,    PARAMETER :: CO2_MIXING_RATIO = 0.                 ! Dummy, unused
#endif

#ifdef WRF_EARTH
   ! Earth, if run using the WRF_PLANET OPTION
   REAL, PARAMETER :: P2SI = 1.E0
   INTEGER , PARAMETER :: PLANET_YEAR      = 365
   REAL ,    PARAMETER :: OBLIQUITY        = 23.45              ! Ref [2]
   REAL ,    PARAMETER :: ECCENTRICITY     = 0.01671022         ! Ref [2]
   REAL ,    PARAMETER :: SEMIMAJOR_AXIS   = 1.0                ! by definition
   ! According to the US Naval Observatory,
   !    <http://aa.usno.navy.mil/data/docs/EarthSeasons.php>
   ! from 2000 to 2020, the average date of perihelion is Jan 3, 18:00, and
   ! the average data of the northern spring equinox is Mar 20, 14:13.
   ! This gives a equinox fraction of (Mar 20 14:13 - Jan 3 18:00)/365.2425
   ! Also, the zero date is (Jan 3 18:00 - Jan 1, 00:00)
   ! These values are only necessary when using Earth as a generalized planet
   ! and the assumption is that solar longitude of 0 starts on Jan 1 at midnight
   ! (and not at the instant of northern vernal equinox)
   REAL ,    PARAMETER :: ZERO_DATE        = 2.75
   REAL ,    PARAMETER :: EQUINOX_FRACTION = 0.21038724
   REAL ,    PARAMETER :: RHO_GROUND       = 2.e3               ! Dummy, unused
!   REAL ,    PARAMETER :: CO2ICE_THRESHOLD = 1.e30              ! Dummy, unused
!   REAL ,    PARAMETER :: H2OICE_THRESHOLD = 1.e30              ! Dummy, unused
   REAL ,    PARAMETER :: CO2_MIXING_RATIO = 0.                 ! Dummy, unused
   REAL ,    PARAMETER :: N2_MIXING_RATIO  = 0.78084            ! Unused. Ref. [1].
#endif

#else

   !  Earth

   !  The value for P2SI *must* be set to 1.0 for Earth
   REAL    , PARAMETER :: P2SI         = 1.0

   !  Orbital constants:

   INTEGER , PARAMETER :: PLANET_YEAR = 365
   REAL , PARAMETER :: OBLIQUITY = 23.5
   REAL , PARAMETER :: ECCENTRICITY = 0.014
   REAL , PARAMETER :: SEMIMAJOR_AXIS = 1.0 ! In AU
   ! Don't know the following values, so we'll fake them for now
   REAL , PARAMETER :: zero_date = 0.0   ! Time of perihelion passage
   !  Fraction into the year (from perhelion) of the
   !  occurrence of the Northern Spring Equinox
   REAL , PARAMETER :: EQUINOX_FRACTION= 0.0
   !  The value of CO2ICE_THRESHOLD, H2OICE_THRESHOLD, RHOGRND can be anything
   !  Just need values to get it to compile
!   REAL , PARAMETER :: CO2ICE_THRESHOLD = 1.e30
!   REAL , PARAMETER :: H2OICE_THRESHOLD  = 1.e30
!   REAL , PARAMETER :: rhogrnd           = 2.e3
   REAL ,    PARAMETER :: CO2_MIXING_RATIO = 0.                 ! Dummy, unused
   REAL ,    PARAMETER :: N2_MIXING_RATIO  = 0.78084            ! Unused. Ref. [1].

#endif

 REAL    , PARAMETER :: DPD       = 360./PLANET_YEAR
 REAL, PARAMETER :: SEMIMAJORAXIS = SEMIMAJOR_AXIS

 CONTAINS
   SUBROUTINE init_module_model_constants
   END SUBROUTINE init_module_model_constants
 END MODULE module_model_constants
