!WRF:DRIVER_LAYER:CONSTANTS
!
!  This MODULE contains all of the constants used in the model.  These
!  are separated by uage within the code.

MODULE module_driver_constants

   !  0. The following tells the rest of the model what data ordering we are
   !     using

   INTEGER , PARAMETER :: DATA_ORDER_XYZ = 1
   INTEGER , PARAMETER :: DATA_ORDER_YXZ = 2
   INTEGER , PARAMETER :: DATA_ORDER_ZXY = 3
   INTEGER , PARAMETER :: DATA_ORDER_ZYX = 4
   INTEGER , PARAMETER :: DATA_ORDER_XZY = 5
   INTEGER , PARAMETER :: DATA_ORDER_YZX = 6
   INTEGER , PARAMETER :: DATA_ORDER_XY = DATA_ORDER_XYZ
   INTEGER , PARAMETER :: DATA_ORDER_YX = DATA_ORDER_YXZ


#include <model_data_order.inc>

   !  1. Following are constants for use in defining maximal values for array
   !     definitions.  
   !

   !  The maximum number of levels in the model is how deeply the domains may
   !  be nested.

   INTEGER , PARAMETER :: max_levels      =  20

   !  The maximum number of nests that can depend on a single parent and other way round

   INTEGER , PARAMETER :: max_nests        =  20

   !  The maximum number of parents that a nest can have (simplified assumption -> one only)

   INTEGER , PARAMETER :: max_parents      =  1

   !  The maximum number of domains is how many grids the model will be running.

   INTEGER , PARAMETER :: max_domains     =   ( MAX_DOMAINS_F - 1 ) / 2 + 1

   !  The maximum number of nest move specifications allowed in a namelist

   INTEGER , PARAMETER :: max_moves       =   50

   !  The maximum number of eta levels

   INTEGER , PARAMETER :: max_eta         =   501

   !  The maximum number of outer iterations (for DA minimisation)

   INTEGER , PARAMETER :: max_outer_iterations = 10

   !  The maximum number of instruments (for radiance DA)

   INTEGER , PARAMETER :: max_instruments =   30

   !  The maximum number of bogus storms

   INTEGER , PARAMETER :: max_bogus =  5

#if ( WRF_PLANET == 1 )
   !  The maximum number of output files (wrfout + auxhist)
   !  Needed to allow declaration of maximum sizes for namelist arrays in the
   !  Registry file that will use this value
   !  As with MAX_DOMAINS_F above, MAX_HISTORY will be replaced by the value
   !  set in arch/preamble_new by the C-pre-processor (cpp or equivalent) when
   !  compiling this file.

   INTEGER , PARAMETER :: max_out_streams = MAX_HISTORY

! also look at values of kdm_n_band_xxx in Registry.EM_PLANET
#if ( WRF_MARS == 1 )
   INTEGER , PARAMETER :: max_n_band_ir = 7
   INTEGER , PARAMETER :: max_n_band_solar = 7
#elif ( WRF_TITAN == 1 )
   INTEGER , PARAMETER :: max_n_band_ir = 7
   INTEGER , PARAMETER :: max_n_band_solar = 4
#else
   INTEGER , PARAMETER :: max_n_band_ir = 10
   INTEGER , PARAMETER :: max_n_band_solar = 10
#endif
#endif

   !  2. Following related to driver leve data structures for DM_PARALLEL communications

#ifdef DM_PARALLEL
   INTEGER , PARAMETER :: max_comms       =   1024
#else
   INTEGER , PARAMETER :: max_comms       =   1
#endif

   !  3. Following is information related to the file I/O.

   !  These are the bounds of the available FORTRAN logical unit numbers for the file I/O.
   !  Only logical unti numbers within these bounds will be chosen for I/O unit numbers.

   INTEGER , PARAMETER :: min_file_unit = 10
   INTEGER , PARAMETER :: max_file_unit = 99

#if ( WRF_PLANET == 1 )
   ! Unfortunately, the below rump of code is no longer accurate nor
   ! necessary.  P2SI does NOT need to be defined here, and the #ifdef's
   ! are outdated anyway, so the code will always default to P2SI=1.0.
   ! The following block of code is also, unfortunately, a part of NCAR's
   ! archive and thus we don't have the power to remove it as obsolete.
   ! Thus, the whole block of code is #if/else'ed out of compilation here.
   ! In particular, even though (or perhaps, especially because) it is
   ! obsolete code, we don't want the defining of P2SI here, lest it conflict
   ! with the variable's true setting in share/module_model_constants.
#else
   !  4. Unfortunately, the following definition is needed here (rather
   !     than the more logical place in share/module_model_constants.F)
   !     for the namelist reads in frame/module_configure.F, and for some
   !     conversions in share/set_timekeeping.F
   !     Actually, using it here will mean that we don't need to set it
   !     in share/module_model_constants.F, since this file will be
   !     included (USEd) in:
   !        frame/module_configure.F
   !     which will be USEd in:
   !        share/module_bc.F
   !     which will be USEd in:
   !        phys/module_radiation_driver.F
   !     which is the other important place for it to be, and where
   !     it is passed as a subroutine parameter to any physics subroutine.
   !
   !     P2SI is the number of SI seconds in an planetary solar day
   !     divided by the number of SI seconds in an earth solar day
#if defined MARS
   !     For Mars, P2SI = 88775.2/86400.
   REAL , PARAMETER :: P2SI = 1.0274907
#elif defined TITAN
   !     For Titan, P2SI = 1378080.0/86400.
   REAL , PARAMETER :: P2SI = 15.95
#else
   !     Default for Earth
   REAL , PARAMETER :: P2SI = 1.0
#endif
#endif
 CONTAINS
   SUBROUTINE init_module_driver_constants
   END SUBROUTINE init_module_driver_constants
 END MODULE module_driver_constants

! routines that external packages can call to get at WRF stuff that isn't available
! through argument lists; since they are external we don't want them using WRF 
! modules unnecessarily (complicates the build even more)
 SUBROUTINE inquire_of_wrf_data_order_xyz( data_order )
   USE module_driver_constants, ONLY : DATA_ORDER_XYZ
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: data_order
   data_order = DATA_ORDER_XYZ
 END SUBROUTINE inquire_of_wrf_data_order_xyz

 SUBROUTINE inquire_of_wrf_data_order_xzy( data_order )
   USE module_driver_constants, ONLY : DATA_ORDER_XZY
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: data_order
   data_order = DATA_ORDER_XZY
 END SUBROUTINE inquire_of_wrf_data_order_xzy

 SUBROUTINE inquire_of_wrf_iwordsize( iwordsz )
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: iwordsz
   iwordsz = IWORDSIZE
 END SUBROUTINE inquire_of_wrf_iwordsize

 SUBROUTINE inquire_of_wrf_rwordsize( rwordsz )
   IMPLICIT NONE
   INTEGER, INTENT(OUT) :: rwordsz
   rwordsz = RWORDSIZE
 END SUBROUTINE inquire_of_wrf_rwordsize

