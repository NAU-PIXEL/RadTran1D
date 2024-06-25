!WRF:DRIVER_LAYER:UTIL
!

MODULE module_wrf_error
  INTEGER           :: wrf_debug_level = 0
  CHARACTER*256     :: wrf_err_message
!$OMP THREADPRIVATE (wrf_err_message)
CONTAINS

SUBROUTINE wrf_message( str )
  IMPLICIT NONE
  CHARACTER*(*) str
  PRINT*, TRIM(str)
END SUBROUTINE wrf_message

SUBROUTINE wrf_error_fatal( str )
  IMPLICIT NONE
  CHARACTER*(*) str
  PRINT*, TRIM(str)
  STOP
END SUBROUTINE wrf_error_fatal

END MODULE module_wrf_error
