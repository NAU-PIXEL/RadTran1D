MODULE module_array_utilities

  IMPLICIT NONE 

  PRIVATE   ! make all variable declarations private unless
            ! explicity made public

  ! Make the subroutines we will actually call public
  PUBLIC :: locate, bilinear

  !----------------------------------------------------------------------
  ! Overloaded routines
  !----------------------------------------------------------------------
  INTERFACE locate
#if ( NO_QUAD == 1 )
     MODULE PROCEDURE locate_s
#else
     MODULE PROCEDURE locate_s, locate_d
#endif
  END INTERFACE

CONTAINS

!-----------------------------------------------------------------------------
!
! BILINEAR.F
!
!-----------------------------------------------------------------------------

  REAL FUNCTION BILINEAR(XA,YA,ZA,X,Y)
    ! Performs bilinear interpolation
    ! Given an array of ZA(j,k), where j varies from 1 to M, and k
    ! varies from 1 to N, and where XA(j) describes the values of the
    ! first index, and YA(k) describes the values of the second index
    ! (both must be monotonically increasing), this subroutine returns
    ! the value of Z as the value of ZA evaluated at the location (X,Y)
    IMPLICIT NONE

    ! Input/output variables
    REAL,    INTENT(IN   )                 :: X,Y
    REAL,    INTENT(IN   ), DIMENSION(:)   :: XA, YA
    REAL,    INTENT(IN   ), DIMENSION(:,:) :: ZA

    ! Local variables
    INTEGER :: M,N               ! Array sizes
    INTEGER :: J,K               ! Indices of locations near X,Y
    REAL    :: T,U,Z1,Z2,Z3,Z4   ! Temporary values

    ! Find M, N, the size of the input arrays
    M = SIZE(XA)
    N = SIZE(YA)

    ! Find j and k such that XA(j) <= X <= XA(j+1) and
    !                        YA(k) <= Y <= YA(k+1)
    j = LOCATE(XA,X)
    k = LOCATE(YA,Y)

    ! Here we are saying if our value lies outside our box range,
    ! just use the value at the edge of the range
    IF (j < 1)   j=1
    IF (j > M-1) j=M-1
    IF (k < 1)   k=1
    IF (k > N-1) k=N-1

    ! Define the values of the four points around X,Y
    Z1=ZA(J  ,K  )
    Z2=ZA(J+1,K  )
    Z3=ZA(J+1,K+1)
    Z4=ZA(J  ,K+1)

    ! Set up variables for convenience
    T =(X-XA(J))/(XA(J+1)-XA(J))
    U =(Y-YA(K))/(YA(K+1)-YA(K))
    IF (T < 0.) T = 0.
    IF (T > 1.) T = 1.
    IF (U < 0.) U = 0.
    IF (U > 1.) U = 1.

    ! Calculate the final value
    bilinear = (1.-T)*(1.-U)*Z1 + &
                   T *(1.-U)*Z2 + &
                   T *    U *Z3 + &
               (1.-T)*    U *Z4

  END FUNCTION BILINEAR

!-----------------------------------------------------------------------------
!
! LOCATE.F
!
!-----------------------------------------------------------------------------

  FUNCTION locate_s(y,x)
    ! Return a value j such that an input value x and input array y between
    ! y(j) and y(j+1).  Input array must either be monotonically increasing
    ! or decreasing.  A return value of 0 or n (where n is the size of array y)
    ! is used to indicate that x lies outised the range of y
    IMPLICIT NONE
    ! Input/output variables
    REAL, INTENT(IN), DIMENSION(:) :: y
    REAL, INTENT(IN) :: x
    INTEGER :: locate_s
    ! Local variables
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd

    n=SIZE(y)
    ascnd = (y(n) >= y(1))
    jl=0
    ju=n+1
    DO
       IF (ju-jl <= 1) EXIT
       jm=(ju+jl)/2
       IF (ascnd .EQV. (x >= y(jm))) THEN
          jl=jm
       ELSE
          ju=jm
       END IF
    END DO
    IF (x == y(1)) THEN
       locate_s=1
    ELSE IF (x == y(n)) THEN
       locate_s=n-1
    ELSE
       locate_s=jl
    END IF
  END FUNCTION locate_s

  FUNCTION locate_d(y,x)
    ! Return a value j such that an input value x and input array y between
    ! y(j) and y(j+1).  Input array must either be monotonically increasing
    ! or decreasing.  A return value of 0 or n (where n is the size of array y)
    ! is used to indicate that x lies outised the range of y
    IMPLICIT NONE
    ! Input/output variables
    REAL(KIND=KIND(0.d0)), INTENT(IN), DIMENSION(:) :: y
    REAL(KIND=KIND(0.d0)), INTENT(IN) :: x
    INTEGER :: locate_d
    ! Local variables
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd

    n=SIZE(y)
    ascnd = (y(n) >= y(1))
    jl=0
    ju=n+1
    DO
       IF (ju-jl <= 1) EXIT
       jm=(ju+jl)/2
       IF (ascnd .EQV. (x >= y(jm))) THEN
          jl=jm
       ELSE
          ju=jm
       END IF
    END DO
    IF (x == y(1)) THEN
       locate_d=1
    ELSE IF (x == y(n)) THEN
       locate_d=n-1
    ELSE
       locate_d=jl
    END IF
  END FUNCTION locate_d

END MODULE module_array_utilities
