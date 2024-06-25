MODULE module_nrutils

  IMPLICIT NONE

  PRIVATE   ! make all variable and routine declarations private unless
            ! explicity made public

  ! Make the subroutines we will actually call public
  PUBLIC :: locate
  PUBLIC :: linear_interpolate
  PUBLIC :: bilinear_interpolate
  PUBLIC :: trilinear_interpolate
  PUBLIC :: quadrilinear_interpolate
  PUBLIC :: four1, fourrow, sort2, indexx
  PUBLIC :: tridag
  PUBLIC :: matinv, mprove, lubksb, ludcmp
  PUBLIC :: period, avevar
  PUBLIC :: spline, splint
  PUBLIC :: tred2, tqli, eigsrt
  PUBLIC :: balanc, elmhes, hqr
  PUBLIC :: polint
  PUBLIC :: zbrak, zroots, rtsafe
  ! The following two added to public so that the NR routines hacked in
  ! dyn_em/module_setup_h_grid.F can use them too.
  PUBLIC :: arth
  PUBLIC :: nrerror
  ! The following nrutil routines added so that other subroutines
  ! can use them as well
  PUBLIC :: iminloc, imaxloc

  ! NRTYPE.F

  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  TYPE sprs2_sp
     INTEGER(I4B) :: n,len
     REAL(SP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  TYPE sprs2_dp
     INTEGER(I4B) :: n,len
     REAL(DP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp

  ! NRUTIL.F variable declarations

  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
  INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
  INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
  INTEGER(I4B), PARAMETER :: NPAR_POLY=8
  INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8

  !----------------------------------------------------------------------
  ! Special interpolation status variables 
  !----------------------------------------------------------------------
  INTEGER, PARAMETER, PUBLIC :: signal_interpolate_no_error                 = 0
  INTEGER, PARAMETER, PUBLIC :: signal_interpolate_interpolant_out_of_range = 1

  !----------------------------------------------------------------------
  ! Overloaded routines
  !----------------------------------------------------------------------

  ! NR.F interfaces

  INTERFACE locate
#if (NO_QUAD == 1)
     MODULE PROCEDURE locate_s
#else
     MODULE PROCEDURE locate_s, locate_d
#endif
  END INTERFACE

  INTERFACE linear_interpolate
#if (NO_QUAD == 1)
     MODULE PROCEDURE linear_interpolate_s, &
                    bilinear_interpolate_s, &
                   trilinear_interpolate_s, &
                quadrilinear_interpolate_s
#else
     MODULE PROCEDURE linear_interpolate_s, linear_interpolate_d, &
                  bilinear_interpolate_s, bilinear_interpolate_d, &
                trilinear_interpolate_s, trilinear_interpolate_d, &
          quadrilinear_interpolate_s, quadrilinear_interpolate_d
#endif
  END INTERFACE

  INTERFACE bilinear_interpolate
#if (NO_QUAD == 1)
     MODULE PROCEDURE bilinear_interpolate_s
#else
     MODULE PROCEDURE bilinear_interpolate_s, bilinear_interpolate_d
#endif
  END INTERFACE

  INTERFACE trilinear_interpolate
#if (NO_QUAD == 1)
     MODULE PROCEDURE trilinear_interpolate_s
#else
     MODULE PROCEDURE trilinear_interpolate_s, trilinear_interpolate_d
#endif
  END INTERFACE

  INTERFACE quadrilinear_interpolate
#if (NO_QUAD == 1)
     MODULE PROCEDURE quadrilinear_interpolate_s
#else
     MODULE PROCEDURE quadrilinear_interpolate_s, quadrilinear_interpolate_d
#endif
  END INTERFACE

  INTERFACE distance
#if (NO_QUAD == 1)
     MODULE PROCEDURE distance_d
#else
     MODULE PROCEDURE distance_s, distance_d
#endif
  END INTERFACE

!  INTERFACE four1
!     SUBROUTINE four1_dp(data,isign)
!       COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
!       INTEGER(I4B), INTENT(IN) :: isign
!     END SUBROUTINE four1_dp
!
!     SUBROUTINE four1_sp(data,isign)
!       COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
!       INTEGER(I4B), INTENT(IN) :: isign
!     END SUBROUTINE four1_sp
!  END INTERFACE

!  INTERFACE fourrow
!     SUBROUTINE fourrow_dp(data,isign)
!       COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
!       INTEGER(I4B), INTENT(IN) :: isign
!     END SUBROUTINE fourrow_dp
!
!     SUBROUTINE fourrow_sp(data,isign)
!       COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
!       INTEGER(I4B), INTENT(IN) :: isign
!     END SUBROUTINE fourrow_sp
!  END INTERFACE

  INTERFACE four1
!#if (RWORDSIZE == DWORDSIZE)
#if ( NO_QUAD == 1 )
     MODULE PROCEDURE four1_sp
#else
     MODULE PROCEDURE four1_dp, four1_sp
#endif
  END INTERFACE

  INTERFACE fourrow
#if ( NO_QUAD == 1)
     MODULE PROCEDURE fourrow_sp
#else
     MODULE PROCEDURE fourrow_dp, fourrow_sp
#endif
  END INTERFACE

  INTERFACE tridag
#if ( NO_QUAD == 1)
     MODULE PROCEDURE tridag_r
#else
     MODULE PROCEDURE tridag_r, tridag_d
#endif
  END INTERFACE

  INTERFACE matinv
#if ( NO_QUAD == 1)
     MODULE PROCEDURE matinv_r
#else
     MODULE PROCEDURE matinv_r, matinv_d
#endif
  END INTERFACE

  INTERFACE lubksb
#if ( NO_QUAD == 1)
     MODULE PROCEDURE lubksb_r
#else
     MODULE PROCEDURE lubksb_r, lubksb_d
#endif
  END INTERFACE

  INTERFACE ludcmp
#if ( NO_QUAD == 1)
     MODULE PROCEDURE ludcmp_r
#else
     MODULE PROCEDURE ludcmp_r, ludcmp_d
#endif
  END INTERFACE

  INTERFACE pythag
#if ( NO_QUAD == 1)
     MODULE PROCEDURE pythag_sp
#else
     MODULE PROCEDURE pythag_dp, pythag_sp
#endif
  END INTERFACE

  INTERFACE tred2
#if ( NO_QUAD == 1)
     MODULE PROCEDURE tred2_sp
#else
     MODULE PROCEDURE tred2_dp, tred2_sp
#endif
  END INTERFACE

  INTERFACE tqli
#if ( NO_QUAD == 1)
     MODULE PROCEDURE tqli_sp
#else
     MODULE PROCEDURE tqli_dp, tqli_sp
#endif
  END INTERFACE

  INTERFACE eigsrt
#if ( NO_QUAD == 1)
     MODULE PROCEDURE eigsrt_sp
#else
     MODULE PROCEDURE eigsrt_dp, eigsrt_sp
#endif
  END INTERFACE

  INTERFACE spline
#if ( NO_QUAD == 1)
     MODULE PROCEDURE spline_r
#else
     MODULE PROCEDURE spline_r, spline_d
#endif
  END INTERFACE

  INTERFACE splint
#if ( NO_QUAD == 1)
     MODULE PROCEDURE splint_r
#else
     MODULE PROCEDURE splint_r, splint_d
#endif
  END INTERFACE

  INTERFACE polint
#if ( NO_QUAD == 1)
     MODULE PROCEDURE polint_r
#else
     MODULE PROCEDURE polint_r, polint_d
#endif
  END INTERFACE

  INTERFACE zbrak
#if ( NO_QUAD == 1)
     MODULE PROCEDURE zbrak_r, zbrak2_r
#else
     MODULE PROCEDURE zbrak_r, zbrak_d, zbrak2_r, zbrak2_d
#endif
  END INTERFACE

  INTERFACE zroots
#if ( NO_QUAD == 1)
     MODULE PROCEDURE zroots_c
#else
     MODULE PROCEDURE zroots_c, zroots_z
#endif
  END INTERFACE

  INTERFACE laguer
#if ( NO_QUAD == 1)
     MODULE PROCEDURE laguer_c
#else
     MODULE PROCEDURE laguer_c, laguer_z
#endif
  END INTERFACE

  INTERFACE rtsafe
#if ( NO_QUAD == 1)
     MODULE PROCEDURE rtsafe_r, rtsafe2_r, rtsafe3_r
#else
     MODULE PROCEDURE rtsafe_r,  rtsafe_d, &
                      rtsafe2_r, rtsafe2_d, &
                      rtsafe3_r, rtsafe3_d
#endif
  END INTERFACE

  ! NRUTIL.F interfaces

  INTERFACE arth
#if ( NO_QUAD == 1)
     MODULE PROCEDURE arth_r, arth_i
#else
     MODULE PROCEDURE arth_r, arth_d, arth_i
#endif
  END INTERFACE

  INTERFACE assert
     MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE

  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE

  INTERFACE icomp_xchg
#if ( NO_QUAD == 1)
     MODULE PROCEDURE icomp_xchg_sp, icomp_xchg_i4b
#else
     MODULE PROCEDURE icomp_xchg_sp, icomp_xchg_dp, icomp_xchg_i4b
#endif
  END INTERFACE

  INTERFACE imaxloc
#if ( NO_QUAD == 1)
     MODULE PROCEDURE imaxloc_r, imaxloc_i
#else
     MODULE PROCEDURE imaxloc_d, imaxloc_r, imaxloc_i
#endif
  END INTERFACE

  INTERFACE iminloc
#if ( NO_QUAD == 1)
     MODULE PROCEDURE iminloc_r, iminloc_i
#else
     MODULE PROCEDURE iminloc_d, iminloc_r, iminloc_i
#endif
  END INTERFACE

  INTERFACE indexx
#if ( NO_QUAD == 1)
     MODULE PROCEDURE indexx_sp, indexx_i4b
#else
     MODULE PROCEDURE indexx_sp, indexx_dp, indexx_i4b
#endif
  END INTERFACE

  INTERFACE outerprod
#if ( NO_QUAD == 1)
     MODULE PROCEDURE outerprod_r
#else
     MODULE PROCEDURE outerprod_r, outerprod_d
#endif
  END INTERFACE

  INTERFACE outerdiff
#if ( NO_QUAD == 1)
     MODULE PROCEDURE outerdiff_r, outerdiff_i
#else
     MODULE PROCEDURE outerdiff_r, outerdiff_d, outerdiff_i
#endif
  END INTERFACE

  INTERFACE swap
#if ( NO_QUAD == 1)
     MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
          swap_cv,swap_cm, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
#else
     MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_d,swap_dv,swap_c, &
          swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
          masked_swap_rs,masked_swap_rv,masked_swap_rm
#endif
  END INTERFACE

  INTERFACE diagadd
     MODULE PROCEDURE diagadd_rv,diagadd_r
  END INTERFACE

  INTERFACE poly
#if ( NO_QUAD == 1)
     MODULE PROCEDURE poly_rr,poly_rrv,&
          poly_rc,poly_cc,poly_msk_rrv
#else
     MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
          poly_rc,poly_cc,poly_zz,poly_msk_rrv,poly_msk_ddv
#endif
  END INTERFACE

  INTERFACE poly_term
#if ( NO_QUAD == 1)
     MODULE PROCEDURE poly_term_rr,poly_term_cc
#else
     MODULE PROCEDURE poly_term_rr,poly_term_cc,poly_term_dd,poly_term_zz
#endif
  END INTERFACE

CONTAINS

!-----------------------------------------------------------------------------
!
! NRUTILS.F
!
!-----------------------------------------------------------------------------
! assert
!-----------------------------------------------------------------------------

  SUBROUTINE assert1(n1,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if (.not. n1) then
       write (*,*) 'nrerror: an assertion failed with this tag:',string
       STOP 'program terminated by assert1'
    end if
  END SUBROUTINE assert1

  SUBROUTINE assert2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2
    if (.not. (n1 .and. n2)) then
       write (*,*) 'nrerror: an assertion failed with this tag:',string
       STOP 'program terminated by assert2'
    end if
  END SUBROUTINE assert2

  SUBROUTINE assert3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
       write (*,*) 'nrerror: an assertion failed with this tag:',string
       STOP 'program terminated by assert3'
    end if
  END SUBROUTINE assert3

  SUBROUTINE assert4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
       write (*,*) 'nrerror: an assertion failed with this tag:',string
       STOP 'program terminated by assert4'
    end if
  END SUBROUTINE assert4

  SUBROUTINE assert_v(n,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, DIMENSION(:), INTENT(IN) :: n
    if (.not. all(n)) then
       write (*,*) 'nrerror: an assertion failed with this tag:',string
       STOP 'program terminated by assert_v'
    end if
  END SUBROUTINE assert_v

!-----------------------------------------------------------------------------
! assert_eq
!-----------------------------------------------------------------------------

  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2

  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3

  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn

!-----------------------------------------------------------------------------
! arth
!-----------------------------------------------------------------------------

  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r

  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d

  FUNCTION arth_i(first,increment,n)
    INTEGER(I4B), INTENT(IN) :: first,increment,n
    INTEGER(I4B), DIMENSION(n) :: arth_i
    INTEGER(I4B) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

!-----------------------------------------------------------------------------
! swap
!-----------------------------------------------------------------------------

  SUBROUTINE swap_i(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i

  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r

  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv

  SUBROUTINE swap_d(a,b)
    REAL(DP), INTENT(INOUT) :: a,b
    REAL(DP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_d

  SUBROUTINE swap_dv(a,b)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(DP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_dv

  SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c

  SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv

  SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm

  SUBROUTINE swap_z(a,b)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z

  SUBROUTINE swap_zv(a,b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv

  SUBROUTINE swap_zm(a,b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm

  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(SP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs

  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv

  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm

!-----------------------------------------------------------------------------
! nrerror
!-----------------------------------------------------------------------------

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

!-----------------------------------------------------------------------------
! outerprod/outerdiff
!-----------------------------------------------------------------------------

  FUNCTION outerprod_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = SPREAD(a,dim=2,ncopies=size(b)) * &
                  SPREAD(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r

  FUNCTION outerprod_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = SPREAD(a,dim=2,ncopies=size(b)) * &
                  SPREAD(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_d

  FUNCTION outerdiff_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
                  spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r

  FUNCTION outerdiff_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
    outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
                  spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_d

  FUNCTION outerdiff_i(a,b)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
    INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
                  spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i

!-----------------------------------------------------------------------------
! imaxloc/iminloc
!-----------------------------------------------------------------------------

  FUNCTION imaxloc_d(arr)
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B) :: imaxloc_d
    INTEGER(I4B), DIMENSION(1) :: imax
    imax=MAXLOC(arr(:))
    imaxloc_d=imax(1)
  END FUNCTION imaxloc_d

  FUNCTION imaxloc_r(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B) :: imaxloc_r
    INTEGER(I4B), DIMENSION(1) :: imax
    imax=MAXLOC(arr(:))
    imaxloc_r=imax(1)
  END FUNCTION imaxloc_r

  FUNCTION imaxloc_i(iarr)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(I4B), DIMENSION(1) :: imax
    INTEGER(I4B) :: imaxloc_i
    imax=MAXLOC(iarr(:))
    imaxloc_i=imax(1)
  END FUNCTION imaxloc_i

  FUNCTION iminloc_d(arr)
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc_d
    imin=MINLOC(arr(:))
    iminloc_d=imin(1)
  END FUNCTION iminloc_d

  FUNCTION iminloc_r(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc_r
    imin=MINLOC(arr(:))
    iminloc_r=imin(1)
  END FUNCTION iminloc_r

  FUNCTION iminloc_i(arr)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc_i
    imin=MINLOC(arr(:))
    iminloc_i=imin(1)
  END FUNCTION iminloc_i

!-----------------------------------------------------------------------------
! diagadd
!-----------------------------------------------------------------------------

  SUBROUTINE diagadd_rv(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
    do j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    end do
  END SUBROUTINE diagadd_rv

  SUBROUTINE diagadd_r(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)+diag
    end do
  END SUBROUTINE diagadd_r

!-----------------------------------------------------------------------------
! upper_triangle
!-----------------------------------------------------------------------------

  FUNCTION upper_triangle(j,k,extra)
    INTEGER(I4B), INTENT(IN) :: j,k
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
    LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
    INTEGER(I4B) :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle

!-----------------------------------------------------------------------------
! poly
!-----------------------------------------------------------------------------

  FUNCTION poly_rr(x,coeffs)
    REAL(SP), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
    REAL(SP) :: poly_rr
    REAL(SP) :: pow
    REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rr=0.0_sp
    else if (n < NPAR_POLY) then
       poly_rr=coeffs(n)
       do i=n-1,1,-1
          poly_rr=x*poly_rr+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rr=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_rr

  FUNCTION poly_dd(x,coeffs)
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
    REAL(DP) :: poly_dd
    REAL(DP) :: pow
    REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_dd=0.0_dp
    else if (n < NPAR_POLY) then
       poly_dd=coeffs(n)
       do i=n-1,1,-1
          poly_dd=x*poly_dd+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_dd=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_dd

  FUNCTION poly_rc(x,coeffs)
    COMPLEX(SPC), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(SPC) :: poly_rc
    COMPLEX(SPC) :: pow
    COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rc=0.0_sp
    else if (n < NPAR_POLY) then
       poly_rc=coeffs(n)
       do i=n-1,1,-1
          poly_rc=x*poly_rc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rc=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_rc

  FUNCTION poly_cc(x,coeffs)
    COMPLEX(SPC), INTENT(IN) :: x
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(SPC) :: poly_cc
    COMPLEX(SPC) :: pow
    COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_cc=0.0_sp
    else if (n < NPAR_POLY) then
       poly_cc=coeffs(n)
       do i=n-1,1,-1
          poly_cc=x*poly_cc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_cc=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_cc

  FUNCTION poly_zz(x,coeffs)
    COMPLEX(DPC), INTENT(IN) :: x
    COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(DPC) :: poly_zz
    COMPLEX(DPC) :: pow
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_zz=0.0_dp
    else if (n < NPAR_POLY) then
       poly_zz=coeffs(n)
       do i=n-1,1,-1
          poly_zz=x*poly_zz+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_zz=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_zz

  FUNCTION poly_rrv(x,coeffs)
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
    REAL(SP), DIMENSION(size(x)) :: poly_rrv
    INTEGER(I4B) :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_rrv=0.0_sp
    else if (m < n .or. m < NPAR_POLY) then
       poly_rrv=coeffs(m)
       do i=m-1,1,-1
          poly_rrv=x*poly_rrv+coeffs(i)
       end do
    else
       do i=1,n
          poly_rrv(i)=poly_rr(x(i),coeffs)
       end do
    end if
  END FUNCTION poly_rrv

  FUNCTION poly_ddv(x,coeffs)
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
    REAL(DP), DIMENSION(size(x)) :: poly_ddv
    INTEGER(I4B) :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_ddv=0.0_dp
    else if (m < n .or. m < NPAR_POLY) then
       poly_ddv=coeffs(m)
       do i=m-1,1,-1
          poly_ddv=x*poly_ddv+coeffs(i)
       end do
    else
       do i=1,n
          poly_ddv(i)=poly_dd(x(i),coeffs)
       end do
    end if
  END FUNCTION poly_ddv

  FUNCTION poly_msk_rrv(x,coeffs,mask)
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
    poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
  END FUNCTION poly_msk_rrv

  FUNCTION poly_msk_ddv(x,coeffs,mask)
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
    poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
  END FUNCTION poly_msk_ddv

!-----------------------------------------------------------------------------
! poly_term
!-----------------------------------------------------------------------------

  RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a
    REAL(SP), INTENT(IN) :: b
    REAL(SP), DIMENSION(size(a)) :: u
    INTEGER(I4B) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_rr

  RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(SPC), INTENT(IN) :: b
    COMPLEX(SPC), DIMENSION(size(a)) :: u
    INTEGER(I4B) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_cc

  RECURSIVE FUNCTION poly_term_dd(a,b) RESULT(u)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a
    REAL(DP), INTENT(IN) :: b
    REAL(DP), DIMENSION(size(a)) :: u
    INTEGER(I4B) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_dd(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_dd

  RECURSIVE FUNCTION poly_term_zz(a,b) RESULT(u)
    COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(DPC), INTENT(IN) :: b
    COMPLEX(DPC), DIMENSION(size(a)) :: u
    INTEGER(I4B) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_zz(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_zz

!-----------------------------------------------------------------------------
!
! LOCATE.F
!
!-----------------------------------------------------------------------------

  INTEGER FUNCTION locate_s(y,x)
    ! Return a value j such that an input value x and input array y between
    ! y(j) and y(j+1).  Input array must either be monotonically increasing
    ! or decreasing.  A return value of 0 or n (where n is the size of array y)
    ! is used to indicate that x lies outised the range of y
    IMPLICIT NONE
    ! Input/output variables
    REAL, INTENT(IN), DIMENSION(:) :: y
    REAL, INTENT(IN) :: x
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


  INTEGER FUNCTION locate_d(y,x)
    ! Return a value j such that an input value x and input array y between
    ! y(j) and y(j+1).  Input array must either be monotonically increasing
    ! or decreasing.  A return value of 0 or n (where n is the size of array y)
    ! is used to indicate that x lies outised the range of y
    IMPLICIT NONE
    ! Input/output variables
    REAL(KIND=KIND(0d0)), INTENT(IN), DIMENSION(:) :: y
    REAL(KIND=KIND(0d0)), INTENT(IN) :: x
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

!-----------------------------------------------------------------------------
!
! N-DIMENSIONAL (N=1,2,3,4) LINEAR INTERPOLATION
! Not a specific library routine, but following the standard methodology
! as also explained in the text.
!
!-----------------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Distance
  ! Bookkeeping routine: fractional distance of interpolant between
  ! nearest independent variable gridded values (from lower indexed
  ! location)
  !----------------------------------------------------------------------
  REAL          FUNCTION distance_s(x0, x, index, log_abscissa, &
       outside_to_edge, extrapolate) &
       RESULT(d)
    IMPLICIT NONE
    REAL, INTENT(IN), DIMENSION(:) :: x
    REAL, INTENT(IN)               :: x0
    INTEGER, INTENT(INOUT)         :: index
    LOGICAL, INTENT(IN)            :: log_abscissa
    LOGICAL, INTENT(IN)            :: outside_to_edge
    LOGICAL, INTENT(IN)            :: extrapolate

    IF (outside_to_edge .AND. (index == 0) ) THEN
       index = 1
       d     = 0.
    ELSE IF (outside_to_edge .AND. (index == SIZE(x)) ) THEN
       index = SIZE(x)-1
       d     = 1.
    ELSE
       IF (extrapolate .AND. (index == 0) ) THEN
          index = 1
       ELSE IF (extrapolate .AND. (index == SIZE(x)) ) THEN
          index = SIZE(x)-1
       END IF
       
       IF (log_abscissa) THEN
          d = LOG(x0/x(index))/LOG(x(index+1)/x(index))
       ELSE
          d = (x0-x(index))/(x(index+1)-x(index))
       END IF
    END IF
  END FUNCTION distance_s

  REAL(KIND=KIND(0.d0)) FUNCTION distance_d(x0, x, index, log_abscissa, &
       outside_to_edge, extrapolate) &
       RESULT(d)
    IMPLICIT NONE
    REAL(KIND=KIND(0.d0)), INTENT(IN), DIMENSION(:) :: x
    REAL(KIND=KIND(0.d0)), INTENT(IN)               :: x0
    INTEGER, INTENT(INOUT)                          :: index
    LOGICAL, INTENT(IN)                             :: log_abscissa
    LOGICAL, INTENT(IN)                             :: outside_to_edge
    LOGICAL, INTENT(IN)                             :: extrapolate

    IF (outside_to_edge .AND. (index == 0) ) THEN
       index = 1
       d     = 0.d0
    ELSE IF (outside_to_edge .AND. (index == SIZE(x)) ) THEN
       index = SIZE(x)-1
       d     = 1.d0
    ELSE
       IF (extrapolate .AND. (index == 0) ) THEN
          index = 1
       ELSE IF (extrapolate .AND. (index == SIZE(x)) ) THEN
          index = SIZE(x)-1
       END IF
       
       IF (log_abscissa) THEN
          d = LOG(x0/x(index))/LOG(x(index+1)/x(index))
       ELSE
          d = (x0-x(index))/(x(index+1)-x(index))
       END IF
    END IF
  END FUNCTION distance_d
  


  !----------------------------------------------------------------------
  ! Linear interpolation
  !----------------------------------------------------------------------
  REAL          FUNCTION linear_interpolate_s(x0, x, y, log_x, bad_flag, &
       out_of_range_use_nearest_edge, linear_extrapolate,                &
       error_status)                                                     &
       RESULT(y0)

    IMPLICIT NONE
    REAL, INTENT(IN), DIMENSION(:) :: x, y
    REAL, INTENT(IN)               :: x0
    LOGICAL, INTENT(IN), OPTIONAL  :: log_x
    REAL, INTENT(IN), OPTIONAL     :: bad_flag
    LOGICAL, INTENT(IN), OPTIONAL  :: out_of_range_use_nearest_edge
    LOGICAL, INTENT(IN), OPTIONAL  :: linear_extrapolate
    INTEGER, INTENT(OUT), OPTIONAL :: error_status
    LOGICAL :: log_abscissa
    LOGICAL :: outside_to_edge
    LOGICAL :: extrapolate
    INTEGER :: i
    REAL :: d

    IF (PRESENT(error_status)) error_status = 0

    IF (PRESENT(linear_extrapolate)) THEN
       extrapolate = linear_extrapolate
    ELSE
       extrapolate = .FALSE.
    END IF

    IF (PRESENT(log_x)) THEN
       log_abscissa = log_x
    ELSE
       log_abscissa = .FALSE.
    END IF

    IF (log_abscissa) THEN
       i = LOCATE(LOG(x), LOG(x0))
    ELSE
       i = LOCATE(x, x0)
    END IF

    IF (PRESENT(out_of_range_use_nearest_edge)) THEN
       outside_to_edge = out_of_range_use_nearest_edge
    ELSE
       outside_to_edge = .FALSE.
    END IF
    IF ( .NOT.outside_to_edge .AND. &
         .NOT.extrapolate     .AND. &
         ((i == 0) .OR. (i == SIZE(x)))) THEN
       IF (PRESENT(error_status)) THEN
          error_status = signal_interpolate_interpolant_out_of_range
          y0 = HUGE(0.)
          RETURN
       ELSE
          WRITE(0,*) 'ERROR: linear_interpolate_s: interpolant outside input range'
          WRITE(0,*) x0, x(1), x(SIZE(x))
          STOP
       END IF
    END IF

    IF (PRESENT(bad_flag)) THEN
       IF ( (y(i)   <= bad_flag) .OR. &
            (y(i+1) <= bad_flag) ) THEN
          WRITE(0,*) 'ERROR: linear_interpolate_s: bad value(s):', &
               y(i), y(i+1)
          STOP
       END IF
    END IF

    d = distance(x0, x, i, log_abscissa, outside_to_edge, extrapolate)

    y0 = y(i) + (y(i+1)-y(i))*d

  END FUNCTION linear_interpolate_s
  
  REAL(KIND=KIND(0.d0)) FUNCTION linear_interpolate_d(x0, x, y, log_x, &
       bad_flag,                                                       &
       out_of_range_use_nearest_edge, linear_extrapolate,              &
       error_status)                                                   &
       RESULT(y0)

    IMPLICIT NONE
    REAL(KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:) :: x, y
    REAL(KIND=KIND(0.D0)), INTENT(IN)               :: x0
    LOGICAL, INTENT(IN), OPTIONAL           :: log_x
    REAL(KIND=KIND(0.D0)), INTENT(IN), OPTIONAL     :: bad_flag
    LOGICAL, INTENT(IN), OPTIONAL           :: out_of_range_use_nearest_edge
    LOGICAL, INTENT(IN), OPTIONAL           :: linear_extrapolate
    INTEGER, INTENT(OUT), OPTIONAL          :: error_status
    LOGICAL :: log_abscissa
    LOGICAL :: outside_to_edge
    LOGICAL :: extrapolate
    INTEGER :: i
    REAL :: d

    IF (PRESENT(error_status)) error_status = 0

    IF (PRESENT(linear_extrapolate)) THEN
       extrapolate = linear_extrapolate
    ELSE
       extrapolate = .FALSE.
    END IF

    IF (PRESENT(log_x)) THEN
       log_abscissa = log_x
    ELSE
       log_abscissa = .FALSE.
    END IF

    IF (log_abscissa) THEN
       i = LOCATE(LOG(x), LOG(x0))
    ELSE
       i = LOCATE(x, x0)
    END IF

    IF (PRESENT(out_of_range_use_nearest_edge)) THEN
       outside_to_edge = out_of_range_use_nearest_edge
    ELSE
       outside_to_edge = .FALSE.
    END IF
    IF ( .NOT.outside_to_edge .AND. &
         .NOT.extrapolate     .AND. &
         ((i == 0) .OR. (i == SIZE(x)))) THEN
       IF (PRESENT(error_status)) THEN
          error_status = signal_interpolate_interpolant_out_of_range
          y0 = HUGE(0.d0)
          RETURN
       ELSE
          WRITE(0,*) 'ERROR: linear_interpolate_d: interpolant outside input range'
          WRITE(0,*) x0, x(1), x(SIZE(x))
          STOP
       END IF
    END IF

    IF (PRESENT(bad_flag)) THEN
       IF ( (y(i)   <= bad_flag) .OR. &
            (y(i+1) <= bad_flag) ) THEN
          WRITE(0,*) 'ERROR: linear_interpolate_d: bad value(s):', &
               y(i), y(i+1)
          STOP
       END IF
    END IF

    d = distance(x0, x, i, log_abscissa, outside_to_edge, extrapolate)

    y0 = y(i) + (y(i+1)-y(i))*d

  END FUNCTION linear_interpolate_d
  


  !----------------------------------------------------------------------
  ! Bilinear interpolation
  !----------------------------------------------------------------------
  REAL          FUNCTION bilinear_interpolate_s(x1_0, x2_0, x1, x2, y, &
       log_x1, log_x2, bad_flag, out_of_range_use_nearest_edge,        &
       linear_extrapolate,                                             &
       error_status)                                                   &
       RESULT(y0)

    IMPLICIT NONE
    REAL, INTENT(IN), DIMENSION(:,:) :: y
    REAL, INTENT(IN), DIMENSION(:)   :: x1, x2
    REAL, INTENT(IN)                 :: x1_0, x2_0
    LOGICAL, INTENT(IN), OPTIONAL    :: log_x1, log_x2
    REAL, INTENT(IN), OPTIONAL       :: bad_flag
    LOGICAL, INTENT(IN), OPTIONAL    :: out_of_range_use_nearest_edge
    LOGICAL, INTENT(IN), OPTIONAL    :: linear_extrapolate
    INTEGER, INTENT(OUT), OPTIONAL   :: error_status
    LOGICAL :: log_abscissa1, log_abscissa2
    LOGICAL :: outside_to_edge
    LOGICAL :: extrapolate
    INTEGER :: i1, i2
    REAL :: d1, d2

    IF (PRESENT(error_status)) error_status = 0

    IF (PRESENT(linear_extrapolate)) THEN
       extrapolate = linear_extrapolate
    ELSE
       extrapolate = .FALSE.
    END IF

    IF (PRESENT(log_x1)) THEN
       log_abscissa1 = log_x1
    ELSE
       log_abscissa1 = .FALSE.
    END IF
    IF (PRESENT(log_x2)) THEN
       log_abscissa2 = log_x2
    ELSE
       log_abscissa2 = .FALSE.
    END IF

    IF (log_abscissa1) THEN
       i1 = LOCATE(LOG(x1), LOG(x1_0))
    ELSE
       i1 = LOCATE(x1, x1_0)
    END IF
    IF (log_abscissa2) THEN
       i2 = LOCATE(LOG(x2), LOG(x2_0))
    ELSE
       i2 = LOCATE(x2, x2_0)
    END IF

    IF (PRESENT(out_of_range_use_nearest_edge)) THEN
       outside_to_edge = out_of_range_use_nearest_edge
    ELSE
       outside_to_edge = .FALSE.
    END IF
    IF (.NOT.outside_to_edge .AND. &
        .NOT.extrapolate     .AND. &
         ( (i1 == 0) .OR. (i1 == SIZE(x1)) .OR. &
           (i2 == 0) .OR. (i2 == SIZE(x2)) ) ) THEN
       IF (PRESENT(error_status)) THEN
          error_status = signal_interpolate_interpolant_out_of_range
          y0 = HUGE(0.)
          RETURN
       ELSE
          WRITE(0,*) 'ERROR: bilinear_interpolate_s: interpolant outside input range'
          WRITE(0,*) x1_0, x1(1), x1(SIZE(x1))
          WRITE(0,*) x2_0, x2(1), x2(SIZE(x2))
          STOP
       END IF
    END IF

    IF (PRESENT(bad_flag)) THEN
       IF ( (y(i1  ,i2  ) <= bad_flag) .OR. &
            (y(i1+1,i2  ) <= bad_flag) .OR. &
            (y(i1  ,i2+1) <= bad_flag) .OR. &
            (y(i1+1,i2+1) <= bad_flag) ) THEN
          WRITE(0,*) 'ERROR: bilinear_interpolate_s: bad value(s):', &
               y(i1  ,i2  ), y(i1+1,i2  ), y(i1  ,i2+1), y(i1+1,i2+1)
          STOP
       END IF
    END IF

    d1 = distance(x1_0, x1, i1, log_abscissa1, outside_to_edge, extrapolate)
    d2 = distance(x2_0, x2, i2, log_abscissa2, outside_to_edge, extrapolate)

    y0 = (1.-d1)*(1.-d2)*y(i1  ,i2  ) + &
         (   d1)*(1.-d2)*y(i1+1,i2  ) + &
         (1.-d1)*(   d2)*y(i1  ,i2+1) + &
         (   d1)*(   d2)*y(i1+1,i2+1)

  END FUNCTION bilinear_interpolate_s

  REAL(KIND=KIND(0.D0)) FUNCTION bilinear_interpolate_d(x1_0, x2_0, x1, x2, y, &
       log_x1, log_x2, bad_flag, out_of_range_use_nearest_edge,        &
       linear_extrapolate,                                             &
       error_status)                                                   &
       RESULT(y0)

    IMPLICIT NONE
    REAL(KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:,:) :: y
    REAL(KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:)   :: x1, x2
    REAL(KIND=KIND(0.D0)), INTENT(IN)                 :: x1_0, x2_0
    LOGICAL, INTENT(IN), OPTIONAL             :: log_x1, log_x2
    REAL(KIND=KIND(0.D0)), INTENT(IN), OPTIONAL       :: bad_flag
    LOGICAL, INTENT(IN), OPTIONAL             :: out_of_range_use_nearest_edge
    LOGICAL, INTENT(IN), OPTIONAL             :: linear_extrapolate
    INTEGER, INTENT(OUT), OPTIONAL            :: error_status
    LOGICAL :: log_abscissa1, log_abscissa2
    LOGICAL :: outside_to_edge
    LOGICAL :: extrapolate
    INTEGER :: i1, i2
    REAL(KIND=KIND(0.D0)) :: d1, d2

    IF (PRESENT(error_status)) error_status = 0

    IF (PRESENT(linear_extrapolate)) THEN
       extrapolate = linear_extrapolate
    ELSE
       extrapolate = .FALSE.
    END IF

    IF (PRESENT(log_x1)) THEN
       log_abscissa1 = log_x1
    ELSE
       log_abscissa1 = .FALSE.
    END IF
    IF (PRESENT(log_x2)) THEN
       log_abscissa2 = log_x2
    ELSE
       log_abscissa2 = .FALSE.
    END IF

    IF (log_abscissa1) THEN
       i1 = LOCATE(LOG(x1), LOG(x1_0))
    ELSE
       i1 = LOCATE(x1, x1_0)
    END IF
    IF (log_abscissa2) THEN
       i2 = LOCATE(LOG(x2), LOG(x2_0))
    ELSE
       i2 = LOCATE(x2, x2_0)
    END IF

    IF (PRESENT(out_of_range_use_nearest_edge)) THEN
       outside_to_edge = out_of_range_use_nearest_edge
    ELSE
       outside_to_edge = .FALSE.
    END IF
    IF (.NOT.outside_to_edge .AND. &
        .NOT.extrapolate     .AND. &
         ( (i1 == 0) .OR. (i1 == SIZE(x1)) .OR. &
           (i2 == 0) .OR. (i2 == SIZE(x2)) ) ) THEN
       IF (PRESENT(error_status)) THEN
          error_status = signal_interpolate_interpolant_out_of_range
          y0 = HUGE(0.d0)
          RETURN
       ELSE
          WRITE(0,*) 'ERROR: bilinear_interpolate_d: interpolant outside input range'
          WRITE(0,*) x1_0, x1(1), x1(SIZE(x1))
          WRITE(0,*) x2_0, x2(1), x2(SIZE(x2))
          STOP
       END IF
    END IF

    IF (PRESENT(bad_flag)) THEN
       IF ( (y(i1  ,i2  ) <= bad_flag) .OR. &
            (y(i1+1,i2  ) <= bad_flag) .OR. &
            (y(i1  ,i2+1) <= bad_flag) .OR. &
            (y(i1+1,i2+1) <= bad_flag) ) THEN
          WRITE(0,*) 'ERROR: bilinear_interpolate_d: bad value(s):', &
               y(i1  ,i2  ), y(i1+1,i2  ), y(i1  ,i2+1), y(i1+1,i2+1)
          STOP
       END IF
    END IF

    d1 = distance(x1_0, x1, i1, log_abscissa1, outside_to_edge, extrapolate)
    d2 = distance(x2_0, x2, i2, log_abscissa2, outside_to_edge, extrapolate)

    y0 = (1.d0-d1)*(1.d0-d2)*y(i1  ,i2  ) + &
         (     d1)*(1.d0-d2)*y(i1+1,i2  ) + &
         (1.d0-d1)*(     d2)*y(i1  ,i2+1) + &
         (     d1)*(     d2)*y(i1+1,i2+1)

  END FUNCTION bilinear_interpolate_d



  !----------------------------------------------------------------------
  ! Trilinear interpolation
  !----------------------------------------------------------------------
  REAL          FUNCTION trilinear_interpolate_s(x1_0, x2_0, x3_0, x1, x2,     &
       x3, y, log_x1, log_x2, log_x3, bad_flag, out_of_range_use_nearest_edge, &
       linear_extrapolate,                                                     &
       error_status)                                                           &
       RESULT(y0)

    IMPLICIT NONE
    REAL, INTENT(IN), DIMENSION(:,:,:) :: y
    REAL, INTENT(IN), DIMENSION(:)     :: x1, x2, x3
    REAL, INTENT(IN)                   :: x1_0, x2_0, x3_0
    LOGICAL, INTENT(IN), OPTIONAL      :: log_x1, log_x2, log_x3
    REAL, INTENT(IN), OPTIONAL         :: bad_flag
    LOGICAL, INTENT(IN), OPTIONAL      :: out_of_range_use_nearest_edge
    LOGICAL, INTENT(IN), OPTIONAL      :: linear_extrapolate
    INTEGER, INTENT(OUT), OPTIONAL     :: error_status
    LOGICAL :: log_abscissa1, log_abscissa2, log_abscissa3
    LOGICAL :: outside_to_edge
    LOGICAL :: extrapolate
    INTEGER :: i1, i2, i3
    REAL :: d1, d2, d3

    IF (PRESENT(error_status)) error_status = 0

    IF (PRESENT(linear_extrapolate)) THEN
       extrapolate = linear_extrapolate
    ELSE
       extrapolate = .FALSE.
    END IF

    IF (PRESENT(log_x1)) THEN
       log_abscissa1 = log_x1
    ELSE
       log_abscissa1 = .FALSE.
    END IF
    IF (PRESENT(log_x2)) THEN
       log_abscissa2 = log_x2
    ELSE
       log_abscissa2 = .FALSE.
    END IF
    IF (PRESENT(log_x3)) THEN
       log_abscissa3 = log_x3
    ELSE
       log_abscissa3 = .FALSE.
    END IF

    IF (log_abscissa1) THEN
       i1 = LOCATE(LOG(x1), LOG(x1_0))
    ELSE
       i1 = LOCATE(x1, x1_0)
    END IF
    IF (log_abscissa2) THEN
       i2 = LOCATE(LOG(x2), LOG(x2_0))
    ELSE
       i2 = LOCATE(x2, x2_0)
    END IF
    IF (log_abscissa3) THEN
       i3 = LOCATE(LOG(x3), LOG(x3_0))
    ELSE
       i3 = LOCATE(x3, x3_0)
    END IF

    IF (PRESENT(out_of_range_use_nearest_edge)) THEN
       outside_to_edge = out_of_range_use_nearest_edge
    ELSE
       outside_to_edge = .FALSE.
    END IF
    IF (.NOT.outside_to_edge .AND. &
        .NOT.extrapolate     .AND. &
         ( (i1 == 0) .OR. (i1 == SIZE(x1)) .OR. &
           (i2 == 0) .OR. (i2 == SIZE(x2)) .OR. &
           (i3 == 0) .OR. (i3 == SIZE(x3)) ) ) THEN
       IF (PRESENT(error_status)) THEN
          error_status = signal_interpolate_interpolant_out_of_range
          y0 = HUGE(0.)
          RETURN
       ELSE
          WRITE(0,*) 'ERROR: trilinear_interpolate_s: interpolant outside input range'
          WRITE(0,*) x1_0, x1(1), x1(SIZE(x1))
          WRITE(0,*) x2_0, x2(1), x2(SIZE(x2))
          WRITE(0,*) x3_0, x3(1), x3(SIZE(x3))
          STOP
       END IF
    END IF

    IF (PRESENT(bad_flag)) THEN
       IF ( (y(i1  ,i2  ,i3  ) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3  ) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3  ) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3  ) <= bad_flag) .OR. &
            (y(i1  ,i2  ,i3+1) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3+1) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3+1) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3+1) <= bad_flag) ) THEN
          WRITE(0,*) 'ERROR: trilinear_interpolate_s: bad value(s):', &
               y(i1  ,i2  ,i3  ), &
               y(i1+1,i2  ,i3  ), &
               y(i1  ,i2+1,i3  ), &
               y(i1+1,i2+1,i3  ), &
               y(i1  ,i2  ,i3+1), &
               y(i1+1,i2  ,i3+1), &
               y(i1  ,i2+1,i3+1), &
               y(i1+1,i2+1,i3+1)
          STOP
       END IF
    END IF

    d1 = distance(x1_0, x1, i1, log_abscissa1, outside_to_edge, extrapolate)
    d2 = distance(x2_0, x2, i2, log_abscissa2, outside_to_edge, extrapolate)
    d3 = distance(x3_0, x3, i3, log_abscissa3, outside_to_edge, extrapolate)

    y0 = (1.-d1)*(1.-d2)*(1.-d3)*y(i1  ,i2  ,i3  ) + &
         (   d1)*(1.-d2)*(1.-d3)*y(i1+1,i2  ,i3  ) + &
         (1.-d1)*(   d2)*(1.-d3)*y(i1  ,i2+1,i3  ) + &
         (   d1)*(   d2)*(1.-d3)*y(i1+1,i2+1,i3  ) + &
         (1.-d1)*(1.-d2)*(   d3)*y(i1  ,i2  ,i3+1) + &
         (   d1)*(1.-d2)*(   d3)*y(i1+1,i2  ,i3+1) + &
         (1.-d1)*(   d2)*(   d3)*y(i1  ,i2+1,i3+1) + &
         (   d1)*(   d2)*(   d3)*y(i1+1,i2+1,i3+1)

  END FUNCTION trilinear_interpolate_s

  REAL(KIND=KIND(0.D0)) FUNCTION trilinear_interpolate_d(x1_0, x2_0, x3_0, x1, x2,     &
       x3, y, log_x1, log_x2, log_x3, bad_flag, out_of_range_use_nearest_edge, &
       linear_extrapolate,                                                     &
       error_status)                                                           &
       RESULT(y0)

    IMPLICIT NONE
    REAL(KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:,:,:) :: y
    REAL(KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:)     :: x1, x2, x3
    REAL(KIND=KIND(0.D0)), INTENT(IN)                   :: x1_0, x2_0, x3_0
    LOGICAL, INTENT(IN), OPTIONAL               :: log_x1, log_x2, log_x3
    REAL(KIND=KIND(0.D0)), INTENT(IN), OPTIONAL         :: bad_flag
    LOGICAL, INTENT(IN), OPTIONAL               :: out_of_range_use_nearest_edge
    LOGICAL, INTENT(IN), OPTIONAL               :: linear_extrapolate
    INTEGER, INTENT(OUT), OPTIONAL              :: error_status
    LOGICAL :: log_abscissa1, log_abscissa2, log_abscissa3
    LOGICAL :: outside_to_edge
    LOGICAL :: extrapolate
    INTEGER :: i1, i2, i3
    REAL(KIND=KIND(0.D0)) :: d1, d2, d3

    IF (PRESENT(error_status)) error_status = 0

    IF (PRESENT(linear_extrapolate)) THEN
       extrapolate = linear_extrapolate
    ELSE
       extrapolate = .FALSE.
    END IF

    IF (PRESENT(log_x1)) THEN
       log_abscissa1 = log_x1
    ELSE
       log_abscissa1 = .FALSE.
    END IF
    IF (PRESENT(log_x2)) THEN
       log_abscissa2 = log_x2
    ELSE
       log_abscissa2 = .FALSE.
    END IF
    IF (PRESENT(log_x3)) THEN
       log_abscissa3 = log_x3
    ELSE
       log_abscissa3 = .FALSE.
    END IF

    IF (log_abscissa1) THEN
       i1 = LOCATE(LOG(x1), LOG(x1_0))
    ELSE
       i1 = LOCATE(x1, x1_0)
    END IF
    IF (log_abscissa2) THEN
       i2 = LOCATE(LOG(x2), LOG(x2_0))
    ELSE
       i2 = LOCATE(x2, x2_0)
    END IF
    IF (log_abscissa3) THEN
       i3 = LOCATE(LOG(x3), LOG(x3_0))
    ELSE
       i3 = LOCATE(x3, x3_0)
    END IF

    IF (PRESENT(out_of_range_use_nearest_edge)) THEN
       outside_to_edge = out_of_range_use_nearest_edge
    ELSE
       outside_to_edge = .FALSE.
    END IF
    IF (.NOT.outside_to_edge .AND. &
        .NOT.extrapolate     .AND. &
         ( (i1 == 0) .OR. (i1 == SIZE(x1)) .OR. &
           (i2 == 0) .OR. (i2 == SIZE(x2)) .OR. &
           (i3 == 0) .OR. (i3 == SIZE(x3)) ) ) THEN
       IF (PRESENT(error_status)) THEN
          error_status = signal_interpolate_interpolant_out_of_range
          y0 = HUGE(0.d0)
          RETURN
       ELSE
          WRITE(0,*) 'ERROR: trilinear_interpolate_d: interpolant outside input range'
          WRITE(0,*) x1_0, x1(1), x1(SIZE(x1))
          WRITE(0,*) x2_0, x2(1), x2(SIZE(x2))
          WRITE(0,*) x3_0, x3(1), x3(SIZE(x3))
          STOP
       END IF
    END IF

    IF (PRESENT(bad_flag)) THEN
       IF ( (y(i1  ,i2  ,i3  ) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3  ) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3  ) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3  ) <= bad_flag) .OR. &
            (y(i1  ,i2  ,i3+1) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3+1) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3+1) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3+1) <= bad_flag) ) THEN
          WRITE(0,*) 'ERROR: trilinear_interpolate_d: bad value(s):', &
               y(i1  ,i2  ,i3  ), &
               y(i1+1,i2  ,i3  ), &
               y(i1  ,i2+1,i3  ), &
               y(i1+1,i2+1,i3  ), &
               y(i1  ,i2  ,i3+1), &
               y(i1+1,i2  ,i3+1), &
               y(i1  ,i2+1,i3+1), &
               y(i1+1,i2+1,i3+1)
          STOP
       END IF
    END IF

    d1 = distance(x1_0, x1, i1, log_abscissa1, outside_to_edge, extrapolate)
    d2 = distance(x2_0, x2, i2, log_abscissa2, outside_to_edge, extrapolate)
    d3 = distance(x3_0, x3, i3, log_abscissa3, outside_to_edge, extrapolate)

    y0 = (1.d0-d1)*(1.d0-d2)*(1.d0-d3)*y(i1  ,i2  ,i3  ) + &
         (     d1)*(1.d0-d2)*(1.d0-d3)*y(i1+1,i2  ,i3  ) + &
         (1.d0-d1)*(     d2)*(1.d0-d3)*y(i1  ,i2+1,i3  ) + &
         (     d1)*(     d2)*(1.d0-d3)*y(i1+1,i2+1,i3  ) + &
         (1.d0-d1)*(1.d0-d2)*(     d3)*y(i1  ,i2  ,i3+1) + &
         (     d1)*(1.d0-d2)*(     d3)*y(i1+1,i2  ,i3+1) + &
         (1.d0-d1)*(     d2)*(     d3)*y(i1  ,i2+1,i3+1) + &
         (     d1)*(     d2)*(     d3)*y(i1+1,i2+1,i3+1)

  END FUNCTION trilinear_interpolate_d



  !----------------------------------------------------------------------
  ! Quadrilinear interpolation
  !----------------------------------------------------------------------
  REAL          FUNCTION quadrilinear_interpolate_s(x1_0, x2_0, x3_0, x4_0,  &
       x1, x2, x3, x4, y, log_x1, log_x2, log_x3, log_x4, bad_flag,          &
       out_of_range_use_nearest_edge, linear_extrapolate,                    &
       error_status)                                                         &
       RESULT(y0)

    IMPLICIT NONE
    REAL, INTENT(IN), DIMENSION(:,:,:,:) :: y
    REAL, INTENT(IN), DIMENSION(:)       :: x1, x2, x3, x4
    REAL, INTENT(IN)                     :: x1_0, x2_0, x3_0, x4_0
    LOGICAL, INTENT(IN), OPTIONAL        :: log_x1, log_x2, log_x3, log_x4
    REAL, INTENT(IN), OPTIONAL           :: bad_flag
    LOGICAL, INTENT(IN), OPTIONAL        :: out_of_range_use_nearest_edge
    LOGICAL, INTENT(IN), OPTIONAL        :: linear_extrapolate
    INTEGER, INTENT(OUT), OPTIONAL       :: error_status
    LOGICAL :: log_abscissa1, log_abscissa2, log_abscissa3, log_abscissa4
    LOGICAL :: outside_to_edge
    LOGICAL :: extrapolate
    INTEGER :: i1, i2, i3, i4
    REAL :: d1, d2, d3, d4

    IF (PRESENT(error_status)) error_status = 0

    IF (PRESENT(linear_extrapolate)) THEN
       extrapolate = linear_extrapolate
    ELSE
       extrapolate = .FALSE.
    END IF

    IF (PRESENT(log_x1)) THEN
       log_abscissa1 = log_x1
    ELSE
       log_abscissa1 = .FALSE.
    END IF
    IF (PRESENT(log_x2)) THEN
       log_abscissa2 = log_x2
    ELSE
       log_abscissa2 = .FALSE.
    END IF
    IF (PRESENT(log_x3)) THEN
       log_abscissa3 = log_x3
    ELSE
       log_abscissa3 = .FALSE.
    END IF
    IF (PRESENT(log_x4)) THEN
       log_abscissa4 = log_x4
    ELSE
       log_abscissa4 = .FALSE.
    END IF

    IF (log_abscissa1) THEN
       i1 = LOCATE(LOG(x1), LOG(x1_0))
    ELSE
       i1 = LOCATE(x1, x1_0)
    END IF
    IF (log_abscissa2) THEN
       i2 = LOCATE(LOG(x2), LOG(x2_0))
    ELSE
       i2 = LOCATE(x2, x2_0)
    END IF
    IF (log_abscissa3) THEN
       i3 = LOCATE(LOG(x3), LOG(x3_0))
    ELSE
       i3 = LOCATE(x3, x3_0)
    END IF
    IF (log_abscissa4) THEN
       i4 = LOCATE(LOG(x4), LOG(x4_0))
    ELSE
       i4 = LOCATE(x4, x4_0)
    END IF

    IF (PRESENT(out_of_range_use_nearest_edge)) THEN
       outside_to_edge = out_of_range_use_nearest_edge
    ELSE
       outside_to_edge = .FALSE.
    END IF
    IF (.NOT.outside_to_edge .AND. &
        .NOT.extrapolate     .AND. &
         ( (i1 == 0) .OR. (i1 == SIZE(x1)) .OR. &
           (i2 == 0) .OR. (i2 == SIZE(x2)) .OR. &
           (i3 == 0) .OR. (i3 == SIZE(x3)) .OR. &
           (i4 == 0) .OR. (i4 == SIZE(x4)) ) ) THEN
       IF (PRESENT(error_status)) THEN
          error_status = signal_interpolate_interpolant_out_of_range
          y0 = HUGE(0.)
          RETURN
       ELSE
          WRITE(0,*) 'ERROR: quadrilinear_interpolate_s: interpolant outside input range'
          WRITE(0,*) x1_0, x1(1), x1(SIZE(x1))
          WRITE(0,*) x2_0, x2(1), x2(SIZE(x2))
          WRITE(0,*) x3_0, x3(1), x3(SIZE(x3))
          WRITE(0,*) x4_0, x4(1), x4(SIZE(x4))
          STOP
       END IF
    END IF

    IF (PRESENT(bad_flag)) THEN
       IF ( (y(i1  ,i2  ,i3  ,i4  ) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3  ,i4  ) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3  ,i4  ) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3  ,i4  ) <= bad_flag) .OR. &
            (y(i1  ,i2  ,i3+1,i4  ) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3+1,i4  ) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3+1,i4  ) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3+1,i4  ) <= bad_flag) .OR. &
            (y(i1  ,i2  ,i3  ,i4+1) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3  ,i4+1) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3  ,i4+1) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3  ,i4+1) <= bad_flag) .OR. &
            (y(i1  ,i2  ,i3+1,i4+1) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3+1,i4+1) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3+1,i4+1) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3+1,i4+1) <= bad_flag) ) THEN
          WRITE(0,*) 'ERROR: quadrilinear_interpolate_s: bad value(s):', &
               y(i1  ,i2  ,i3  ,i4  ), &
               y(i1+1,i2  ,i3  ,i4  ), &
               y(i1  ,i2+1,i3  ,i4  ), &
               y(i1+1,i2+1,i3  ,i4  ), &
               y(i1  ,i2  ,i3+1,i4  ), &
               y(i1+1,i2  ,i3+1,i4  ), &
               y(i1  ,i2+1,i3+1,i4  ), &
               y(i1+1,i2+1,i3+1,i4  ), &
               y(i1  ,i2  ,i3  ,i4+1), &
               y(i1+1,i2  ,i3  ,i4+1), &
               y(i1  ,i2+1,i3  ,i4+1), &
               y(i1+1,i2+1,i3  ,i4+1), &
               y(i1  ,i2  ,i3+1,i4+1), &
               y(i1+1,i2  ,i3+1,i4+1), &
               y(i1  ,i2+1,i3+1,i4+1), &
               y(i1+1,i2+1,i3+1,i4+1)
          STOP
       END IF
    END IF

    d1 = distance(x1_0, x1, i1, log_abscissa1, outside_to_edge, extrapolate)
    d2 = distance(x2_0, x2, i2, log_abscissa2, outside_to_edge, extrapolate)
    d3 = distance(x3_0, x3, i3, log_abscissa3, outside_to_edge, extrapolate)
    d4 = distance(x4_0, x4, i4, log_abscissa4, outside_to_edge, extrapolate)

    y0 = (1.-d1)*(1.-d2)*(1.-d3)*(1.-d4)*y(i1  ,i2  ,i3  ,i4  ) + &
         (   d1)*(1.-d2)*(1.-d3)*(1.-d4)*y(i1+1,i2  ,i3  ,i4  ) + &
         (1.-d1)*(   d2)*(1.-d3)*(1.-d4)*y(i1  ,i2+1,i3  ,i4  ) + &
         (   d1)*(   d2)*(1.-d3)*(1.-d4)*y(i1+1,i2+1,i3  ,i4  ) + &
         (1.-d1)*(1.-d2)*(   d3)*(1.-d4)*y(i1  ,i2  ,i3+1,i4  ) + &
         (   d1)*(1.-d2)*(   d3)*(1.-d4)*y(i1+1,i2  ,i3+1,i4  ) + &
         (1.-d1)*(   d2)*(   d3)*(1.-d4)*y(i1  ,i2+1,i3+1,i4  ) + &
         (   d1)*(   d2)*(   d3)*(1.-d4)*y(i1+1,i2+1,i3+1,i4  ) + &
         (1.-d1)*(1.-d2)*(1.-d3)*(   d4)*y(i1  ,i2  ,i3  ,i4+1) + &
         (   d1)*(1.-d2)*(1.-d3)*(   d4)*y(i1+1,i2  ,i3  ,i4+1) + &
         (1.-d1)*(   d2)*(1.-d3)*(   d4)*y(i1  ,i2+1,i3  ,i4+1) + &
         (   d1)*(   d2)*(1.-d3)*(   d4)*y(i1+1,i2+1,i3  ,i4+1) + &
         (1.-d1)*(1.-d2)*(   d3)*(   d4)*y(i1  ,i2  ,i3+1,i4+1) + &
         (   d1)*(1.-d2)*(   d3)*(   d4)*y(i1+1,i2  ,i3+1,i4+1) + &
         (1.-d1)*(   d2)*(   d3)*(   d4)*y(i1  ,i2+1,i3+1,i4+1) + &
         (   d1)*(   d2)*(   d3)*(   d4)*y(i1+1,i2+1,i3+1,i4+1)

  END FUNCTION quadrilinear_interpolate_s

  REAL(KIND=KIND(0.D0)) FUNCTION quadrilinear_interpolate_d(x1_0, x2_0, x3_0, x4_0,  &
       x1, x2, x3, x4, y, log_x1, log_x2, log_x3, log_x4, bad_flag,          &
       out_of_range_use_nearest_edge, linear_extrapolate,                    &
       error_status)                                                         &
       RESULT(y0)

    IMPLICIT NONE
    REAL(KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:,:,:,:) :: y
    REAL(KIND=KIND(0.D0)), INTENT(IN), DIMENSION(:)       :: x1, x2, x3, x4
    REAL(KIND=KIND(0.D0)), INTENT(IN)                     :: x1_0, x2_0, x3_0, x4_0
    LOGICAL, INTENT(IN), OPTIONAL                 :: log_x1, log_x2, log_x3, log_x4
    REAL(KIND=KIND(0.D0)), INTENT(IN), OPTIONAL           :: bad_flag
    LOGICAL, INTENT(IN), OPTIONAL                 :: out_of_range_use_nearest_edge
    LOGICAL, INTENT(IN), OPTIONAL                 :: linear_extrapolate
    INTEGER, INTENT(OUT), OPTIONAL                :: error_status
    LOGICAL :: log_abscissa1, log_abscissa2, log_abscissa3, log_abscissa4
    LOGICAL :: outside_to_edge
    LOGICAL :: extrapolate
    INTEGER :: i1, i2, i3, i4
    REAL(KIND=KIND(0.D0)) :: d1, d2, d3, d4

    IF (PRESENT(error_status)) error_status = 0

    IF (PRESENT(linear_extrapolate)) THEN
       extrapolate = linear_extrapolate
    ELSE
       extrapolate = .FALSE.
    END IF

    IF (PRESENT(log_x1)) THEN
       log_abscissa1 = log_x1
    ELSE
       log_abscissa1 = .FALSE.
    END IF
    IF (PRESENT(log_x2)) THEN
       log_abscissa2 = log_x2
    ELSE
       log_abscissa2 = .FALSE.
    END IF
    IF (PRESENT(log_x3)) THEN
       log_abscissa3 = log_x3
    ELSE
       log_abscissa3 = .FALSE.
    END IF
    IF (PRESENT(log_x4)) THEN
       log_abscissa4 = log_x4
    ELSE
       log_abscissa4 = .FALSE.
    END IF

    IF (log_abscissa1) THEN
       i1 = LOCATE(LOG(x1), LOG(x1_0))
    ELSE
       i1 = LOCATE(x1, x1_0)
    END IF
    IF (log_abscissa2) THEN
       i2 = LOCATE(LOG(x2), LOG(x2_0))
    ELSE
       i2 = LOCATE(x2, x2_0)
    END IF
    IF (log_abscissa3) THEN
       i3 = LOCATE(LOG(x3), LOG(x3_0))
    ELSE
       i3 = LOCATE(x3, x3_0)
    END IF
    IF (log_abscissa4) THEN
       i4 = LOCATE(LOG(x4), LOG(x4_0))
    ELSE
       i4 = LOCATE(x4, x4_0)
    END IF

    IF (PRESENT(out_of_range_use_nearest_edge)) THEN
       outside_to_edge = out_of_range_use_nearest_edge
    ELSE
       outside_to_edge = .FALSE.
    END IF
    IF (.NOT.outside_to_edge .AND. &
        .NOT.extrapolate     .AND. &
         ( (i1 == 0) .OR. (i1 == SIZE(x1)) .OR. &
           (i2 == 0) .OR. (i2 == SIZE(x2)) .OR. &
           (i3 == 0) .OR. (i3 == SIZE(x3)) .OR. &
           (i4 == 0) .OR. (i4 == SIZE(x4)) ) ) THEN
       IF (PRESENT(error_status)) THEN
          error_status = signal_interpolate_interpolant_out_of_range
          y0 = HUGE(0.d0)
          RETURN
       ELSE
          WRITE(0,*) 'ERROR: quadrilinear_interpolate_s: interpolant outside input range'
          WRITE(0,*) x1_0, x1(1), x1(SIZE(x1))
          WRITE(0,*) x2_0, x2(1), x2(SIZE(x2))
          WRITE(0,*) x3_0, x3(1), x3(SIZE(x3))
          WRITE(0,*) x4_0, x4(1), x4(SIZE(x4))
          STOP
       END IF
    END IF

    IF (PRESENT(bad_flag)) THEN
       IF ( (y(i1  ,i2  ,i3  ,i4  ) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3  ,i4  ) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3  ,i4  ) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3  ,i4  ) <= bad_flag) .OR. &
            (y(i1  ,i2  ,i3+1,i4  ) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3+1,i4  ) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3+1,i4  ) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3+1,i4  ) <= bad_flag) .OR. &
            (y(i1  ,i2  ,i3  ,i4+1) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3  ,i4+1) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3  ,i4+1) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3  ,i4+1) <= bad_flag) .OR. &
            (y(i1  ,i2  ,i3+1,i4+1) <= bad_flag) .OR. &
            (y(i1+1,i2  ,i3+1,i4+1) <= bad_flag) .OR. &
            (y(i1  ,i2+1,i3+1,i4+1) <= bad_flag) .OR. &
            (y(i1+1,i2+1,i3+1,i4+1) <= bad_flag) ) THEN
          WRITE(0,*) 'ERROR: quadrilinear_interpolate_s: bad value(s):', &
               y(i1  ,i2  ,i3  ,i4  ), &
               y(i1+1,i2  ,i3  ,i4  ), &
               y(i1  ,i2+1,i3  ,i4  ), &
               y(i1+1,i2+1,i3  ,i4  ), &
               y(i1  ,i2  ,i3+1,i4  ), &
               y(i1+1,i2  ,i3+1,i4  ), &
               y(i1  ,i2+1,i3+1,i4  ), &
               y(i1+1,i2+1,i3+1,i4  ), &
               y(i1  ,i2  ,i3  ,i4+1), &
               y(i1+1,i2  ,i3  ,i4+1), &
               y(i1  ,i2+1,i3  ,i4+1), &
               y(i1+1,i2+1,i3  ,i4+1), &
               y(i1  ,i2  ,i3+1,i4+1), &
               y(i1+1,i2  ,i3+1,i4+1), &
               y(i1  ,i2+1,i3+1,i4+1), &
               y(i1+1,i2+1,i3+1,i4+1)
          STOP
       END IF
    END IF

    d1 = distance(x1_0, x1, i1, log_abscissa1, outside_to_edge, extrapolate)
    d2 = distance(x2_0, x2, i2, log_abscissa2, outside_to_edge, extrapolate)
    d3 = distance(x3_0, x3, i3, log_abscissa3, outside_to_edge, extrapolate)
    d4 = distance(x4_0, x4, i4, log_abscissa4, outside_to_edge, extrapolate)

    y0 = (1.d0-d1)*(1.d0-d2)*(1.d0-d3)*(1.d0-d4)*y(i1  ,i2  ,i3  ,i4  ) + &
         (     d1)*(1.d0-d2)*(1.d0-d3)*(1.d0-d4)*y(i1+1,i2  ,i3  ,i4  ) + &
         (1.d0-d1)*(     d2)*(1.d0-d3)*(1.d0-d4)*y(i1  ,i2+1,i3  ,i4  ) + &
         (     d1)*(     d2)*(1.d0-d3)*(1.d0-d4)*y(i1+1,i2+1,i3  ,i4  ) + &
         (1.d0-d1)*(1.d0-d2)*(     d3)*(1.d0-d4)*y(i1  ,i2  ,i3+1,i4  ) + &
         (     d1)*(1.d0-d2)*(     d3)*(1.d0-d4)*y(i1+1,i2  ,i3+1,i4  ) + &
         (1.d0-d1)*(     d2)*(     d3)*(1.d0-d4)*y(i1  ,i2+1,i3+1,i4  ) + &
         (     d1)*(     d2)*(     d3)*(1.d0-d4)*y(i1+1,i2+1,i3+1,i4  ) + &
         (1.d0-d1)*(1.d0-d2)*(1.d0-d3)*(     d4)*y(i1  ,i2  ,i3  ,i4+1) + &
         (     d1)*(1.d0-d2)*(1.d0-d3)*(     d4)*y(i1+1,i2  ,i3  ,i4+1) + &
         (1.d0-d1)*(     d2)*(1.d0-d3)*(     d4)*y(i1  ,i2+1,i3  ,i4+1) + &
         (     d1)*(     d2)*(1.d0-d3)*(     d4)*y(i1+1,i2+1,i3  ,i4+1) + &
         (1.d0-d1)*(1.d0-d2)*(     d3)*(     d4)*y(i1  ,i2  ,i3+1,i4+1) + &
         (     d1)*(1.d0-d2)*(     d3)*(     d4)*y(i1+1,i2  ,i3+1,i4+1) + &
         (1.d0-d1)*(     d2)*(     d3)*(     d4)*y(i1  ,i2+1,i3+1,i4+1) + &
         (     d1)*(     d2)*(     d3)*(     d4)*y(i1+1,i2+1,i3+1,i4+1)

  END FUNCTION quadrilinear_interpolate_d

!-----------------------------------------------------------------------------
!
! FOUR1.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE four1_sp(data,isign)
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER(I4B) :: n,m1,m2,j
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_sp')
    m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
    m2=n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat=reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta=arth(0,isign,m1)*TWOPI_D/n
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2,m2
       w=w*wp+w
       dat(:,j)=dat(:,j)*w
    end do
    temp=transpose(dat)
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)
  END SUBROUTINE four1_sp

  SUBROUTINE four1_dp(data,isign)
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER(I4B) :: n,m1,m2,j
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_dp')
    m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
    m2=n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat=reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta=arth(0,isign,m1)*TWOPI_D/n
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2,m2
       w=w*wp+w
       dat(:,j)=dat(:,j)*w
    end do
    temp=transpose(dat)
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)
  END SUBROUTINE four1_dp

!-----------------------------------------------------------------------------
!
! FOURROW.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE fourrow_sp(data,isign)
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp
    COMPLEX(SPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
    n2=n/2
    j=n2
    do i=1,n-2
       if (j > i) call swap(data(:,j+1),data(:,i+1))
       m=n2
       do
          if (m < 2 .or. j < m) exit
          j=j-m
          m=m/2
       end do
       j=j+m
    end do
    mmax=1
    do
       if (n <= mmax) exit
       istep=2*mmax
       theta=PI_D/(isign*mmax)
       wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
       w=cmplx(1.0_dp,0.0_dp,kind=dpc)
       do m=1,mmax
          ws=w
          do i=m,n,istep
             j=i+mmax
             temp=ws*data(:,j)
             data(:,j)=data(:,i)-temp
             data(:,i)=data(:,i)+temp
          end do
          w=w*wp+w
       end do
       mmax=istep
    end do
  END SUBROUTINE fourrow_sp

  SUBROUTINE fourrow_dp(data,isign)
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(DPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp
    COMPLEX(DPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
    n2=n/2
    j=n2
    do i=1,n-2
       if (j > i) call swap(data(:,j+1),data(:,i+1))
       m=n2
       do
          if (m < 2 .or. j < m) exit
          j=j-m
          m=m/2
       end do
       j=j+m
    end do
    mmax=1
    do
       if (n <= mmax) exit
       istep=2*mmax
       theta=PI_D/(isign*mmax)
       wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
       w=cmplx(1.0_dp,0.0_dp,kind=dpc)
       do m=1,mmax
          ws=w
          do i=m,n,istep
             j=i+mmax
             temp=ws*data(:,j)
             data(:,j)=data(:,i)-temp
             data(:,i)=data(:,i)+temp
          end do
          w=w*wp+w
       end do
       mmax=istep
    end do
  END SUBROUTINE fourrow_dp

!-----------------------------------------------------------------------------
!
! SORT2.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE sort2(arr,slave)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
    INTEGER(I4B) :: ndum
    INTEGER(I4B), DIMENSION(size(arr)) :: index
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexx(arr,index)
    arr=arr(index)
    slave=slave(index)
  END SUBROUTINE sort2

!-----------------------------------------------------------------------------
!
! INDEXX.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE indexx_sp(arr,index)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
    INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
    REAL(SP) :: a
    INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
    INTEGER(I4B), DIMENSION(NSTACK) :: istack
    n=assert_eq(size(index),size(arr),'indexx_sp')
    index=arth(1,1,n)
    jstack=0
    l=1
    r=n
    do
       if (r-l < NN) then
          do j=l+1,r
             indext=index(j)
             a=arr(indext)
             do i=j-1,l,-1
                if (arr(index(i)) <= a) exit
                index(i+1)=index(i)
             end do
             index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
       else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r),arr(index(l)),arr(index(r)))
          call icomp_xchg(index(l+1),index(r),arr(index(l+1)),arr(index(r)))
          call icomp_xchg(index(l),index(l+1),arr(index(l)),arr(index(l+1)))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
             do
                i=i+1
                if (arr(index(i)) >= a) exit
             end do
             do
                j=j-1
                if (arr(index(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
          if (r-i+1 >= j-l) then
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
          else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
          end if
       end if
    end do
  END SUBROUTINE indexx_sp

  SUBROUTINE indexx_dp(arr,index)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
    INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
    REAL(DP) :: a
    INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
    INTEGER(I4B), DIMENSION(NSTACK) :: istack
    n=assert_eq(size(index),size(arr),'indexx_sp')
    index=arth(1,1,n)
    jstack=0
    l=1
    r=n
    do
       if (r-l < NN) then
          do j=l+1,r
             indext=index(j)
             a=arr(indext)
             do i=j-1,l,-1
                if (arr(index(i)) <= a) exit
                index(i+1)=index(i)
             end do
             index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
       else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r),arr(index(l)),arr(index(r)))
          call icomp_xchg(index(l+1),index(r),arr(index(l+1)),arr(index(r)))
          call icomp_xchg(index(l),index(l+1),arr(index(l)),arr(index(l+1)))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
             do
                i=i+1
                if (arr(index(i)) >= a) exit
             end do
             do
                j=j-1
                if (arr(index(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
          if (r-i+1 >= j-l) then
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
          else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
          end if
       end if
    end do
  END SUBROUTINE indexx_dp

  SUBROUTINE indexx_i4b(iarr,index)
    IMPLICIT NONE
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
    INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
    INTEGER(I4B) :: a
    INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
    INTEGER(I4B), DIMENSION(NSTACK) :: istack
    n=assert_eq(size(index),size(iarr),'indexx_sp')
    index=arth(1,1,n)
    jstack=0
    l=1
    r=n
    do
       if (r-l < NN) then
          do j=l+1,r
             indext=index(j)
             a=iarr(indext)
             do i=j-1,l,-1
                if (iarr(index(i)) <= a) exit
                index(i+1)=index(i)
             end do
             index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
       else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r),iarr(index(l)),iarr(index(r)))
          call icomp_xchg(index(l+1),index(r),iarr(index(l+1)),iarr(index(r)))
          call icomp_xchg(index(l),index(l+1),iarr(index(l)),iarr(index(l+1)))
          i=l+1
          j=r
          indext=index(l+1)
          a=iarr(indext)
          do
             do
                i=i+1
                if (iarr(index(i)) >= a) exit
             end do
             do
                j=j-1
                if (iarr(index(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
          if (r-i+1 >= j-l) then
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
          else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
          end if
       end if
    end do
  END SUBROUTINE indexx_i4b

!-----------------------------------------------------------------------------
! icomp_xchg
!-----------------------------------------------------------------------------

  SUBROUTINE icomp_xchg_sp(i,j,ai,aj)
    REAL(SP), INTENT(IN) :: ai, aj
    INTEGER(I4B), INTENT(INOUT) :: i,j
    INTEGER(I4B) :: swp
    if (aj < ai) then
       swp=i
       i=j
       j=swp
    end if
  END SUBROUTINE icomp_xchg_sp

  SUBROUTINE icomp_xchg_dp(i,j,ai,aj)
    REAL(DP), INTENT(IN) :: ai, aj
    INTEGER(I4B), INTENT(INOUT) :: i,j
    INTEGER(I4B) :: swp
    if (aj < ai) then
       swp=i
       i=j
       j=swp
    end if
  END SUBROUTINE icomp_xchg_dp

  SUBROUTINE icomp_xchg_i4b(i,j,ai,aj)
    INTEGER(I4B), INTENT(IN) :: ai, aj
    INTEGER(I4B), INTENT(INOUT) :: i,j
    INTEGER(I4B) :: swp
    if (aj < ai) then
       swp=i
       i=j
       j=swp
    end if
  END SUBROUTINE icomp_xchg_i4b

!-----------------------------------------------------------------------------
!
! TRIDAG.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE tridag_r(a,b,c,r,u)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(SP), DIMENSION(:), INTENT(OUT) :: u
    REAL(SP), DIMENSION(size(b)) :: gam
    INTEGER(I4B) :: n,j
    REAL(SP) :: bet
    n=ASSERT_EQ((/SIZE(a)+1,SIZE(b),SIZE(c)+1,SIZE(r),SIZE(u)/),'tridag_ser')
    bet=b(1)
    IF (bet == 0.0) CALL NRERROR('tridag_ser: Error at code stage 1')
    u(1)=r(1)/bet
    DO j=2,n
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j-1)*gam(j)
       IF (bet == 0.0) CALL NRERROR('tridag_ser: Error at code stage 2')
       u(j)=(r(j)-a(j-1)*u(j-1))/bet
    END DO
    DO j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    END DO
  END SUBROUTINE tridag_r

  SUBROUTINE tridag_d(a,b,c,r,u)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(DP), DIMENSION(:), INTENT(OUT) :: u
    REAL(DP), DIMENSION(size(b)) :: gam
    INTEGER(I4B) :: n,j
    REAL(DP) :: bet
    n=ASSERT_EQ((/SIZE(a)+1,SIZE(b),SIZE(c)+1,SIZE(r),SIZE(u)/),'tridag_ser')
    bet=b(1)
    IF (bet == 0.d0) CALL NRERROR('tridag_ser: Error at code stage 1')
    u(1)=r(1)/bet
    DO j=2,n
       gam(j)=c(j-1)/bet
       bet=b(j)-a(j-1)*gam(j)
       IF (bet == 0.d0) CALL NRERROR('tridag_ser: Error at code stage 2')
       u(j)=(r(j)-a(j-1)*u(j-1))/bet
    END DO
    DO j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    END DO
  END SUBROUTINE tridag_d

!-----------------------------------------------------------------------------
!
! MATINV.F (not in package, but in book)
!
!-----------------------------------------------------------------------------
  FUNCTION matinv_r(a) RESULT(ainv)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
    REAL(SP), DIMENSION(SIZE(a,1),SIZE(a,2)) :: ainv, acopy
    INTEGER(I4B), DIMENSION(SIZE(a,1)) :: indx
    REAL(SP) :: d
    INTEGER :: i, n
    n = ASSERT_EQ((/SIZE(a,1),SIZE(a,2)/),'matinv_r')
    ainv = 0.
    DO i = 1, n
       ainv(i,i) = 1.
    END DO
    acopy = a
    CALL ludcmp(acopy,indx,d)
    DO i = 1, n
       CALL lubksb(acopy,indx,ainv(:,i))
    END DO
  END FUNCTION matinv_r

  FUNCTION matinv_d(a) RESULT(ainv)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
    REAL(DP), DIMENSION(SIZE(a,1),SIZE(a,2)) :: ainv, acopy
    INTEGER(I4B), DIMENSION(SIZE(a,1)) :: indx
    REAL(DP) :: d
    INTEGER :: i, n
    n = ASSERT_EQ((/SIZE(a,1),SIZE(a,2)/),'matinv_d')
    ainv = 0.
    DO i = 1, n
       ainv(i,i) = 1.
    END DO
    acopy = a
    CALL ludcmp(acopy,indx,d)
    DO i = 1, n
       CALL lubksb(acopy,indx,ainv(:,i))
    END DO
  END FUNCTION matinv_d

!-----------------------------------------------------------------------------
!
! MPROVE.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE mprove(a,alud,indx,b,x)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,alud
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(SP), DIMENSION(:), INTENT(IN) :: b
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
    INTEGER(I4B) :: ndum
    REAL(SP), DIMENSION(SIZE(a,1)) :: r
    ndum=ASSERT_EQ((/SIZE(a,1),SIZE(a,2),SIZE(alud,1),SIZE(alud,2),SIZE(b),&
                     SIZE(x),SIZE(indx)/),'mprove')
    r=MATMUL(REAL(a,dp),REAL(x,dp))-REAL(b,dp)
    CALL LUBKSB(alud,indx,r)
    x=x-r
  END SUBROUTINE mprove

!-----------------------------------------------------------------------------
!
! LUBKSB.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE lubksb_r(a,indx,b)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER(I4B) :: i,n,ii,ll
    REAL(SP) :: summ
    n=ASSERT_EQ(size(a,1),size(a,2),size(indx),'lubksb_r')
    ii=0
    DO i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       IF (ii /= 0) THEN
          summ=summ-DOT_PRODUCT(a(i,ii:i-1),b(ii:i-1))
       ELSE IF (summ /= 0.0) THEN
          ii=i
       END IF
       b(i)=summ
    END DO
    DO i=n,1,-1
       b(i) = (b(i)-DOT_PRODUCT(a(i,i+1:n),b(i+1:n)))/a(i,i)
    END DO
  END SUBROUTINE lubksb_r

  SUBROUTINE lubksb_d(a,indx,b)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER(I4B) :: i,n,ii,ll
    REAL(DP) :: summ
    n=ASSERT_EQ(size(a,1),size(a,2),size(indx),'lubksb_d')
    ii=0
    DO i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       IF (ii /= 0) THEN
          summ=summ-DOT_PRODUCT(a(i,ii:i-1),b(ii:i-1))
       ELSE IF (summ /= 0.0) THEN
          ii=i
       END IF
       b(i)=summ
    END DO
    DO i=n,1,-1
       b(i) = (b(i)-DOT_PRODUCT(a(i,i+1:n),b(i+1:n)))/a(i,i)
    END DO
  END SUBROUTINE lubksb_d

!-----------------------------------------------------------------------------
!
! LUDCMP.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE ludcmp_r(a,indx,d)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
    REAL(SP), INTENT(OUT) :: d
    REAL(SP), DIMENSION(size(a,1)) :: vv
    REAL(SP), PARAMETER :: TINY=1.0e-20_sp
    INTEGER(I4B) :: j,n,imax
    n=ASSERT_EQ(SIZE(a,1),SIZE(a,2),SIZE(indx),'ludcmp_r')
    d=1.0
    vv=MAXVAL(ABS(a),dim=2)
    IF (ANY(vv == 0.0)) CALL NRERROR('singular matrix in ludcmp_r')
    vv=1.0_sp/vv
    DO j=1,n
       imax=(j-1)+IMAXLOC(vv(j:n)*ABS(a(j:n,j)))
       IF (j /= imax) THEN
          CALL SWAP(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       END IF
       indx(j)=imax
       IF (a(j,j) == 0.0) a(j,j)=TINY
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-OUTERPROD(a(j+1:n,j),a(j,j+1:n))
    END DO
  END SUBROUTINE ludcmp_r

  SUBROUTINE ludcmp_d(a,indx,d)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
    REAL(DP), INTENT(OUT) :: d
    REAL(DP), DIMENSION(size(a,1)) :: vv
    REAL(DP), PARAMETER :: TINY=1.0d-20
    INTEGER(I4B) :: j,n,imax
    n=ASSERT_EQ(SIZE(a,1),SIZE(a,2),SIZE(indx),'ludcmp_d')
    d=1.0
    vv=MAXVAL(ABS(a),dim=2)
    IF (ANY(vv == 0.0)) CALL NRERROR('singular matrix in ludcmp_d')
    vv=1.0_dp/vv
    DO j=1,n
       imax=(j-1)+IMAXLOC(vv(j:n)*ABS(a(j:n,j)))
       IF (j /= imax) THEN
          CALL SWAP(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       END IF
       indx(j)=imax
       IF (a(j,j) == 0.0) a(j,j)=TINY
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-OUTERPROD(a(j+1:n,j),a(j,j+1:n))
    END DO
  END SUBROUTINE ludcmp_d

!-----------------------------------------------------------------------------
!
! TRED2.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE tred2_sp(a,d,e,novectors)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
    LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
    INTEGER(I4B) :: i,j,l,n
    REAL(SP) :: f,g,h,hh,scale
    REAL(SP), DIMENSION(size(a,1)) :: gg
    LOGICAL(LGT) :: yesvec
    n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2')
    if (present(novectors)) then
       yesvec=.not. novectors
    else
       yesvec=.true.
    end if
    do i=n,2,-1
       l=i-1
       h=0.0
       if (l > 1) then
          scale=sum(abs(a(i,1:l)))
          if (scale == 0.0) then
             e(i)=a(i,l)
          else
             a(i,1:l)=a(i,1:l)/scale
             h=sum(a(i,1:l)**2)
             f=a(i,l)
             g=-sign(sqrt(h),f)
             e(i)=scale*g
             h=h-f*g
             a(i,l)=f-g
             if (yesvec) a(1:l,i)=a(i,1:l)/h
             do j=1,l
                e(j)=(dot_product(a(j,1:j),a(i,1:j)) &
                     +dot_product(a(j+1:l,j),a(i,j+1:l)))/h
             end do
             f=dot_product(e(1:l),a(i,1:l))
             hh=f/(h+h)
             e(1:l)=e(1:l)-hh*a(i,1:l)
             do j=1,l
                a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
             end do
          end if
       else
          e(i)=a(i,l)
       end if
       d(i)=h
    end do
    if (yesvec) d(1)=0.0
    e(1)=0.0
    do i=1,n
       if (yesvec) then
          l=i-1
          if (d(i) /= 0.0) then
             gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
             a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
          end if
          d(i)=a(i,i)
          a(i,i)=1.0
          a(i,1:l)=0.0
          a(1:l,i)=0.0
       else
          d(i)=a(i,i)
       end if
    end do
  END SUBROUTINE tred2_sp

  SUBROUTINE tred2_dp(a,d,e,novectors)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(DP), DIMENSION(:), INTENT(OUT) :: d,e
    LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
    INTEGER(I4B) :: i,j,l,n
    REAL(DP) :: f,g,h,hh,scale
    REAL(DP), DIMENSION(size(a,1)) :: gg
    LOGICAL(LGT) :: yesvec
    n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2')
    if (present(novectors)) then
       yesvec=.not. novectors
    else
       yesvec=.true.
    end if
    do i=n,2,-1
       l=i-1
       h=0.0
       if (l > 1) then
          scale=sum(abs(a(i,1:l)))
          if (scale == 0.0) then
             e(i)=a(i,l)
          else
             a(i,1:l)=a(i,1:l)/scale
             h=sum(a(i,1:l)**2)
             f=a(i,l)
             g=-sign(sqrt(h),f)
             e(i)=scale*g
             h=h-f*g
             a(i,l)=f-g
             if (yesvec) a(1:l,i)=a(i,1:l)/h
             do j=1,l
                e(j)=(dot_product(a(j,1:j),a(i,1:j)) &
                     +dot_product(a(j+1:l,j),a(i,j+1:l)))/h
             end do
             f=dot_product(e(1:l),a(i,1:l))
             hh=f/(h+h)
             e(1:l)=e(1:l)-hh*a(i,1:l)
             do j=1,l
                a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
             end do
          end if
       else
          e(i)=a(i,l)
       end if
       d(i)=h
    end do
    if (yesvec) d(1)=0.0
    e(1)=0.0
    do i=1,n
       if (yesvec) then
          l=i-1
          if (d(i) /= 0.0) then
             gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
             a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
          end if
          d(i)=a(i,i)
          a(i,i)=1.0
          a(i,1:l)=0.0
          a(1:l,i)=0.0
       else
          d(i)=a(i,i)
       end if
    end do
  END SUBROUTINE tred2_dp

!-----------------------------------------------------------------------------
!
! TQLI.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE tqli_sp(d,e,z)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
    REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
    INTEGER(I4B) :: i,iter,l,m,n,ndum
    REAL(SP) :: b,c,dd,f,g,p,r,s
    REAL(SP), DIMENSION(size(e)) :: ff
    n=assert_eq(size(d),size(e),'tqli: n')
    if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
    e(:)=eoshift(e(:),1)
    do l=1,n
       iter=0
       iterate: do
          do m=l,n-1
             dd=abs(d(m))+abs(d(m+1))
             if (abs(e(m))+dd == dd) exit
          end do
          if (m == l) exit iterate
          if (iter == 30) call nrerror('too many iterations in tqli')
          iter=iter+1
          g=(d(l+1)-d(l))/(2.0_sp*e(l))
          r=pythag(g,1.0_sp)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0
          c=1.0
          p=0.0
          do i=m-1,l,-1
             f=s*e(i)
             b=c*e(i)
             r=pythag(f,g)
             e(i+1)=r
             if (r == 0.0) then
                d(i+1)=d(i+1)-p
                e(m)=0.0
                cycle iterate
             end if
             s=f/r
             c=g/r
             g=d(i+1)-p
             r=(d(i)-g)*s+2.0_sp*c*b
             p=s*r
             d(i+1)=g+p
             g=c*r-b
             if (present(z)) then
                ff(1:n)=z(1:n,i+1)
                z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
                z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
             end if
          end do
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0
       end do iterate
    end do
  END SUBROUTINE tqli_sp

  SUBROUTINE tqli_dp(d,e,z)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: d,e
    REAL(DP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
    INTEGER(I4B) :: i,iter,l,m,n,ndum
    REAL(DP) :: b,c,dd,f,g,p,r,s
    REAL(DP), DIMENSION(size(e)) :: ff
    n=assert_eq(size(d),size(e),'tqli: n')
    if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
    e(:)=eoshift(e(:),1)
    do l=1,n
       iter=0
       iterate: do
          do m=l,n-1
             dd=abs(d(m))+abs(d(m+1))
             if (abs(e(m))+dd == dd) exit
          end do
          if (m == l) exit iterate
          if (iter == 30) call nrerror('too many iterations in tqli')
          iter=iter+1
          g=(d(l+1)-d(l))/(2.0_dp*e(l))
          r=pythag(g,1.0_dp)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0
          c=1.0
          p=0.0
          do i=m-1,l,-1
             f=s*e(i)
             b=c*e(i)
             r=pythag(f,g)
             e(i+1)=r
             if (r == 0.0) then
                d(i+1)=d(i+1)-p
                e(m)=0.0
                cycle iterate
             end if
             s=f/r
             c=g/r
             g=d(i+1)-p
             r=(d(i)-g)*s+2.0_dp*c*b
             p=s*r
             d(i+1)=g+p
             g=c*r-b
             if (present(z)) then
                ff(1:n)=z(1:n,i+1)
                z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
                z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
             end if
          end do
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0
       end do iterate
    end do
  END SUBROUTINE tqli_dp

!-----------------------------------------------------------------------------
!
! PYTHAG.F
!
!-----------------------------------------------------------------------------

  FUNCTION pythag_sp(a,b)
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: a,b
    REAL(SP) :: pythag_sp
    REAL(SP) :: absa,absb
    absa=abs(a)
    absb=abs(b)
    if (absa > absb) then
       pythag_sp=absa*sqrt(1.0_sp+(absb/absa)**2)
    else
       if (absb == 0.0) then
          pythag_sp=0.0
       else
          pythag_sp=absb*sqrt(1.0_sp+(absa/absb)**2)
       end if
    end if
  END FUNCTION pythag_sp

  FUNCTION pythag_dp(a,b)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: a,b
    REAL(DP) :: pythag_dp
    REAL(DP) :: absa,absb
    absa=abs(a)
    absb=abs(b)
    if (absa > absb) then
       pythag_dp=absa*sqrt(1.0_dp+(absb/absa)**2)
    else
       if (absb == 0.0) then
          pythag_dp=0.0
       else
          pythag_dp=absb*sqrt(1.0_dp+(absa/absb)**2)
       end if
    end if
  END FUNCTION pythag_dp

!-----------------------------------------------------------------------------
!
! EIGSRT.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE eigsrt_sp(d,v)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: v
    INTEGER(I4B) :: i,j,n
    n=assert_eq(size(d),size(v,1),size(v,2),'eigsrt')
    do i=1,n-1
       j=imaxloc(d(i:n))+i-1
       if (j /= i) then
          call swap(d(i),d(j))
          call swap(v(:,i),v(:,j))
       end if
    end do
  END SUBROUTINE eigsrt_sp

  SUBROUTINE eigsrt_dp(d,v)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: d
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: v
    INTEGER(I4B) :: i,j,n
    n=assert_eq(size(d),size(v,1),size(v,2),'eigsrt')
    do i=1,n-1
       j=imaxloc(d(i:n))+i-1
       if (j /= i) then
          call swap(d(i),d(j))
          call swap(v(:,i),v(:,j))
       end if
    end do
  END SUBROUTINE eigsrt_dp

!-----------------------------------------------------------------------------
!
! BALANC.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE balanc(a)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(SP), PARAMETER :: RADX=radix(a),SQRADX=RADX**2
    INTEGER(I4B) :: i,last,ndum
    REAL(SP) :: c,f,g,r,s
    ndum=assert_eq(size(a,1),size(a,2),'balanc')
    do
       last=1
       do i=1,size(a,1)
          c=sum(abs(a(:,i)))-a(i,i)
          r=sum(abs(a(i,:)))-a(i,i)
          if (c /= 0.0 .and. r /= 0.0) then
             g=r/RADX
             f=1.0
             s=c+r
             do
                if (c >= g) exit
                f=f*RADX
                c=c*SQRADX
             end do
             g=r*RADX
             do
                if (c <= g) exit
                f=f/RADX
                c=c/SQRADX
             end do
             if ((c+r)/f < 0.95_sp*s) then
                last=0
                g=1.0_sp/f
                a(i,:)=a(i,:)*g
                a(:,i)=a(:,i)*f
             end if
          end if
       end do
       if (last /= 0) exit
    end do
  END SUBROUTINE balanc

!-----------------------------------------------------------------------------
!
! ELMHES.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE elmhes(a)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER(I4B) :: i,m,n
    REAL(SP) :: x
    REAL(SP), DIMENSION(size(a,1)) :: y
    n=assert_eq(size(a,1),size(a,2),'elmhes')
    do m=2,n-1
       i=imaxloc(abs(a(m:n,m-1)))+m-1
       x=a(i,m-1)
       if (i /= m) then
          call swap(a(i,m-1:n),a(m,m-1:n))
          call swap(a(:,i),a(:,m))
       end if
       if (x /= 0.0) then
          y(m+1:n)=a(m+1:n,m-1)/x
          a(m+1:n,m-1)=y(m+1:n)
          a(m+1:n,m:n)=a(m+1:n,m:n)-outerprod(y(m+1:n),a(m,m:n))
          a(:,m)=a(:,m)+matmul(a(:,m+1:n),y(m+1:n))
       end if
    end do
  END SUBROUTINE elmhes

!-----------------------------------------------------------------------------
!
! HQR.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE hqr(a,wr,wi)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(OUT) :: wr,wi
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER(I4B) :: i,its,k,l,m,n,nn,mnnk
    REAL(SP) :: anorm,p,q,r,s,t,u,v,w,x,y,z
    REAL(SP), DIMENSION(size(a,1)) :: pp
    n=assert_eq(size(a,1),size(a,2),size(wr),size(wi),'hqr')
    anorm=sum(abs(a),mask=upper_triangle(n,n,extra=2))
    nn=n
    t=0.0
    do
       if (nn < 1) exit
       its=0
       iterate: do
          small: do l=nn,2,-1
             s=abs(a(l-1,l-1))+abs(a(l,l))
             if (s == 0.0) s=anorm
             if (abs(a(l,l-1))+s == s) then
                a(l,l-1)=0.0
                exit small
             end if
          end do small
          x=a(nn,nn)
          if (l == nn) then
             wr(nn)=x+t
             wi(nn)=0.0
             nn=nn-1
             exit iterate
          end if
          y=a(nn-1,nn-1)
          w=a(nn,nn-1)*a(nn-1,nn)
          if (l == nn-1) then
             p=0.5_sp*(y-x)
             q=p**2+w
             z=sqrt(abs(q))
             x=x+t
             if (q >= 0.0) then
                z=p+sign(z,p)
                wr(nn)=x+z
                wr(nn-1)=wr(nn)
                if (z /= 0.0) wr(nn)=x-w/z
                wi(nn)=0.0
                wi(nn-1)=0.0
             else
                wr(nn)=x+p
                wr(nn-1)=wr(nn)
                wi(nn)=z
                wi(nn-1)=-z
             end if
             nn=nn-2
             exit iterate
          end if
          if (its == 30) call nrerror('too many iterations in hqr')
          if (its == 10 .or. its == 20) then
             t=t+x
             call diagadd(a(1:nn,1:nn),-x)
             s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
             x=0.75_sp*s
             y=x
             w=-0.4375_sp*s**2
          end if
          its=its+1
          do m=nn-2,l,-1
             z=a(m,m)
             r=x-z
             s=y-z
             p=(r*s-w)/a(m+1,m)+a(m,m+1)
             q=a(m+1,m+1)-z-r-s
             r=a(m+2,m+1)
             s=abs(p)+abs(q)+abs(r)
             p=p/s
             q=q/s
             r=r/s
             if (m == l) exit
             u=abs(a(m,m-1))*(abs(q)+abs(r))
             v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
             if (u+v == v) exit
          end do
          do i=m+2,nn
             a(i,i-2)=0.0
             if (i /= m+2) a(i,i-3)=0.0
          end do
          do k=m,nn-1
             if (k /= m) then
                p=a(k,k-1)
                q=a(k+1,k-1)
                r=0.0
                if (k /= nn-1) r=a(k+2,k-1)
                x=abs(p)+abs(q)+abs(r)
                if (x /= 0.0) then
                   p=p/x
                   q=q/x
                   r=r/x
                end if
             end if
             s=sign(sqrt(p**2+q**2+r**2),p)
             if (s /= 0.0) then
                if (k == m) then
                   if (l /= m) a(k,k-1)=-a(k,k-1)
                else
                   a(k,k-1)=-s*x
                end if
                p=p+s
                x=p/s
                y=q/s
                z=r/s
                q=q/p
                r=r/p
                pp(k:nn)=a(k,k:nn)+q*a(k+1,k:nn)
                if (k /= nn-1) then
                   pp(k:nn)=pp(k:nn)+r*a(k+2,k:nn)
                   a(k+2,k:nn)=a(k+2,k:nn)-pp(k:nn)*z
                end if
                a(k+1,k:nn)=a(k+1,k:nn)-pp(k:nn)*y
                a(k,k:nn)=a(k,k:nn)-pp(k:nn)*x
                mnnk=min(nn,k+3)
                pp(l:mnnk)=x*a(l:mnnk,k)+y*a(l:mnnk,k+1)
                if (k /= nn-1) then
                   pp(l:mnnk)=pp(l:mnnk)+z*a(l:mnnk,k+2)
                   a(l:mnnk,k+2)=a(l:mnnk,k+2)-pp(l:mnnk)*r
                end if
                a(l:mnnk,k+1)=a(l:mnnk,k+1)-pp(l:mnnk)*q
                a(l:mnnk,k)=a(l:mnnk,k)-pp(l:mnnk)
             end if
          end do
       end do iterate
    end do
  END SUBROUTINE hqr

!-----------------------------------------------------------------------------
!
! PERIOD.F (with modifications for phase information)
!
!-----------------------------------------------------------------------------

  SUBROUTINE period(x,y,ofac,hifac,px,py,jmax,prob,prl,pim)
    IMPLICIT NONE
    INTEGER(I4B), INTENT(OUT) :: jmax
    REAL(SP), INTENT(IN) :: ofac,hifac
    REAL(SP), INTENT(OUT) :: prob
    REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
    REAL(SP), DIMENSION(:), POINTER :: px,py
    ! ADT changes start here
    REAL(SP), DIMENSION(:), POINTER :: prl,pim
    REAL(SP), DIMENSION(:), ALLOCATABLE :: ph1
    INTEGER(I4B) :: nlen
    REAL(SP) :: ry,iy,phLS,arg0,arg1,amp
    ! End ADT changes
    INTEGER(I4B) :: i,n,nout
    REAL(SP) :: ave,cwtau,effm,expy,pnow,sumc,sumcy,&
            sums,sumsh,sumsy,swtau,var,wtau,xave,xdif,xmax,xmin
    REAL(DP), DIMENSION(size(x)) :: tmp1,tmp2,wi,wpi,wpr,wr
    LOGICAL(LGT), SAVE :: init=.true.
    n=assert_eq(size(x),size(y),'period')
    if (init) then
       init=.false.
       nullify(px,py)
       ! ADT changes start here
       nullify(prl,pim)
       ! End ADT changes
    else
       if (associated(px)) deallocate(px)
       if (associated(py)) deallocate(py)
       ! ADT changes start here
       if (associated(prl)) deallocate(prl)
       if (associated(pim)) deallocate(pim)
       ! End ADT changes
    end if
    nout=0.5_sp*ofac*hifac*n
    allocate(px(nout),py(nout))
    ! ADT changes start here
    allocate(prl(nout),pim(nout))
    allocate(ph1(nout))
    ! End ADT changes
    call avevar(y(:),ave,var)
    ! ADT changes start here
    ! Handle this case below ... we are probably trying to do a transform
    ! of a flat function, so just set all py's to 0 below
    !if (var == 0.0) call nrerror('zero variance in period')
    ! End ADT changes
    xmax=maxval(x(:))
    xmin=minval(x(:))
    xdif=xmax-xmin
    xave=0.5_sp*(xmax+xmin)
    pnow=1.0_sp/(xdif*ofac)
    tmp1(:)=TWOPI_D*((x(:)-xave)*pnow)
    wpr(:)=-2.0_dp*sin(0.5_dp*tmp1)**2
    wpi(:)=sin(tmp1(:))
    wr(:)=cos(tmp1(:))
    wi(:)=wpi(:)
    do i=1,nout
       px(i)=pnow
       sumsh=dot_product(wi,wr)
       sumc=dot_product(wr(:)-wi(:),wr(:)+wi(:))
       wtau=0.5_sp*atan2(2.0_sp*sumsh,sumc)
       swtau=sin(wtau)
       cwtau=cos(wtau)
       tmp1(:)=wi(:)*cwtau-wr(:)*swtau
       tmp2(:)=wr(:)*cwtau+wi(:)*swtau
       sums=dot_product(tmp1,tmp1)
       sumc=dot_product(tmp2,tmp2)
       sumsy=dot_product(y(:)-ave,tmp1)
       sumcy=dot_product(y(:)-ave,tmp2)
       tmp1(:)=wr(:)
       wr(:)=(wr(:)*wpr(:)-wi(:)*wpi(:))+wr(:)
       wi(:)=(wi(:)*wpr(:)+tmp1(:)*wpi(:))+wi(:)
       ! ADT changes start here
       !py(i)=0.5_sp*(sumcy**2/sumc+sumsy**2/sums)/var
       if (var == 0.0) then
          py(i) = 0.0
       else
          py(i)=0.5_sp*(sumcy**2/sumc+sumsy**2/sums)/var
       end if
       iy = sumsy/sqrt(sums)
       ry = sumcy/sqrt(sumc)
       phLS=atan2(iy,ry)
       !arg0=twopi*(xave+xmin)*pnow+wtau
       arg1=twopi*xave*pnow+wtau
       !ph(i) =mod(phLS+arg0,twopi)
       ph1(i)=mod(phLS+arg1,twopi)
       ! End ADT changes
       pnow=pnow+1.0_sp/(ofac*xdif)
    end do
    jmax=imaxloc(py(1:nout))
    expy=exp(-py(jmax))
    effm=2.0_sp*nout/ofac
    prob=effm*expy
    if (prob > 0.01_sp) prob=1.0_sp-(1.0_sp-expy)**effm
    !ADT changes start here
    nlen=2*nout+1
    do i=1,nout
       amp=sqrt(var*nlen/2.)*sqrt(py(i))
       prl(i)=amp*cos(ph1(i))
       pim(i)=amp*sin(ph1(i))
    end do
    deallocate(ph1)
    ! End ADT changes
  END SUBROUTINE period

!-----------------------------------------------------------------------------
!
! AVEVAR.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE avevar(data,ave,var)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: data
    REAL(SP), INTENT(OUT) :: ave,var
    INTEGER(I4B) :: n
    REAL(SP), DIMENSION(size(data)) :: s
    n=size(data)
    ave=sum(data(:))/n
    s(:)=data(:)-ave
    var=dot_product(s,s)
    var=(var-sum(s)**2/n)/(n-1)
  END SUBROUTINE avevar

!-----------------------------------------------------------------------------
!
! SPLINE/SPLINT.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE spline_r(x,y,yp1,ypn,y2)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
    REAL(SP), INTENT(IN) :: yp1,ypn
    REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER(I4B) :: n
    REAL(SP), DIMENSION(size(x)) :: a,b,c,r
    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99e30_sp) then
       r(1)=0.0
       c(1)=0.0
    else
       r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5
    end if
    if (ypn > 0.99e30_sp) then
       r(n)=0.0
       a(n)=0.0
    else
       r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  END SUBROUTINE spline_r

  SUBROUTINE spline_d(x,y,yp1,ypn,y2)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
    REAL(DP), INTENT(IN) :: yp1,ypn
    REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER(I4B) :: n
    REAL(DP), DIMENSION(size(x)) :: a,b,c,r
    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99e30_sp) then
       r(1)=0.0
       c(1)=0.0
    else
       r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5
    end if
    if (ypn > 0.99e30_sp) then
       r(n)=0.0
       a(n)=0.0
    else
       r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  END SUBROUTINE spline_d

  FUNCTION splint_r(xa,ya,y2a,x)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: splint_r
    INTEGER(I4B) :: khi,klo,n
    REAL(SP) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_r=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
  END FUNCTION splint_r

  FUNCTION splint_d(xa,ya,y2a,x)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: splint_d
    INTEGER(I4B) :: khi,klo,n
    REAL(DP) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_d=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
  END FUNCTION splint_d

!-----------------------------------------------------------------------------
!
! POLINT.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE polint_r(xa,ya,x,y,dy)
    !USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
    REAL(SP), INTENT(IN) :: x
    REAL(SP), INTENT(OUT) :: y,dy
    INTEGER(I4B) :: m,n,ns
    REAL(SP), DIMENSION(SIZE(xa)) :: c,d,den,ho
    n=ASSERT_EQ(SIZE(xa),SIZE(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=IMINLOC(ABS(x-xa))
    y=ya(ns)
    ns=ns-1
    DO m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       IF (ANY(den(1:n-m) == 0.0)) &
            CALL NRERROR('polint: calculation failure')
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       IF (2*ns < n-m) THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       END IF
       y=y+dy
    END DO
  END SUBROUTINE polint_r

  SUBROUTINE polint_d(xa,ya,x,y,dy)
    !USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(OUT) :: y,dy
    INTEGER(I4B) :: m,n,ns
    REAL(DP), DIMENSION(SIZE(xa)) :: c,d,den,ho
    n=ASSERT_EQ(SIZE(xa),SIZE(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=IMINLOC(ABS(x-xa))
    y=ya(ns)
    ns=ns-1
    DO m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       IF (ANY(den(1:n-m) == 0.0)) &
            CALL NRERROR('polint: calculation failure')
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       IF (2*ns < n-m) THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       END IF
       y=y+dy
    END DO
  END SUBROUTINE polint_d

!-----------------------------------------------------------------------------
!
! ZROOTS.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE zroots_c(a,roots,polish)
    !USE nrtype; USE nrutil, ONLY : assert_eq,poly_term
    !USE nr, ONLY : laguer,indexx
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: roots
    LOGICAL(LGT), INTENT(IN) :: polish
    REAL(SP), PARAMETER :: EPS=1.0e-6_sp
    INTEGER(I4B) :: j,its,m
    INTEGER(I4B), DIMENSION(size(roots)) :: indx
    COMPLEX(SPC) :: x
    COMPLEX(SPC), DIMENSION(size(a)) :: ad
    m=assert_eq(size(roots),size(a)-1,'zroots')
    ad(:)=a(:)
    do j=m,1,-1
       x=cmplx(0.0_sp,kind=spc)
       call laguer(ad(1:j+1),x,its)
       if (abs(aimag(x)) <= 2.0_sp*EPS**2*abs(real(x))) &
            x=cmplx(real(x),kind=spc)
       roots(j)=x
       ad(j:1:-1)=poly_term(ad(j+1:2:-1),x)
    end do
    if (polish) then
       do j=1,m
          call laguer(a(:),roots(j),its)
       end do
    end if
    call indexx(real(roots),indx)
    roots=roots(indx)
  END SUBROUTINE zroots_c

  SUBROUTINE zroots_z(a,roots,polish)
    !USE nrtype; USE nrutil, ONLY : assert_eq,poly_term
    !USE nr, ONLY : laguer,indexx
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(DPC), DIMENSION(:), INTENT(OUT) :: roots
    LOGICAL(LGT), INTENT(IN) :: polish
    REAL(DP), PARAMETER :: EPS=1.0e-6_dp
    INTEGER(I4B) :: j,its,m
    INTEGER(I4B), DIMENSION(size(roots)) :: indx
    COMPLEX(DPC) :: x
    COMPLEX(DPC), DIMENSION(size(a)) :: ad
    m=assert_eq(size(roots),size(a)-1,'zroots')
    ad(:)=a(:)
    do j=m,1,-1
       x=cmplx(0.0_dp,kind=dpc)
       call laguer(ad(1:j+1),x,its)
       if (abs(aimag(x)) <= 2.0_dp*EPS**2*abs(real(x))) &
            x=cmplx(real(x),kind=dpc)
       roots(j)=x
       ad(j:1:-1)=poly_term(ad(j+1:2:-1),x)
    end do
    if (polish) then
       do j=1,m
          call laguer(a(:),roots(j),its)
       end do
    end if
    call indexx(real(roots),indx)
    roots=roots(indx)
  END SUBROUTINE zroots_z

!-----------------------------------------------------------------------------
!
! LAGUER.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE laguer_c(a,x,its)
    !USE nrtype; USE nrutil, ONLY : nrerror,poly,poly_term
    IMPLICIT NONE
    INTEGER(I4B), INTENT(OUT) :: its
    COMPLEX(SPC), INTENT(INOUT) :: x
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
    REAL(SP), PARAMETER :: EPS=epsilon(1.0_sp)
    INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
    INTEGER(I4B) :: iter,m
    REAL(SP) :: abx,abp,abm,err
    COMPLEX(SPC) :: dx,x1,f,g,h,sq,gp,gm,g2
    COMPLEX(SPC), DIMENSION(size(a)) :: b,d
    REAL(SP), DIMENSION(MR) :: frac = &
         (/ 0.5_sp,0.25_sp,0.75_sp,0.13_sp,0.38_sp,0.62_sp,0.88_sp,1.0_sp /)
    m=size(a)-1
    do iter=1,MAXIT
       its=iter
       abx=abs(x)
       b(m+1:1:-1)=poly_term(a(m+1:1:-1),x)
       d(m:1:-1)=poly_term(b(m+1:2:-1),x)
       f=poly(x,d(2:m))
       err=EPS*poly(abx,abs(b(1:m+1)))
       if (abs(b(1)) <= err) RETURN
       g=d(1)/b(1)
       g2=g*g
       h=g2-2.0_sp*f/b(1)
       sq=sqrt((m-1)*(m*h-g2))
       gp=g+sq
       gm=g-sq
       abp=abs(gp)
       abm=abs(gm)
       if (abp < abm) gp=gm
       if (max(abp,abm) > 0.0) then
          dx=m/gp
       else
          dx=exp(cmplx(log(1.0_sp+abx),iter,kind=spc))
       end if
       x1=x-dx
       if (x == x1) RETURN
       if (mod(iter,MT) /= 0) then
          x=x1
       else
          x=x-dx*frac(iter/MT)
       end if
    end do
    call nrerror('laguer: too many iterations')
  END SUBROUTINE laguer_c

  SUBROUTINE laguer_z(a,x,its)
    !USE nrtype; USE nrutil, ONLY : nrerror,poly,poly_term
    IMPLICIT NONE
    INTEGER(I4B), INTENT(OUT) :: its
    COMPLEX(DPC), INTENT(INOUT) :: x
    COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
    REAL(DP), PARAMETER :: EPS=epsilon(1.0_dp)
    INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
    INTEGER(I4B) :: iter,m
    REAL(DP) :: abx,abp,abm,err
    COMPLEX(DPC) :: dx,x1,f,g,h,sq,gp,gm,g2
    COMPLEX(DPC), DIMENSION(size(a)) :: b,d
    REAL(DP), DIMENSION(MR) :: frac = &
         (/ 0.5_dp,0.25_dp,0.75_dp,0.13_dp,0.38_dp,0.62_dp,0.88_dp,1.0_dp /)
    m=size(a)-1
    do iter=1,MAXIT
       its=iter
       abx=abs(x)
       b(m+1:1:-1)=poly_term(a(m+1:1:-1),x)
       d(m:1:-1)=poly_term(b(m+1:2:-1),x)
       f=poly(x,d(2:m))
       err=EPS*poly(abx,abs(b(1:m+1)))
       if (abs(b(1)) <= err) RETURN
       g=d(1)/b(1)
       g2=g*g
       h=g2-2.0_dp*f/b(1)
       sq=sqrt((m-1)*(m*h-g2))
       gp=g+sq
       gm=g-sq
       abp=abs(gp)
       abm=abs(gm)
       if (abp < abm) gp=gm
       if (max(abp,abm) > 0.0) then
          dx=m/gp
       else
          dx=exp(cmplx(log(1.0_dp+abx),iter,kind=dpc))
       end if
       x1=x-dx
       if (x == x1) RETURN
       if (mod(iter,MT) /= 0) then
          x=x1
       else
          x=x-dx*frac(iter/MT)
       end if
    end do
    call nrerror('laguer: too many iterations')
  END SUBROUTINE laguer_z

!-----------------------------------------------------------------------------
!
! ZBRAK.F
!
!-----------------------------------------------------------------------------

  SUBROUTINE zbrak_r(func,x1,x2,n,xb1,xb2,nb)
    !USE nrtype; USE nrutil, ONLY : arth
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B), INTENT(OUT) :: nb
    REAL(SP), INTENT(IN) :: x1,x2
    REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
    INTERFACE
       FUNCTION func(x)
         !USE nrtype
         IMPLICIT NONE
         !REAL(SP), INTENT(IN) :: x
         !REAL(SP) :: func
         REAL, INTENT(IN) :: x
         REAL :: func
       END FUNCTION func
    END INTERFACE
    INTEGER(I4B) :: i
    REAL(SP) :: dx
    REAL(SP), DIMENSION(0:n) :: f,x
    LOGICAL(LGT), DIMENSION(1:n) :: mask
    LOGICAL(LGT), SAVE :: init=.true.
    if (init) then
       init=.false.
       nullify(xb1,xb2)
    end if
    if (associated(xb1)) deallocate(xb1)
    if (associated(xb2)) deallocate(xb2)
    dx=(x2-x1)/n
    x=x1+dx*arth(0,1,n+1)
    do i=0,n
       f(i)=func(x(i))
    end do
    mask=f(1:n)*f(0:n-1) <= 0.0
    nb=count(mask)
    allocate(xb1(nb),xb2(nb))
    xb1(1:nb)=pack(x(0:n-1),mask)
    xb2(1:nb)=pack(x(1:n),mask)
  END SUBROUTINE zbrak_r

  SUBROUTINE zbrak_d(func,x1,x2,n,xb1,xb2,nb)
    !USE nrtype; USE nrutil, ONLY : arth
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B), INTENT(OUT) :: nb
    REAL(DP), INTENT(IN) :: x1,x2
    REAL(DP), DIMENSION(:), POINTER :: xb1,xb2
    INTERFACE
       FUNCTION func(x)
         !USE nrtype
         IMPLICIT NONE
         !REAL(DP), INTENT(IN) :: x
         !REAL(DP) :: func
         REAL(KIND(0d0)), INTENT(IN) :: x
         REAL(KIND(0d0)) :: func
       END FUNCTION func
    END INTERFACE
    INTEGER(I4B) :: i
    REAL(DP) :: dx
    REAL(DP), DIMENSION(0:n) :: f,x
    LOGICAL(LGT), DIMENSION(1:n) :: mask
    LOGICAL(LGT), SAVE :: init=.true.
    if (init) then
       init=.false.
       nullify(xb1,xb2)
    end if
    if (associated(xb1)) deallocate(xb1)
    if (associated(xb2)) deallocate(xb2)
    dx=(x2-x1)/n
    x=x1+dx*arth(0,1,n+1)
    do i=0,n
       f(i)=func(x(i))
    end do
    mask=f(1:n)*f(0:n-1) <= 0.d0
    nb=count(mask)
    allocate(xb1(nb),xb2(nb))
    xb1(1:nb)=pack(x(0:n-1),mask)
    xb2(1:nb)=pack(x(1:n),mask)
  END SUBROUTINE zbrak_d

  SUBROUTINE zbrak2_r(func,x1,x2,n,y,xb1,xb2,nb)
    !USE nrtype; USE nrutil, ONLY : arth
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B), INTENT(OUT) :: nb
    REAL(SP), INTENT(IN) :: x1,x2,y
    REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
    INTERFACE
       FUNCTION func(x,y)
         !USE nrtype
         IMPLICIT NONE
         !REAL(SP), INTENT(IN) :: x,y
         !REAL(SP) :: func
         REAL, INTENT(IN) :: x,y
         REAL :: func
       END FUNCTION func
    END INTERFACE
    INTEGER(I4B) :: i
    REAL(SP) :: dx
    REAL(SP), DIMENSION(0:n) :: f,x
    LOGICAL(LGT), DIMENSION(1:n) :: mask
    LOGICAL(LGT), SAVE :: init=.true.
    if (init) then
       init=.false.
       nullify(xb1,xb2)
    end if
    if (associated(xb1)) deallocate(xb1)
    if (associated(xb2)) deallocate(xb2)
    dx=(x2-x1)/n
    x=x1+dx*arth(0,1,n+1)
    do i=0,n
       f(i)=func(x(i),y)
    end do
    mask=f(1:n)*f(0:n-1) <= 0.0
    nb=count(mask)
    allocate(xb1(nb),xb2(nb))
    xb1(1:nb)=pack(x(0:n-1),mask)
    xb2(1:nb)=pack(x(1:n),mask)
  END SUBROUTINE zbrak2_r

  SUBROUTINE zbrak2_d(func,x1,x2,n,y,xb1,xb2,nb)
    !USE nrtype; USE nrutil, ONLY : arth
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B), INTENT(OUT) :: nb
    REAL(DP), INTENT(IN) :: x1,x2,y
    REAL(DP), DIMENSION(:), POINTER :: xb1,xb2
    INTERFACE
       FUNCTION func(x,y)
         !USE nrtype
         IMPLICIT NONE
         !REAL(DP), INTENT(IN) :: x,y
         !REAL(DP) :: func
         REAL(KIND(0d0)), INTENT(IN) :: x,y
         REAL(KIND(0d0)) :: func
       END FUNCTION func
    END INTERFACE
    INTEGER(I4B) :: i
    REAL(DP) :: dx
    REAL(DP), DIMENSION(0:n) :: f,x
    LOGICAL(LGT), DIMENSION(1:n) :: mask
    LOGICAL(LGT), SAVE :: init=.true.
    if (init) then
       init=.false.
       nullify(xb1,xb2)
    end if
    if (associated(xb1)) deallocate(xb1)
    if (associated(xb2)) deallocate(xb2)
    dx=(x2-x1)/n
    x=x1+dx*arth(0,1,n+1)
    do i=0,n
       f(i)=func(x(i),y)
    end do
    mask=f(1:n)*f(0:n-1) <= 0.d0
    nb=count(mask)
    allocate(xb1(nb),xb2(nb))
    xb1(1:nb)=pack(x(0:n-1),mask)
    xb2(1:nb)=pack(x(1:n),mask)
  END SUBROUTINE zbrak2_d

!-----------------------------------------------------------------------------
!
! RTSAFE.F
!
!-----------------------------------------------------------------------------

  FUNCTION rtsafe_r(funcd,x1,x2,xacc)
    !USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: x1,x2,xacc
    REAL(SP) :: rtsafe_r
    INTERFACE
       SUBROUTINE funcd(x,fval,fderiv)
         !USE nrtype
         IMPLICIT NONE
         !REAL(SP), INTENT(IN) :: x
         !REAL(SP), INTENT(OUT) :: fval,fderiv
         REAL, INTENT(IN) :: x
         REAL, INTENT(OUT) :: fval,fderiv
       END SUBROUTINE funcd
    END INTERFACE
    INTEGER(I4B), PARAMETER :: MAXIT=100
    INTEGER(I4B) :: j
    REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
    call funcd(x1,fl,df)
    call funcd(x2,fh,df)
    if ( (fl > 0.0 .and. fh > 0.0) .or. &
         (fl < 0.0 .and. fh < 0.0) ) &
         call nrerror('root must be bracketed in rtsafe')
    if (fl == 0.0) then
       rtsafe_r=x1
       RETURN
    else if (fh == 0.0) then
       rtsafe_r=x2
       RETURN
    else if (fl < 0.0) then
       xl=x1
       xh=x2
    else
       xh=x1
       xl=x2
    end if
    rtsafe_r=0.5_sp*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(rtsafe_r,f,df)
    do j=1,MAXIT
       if (((rtsafe_r-xh)*df-f)*((rtsafe_r-xl)*df-f) > 0.0 .or. &
            abs(2.0_sp*f) > abs(dxold*df) ) then
          dxold=dx
          dx=0.5_sp*(xh-xl)
          rtsafe_r=xl+dx
          if (xl == rtsafe_r) RETURN
       else
          dxold=dx
          dx=f/df
          temp=rtsafe_r
          rtsafe_r=rtsafe_r-dx
          if (temp == rtsafe_r) RETURN
       end if
       if (abs(dx) < xacc) RETURN
       call funcd(rtsafe_r,f,df)
       if (f < 0.0) then
          xl=rtsafe_r
       else
          xh=rtsafe_r
       end if
    end do
    call nrerror('rtsafe: exceeded maximum iterations')
  END FUNCTION rtsafe_r

  FUNCTION rtsafe_d(funcd,x1,x2,xacc)
    !USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,xacc
    REAL(DP) :: rtsafe_d
    INTERFACE
       SUBROUTINE funcd(x,fval,fderiv)
         !USE nrtype
         IMPLICIT NONE
         !REAL(DP), INTENT(IN) :: x
         !REAL(DP), INTENT(OUT) :: fval,fderiv
         REAL(KIND(0d0)), INTENT(IN) :: x
         REAL(KIND(0d0)), INTENT(OUT) :: fval,fderiv
       END SUBROUTINE funcd
    END INTERFACE
    INTEGER(I4B), PARAMETER :: MAXIT=100
    INTEGER(I4B) :: j
    REAL(DP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
    call funcd(x1,fl,df)
    call funcd(x2,fh,df)
    if ( (fl > 0.0 .and. fh > 0.0) .or. &
         (fl < 0.0 .and. fh < 0.0) ) &
         call nrerror('root must be bracketed in rtsafe')
    if (fl == 0.0) then
       rtsafe_d=x1
       RETURN
    else if (fh == 0.0) then
       rtsafe_d=x2
       RETURN
    else if (fl < 0.0) then
       xl=x1
       xh=x2
    else
       xh=x1
       xl=x2
    end if
    rtsafe_d=0.5_dp*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(rtsafe_d,f,df)
    do j=1,MAXIT
       if (((rtsafe_d-xh)*df-f)*((rtsafe_d-xl)*df-f) > 0.0 .or. &
            abs(2.0_dp*f) > abs(dxold*df) ) then
          dxold=dx
          dx=0.5_dp*(xh-xl)
          rtsafe_d=xl+dx
          if (xl == rtsafe_d) RETURN
       else
          dxold=dx
          dx=f/df
          temp=rtsafe_d
          rtsafe_d=rtsafe_d-dx
          if (temp == rtsafe_d) RETURN
       end if
       if (abs(dx) < xacc) RETURN
       call funcd(rtsafe_d,f,df)
       if (f < 0.0) then
          xl=rtsafe_d
       else
          xh=rtsafe_d
       end if
    end do
    call nrerror('rtsafe: exceeded maximum iterations')
  END FUNCTION rtsafe_d

  FUNCTION rtsafe2_r(funcd,x1,x2,xacc,y)
    !USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: x1,x2,xacc,y
    REAL(SP) :: rtsafe2_r
    INTERFACE
       SUBROUTINE funcd(x,y,fval,fderiv)
         !USE nrtype
         IMPLICIT NONE
         !REAL(SP), INTENT(IN) :: x, y
         !REAL(SP), INTENT(OUT) :: fval,fderiv
         REAL, INTENT(IN) :: x, y
         REAL, INTENT(OUT) :: fval,fderiv
       END SUBROUTINE funcd
    END INTERFACE
    INTEGER(I4B), PARAMETER :: MAXIT=100
    INTEGER(I4B) :: j
    REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
    call funcd(x1,y,fl,df)
    call funcd(x2,y,fh,df)
    if ((fl > 0.0 .and. fh > 0.0) .or. &
         (fl < 0.0 .and. fh < 0.0)) &
         call nrerror('root must be bracketed in rtsafe')
    if (fl == 0.0) then
       rtsafe2_r=x1
       RETURN
    else if (fh == 0.0) then
       rtsafe2_r=x2
       RETURN
    else if (fl < 0.0) then
       xl=x1
       xh=x2
    else
       xh=x1
       xl=x2
    end if
    rtsafe2_r=0.5_sp*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(rtsafe2_r,y,f,df)
    do j=1,MAXIT
       if (((rtsafe2_r-xh)*df-f)*((rtsafe2_r-xl)*df-f) > 0.0 .or. &
            abs(2.0_sp*f) > abs(dxold*df) ) then
          dxold=dx
          dx=0.5_sp*(xh-xl)
          rtsafe2_r=xl+dx
          if (xl == rtsafe2_r) RETURN
       else
          dxold=dx
          dx=f/df
          temp=rtsafe2_r
          rtsafe2_r=rtsafe2_r-dx
          if (temp == rtsafe2_r) RETURN
       end if
       if (abs(dx) < xacc) RETURN
       call funcd(rtsafe2_r,y,f,df)
       if (f < 0.0) then
          xl=rtsafe2_r
       else
          xh=rtsafe2_r
       end if
    end do
    call nrerror('rtsafe: exceeded maximum iterations')
  END FUNCTION rtsafe2_r

  FUNCTION rtsafe2_d(funcd,x1,x2,xacc,y)
    !USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,xacc,y
    REAL(DP) :: rtsafe2_d
    INTERFACE
       SUBROUTINE funcd(x,y,fval,fderiv)
         !USE nrtype
         IMPLICIT NONE
         !REAL(DP), INTENT(IN) :: x, y
         !REAL(DP), INTENT(OUT) :: fval,fderiv
         REAL(KIND(0d0)), INTENT(IN) :: x, y
         REAL(KIND(0d0)), INTENT(OUT) :: fval,fderiv
       END SUBROUTINE funcd
    END INTERFACE
    INTEGER(I4B), PARAMETER :: MAXIT=100
    INTEGER(I4B) :: j
    REAL(DP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
    call funcd(x1,y,fl,df)
    call funcd(x2,y,fh,df)
    if ((fl > 0.0 .and. fh > 0.0) .or. &
         (fl < 0.0 .and. fh < 0.0)) &
         call nrerror('root must be bracketed in rtsafe')
    if (fl == 0.0) then
       rtsafe2_d=x1
       RETURN
    else if (fh == 0.0) then
       rtsafe2_d=x2
       RETURN
    else if (fl < 0.0) then
       xl=x1
       xh=x2
    else
       xh=x1
       xl=x2
    end if
    rtsafe2_d=0.5_dp*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(rtsafe2_d,y,f,df)
    do j=1,MAXIT
       if (((rtsafe2_d-xh)*df-f)*((rtsafe2_d-xl)*df-f) > 0.0 .or. &
            abs(2.0_dp*f) > abs(dxold*df) ) then
          dxold=dx
          dx=0.5_dp*(xh-xl)
          rtsafe2_d=xl+dx
          if (xl == rtsafe2_d) RETURN
       else
          dxold=dx
          dx=f/df
          temp=rtsafe2_d
          rtsafe2_d=rtsafe2_d-dx
          if (temp == rtsafe2_d) RETURN
       end if
       if (abs(dx) < xacc) RETURN
       call funcd(rtsafe2_d,y,f,df)
       if (f < 0.0) then
          xl=rtsafe2_d
       else
          xh=rtsafe2_d
       end if
    end do
    call nrerror('rtsafe: exceeded maximum iterations')
  END FUNCTION rtsafe2_d

  FUNCTION rtsafe3_r(funcd,x1,x2,xacc,y,z)
    !USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: x1,x2,xacc,y,z
    REAL(SP) :: rtsafe3_r
    INTERFACE
       SUBROUTINE funcd(x,y,z,fval,fderiv)
         !USE nrtype
         IMPLICIT NONE
         !REAL(SP), INTENT(IN) :: x,y,z
         !REAL(SP), INTENT(OUT) :: fval,fderiv
         REAL, INTENT(IN) :: x,y,z
         REAL, INTENT(OUT) :: fval,fderiv
       END SUBROUTINE funcd
    END INTERFACE
    INTEGER(I4B), PARAMETER :: MAXIT=100
    INTEGER(I4B) :: j
    REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
    call funcd(x1,y,z,fl,df)
    call funcd(x2,y,z,fh,df)
    if ((fl > 0.0 .and. fh > 0.0) .or. &
         (fl < 0.0 .and. fh < 0.0)) &
         call nrerror('root must be bracketed in rtsafe')
    if (fl == 0.0) then
       rtsafe3_r=x1
       RETURN
    else if (fh == 0.0) then
       rtsafe3_r=x2
       RETURN
    else if (fl < 0.0) then
       xl=x1
       xh=x2
    else
       xh=x1
       xl=x2
    end if
    rtsafe3_r=0.5_sp*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(rtsafe3_r,y,z,f,df)
    do j=1,MAXIT
       if (((rtsafe3_r-xh)*df-f)*((rtsafe3_r-xl)*df-f) > 0.0 .or. &
            abs(2.0_sp*f) > abs(dxold*df) ) then
          dxold=dx
          dx=0.5_sp*(xh-xl)
          rtsafe3_r=xl+dx
          if (xl == rtsafe3_r) RETURN
       else
          dxold=dx
          dx=f/df
          temp=rtsafe3_r
          rtsafe3_r=rtsafe3_r-dx
          if (temp == rtsafe3_r) RETURN
       end if
       if (abs(dx) < xacc) RETURN
       call funcd(rtsafe3_r,y,z,f,df)
       if (f < 0.0) then
          xl=rtsafe3_r
       else
          xh=rtsafe3_r
       end if
    end do
    call nrerror('rtsafe: exceeded maximum iterations')
  END FUNCTION rtsafe3_r

  FUNCTION rtsafe3_d(funcd,x1,x2,xacc,y,z)
    !USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,xacc,y,z
    REAL(DP) :: rtsafe3_d
    INTERFACE
       SUBROUTINE funcd(x,y,z,fval,fderiv)
         !USE nrtype
         IMPLICIT NONE
         !REAL(DP), INTENT(IN) :: x,y,z
         !REAL(DP), INTENT(OUT) :: fval,fderiv
         REAL(KIND(0d0)), INTENT(IN) :: x,y,z
         REAL(KIND(0d0)), INTENT(OUT) :: fval,fderiv
       END SUBROUTINE funcd
    END INTERFACE
    INTEGER(I4B), PARAMETER :: MAXIT=100
    INTEGER(I4B) :: j
    REAL(DP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
    call funcd(x1,y,z,fl,df)
    call funcd(x2,y,z,fh,df)
    if ((fl > 0.0 .and. fh > 0.0) .or. &
         (fl < 0.0 .and. fh < 0.0)) &
         call nrerror('root must be bracketed in rtsafe')
    if (fl == 0.0) then
       rtsafe3_d=x1
       RETURN
    else if (fh == 0.0) then
       rtsafe3_d=x2
       RETURN
    else if (fl < 0.0) then
       xl=x1
       xh=x2
    else
       xh=x1
       xl=x2
    end if
    rtsafe3_d=0.5_dp*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(rtsafe3_d,y,z,f,df)
    do j=1,MAXIT
       if (((rtsafe3_d-xh)*df-f)*((rtsafe3_d-xl)*df-f) > 0.0 .or. &
            abs(2.0_dp*f) > abs(dxold*df) ) then
          dxold=dx
          dx=0.5_dp*(xh-xl)
          rtsafe3_d=xl+dx
          if (xl == rtsafe3_d) RETURN
       else
          dxold=dx
          dx=f/df
          temp=rtsafe3_d
          rtsafe3_d=rtsafe3_d-dx
          if (temp == rtsafe3_d) RETURN
       end if
       if (abs(dx) < xacc) RETURN
       call funcd(rtsafe3_d,y,z,f,df)
       if (f < 0.0) then
          xl=rtsafe3_d
       else
          xh=rtsafe3_d
       end if
    end do
    call nrerror('rtsafe: exceeded maximum iterations')
  END FUNCTION rtsafe3_d

END MODULE module_nrutils
