!c-mm------------------------------------------------------------------
!c-mm Using this routine to calculate the Chapman function, as an alternative
!c-mm    to dividing by cos(SZA).  Should be better higher up.
!c-mm------------------------------------------------------------------
      MODULE module_ra_chapman

      CONTAINS
!C*********************************************************************C
!C*                                                                   *C
!C*  chapman.for                                                      *C
!C*                                                                   *C
!C*  Written by:  David L. Huestis, Molecular Physics Laboratory      *C
!C*                                                                   *C
!C*  Copyright (c) 2000  SRI International                            *C
!C*  All Rights Reserved                                              *C
!C*                                                                   *C
!C*  This software is provided on an as is basis; without any         *C
!C*  warranty; without the implied warranty of merchantability or     *C
!C*  fitness for a particular purpose.                                *C
!C*                                                                   *C
!C*********************************************************************C
!C*
!C*	To calculate the Chapman Function, Ch(X,chi0), the column 
!C*	depth of an exponential atmosphere integrated along a line 
!C*	from a given point to the sun, divided by the column depth for 
!C*	a vertical sun.
!C*
!C*  USAGE:
!C*
!C*	  z = altitude above the surface
!C*	  R = radius of the planet
!C*	  H = atmospheric scale height
!C*
!C*	  X = (R+z)/H
!C*	  chi0 = solar zenith angle (in degrees)
!C*
!C*	  IMPLICIT real*4(a-h,o-z)
!C*	  depth = atm_chapman(X,chi0)	! analytical
!C*	  depth = atm_chap_num(X,chi0)	! numerical (chi0 <= 90)
!C*
!C*	  IMPLICIT real*8(a-h,o-z)
!C*	  depth = atm8_chapman(X,chi0)	! analytical
!C*	  depth = atm8_chap_num(X,chi0)	! numerical (chi0 <= 90)
!C*
!C*  PERFORMANCE:
!C*
!C*	Compiled and linked using Microsoft FORTRAN 5.1, and executed 
!C*	in MS-DOS mode under Windows 95 on a 160 MHz PC.
!C*
!C*    TIMING (in microseconds, typical)
!C*
!C*	  120	atm_chapman and atm8_chapman for X < 36
!C*	   25	atm_chapman and atm8_chapman for X >= 36
!C*	  500	atm_chap_num
!C*	 5000	atm8_chap_num
!C*
!C*    ACCURACY (maximum relative error, 0<=chi0<=90, 1<=X<=820)
!C*
!C*	6.0E-7	atm_chapman and atm8_chapman for X < 60
!C*	1.5E-7	atm_chapman and atm8_chapman for X >= 60
!C*	6.0E-8	atm_chap_num
!C*	1.E-15	atm8_chap_num (convergence test)
!C*
!C*    CODING
!C*
!C*	No claims are made that the code is optimized for speed, 
!C*	accuracy, or compactness.  The principal objectives were 
!C*
!C*	  (1) Robustness with respect to argument values
!C*	  (2) Rigorous mathematical derivation and error control
!C*	  (3) Maximal use of "well known" mathematical functions
!C*	  (4) Ease of readability and mapping of theory to coding
!C*
!C*	The real*8 accuracy could be improved with more accurate 
!C*	representations of E1(), erfc(), I0(), I1(), K0(), K1().
!C*
!C*	In the course of development, many representations and 
!C*	approximations of the Chapman Function were attempted that 
!C*	failed to be robustly extendable to machine-precision.
!C*
!C*  INTERNET ACCESS:
!C*
!C*	Source: http://www-mpl.sri.com/software/chapman/chapman.html
!C*	Author: mailto:david.huestis@sri.com
!C*	        http://www-mpl.sri.com/bios/Huestis-DL.html
!C*
!C*  EDIT HISTORY:
!C*
!C*	01/22/2000 DLH	First complete documentation
!C*
!C*	01/15/2000 DLH	First complete version of chapman.for
!C*
!C**********************************************************************
!C*
!C*  THEORY:
!C*
!C*    INTRODUCTION
!C*
!C*	    This computer code models the absorption of solar radiation 
!C*	by an atmosphere that depends exponentionally on altitude.  In 
!C*	specific we calculate the effective column depth of a species 
!C*	of local density, n(z), from a point at a given altitude, z0, 
!C*	to the sun at a given solar zenith angle, chi0.  Following Rees 
!C*	[Re89, Section 2.2] we write the column depth for chi0 <= 90 
!C*	degrees as
!C*
!C*   (A)  N(z0,chi0) = int{z=z0,infinity} 
!C*	     [ n(z)/sqrt( 1 - ( sin(chi0) * (R+z0) / (R+z) ) **2 ) dz ]
!C*
!C*	where R is the radius of the solid planet (e.g. Earth).  For 
!C*	chi0 .gt. 90 degrees we write
!C*
!C*	  N(z0,chi0) = 2*N(zs,90) - N(z0,180-chi0)
!C*
!C*	where zs = (R+z0)*sin(chi0)-R is the tangent height.
!C*
!C*	    For an exponential atmosphere, with
!C*
!C*	  n(z) = n(z0) * exp(-(z-z0)/H)
!C*
!C*	with a constant scale height, H, the column depth can be 
!C*	represented by the Chapman function, Ch(X,chi0), named after 
!C*	the author of the first quantitative mathematical investigation 
!C*	[Ch31b] trough the relation
!C*
!C*	  N(z0,chi0) = H * n(z0) * Ch(X,chi0)
!C*
!C*	where X = (R+z0)/H is a DIMENSIONless measure of the radius 
!C*	of curvature, with values from about 300 to 1300 on Earth.
!C*
!C*
!C*    APPROACH
!C*
!C*	    We provide function entry points for very stable and 
!C*	reasonably efficient evaluation of Ch(X,chi0) with full 
!C*	single-precision accuracy (<= 6.0E-7 relative) for a wide 
!C*	range of PARAMETERs.  A 15-digit-accurate double precision 
!C*	numerical integration routine is also provided.
!C*
!C*	    Below we will develop (1) a compact asymptotic expansion of 
!C*	good accuracy for moderately large values of X (.gt. 36) and all 
!C*	values of chi0, (2) an efficient numerical integral for 
!C*	all values of X and chi0, and (3) an explicit analytical 
!C*	representation, valid for all values of X and chi0, based 
!C*	the differential equation satisfied by Ch(X,chi0).
!C*
!C*	    All three of these represent new research results as well 
!C*	as significant computational improvements over the previous 
!C*	literature, much of which is cited below.
!C*
!C*
!C*    CHANGES OF THE VARIABLE OF INTEGRATION
!C*
!C*	Substituting y = (R+z)/(R+z0) - 1 we find
!C*
!C*   (B)  Ch(X,chi0) = X * int{y=0,infinity}
!C*	     [ exp(-X*y) / sqrt( 1 - ( sin(chi0) / (1+y) )**2 ) dy ]
!C*
!C*	The futher substitutions s = (1+y)/sin(chi0), s0 = 1/sin(chi0) 
!C*	give
!C*
!C*   (C)  Ch(X,chi0) = X*sin(chi0) * int{s=s0,infinity}
!C*	     [ exp(X*(1-sin(chi0)*s)) * s / sqrt(s**2-1) ds ]
!C*
!C*	From this equation we can establish that
!C*
!C*	  Ch(X,90) = X*exp(X)*K1(X)
!C*
!C*	[AS64, Equations 9.6.23 and 9.6.27].  If we now substitute
!C*	s = 1/sin(lambda) we obtain
!C*
!C*   (D)  Ch(X,chi0) = X*sin(chi0) * int{lambda=0,chi0} 
!C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) * csc(lambda)**2 dlambda]
!C*
!C*	which is the same as Chapman's original formulation [Ch31b, p486,
!C*	eqn (10)].  If we first expand the square root in (B)
!C*
!C*	  1/sqrt(1-q) = 1 + q/( sqrt(1-q)*(1+sqrt(1-q)) )
!C*
!C*	with q = ( sin(chi0) / (1+y) )**2 = sin(lambda)**2, we obtain 
!C*	a new form of (D) without numerical sigularities and simple 
!C*	convergence to Ch(0,chi0) = Ch(X,0) = 1
!C*
!C*   (E)  Ch(X,chi0) = 1 + X*sin(chi0) * int{lambda=0,chi0} 
!C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) 
!C*		/ (1 + cos(lambda) ) dlambda ]
!C*
!C*	Alternatively, we may substitute t**2 = y + t0**2, 
!C*	into Equation (B), with t0**2 = 1-sin(chi0), finding
!C*
!C*   (F)  Ch(X,chi0) = X * int{s=t0,infinity} 
!C*	    [ exp(-X*(t**2-t0**2)) * f(t,chi0) dt ]
!C* 
!C*	where
!C*
!C*	  f(t,chi0) = (t**2 + sin(chi0)) / sqrt(t**2+2*sin(chi0))
!C*
!C*	  f(t,chi0) = (t**2-t0**2+1)/sqrt(t**2-t0**2+1+sin(chi0))
!C*
!C*	    Below we will use Equation (F) above to develop a
!C*	compact asymptotic expansion of good accuracy for moderately 
!C*	large values of X (.gt. 36) and all values of chi0, Equation (E) 
!C*	to develop an efficient numerical integral for Ch(X,chi0) for 
!C*	all values of X and chi0, and Equation (C) to derive an explicit 
!C*	analytical representation, valid for all values of X and chi0,  
!C*	based on the differential equation satisfied by Ch(X,chi0).
!C*
!C*    atm_chapman(X,chi0) and atm8_chapman(X,chi0)
!C*
!C*	These routines RETURN real*4 and real*8 values of Ch(X,chi0)
!C*	selecting the asymptotic expansion or differential equation 
!C*	approaches, depending on the value of X.  These routines also 
!C*	handle the case of chi0 .gt. 90 degrees.
!C*
!C*    atm_chap_num(X,chi0) and atm8_chap_num(X,chi0)
!C*
!C*	These routines RETURN real*4 and real*8 values of Ch(X,chi0)
!C*	evaluated numerically.  They are both more accurate than the 
!C*	corresponding atm*_chapman() functions, but take significantly 
!C*	more CPU time.
!C*
!C*
!C*    ASYMPTOTIC EXPANSION
!C*
!C*	From Equation (F) we expand, with t0**2 = 1-sin(chi0), 
!C*
!C*	  f(t,chi0) = sum{n=0,3} [ C(n,chi0) * (t**2-t0**2)**n ]
!C*
!C*	The function atm8_chap_asy(X,chi0) evaluates integrals of the 
!C*	form
!C*
!C*	  int{t=t0,infinity} [exp(-X*(t**2-t0**2))*(t**2-t0**2)**n dt]
!C*
!C*	in terms of incomplete gamma functions, and sums them to 
!C*	compute Ch(X,chi0).  For large values of X, this results in an 
!C*	asymptotic expansion in negative powers of X, with coefficients 
!C*	that are stable for all values of chi0.
!C*
!C*	In contrast, the asymptotic expansions of Chapman [Ch31b, 
!C*	p488, Equation (22) and p490, Equation (38)], Hulburt [He39], 
!C*	and Swider [Sw64, p777, Equation (43)] use negative powers of 
!C*	X*cos(chi0)**2 or X*sin(chi0), and are accurate only for 
!C*	small values or large values of chi0, respectively.
!C*
!C*	Taking only the first term in the present expansion gives the 
!C*	simple formula
!C*
!C*	  Ch(X,chi0) = sqrt(pi*X/(1+sin(chi0))) * exp(X*(1-sin(chi0)))
!C*		* erfc( sqrt(X*(1-sin(chi0))) )
!C*
!C*	This is slightly more accurate than the semiempirical 
!C*	formula of Fitzmaurice [Fi64, Equation (3)], and sightly less 
!C*	accurate than that of Swider [Sw64, p780, Equation (52), 
!C*	corrected in SG69].
!C*
!C*
!C*    NUMERICAL INTEGRATION
!C*
!C*	We are integrating
!C*
!C*   (E)  Ch(X,chi0) = 1 + X*sin(chi0) * int{lambda=0,chi0} 
!C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) 
!C*		/ ( 1 + cos(lambda) ) dlambda ]
!C*
!C*	The integrand is numerically very smooth, and rapidly varying 
!C*	only near lambda = 0.  For X .ne. 0 we choose the lower limit 
!C*	of numerical integration such that the integrand is 
!C*	exponentially small, 7.0E-13 (3.0E-20 for real*8).  The domain 
!C*	of integration is divided into 64 equal intervals (6000 for 
!C*	real*8), and integrated numerically using the 9-point closed 
!C*	Newton-Cotes formula from Hildebrand [Hi56a, page 75, Equation
!C*	(3.5.17)].
!C*
!C*
!C*    INHOMOGENOUS DIFFERENTIAL EQUATION
!C*
!C*	    The function atm8_chap_deq(X,chi0) calculates Ch(X,chi0), 
!C*	based on Equation (C) above, using the inhomogeneous 
!C*	Bessel's equation as described below.  Consider the function 
!C*
!C*	  Z(Q) = int{s=s0,infinity} [ exp(-Q*s) / sqrt(s**2-1) ds ]
!C*
!C*	Differentiating with respect to Q we find that 
!C*
!C*	  Ch(X,chi0) = - Q * exp(X) * d/dQ [ Z(Q) ]
!C*
!C*	with Q = X*sin(chi0), s0 = 1/sin(chi0).  Differentiating 
!C*	inside the integral, we find that
!C*
!C*	  Z"(Q) + Z'(Q)/Q - Z(Q) = sqrt(s0**2-1) * exp(-Q*s0) / Q
!C*
!C*	giving us an inhomogeneous modified Bessel's equation of order 
!C*	zero.  Following Rabenstein [Ra66, pp43-45,149] the solution 
!C*	of this equation can be written as
!C*
!C*	  Z(Q) = A*I0(Q) + B*K0(Q) - sqrt(s0**2-1) 
!C*	         * int{t=Q,infinity} [ exp(-t*s0) 
!C*		   * ( I0(Q)*K0(t) - I0(t)*K0(Q) ) dt ] 
!C*
!C*	with coefficients A and B to be determined by matching 
!C*	boundary conditions.
!C*
!C*	    Differentiating with respect to Q we obtain
!C*
!C*	  Ch(X,chi0) = X*sin(chi0)*exp(X)*( 
!C*		- A*I1(X*sin(chi0)) + B*K1(X*sin(chi0)) 
!C*		+ cos(chi0) * int{y=X,infinity} [ exp(-y) 
!C*		  * ( I1(X*sin(chi0))*K0(y*sin(chi0))
!C*		    + K1(X*sin(chi0))*I0(y*sin(chi0)) ) dy ] )
!C*
!C*	Applying the boundary condition Ch(X,0) = 1 requires that 
!C*	B = 0.  Similarly, the requirement that Ch(X,chi0) approach 
!C*	the finite value of sec(chi0) as X approaches infinity [Ch31b, 
!C*	p486, Equation (12)] implies A = 0.  Thus we have
!C*
!C*	  Ch(X,chi0) = X*sin(chi0)*cos(chi0)*exp(X)*
!C*		int{y=X,infinity} [ exp(-y) 
!C*		  * ( I1(X*sin(chi0))*K0(y*sin(chi0))
!C*		    + K1(X*sin(chi0))*I0(y*sin(chi0)) ) dy ]
!C*
!C*	The function atm8_chap_deq(X,chi0) evaluates this expression.
!C*	Since explicit approximations are available for I1(z) and K1(z),
!C*	the remaining challenge is evaluation of the integrals
!C*
!C*	  int{y=X,infinity} [ exp(-y) I0(y*sin(chi0)) dy ]
!C*
!C*	and
!C*
!C*	  int{y=X,infinity} [ exp(-y) K0(y*sin(chi0)) dy ]
!C*
!C*	which are accomplished by term-by-term integration of ascending
!C*	and descending power series expansions of I0(z) and K0(z).
!C*
!C*  REFERENCES:
!C*
!C*	AS64	M. Abramowitz and I. A. Stegun, "Handbook of 
!C*		Mathematical Functions," NBS AMS 55 (USGPO, 
!C*		Washington, DC, June 1964, 9th printing, November 1970).
!C*
!C*	Ch31b	S. Chapman, "The Absorption and Dissociative or
!C*		Ionizing Effect of Monochromatic Radiation in an
!C*		Atmosphere on a Rotating Earth: Part II. Grazing
!C*		Incidence," Proc. Phys. Soc. (London), _43_, 483-501 
!C*		(1931).
!C*
!C*	Fi64	J. A. Fitzmaurice, "Simplfication of the Chapman
!C*		Function for Atmospheric Attenuation," Appl. Opt. _3_,
!C*		640 (1964).
!C*
!C*	Hi56a	F. B. Hildebrand, "Introduction to Numerical
!C*		Analysis," (McGraw-Hill, New York, 1956).
!C*
!C*	Hu39	E. O. Hulburt, "The E Region of the Ionosphere," 
!C*		Phys. Rev. _55_, 639-645 (1939).
!C*
!C*	PFT86	W. H. Press, B. P. Flannery, S. A. Teukolsky, and 
!C*		W. T. Vetterling, "Numerical Recipes," (Cambridge, 
!C*		1986).
!C*
!C*	Ra66	A. L. Rabenstein, "Introduction to Ordinary
!C*		Differential Equations," (Academic, NY, 1966).
!C*
!C*	Re89	M. H. Rees, "Physics and Chemistry of the Upper
!C*		Atmosphere," (Cambridge, 1989).
!C*
!C*	SG69	W. Swider, Jr., and M. E. Gardner, "On the Accuracy 
!C*		of Chapman FUNCTION Approximations," Appl. Opt. _8_,
!C*		725 (1969).
!C*
!C*	Sw64	W. Swider, Jr., "The Determination of the Optical 
!C*		Depth at Large Solar Zenith Angles," Planet. Space 
!C*		Sci. _12_, 761-782 (1964).
!C
!C  ####################################################################
!C
!C	Chapman FUNCTION calculated by various methods
!C
!C	  Ch(X,chi0) = atm_chapman(X,chi0)   : REAL*4 entry
!C	  Ch(X,chi0) = atm8_chapman(X,chi0)  : REAL*8 entry
!C
!C	Internal service routines - user should not call, except for
!C	testing.
!C
!C	  Ch(X,chi0) = atm8_chap_asy(X,chi0) : asymptotic expansion
!C	  Ch(X,chi0) = atm8_chap_deq(X,chi0) : differential equation
!C	  Ch(X,chi0) = atm_chap_num(X,chi0)  : REAL*4 numerical integral
!C	  Ch(X,chi0) = atm8_chap_num(X,chi0) : REAL*8 numerical integral
!C
!C  ####################################################################

!C  ====================================================================
!C
!C	These are the entries for the user to call.
!C
!C	chi0 can range from 0 to 180 in degrees.  For chi0 .gt. 90, the 
!C	product X*(1-sin(chi0)) must not be too large, otherwise we 
!C	will get an exponential overflow.
!C
!C	For chi0 <= 90 degrees, X can range from 0 to thousands 
!C	without overflow.
!C
!C  ====================================================================

	REAL FUNCTION atm_chapman( X, chi0 )
	atm_chapman = atm8_chapman( DBLE(X), DBLE(chi0) )
	RETURN
END FUNCTION atm_chapman

!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chapman( X, chi0 )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (rad=57.2957795130823208768d0)

	IF ( (X <= 0) .OR. (chi0 <= 0) .OR. (chi0 >= 180) ) THEN
	  atm8_chapman = 1
	  RETURN
	ENDIF

	IF ( chi0 > 90 ) THEN
	  chi = 180 - chi0
	ELSE
	  chi = chi0
	ENDIF

	IF ( X < 36 ) THEN
	  atm8_chapman = atm8_chap_deq(X,chi)
	ELSE
	  atm8_chapman = atm8_chap_asy(X,chi)
	ENDIF

	IF ( chi0 > 90 ) THEN
	  atm8_chapman = 2*EXP(X*2*SIN((90-chi)/(2*rad))**2) &
     		* atm8_chap_xK1(X*SIN(chi/rad)) - atm8_chapman
	ENDIF

	RETURN
END FUNCTION atm8_chapman

!C  ====================================================================
!C
!C	This Chapman FUNCTION routine calculates
!C
!C	  Ch(X,chi0) = atm8_chap_asy(X,chi0)
!C		     = sum{n=0,3} [C(n) * int{t=t0,infinity} 
!C			[ EXP(-X*(t**2-t0**2) * (t**2-t0**2)**n dy ] ]
!C
!C	with t0**2 = 1 - SIN(chi0)
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_asy( X, chi0 )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (rad=57.2957795130823208768d0)
	DIMENSION C(0:3), XI(0:3), Dn(0:3)
	COMMON/atm8_chap_cm/Fn(0:3)

	IF ( (X <= 0) .OR. (chi0 <= 0) ) THEN
	  DO i=0,3
	    Fn(i) = 1
	  ENDDO
        ELSE
   	   sinchi = SIN(chi0/rad)
	   s1 = 1 + sinchi
	   rx = sqrt(X)
	   Y0 = rx * sqrt( 2*SIN( (90-chi0)/(2*rad) )**2 )

	   C(0) = 1/sqrt(s1)
	   fact = C(0)/s1
	   C(1) = fact * (0.5d0+sinchi)
	   fact = fact/s1
	   C(2) = - fact * (0.125d0+0.5d0*sinchi)
	   fact = fact/s1
	   C(3) = fact * (0.0625d0+0.375d0*sinchi)

	   CALL atm8_chap_gd3( Y0, Dn )
	   fact = 2*rx
	   DO n=0,3
	     XI(n) = fact * Dn(n)
	     fact = fact/X
	   ENDDO

	   Fn(0) = C(0) * XI(0)
	   DO i=1,3
	     Fn(i) = Fn(i-1) + C(i)*XI(i)
	   ENDDO
        ENDIF
	atm8_chap_asy = Fn(3)
	RETURN
END FUNCTION atm8_chap_asy

!C  ====================================================================
!C
!C	This Chapman FUNCTION routine calculates
!C
!C	  Ch(X,chi0) = atm8_chap_deq(X,chi0)
!C		     = X * SIN(chi0) * cos(chi0) * EXP(X*SIN(chi0))
!C		       * int{y=X,infinity} [ EXP(-y)*( 
!C			 I1(X*SIN(chi0))*K0(y*SIN(chi0)) 
!C			 + K1(X*SIN(chi0))*I0(y*SIN(chi0)) ) dy ]
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_deq( X, chi0 )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (rad=57.2957795130823208768d0)
	COMMON/atm8_chap_cm/xI1,xK1,yI0,yK0

	IF ( (X <= 0) .OR. (chi0 <= 0) ) THEN 
           atm8_chap_deq = 1
        ELSE
	   alpha = X * SIN(chi0/rad)

!C  --------------------------------------------------------------------
!C
!C	This code fragment calculates
!C
!C	  yI0 = EXP(x*(1-SIN(chi0))) * cos(chi0) * 
!C		int{y=x,infinity} [ EXP(-y) * I0(y*SIN(chi0)) dy ]
!C
!C  --------------------------------------------------------------------

	   yI0 = atm8_chap_yI0( X, chi0 )

!C  --------------------------------------------------------------------
!C
!C	This code fragment calculates
!C
!C	  yK0 = EXP(x*(1+SIN(chi0))) x * SIN(chi0) * cos(chi0) * 
!C		int{y=x,infinity} [ EXP(-y) * K0(y*SIN(chi0)) dy ]
!C
!C  --------------------------------------------------------------------

	   yK0 = atm8_chap_yK0( X, chi0 )

!C  --------------------------------------------------------------------
!C
!C	This code fragment calculates
!C
!C	  xI1 = EXP(-x*SIN(chi0)) * I1(x*SIN(chi0))
!C
!C  --------------------------------------------------------------------

	   xI1 = atm8_chap_xI1( alpha )

!C  --------------------------------------------------------------------
!C
!C	This code fragment calculates
!C
!C	  xK1 = x*SIN(chi0) * EXP(x*SIN(chi0)) * K1(x*SIN(chi0))
!C
!C  --------------------------------------------------------------------

	   xK1 = atm8_chap_xK1( alpha )

!C  --------------------------------------------------------------------
!C
!C	Combine the terms
!C
!C  --------------------------------------------------------------------

	   atm8_chap_deq = xI1*yK0 + xK1*yI0
        ENDIF

	RETURN
END FUNCTION atm8_chap_deq

!C  ====================================================================
!C
!C	This Chapman FUNCTION routine calculates
!C
!C	  Ch(X,chi0) = atm_chap_num(X,chi0) = numerical integral
!C
!C  ====================================================================

	REAL*4 FUNCTION atm_chap_num(X,chi0)
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	REAL*4 X, chi0
	PARAMETER (rad=57.2957795130823208768D0)
	PARAMETER (n=65,nfact=8)
	DIMENSION factor(0:nfact)
	DATA factor/14175.0D0, 23552.0D0, -3712.0D0, 41984.0D0, &
     	  -18160.0D0, 41984.0D0, -3712.0D0, 23552.0D0, 7912.0D0/

	IF ( (chi0 <= 0) .OR. (chi0 > 90) .OR. (X <= 0) ) THEN
	  atm_chap_num = 1
	  RETURN
	ENDIF

	X8 = X
	chi0rad = chi0/rad
	sinchi = SIN(chi0rad)

	alpha0 = aSIN( (X8/(X8+28)) * sinchi )
	delta = (chi0rad - alpha0)/(n-1)

	sum_total = 0

	DO i=1,n
	  alpha = -(i-1)*delta + chi0rad

	  IF ( (i .EQ. 1) .OR. (X <= 0) ) THEN
	    f = 1/(1+COS(alpha))
	  ELSE IF ( alpha <= 0 ) THEN
	    f = 0
	  ELSE
	    f = EXP(-X8*(sinchi/SIN(alpha)-1) ) /(1+COS(alpha))
	  ENDIF

	  IF ( (i.EQ.1) .OR. (i.EQ.n) ) THEN
	    fact = factor(nfact)/2
	  ELSE
	    fact = factor( mod(i-2,nfact)+1 )
	  ENDIF

	  sum_total = sum_total + fact*f
	ENDDO

	atm_chap_num = 1 + X8*sinchi*sum_total*delta/factor(0)
	RETURN
END FUNCTION atm_chap_num

!C  ====================================================================
!C
!C	This Chapman FUNCTION routine calculates
!C
!C	  Ch(X,chi0) = atm8_chap_num(X,chi0) = numerical integral
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_num(X,chi0)
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (rad=57.2957795130823208768D0)
	PARAMETER (n=601,nfact=8)
	DIMENSION factor(0:nfact)
	DATA factor/14175.0D0, 23552.0D0, -3712.0D0, 41984.0D0, &
     	  -18160.0D0, 41984.0D0, -3712.0D0, 23552.0D0, 7912.0D0/

	IF ( (chi0 <= 0) .OR. (chi0 > 90) .OR. (X <= 0) ) THEN
	  atm8_chap_num = 1
	  RETURN
	ENDIF

	chi0rad = chi0/rad
	sinchi = SIN(chi0rad)

	alpha0 = aSIN( (X/(X+45)) * sinchi )
	delta = (chi0rad - alpha0)/(n-1)

	sum_total = 0

	DO i=1,n
	  alpha = -(i-1)*delta + chi0rad

	  IF ( (i .EQ. 1) .OR. (X <= 0) ) THEN
	    f = 1/(1+COS(alpha))
	  ELSE IF ( alpha <= 0 ) THEN
	    f = 0
	  ELSE
	    f = EXP(-X*(sinchi/SIN(alpha)-1) ) /(1+COS(alpha))
	  ENDIF

	  IF ( (i.EQ.1) .OR. (i.EQ.n) ) THEN
	    fact = factor(nfact)/2
	  ELSE
	    fact = factor( mod(i-2,nfact)+1 )
	  ENDIF

	  sum_total = sum_total + fact*f
	ENDDO

	atm8_chap_num = 1 + X*sinchi*sum_total*delta/factor(0)
	RETURN
END FUNCTION atm8_chap_num

!C  ####################################################################
!C
!C	The following "Bessel integral" routines RETURN various 
!C	combinations of integrals of Bessel FUNCTIONs, powers, 
!C	and EXPonentials, involving trigonometric FUNCTIONs of chi0.
!C
!C	For small values of z = X*SIN(chi0) we EXPand
!C
!C	  I0(z) = sum{n=0,6} [ aI0(n) * z**(2*n) ]
!C	  K0(z) = -log(z)*I0(z) + sum{n=0,6} [ aK0(n) * z**(2*n) ]
!C
!C	For large values of z we EXPand in reciprocal powers
!C
!C	  I0(z) = EXP(z) * sum{n=0,8} [ bI0(n) * z**(-n-0.5) ]
!C	  K0(z) = EXP(-z) * sum{n=0,6} [ bK0(n) * z**(-n-0.5) ]
!C
!C	The EXPansion coefficients are calculated from those given 
!C	by Abramowitz and Stegun [AS64, pp378-9, Section 9.8] and
!C	Press et al. [PFT86, pp177-8, BESSI0.FOR, BESSK0.FOR].
!C
!C	For small values of X*SIN(chi0) we break the integral
!C	into two parts (with F(z) = I0(z) or K0(z)):
!C
!C	  int{y=X,infinity} [ EXP(-y) * F(y*SIN(chi0)) dy ]
!C
!C	    = int{y=X,x1} [ EXP(-y) * F(y*SIN(chi0)) dy ]
!C	      + int{y=x1,infinity} [ EXP(-y) * F(y*SIN(chi0)) dy ]
!C
!C	where x1 = 3.75/SIN(chi0) for I0 and 2/SIN(chi0) for K0.
!C
!C	In the range y=X,x1 we integrate the term-by-term uSINg
!C
!C	  int{z=a,b} [ EXP(-z) * z**(2*n) dz ]
!C	    = Gamma(2*n+1,a) - Gamma(2*n+1,b)
!C
!C	and a similar but more complicated formula for
!C
!C	  int{z=a,b} [ log(z) * EXP(-z) * z**(2*n) dz ]
!C
!C	In the range y=x1,infinity we use
!C
!C	  int{z=b,infinity} [ EXP(-z) * z**(-n-0.5) dz]
!C	    = Gamma(-n+0.5,b)
!C
!C  ####################################################################

!C  ====================================================================
!C
!C	This Bessel integral routine calculates
!C
!C	  yI0 = EXP(X*(1-SIN(chi0))) * cos(chi0) * 
!C		int{y=X,infinity} [ EXP(-y) * I0(y*SIN(chi0)) dy ]
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_yI0( X, chi0 )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (rad=57.2957795130823208768d0)
	DIMENSION qbeta(0:8), gg(0:6)
	DIMENSION aI0(0:6), bI0(0:8)

        DATA aI0/ 1.0000000D+00, 2.4999985D-01, 1.5625190D-02, &
           4.3393973D-04, 6.8012343D-06, 6.5601736D-08,        &
           5.9239791D-10/
        DATA bI0/ 3.9894228D-01, 4.9822200D-02, 3.1685484D-02, &
          -8.3090918D-02, 1.8119815D+00,-1.5259477D+01,        &
           7.3292025D+01,-1.7182223D+02, 1.5344533D+02/

	theta = (90-chi0)/(2*rad)
	sint = SIN(theta)
	cost = COS(theta)
	sinchi = SIN(chi0/rad)
	coschi = COS(chi0/rad)
	sc1m = 2*sint**2	! = (1-sinchi)

	alpha = X * sinchi

	IF ( alpha <= 0 ) THEN
	  atm8_chap_yI0 = 1
	ELSE IF ( alpha < 3.75d0 ) THEN
	  x1 = 3.75d0/sinchi
	  CALL atm8_chap_gg06( X, x1, gg )
	  IF ( X <= 1 ) THEN
	    rho = 1
	  ELSE
	    rho = 1/X
	  ENDIF
	  f = (sinchi/rho)**2
	  sum_total = aI0(6)*gg(6)
	  DO i=5,0,-1
	    sum_total = sum_total*f + aI0(i)*gg(i)
!C	    write(*,1900)i,sum,gg(i)
!C1900	format(i5,1p5d14.6)
	  ENDDO
	  CALL atm8_chap_gq85( x1*sc1m, qbeta )
	  sum2 = bI0(8) * qbeta(8)
	  DO n=7,0,-1
	    sum2 = sum2/3.75d0 + bI0(n)*qbeta(n)
	  ENDDO
	  atm8_chap_yI0 = EXP(-alpha)*coschi*sum_total &
     		+ EXP((X-x1)*sc1m)*sum2*cost*sqrt(2/sinchi)
	ELSE
	  CALL atm8_chap_gq85( X*sc1m, qbeta )
	  sum_total = bI0(8) * qbeta(8)
	  DO n=7,0,-1
	    sum_total = sum_total/alpha + bI0(n)*qbeta(n)
	  ENDDO
	  atm8_chap_yI0 = sum_total * cost * sqrt( 2 / sinchi )
	ENDIF
	RETURN
END FUNCTION atm8_chap_yI0

!C  ====================================================================
!C
!C	This Bessel integral routine calculates
!C
!C	  yK0 = EXP(x*(1+SIN(chi0))) x * SIN(chi0) * cos(chi0) * 
!C		int{y=x,infinity} [ EXP(-y) * K0(y*SIN(chi0)) dy ]
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_yK0( x, chi0 )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (rad=57.2957795130823208768d0)
	DIMENSION aI0(0:6), aK0(0:6), bK0(0:6)
	DIMENSION gf(0:6), gg(0:6), qgamma(0:8)
	
        DATA aI0/ 1.0000000D+00, 2.4999985D-01, 1.5625190D-02, &
           4.3393973D-04, 6.8012343D-06, 6.5601736D-08,        &
           5.9239791D-10/
        DATA aK0/ 1.1593152D-01, 2.7898274D-01, 2.5249154D-02, &
           8.4587629D-04, 1.4975897D-05, 1.5045213D-07,        &
           2.2172596D-09/
        DATA bK0/ 1.2533141D+00,-1.5664716D-01, 8.7582720D-02, &
          -8.4995680D-02, 9.4059520D-02,-8.0492800D-02,        &
           3.4053120D-02/

	theta = (90-chi0)/(2*rad)
	sint = SIN(theta)
	cost = COS(theta)
	sinchi = SIN(chi0/rad)
	sc1 = 1+sinchi
	coschi = SIN(2*theta)

	alpha = X * sinchi
	gamma = X * sc1

	IF ( alpha <= 0 ) THEN
	  atm8_chap_yK0 = 0
	ELSE IF ( alpha < 2 ) THEN
	  x1 = 2/sinchi
	  CALL atm8_chap_gfg06( X, x1, gf, gg )
	  IF ( x <= 1 ) THEN
	    rho = 1
	  ELSE
	    rho = 1/X
	  ENDIF
	  sl = LOG(sinchi)
	  f = (sinchi/rho)**2
	  sum_total = -aI0(6)*gf(6) + (-sl*aI0(6)+aK0(6))*gg(6)
	  DO i=5,0,-1
	    sum_total = sum_total*f - aI0(i)*gf(i) + (-sl*aI0(i)+aK0(i))*gg(i)
!C	    write(*,1900)i,sum,gf(i),gg(i)
!C1900	format(i5,1p5d14.6)
	  ENDDO
	  CALL atm8_chap_gq85( x1*sc1, qgamma )
	  sum2 = bK0(6)*qgamma(6)
	  DO i=5,0,-1
	    sum2 = sum2*0.5d0 + bK0(i)*qgamma(i)
!C	    write(*,1900)i,sum2,bK0(i),qgamma(i)
	  ENDDO
	  sum_total = sum_total + EXP(X-x1-2)*sum2/sqrt(sinchi*sc1)
	  atm8_chap_yK0 = sum_total * EXP(alpha) * alpha * coschi
	ELSE
	  CALL atm8_chap_gq85( gamma, qgamma )
	  sum_total = bK0(6) * qgamma(6)
	  DO i=5,0,-1
	    sum_total = sum_total/alpha + bK0(i)*qgamma(i)
	  ENDDO
	  atm8_chap_yK0 = sum_total * sint * sqrt( 2 * sinchi ) * X
	ENDIF

	RETURN
END FUNCTION atm8_chap_yK0

!C  ####################################################################
!C
!C	The following "pure math" routines RETURN various combinations
!C	of Bessel FUNCTIONs, powers, and EXPonentials.
!C
!C  ####################################################################

!C  ====================================================================
!C
!C	This Bessel FUNCTION math routine RETURNs
!C
!C	  xI1 = EXP(-|z|) * I1(z)
!C
!C	Following Press et al [PFT86, page 178, BESSI1.FOR] and 
!C	Abrahamson and Stegun [AS64, page 378, 9.8.3, 9.8.4].
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_xI1( z )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
        DIMENSION aI1(0:6), bI1(0:8)

        DATA aI1/ 5.00000000D-01, 6.2499978D-02, 2.6041897D-03, &
           5.4244512D-05, 6.7986797D-07, 5.4830314D-09,         &
           4.1909957D-11/
        DATA bI1/ 3.98942280D-01,-1.4955090D-01,-5.0908781D-02, &
           8.6379434D-02,-2.0399403D+00, 1.6929962D+01,         &
          -8.0516146D+01, 1.8642422D+02,-1.6427082D+02/

	IF ( z < 0 ) THEN
	  az = -z
	ELSE IF ( z .EQ. 0 ) THEN
	  atm8_chap_xI1 = 0
	  RETURN
	ELSE
	  az = z
	ENDIF
	IF ( az < 3.75d0 ) THEN
	  z2 = z*z
	  sum_total = aI1(6)
	  DO i=5,0,-1
	    sum_total = sum_total*z2 + aI1(i)
	  ENDDO
	  atm8_chap_xI1 = z*EXP(-az) * sum_total
	ELSE
	  sum_total = bI1(8)
	  DO i=7,0,-1
	    sum_total = sum_total/az + bI1(i)
	  ENDDO
	  atm8_chap_xI1 = sum_total*sqrt(az)/z
	ENDIF
	RETURN
END FUNCTION atm8_chap_xI1

!C  ====================================================================
!C
!C	This Bessel FUNCTION math routine RETURNs
!C
!C	  xK1 = z * EXP(+z) * K1(z)
!C
!C	Following Press et al [PFT86, page 179, BESSK1.FOR] and 
!C	Abrahamson and Stegun [AS64, page 379, 9.8.7, 9.8.8].
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_xK1( z )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
        DIMENSION aK1(0:6), bK1(0:6)

        DATA aK1/ 1.00000000D+00, 3.8607860D-02,-4.2049112D-02, &
          -2.8370152D-03,-7.4976641D-05,-1.0781641D-06,         &
          -1.1440430D-08/
        DATA bK1/ 1.25331414D+00, 4.6997238D-01,-1.4622480D-01, &
           1.2034144D-01,-1.2485648D-01, 1.0419648D-01,         &
          -4.3676800D-02/

	IF ( z <= 0 ) THEN
	  atm8_chap_xK1 = 1
	ELSE IF ( z < 2 ) THEN
	  xz = EXP(z)
	  z2 = z*z
	  sum_total = aK1(6)
	  DO i=5,0,-1
	    sum_total = sum_total*z2 + aK1(i)
	  ENDDO
	  atm8_chap_xK1 = xz * ( sum_total &
     		+ z*LOG(z/2)*atm8_chap_xI1(z)*xz )
	ELSE
	  sum_total = bk1(6)
	  DO i=5,0,-1
	    sum_total = sum_total/z + bK1(i)
	  ENDDO
	  atm8_chap_xK1 = sum_total*sqrt(z)
	ENDIF

	RETURN
END FUNCTION atm8_chap_xK1

!C  ####################################################################
!C
!C	The following "pure math" routines RETURN various combinations
!C	of the Error FUNCTION, powers, and EXPonentials.
!C
!C  ####################################################################

!C  ====================================================================
!C
!C	This Error FUNCTION math routine RETURNs
!C
!C	  xerfc(x) = EXP(x**2)*erfc(x)
!C
!C	following Press et al. [PFT86, p164, ERFCC.FOR]
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_xerfc(x)
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
        T=1.0D0/(1.0D0+0.5D0*x)
	atm8_chap_xerfc = &
     	  T*EXP( -1.26551223D0 +T*(1.00002368D0 +T*( .37409196D0 &
            +T*(  .09678418D0 +T*(-.18628806D0 +T*( .27886807D0 &
     	     +T*(-1.13520398D0 +T*(1.48851587D0 +T*(-.82215223D0 &
     	     +T*   .17087277D0) ))))))))
        RETURN
      END FUNCTION atm8_chap_xerfc

!C  ####################################################################
!C
!C	The following "pure math" routines RETURN various combinations
!C	of EXPonential integrals, powers, and EXPonentials.
!C
!C  ####################################################################

!C  ====================================================================
!C
!C	This EXPonential math routine evaluates
!C
!C	  zxE1(x) = x*EXP(x) int{y=1,infinity} [ EXP(-x*y)/y dy ]
!C
!C	following Abramowitz and Stegun [AS64, p229;231, equations
!C	5.1.11 and 5.1.56]
!C
!C  ====================================================================

	REAL(kind(1.d0)) FUNCTION atm8_chap_zxE1(x)
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (gamma = 0.5772156649015328606d0)
	DIMENSION aE1(0:4), bE1(0:4), cEin(1:10)

	DATA aE1/1.0d0, 8.5733287401d0, 18.0590169730d0, &
     	    8.6347608925d0, 0.2677737343d0 /
	DATA bE1/1.0d0, 9.5733223454d0, 25.6329561486d0, &
     	    21.0996530827d0, 3.9584969228d0/
        DATA cEin/ 1.00000000D+00,-2.50000000D-01, 5.55555556D-02, &
         -1.0416666667D-02, 1.6666666667D-03,-2.3148148148D-04,    &
          2.8344671202D-05,-3.1001984127D-06, 3.0619243582D-07,    &
         -2.7557319224D-08/

	IF ( x <= 0 ) THEN
	  atm8_chap_zxE1 = 0
	ELSE IF ( x <= 1 ) THEN
	  sum_total = cEin(10)
	  DO i=9,1,-1
	    sum_total = sum_total*x + cEin(i)
	  ENDDO
	  atm8_chap_zxE1 = x*EXP(x)*( x * sum_total - LOG(x) - gamma )
	ELSE
	  top = aE1(4)
	  bot = bE1(4)
	  DO i=3,0,-1
	    top = top/x + aE1(i)
	    bot = bot/x + bE1(i)
	  ENDDO
	  atm8_chap_zxE1 = top/bot
	ENDIF
	RETURN
END FUNCTION atm8_chap_zxE1

!C  ####################################################################
!C
!C	The following "pure math" routines RETURN various combinations
!C	of incomplete gamma FUNCTIONs, powers, and EXPonentials.
!C
!C  ####################################################################

!C  ====================================================================
!C
!C	This gamma FUNCTION math routine calculates
!C
!C	Dn(n) = int{t=z,infinity}
!C		[ EXP( -(t**2-z**2) ) * (t**2-z**2)**n dt ]
!C
!C  ====================================================================

	SUBROUTINE atm8_chap_gd3( z, Dn )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (rpi=1.7724538509055160273d0)
	DIMENSION Dn(0:3), xg(0:3)

	IF ( z <= 0 ) THEN
	  Dn(0) = rpi/2
	  DO i=1,3
	    Dn(i) = (i-0.5d0)*Dn(i-1)
	  ENDDO
	  RETURN
	ENDIF

	z2 = z*z
	IF ( z >= 7 ) r = 1/z2

	IF ( z < 14 ) THEN
	  z4 = z2*z2
	  xg(0) = rpi * atm8_chap_xerfc(z)
	  xg(1) = 0.5d0*xg(0) + z
	  xg(2) = 1.5d0*xg(1) + z*z2
	  Dn(0) = 0.5d0*xg(0)
	  Dn(1) = 0.5d0*(xg(1)-z2*xg(0))
	  Dn(2) = 0.5d0*(xg(2)-2*z2*xg(1)+z4*xg(0))
	ELSE
	  Dn(0) = ( 1 + r*(-0.5d0 +r*(0.75d0 +r*(-1.875d0 &
     		+r*6.5625d0) ) ) )/(2*z)
	  Dn(1) = ( 1 + r*(-1.0d0 +r*(2.25d0 +r*(-7.5d0   &
     		+r*32.8125d0) ) ) )/(2*z)
	  Dn(2) = ( 2 + r*(-3.0d0 +r*(9.00d0 +r*(-37.5d0  &
     		+r*196.875d0) ) ) )/(2*z)
	ENDIF

	IF ( z < 7 ) THEN
	  z6 = z4*z2
	  xg(3) = 2.5d0*xg(2) + z*z4
	  Dn(3) = 0.5d0*(xg(3)-3*z2*xg(2)+3*z4*xg(1)-z6*xg(0))
	ELSE
	  Dn(3) = ( 6 + r*(-12.0d0 +r*(45.0d0 +r*(-225.0d0 &
     		+r*1378.125d0) ) ) )/(2*z)
	ENDIF

	RETURN
END SUBROUTINE atm8_chap_gd3

!C  ====================================================================
!C
!C	This Gamma FUNCTION math routine calculates
!C
!C	  gf06(n) = g(n,x) * int{y=x,z} [LOG(y) * EXP(-y) * y**(2*n) dy]
!C
!C	and
!C
!C	  gg06(n) = g(n,x) * int{y=x,z} [ EXP(-y) * y**(2*n) dy ]
!C	          = g(n,x) * ( Gamma(2*n+1,x) - Gamma(2*n+1,z) )
!C
!C	for n=0,6, with g(n,x) = EXP(x) * max(1,x)**(-2*n)
!C
!!C  ====================================================================

	SUBROUTINE atm8_chap_gfg06( x, z, gf06, gg06 )
	IMPLICIT REAL(kind(1.d0)) (a-h,o-z)
	PARAMETER (gamma = 0.5772156649015328606d0)
	DIMENSION gf06(0:6), gg06(0:6)
	DIMENSION gh13x(13), gh13z(13), rgn(13), delta(13)
	CALL atm8_chap_gh13( x, x, gh13x )
	CALL atm8_chap_gh13( x, z, gh13z )
	IF ( x <= 1 ) THEN
	  rho = 1
	ELSE
	  rho = 1/x
	ENDIF

	delta(1) = 0
	delta(2) = ( gh13x(1) - gh13z(1) ) * rho
	rgn(1) = 1
	rgn(2) = rho
	DO n=2,12
	  delta(n+1) = rho*( n*delta(n) + gh13x(n) - gh13z(n) )
	  rgn(n+1) = (n*rho)*rgn(n)
	ENDDO

	IF ( x > 0 ) THEN
	  xE1_x = atm8_chap_zxE1(x)/x
	  xLOG = LOG(x)
	ENDIF
	IF ( z > 0 ) THEN
	  xE1_z = EXP(x-z)*atm8_chap_zxE1(z)/z
	  zLOG = LOG(z)
	ENDIF

	DO k=0,6
	  n = 2*k+1
	  IF ( x <= 0 ) THEN
	    gf06(k) = -gamma*rgn(n) + delta(n)
	  ELSE
	    gf06(k) = xLOG*gh13x(n) + rgn(n)*xE1_x + delta(n)
	  ENDIF
	  IF ( z <= 0 ) THEN
	    gf06(k) = gf06(k) + gamma*rgn(n)
	  ELSE
	    gf06(k) = gf06(k) - (zLOG*gh13z(n) + rgn(n)*xE1_z)
	  ENDIF
	  gg06(k) = gh13x(n) - gh13z(n)
	ENDDO

	RETURN
END SUBROUTINE atm8_chap_gfg06

!C  ====================================================================
!C
!C	This Gamma FUNCTION math routine calculates
!C
!C	  gg06(n) = g(n,x) * int{y=x,z} [ EXP(-y) * y**(2*n) dy ]
!C	          = g(n,x) * ( Gamma(2*n+1,x) - Gamma(2*n+1,z) )
!C
!C	for n=0,6, with g(n,x) = EXP(x) * max(1,x)**(-2*n)
!C
!C  ====================================================================

	SUBROUTINE atm8_chap_gg06( x, z, gg06 )
	IMPLICIT REAL(kind(1.d0)) (a-h,o-z)
	DIMENSION gg06(0:6), gh13x(13), gh13z(13)
	CALL atm8_chap_gh13( x, x, gh13x )
	CALL atm8_chap_gh13( x, z, gh13z )
	DO n=0,6
	  gg06(n) = gh13x(2*n+1) - gh13z(2*n+1)
	ENDDO
	RETURN
END SUBROUTINE atm8_chap_gg06

!C  ====================================================================
!C
!C	This Gamma FUNCTION math routine calculates
!C
!C	  gh13(n) = f(n,x) * int{y=z,infinity} [EXP(-y) * y**(n-1) dy]
!C	          = f(n,x) * Gamma(n,z)
!C
!C	for n=1,13, with f(n,x) = EXP(x) * max(1,x)**(-n+1)
!C
!C  ====================================================================

	SUBROUTINE atm8_chap_gh13( x, z, gh13 )
	IMPLICIT REAL(kind(1.d0)) (a-h,o-z)
	DIMENSION gh13(13), Tab(12)

	IF ( z <= 0 ) THEN
	  gh13(1) = 1
	  DO n=1,12
	    gh13(n+1) = n*gh13(n)
	  ENDDO
	  RETURN
	ENDIF

	IF ( x <= 1 ) THEN
	  rho = 1
	ELSE
	  rho = 1/x
	ENDIF
	rhoz = rho * z
	exz = EXP(x-z)
	Tab(12) = EXP( (x-z) + 12*LOG(rhoz) )
	DO n=11,1,-1
	  Tab(n) = Tab(n+1)/rhoz
	ENDDO
	gh13(1) = exz
	DO n=1,12
	  gh13(n+1) = rho*n*gh13(n) + Tab(n)
	ENDDO
	RETURN
END SUBROUTINE atm8_chap_gh13

!C  ====================================================================
!C
!C	This Gamma FUNCTION math SUBROUTINE calculates
!C
!C	  Qn(x) = x**n * EXP(x) * Gamma(-n+0.5,x), n=0,8
!C	    = x**n * EXP(x) * int{y=x,infinity} [EXP(-y)*y**(-n-0.5)dy]
!C
!C	For x < 2 we first calculate
!C
!C	  Q0(x) = sqrt(pi)*EXP(x)*erfc(sqrt(x)) = EXP(x)*Gamma(0.5,x)
!C
!C	and use upward recursion.  ELSE, we first calculate
!C
!C	  Q8(x) = x**8 * EXP(x) * Gamma(-7.5,x)
!C
!C	following Press et al. [PFT86, pp162-63, GCF.FOR] and THEN
!C	recur downward.  Also see Abramowitz and Stegun [AS64, 6.5].
!C
!C  ====================================================================

	SUBROUTINE atm8_chap_gq85( x, qn )
	IMPLICIT REAL(kind(1.d0))(a-h,o-z)
	PARAMETER (rpi=1.7724538509055160273d0)
	PARAMETER (itmax=100,eps=3.0d-9)
	DIMENSION qn(0:8)

	IF ( x <= 0 ) THEN
	  qn(0) = rpi
	  DO i=1,8
	    qn(i) = 0
	  ENDDO
	  RETURN
	ENDIF

	rx = sqrt(x)

	IF ( x < 2 ) THEN
	  qn(0) = rpi * atm8_chap_xerfc( rx )
	  DO n=1,8
	    qn(n) = ( -rx*qn(n-1) + 1 ) * rx / ( n - 0.5d0 )
	  ENDDO
	ELSE
          GOLD=0.0d0
	  A0=1.0d0
	  A1=x
	  B0=0.0d0
	  B1=1.0d0
	  FAC=1.0d0
	  DO N=1,ITMAX
	    AN= (N)
	    ANA=AN + 7.5d0
	    A0=(A1+A0*ANA)*FAC
	    B0=(B1+B0*ANA)*FAC
	    ANF=AN*FAC
	    A1=x*A0+ANF*A1
	    B1=x*B0+ANF*B1
	    FAC=1./A1
            G=B1*FAC
	    test = G*eps
	    del = G - Gold
	    IF ( test < 0 ) test = - test
	    IF ( (del >= -test) .AND. (del <= test) ) EXIT
	    GOLD=G
          ENDDO
	  qn(8) = G * rx
	  DO n=8,1,-1
	    qn(n-1) = ( (-n+0.5d0)*qn(n)/rx + 1 ) / rx
	  ENDDO
	ENDIF

	RETURN
END SUBROUTINE atm8_chap_gq85

      END MODULE module_ra_chapman
