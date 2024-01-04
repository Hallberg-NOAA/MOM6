!> A module with intrinsic functions that are used by MOM but are not supported
!!  by some compilers.
module MOM_intrinsic_functions

! This file is part of MOM6. See LICENSE.md for the license.

use iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit

implicit none ; private

public :: invcosh, cuberoot
public :: intrinsic_functions_unit_tests

contains

!> Evaluate the inverse cosh, either using a math library or an
!! equivalent expression
function invcosh(x)
  real, intent(in) :: x !< The argument of the inverse of cosh [nondim].  NaNs will
                        !! occur if x<1, but there is no error checking
  real :: invcosh  ! The inverse of cosh of x [nondim]

#ifdef __INTEL_COMPILER
  invcosh = acosh(x)
#else
  invcosh = log(x+sqrt(x*x-1))
#endif

end function invcosh

!> Returns the cube root of a real argument at roundoff accuracy, in a form that works properly with
!! rescaling of the argument by integer powers of 8.  If the argument is a NaN, a NaN is returned.
elemental function cuberoot(x) result(root)
  real, intent(in) :: x !< The argument of cuberoot in arbitrary units cubed [A3]
  real :: root !< The real cube root of x in arbitrary units [A]

  real :: asx ! The absolute value of x rescaled by an integer power of 8 to put it into
              ! the range from 0.125 < a <= 1.0, in ambiguous units cubed [B3]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [B C]
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [C]
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]
  integer :: ex_3 ! One third of the exponent part of x, used to rescale x to get a.
  integer :: itt

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    ex_3 = exponent(x) / 3
    asx = scale(abs(x), -3*ex_3)
    if (asx > 1.0) then
      asx = 0.125 * asx ; ex_3 = ex_3 + 1
!    elseif (asx <= 0.125) then
!      asx = 8.0 * asx ; ex_3 = ex_3 - 1
    endif

    num = (2.0 + asx)
    den = 3.0
    if (asx < 1.0) then
      ! Iteratively determine R = asx**1/3 using Newton's method, noting that in this case Newton's
      ! method converges monotonically from above and needs no bounding.  For the range of asx from
      ! 0.125 to 1.0 with a first guess of 1.0, 6 iterations suffice to converge to roundoff.
      do itt=1,9
        ! Newton's method iterates estimates as Root = Root - (Root**3 - asx) / (3.0 * Root**2).
        ! Keeping the estimates in a fractional form Root = num / den allows this calculation with
        ! fewer (or no) real divisions during the iterations before doing a single real division
        ! at the end, and it is therefore more computationally efficienty.

        ! Because successive estimates of the numerator and denominator tend to be the cube of their
        ! predecessors, the numerator and denominator need to be rescaled by division when they get
        ! too large or small to avoid overflow or underflow.  Were this being done for 32 bit reals,
        ! the values to compare with would be about 1.0e9 and 1.0e-9.
        if ((den > 1.0e75) .or. (den < 1.0e-75)) then
          num = num / den ; den = 1.0
        endif

        num_prev = num ; den_prev = den
        num = (2.0 * num_prev**3 + asx * den_prev**3)
        den = 3.0 * (den_prev * num_prev**2)

        if (num * den_prev == num_prev * den) &
          exit  ! This is an exact or converged solution, so stop.

        if (itt > 1) then
          if ((abs(num*den_prev - num_prev*den) <= 3.0*epsilon(num)*num*den_prev) .and. &
              (num * den_prev >= num_prev * den)) then
            ! Because Newton's method converges monotonically from above after the first iteration,
            ! if the answer has increased slightly (at the last bit) from the previous iteration, the
            ! answers may have converged to a roundoff-level limit cycle around an irrational solution,
            ! so take the previous estimate and stop iterating.  (The test above for whether the two
            ! solutions are within 1e-15 of each other might not actually be needed.)
            num = num_prev ; den = den_prev
            exit
          endif
        endif
      enddo
    endif
    root = sign(scale(num / den, ex_3), x)
  endif

end function cuberoot

!> Returns true if any unit test of intrinsic_functions fails, or false if they all pass.
logical function intrinsic_functions_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout

  ! Local variables
  real :: testval  ! A test value for self-consistency testing [nondim]
  logical :: fail, v

  fail = .false.
  v = verbose
  write(stdout,*) '==== MOM_intrinsic_functions: intrinsic_functions_unit_tests ==='

  fail = fail .or. Test_cuberoot(v, 1.2345678901234e9)
  fail = fail .or. Test_cuberoot(v, -9.8765432109876e-21)
  fail = fail .or. Test_cuberoot(v, 64.0)
  fail = fail .or. Test_cuberoot(v, -0.5000000000001)
  fail = fail .or. Test_cuberoot(v, 0.0)

end function intrinsic_functions_unit_tests

!> True if the cube of cuberoot(val) does not closely match val. False otherwise.
logical function Test_cuberoot(verbose, val)
  logical, intent(in) :: verbose !< If true, write results to stdout
  real, intent(in) :: val  !< The real value to test, in arbitrary units [A]
  ! Local variables
  real :: diff ! The difference between val and the cube root of its cube.

  diff = val - cuberoot(val**3)
  Test_cuberoot = (abs(diff) > 2.0e-15*abs(val))

  if (Test_cuberoot) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2," <-- FAIL")') val, diff
  elseif (verbose) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2)') val, diff

  endif
end function Test_cuberoot

end module MOM_intrinsic_functions
