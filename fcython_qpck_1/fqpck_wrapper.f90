module fqpck_wrapper

  use types, only: dp
  use iso_c_binding, only: c_double, c_int
  use f90_vs_numpy_arrays_5_qpck only: qags

  implicit none

contains

  subroutine c_quad(f, a, b, epsabs, epsrel, result, abserr, neval, ier) bind(c)
    real(c_double), external :: f
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_double), parameter, intent(in) :: epsabs = 0.0_dp
    real(c_double), parameter, intent(in) :: epsrel = 0.001_dp
    real(c_double), intent(out) :: result
    real(c_double), intent(in) :: abserr
    integer(c_int), intent(in) :: neval
    integer(c_int), intent(in) :: ier
    call qags ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )
  end subroutine c_quad

end module fqpck_wrapper

