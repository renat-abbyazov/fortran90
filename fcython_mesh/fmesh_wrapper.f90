module fmesh_wrapper

  use iso_c_binding, only: c_double, c_int
  use fmesh, only: meshexp

  implicit none

contains

  subroutine c_meshexp(r_min, r_max, a, N, mesh) bind(c)
    real(c_double), intent(in) :: r_min
    real(c_double), intent(in) :: r_max
    real(c_double), intent(in) :: a
    integer(c_int), intent(in) :: N
    real(c_double), intent(out) :: mesh(N)
    call meshexp(r_min, r_max, a, N, mesh)
  end subroutine c_meshexp

  ! wrap more functions here
  ! ...

end module fmesh_wrapper

