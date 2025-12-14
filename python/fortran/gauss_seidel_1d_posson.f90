program gauss_seidel_1d_poisson
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: n, i, k, max_iter
  real(dp) :: a, b, h, pi, tol, relres
  real(dp), allocatable :: x(:), rhs(:), r(:)

  n = 200
  max_iter = 20000
  tol = 1.0d-10
  a = 0.0_dp
  b = 1.0_dp
  pi = acos(-1.0_dp)

  allocate(x(n), rhs(n), r(n))

  h = (b - a) / real(n + 1, dp)

  do i = 1, n
     rhs(i) = (pi*pi) * sin(pi * (a + real(i, dp)*h))
  end do

  x = 0.0_dp

  do k = 1, max_iter
     call gauss_seidel_step_1d_poisson(x, rhs, h, n)
     call residual_1d_poisson(x, rhs, h, n, r)
     relres = norm2(r) / max(1.0_dp, norm2(rhs))
     if (relres < tol) exit
  end do

  print *, "Gauss-Seidel - 1D Poisson (-u''=f), Dirichlet 0"
  print *, "N =", n
  print *, "Iterations =", k
  print *, "Relative residual =", relres

contains

  subroutine gauss_seidel_step_1d_poisson(x, rhs, h, n)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: h
    real(dp), intent(in) :: rhs(n)
    real(dp), intent(inout) :: x(n)
    integer :: i

    do i = 1, n
       if (i == 1) then
          x(i) = 0.5_dp * (0.0_dp + x(i+1) + h*h*rhs(i))
       else if (i == n) then
          x(i) = 0.5_dp * (x(i-1) + 0.0_dp + h*h*rhs(i))
       else
          x(i) = 0.5_dp * (x(i-1) + x(i+1) + h*h*rhs(i))
       end if
    end do
  end subroutine gauss_seidel_step_1d_poisson

  subroutine residual_1d_poisson(x, rhs, h, n, r)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: h
    real(dp), intent(in) :: x(n), rhs(n)
    real(dp), intent(out) :: r(n)
    integer :: i
    real(dp) :: ax

    do i = 1, n
       if (i == 1) then
          ax = (2.0_dp*x(i) - 0.0_dp - x(i+1)) / (h*h)
       else if (i == n) then
          ax = (2.0_dp*x(i) - x(i-1) - 0.0_dp) / (h*h)
       else
          ax = (2.0_dp*x(i) - x(i-1) - x(i+1)) / (h*h)
       end if
       r(i) = rhs(i) - ax
    end do
  end subroutine residual_1d_poisson

  pure real(dp) function norm2(v)
    implicit none
    real(dp), intent(in) :: v(:)
    integer :: i
    norm2 = 0.0_dp
    do i = 1, size(v)
       norm2 = norm2 + v(i)*v(i)
    end do
    norm2 = sqrt(norm2)
  end function norm2

end program gauss_seidel_1d_poisson
