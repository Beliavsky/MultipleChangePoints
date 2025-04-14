program xsegmentation
!> Main program for simulating a zero-mean time series with
!> discrete changes in standard deviation and testing the change
!> point detection routine. The program simulates data with known
!> change points, then computes the maximum likelihood models with
!> k = 0 to k_max change points and prints a table of the log-
!> likelihood values and the corresponding detected change point
!> indices.
  use kind_mod        , only: dp
  use segmentation_mod, only: find_change_points
  use random_mod      , only: random_normal
  implicit none

  integer, parameter :: n = 400    ! # of observations
  integer, parameter :: k_true = 2 ! True number of change points in the simulation
  integer, parameter :: k_max = 3  ! Maximum number of change points to test
  real(kind=dp) :: x(n), xsq(n)
  integer :: i, k, nseg
  integer, parameter :: cp_true(k_true) = [30, 70]
  real(kind=dp), parameter :: sigma(k_true+1) = [1.0_dp, 5.0_dp, 0.2_dp]
  real(kind=dp), allocatable :: times(:)
  ! For storing the estimated change points for each model with k>0
  integer, allocatable :: cp_est(:)
  real(kind=dp) :: best_ll, ll, sumsq
  logical, parameter :: print_times = .true.

  do i = 1, n
     if (i <= cp_true(1)) then
        x(i) = sigma(1) * random_normal()
     else if (i <= cp_true(2)) then
        x(i) = sigma(2) * random_normal()
     else
        x(i) = sigma(3) * random_normal()
     end if
  end do
  xsq = x**2
  print "(/,a,*(1x,i0))", "True Change Points (indices):", cp_true
  print "(/,a)", " k   Log-Likelihood   Detected Change Points"
  ! Loop over models with k = 0 to k_max change points.
  if (print_times) allocate (times(0:k_max+1))
  allocate (cp_est(k_max))
  do k = 0, k_max
     if (print_times) call cpu_time(times(k))
     if (k == 0) then
        ! No change points: the model has one segment spanning the whole data.
        nseg = n
        sumsq = sum(xsq)
        if (sumsq <= 0.0_dp) sumsq = 1.0e-10_dp
        ll = -0.5_dp * nseg * (1.0_dp + log(sumsq/nseg))
        write(*,"(I2,3X,F12.4)") k, ll
     else
        call find_change_points(xsq, cp_est(:k), best_ll)
        write(*,"(I2,3X,F12.4,3X,*(i4))") k, best_ll, cp_est(:k)
     end if
  end do
  if (print_times) then
     call cpu_time(times(k))
     print "(/,a20, *(i8))", "#changepoints", (k, k=0, k_max)
     print "(a20,*(f8.4))", "times", times(1:) - times(:k_max)
  end if
end program xsegmentation
