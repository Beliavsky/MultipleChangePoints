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

  integer, parameter :: n = 100
  integer, parameter :: k_true = 2       ! True number of change points in the simulation
  integer, parameter :: k_max = 4        ! Maximum number of change points to test
  real(kind=dp) :: x(n)
  integer :: i, k

  integer :: cp_true(k_true) = [30, 70]
  real(kind=dp), allocatable :: sigma(:)
  ! For storing the estimated change points for each model with k>0
  integer, allocatable :: cp_est(:)
  real(kind=dp) :: best_ll, ll, sumsq
  integer :: nseg
  sigma = [1.0_dp, 5.0_dp, 0.2_dp]

  do i = 1, n
     if (i <= cp_true(1)) then
        x(i) = sigma(1) * random_normal()
     else if (i <= cp_true(2)) then
        x(i) = sigma(2) * random_normal()
     else
        x(i) = sigma(3) * random_normal()
     end if
  end do

  ! Print the true change points.
  print *
  print *, "True Change Points (indices):"
  do i = 1, k_true
     print *, cp_true(i)
  end do
  print *

  ! Print table header.
  print *, " k   Log-Likelihood     Detected Change Points"
  print *, "-----------------------------------------------------"

  ! Loop over models with k = 0 to k_max change points.
  do k = 0, k_max
     if (k == 0) then
        ! No change points: the model has one segment spanning the whole data.
        nseg = n
        sumsq = sum(x(1:n)**2)
        if (sumsq <= 0.0_dp) sumsq = 1.0e-10_dp
        ll = -0.5_dp * nseg * (1.0_dp + log(sumsq/nseg))
        write(*,'(I2,3X,F12.4,5X,A)') k, ll, "None"
     else
        allocate(cp_est(k))
        call find_change_points(x, cp_est, best_ll)
        write(*,"(I2,3X,F12.4,3X)", advance="no") k, best_ll
        ! Print each detected change point in the candidate.
        do i = 1, k
           write(*,"(I4)", advance="no") cp_est(i)
        end do
        write(*,*)
        deallocate(cp_est)
     end if
  end do

end program xsegmentation
