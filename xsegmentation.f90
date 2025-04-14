program xsegmentation
  use kind_mod, only: dp
  use segmentation_mod, only: find_change_points
  use random_mod,      only: random_normal
  use info_crit_mod,   only: aic, aicc, bic
  implicit none

  integer, parameter :: n = 100       ! Number of observations.
  integer, parameter :: k_true = 3    ! True number of change points.
  integer, parameter :: k_max = 4     ! Maximum number of change points to test.
  integer, parameter :: niter = 3
  real(kind=dp) :: x(n), xsq(n)
  integer :: i, isig, iter, j, k, nseg, p
  integer, parameter :: cp_true(k_true) = [25, 50, 75]
  real(kind=dp), parameter :: sigma(k_true+1) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
  real(kind=dp), allocatable :: times(:)
  ! For storing estimated change points for models with k > 0.
  integer, allocatable :: cp_est(:)
  real(kind=dp) :: best_ll, ll, sumsq
  logical, parameter :: print_times = .true.
  
  ! Arrays for storing information criteria for each model.
  real(kind=dp), dimension(0:k_max) :: ll_arr, aicc_arr, bic_arr
  integer :: k_aicc, k_bic
  real(kind=dp) :: aic_val, aicc_val, bic_val
  real(kind=dp) :: min_aicc, min_bic

  if (print_times) allocate(times(0:k_max+1))
  allocate(cp_est(k_max))

  print "(/,a30,*(i6))", "True ChangePoints (indices):", cp_true
  print "(a30,*(f6.2))", "True SD:", sigma

  ! Simulate the data.
  print "('#obs: ', i0)", n
do iter=1, niter
  do i = 1, n
     isig = k_true + 1
     do j = 1, k_true
        if (i <= cp_true(j)) then
           isig = j
           exit
        end if
     end do
     x(i) = sigma(isig) * random_normal()
  end do
  xsq = x**2
  
  print "(/,a)", " k   Log-Likelihood     AICC         BIC        Estimated-ChangePoints"
  
  ! Loop over models with k = 0 to k_max change points.
  do k = 0, k_max
     if (print_times) call cpu_time(times(k))
     if (k == 0) then
        ! Model with no change points: one segment.
        nseg = n
        sumsq = sum(xsq)
        if (sumsq <= 0.0_dp) sumsq = 1.0e-10_dp
        ll = -0.5_dp * nseg * (1.0_dp + log(sumsq/nseg))
     else
        call find_change_points(xsq, cp_est(:k), best_ll)
        ll = best_ll
     end if

     ll_arr(k) = ll
     p = k + 1   ! p is the number of free parameters (one per segment).
     aic_val  = aic(ll, p)
     aicc_val = aicc(ll, p, n)
     bic_val  = bic(ll, p, n)
     aicc_arr(k) = aicc_val
     bic_arr(k) = bic_val

     if (k == 0) then
        write(*, "(I2,3X,F12.4,3X,F10.4,3X,F10.4)") k, ll, aicc_val, bic_val
     else
        write(*, "(I2,3X,F12.4,3X,F10.4,3X,F10.4,3X,*(i4))") k, ll, aicc_val, bic_val, cp_est(:k)
     end if
  end do

  ! Determine the best model by AICC and BIC.
  min_aicc = aicc_arr(0)
  min_bic  = bic_arr(0)
  k_aicc = 0
  k_bic  = 0
  do k = 1, k_max
     if (aicc_arr(k) < min_aicc) then
        min_aicc = aicc_arr(k)
        k_aicc = k
     end if
     if (bic_arr(k) < min_bic) then
        min_bic = bic_arr(k)
        k_bic = k
     end if
  end do

  print *
  print *, "Model selected by AICC: k = ", k_aicc, " with AICC = ", min_aicc
  print *, "Model selected by BIC:  k = ", k_bic, " with BIC  = ", min_bic

  if (print_times) then
     call cpu_time(times(k))
     print "(/,a20, *(i8))", "#changepoints", (k, k=0, k_max)
     print "(a20,*(f8.4))", "times", times(1:) - times(:k_max)
  end if
end do
end program xsegmentation
