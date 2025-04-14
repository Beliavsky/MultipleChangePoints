program xsegmentation
!***********************************************************************
! Simulates multiple datasets of a zero-mean time series with a constant
! standard deviation and applies a change-point detection routine to
! estimate segmentation models with 0 to k_max change points. For each
! candidate model, the log-likelihood and information criteria (AICC, BIC,
! and CPIC) are computed, and the model with the minimum criterion value
! is selected. The program then summarizes how frequently each candidate
! number of change points is chosen by each criterion over all simulations.
!***********************************************************************
  use kind_mod,         only: dp
  use segmentation_mod, only: find_change_points
  use random_mod,       only: random_normal
  use info_crit_mod,    only: aic, aicc, bic, cpic
  implicit none

  integer, parameter :: n = 100           ! n: Number of observations.
  integer, parameter :: k_true = 2        ! k_true: True number of change points.
  integer, parameter :: k_max = 4         ! k_max: Maximum number of change points to test.
  integer, parameter :: niter = 100       ! niter: Number of simulations.
  real(kind=dp) :: x(n), xsq(n)
  integer :: i, isig, iter, j, k, nseg, p
  integer, parameter :: cp_true(k_true) = [50, 75] ! True change point indices.
  real(kind=dp), parameter :: sigma(k_true+1) = [1.0_dp, 2.0_dp, 1.0_dp]
  real(kind=dp), allocatable :: times(:)
  integer, allocatable :: cp_est(:)
  real(kind=dp) :: best_ll, ll, sumsq
  logical, parameter :: print_times = .true.

  ! Arrays to store criteria for each candidate model.
  real(kind=dp), dimension(0:k_max) :: ll_arr, aicc_arr, bic_arr, cpic_arr
  integer :: k_aicc, k_bic, k_cpic
  real(kind=dp) :: aic_val, aicc_val, bic_val, cpic_val
  real(kind=dp) :: min_aicc, min_bic, min_cpic

  ! Counters: number of times each candidate model is selected.
  integer, dimension(0:k_max) :: count_aicc, count_bic, count_cpic
  count_aicc = 0
  count_bic  = 0
  count_cpic = 0

  if (print_times) allocate(times(0:k_max+1))
  allocate(cp_est(k_max))

  print "(/,a30,*(i6))", "True ChangePoints (indices):", cp_true
  print "(a30,*(f6.2))", "True SD:", sigma
  print "('#obs: ', i0)", n
  if (any(cp_true < 1) .or. any(cp_true > n)) error stop "need 0 < cp_true <= n"
  do iter = 1, niter
     ! Simulate data.
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

     print "(/,a, I0)", "Iteration ", iter
     print*,"k   Log-Likelihood     AICC         BIC         CPIC        Estimated-ChangePoints"

     do k = 0, k_max
        if (print_times) call cpu_time(times(k))
        if (k == 0) then
           ! Model with no change points: one segment.
           nseg = n
           sumsq = sum(xsq)
           if (sumsq <= 0.0_dp) sumsq = 1.0e-10_dp
           ll = -0.5_dp * nseg * (1.0_dp + log(sumsq / nseg))
        else
           call find_change_points(xsq, cp_est(:k), best_ll)
           ll = best_ll
        end if

        ll_arr(k) = ll
        p = k + 1  ! p: number of free parameters (one per segment).

        aic_val  = aic(ll, p)
        aicc_val = aicc(ll, p, n)
        bic_val  = bic(ll, p, n)
        cpic_val = cpic(ll, p, n)

        aicc_arr(k) = aicc_val
        bic_arr(k)  = bic_val
        cpic_arr(k) = cpic_val

        if (k == 0) then
           write(*, "(I2,3X,F12.4,3X,F10.4,3X,F10.4,3X,F10.4)") &
                k, ll, aicc_val, bic_val, cpic_val
        else
           write(*, "(I2,3X,F12.4,3X,F10.4,3X,F10.4,3X,F10.4,3X,*(i4))") &
                k, ll, aicc_val, bic_val, cpic_val, cp_est(:k)
        end if
     end do

     ! Select best model by AICC, BIC, and CPIC.
     min_aicc = aicc_arr(0)
     min_bic  = bic_arr(0)
     min_cpic = cpic_arr(0)
     k_aicc = 0
     k_bic  = 0
     k_cpic = 0
     do k = 1, k_max
        if (aicc_arr(k) < min_aicc) then
           min_aicc = aicc_arr(k)
           k_aicc = k
        end if
        if (bic_arr(k) < min_bic) then
           min_bic = bic_arr(k)
           k_bic = k
        end if
        if (cpic_arr(k) < min_cpic) then
           min_cpic = cpic_arr(k)
           k_cpic = k
        end if
     end do

     print *
     print *, "Iteration ", iter, ": Model selected by AICC: k = ", k_aicc, " with AICC = ", min_aicc
     print *, "Iteration ", iter, ": Model selected by BIC:  k = ", k_bic,  " with BIC  = ", min_bic
     print *, "Iteration ", iter, ": Model selected by CPIC: k = ", k_cpic, " with CPIC = ", min_cpic

     count_aicc(k_aicc) = count_aicc(k_aicc) + 1
     count_bic(k_bic)   = count_bic(k_bic) + 1
     count_cpic(k_cpic) = count_cpic(k_cpic) + 1

     if (print_times) then
        call cpu_time(times(k))
        print "(/,a20, *(i8))", "#changepoints", (k, k=0, k_max)
        print "(a20,*(f8.4))", "times", times(1:) - times(:k_max)
     end if
  end do

  ! Summary table.
  print "(/,a)", "Summary of model selection over all iterations:"
  print "(*(a15))", "#changepoints", "#AICC", "#BIC", "#CPIC"
  do k = 0, k_max
     write(*, "(*(i15))") k, count_aicc(k), count_bic(k), count_cpic(k)
  end do
end program xsegmentation
