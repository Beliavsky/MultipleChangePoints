program test_split_data
  !------------------------------------------------------------------------------
  ! Program: test_split_data
  !
  ! Purpose:
  !   Tests the generalized find_multiple_change_points and compute_test_stats
  !   subroutines by performing experiments with data having multiple change points.
  !
  ! Description:
  !   - Initializes the random number generator.
  !   - Generates random data with multiple segments, each with its own mean.
  !   - Finds the best set of change points and computes test statistics.
  !   - Outputs the results, including change points, means, total SSE, and p-value.
  !------------------------------------------------------------------------------
  use split_data_mod, only: find_multiple_change_points, compute_test_stats_mcp
  use random, only: random_normal, random_seed_init
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n = 1000, nchanges = 3, nregions = nchanges + 1, &
                        ndata_sets = 5
  integer, parameter :: num_random_sets = 10**5  ! Number of random sets to try
  integer, parameter :: n_perm = 0  ! Number of permutations
  real(kind=dp), allocatable :: x(:)
  integer, allocatable :: true_change_points(:) ! , best_change_points(:)
  real(kind=dp), allocatable :: segment_means(:)
  real(kind=dp) :: total_sse, p_value
  real(kind=dp), parameter :: mean_chg = 1.0_dp
  integer :: i, iset, start_idx, end_idx
  logical, parameter :: call_random_seed_init = .true., uniform_change_points = .false.
  integer, allocatable :: best_change_points(:), segment_sizes(:)
  if (call_random_seed_init) call random_seed_init(123)
  allocate(x(n), true_change_points(nchanges), segment_means(nregions), &
     segment_sizes(nregions), best_change_points(nchanges))
  ! Generate true change points
  if (uniform_change_points) then
     do i = 1, nchanges
        true_change_points(i) = i * n / (nchanges + 1)
     end do
  else
     true_change_points = [100, 200, 700]
  end if
  ! Generate segment means
  do i = 1, nchanges + 1
    segment_means(i) = 5.0_dp + mean_chg * (i - 1)
  end do
  print*,"#obs:", n
  print*, "True Change Points:", true_change_points
  print "(a20, *(1x,f0.4))", "Segment Means:", segment_means
  do iset=1, ndata_sets
     ! Generate data with multiple change points
     start_idx = 1
     do i = 1, nchanges + 1
       if (i <= nchanges) then
         end_idx = true_change_points(i)
       else
         end_idx = n
       end if
       x(start_idx:end_idx) = random_normal(end_idx - start_idx + 1) + segment_means(i)
       start_idx = end_idx + 1
     end do
     ! Find multiple change points
     call find_multiple_change_points(x, num_random_sets, best_change_points, total_sse)
     ! Output results
     print*
     print*, "Estimated Change Points:", best_change_points
     print*, "Total SSE:", total_sse
     ! Perform permutation test
     if (n_perm > 0) then
        call compute_test_stats_mcp(x, nchanges, total_sse, p_value, n_perm)
        print*, "P-value from permutation test:", p_value
     end if
  end do
end program test_split_data
