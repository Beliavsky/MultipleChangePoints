module split_data_mod
  use kind_mod, only: dp
  use random, only: random_shuffle
use change_point_util, only: change_points_neighbors
  implicit none
  private
  public :: find_multiple_change_points, compute_test_stats_mcp, optimize_change_points
contains
  !------------------------------------------------------------------------------
  ! Subroutine: total_sse
  !
  ! Purpose:
  !   Computes the total SSE given the data and a set of change points.
  !
  ! Arguments:
  !   x          - Input data array (INTENT(IN))
  !   change_cps - Array of change point locations (INTENT(IN))
  !   sse_total  - Total SSE for the given partition (INTENT(OUT))
  !------------------------------------------------------------------------------
  subroutine total_sse(x, change_cps, sse_total)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: change_cps(:)
    real(dp), intent(out) :: sse_total
    integer :: i, K, n, start_cp, end_cp
    real(dp) :: segment_mean, segment_sse

    n = size(x)
    K = size(change_cps)
    sse_total = 0.0_dp

    ! First segment: from 1 to the first change point
    start_cp = 1
    end_cp = change_cps(1)
    call mean_and_sse(x(start_cp:end_cp), segment_mean, segment_sse)
    sse_total = sse_total + segment_sse

    ! Intermediate segments
    do i = 2, K
      start_cp = change_cps(i - 1) + 1
      end_cp = change_cps(i)
      call mean_and_sse(x(start_cp:end_cp), segment_mean, segment_sse)
      sse_total = sse_total + segment_sse
    end do

    ! Last segment: from last change point to the end
    start_cp = change_cps(K) + 1
    if (start_cp <= n) then
      call mean_and_sse(x(start_cp:n), segment_mean, segment_sse)
      sse_total = sse_total + segment_sse
    end if
  end subroutine total_sse
  !------------------------------------------------------------------------------
  ! Subroutine: mean_and_sse
  !
  ! Purpose:
  !   Computes the mean and sum of squared errors (SSE) for the input array.
  !
  ! Arguments:
  !   x    - Input data array (INTENT(IN))
  !   mean - Calculated mean of the array (INTENT(OUT))
  !   sse  - Sum of squared errors of the array (INTENT(OUT))
  !------------------------------------------------------------------------------
  subroutine mean_and_sse(x, mean, sse)
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp), intent(out) :: mean, sse
    integer :: n
    n = size(x)
    if (n > 0) then
      mean = sum(x) / n
      sse = sum((x - mean)**2)
    else
      mean = 0.0_dp
      sse = 0.0_dp
    end if
  end subroutine mean_and_sse
  !------------------------------------------------------------------------------
  ! Subroutine: find_multiple_change_points
  !
  ! Purpose:
  !   Finds the best set of change points by trying a specified number of random sets
  !   of ascending integers and selecting the set with the minimum total SSE.
  !
  ! Arguments:
  !   x                  - Input data array (INTENT(IN))
  !   num_change_points  - Number of change points to find (INTENT(IN))
  !   num_random_sets    - Number of random sets to try (INTENT(IN))
  !   best_change_points - Best set of change points found (INTENT(OUT))
  !   total_sse          - Total sum of squared errors for the best set (INTENT(OUT))
  !------------------------------------------------------------------------------
  subroutine find_multiple_change_points(x, num_random_sets, best_change_points, total_sse)
    real(kind=dp), intent(in)  :: x(:)
    integer, intent(in)        :: num_random_sets
    integer, intent(out) :: best_change_points(:)
    real(kind=dp), intent(out) :: total_sse
    integer :: n, set_idx, cp_idx, num_change_points
    integer, allocatable :: candidate_change_points(:)
    real(kind=dp) :: sse, min_total_sse, segment_mean, segment_sse
    real(kind=dp), allocatable :: segment(:)
    num_change_points = size(best_change_points)
    n = size(x)
    min_total_sse = huge(0.0_dp)
    do set_idx = 1, num_random_sets
      ! Generate a set of random change points
      call generate_random_change_points(n, num_change_points, candidate_change_points)
      sse = 0.0_dp
      ! Compute SSE for each segment
      do cp_idx = 1, num_change_points + 1
        if (cp_idx == 1) then
          segment = x(1 : candidate_change_points(1))
        else if (cp_idx == num_change_points + 1) then
          segment = x(candidate_change_points(num_change_points) + 1 : n)
        else
          segment = x(candidate_change_points(cp_idx - 1) + 1 : candidate_change_points(cp_idx))
        end if
        call mean_and_sse(segment, segment_mean, segment_sse)
        sse = sse + segment_sse
      end do
      ! Check if this set has the minimum total SSE
      if (sse < min_total_sse) then
        min_total_sse = sse
        best_change_points = candidate_change_points
      end if
      deallocate(candidate_change_points)
    end do
    total_sse = min_total_sse
  end subroutine find_multiple_change_points

  !------------------------------------------------------------------------------
  ! Subroutine: generate_random_change_points
  !
  ! Purpose:
  !   Generates a sorted array of ascending unique random integers as change points.
  !
  ! Arguments:
  !   n                 - Size of the data array (INTENT(IN))
  !   num_change_points - Number of change points to generate (INTENT(IN))
  !   change_points     - Sorted array of change points (INTENT(OUT))
  !------------------------------------------------------------------------------
  subroutine generate_random_change_points(n, num_change_points, change_points)
    integer, intent(in)  :: n, num_change_points
    integer, intent(out), allocatable :: change_points(:)
    integer, allocatable :: temp_points(:)
    integer :: i, total_points
    total_points = n - 2
    allocate(temp_points(total_points))
    ! Generate positions from 2 to n - 1
    temp_points = [(i, i = 2, n - 1)]
    ! Shuffle temp_points
    call random_shuffle_int(temp_points)
    ! Select the first num_change_points positions as change points
    change_points = temp_points(1:num_change_points)
    ! Sort the change points in ascending order
    call sort_integers(change_points)
    deallocate(temp_points)
  end subroutine generate_random_change_points

  !------------------------------------------------------------------------------
  ! Subroutine: random_shuffle_int
  !
  ! Purpose:
  !   Randomly shuffles an integer array using the Fisher-Yates algorithm.
  !
  ! Arguments:
  !   x - The integer array to shuffle (INTENT(INOUT))
  !------------------------------------------------------------------------------
  subroutine random_shuffle_int(x)
    integer, intent(inout) :: x(:)
    integer :: i, j, temp
    real(kind=dp) :: rand_val

    do i = size(x), 2, -1
      call random_number(rand_val)
      j = int(rand_val * i) + 1
      if (j /= i) then
        temp = x(j)
        x(j) = x(i)
        x(i) = temp
      end if
    end do
  end subroutine random_shuffle_int

  !------------------------------------------------------------------------------
  ! Subroutine: sort_integers
  !
  ! Purpose:
  !   Sorts an integer array in ascending order using the insertion sort algorithm.
  !
  ! Arguments:
  !   array - The integer array to sort (INTENT(INOUT))
  !------------------------------------------------------------------------------
subroutine sort_integers(array)
  integer, intent(inout) :: array(:)
  integer :: i, j, key, n
  n = size(array)
  do i = 2, n
    key = array(i)
    j = i - 1
    do
      if (j > 0) then
        if (array(j) > key) then
          array(j + 1) = array(j)
          j = j - 1
        else
          exit
        end if
      else
        exit
      end if
    end do
    array(j + 1) = key
  end do
end subroutine sort_integers

  !------------------------------------------------------------------------------
  ! Subroutine: compute_test_stats_mcp
  !
  ! Purpose:
  !   Calculates the total SSE for the observed data and performs a permutation test
  !   to determine the p-value, accounting for multiple testing.
  !
  ! Arguments:
  !   x                  - Original data array (INTENT(IN))
  !   num_change_points  - Number of change points (INTENT(IN))
  !   total_sse          - Total SSE for the best set (INTENT(IN))
  !   p_value            - P-value from permutation test (INTENT(OUT))
  !   n_perm             - Number of permutations to perform (INTENT(IN))
  !------------------------------------------------------------------------------
  subroutine compute_test_stats_mcp(x, num_change_points, total_sse, p_value, n_perm)
    real(kind=dp), intent(in)  :: x(:)
    integer, intent(in)        :: num_change_points
    real(kind=dp), intent(in)  :: total_sse
    integer, intent(in)        :: n_perm
    real(kind=dp), intent(out) :: p_value
    real(kind=dp), allocatable :: permuted_x(:)
    integer, allocatable :: perm_change_points(:)
    real(kind=dp) :: perm_sse
    integer :: i, count, n
    n = size(x)
    count = 0
    allocate(permuted_x(n), perm_change_points(num_change_points))
    do i = 1, n_perm
      permuted_x = x
      call random_shuffle(permuted_x)
      call find_multiple_change_points(permuted_x, 1, perm_change_points, perm_sse)
      if (perm_sse <= total_sse) then
        count = count + 1
      end if
    end do
    p_value = real(count + 1, kind=dp) / real(n_perm + 1, kind=dp)
  end subroutine compute_test_stats_mcp
!
  !------------------------------------------------------------------------------
  ! Subroutine: optimize_change_points
  !
  ! Purpose:
  !   Iteratively finds the best neighbor in terms of total SSE, starting from an
  !   initial guess, and stops when no neighboring set of change points provides
  !   a lower SSE.
  !
  ! Arguments:
  !   x               - Input data array (INTENT(IN))
  !   initial_cps      - Initial guess of change points (INTENT(IN))
  !   optimized_cps    - Optimized change points with lowest SSE (INTENT(OUT))
  !------------------------------------------------------------------------------
  subroutine optimize_change_points(x, initial_cps, optimized_cps)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: initial_cps(:)
    integer, intent(out) :: optimized_cps(:)
    integer, allocatable :: neighbor_cps(:,:)
    integer :: num_neighbors, i
    real(dp) :: current_sse, neighbor_sse
    integer :: best_neighbor_idx
    logical :: improvement
    ! Initialize with the initial guess of change points
    optimized_cps = initial_cps
    call total_sse(x, optimized_cps, current_sse)
    improvement = .true.
    ! Iterate until no better neighbor is found
    do while (improvement)
      improvement = .false.
      ! Generate neighboring sets of change points
      neighbor_cps = change_points_neighbors(optimized_cps, size(x))
      num_neighbors = size(neighbor_cps, 1)
      ! Check all neighbors to find the one with the lowest SSE
      do i = 1, num_neighbors
        call total_sse(x, neighbor_cps(i,:), neighbor_sse)
        if (neighbor_sse < current_sse) then
          current_sse = neighbor_sse
          best_neighbor_idx = i
          improvement = .true.
        end if
      end do
      ! Update the optimized change points if a better neighbor was found
      if (improvement) then
        optimized_cps = neighbor_cps(best_neighbor_idx, :)
      end if
      if (allocated(neighbor_cps)) deallocate(neighbor_cps)
    end do
  end subroutine optimize_change_points
end module split_data_mod
