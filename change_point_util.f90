module change_point_util
  implicit none
  private
  public :: change_points_neighbors

contains

  !------------------------------------------------------------------------------
  ! Function: generate_neighbor_change_points
  !
  ! Purpose:
  !   Generates a matrix of neighboring change point sets by moving each change point
  !   one position to the left or right, ensuring that change points remain within
  !   valid bounds and maintain strict ascending order.
  !
  ! Arguments:
  !   current_cps - Array of current change point locations (INTENT(IN))
  !   N           - Total number of observations (INTENT(IN))
  !
  ! Returns:
  !   neighbor_cps - Allocatable matrix where each row is a new set of change points
  !                  with one change point moved left or right (INTENT(OUT))
  !
  ! Example:
  !   current_cps = [4, 7, 10], N = 10
  !   neighbor_cps =
  !     [3, 7, 10]
  !     [5, 7, 10]
  !     [4, 6, 10]
  !     [4, 8, 10]
  !     [4, 7, 9]
  !------------------------------------------------------------------------------
  function change_points_neighbors(current_cps, N) result(neighbor_cps)
    integer, intent(in) :: current_cps(:)
    integer, intent(in) :: N
    integer, allocatable :: neighbor_cps(:,:)
    integer :: K
    integer, allocatable :: temp_cps(:,:)
    integer :: i
    integer :: new_cp
    logical :: can_move_left, can_move_right
    integer :: prev_cp, next_cp
    integer, allocatable :: new_set(:)
    integer :: max_neighbors, num_neighbors

    ! Determine the number of change points
    K = size(current_cps)

    ! Maximum possible neighbors is 2*K (each change point can move left and right)
    max_neighbors = 2*K
    allocate(temp_cps(max_neighbors, K))
    num_neighbors = 0

    do i = 1, K
      ! Current change point
      new_cp = current_cps(i)

      ! Attempt to move the i-th change point to the left
      can_move_left = .true.
      if (new_cp > 2) then
        if (i > 1) then
          prev_cp = current_cps(i-1)
          if ((new_cp - 1) <= prev_cp) then
            can_move_left = .false.
          end if
        end if
      else
        can_move_left = .false.
      end if

      if (can_move_left) then
        num_neighbors = num_neighbors + 1
        allocate(new_set(K))
        new_set = current_cps
        new_set(i) = new_cp - 1
        temp_cps(num_neighbors, :) = new_set
        deallocate(new_set)
      end if

      ! Attempt to move the i-th change point to the right
      can_move_right = .true.
      if (new_cp < N) then
        if (i < K) then
          next_cp = current_cps(i+1)
          if ((new_cp + 1) >= next_cp) then
            can_move_right = .false.
          end if
        end if
      else
        can_move_right = .false.
      end if

      if (can_move_right) then
        num_neighbors = num_neighbors + 1
        allocate(new_set(K))
        new_set = current_cps
        new_set(i) = new_cp + 1
        temp_cps(num_neighbors, :) = new_set
        deallocate(new_set)
      end if
    end do
    ! Allocate the neighbor_cps matrix with the exact number of neighbors found
    if (num_neighbors > 0) then
      allocate(neighbor_cps(num_neighbors, K))
      neighbor_cps = temp_cps(1:num_neighbors, :)
    else
      ! No valid neighbors found
      allocate(neighbor_cps(0, K))
    end if
  end function change_points_neighbors
end module change_point_util
