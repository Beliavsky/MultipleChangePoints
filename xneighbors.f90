program xneighbors
  use change_point_util, only: change_points_neighbors
  implicit none
  integer, parameter :: K = 3     ! Number of change points
  integer, parameter :: N = 10    ! Total number of observations
  integer, dimension(K) :: current_cps
  integer, allocatable :: neighbor_cps(:,:)
  integer :: i, num_neighbors
  ! Initialize current change points
  current_cps = [4, 7, 10]
  ! Call the function to generate neighboring change points
  neighbor_cps = change_points_neighbors(current_cps, N)
  ! Determine the number of neighbors generated
  num_neighbors = size(neighbor_cps, 1)
  ! Print the current change points
  print *, "Current Change Points:", current_cps
  ! Print the neighboring change points
  if (num_neighbors > 0) then
    print *, "Neighboring Change Points:"
    do i = 1, num_neighbors
      write(*, '(3I5)') neighbor_cps(i, :)
    end do
  else
    print *, "No valid neighboring change points found."
  end if
end program xneighbors