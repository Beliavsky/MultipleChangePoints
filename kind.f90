module kind_mod
implicit none
private
public :: dp
integer, parameter :: dp = selected_real_kind(15, 307) ! double precision
end module kind_mod
