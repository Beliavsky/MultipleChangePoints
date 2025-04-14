module kind_mod
implicit none
private
public :: sp, dp
integer, parameter :: sp = kind(1.0), dp = selected_real_kind(15, 307) ! double precision
end module kind_mod
