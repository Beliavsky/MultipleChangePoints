program test_optimize_change_points
  use random, only: random_normal
  use split_data_mod, only: optimize_change_points
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 150, ncp = 2, &
                        true_cps(ncp) = [61, 101], nsets = 5
  real(dp), parameter :: mean_chg = 2.0_dp
  real(dp) :: xtrue(n), x(n)
  integer :: i, initial_cps(ncp), optimized_cps(ncp), iset
  call random_seed()
  print *, "     True Change Points:", true_cps
  do i=1, n
     if (i < true_cps(1)) then
        xtrue(i) = 0.0_dp
     else if (i < true_cps(2)) then
        xtrue(i) = mean_chg
     else
        xtrue(i) = 2*mean_chg
     end if
  end do
  initial_cps = [40, 110]
  do iset=1, nsets
     x = xtrue + random_normal(n)
     call optimize_change_points(x, initial_cps, optimized_cps)
     ! Output the initial and optimized change points
     print*
     print*, "  Initial Change Points:", initial_cps
     print*, "Optimized Change Points:", optimized_cps
  end do
end program test_optimize_change_points
