module segmentation_mod
  use kind_mod, only: dp
  implicit none
  private
  public :: find_change_points
contains

  pure subroutine find_change_points(x, best_cp, best_ll)
    !> Find change points in a zero-mean time series assuming that the
    !> standard deviation changes at k discrete points.
    real(kind=dp), intent(in)  :: x(:)       ! data array
    integer      , intent(out) :: best_cp(:) ! best change points.
    real(kind=dp), intent(out) :: best_ll    ! maximum log-likelihood value.
    integer                    :: n, k, cp(size(best_cp))
    n = size(x)       ! # of observations
    k = size(best_cp) ! # of change points
    ! Initialize best_ll to a very low value and clear candidate change points.
    best_ll = -huge(1.0_dp)
    best_cp = 0
    cp = 0
    call search_cp(x, 1, cp, best_ll, best_cp)
  end subroutine find_change_points

  pure recursive subroutine search_cp(x, start, cp, best_ll, best_cp)
    !> Recursively loop over all combinations of k change points.
    integer      , intent(in)     :: start       ! starting index for the current recursion
    real(kind=dp), intent(in)     :: x(:)        ! data array
    integer      , intent(in out) :: cp(:)       ! current candidate array for change points
    real(kind=dp), intent(in out) :: best_ll     ! maximum log-likelihood found so far
    integer      , intent(in out) :: best_cp(:)  ! best candidate change points.
    integer :: i, k, n, pos
    n = size(x)
    k = size(cp)
    if (size(best_cp) /= k) error stop "need size(cp) == size(best_cp)"
    ! Determine the next position (1 <= pos <= k) to set in cp.
    pos = 0
    do i = 1, k
       if (cp(i) == 0) then
          pos = i
          exit
       end if
    end do
    ! If pos remains zero, all change points have been set;
    ! evaluate this candidate segmentation.
    if (pos == 0) then
      call evaluate_segmentation(x, cp, best_ll, best_cp)
      return
    end if
    ! The upper bound for a valid candidate index is n - (k - pos)
    do i = start, n - (k - pos)
      cp(pos) = i
      call search_cp(x, i+1, cp, best_ll, best_cp)
      cp(pos) = 0  ! reset for the next candidate
    end do
  end subroutine search_cp

  pure subroutine evaluate_segmentation(x, cp, best_ll, best_cp)
    !> Evaluate the log-likelihood of a candidate segmentation.
    !>
    !> The data are segmented into k+1 parts. For segment i the log-likelihood is:
    !>   L_seg = -0.5 * n_seg * [ 1 + log(sum(x**2)/n_seg) ]
    real(kind=dp), intent(in)     :: x(:)       ! data array
    integer      , intent(in)     :: cp(:)      ! candidate change points (indices marking the end of segments)
    real(kind=dp), intent(in out) :: best_ll    ! best log-likelihood so far
    integer      , intent(in out) :: best_cp(:) ! best candidate change points
    real(kind=dp) :: ll, ll_seg, sumsq
    integer :: seg_start, seg_end, nseg, i, n, k
    n = size(x)
    k = size(cp)
    if (size(best_cp) /= k) error stop "need size(cp) == size(best_cp)"
    ll = 0.0_dp
    seg_start = 1
    ! Loop over the first k segments.
    do i = 1, k
       seg_end = cp(i)
       nseg = seg_end - seg_start + 1
       if (nseg <= 0) return  ! invalid segmentation (should not happen if indices are ascending)
       sumsq = sum(x(seg_start:seg_end)**2)
       if (sumsq <= 0.0_dp) sumsq = 1.0e-10_dp  ! avoid log(0)
       ll_seg = -0.5_dp * nseg * (1.0_dp + log(sumsq/nseg))
       ll = ll + ll_seg
       seg_start = seg_end + 1
    end do
    ! Last segment: from cp(k)+1 to n.
    seg_end = n
    nseg = seg_end - cp(k)
    if (nseg <= 0) return  ! invalid segmentation
    sumsq = sum(x(cp(k)+1:n)**2)
    if (sumsq <= 0.0_dp) then
       sumsq = 1.0e-10_dp
    end if
    ll_seg = -0.5_dp * nseg * (1.0_dp + log(sumsq/nseg))
    ll = ll + ll_seg
    ! Save this candidate if its log-likelihood is higher than any found before.
    if (ll > best_ll) then
       best_ll = ll
       best_cp = cp
    end if
  end subroutine evaluate_segmentation

end module segmentation_mod
