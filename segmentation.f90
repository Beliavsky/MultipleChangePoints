module segmentation_mod
  use kind_mod, only: dp
  implicit none
  private
  public :: find_change_points
contains

  pure subroutine find_change_points(xsq, best_cp, best_ll)
    !> Find change points in a zero-mean time series assuming that the
    !> standard deviation changes at k discrete points.
    !>
    !> xsq is the data array squared.
    !> The routine precomputes the cumulative sum of xsq (cumulsq) and passes
    !> it to the search and evaluation subroutines.
    real(kind=dp), intent(in)  :: xsq(:)     ! data array squared
    integer      , intent(out) :: best_cp(:) ! best change points.
    real(kind=dp), intent(out) :: best_ll    ! maximum log-likelihood value.
    integer                    :: n, k, cp(size(best_cp)), i
    real(kind=dp), allocatable :: cumulsq(:)

    n = size(xsq)
    k = size(best_cp)

    ! Initialize best_ll to a very low value and clear candidate change points.
    best_ll = -huge(1.0_dp)
    best_cp = 0
    cp = 0

    ! Compute the cumulative sum of xsq:
    allocate(cumulsq(n))
    cumulsq(1) = xsq(1)
    do i = 2, n
       cumulsq(i) = cumulsq(i-1) + xsq(i)
    end do

    call search_cp(cumulsq, 1, cp, best_ll, best_cp)
  end subroutine find_change_points

  pure recursive subroutine search_cp(cumulsq, start, cp, best_ll, best_cp)
    !> Recursively loop over all combinations of k change points using the
    !> precomputed cumulative sum array.
    !
    !> Input:
    !>   cumulsq  - cumulative sum of the squared data values.
    !>   start    - starting index for the current recursion.
    !>   cp       - candidate change points (indices marking the end of segments).
    !> In/Out:
    !>   best_ll  - best (highest) log-likelihood found so far.
    !>   best_cp  - best candidate change points corresponding to best_ll.
    integer, intent(in)     :: start       ! starting index for current recursion
    real(kind=dp), intent(in) :: cumulsq(:)  ! cumulative sum array (size n)
    integer, intent(in out) :: cp(:)         ! candidate change points array (size k)
    real(kind=dp), intent(in out) :: best_ll   ! best log-likelihood so far
    integer, intent(in out) :: best_cp(:)      ! best candidate change points (size k)
    integer :: i, k, n, pos

    n = size(cumulsq)
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

    ! If all change points have been set, evaluate this candidate segmentation.
    if (pos == 0) then
      call evaluate_segmentation(cumulsq, cp, best_ll, best_cp)
      return
    end if

    ! The upper bound for a valid candidate index is n - (k - pos)
    do i = start, n - (k - pos)
      cp(pos) = i
      call search_cp(cumulsq, i+1, cp, best_ll, best_cp)
      cp(pos) = 0  ! reset for the next candidate
    end do
  end subroutine search_cp

  pure subroutine evaluate_segmentation(cumulsq, cp, best_ll, best_cp)
    !> Evaluate the log-likelihood of a candidate segmentation using the cumulative
    !> sum array. The data are segmented into k+1 parts. For segment i the log-likelihood is:
    !>
    !>   L_seg = -0.5 * n_seg * [ 1 + log( sum(x**2)/n_seg ) ]
    !>
    !> The cumulative sum array (cumulsq) is assumed to satisfy:
    !>
    !>   cumulsq(i) = xsq(1) + xsq(2) + ... + xsq(i)
    !>
    !> To compute the segment sum over indices seg_start to seg_end:
    !>
    !>   if (seg_start == 1) then
    !>       sumsq = cumulsq(seg_end)
    !>   else
    !>       sumsq = cumulsq(seg_end) - cumulsq(seg_start - 1)
    !>
    !> Input/Output:
    !>   cumulsq  - cumulative sum of the squared data.
    !>   cp       - candidate change points (indices marking the end of segments).
    !>   best_ll  - best log-likelihood so far.
    !>   best_cp  - candidate change points corresponding to best_ll.
    real(kind=dp), intent(in) :: cumulsq(:)  ! cumulative sum array
    integer, intent(in)       :: cp(:)       ! candidate change points indices (size k)
    real(kind=dp), intent(in out) :: best_ll   ! best log-likelihood so far
    integer, intent(in out)   :: best_cp(:)    ! best candidate change points (size k)
    real(kind=dp) :: ll, ll_seg, sumsq
    integer :: seg_start, seg_end, nseg, i, n, k

    n = size(cumulsq)
    k = size(cp)
    if (size(best_cp) /= k) error stop "need size(cp) == size(best_cp)"

    ll = 0.0_dp
    seg_start = 1
    ! Loop over the first k segments.
    do i = 1, k
       seg_end = cp(i)
       nseg = seg_end - seg_start + 1
       if (nseg <= 0) return  ! invalid segmentation; indices should be ascending
       if (seg_start == 1) then
          sumsq = cumulsq(seg_end)
       else
          sumsq = cumulsq(seg_end) - cumulsq(seg_start - 1)
       end if
       ! Protect against taking log of zero.
       if (sumsq <= 0.0_dp) sumsq = 1.0e-10_dp
       ll_seg = -0.5_dp * nseg * (1.0_dp + log(sumsq/nseg))
       ll = ll + ll_seg
       seg_start = seg_end + 1
    end do

    ! Process the last segment from cp(k)+1 to n.
    seg_end = n
    nseg = seg_end - cp(k)
    if (nseg <= 0) return  ! invalid segmentation
    sumsq = cumulsq(seg_end) - cumulsq(cp(k))
    if (sumsq <= 0.0_dp) sumsq = 1.0e-10_dp
    ll_seg = -0.5_dp * nseg * (1.0_dp + log(sumsq/nseg))
    ll = ll + ll_seg

    ! Update the best candidate segmentation if this candidate's log-likelihood is higher.
    if (ll > best_ll) then
       best_ll = ll
       best_cp = cp
    end if
  end subroutine evaluate_segmentation

end module segmentation_mod
