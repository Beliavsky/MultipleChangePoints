module info_crit_mod
  use kind_mod, only: dp  ! Using double-precision kind dp.
  implicit none
  private
  public :: aic, aicc, bic, cpic
  contains

  elemental pure function aic(ll, p) result(aic_val)
    ! Compute AIC = -2*ll + 2*p
    real(kind=dp), intent(in) :: ll  ! ll: log-likelihood of the model.
    integer, intent(in)       :: p   ! p: number of parameters.
    real(kind=dp)             :: aic_val

    aic_val = -2.0_dp * ll + 2.0_dp * real(p, dp)
  end function aic

  elemental pure function aicc(ll, p, n) result(aicc_val)
    ! Compute AICC = AIC + [2*p*(p+1)]/(n - p - 1)
    real(kind=dp), intent(in) :: ll  ! ll: log-likelihood of the model.
    integer, intent(in)       :: p   ! p: number of parameters.
    integer, intent(in)       :: n   ! n: sample size.
    real(kind=dp)             :: aicc_val, aic_val

    aic_val = aic(ll, p)
    aicc_val = aic_val + 2.0_dp * real(p, dp) * (real(p, dp) + 1.0_dp) / ( real(n, dp) - real(p, dp) - 1.0_dp )
  end function aicc

  elemental pure function bic(ll, p, n) result(bic_val)
    ! Compute BIC = -2*ll + p*log(n)
    real(kind=dp), intent(in) :: ll  ! ll: log-likelihood of the model.
    integer, intent(in)       :: p   ! p: number of parameters.
    integer, intent(in)       :: n   ! n: sample size.
    real(kind=dp)             :: bic_val

    bic_val = -2.0_dp * ll + real(p, dp) * log(real(n, dp))
  end function bic

  elemental pure function cpic(ll, p, n) result(cpic_val)
    ! Compute the alternative change-point information criterion:
    ! CPIC = -2*ll + p*log(n) + 2*(p-1)*log(log(n))
    ! ll: log-likelihood of the model.
    ! p: number of parameters (p = k+1, so p-1 is the number of change points).
    ! n: sample size.
    real(kind=dp), intent(in) :: ll  ! ll: log-likelihood.
    integer, intent(in)       :: p   ! p: number of parameters.
    integer, intent(in)       :: n   ! n: sample size.
    real(kind=dp)             :: cpic_val

    cpic_val = -2.0_dp * ll + real(p, dp) * log(real(n, dp)) + 2.0_dp * real(p-1, dp) * log(log(real(n, dp)))
  end function cpic

end module info_crit_mod
