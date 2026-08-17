module prt_mod_digitiser
#if defined(_PARTICLE)
!#if defined(_EULER)
  use mod_types
  use mod_param         , only: dl
  use prt_mod_param     , only: eps_sol
  !
  implicit none
  private
  public digitiser
  !
  contains
  !
  subroutine digitiser(ds, n, alpha)
    !
    implicit none
    real(rp), intent(in) :: ds
    real(rp), dimension(3), intent(in) :: n
    real(rp) :: delta, sigma, lambda
    real(rp), intent(out) :: alpha
    !
    delta=((dl(1))*(dl(2))*(dl(3)))**(1.0_rp/3)
    lambda=abs(n(1))+abs(n(2))+abs(n(3))
    sigma=0.05_rp*(1._rp-lambda**2)+0.3_rp
    if (ds<=eps_sol*delta) then
      alpha = 0.5*(1-tanh(ds/(sigma*lambda*delta)))
    else
      alpha = 0
    endif
    !
    return
  end subroutine digitiser
  !
!#endif
#endif
end module prt_mod_digitiser
