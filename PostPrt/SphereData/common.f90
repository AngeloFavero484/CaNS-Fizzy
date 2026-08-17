module mod_common
  !
  use mod_param       , only: rp,nl
  !
  implicit none
  public
  !
  real(rp), dimension(nl) :: lpx,lpy,lpz
  real(rp), dimension(nl) :: nx,ny,nz
  real(rp), dimension(nl) :: lptheta,lpphi
  real(rp) :: rad
  real(rp) :: lpA
  !
end module mod_common
