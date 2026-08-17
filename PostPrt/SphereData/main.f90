program sphere
  !
  use mod_param     , only: rp,nl,ng,radius,retraction,dl
  use mod_loadd     , only: flpcoord,loadfld 
  use mod_common    , only: lpx,lpy,lpz,rad,lptheta,lpphi
  use mod_forces    , only: CmptDNSForces
  !
  implicit none
  !
  real(rp), allocatable :: u(:,:,:),v(:,:,:),w(:,:,:),p(:,:,:)
  real(rp) :: fx,fy,fz
  !
  allocate(u(1:ng(1),1:ng(2),1:ng(3)), &
           v(1:ng(1),1:ng(2),1:ng(3)), &
           w(1:ng(1),1:ng(2),1:ng(3)), &
           p(1:ng(1),1:ng(2),1:ng(3)))
  !
  call loadfld(u,v,w,p)
  !
  rad=radius
  call flpcoord(lptheta,lpphi,rad)
  !
  call CmptDNSForces(ng,dl,u,v,w,p,fx,fy,fz)
  !
end program sphere
