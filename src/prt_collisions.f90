module mod_collisions !VER O QUE ACONTECE SE TIVER NDT=40 E ITERMAX = 1
#if defined(_PARTICLE)
!  use mpi
  use mod_types
!  use mod_common
!  use mod_common_mpi
  use prt_mod_common       , only: ep,op,coeff_f,coeff_t
  use prt_mod_param        , only: a11_ini_pp,a11_ini_pw,a11_sat_pp,a11_sat_pw, &
                                   a11a_ini_pp,a11a_sat_pp,a11b_ini_pp,a11b_sat_pp, &
                                   a22_ini_pp,a22_ini_pw,a22_sat_pp,a22_sat_pw, &
                                   a22a_ini_pp,a22a_sat_pp,a22b_ini_pp,a22b_sat_pp, &
                                   a33_ini_pw,a33_sat_pw, &
                                   a33a_ini_pp,a33a_sat_pp,a33b_ini_pp,a33b_sat_pp, &
                                   b23_ini_pp,b23_ini_pw,b23_sat_pp,b23_sat_pw, &
                                   b23a_ini_pp,b23a_sat_pp,b23b_ini_pp,b23b_sat_pp, &
                                   b32_ini_pw,b32_sat_pw, &
                                   b32a_ini_pp,b32a_sat_pp,b32b_ini_pp,b32b_sat_pp, &
                                   c23_ini_pw,c23_sat_pw, &
                                   c23a_ini_pp,c23a_sat_pp,c23b_ini_pp,c23b_sat_pp, &
                                   c32_ini_pw,c32_sat_pw, &
                                   c32a_ini_pp,c32a_sat_pp,c32b_ini_pp,c32b_sat_pp, &
                                   d11_ini_pp,d11_ini_pw,d11_sat_pp,d11_sat_pw, &
                                   d11a_ini_pp,d11a_sat_pp,d11b_ini_pp,d11b_sat_pp, &
                                   d22_ini_pw,d22_sat_pw, &
                                   d22a_ini_pp,d22a_sat_pp,d22b_ini_pp,d22b_sat_pp, &
                                   d33_ini_pw,d33_sat_pw, &
                                   d33a_ini_pp,d33a_sat_pp,d33b_ini_pp,d33b_sat_pp, &
                                   eps_sat_pp,eps_sat_pw, &
                                   etan_ss,etan_sw,etat_ss,etat_sw, &
                                   kn_ss,kn_sw,kt_ss,kt_sw, &
                                   meffn_ss,meffn_sw,mefft_ss,mefft_sw, &
                                   muc_ss,muc_sw, &
                                   psi_crit_ss,psi_crit_sw, &
                                   np,radius
  implicit none
  private
  public collisions, lubrication
  !
  contains
  !
  subroutine collisions(rkpar,p,qq,idq,dist,deltax,deltay,deltaz, & !rkpar,p,q,qq,idp,idq,dist,deltax,deltay,deltaz, &
                        u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb,dtsub)
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer, intent(in) :: p,qq,idq !,q,idq
    real(rp), intent(in) :: dist,deltax,deltay,deltaz, &
                            u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb, &
                            dtsub
    !integer :: i
    real(rp) :: distold,dtabs
    real(rp) :: nx,ny,nz,nxold,nyold,nzold,tx,ty,tz
    !real(rp) :: dxn,dyn,dzn
    real(rp) :: dxt,dyt,dzt
    real(rp) :: dun,dvn,dwn,dut,dvt,dwt,vabs,vtabs
    real(rp) :: fxn,fyn,fzn,fnabs,fxt,fyt,fzt,ftabs
    real(rp) :: psi,psi_crit
    real(rp) :: deltan,kn,etan,kt,etat,muc
    real(rp) :: torqx,torqy,torqz
    real(rp) :: hx,hy,hz,aa,bb,cc
    real(rp) :: h11,h12,h13, &
                h21,h22,h23, &
                h31,h32,h33
    real(rp) :: meffn,mefft
    real(rp) :: rkcoeffab
    !
    rkcoeffab=rkpar(1)+rkpar(2)
    !
    ! vector (nx,ny,nz) is normal unit vector, pointing from particle p towards particle q
    !
    nx = deltax/dist
    ny = deltay/dist
    nz = deltaz/dist
    !
    if (idq <= np) then
      deltan = 2*radius-dist
      meffn = meffn_ss
      mefft = mefft_ss
      kn = kn_ss
      kt = kt_ss
      etan = etan_ss
      etat = etat_ss
      muc = muc_ss
      psi_crit = psi_crit_ss
    else
      deltan = radius-dist
      meffn = meffn_sw
      mefft = mefft_sw
      kn = kn_sw
      kt = kt_sw
      etan = etan_sw
      etat = etat_sw
      muc = muc_sw
      psi_crit = psi_crit_sw
    endif
    !
    ! relative velocity between particle p and q (Nb: p - q !)
    !
    dun = (ep(p)%u+(ep(p)%omy*radius*nz) - (ep(p)%omz*radius*ny)) - &
          (u_nb + (omy_nb*radius*(-nz)) - (omz_nb*radius*(-ny)))
    dvn = (ep(p)%v-(ep(p)%omx*radius*nz) + (ep(p)%omz*radius*nx)) - &
          (v_nb - (omx_nb*radius*(-nz)) + (omz_nb*radius*(-nx)))
    dwn = (ep(p)%w+(ep(p)%omx*radius*ny) - (ep(p)%omy*radius*nx)) - &
          (w_nb + (omx_nb*radius*(-ny)) - (omy_nb*radius*(-nx)))
    vabs = dun*nx + dvn*ny + dwn*nz
    !
    ! computation of contact forces based on soft-sphere model
    !
    fxn = - kn*(deltan*nx) - etan*(vabs*nx)
    fyn = - kn*(deltan*ny) - etan*(vabs*ny)
    fzn = - kn*(deltan*nz) - etan*(vabs*nz)
    !
    ! relative tangential velocity
    !
    dut = dun - vabs*nx
    dvt = dvn - vabs*ny
    dwt = dwn - vabs*nz
    !
    vtabs = sqrt(dut**2+dvt**2+dwt**2)
    !
    ! prediction of contact instant skipped for now for simplicity's sake
    !
    distold = sqrt(op(p)%dx(qq)**2+op(p)%dy(qq)**2+op(p)%dz(qq)**2)
    if(distold == 0) then
      nxold = nx
      nyold = ny
      nzold = nz
    else
      nxold = op(p)%dx(qq)/distold
      nyold = op(p)%dy(qq)/distold
      nzold = op(p)%dz(qq)/distold
    endif
    !
    ! computation of rotation matrix hij
    !
    if(ep(p)%firstc(qq) == idq) then
      psi = ep(p)%psi(qq)
      hx = ny*nzold-nyold*nz
      hy = nz*nxold-nzold*nx
      hz = nx*nyold-nxold*ny
      aa = sqrt(hx**2 + hy**2 + hz **2)
      cc = cos(asin(aa))
      bb = 1.0_rp - cc
      if(aa == 0.0_rp) then
        hx = 0.0_rp
        hy = 0.0_rp
        hz = 0.0_rp
      else
        hx = hx/aa
        hy = hy/aa
        hz = hz/aa
      endif
      h11 = bb*hx**2 + cc
      h12 = bb*hx*hy - aa*hz
      h13 = bb*hx*hz + aa*hy
      h21 = bb*hx*hy + aa*hz
      h22 = bb*hy**2 + cc
      h23 = bb*hy*hz - aa*hx
      h31 = bb*hx*hz - aa*hy
      h32 = bb*hy*hz + aa*hx
      h33 = bb*hz**2 + cc
      !  
      ! tangential displacement evolved in time
      !
      dxt = op(p)%dxt(qq)*h11 + op(p)%dyt(qq)*h21 + op(p)%dzt(qq)*h31 + &
            0.5_rp*dtsub*rkcoeffab*(dut + &
            op(p)%dut(qq)*h11 + op(p)%dvt(qq)*h21 + op(p)%dwt(qq)*h31)
      dyt = op(p)%dxt(qq)*h12 + op(p)%dyt(qq)*h22 + op(p)%dzt(qq)*h32 + &
            0.5_rp*dtsub*rkcoeffab*(dvt + &
            op(p)%dut(qq)*h12 + op(p)%dvt(qq)*h22 + op(p)%dwt(qq)*h32)
      dzt = op(p)%dxt(qq)*h13 + op(p)%dyt(qq)*h23 + op(p)%dzt(qq)*h33 + &
            0.5_rp*dtsub*rkcoeffab*(dwt + &
            op(p)%dut(qq)*h13 + op(p)%dvt(qq)*h23 + op(p)%dwt(qq)*h33)
    else
      !  ap(p)%firstc(qq) = idq
      psi = vtabs/abs(vabs)
      ep(p)%psi(qq) = psi
      dxt = 0.0_rp
      dyt = 0.0_rp
      dzt = 0.0_rp
    endif
    dtabs = sqrt(dxt**2+dyt**2+dzt**2)
    !
    ! computation of tangential force
    !
    !if(psi.lt.psi_crit) then ! stick
    !  fxt = - (kt*dxt) - (etat*dut)
    !  fyt = - (kt*dyt) - (etat*dvt)
    !  fzt = - (kt*dzt) - (etat*dwt)
    !  print*,'stick'
    !else ! slip
    !  fnabs = sqrt(fxn**2.+fyn**2.+fzn**2.)
    !  fxt = - muc*fnabs*tx
    !  fyt = - muc*fnabs*ty
    !  fzt = - muc*fnabs*tz
    !  print*,'slip'
    !endif
    !
    fxt = - (kt*dxt) - (etat*dut)
    fyt = - (kt*dyt) - (etat*dvt)
    fzt = - (kt*dzt) - (etat*dwt)
    fnabs = sqrt(fxn**2+fyn**2+fzn**2)
    ftabs = sqrt(fxt**2+fyt**2+fzt**2)
    tx = - fxt/ftabs
    ty = - fyt/ftabs
    tz = - fzt/ftabs
    if(ftabs <= muc*fnabs) then
      !  write(*,*) 'Sticking'
    else
      dxt = muc*fnabs*tx/kt
      dyt = muc*fnabs*ty/kt
      dzt = muc*fnabs*tz/kt
      fxt = - (kt*dxt) - (etat*dut)
      fyt = - (kt*dyt) - (etat*dvt)
      fzt = - (kt*dzt) - (etat*dwt)
      ftabs = sqrt(fxt**2+fyt**2+fzt**2)
      if(ftabs <= muc*fnabs) then
        !we were lucky
        !    write(*,*) 'frontier between stick and slip'
      else
        fxt = -muc*fnabs*tx
        fyt = -muc*fnabs*ty
        fzt = -muc*fnabs*tz
        !    write(*,*) 'Slipping'
      endif
    endif
    !
    ! computation of collision torque
    !
    torqx = radius*(ny*fzt-nz*fyt)
    torqy = radius*(nz*fxt-nx*fzt)
    torqz = radius*(nx*fyt-ny*fxt)
    !
    ! add contribution from this contact to the total contact forces/torques
    !
    ep(p)%colfx = ep(p)%colfx + fxn! + fxt*0.
    ep(p)%colfy = ep(p)%colfy + fyn! + fyt*0.
    ep(p)%colfz = ep(p)%colfz + fzn! + fzt*0.
    !
    ep(p)%coltx = ep(p)%coltx! + torqx*0.
    ep(p)%colty = ep(p)%colty! + torqy*0.
    ep(p)%coltz = ep(p)%coltz! + torqz*0.
    !
    ! update variables needed for integrating the tangential displacement
    !
    ep(p)%dx(qq) = deltax
    ep(p)%dy(qq) = deltay
    ep(p)%dz(qq) = deltaz
    ep(p)%dxt(qq) = dxt
    ep(p)%dyt(qq) = dyt
    ep(p)%dzt(qq) = dzt
    ep(p)%dut(qq) = dut
    ep(p)%dvt(qq) = dvt
    ep(p)%dwt(qq) = dwt
    !
    return
  end subroutine collisions
  !
  subroutine lubrication(p,idq,nx,ny,nz,eps,u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb) !(p,q,idp,idq,nx,ny,nz,eps,u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb)
    implicit none
    integer, intent(in) :: p,idq !,q,idp
    real(rp), intent(in) :: nx,ny,nz,eps, &
                            u_nb,v_nb,w_nb, &
                            omx_nb,omy_nb,omz_nb
    !
    ! 1-> squeezing direction; 2 and 3 -> shearing directions
    ! see paper by dance and Maxey
    !
    real(rp) :: a11,a22,a33,b23,b32,c23,c32,d11,d22,d33           
    real(rp) :: a11a,a22a,a33a,b23a,b32a,c23a,c32a,d11a,d22a,d33a 
    real(rp) :: a11b,a22b,a33b,b23b,b32b,c23b,c32b,d11b,d22b,d33b 
    real(rp) :: u1a,u2a,u3a,omg1a,omg2a,omg3a 
    real(rp) :: u1b,u2b,u3b,omg1b,omg2b,omg3b
    real(rp) :: t1x,t1y,t1z,t2x,t2y,t2z
    real(rp) :: f1,f2,f3,t1,t2,t3
    real(rp) :: dun,dvn,dwn,dut,dvt,dwt,vabs,vtabs
    !
    ! determine the unit vectors t1x and t2x from nx,ny,nz
    ! and the velocity differences
    !
    dun = (ep(p)%u+(ep(p)%omy*radius*nz) - (ep(p)%omz*radius*ny)) - &
          (u_nb + (omy_nb*radius*(-nz)) - (omz_nb*radius*(-ny)))
    dvn = (ep(p)%v-(ep(p)%omx*radius*nz) + (ep(p)%omz*radius*nx)) - &
          (v_nb - (omx_nb*radius*(-nz)) + (omz_nb*radius*(-nx)))
    dwn = (ep(p)%w+(ep(p)%omx*radius*ny) - (ep(p)%omy*radius*nx)) - &
          (w_nb + (omx_nb*radius*(-ny)) - (omy_nb*radius*(-nx)))
    vabs = dun*nx + dvn*ny + dwn*nz
    dut = dun - vabs*nx
    dvt = dvn - vabs*ny
    dwt = dwn - vabs*nz
    !
    vtabs = sqrt(dut**2+dvt**2+dwt**2)
    !
    ! t1i has the direction of the tangential velocity at contact point
    !
    t1x = dut/vtabs
    t1y = dvt/vtabs
    t1z = dwt/vtabs
    !
    ! t2i results from the cross product of ni and t1i
    !
    t2x = ny*t1z-nz*t1y
    t2y = nz*t1x-nx*t1z
    t2z = nx*t1y-ny*t1x
    !
    ! compute velocity differences in local reference frame
    !
    u1a = ep(p)%u*nx + ep(p)%v*ny + ep(p)%w*nz
    u1b = u_nb*nx + v_nb*ny + w_nb*nz
    u2a = ep(p)%u*t1x + ep(p)%v*t1y + ep(p)%w*t1z
    u2b = u_nb*t1x + v_nb*t1y + w_nb*t1z
    u3a = ep(p)%u*t2x + ep(p)%v*t2y + ep(p)%w*t2z
    u3b = u_nb*t2x + v_nb*t2y + w_nb*t2z
    omg1a = ep(p)%omx*nx + ep(p)%omy*ny+ep(p)%omz*nz
    omg1b = omx_nb*nx + omy_nb*ny+omz_nb*nz
    omg2a = ep(p)%omx*t1x + ep(p)%omy*t1y+ep(p)%omz*t1z
    omg2b = omx_nb*t1x + omy_nb*t1y+omz_nb*t1z
    omg3a = ep(p)%omx*t2x + ep(p)%omy*t2y+ep(p)%omz*t2z
    omg3b = omx_nb*t2x + omy_nb*t2y+omz_nb*t2z
    !
    ! compute lubrication forces
    !
    if(idq > np) then ! particle-wall collision
      if(eps.lt.eps_sat_pw) then
        a11 = a11_sat_pw
        a22 = a22_sat_pw
        a33 = a33_sat_pw
        b23 = b23_sat_pw
        b32 = b32_sat_pw
        c23 = c23_sat_pw
        c32 = c32_sat_pw
        d11 = d11_sat_pw
        d22 = d22_sat_pw
        d33 = d33_sat_pw
      else
        a11 = -1.0_rp/eps+1.0_rp/5.0_rp*log(eps)+1.0_rp/21.0_rp*eps*log(eps)-0.9713_rp
        a22 = 8.0_rp/15.0_rp*log(eps)+64.0_rp/375.0_rp*eps*log(eps)-0.952_rp
        a33 = a22
        b23 = -2.0_rp/15.0_rp*log(eps)-86.0_rp/375.0_rp*eps*log(eps)-0.257_rp
        b32 = -b23
        c23 = b32
        c32 = b23
        d11 = 1.0_rp/2.0_rp*eps*log(eps)-1.202_rp
        d22 = 2.0_rp/5.0_rp*log(eps)+66.0_rp/125.0_rp*eps*log(eps)-0.371_rp
        d33 = 2.0_rp/5.0_rp*log(eps)+66.0_rp/125.0_rp*eps*log(eps)-0.371_rp
      endif
      f1 = (a11*u1a) ! squeezing
      f2 = (a22*u2a + radius*b23*omg3a) ! translational + rotational shearing
      f3 = (a33*u3a + radius*b32*omg2a) ! translational + rotational shearing
      t1 = (radius*d11*omg1a) ! torsoidal shearing
      t2 = (c23*u3a+radius*d22*omg2a) ! translational + rotational shearing
      t3 = (c32*u2a+radius*d33*omg3a) ! translational + rotational shearing
    else
      if(eps < eps_sat_pp) then
        a11  = a11_sat_pp
        a11a = a11a_sat_pp
        a11b = a11b_sat_pp
        a22  = a22_sat_pp
        a22a = a22a_sat_pp
        a22b = a22b_sat_pp
        a33a = a33a_sat_pp
        a33b = a33b_sat_pp
        b23  = b23_sat_pp
        b23a = b23a_sat_pp
        b23b = b23b_sat_pp
        b32a = b32a_sat_pp
        b32b = b32b_sat_pp
        c23a = c23a_sat_pp
        c23b = c23b_sat_pp
        c32a = c32a_sat_pp
        c32b = c32b_sat_pp
        d11  = d11_sat_pp
        d11a = d11a_sat_pp
        d11a = d11b_sat_pp
        d22a = d22a_sat_pp
        d22b = d22b_sat_pp
        d33a = d33a_sat_pp
        d33b = d33b_sat_pp
      else
        a11  = -1.0_rp/4.0_rp*eps**(-1.0_rp)+9.0_rp/40.0_rp*log(eps)+3.0_rp/112.0_rp*eps*log(eps)
        a11a = a11-0.995_rp
        a11b = -a11+0.350_rp
        a22  = 1.0_rp/6.0_rp*log(eps)
        a22a = a22-0.998_rp
        a22b = -a22+0.274_rp
        a33a = a22a
        a33b = a22b
        b23  = -1.0_rp/6.0_rp*log(eps)-1.0_rp/12.0_rp*eps*log(eps)
        b23a = b23-0.159_rp
        b23b = -b23+0.001_rp
        b32a = -b23a
        b32b = -b23b
        c23a = b32a
        c23b = b32b
        c32a = b23a
        c32b = b23b
        d11  = 1.0_rp/8.0_rp*eps*log(eps)
        d11a = d11
        d11b = -d11
        d22a = 1.0_rp/5.0_rp*log(eps)+47.0_rp/250.0_rp*eps*log(eps)-0.703_rp
        d22b = -1.0_rp/20.0_rp*log(eps)+31.0_rp/500.0_rp*eps*log(eps)-0.027_rp
        d33a = 1.0_rp/5.0_rp*log(eps)+47.0_rp/250.0_rp*eps*log(eps)-0.703_rp
        d33b = -1.0_rp/20.0_rp*log(eps)+31.0_rp/500.0_rp*eps*log(eps)-0.027_rp
      endif
      f1 = (a11a*u1a+a11b*u1b) ! squeezing 
      f2 = (a22a*u2a+a22b*u2b+radius*(b23a*omg3a-b23b*omg3b)) ! translational + rotational shearing
      f3 = (a33a*u3a+a33b*u3b+radius*(b32a*omg2a-b32b*omg2b)) ! translational + rotational shearing
      t1 = (radius*(d11a*omg1a-d11b*omg1b)) ! torsoidal shearing
      t2 = (c23a*u3a+c23b*u3b+radius*(d22a*omg2a-d22b*omg2b)) ! translational + rotational shearing
      t3 = (c32a*u2a+c32b*u2b+radius*(d33a*omg3a-d33b*omg3b)) ! translational + rotational shearing
    endif
    !
    if(idq > np) then ! particle-wall collision
      a11 = a11_ini_pw
      a22 = a22_ini_pw
      a33 = a33_ini_pw
      b23 = b23_ini_pw
      b32 = b32_ini_pw
      c23 = c23_ini_pw
      c32 = c32_ini_pw
      d11 = d11_ini_pw
      d22 = d22_ini_pw
      d33 = d33_ini_pw
      f1 = coeff_f*(f1-(a11*u1a))
      f2 = coeff_f*(f2-(a22*u2a + radius*b23*omg3a))
      f3 = coeff_f*(f3-(a33*u3a + radius*b32*omg2a))
      t1 = coeff_t*(t1-(radius*d11*omg1a))
      t2 = coeff_t*(t2-(c23*u3a+radius*d22*omg2a))
      t3 = coeff_t*(t3-(c32*u2a+radius*d33*omg3a))
    else
      a11  = a11_ini_pp
      a11a = a11a_ini_pp
      a11b = a11b_ini_pp
      a22  = a22_ini_pp
      a22a = a22a_ini_pp
      a22b = a22b_ini_pp
      a33a = a33a_ini_pp
      a33b = a33b_ini_pp
      b23  = b23_ini_pp
      b23a = b23a_ini_pp
      b23b = b23b_ini_pp
      b32a = b32a_ini_pp
      b32b = b32b_ini_pp
      c23a = c23a_ini_pp
      c23b = c23b_ini_pp
      c32a = c32a_ini_pp
      c32b = c32b_ini_pp
      d11  = d11_ini_pp
      d11a  = d11a_ini_pp
      d11b  = d11b_ini_pp
      d22a = d22a_ini_pp
      d22b = d22b_ini_pp
      d33a = d33a_ini_pp
      d33b = d33b_ini_pp
      f1 = coeff_f*(f1-(a11a*u1a+a11b*u1b))
      f2 = coeff_f*(f2-(a22a*u2a+a22b*u2b+radius*(b23a*omg3a-b23b*omg3b)))
      f3 = coeff_f*(f3-(a33a*u3a+a33b*u3b+radius*(b32a*omg2a-b32b*omg2b)))
      t1 = coeff_t*(t1-(radius*(d11a*omg1a-d11b*omg1b)))
      t2 = coeff_t*(t2-(c23a*u3a+c23b*u3b+radius*(d22a*omg2a-d22b*omg2b)))
      t3 = coeff_t*(t3-(c32a*u2a+c32b*u2b+radius*(d33a*omg3a-d33b*omg3b)))
      !
    endif
    !
    ! project forces in global reference frame
    ! non-normal interactions neglected for now
    !
    ep(p)%colfx = ep(p)%colfx + f1*nx! + f2*t1x + f3*t2x
    ep(p)%colfy = ep(p)%colfy + f1*ny! + f2*t1y + f3*t2y
    ep(p)%colfz = ep(p)%colfz + f1*nz! + f2*t1z + f3*t2z
    !ap(p)%coltx = ap(p)%coltx + t1*nx + t2*t1x + t3*t2x
    !ap(p)%colty = ap(p)%coltx + t1*ny + t2*t1y + t3*t2y
    !ap(p)%coltz = ap(p)%coltz + t1*nz + t2*t1z + t3*t2z
    !
    return
  end subroutine lubrication
  !
#endif
end module mod_collisions
