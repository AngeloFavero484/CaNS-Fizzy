module prt_mod_eulint
#if defined(_PARTICLE)
#if defined(_EULER)
  use mod_types
  use mpi
  use mod_common_mpi    , only: boundfrontmyid,boundleftmyid,prt_comm_cart,ierr,myid
  use mod_param         , only: l,dli,dims,nh_wide,rho12,sigma
  use prt_mod_param     , only: radius,eps_sol
  use prt_mod_common    , only: ep,npmax,pmax,dVeul,neighbor,coords,retrac, &
                                norm_partx,norm_party,norm_partz, &
                                fx_tot,fy_tot,fz_tot, &
                                tx_tot,ty_tot,tz_tot
  use prt_mod_digitiser , only: digitiser
  !
  implicit none
  !
  private
  public :: eulint
  !
  contains
  !
  subroutine eulint(rkpar,dttot,rkit,n,psi,psio,kappa,fx_old,fy_old,fz_old,unew,vnew,wnew)
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    real(rp), intent(in) :: dttot
    integer, intent(in) :: rkit
    integer, dimension(3), intent(in) :: n
    real(rp), dimension(1-nh_wide:,1-nh_wide:,1-nh_wide:), intent(in) :: psi
    real(rp), dimension(1-nh_wide:,1-nh_wide:,1-nh_wide:), intent(in) :: psio
    real(rp), dimension(1-nh_wide:,1-nh_wide:,1-nh_wide:), intent(in) :: kappa
    real(rp), dimension(1-nh_wide:,1-nh_wide:,1-nh_wide:), intent(inout) :: fx_old,fy_old,fz_old,unew,vnew,wnew
    type pneighbor
       real(rp) :: x,y,z, &
                   intu,intv,intw, &
                   intomx, intomy, intomz
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb ! can be pmax because it is a subroutine-specific array!
    real(rp), dimension(1:n(1),1:n(2),1:n(3)) :: ustar,vstar,wstar
    integer :: i,j,k,p,s
    integer :: ilow,ihigh,jlow,jhigh,klow,khigh
    real(rp) :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    real(rp) :: coorxc,cooryc,coorzc
    real(rp) :: coorx_cent,coorx_stag,coory_cent,coory_stag,coorz_cent,coorz_stag
    real(rp) :: radin,radin2,dist2x,dist2y,dist2z,dist2c
    real(rp) :: dx,dy,dz
    real(rp) :: sum1,sum2
    integer :: nb,nbsend,nbrecv
    integer :: nrrequests
    integer :: arrayrequests(1:3)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    real(rp) :: sgndist(1:8)
    real(rp) :: rx,ry,rz
    integer :: idp,tag
!    real(rp) :: dxm2,dxp2,dym2,dyp2,dzm2,dzp2
    real(rp) :: dxp2_cent,dxp2_stag,dyp2_cent,dyp2_stag,dzp2_cent,dzp2_stag
    real(rp) :: fx,fy,fz
    real(rp) :: fsurf_x,fsurf_y,fsurf_z
    real(rp) :: ul,vl,wl
    real(rp) :: auxu,auxv,auxw,auxomx,auxomy,auxomz
!    real(rp) :: intu,intv,intw
    real(rp) :: alpha_eulx, alpha_euly, alpha_eulz,alpha_eulc
    real(rp) :: deltasx,deltasy,deltasz,deltasc
    real(rp) :: normalx(1:3), normaly(1:3), normalz(1:3), normalc(1:3)
    real(rp) :: dt,dti,dtitot
    real(rp) :: rho,drho,rhox,rhoy,rhoz,rhox_old,rhoy_old,rhoz_old
!    real(rp) :: buf_send(3),buf_recv(3)
!    real(rp), dimension(1:3,0:8,1:npmax) :: buf_send,buf_recv
    !
    dx       = 1./dli(1)
    dy       = 1./dli(2)
    dz       = 1./dli(3)
    rho=rho12(2)
    drho=rho12(1)-rho12(2)
    radin    = (radius-retrac)-eps_sol*dx
    radin2   = radin**2
    dt=(rkpar(1)+rkpar(2))*dttot
    dti=1./dt
    dtitot=1./dttot
!    alphac(:,:,:)=0._rp
    !
    ! for r < radius - sqrt( (dx)**2 + (dy)**2 + (dz)**2 ) = 
    !        radius - 1.732050808*dx : all cell corner points within sphere
    !
    ! first step: slaves need from master x,..,xfp,.. etc
    !
    !$acc enter data create(anb,sgndist) async(1)
    !
    !$acc parallel loop default(present) async(1)
    do p = 1,npmax
      anb(0,p)%x = ep(p)%x
      anb(0,p)%y = ep(p)%y
      anb(0,p)%z = ep(p)%z
    enddo
    !$acc end parallel
    !
    !$acc parallel loop collapse(2) default(present) async(1)
    do p=1,npmax
      do i=1,8
        anb(i,p)%x = 0._rp
        anb(i,p)%y = 0._rp
        anb(i,p)%z = 0._rp
      enddo
    enddo
    !$acc end parallel
    !
    !$acc wait
    !$acc update self(anb)
    !
    do p=1,pmax
      nrrequests = 0
      !
      do nb=1,8
        nbsend = nb    ! neighbor(nbsend) = rank of process to which data is send
        idp = abs(ep(p)%mslv)
        tag = idp*10+nbsend-idp*10
        nbrecv  = nb+4  ! neighbor(nbrecv)  = rank of process from which data is received
        !
        if (nbrecv > 8) nbrecv = nbrecv - 8
        !
        if (ep(p)%mslv > 0) then 
          ! myid is master of particle ap(p)%mslv
          !
          if (ep(p)%nb(nbsend) == 1) then
            !
            if (neighbor(nbsend) == myid) then
              ! process might be both master and slave of same particle due to periodic b.c.'s
              anb(nbrecv,p)%x = anb(0,p)%x
              anb(nbrecv,p)%y = anb(0,p)%y
              anb(nbrecv,p)%z = anb(0,p)%z
              ! recompute particle positions due to periodic b.c.'s
              if (nbrecv == 1) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
              endif
              !
              if (nbrecv == 2) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
              endif
              !
              if (nbrecv == 3) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
              endif
              !
              if (nbrecv == 4) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
              endif
              !
              if (nbrecv == 5) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
              endif
              !
              if (nbrecv == 6) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
                anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
              endif
              !
              if (nbrecv == 7) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
              endif
              !
              if (nbrecv == 8) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
                anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
              endif
            else
              ! neighbor(nbsend) is rank of slave for particle ep(p)%mslv
              nrrequests = nrrequests + 1
              call MPI_ISEND(anb(0,p)%x,3,MPI_REAL_RP,neighbor(nbsend), &
                             tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
              ! send x,y,z -> 3 contiguous info
              ! (see definition of type pneighbor in the begining of the subroutine)
            endif
          endif
        endif
        if (ep(p)%mslv < 0) then
          ! myid is slave of particle -ap(p)%mslv
          !
          if (ep(p)%nb(nbrecv) == 1) then
            ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
            nrrequests = nrrequests + 1
            call MPI_IRECV(anb(nbrecv,p)%x,3,MPI_REAL_RP,neighbor(nbrecv), &
                           tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
          endif
        endif
      enddo ! do nb=
      call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
    enddo
    !
    ! second step: recompute particle positions for slaves due to periodic b.c.'s.
    ! Required: (part of) particle within domain bounds of slave process.
    !
    !$omp parallel default(none) &
    !$omp shared(ep,anb,pmax,coords,l,dims)  &
    !$omp private(p,nbrecv,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)
    !$omp do
    do p=1,pmax
      !
      if (ep(p)%mslv < 0) then
        ! myid is slave of particle -ap(p)%mslv
        nbrecv=1
        if (ep(p)%nb(nbrecv) > 0) then
          ! neighbor(nbrecv) is rank of master of particle -ep(p)%mslv
          boundleftnb  = (coords(1)+1)*l(1)/(1._rp*dims(1)) ! left boundary of neighbor nb
          if (anb(nbrecv,p)%x < boundleftnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
          endif
        endif
        !
        nbrecv=2
        if (ep(p)%nb(nbrecv) > 0) then
          boundleftnb  = (coords(1)+1)*l(1)/(1._rp*dims(1)) ! left boundary of neighbor nb
          boundbacknb  = (coords(2))*l(2)/(1._rp*dims(2)) ! back boundary of neighbor nb
          if (anb(nbrecv,p)%x < boundleftnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
          endif
          if (anb(nbrecv,p)%y > boundbacknb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
          endif
        endif
        !
        nbrecv=3
        if (ep(p)%nb(nbrecv) > 0) then
          boundbacknb  = (coords(2))*l(2)/(1._rp*dims(2)) ! back boundary of neighbor nb
          if (anb(nbrecv,p)%y > boundbacknb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
          endif
        endif
        !
        nbrecv=4
        if (ep(p)%nb(nbrecv) > 0) then
          boundrightnb = (coords(1))*l(1)/(1._rp*dims(1)) ! right boundary of neighbor nb
          boundbacknb  = (coords(2))*l(2)/(1._rp*dims(2)) ! back  boundary of neighbor nb
          if (anb(nbrecv,p)%x > boundrightnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
          endif
          if (anb(nbrecv,p)%y > boundbacknb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
          endif
        endif
        !
        nbrecv=5
        if (ep(p)%nb(nbrecv) > 0) then
          boundrightnb = (coords(1))*l(1)/(1._rp*dims(1)) ! right boundary of neighbor nb
          if (anb(nbrecv,p)%x > boundrightnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
          endif
        endif
        !
        nbrecv=6
        if (ep(p)%nb(nbrecv) > 0) then
          boundrightnb = (coords(1))*l(1)/(1._rp*dims(1)) ! right boundary of neighbor nb
          boundfrontnb = (coords(2)+1)*l(2)/(1._rp*dims(2)) ! front boundary of neighbor nb
          if (anb(nbrecv,p)%x > boundrightnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
          endif
          if (anb(nbrecv,p)%y < boundfrontnb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
          endif
        endif
        !
        nbrecv=7
        if (ep(p)%nb(nbrecv) > 0) then
          boundfrontnb = (coords(2)+1)*l(2)/(1._rp*dims(2)) ! front boundary of neighbor nb
          if (anb(nbrecv,p)%y < boundfrontnb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
          endif
        endif
        !
        nbrecv=8
        if (ep(p)%nb(nbrecv) > 0) then
          boundleftnb  = (coords(1)+1)*l(1)/(1._rp*dims(1)) ! left  boundary of neighbor nb
          boundfrontnb = (coords(2)+1)*l(2)/(1._rp*dims(2)) ! front boundary of neighbor nb
          if (anb(nbrecv,p)%x < boundleftnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
          endif
          if (anb(nbrecv,p)%y < boundfrontnb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
          endif
        endif
        !
      endif
    enddo
    !$omp end parallel
    !
    ! second step: perform integration.
    !
    !$acc update device(anb) async(1)
    !
    do p=1,npmax
      do i=0,8
        anb(i,p)%intu = 0._rp
        anb(i,p)%intv = 0._rp
        anb(i,p)%intw = 0._rp
      enddo
    enddo
    do p=1,pmax
      ep(p)%fxltot = 0._rp
      ep(p)%fyltot = 0._rp
      ep(p)%fzltot = 0._rp
      ep(p)%torqxltot = 0._rp
      ep(p)%torqyltot = 0._rp
      ep(p)%torqzltot = 0._rp
      ! fcap* accumulates with '+' below and is scaled by dVeul at the end of the
      ! nb loop, so without this reset every step inherited dVeul*(previous step),
      ! i.e. fcap_n = dVeul*(fcap_{n-1} + sum_n). The spurious memory term is only
      ! ~dVeul = dx**3 in relative size (2.4e-4 at D = 32), but it grows as dx**3
      ! on coarser grids and made fcap* -- and hence the F_cap_ibm column of
      ! forces_data.csv -- carry a grid-dependent artefact.
      ep(p)%fcapx = 0._rp
      ep(p)%fcapy = 0._rp
      ep(p)%fcapz = 0._rp
    end do
    !
    !$omp parallel default(none) &
    !$omp shared(ep,anb,pmax) & 
    !$omp shared(unew,vnew,wnew) &
    !$omp shared(cas,dx,dy,dz,dli,n,radin2,dveul) &
    !$omp shared(myid,neighbor,boundleftmyid,boundfrontmyid) &
    !$omp private(p,nb,nbrecv,coorxc,cooryc,coorzc,ilow,ihigh,jlow,jhigh,klow,khigh) &
    !$omp private(coorxmin,coorxplus,coorymin,cooryplus,coorzmin,coorzplus) &
    !$omp private(dxm2,dxp2,dym2,dyp2,dzm2,dzp2,dist2,rx,ry,rz,sgndist,sum1,sum2) &
    !$omp private(i,j,k) 
    !$omp do
    !
    do p=1,pmax
      if (ep(p)%mslv /= 0) then
        ! myid is master or slave of particle abs(ap(p)%mslv)
        !
        do nb=0,8
          !
          if(( (nb > 0)  .and. (ep(p)%mslv < 0) .and. (ep(p)%nb(nb) == 1)) .or. & !slave
            (  (nb == 0) .and. (ep(p)%mslv > 0)) .or. & !pure master
            (  (nb > 0)  .and. (ep(p)%mslv > 0) .and. (ep(p)%nb(nb) == 1) .and. (neighbor(nb) == myid))) then
            ! master that looks like a slave due to periodic bcs
            nbrecv = nb
            ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
            !
            if ( neighbor(nb) == myid .and. nb > 0 .and. ep(p)%mslv > 0) then
              nbrecv = nb + 4
              if (nbrecv > 8) nbrecv = nbrecv-8
            endif
            !
            coorxc = anb(nbrecv,p)%x-boundleftmyid
            cooryc = anb(nbrecv,p)%y-boundfrontmyid
            coorzc = anb(nbrecv,p)%z
            !
            ilow  = nint( (coorxc-(radius-retrac))*dli(1) - eps_sol )
            ihigh = nint( (coorxc+(radius-retrac))*dli(1) + eps_sol )
            jlow  = nint( (cooryc-(radius-retrac))*dli(2) - eps_sol )
            jhigh = nint( (cooryc+(radius-retrac))*dli(2) + eps_sol )
            klow  = nint( (coorzc-(radius-retrac))*dli(3) - eps_sol )
            khigh = nint( (coorzc+(radius-retrac))*dli(3) + eps_sol )
            !
            if (ilow < 1) ilow = 1
            if (jlow < 1) jlow = 1
            if (klow < 1) klow = 1
            if (ihigh > n(1)) ihigh = n(1)
            if (jhigh > n(2)) jhigh = n(2)
            if (khigh > n(3)) khigh = n(3)
            !
            !$acc update device(intu,intv,intw) async(1)
            !
            !$acc parallel loop collapse(3) async(1) &
            !$acc reduction(+:intu,intv,intw) &
            !$acc default(present) &
            !$acc private(coorxmin,coorxplus,dxp2) &
            !$acc private(coorymin,cooryplus,dyp2) &
            !$acc private(coorzmin,coorzplus,dzp2) &
            !$acc private(dist2,ry,rz,dxm2,dym2,dzm2) &
            !$acc private(sgndist,sum1,sum2)            
            do k=klow,khigh
              do j=jlow,jhigh
                do i=ilow,ihigh
                  coorx_cent = boundleftmyid +(i-0.5_rp)*dx
                  coorx_stag = boundleftmyid + i*dx
                  dxp2_cent = (coorx_cent - anb(nbrecv,p)%x)**2
                  dxp2_stag = (coorx_stag - anb(nbrecv,p)%x)**2
                  !
                  coory_cent = boundfrontmyid + (j-0.5_rp)*dy
                  coory_stag= boundfrontmyid + j*dy
                  dyp2_cent = (coory_cent - anb(nbrecv,p)%y)**2
                  dyp2_stag = (coory_stag - anb(nbrecv,p)%y)**2
                  !
                  coorz_cent = (k-0.5_rp)*dz
                  coorz_stag = k*dz
                  dzp2_cent = (coorz_cent - anb(nbrecv,p)%z)**2
                  dzp2_stag = (coorz_stag - anb(nbrecv,p)%z)**2
                  !
                  dist2x = dxp2_stag+dyp2_cent+dzp2_cent
                  dist2y = dxp2_cent+dyp2_stag+dzp2_cent
                  dist2z = dxp2_cent+dyp2_cent+dzp2_stag
                  dist2c = dxp2_cent+dyp2_cent+dzp2_cent
                  !
                  normalx = [coorx_stag - anb(nbrecv,p)%x,  &
                             coory_cent - anb(nbrecv,p)%y,  &
                             coorz_cent - anb(nbrecv,p)%z    ]
                  normalx = normalx / sqrt(dist2x)
                  deltasx = sqrt(dist2x)-(radius-retrac)
                  call digitiser(deltasx,normalx,alpha_eulx)
                  !
                  normaly = [coorx_cent - anb(nbrecv,p)%x,  &
                             coory_stag - anb(nbrecv,p)%y,  &
                             coorz_cent - anb(nbrecv,p)%z    ]
                  normaly = normaly / sqrt(dist2y)
                  deltasy = sqrt(dist2y)-(radius-retrac)
                  call digitiser(deltasy,normaly,alpha_euly)
                  !
                  normalz = [coorx_cent - anb(nbrecv,p)%x,  &
                             coory_cent - anb(nbrecv,p)%y,  &
                             coorz_stag - anb(nbrecv,p)%z    ]
                  normalz = normalz / sqrt(dist2z)
                  deltasz = sqrt(dist2z)-(radius-retrac)
                  call digitiser(deltasz,normalz,alpha_eulz)
                  normalc = [coorx_cent - anb(nbrecv,p)%x,  &
                             coory_cent - anb(nbrecv,p)%y,  &
                             coorz_cent - anb(nbrecv,p)%z    ]
                  normalc = normalc / sqrt(dist2c)
                  norm_partx(i,j,k)=normalc(1)
                  norm_party(i,j,k)=normalc(2)
                  norm_partz(i,j,k)=normalc(3)
!                  deltasc = sqrt(dist2c)-(radius-retrac)
!                  call digitiser(deltasc,normalc,alpha_eulc)
!                  if (dist2c<radin2) then
!                    alphac(i,j,k) = 1._rp
!                  else
!                    alphac(i,j,k)=alpha_eulc
!                  end if
                  !
                  ul = ep(p)%u + ep(p)%omy*(coorz_cent-anb(nbrecv,p)%z) &
                               - ep(p)%omz*(coory_cent-anb(nbrecv,p)%y)
                  vl = ep(p)%v + ep(p)%omz*(coorx_cent-anb(nbrecv,p)%x) &
                               - ep(p)%omx*(coorz_cent-anb(nbrecv,p)%z)
                  wl = ep(p)%w + ep(p)%omx*(coory_cent-anb(nbrecv,p)%y) &
                               - ep(p)%omy*(coorx_cent-anb(nbrecv,p)%x)
                  !
                  rhox = rho + drho*0.5*(psi( i,j,k)+psi( i+1,j,k))
                  rhoy = rho + drho*0.5*(psi( i,j,k)+psi( i,j+1,k))
                  rhoz = rho + drho*0.5*(psi( i,j,k)+psi( i,j,k+1))
                  !
                  rhox_old = rho + drho*0.5*(psio( i,j,k)+psio( i+1,j,k))
                  rhoy_old = rho + drho*0.5*(psio( i,j,k)+psio( i,j+1,k))
                  rhoz_old = rho + drho*0.5*(psio( i,j,k)+psio( i,j,k+1))
                  !
                  fsurf_x = alpha_eulx*0.5_rp*(kappa(i+1,j,k)+kappa(i,j,k))* &
                            sigma*(2._rp/(rho12(1)+rho12(2)))*rhox*(psi(i+1,j,k)-psi(i,j,k))/dx
                  fsurf_y = alpha_euly*0.5_rp*(kappa(i,j+1,k)+kappa(i,j,k))* &
                            sigma*(2._rp/(rho12(1)+rho12(2)))*rhoy*(psi(i,j+1,k)-psi(i,j,k))/dy
                  fsurf_z = alpha_eulz*0.5_rp*(kappa(i,j,k+1)+kappa(i,j,k))* &
                            sigma*(2._rp/(rho12(1)+rho12(2)))*rhoz*(psi(i,j,k+1)-psi(i,j,k))/dz
                  !
                  if (abs(deltasx)>2*eps_sol*dx) then
                    alpha_eulx=0
                  end if
                  !
                  if (abs(deltasy)>2*eps_sol*dy) then
                    alpha_euly=0
                  end if
                  !
                  if (abs(deltasz)>2*eps_sol*dz) then
                    alpha_eulz=0
                  end if
                  !
                  fx = alpha_eulx*rhox*(ul-unew(i,j,k))*dti
                  ustar(i,j,k) = unew(i,j,k) + dt*(fx/rhox)
                  !
                  fy = alpha_euly*rhoy*(vl-vnew(i,j,k))*dti
                  vstar(i,j,k) = vnew(i,j,k) + dt*(fy/rhoy)
                  !
                  fz = alpha_eulz*rhoz*(wl-wnew(i,j,k))*dti
                  wstar(i,j,k) = wnew(i,j,k) + dt*(fz/rhoz)
                  !
                  !
!                  do s=1,3
!                  fx = fx + alpha_eulx*rhox*(ul-ustar(i,j,k))*dti
!                  fy = fy + alpha_euly*rhoy*(vl-vstar(i,j,k))*dti
!                  fz = fz + alpha_eulz*rhoz*(wl-wstar(i,j,k))*dti
!                  !
!                  ustar(i,j,k) = unew(i,j,k) + fx*dt/rhox
!                  vstar(i,j,k) = vnew(i,j,k) + fy*dt/rhoy
!                  wstar(i,j,k) = wnew(i,j,k) + fz*dt/rhoz
!                  end do
                  !
                  unew(i,j,k) = ustar(i,j,k)
                  vnew(i,j,k) = vstar(i,j,k)
                  wnew(i,j,k) = wstar(i,j,k)
                  !
                  ep(p)%fxltot = ep(p)%fxltot + fx
                  ep(p)%fyltot = ep(p)%fyltot + fy
                  ep(p)%fzltot = ep(p)%fzltot + fz
                  !
                  ep(p)%torqxltot = ep(p)%torqxltot + (coory_cent-anb(nbrecv,p)%y)*fz - &
                                                      (coorz_cent-anb(nbrecv,p)%z)*fy
                  ep(p)%torqyltot = ep(p)%torqyltot + (coorz_cent-anb(nbrecv,p)%z)*fx - &
                                                      (coorx_cent-anb(nbrecv,p)%x)*fz
                  ep(p)%torqzltot = ep(p)%torqzltot + (coorx_cent-anb(nbrecv,p)%x)*fy - &
                                                      (coory_cent-anb(nbrecv,p)%y)*fx
                  !
                  ep(p)%fcapx = ep(p)%fcapx + fsurf_x
                  ep(p)%fcapy = ep(p)%fcapy + fsurf_y
                  ep(p)%fcapz = ep(p)%fcapz + fsurf_z
                  !
                  fx_old(i,j,k)=fx
                  fy_old(i,j,k)=fy
                  fz_old(i,j,k)=fz
!                  if (myid ==4 .and. alpha_euly>0) then
!                    PRINT *, "alpha_eulx", alpha_eulx
!                    PRINT *, "alpha_euly", alpha_euly
!                    PRINT *, "alpha_eulz", alpha_eulz
!                    PRINT *, "fx", fx
!                    PRINT *, "fy", fy
!                    PRINT *, "fz", fz
!                  endif
                enddo ! do i=
                !
              enddo ! do j=
              !
            enddo ! do k=
            !$acc end parallel
            !
            !$acc wait
            ep(p)%fxltot = ep(p)%fxltot*dVeul
            ep(p)%fyltot = ep(p)%fyltot*dVeul
            ep(p)%fzltot = ep(p)%fzltot*dVeul
            ep(p)%torqxltot = ep(p)%torqxltot*dVeul
            ep(p)%torqyltot = ep(p)%torqyltot*dVeul
            ep(p)%torqzltot = ep(p)%torqzltot*dVeul
            ep(p)%torqtheta = (ep(p)%torqyltot*cos(ep(p)%phi)) - &
                              (ep(p)%torqxltot*sin(ep(p)%phi))
            ep(p)%fcapx = ep(p)%fcapx*dVeul
            ep(p)%fcapy = ep(p)%fcapy*dVeul
            ep(p)%fcapz = ep(p)%fcapz*dVeul
            !
          endif
          !
        enddo 
        !
      endif
      !
    enddo
    !$omp end parallel
    !
    !$acc exit data delete(anb)
    !
    ! third step: communicate data of slaves to their masters
    !
     do p=1,npmax
      do nb=1,8
        anb(nb,p)%intu   = 0._rp
        anb(nb,p)%intv   = 0._rp
        anb(nb,p)%intw   = 0._rp
        anb(nb,p)%intomx = 0._rp
        anb(nb,p)%intomy = 0._rp
        anb(nb,p)%intomz = 0._rp
      enddo
    enddo
    !
    anb(0,1:npmax)%intu   = ep(1:npmax)%fxltot 
    anb(0,1:npmax)%intv   = ep(1:npmax)%fyltot
    anb(0,1:npmax)%intw   = ep(1:npmax)%fzltot
    anb(0,1:npmax)%intomx = ep(1:npmax)%torqxltot
    anb(0,1:npmax)%intomy = ep(1:npmax)%torqyltot
    anb(0,1:npmax)%intomz = ep(1:npmax)%torqzltot
    !
    do p=1,pmax
      nrrequests = 0
      do nb=1,8
        nbsend = nb    ! rank of process which sends data ('data is received from neighbor nbsend')
        idp = abs(ep(p)%mslv)
        !tag = idp*10+nbsend
        tag = idp*10+nbsend-idp*10
        nbrecv  = nb+4  ! rank of process which receives data ('data is send to neighbor nbrecv')
        if (nbrecv > 8) nbrecv = nbrecv - 8
        if (ep(p)%mslv > 0) then
          ! myid is master of particle ep(p)%mslv
          if (ep(p)%nb(nbsend) == 1) then
            ! neighbor(nbsend) is rank of slave for particle ep(p)%mslv
            if ( neighbor(nbsend) /= myid ) then
              nrrequests = nrrequests + 1
              call MPI_IRECV(anb(nbrecv,p)%intu,6,MPI_REAL_RP,neighbor(nbsend), &
                             tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
            endif
          endif
        endif
        if (ep(p)%mslv < 0) then
          ! myid is slave of particle -ep(p)%mslv
          if (ep(p)%nb(nbrecv) == 1) then
            ! neighbor(nbrecv) is rank of master of particle -ep(p)%mslv
            nrrequests = nrrequests + 1
            call MPI_ISEND(anb(0,p)%intu,6,MPI_REAL_RP,neighbor(nbrecv), &
                           tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
          endif
        endif
      enddo ! do nb=
      call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
    enddo
    !
    if (rkit==1) then
      fx_tot(:) = 0._rp
      fy_tot(:) = 0._rp
      fz_tot(:) = 0._rp
      tx_tot(:) = 0._rp
      ty_tot(:) = 0._rp
      tz_tot(:) = 0._rp
    endif
    !$omp parallel default(none) &
    !$omp shared(ep,anb,pmax)    &
    !$omp private(p,nb) reduction(+:auxu,auxv,auxw,auxomx,auxomy,auxomz)
    !$omp do 
    do p=1,pmax
      ep(p)%fxltot    = 0._rp
      ep(p)%fyltot    = 0._rp
      ep(p)%fzltot    = 0._rp
      ep(p)%torqxltot = 0._rp
      ep(p)%torqyltot = 0._rp
      ep(p)%torqzltot = 0._rp
      !
      if (ep(p)%mslv > 0) then
        auxu   = 0._rp
        auxv   = 0._rp
        auxw   = 0._rp
        auxomx = 0._rp
        auxomy = 0._rp
        auxomz = 0._rp
        !
        do nb=0,8
          auxu   = auxu   + anb(nb,p)%intu
          auxv   = auxv   + anb(nb,p)%intv
          auxw   = auxw   + anb(nb,p)%intw
          auxomx = auxomx + anb(nb,p)%intomx
          auxomy = auxomy + anb(nb,p)%intomy
          auxomz = auxomz + anb(nb,p)%intomz
        enddo
        !
        ep(p)%fxltot    = auxu
        ep(p)%fyltot    = auxv
        ep(p)%fzltot    = auxw
        ep(p)%torqxltot = auxomx
        ep(p)%torqyltot = auxomy
        ep(p)%torqzltot = auxomz
        !
        fx_tot(p) = fx_tot(p) + ep(p)%fxltot*dt*dtitot
        fy_tot(p) = fy_tot(p) + ep(p)%fyltot*dt*dtitot
        fz_tot(p) = fz_tot(p) + ep(p)%fzltot*dt*dtitot
        tx_tot(p) = tx_tot(p) + ep(p)%torqxltot*dt*dtitot
        ty_tot(p) = ty_tot(p) + ep(p)%torqyltot*dt*dtitot
        tz_tot(p) = tz_tot(p) + ep(p)%torqzltot*dt*dtitot
!        PRINT *, "fx_tot", fx_tot(p)
!        PRINT *, "fy_tot", fy_tot(p)
!        PRINT *, "fz_tot", fz_tot(p)
!        PRINT *, "tx_tot", tx_tot(p)
!        PRINT *, "ty_tot", ty_tot(p)
!        PRINT *, "tz_tot", tz_tot(p)
      endif
    enddo
    !
    return
  end subroutine eulint
  !
#endif
#endif
end module prt_mod_eulint
