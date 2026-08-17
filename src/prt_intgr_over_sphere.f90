module prt_mod_intgr_over_sphere
#if defined(_PARTICLE)
  use mod_types
  use mpi
  use mod_common_mpi    , only: boundfrontmyid,boundleftmyid,prt_comm_cart,ierr,myid
  use mod_param         , only: l,dli,dims,rho12
  use prt_mod_param     , only: volp,radius,eps_sol
  use prt_mod_common    , only: ep,npmax,pmax,dVeul,neighbor,coords
  use prt_mod_digitiser , only: digitiser
  !
  implicit none
  !
  private
  public :: intgr_over_sphere
  !
  contains
  !
  subroutine intgr_over_sphere(cas,n,psi,unew,vnew,wnew) 
    implicit none
    integer, intent(in) :: cas
    integer, dimension(3), intent(in) :: n
    real(rp), dimension(0:,0:,0:), intent(in) :: psi
    real(rp), dimension(0:,0:,0:), intent(in) :: unew,vnew,wnew
    type pneighbor
       real(rp) :: x,y,z, &
                   intu,intv,intw
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb ! can be pmax because it is a subroutine-specific array!
    integer :: i,j,k,p,cp
    integer :: ilow,ihigh,jlow,jhigh,klow,khigh
    real(rp) :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    real(rp) :: coorxc,cooryc,coorzc
    real(rp) :: coorxmin,coorxplus,coorymin,cooryplus,coorzmin,coorzplus
    real(rp) :: radin,radin2,dist2
    real(rp) :: dx,dy,dz
    real(rp) :: sum1,sum2
    integer :: nb,nbsend,nbrecv
    integer :: nrrequests
    integer :: arrayrequests(1:3) ! 3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    real(rp) :: sgndist(1:8)
    real(rp) :: rx,ry,rz
    integer :: idp,tag
    character(len=5) rankpr
    real(rp) :: dxm2,dxp2,dym2,dyp2,dzm2,dzp2
    real(rp) :: auxu,auxv,auxw
    real(rp) :: intu,intv,intw
    real(rp) :: alpha_eul
    real(rp) :: deltas
    real(rp) :: normal(1:3)
    real(rp) :: rho,drho,rhox,rhoy,rhoz
    character(len=64) :: csv_filename
    character(len=5)  :: str_rank
    integer, parameter :: csv_unit = 101
!    real(rp) :: buf_send(3),buf_recv(3)
!    real(rp), dimension(1:3,0:8,1:npmax) :: buf_send,buf_recv
    !
    dx       = 1./dli(1)
    dy       = 1./dli(2)
    dz       = 1./dli(3)
    radin    = radius-eps_sol*dx
    radin2   = radin**2
    rho=rho12(2)
    drho=rho12(1)-rho12(2)
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
    !$acc enter data create(intu,intv,intw) async(1)
    !
    do p=1,npmax
      do i=0,8
        anb(i,p)%intu = 0._rp
        anb(i,p)%intv = 0._rp
        anb(i,p)%intw = 0._rp
      enddo
    enddo
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
            ilow  = nint( (coorxc-radius)*dli(1) - eps_sol )
            ihigh = nint( (coorxc+radius)*dli(1) + eps_sol )
            jlow  = nint( (cooryc-radius)*dli(2) - eps_sol )
            jhigh = nint( (cooryc+radius)*dli(2) + eps_sol )
            klow  = nint( (coorzc-radius)*dli(3) - eps_sol )
            khigh = nint( (coorzc+radius)*dli(3) + eps_sol )
            !
            if (ilow < 1) ilow = 1
            if (jlow < 1) jlow = 1
            if (klow < 1) klow = 1
            if (ihigh > n(1)) ihigh = n(1)
            if (jhigh > n(2)) jhigh = n(2)
            if (khigh > n(2)) khigh = n(3)
            !
            intu = 0.0_rp
            intv = 0.0_rp
            intw = 0.0_rp
            !$acc update device(intu,intv,intw) async(1)
            !
            ! u-velocity
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
#if !defined(_EULER)
                  coorzmin  = (k-1)*dz
                  coorzplus = k*dz
#else
                  coorzplus = (k-0.5_rp)*dz
#endif
                  dzp2 = (coorzplus - anb(nbrecv,p)%z)**2
                  !
#if !defined(_EULER)
                  coorymin  = boundfrontmyid + (j-1)*dy
                  cooryplus = boundfrontmyid + j*dy
#else
                  cooryplus= boundfrontmyid + (j-0.5_rp)*dy
#endif
                  dyp2 = (cooryplus - anb(nbrecv,p)%y)**2
                  !
#if !defined(_EULER)
                  coorxmin  = boundleftmyid + (i-0.5_rp)*dx
                  coorxplus = boundleftmyid + (i+0.5_rp)*dx
#else
                  coorxplus = boundleftmyid + i*dx
#endif
                  dxp2 = (coorxplus - anb(nbrecv,p)%x)**2
                  dist2 = dxp2+dyp2+dzp2 
                  rhox = rho + drho*0.5*(psi( i,j,k)+psi( i+1,j,k))
                  !
                  if (dist2 < radin2) then
                    if (cas == 1) then
                      intu = intu + dVeul
                    elseif(cas == 2) then
                      intu = intu + rhox*unew(i,j,k)*dVeul
                    elseif(cas == 3) then
                      ry = boundfrontmyid + (j-0.5_rp)*dy - anb(nbrecv,p)%y
                      rz = (k-0.5_rp)*dz                  - anb(nbrecv,p)%z
                      intv = intv + rz*rhox*unew(i,j,k)*dVeul
                      intw = intw - ry*rhox*unew(i,j,k)*dVeul
                    elseif(cas ==4) then
                      intu = intu + rhox*dVeul
                    endif
                  else
#if !defined(_EULER)
                    dxm2 = (coorxmin - anb(nbrecv,p)%x)**2
                    dym2 = (coorymin - anb(nbrecv,p)%y)**2
                    dzm2 = (coorzmin - anb(nbrecv,p)%z)**2
                    sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                    sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                    sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                    sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                    sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                    sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                    sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                    sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                    sum1 = 0._rp
                    sum2 = 0._rp
                    !
                    do cp=1,8 ! cp = corner point
                      if (sgndist(cp) > 0._rp) sum1 = sum1 + sgndist(cp)
                      sum2 = sum2 + abs( sgndist(cp) )
                    enddo
#else
                    normal = [coorxplus - anb(nbrecv,p)%x,  &
                              cooryplus - anb(nbrecv,p)%y,  &
                              coorzplus - anb(nbrecv,p)%z    ]
                    normal = normal / sqrt(dist2)
                    deltas = sqrt(dist2)-radius
                    call digitiser(deltas,normal,alpha_eul)
#endif
                    if (cas == 1) then
#if defined(_EULER)
                      intu = intu + alpha_eul*dVeul
#else
                      intu = intu + (sum1/sum2)*dVeul
#endif
                    elseif(cas == 2) then
#if defined(_EULER)
                      intu = intu + rhox*unew(i,j,k)*alpha_eul*dVeul
#else
                      intu = intu + rhox*unew(i,j,k)*(sum1/sum2)*dVeul
#endif
                    elseif(cas == 3) then
#if defined(_EULER)
                      ry = cooryplus - anb(nbrecv,p)%y
                      rz = coorzplus - anb(nbrecv,p)%z
                      intv = intv + rz*rhox*unew(i,j,k)*alpha_eul*dVeul
                      intw = intw - ry*rhox*unew(i,j,k)*alpha_eul*dVeul
#else
                      ry = boundfrontmyid + (j-0.5_rp)*dy - anb(nbrecv,p)%y
                      rz = (k-0.5_rp)*dz                  - anb(nbrecv,p)%z
                      intv = intv + rz*rhox*unew(i,j,k)*(sum1/sum2)*dVeul
                      intw = intw - ry*rhox*unew(i,j,k)*(sum1/sum2)*dVeul
#endif
                    elseif (cas == 4) then
#if defined(_EULER)
                      intu = intu + rhox*alpha_eul*dVeul
#else
                      intu = intu + rhox*(sum1/sum2)*dVeul
#endif
                    endif
                    !
                  endif
                  !
                enddo ! do i=
                !
              enddo ! do j=
              !
            enddo ! do k=
            !$acc end parallel
            !
            !
            ! v-velocity
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
#if !defined(_EULER)
                  coorzmin  = (k-1)*dz
                  coorzplus = k*dz
#else
                  coorzplus=(k-0.5_rp)*dz
#endif
                  dzp2 = (coorzplus - anb(nbrecv,p)%z)**2
                  !
#if !defined(_EULER)
                  coorymin  = boundfrontmyid + (j-0.5_rp)*dy
                  cooryplus = boundfrontmyid + (j+0.5_rp)*dy
#else
                  cooryplus = boundfrontmyid + j*dy
#endif
                  dyp2 = (cooryplus - anb(nbrecv,p)%y)**2
                  !
#if !defined(_EULER)
                  coorxmin  = boundleftmyid + (i-1)*dx
                  coorxplus = boundleftmyid + i*dx
#else
                  coorxplus = boundleftmyid + (i-0.5_rp)*dx
#endif
                  dxp2 = (coorxplus - anb(nbrecv,p)%x)**2
                  dist2 = dxp2+dyp2+dzp2
                  rhoy = rho + drho*0.5*(psi( i,j,k)+psi( i,j+1,k))
                  !
                  if (dist2 < radin2) then
                    if (cas == 1) then
                      intv = intv + dVeul
                    elseif(cas == 2) then
                      intv = intv + rhoy*vnew(i,j,k)*dVeul
                    elseif(cas == 3) then
                      rx = boundleftmyid + (i-0.5_rp)*dx - anb(nbrecv,p)%x
                      rz = (k-0.5_rp)*dz                 - anb(nbrecv,p)%z
                      intu = intu - rz*rhoy*vnew(i,j,k)*dVeul
                      intw = intw + rx*rhoy*vnew(i,j,k)*dVeul
                    elseif(cas == 4) then
                      intv = intv + rhoy*dVeul
                    endif
                  else
#if !defined(_EULER)
                    dxm2 = (coorxmin - anb(nbrecv,p)%x)**2
                    dym2 = (coorymin - anb(nbrecv,p)%y)**2
                    dzm2 = (coorzmin - anb(nbrecv,p)%z)**2
                    sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                    sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                    sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                    sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                    sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                    sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                    sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                    sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                    sum1 = 0._rp
                    sum2 = 0._rp
                    !
                    !$acc loop seq
                    do cp=1,8 ! cp = corner point
                      if (sgndist(cp) > 0._rp) sum1 = sum1 + sgndist(cp)
                      sum2 = sum2 + abs( sgndist(cp) )
                    enddo
#else
                    normal = [coorxplus - anb(nbrecv,p)%x,  &
                              cooryplus - anb(nbrecv,p)%y,  &
                              coorzplus - anb(nbrecv,p)%z    ]
                    normal = normal / sqrt(dist2)
                    deltas = sqrt(dist2)-radius
                    call digitiser(deltas,normal,alpha_eul)
#endif
                    if (cas == 1) then
#if defined(_EULER)
                      intv = intv + alpha_eul*dVeul
#else
                      intv = intv + (sum1/sum2)*dVeul
#endif
                    elseif(cas == 2) then
#if defined(_EULER)
                      intv = intv + rhoy*vnew(i,j,k)*alpha_eul*dVeul
#else
                      intv = intv + rhoy*vnew(i,j,k)*(sum1/sum2)*dVeul
#endif
                    elseif(cas == 3) then
#if defined(_EULER)
                      rx = coorxplus - anb(nbrecv,p)%x
                      rz = coorzplus - anb(nbrecv,p)%z 
                      intu = intu - rz*rhoy*vnew(i,j,k)*alpha_eul*dVeul
                      intw = intw + rx*rhoy*vnew(i,j,k)*alpha_eul*dVeul
#else
                      rx = boundleftmyid + (i-0.5)*dx - anb(nbrecv,p)%x
                      rz = (k-0.5)*dz                 - anb(nbrecv,p)%z
                      intu = intu - rz*rhoy*vnew(i,j,k)*(sum1/sum2)*dVeul
                      intw = intw + rx*rhoy*vnew(i,j,k)*(sum1/sum2)*dVeul
#endif
                    elseif (cas == 4) then
#if defined(_EULER)
                      intv = intv + rhoy*alpha_eul*dVeul
#else
                      intv = intv + rhoy*(sum1/sum2)*dVeul
#endif
                    endif
                    !
                  endif
                  !
                enddo ! do i=
                !
              enddo ! do j=
              !
            enddo ! do k=
            !$acc end parallel
            !
            ! w-velocity
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
#if !defined(_EULER)
                  coorzmin  = (k-0.5_rp)*dz
                  coorzplus = (k+0.5_rp)*dz
#else
                  coorzplus = k*dz
#endif
                  dzp2 = (coorzplus - anb(nbrecv,p)%z)**2
#if !defined(_EULER)
                  coorymin  = boundfrontmyid + (j-1)*dy
                  cooryplus = boundfrontmyid + j*dy
#else
                  cooryplus = boundfrontmyid + (j-0.5_rp)*dy
#endif
                  dyp2 = (cooryplus - anb(nbrecv,p)%y)**2
                  
#if !defined(_EULER)
                  coorxmin  = boundleftmyid + (i-1)*dx
                  coorxplus = boundleftmyid + i*dx
#else
                  coorxplus = boundleftmyid + (i-0.5_rp)*dx
#endif
                  dxp2 = (coorxplus - anb(nbrecv,p)%x)**2
                  dist2 = dxp2+dyp2+dzp2
                  rhoz = rho + drho*0.5*(psi( i,j,k)+psi( i,j,k+1))
                  !
                  if (dist2 < radin2) then
                    if (cas == 1) then
                      intw = intw + dVeul
                    elseif(cas == 2) then
                      intw = intw + rhoz*wnew(i,j,k)*dVeul
                    elseif(cas == 3) then
                      rx = boundleftmyid  + (i-0.5)*dx - anb(nbrecv,p)%x
                      ry = boundfrontmyid + (j-0.5)*dy - anb(nbrecv,p)%y
                      intu = intu + ry*rhoz*wnew(i,j,k)*dVeul
                      intv = intv - rx*rhoz*wnew(i,j,k)*dVeul
                    elseif(cas == 4) then
                      intw = intw + rhoz*dVeul
                    endif
                  else
#if !defined(_EULER)
                    dxm2 = (coorxmin - anb(nbrecv,p)%x)**2
                    dym2 = (coorymin - anb(nbrecv,p)%y)**2
                    dzm2 = (coorzmin - anb(nbrecv,p)%z)**2
                    sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                    sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                    sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                    sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                    sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                    sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                    sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                    sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                    sum1 = 0._rp
                    sum2 = 0._rp
                    !
                    !$acc loop seq
                    do cp=1,8 ! cp = corner point
                      if (sgndist(cp) > 0._rp) sum1 = sum1 + sgndist(cp)
                      sum2 = sum2 + abs( sgndist(cp) )
                    enddo
#else
                    normal = [coorxplus - anb(nbrecv,p)%x,  &
                              cooryplus - anb(nbrecv,p)%y,  &
                              coorzplus - anb(nbrecv,p)%z    ]
                    normal = normal / sqrt(dist2)
                    deltas = sqrt(dist2)-radius
                    call digitiser(deltas,normal,alpha_eul)
#endif
                    if (cas == 1) then
#if defined(_EULER)
                      intw = intw + alpha_eul*dVeul
#else
                      intw = intw + (sum1/sum2)*dVeul
#endif
                    elseif(cas == 2) then
#if defined(_EULER)
                      intw = intw + rhoz*wnew(i,j,k)*alpha_eul*dVeul
#else
                      intw = intw + rhoz*wnew(i,j,k)*(sum1/sum2)*dVeul
#endif
                    elseif(cas == 3) then
#if defined(_EULER)
                      rx = coorxplus - anb(nbrecv,p)%x
                      ry = cooryplus - anb(nbrecv,p)%y 
                      intu = intu + ry*rhoz*wnew(i,j,k)*alpha_eul*dVeul
                      intv = intv - rx*rhoz*wnew(i,j,k)*alpha_eul*dVeul
#else
                      rx = boundleftmyid  + (i-0.5_rp)*dx - anb(nbrecv,p)%x
                      ry = boundfrontmyid + (j-0.5_rp)*dy - anb(nbrecv,p)%y
                      intu = intu + ry*rhoz*wnew(i,j,k)*(sum1/sum2)*dVeul
                      intv = intv - rx*rhoz*wnew(i,j,k)*(sum1/sum2)*dVeul
#endif
                    elseif (cas == 4) then
#if defined(_EULER)
                      intw = intw + rhoz*alpha_eul*dVeul
#else
                      intw = intw + rhoz*(sum1/sum2)*dVeul
#endif
                    endif
                    !
                  endif
                  !
                enddo ! do i=
                !
              enddo ! do j=
              !
            enddo ! do k=
            !$acc end parallel
            !
            !$acc wait
            !$acc update self(intu,intv,intw)
            anb(nbrecv,p)%intu = intu
            anb(nbrecv,p)%intv = intv
            anb(nbrecv,p)%intw = intw
            !
          endif
          !
        enddo ! do nbrecv=
        !
      endif
      !
    enddo
    !$omp end parallel
    !
    !
    !$acc exit data delete(intu,intv,intw) async(1)
    !
    ! third step: communicate data of slaves to their masters
    !
    do p=1,pmax
      nrrequests = 0
      !
      do nb=1,8
        nbsend = nb    ! rank of process which sends data ('data is received from neighbor nbsend')
        idp = abs(ep(p)%mslv)
        !tag = idp*10+nbsend
        tag = idp*10+nbsend-idp*10
        nbrecv  = nb+4  ! rank of process which receives data ('data is send to neighbor nbrecv')
        if (nbrecv > 8) nbrecv = nbrecv - 8
        !
        if (ep(p)%mslv > 0) then
          ! myid is master of particle ep(p)%mslv
          !
          if (ep(p)%nb(nbsend) == 1) then
            ! neighbor(nbsend) is rank of slave for particle ep(p)%mslv
            !
            if ( neighbor(nbsend) /= myid ) then
              nrrequests = nrrequests + 1
              call MPI_IRECV(anb(nbsend,p)%intu,3,MPI_REAL_RP,neighbor(nbsend), &
                             tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
            endif
            !
          endif
          !
        endif
        !
        if (ep(p)%mslv < 0) then
          ! myid is slave of particle -ep(p)%mslv
          !
          if (ep(p)%nb(nbrecv) == 1) then
            ! neighbor(nbrecv) is rank of master of particle -ep(p)%mslv
            nrrequests = nrrequests + 1
            call MPI_ISEND(anb(nbrecv,p)%intu,3,MPI_REAL_RP,neighbor(nbrecv), &
                           tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
            ! send intu,intv,intwx,y,z -> 3 contiguous info
            ! (see definition of type pneighbor in the begining of the subroutine)
          endif
          !
        endif
        !
      enddo ! do nb=
      call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
      !
    enddo
    !
    !$acc update device(anb)
    !$acc wait
    !
    ! Sum all contributions together.
    !
    !$omp parallel default(none) &
    !$omp shared(ep,anb,pmax,cas) &
    !$omp private(p,nb) reduction(+:auxu,auxv,auxw)
    SELECT CASE (cas)
    CASE (1,2)
      !$omp do
      !$acc parallel loop async(1) &
      !$acc default(present) &
      !$acc private(auxu,auxv,auxw)
      do p=1,pmax
        if (ep(p)%mslv > 0) then
          auxu = 0._rp
          auxv = 0._rp
          auxw = 0._rp
          !
          !$acc loop &
          !$acc reduction(+:auxu,auxv,auxw)
          do nb=0,8
            auxu = auxu + anb(nb,p)%intu
            auxv = auxv + anb(nb,p)%intv
            auxw = auxw + anb(nb,p)%intw
          enddo
          !
          ep(p)%intu = auxu
          ep(p)%intv = auxv
          ep(p)%intw = auxw
        else
          ep(p)%intu = 0._rp
          ep(p)%intv = 0._rp
          ep(p)%intw = 0._rp
        end if
      enddo
      !$acc end parallel
    CASE (3)
      !$omp do
      !$acc parallel loop async(1) &
      !$acc default(present) &
      !$acc private(auxu,auxv,auxw)
      do p=1,pmax
        if (ep(p)%mslv > 0) then
          auxu = 0._rp
          auxv = 0._rp
          auxw = 0._rp
          !
          !$acc loop &
          !$acc reduction(+:auxu,auxv,auxw)
          do nb=0,8
            auxu = auxu + anb(nb,p)%intu
            auxv = auxv + anb(nb,p)%intv
            auxw = auxw + anb(nb,p)%intw
          enddo
          !
          ep(p)%intomx = auxu
          ep(p)%intomy = auxv
          ep(p)%intomz = auxw
        else
          ep(p)%intomx = 0._rp
          ep(p)%intomy = 0._rp
          ep(p)%intomz = 0._rp
        end if
      enddo
      CASE (4)
      !$omp do
      !$acc parallel loop async(1) &
      !$acc default(present) &
      !$acc private(auxu,auxv,auxw)
      do p=1,pmax
        if (ep(p)%mslv > 0) then
          auxu = 0._rp
          auxv = 0._rp
          auxw = 0._rp
          !
          !$acc loop &
          !$acc reduction(+:auxu,auxv,auxw)
          do nb=0,8
            auxu = auxu + anb(nb,p)%intu
            auxv = auxv + anb(nb,p)%intv
            auxw = auxw + anb(nb,p)%intw
          enddo
          !
          ep(p)%intrhox = auxu
          ep(p)%intrhoy = auxv
          ep(p)%intrhoz = auxw
        else
          ep(p)%intrhox = 0._rp
          ep(p)%intrhoy = 0._rp
          ep(p)%intrhoz = 0._rp
        end if
      enddo
      !$acc end parallel
    END SELECT
    !$omp end parallel
    !
    if (cas == 1) then
      !write(rankpr,'(i5.5)') myid
      do p=1,pmax
        if (ep(p)%mslv > 0) then
          idp = ep(p)%mslv
          !    open(22,file=datadir//'volsphr'//rankpr//'.txt',position='append')
          !    write(22,'(2I8,4E16.8)') myid,idp,intu(p),intv(p),intw(p),Volp
          !    close(22)
          if (idp == 1) then
            write(6,'(A31,E16.8)') 'Calc. value of volume sphere = ',ep(p)%intu
            write(6,'(A31,E16.8)') 'Exact value of volume sphere = ',Volp
            write(6,'(A41,E16.8)') 'Error in calc. of volume sphere (in %) = ',100._rp*(ep(p)%intu-Volp)/Volp
          endif
        endif
      enddo
    endif
!    if (cas == 4) then
!      !write(rankpr,'(i5.5)') myid
!      do p=1,pmax
!        if (ep(p)%mslv > 0) then
!          idp = ep(p)%mslv
!          !    open(22,file=datadir//'volsphr'//rankpr//'.txt',position='append')
!          !    write(22,'(2I8,4E16.8)') myid,idp,intu(p),intv(p),intw(p),Volp
!          !    close(22)
!          if (idp == 1) then
!            PRINT *, 'Mass of the fluid occupied by the particle = ',ep(p)%intrhox,ep(p)%intrhoy,ep(p)%intrhoz
!!            write(6,'(A31,E16.8)') 'Exact value of volume sphere = ',Volp
!!            write(6,'(A41,E16.8)') 'Error in calc. of volume sphere (in %) = ',100._rp*(ep(p)%intu-Volp)/Volp
!          endif
!        endif
!      enddo
!    endif
    !
    return
  end subroutine intgr_over_sphere
  !
#endif
end module prt_mod_intgr_over_sphere
