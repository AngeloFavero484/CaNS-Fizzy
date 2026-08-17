module prt_mod_phase_indicator
#if defined(_PARTICLE)
  use mpi
  use mod_types 
  use mod_common_mpi      , only: prt_comm_cart,ierr,myid, &
                                  boundfrontmyid,boundleftmyid
  use mod_param           , only: bcvel,bcpre,cbcvel,cbcpre,is_bound, &
                                  dims,l,dl,dli,nb
  use mod_bound           , only: boundp,bounduvw
  use prt_mod_param       , only: radius
  use prt_mod_common      , only: npmax,pmax,ep,neighbor,coords
  !
  implicit none
  private
  public phase_indicator
  !
  contains
  !
  subroutine phase_indicator(n,dzc,dzf,itype,gamu,gamv,gamw,gamp,upart,vpart,wpart)
    implicit none
    integer, intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    integer, intent(in) :: itype ! 1 -> impose rigid body motion -> translation + rotation
                                 ! 0 -> impose rigid body motion -> rotation
    real(rp), intent(out), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: gamu,gamv,gamw,gamp,upart,vpart,wpart
    type pneighbor
      real(rp) :: x,y,z,u,v,w,omx,omy,omz
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb ! can be pmax because it is a subroutine-specific array!
    integer :: i,j,k,p,cp
    integer :: ilow,ihigh,jlow,jhigh,klow,khigh
    real(rp) :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    real(rp) :: coorxc,cooryc,coorzc
    real(rp) :: coorxmin,coorxplus,coorymin,cooryplus,coorzmin,coorzplus
    real(rp) :: radin,radin2,dist2
    !real(rp) :: dx,dy,dz
    real(rp) :: sum1,sum2
    integer :: nb_ibm,nbsend,nbrecv
    integer :: nrrequests
    integer :: arrayrequests(1:3) ! p=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    real(rp)  :: sgndist(1:8)
    !real(rp) :: rx,ry,rz
    integer :: idp,tag
    !character(len=3) :: rankpr
    real(rp) :: dxm2,dxp2,dym2,dyp2,dzm2,dzp2
    !real(rp) :: aux,aux_all 
    !
!    dx       = 1./dxi
!    dy       = 1./dyi
!    dz       = 1./dzi
    radin    = radius-2.0_rp*dl(1)
    radin2   = radin**2
    !
    ! for r < radius - sqrt( (dx)**2 + (dy)**2 + (dz)**2 ) = 
    !        radius - 1.732050808*dx : all cell corner points within sphere
    !
    ! first step: slaves need from master x,..,xfp,.. etc
    !
    !$omp workshare
    anb(0,1:pmax)%x = ep(1:pmax)%x
    anb(0,1:pmax)%y = ep(1:pmax)%y
    anb(0,1:pmax)%z = ep(1:pmax)%z
    anb(0,1:pmax)%u = ep(1:pmax)%u
    anb(0,1:pmax)%v = ep(1:pmax)%v
    anb(0,1:pmax)%w = ep(1:pmax)%w
    anb(0,1:pmax)%omx = ep(1:pmax)%omx
    anb(0,1:pmax)%omy = ep(1:pmax)%omy
    anb(0,1:pmax)%omz = ep(1:pmax)%omz
    do p=1,pmax
       anb(1:8,p)%x = 0.0_rp
       anb(1:8,p)%y = 0.0_rp
       anb(1:8,p)%z = 0.0_rp
       anb(1:8,p)%u = 0.0_rp
       anb(1:8,p)%v = 0.0_rp
       anb(1:8,p)%w = 0.0_rp
       anb(1:8,p)%omx = 0.0_rp
       anb(1:8,p)%omy = 0.0_rp
       anb(1:8,p)%omz = 0.0_rp
    enddo
    !$omp end workshare
    !
    do p=1,pmax
      nrrequests = 0
      do nb_ibm=1,8
        nbsend = nb_ibm    ! neighbor(nbsend) = rank of process to which data is send
        idp = abs(ep(p)%mslv)
        tag = idp
        nbrecv  = nb_ibm+4  ! neighbor(nbrecv)  = rank of process from which data is received
        if (nbrecv > 8) nbrecv = nbrecv - 8
        if (ep(p)%mslv > 0) then 
          ! myid is master of particle ep(p)%mslv
          if (ep(p)%nb(nbsend) == 1) then
            if (neighbor(nbsend) == myid) then
              ! process might be both master and slave of same particle due to periodic b.c.'s
              anb(nbrecv,p)%x = anb(0,p)%x
              anb(nbrecv,p)%y = anb(0,p)%y
              anb(nbrecv,p)%z = anb(0,p)%z
              ! recompute particle positions due to periodic b.c.'s
              if (nbrecv == 1) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
              endif
              if (nbrecv == 2) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
              endif
              if (nbrecv == 3) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
              endif
              if (nbrecv == 4) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
              endif
              if (nbrecv == 5) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
              endif
              if (nbrecv == 6) then
                 anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
                 anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
              endif
              if (nbrecv == 7) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
              endif
              if (nbrecv == 8) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
                anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
              endif
            else
              ! neighbor(nbsend) is rank of slave for particle ep(p)%mslv
              nrrequests = nrrequests + 1
              call MPI_ISEND(anb(0,p)%x,9,MPI_REAL_RP,neighbor(nbsend),tag, &
                             prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
              ! send x,y,z -> 3 contiguous info
              ! (see definition of type pneighbor in the begining of the subroutine)
            endif
          endif
        endif
        if (ep(p)%mslv < 0) then 
          ! myid is slave of particle -ep(p)%mslv
          if (ep(p)%nb(nbrecv) == 1) then 
            ! neighbor(nbrecv) is rank of master of particle -ep(p)%mslv
            nrrequests = nrrequests + 1
            call MPI_IRECV(anb(nbrecv,p)%x,9,MPI_REAL_RP,neighbor(nbrecv),tag, &
                           prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
            ! recv x,y,z -> 3 contiguous info
            ! (see definition of type pneighbor in the begining of the subroutine)
          endif
        endif
      enddo ! do nb=
      call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
    enddo
    !
    ! second step: recompute particle positions for slaves due to periodic b.c.'s.
    ! Required: (part of) particle within domain bounds of slave process.
    !
    !$omp parallel default(shared) &
    !$omp private(p,nbrecv,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)
    !$omp do
    do p=1,pmax
      if (ep(p)%mslv < 0) then
        ! myid is slave of particle -ep(p)%mslv
        nbrecv=1
        if (ep(p)%nb(nbrecv) == 1) then
          ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
          boundleftnb  = (coords(1)+1)*l(1)/(1.0*dims(1)) ! left boundary of neighbor nb
          if (anb(nbrecv,p)%x < boundleftnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
          endif
        endif
        !
        nbrecv=2
        if (ep(p)%nb(nbrecv) > 0) then
          boundleftnb  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! left boundary of neighbor nb
          boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back boundary of neighbor nb
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
          boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back boundary of neighbor nb
          if (anb(nbrecv,p)%y > boundbacknb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - l(2)
          endif
        endif
        !
        nbrecv=4
        if (ep(p)%nb(nbrecv) > 0) then
          boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
          boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back  boundary of neighbor nb
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
           boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
           if (anb(nbrecv,p)%x > boundrightnb) then
              anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
           endif
        endif
        !
        nbrecv=6
        if (ep(p)%nb(nbrecv) > 0) then
          boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
          boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
          if (anb(nbrecv,p)%x > boundrightnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - l(1)
          endif
          if (anb(nbrecv,p)%y > boundfrontnb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
          endif
        endif
        !
        nbrecv=7
        if (ep(p)%nb(nbrecv) > 0) then
          boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
          if (anb(nbrecv,p)%y < boundfrontnb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
          endif
        endif
        !
        nbrecv=8
        if (ep(p)%nb(nbrecv) > 0) then
          boundleftnb  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! left  boundary of neighbor nb
          boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
          if (anb(nbrecv,p)%x < boundleftnb) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + l(1)
          endif
          if (anb(nbrecv,p)%y < boundfrontnb) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + l(2)
          endif
        endif
      endif
    enddo
    !$omp end parallel
    !
    ! third step: perform integration.
    !
    gamp(:,:,:) = 0.0_rp
    gamu(:,:,:) = 0.0_rp
    gamv(:,:,:) = 0.0_rp
    gamw(:,:,:) = 0.0_rp
    upart(:,:,:) = 0.0_rp
    vpart(:,:,:) = 0.0_rp
    wpart(:,:,:) = 0.0_rp
    !$omp parallel default(shared) &
    !$omp private(p,nb,nbrecv,coorxc,cooryc,coorzc,ilow,ihigh,jlow,jhigh,klow,khigh) &
    !$omp private(coorxmin,coorxplus,coorymin,cooryplus,coorzmin,coorzplus) &
    !$omp private(dxm2,dxp2,dym2,dyp2,dzm2,dzp2,dist2,rx,ry,rz,sgndist,sum1,sum2) &
    !$omp private(i,j,k) 
    !$omp do 
    do p=1,pmax
      if (ep(p)%mslv /= 0) then
        ! myid is master or slave of particle abs(ep(p)%mslv)
        do nb_ibm=0,8
          if ((nb_ibm > 0 .and. ep(p)%mslv < 0 .and. ep(p)%nb(nb_ibm) == 1) .or. & !slave
             (nb_ibm == 0 .and. ep(p)%mslv> 0) .or. & !pure master
             (nb_ibm > 0 .and. ep(p)%mslv > 0 .and. ep(p)%nb(nb_ibm) == 1 .and. neighbor(nb_ibm) == myid)) then
            ! master that looks like a slave due to periodic bcs
            nbrecv = nb_ibm
            ! neighbor(nbrecv) is rank of master of particle -ep(p)%mslv
            if ( neighbor(nb_ibm) == myid .and. nb_ibm > 0 .and. ep(p)%mslv > 0) then
              nbrecv = nb_ibm + 4
              if (nbrecv > 8) nbrecv = nbrecv-8
            endif
            coorxc = anb(nbrecv,p)%x-boundleftmyid
            cooryc = anb(nbrecv,p)%y-boundfrontmyid
            coorzc = anb(nbrecv,p)%z
            ilow  = nint( (coorxc-radius)*dli(1) - 2.0_rp)
            ihigh = nint( (coorxc+radius)*dli(1) + 2.0_rp)
            jlow  = nint( (cooryc-radius)*dli(2) - 2.0_rp)
            jhigh = nint( (cooryc+radius)*dli(2) + 2.0_rp)
            klow  = nint( (coorzc-radius)*dli(3) - 2.0_rp)
            khigh = nint( (coorzc+radius)*dli(3) + 2.0_rp)
            if (ilow < 1) ilow = 1
            if (jlow < 1) jlow = 1
            if (klow < 1) klow = 1
            if (ihigh > n(1)) ihigh = n(1)
            if (jhigh > n(2)) jhigh = n(2)
            if (khigh > n(3)) khigh = n(3)
            !
            ! u-velocity
            !
            do k=klow,khigh
              coorzmin  = (k-1)*dl(3)
              coorzplus = k*dl(3)
              do j=jlow,jhigh
                coorymin  = boundfrontmyid + (j-1)*dl(2)
                cooryplus = boundfrontmyid + j*dl(2)
                do i=ilow,ihigh
                  coorxmin  = boundleftmyid + (i-0.5_rp)*dl(1)
                  coorxplus = boundleftmyid + (i+0.5_rp)*dl(1)
                  dxp2 = (coorxplus - anb(nbrecv,p)%x)**2.0_rp
                  dyp2 = (cooryplus - anb(nbrecv,p)%y)**2.0_rp
                  dzp2 = (coorzplus - anb(nbrecv,p)%z)**2.0_rp
                  dist2 = dxp2+dyp2+dzp2 
                  if (dist2 < radin2 .and. itype == 1) then
                    gamu(i,j,k) = 1.0_rp
                    upart(i,j,k) = anb(nbrecv,p)%u + &
                                   anb(nbrecv,p)%omy*((k-0.5_rp)*dl(3)-anb(nbrecv,p)%z) - &
                                   anb(nbrecv,p)%omz*((j-0.5_rp)*dl(2)-anb(nbrecv,p)%y + boundfrontmyid)
                  else
                    dxm2 = (coorxmin - anb(nbrecv,p)%x)**2.0_rp
                    dym2 = (coorymin - anb(nbrecv,p)%y)**2.0_rp
                    dzm2 = (coorzmin - anb(nbrecv,p)%z)**2.0_rp
                    sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                    sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                    sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                    sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                    sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                    sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                    sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                    sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                    sum1 = 0.0_rp
                    sum2 = 0.0_rp
                    do cp=1,8 ! cp = corner point
                      if (sgndist(cp) > 0.0_rp) sum1 = sum1 + sgndist(cp)
                      sum2 = sum2 + abs( sgndist(cp) )
                    enddo
                    if(sum1/sum2 > gamu(i,j,k) .and. itype == 1) then
                      gamu(i,j,k) = sum1/sum2
                      upart(i,j,k) = anb(nbrecv,p)%u + &
                                     anb(nbrecv,p)%omy*((k-0.5_rp)*dl(3)-anb(nbrecv,p)%z) - &
                                     anb(nbrecv,p)%omz*((j-0.5_rp)*dl(2)-anb(nbrecv,p)%y + boundfrontmyid)
                    endif
                  endif
                enddo ! do i=
              enddo ! do j=
            enddo ! do k=
            !
            ! v-velocity
            !
            do k=klow,khigh
              coorzmin  = (k-1)*dl(3)
              coorzplus = k*dl(3)
              do j=jlow,jhigh
                coorymin  = boundfrontmyid + (j-0.5_rp)*dl(2)
                cooryplus = boundfrontmyid + (j+0.5_rp)*dl(2)
                do i=ilow,ihigh
                  coorxmin  = boundleftmyid + (i-1)*dl(1)
                  coorxplus = boundleftmyid + i*dl(1)
                  dxp2 = (coorxplus - anb(nbrecv,p)%x)**2.0_rp
                  dyp2 = (cooryplus - anb(nbrecv,p)%y)**2.0_rp
                  dzp2 = (coorzplus - anb(nbrecv,p)%z)**2.0_rp
                  dist2 = dxp2+dyp2+dzp2 
                  if (dist2 < radin2 .and. itype == 1) then
                    gamv(i,j,k) = 1.0_rp
                    vpart(i,j,k) = anb(nbrecv,p)%v + &
                                   anb(nbrecv,p)%omz*((i-0.5_rp)*dl(1)-anb(nbrecv,p)%x + boundleftmyid) - &
                                   anb(nbrecv,p)%omx*((k-0.5_rp)*dl(3)-anb(nbrecv,p)%z)
                  else
                    dxm2 = (coorxmin - anb(nbrecv,p)%x)**2.0_rp
                    dym2 = (coorymin - anb(nbrecv,p)%y)**2.0_rp
                    dzm2 = (coorzmin - anb(nbrecv,p)%z)**2.0_rp
                    sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                    sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                    sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                    sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                    sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                    sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                    sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                    sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                    sum1 = 0.0_rp
                    sum2 = 0.0_rp
                    do cp=1,8 ! cp = corner point
                      if (sgndist(cp) > 0.0_rp) sum1 = sum1 + sgndist(cp)
                      sum2 = sum2 + abs( sgndist(cp) )
                    enddo
                    if(sum1/sum2 > gamv(i,j,k) .and. itype == 1) then
                      gamv(i,j,k) = sum1/sum2
                      vpart(i,j,k) = anb(nbrecv,p)%v + &
                                     anb(nbrecv,p)%omz*((i-0.5_rp)*dl(1)-anb(nbrecv,p)%x + boundleftmyid) - &
                                     anb(nbrecv,p)%omx*((k-0.5_rp)*dl(2)-anb(nbrecv,p)%z)
                    endif
                  endif
                enddo ! do i=
              enddo ! do j=
            enddo ! do k=
            !
            ! w-velocity
            !
            do k=klow,khigh
              coorzmin  = (k-0.5_rp)*dl(3)
              coorzplus = (k+0.5_rp)*dl(3)
              do j=jlow,jhigh
                coorymin  = boundfrontmyid + (j-1)*dl(2)
                cooryplus = boundfrontmyid + j*dl(2)
                do i=ilow,ihigh
                  coorxmin  = boundleftmyid + (i-1)*dl(1)
                  coorxplus = boundleftmyid + i*dl(1)
                  dxp2 = (coorxplus - anb(nbrecv,p)%x)**2.0_rp
                  dyp2 = (cooryplus - anb(nbrecv,p)%y)**2.0_rp
                  dzp2 = (coorzplus - anb(nbrecv,p)%z)**2.0_rp
                  dist2 = dxp2+dyp2+dzp2
                  if (dist2 < radin2 .and. itype == 1) then
                    gamw(i,j,k) = 1.0_rp
                    wpart(i,j,k) = anb(nbrecv,p)%w + &
                                   anb(nbrecv,p)%omx*((j-0.5_rp)*dl(2)-anb(nbrecv,p)%y + boundfrontmyid) - &
                                   anb(nbrecv,p)%omy*((i-0.5_rp)*dl(1)-anb(nbrecv,p)%x + boundleftmyid)
                  else
                    dxm2 = (coorxmin - anb(nbrecv,p)%x)**2.0_rp
                    dym2 = (coorymin - anb(nbrecv,p)%y)**2.0_rp
                    dzm2 = (coorzmin - anb(nbrecv,p)%z)**2.0_rp
                    sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                    sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                    sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                    sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                    sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                    sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                    sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                    sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                    sum1 = 0.0_rp
                    sum2 = 0.0_rp
                    do cp=1,8 ! cp = corner point
                      if (sgndist(cp) > 0.0_rp) sum1 = sum1 + sgndist(cp)
                      sum2 = sum2 + abs( sgndist(cp) )
                    enddo
                    if(sum1/sum2 > gamw(i,j,k) .and. itype == 1) then
                      gamw(i,j,k) = sum1/sum2
                      wpart(i,j,k) = anb(nbrecv,p)%w + &
                                     anb(nbrecv,p)%omx*((j-0.5_rp)*dl(2)-anb(nbrecv,p)%y + boundfrontmyid) - &
                                     anb(nbrecv,p)%omy*((i-0.5_rp)*dl(1)-anb(nbrecv,p)%x + boundleftmyid)
                    endif
                  endif
                enddo ! do i=
              enddo ! do j=
            enddo ! do k=
            !
            ! pressure
            !
            do k=klow,khigh
              coorzmin  = (k-1)*dl(3)
              coorzplus = (k)*dl(3)
              do j=jlow,jhigh
                coorymin  = boundfrontmyid + (j-1)*dl(2)
                cooryplus = boundfrontmyid + j*dl(2)
                do i=ilow,ihigh
                  coorxmin  = boundleftmyid + (i-1)*dl(1)
                  coorxplus = boundleftmyid + i*dl(1)
                  dxp2 = (coorxplus - anb(nbrecv,p)%x)**2.0_rp
                  dyp2 = (cooryplus - anb(nbrecv,p)%y)**2.0_rp
                  dzp2 = (coorzplus - anb(nbrecv,p)%z)**2.0_rp
                  dist2 = dxp2+dyp2+dzp2 
                  if (dist2 < radin2 .and. itype == 1) then
                    gamp(i,j,k) = 1.0_rp
                  else
                    dxm2 = (coorxmin - anb(nbrecv,p)%x)**2.0_rp
                    dym2 = (coorymin - anb(nbrecv,p)%y)**2.0_rp
                    dzm2 = (coorzmin - anb(nbrecv,p)%z)**2.0_rp
                    sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                    sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                    sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                    sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                    sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                    sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                    sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                    sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                    sum1 = 0.0_rp
                    sum2 = 0.0_rp
                    do cp=1,8 ! cp = corner point
                      if (sgndist(cp) > 0.0_rp) sum1 = sum1 + sgndist(cp)
                      sum2 = sum2 + abs( sgndist(cp) )
                    enddo
                    if(itype == 1) gamp(i,j,k) = max(sum1/sum2,gamp(i,j,k)) 
                  endif
                  if (dist2 < radin2 .and. itype == 0) then
                    gamp(i,j,k) = 1.0_rp
                    upart(i,j,k) = anb(nbrecv,p)%omx
                    vpart(i,j,k) = anb(nbrecv,p)%omy
                    wpart(i,j,k) = anb(nbrecv,p)%omz
                  else
                    dxm2 = (coorxmin - anb(nbrecv,p)%x)**2.0_rp
                    dym2 = (coorymin - anb(nbrecv,p)%y)**2.0_rp
                    dzm2 = (coorzmin - anb(nbrecv,p)%z)**2.0_rp
                    sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                    sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                    sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                    sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                    sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                    sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                    sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                    sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                    sum1 = 0.0_rp
                    sum2 = 0.0_rp
                    do cp=1,8 ! cp = corner point
                      if (sgndist(cp) > 0.0_rp) sum1 = sum1 + sgndist(cp)
                      sum2 = sum2 + abs( sgndist(cp) )
                    enddo
                    if(sum1/sum2 > gamp(i,j,k) .and. itype == 0) then
                      gamp(i,j,k) = sum1/sum2
                      upart(i,j,k) = anb(nbrecv,p)%omx
                      vpart(i,j,k) = anb(nbrecv,p)%omy
                      wpart(i,j,k) = anb(nbrecv,p)%omz
                    endif
                  endif
                enddo ! do i=
              enddo ! do j=
            enddo ! do k=
          endif
        enddo ! do nbrecv=
      endif
    enddo
    !$omp end parallel
    gamp(:,:,:) = 1.0_rp - gamp(:,:,:)
    gamu(:,:,:) = 1.0_rp - gamu(:,:,:)
    gamv(:,:,:) = 1.0_rp - gamv(:,:,:)
    gamw(:,:,:) = 1.0_rp - gamw(:,:,:)
    !
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,gamu,gamv,gamw)
    call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,gamp)
    !
    call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,upart,vpart,wpart)
!    call updthalos(gamp,1)
!    call updthalos(gamu,1)
!    call updthalos(gamv,1)
!    call updthalos(gamw,1)
!    call updthalos(gamp,2)
!    call updthalos(gamu,2)
!    call updthalos(gamv,2)
!    call updthalos(gamw,2)
!    call updthalos(upart,1)
!    call updthalos(vpart,1)
!    call updthalos(wpart,1)
!    call updthalos(upart,2)
!    call updthalos(vpart,2)
!    call updthalos(wpart,2)
!    !
!    do j=0,j1
!       do i=0,i1
!          gamu(i,j,0)  = gamu(i,j,1)
!          gamv(i,j,0)  = gamv(i,j,1)
!          gamw(i,j,0)  = gamw(i,j,1)
!          gamp(i,j,0)  = gamp(i,j,1)
!          gamu(i,j,k1) = gamu(i,j,kmax)
!          gamv(i,j,k1) = gamv(i,j,kmax)
!          gamw(i,j,k1) = gamw(i,j,kmax)
!          gamp(i,j,k1) = gamp(i,j,kmax)
!       enddo
!    enddo
!    do j=0,j1
!       do i=0,i1
!          upart(i,j,0)    = -upart(i,j,1)     ! no-slip
!          vpart(i,j,0)    = -vpart(i,j,1)     ! no-slip
!          wpart(i,j,0)    = 0.                ! no-penetration
!          upart(i,j,k1)   = -upart(i,j,kmax)  ! no-slip
!          vpart(i,j,k1)   = -vpart(i,j,kmax)  ! no-slip
!          wpart(i,j,kmax) = 0.                ! no-penetration
!          wpart(i,j,k1)   = wpart(i,j,kmax-1) ! dw/dz=0 at wall (not used) 
!       enddo
!    enddo
    !
    !aux = sum(1.-gamp(1:imax,1:jmax,1:kmax))
    !call mpi_allreduce(aux,aux_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    !if(myid.eq.0) print*,'Check Phase Indicator p:',np*volp,aux_all*dveul, &
    !                                               (np*volp-aux_all*dveul)/(np*volp)
    !aux = sum(1.-gamu(1:imax,1:jmax,1:kmax))
    !call mpi_allreduce(aux,aux_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    !if(myid.eq.0) print*,'Check Phase Indicator u:',np*volp,aux_all*dveul, &
    !                                               (np*volp-aux_all*dveul)/(np*volp)
    !aux = sum(1.-gamv(1:imax,1:jmax,1:kmax))
    !call mpi_allreduce(aux,aux_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    !if(myid.eq.0) print*,'Check Phase Indicator v:',np*volp,aux_all*dveul, &
    !                                               (np*volp-aux_all*dveul)/(np*volp)
    !aux = sum(1.-gamw(1:imax,1:jmax,1:kmax))
    !call mpi_allreduce(aux,aux_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    !if(myid.eq.0) print*,'Check Phase Indicator w:',np*volp,aux_all*dveul, &
    !                                               (np*volp-aux_all*dveul)/(np*volp)
    !
    return
  end subroutine phase_indicator
  !
#endif
end module prt_mod_phase_indicator
