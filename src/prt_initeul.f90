module prt_mod_initeul
#if defined(_PARTICLE)
#if defined(_EULER)
  use mod_types
  use mpi
  use mod_common_mpi    , only: boundfrontmyid,boundleftmyid,prt_comm_cart,ierr,myid
  use mod_param         , only: l,dli,dims,nh_wide
  use prt_mod_param     , only: radius,eps_sol
  use prt_mod_common    , only: ep,npmax,pmax,dVeul,neighbor,coords,retrac,alphac, &
                                norm_partx,norm_party,norm_partz
  use prt_mod_digitiser , only: digitiser
  !
  implicit none
  !
  private
  public :: initeul
  !
  contains
  !
  subroutine initeul(n)
    implicit none
    integer, dimension(3), intent(in) :: n
    type pneighbor
       real(rp) :: x,y,z
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb 
    integer  ::  i,j,k,p,s
    integer  ::  ilow,ihigh,jlow,jhigh,klow,khigh
    real(rp) :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    real(rp) :: coorxc,cooryc,coorzc
    real(rp) :: coorx_cent,coory_cent,coorz_cent
    real(rp) :: radin,radin2,dist2c
    real(rp) :: dx,dy,dz
    integer  :: nb,nbsend,nbrecv
    integer  :: nrrequests
    integer  :: arrayrequests(1:3)
    integer  :: arraystatuses(MPI_STATUS_SIZE,1:3)
    real(rp) :: sgndist(1:8)
    integer  :: idp,tag
    real(rp) :: dxp2_cent,dyp2_cent,dzp2_cent
    real(rp) :: alpha_eulc
    real(rp) :: deltasc
    real(rp) :: normalc(1:3)
    !
    dx       = 1./dli(1)
    dy       = 1./dli(2)
    dz       = 1./dli(3)
    radin    = (radius-retrac)-eps_sol*dx
    radin2   = radin**2
    alphac(:,:,:)=0._rp
    !
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
                  dxp2_cent = (coorx_cent - anb(nbrecv,p)%x)**2
                  !
                  coory_cent = boundfrontmyid + (j-0.5_rp)*dy
                  dyp2_cent = (coory_cent - anb(nbrecv,p)%y)**2
                  !
                  coorz_cent = (k-0.5_rp)*dz
                  dzp2_cent = (coorz_cent - anb(nbrecv,p)%z)**2
                  !
                  dist2c = dxp2_cent+dyp2_cent+dzp2_cent
                  !
                  normalc = [coorx_cent - anb(nbrecv,p)%x,  &
                             coory_cent - anb(nbrecv,p)%y,  &
                             coorz_cent - anb(nbrecv,p)%z    ]
                  normalc = normalc / sqrt(dist2c)
                  !
                  norm_partx(i,j,k)=normalc(1)
                  norm_party(i,j,k)=normalc(2)
                  norm_partz(i,j,k)=normalc(3)
                  deltasc = sqrt(dist2c)-(radius-retrac)
                  call digitiser(deltasc,normalc,alpha_eulc)
                  if (dist2c<radin2) then
                    alphac(i,j,k) = 1._rp
                  else
                    alphac(i,j,k)=alphac(i,j,k)+alpha_eulc
                  end if
                  if (alphac(i,j,k)>1._rp) then
                    alphac(i,j,k)=1._rp
                  endif
                  !
                enddo ! do i=
                !
              enddo ! do j=
              !
            enddo ! do k=
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
    return
  end subroutine initeul
  !
#endif
#endif
end module prt_mod_initeul
