module prt_mod_forcing
#if defined(_PARTICLE)
#if !defined(_EULER)
  use mpi
  use mod_types
  use mod_common_mpi  , only: prt_comm_cart,myid,ierr
  use prt_mod_common  , only: ep,npmax,dVlagr,pmax,nla,neighbor, &
                              fx_tot,fy_tot,fz_tot, &
                              tx_tot,ty_tot,tz_tot
  !
  implicit none
  private
  public complagrforces,updtlagrforces  !,updtintermediatevel
  !
  contains
  !
  subroutine complagrforces(dti,rkit)
    implicit none
    real(rp), intent(in) :: dti
    integer, intent(in) :: rkit
    integer :: ll,p
    type pneighbor
       real(rp) :: intu,intv,intw, &
                   intomx,intomy,intomz
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb
    integer :: nb,nbsend,nbrecv
    integer :: nrrequests
    integer :: arrayrequests(1:3) ! 3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    integer :: idp,tag
    real(rp) :: auxu,auxv,auxw,auxomx,auxomy,auxomz
    !real(rp) :: isperiodx,isperiody
    !real(rp) :: coorxfp,cooryfp,coorzfp
    !logical :: isout
    !real(rp) ::buf_send(6),buf_recv(6)
    !
    !$omp parallel default(none)         &
    !$omp shared(ep,pmax,nla,dvlagr,dti) &
    !$omp private(p,ll)
    !$omp do
    do p=1,pmax
      if(ep(p)%mslv /= 0) then
        !ap(p)%fxltot_old  = ap(p)%fxltot ***ask wp why this was here before***
        !ap(p)%fyltot_old  = ap(p)%fyltot
        !ap(p)%fzltot_old  = ap(p)%fzltot
        !ap(p)%torqxltotold  = ap(p)%torqxltot
        !ap(p)%torqyltotold  = ap(p)%torqyltot
        !ap(p)%torqzltotold  = ap(p)%torqzltot
        !ap(p)%torqthetaold = ap(p)%torqtheta
        ep(p)%fxltot  = 0._rp
        ep(p)%fyltot  = 0._rp
        ep(p)%fzltot  = 0._rp
        ep(p)%torqxltot  = 0._rp
        ep(p)%torqyltot  = 0._rp
        ep(p)%torqzltot  = 0._rp
        !
        do ll=1,nla(p)
          ep(p)%fxl(ll) = (ep(p)%ul(ll)-ep(p)%dudtl(ll))*dti
          ep(p)%fyl(ll) = (ep(p)%vl(ll)-ep(p)%dvdtl(ll))*dti
          ep(p)%fzl(ll) = (ep(p)%wl(ll)-ep(p)%dwdtl(ll))*dti
          ep(p)%fxltot = ep(p)%fxltot + ep(p)%fxl(ll)
          ep(p)%fyltot = ep(p)%fyltot + ep(p)%fyl(ll)
          ep(p)%fzltot = ep(p)%fzltot + ep(p)%fzl(ll)
          ep(p)%torqxltot = ep(p)%torqxltot + (ep(p)%yfp(ll)-ep(p)%y)*ep(p)%fzl(ll) - &
                                              (ep(p)%zfp(ll)-ep(p)%z)*ep(p)%fyl(ll)
          ep(p)%torqyltot = ep(p)%torqyltot + (ep(p)%zfp(ll)-ep(p)%z)*ep(p)%fxl(ll) - &
                                              (ep(p)%xfp(ll)-ep(p)%x)*ep(p)%fzl(ll)
          ep(p)%torqzltot = ep(p)%torqzltot + (ep(p)%xfp(ll)-ep(p)%x)*ep(p)%fyl(ll) - &
                                              (ep(p)%yfp(ll)-ep(p)%y)*ep(p)%fxl(ll)
!             endif
        enddo
        ep(p)%fxltot = ep(p)%fxltot*dVlagr
        ep(p)%fyltot = ep(p)%fyltot*dVlagr
        ep(p)%fzltot = ep(p)%fzltot*dVlagr
        ep(p)%torqxltot = ep(p)%torqxltot*dVlagr
        ep(p)%torqyltot = ep(p)%torqyltot*dVlagr
        ep(p)%torqzltot = ep(p)%torqzltot*dVlagr
        ! Torque working on angle theta:
        ! phi = 0      --> torqtheta =  torqyltot
        ! phi = pi/2   --> torqtheta = -torqxltot
        ! phi = pi     --> torqtheta = -torqyltot
        ! phi = 3*pi/2 --> torqtheta =  torqxltot
        ep(p)%torqtheta = (ep(p)%torqyltot*cos(ep(p)%phi)) - &
                          (ep(p)%torqxltot*sin(ep(p)%phi))
        ! phi is value of phic at time step n (weak coupling)
      endif
    enddo
    !$omp end parallel
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
       !       call MPI_IRECV(buf_recv,6,MPI_REAL_RP,neighbor(nbsend), &
       !                      tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
       !       anb(nbrecv,p)%intu=buf_recv(1)
       !       anb(nbrecv,p)%intv=buf_recv(2)
       !       anb(nbrecv,p)%intw=buf_recv(3)
       !       anb(nbrecv,p)%intomx=buf_recv(4)
       !       anb(nbrecv,p)%intomy=buf_recv(5)
       !       anb(nbrecv,p)%intomz=buf_recv(6)
              ! recv intu,intv,intw (...) -> 6 contiguous info
              ! (see definition of type pneighbor in the begining of the subroutine)
            endif
          endif
        endif
        if (ep(p)%mslv < 0) then
          ! myid is slave of particle -ep(p)%mslv
          if (ep(p)%nb(nbrecv) == 1) then
            ! neighbor(nbrecv) is rank of master of particle -ep(p)%mslv
            nrrequests = nrrequests + 1
    !        buf_send=[anb(0,p)%intu,anb(0,p)%intv,anb(0,p)%intw,anb(0,p)%intomx,anb(0,p)%intomy,anb(0,p)%intomz]
    !        call MPI_ISEND(buf_send,6,MPI_REAL_RP,neighbor(nbrecv), &
    !                       tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
            call MPI_ISEND(anb(0,p)%intu,6,MPI_REAL_RP,neighbor(nbrecv), &
                           tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
            ! send intu,intv,intw (...) -> 6 contiguous info
            ! (see definition of type pneighbor in the begining of the subroutine)
          endif
        endif
      enddo ! do nb=
      call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
    enddo
    !
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
        if (rkit == 1) then
          fx_tot(p) = 0.0_rp
          fy_tot(p) = 0.0_rp
          fz_tot(p) = 0.0_rp
          tx_tot(p) = 0.0_rp
          ty_tot(p) = 0.0_rp
          tz_tot(p) = 0.0_rp
        endif
      endif
    enddo
    !$omp end parallel
    !
    return
  end subroutine complagrforces
  !
  subroutine updtlagrforces(dti,dtitot,ibmit)  ! ,maxerror_all,maxavererror_all)
    implicit none
    real(rp), intent(in) :: dti
    real(rp), intent(in) :: dtitot
    integer, intent(in) :: ibmit
!    real(rp), intent(out) :: maxavererror_all,maxerror_all
    integer :: ll,p
!    real(rp) :: maxavererror,maxerror
    real(rp) :: avererror,averlagvel  !,lagvel,relerror,averlagvel_all
    !real(rp) :: errdistr
    type pneighbor
       real(rp) :: intu,intv,intw, &
                   intomx,intomy,intomz
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb
    integer :: nb,nbsend,nbrecv  !,i
    integer :: nrrequests
    integer :: arrayrequests(1:3) ! 3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    integer :: idp,tag
    real(rp) :: auxu,auxv,auxw,auxomx,auxomy,auxomz
    real(rp) :: dt
    !real(rp) :: isperiodx,isperiody
    !real(rp) :: coorxfp,cooryfp,coorzfp
    !logical :: isout
    !real(rp) ::buf_send(6),buf_recv(6)
    !
    !maxerror = 0._rp
    !maxavererror = 0._rp
    dt=1.0_rp/dti
    !
    !$omp parallel default(none)         &
    !$omp shared(ep,pmax,nla,dvlagr,dti) &
    !$omp private(p,ll,lagvel,relerror,avererror,averlagvel)
!    !$omp reduction(max:maxerror,maxavererror)
    !$omp do
    do p=1,pmax
      if(ep(p)%mslv.ne.0) then
        ep(p)%fxltot = 0._rp
        ep(p)%fyltot = 0._rp
        ep(p)%fzltot = 0._rp
        ep(p)%torqxltot = 0._rp
        ep(p)%torqyltot = 0._rp
        ep(p)%torqzltot = 0._rp
        avererror  = 0._rp
        averlagvel = 0._rp
        !
        do ll=1,nla(p)
!          isperiodx = 0.
!          isperiody = 0.
!          if (ap(p)%xfp(l).lt.0.+0.5*dx) isperiodx =  1.
!          if (ap(p)%xfp(l).ge.lx+0.5*dx) isperiodx = -1.
!          if (ap(p)%yfp(l).lt.0.+0.5*dy) isperiody =  1.
!          if (ap(p)%yfp(l).ge.ly+0.5*dy) isperiody = -1.
!          isout = .false.
!          coorxfp = (ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi
!          if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
!          cooryfp = (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi
!          if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
!          if (.not.isout) then
!             coorzfp =  ap(p)%zfp(l)*dzi
                !      if (abs(ap(p)%mslv) .eq. 1) then
                !        lagvel        = sqrt(ap(p)%ul(l)**2 + ap(p)%vl(l)**2 + ap(p)%wl(l)**2)
                !        errdistr      = sqrt(ap(p)%dudtl(l)**2 + ap(p)%dvdtl(l)**2 + ap(p)%dwdtl(l)**2) - &
                !                        lagvel
                !        avererror     = avererror + errdistr
                !        averlagvel    = averlagvel + lagvel
                !        relerror      = 100.*errdistr/(1.e-12 + lagvel)
                !        !if ( maxerror .lt. abs(relerror)) maxerror = abs(relerror)
               !        maxerror=max(maxerror,abs(relerror))
                !      endif
          ep(p)%fxl(ll) = ep(p)%fxl(ll) + (ep(p)%ul(ll)-ep(p)%dudtl(ll))*dti ! update force
          ep(p)%fyl(ll) = ep(p)%fyl(ll) + (ep(p)%vl(ll)-ep(p)%dvdtl(ll))*dti ! update force
          ep(p)%fzl(ll) = ep(p)%fzl(ll) + (ep(p)%wl(ll)-ep(p)%dwdtl(ll))*dti ! update force
          ep(p)%fxltot = ep(p)%fxltot + ep(p)%fxl(ll)
          ep(p)%fyltot = ep(p)%fyltot + ep(p)%fyl(ll)
          ep(p)%fzltot = ep(p)%fzltot + ep(p)%fzl(ll)
          ep(p)%torqxltot = ep(p)%torqxltot + (ep(p)%yfp(ll)-ep(p)%y)*ep(p)%fzl(ll) - &
                                              (ep(p)%zfp(ll)-ep(p)%z)*ep(p)%fyl(ll)
          ep(p)%torqyltot = ep(p)%torqyltot + (ep(p)%zfp(ll)-ep(p)%z)*ep(p)%fxl(ll) - &
                                              (ep(p)%xfp(ll)-ep(p)%x)*ep(p)%fzl(ll)
          ep(p)%torqzltot = ep(p)%torqzltot + (ep(p)%xfp(ll)-ep(p)%x)*ep(p)%fyl(ll) - &
                                              (ep(p)%yfp(ll)-ep(p)%y)*ep(p)%fxl(ll)
!             endif
        enddo
          !  if (abs(ap(p)%mslv) .eq. 1) then
          !    averlagvel = averlagvel/(1.*NL)
          !    avererror  = 100.*(avererror/(1.*NL))/(1.e-12+averlagvel) ! incorrect now!
          !  endif
          !!
          !  maxavererror=max(maxavererror,abs(avererror))
        ep(p)%fxltot = ep(p)%fxltot*dVlagr
        ep(p)%fyltot = ep(p)%fyltot*dVlagr
        ep(p)%fzltot = ep(p)%fzltot*dVlagr
        ep(p)%torqxltot = ep(p)%torqxltot*dVlagr
        ep(p)%torqyltot = ep(p)%torqyltot*dVlagr
        ep(p)%torqzltot = ep(p)%torqzltot*dVlagr
        ! Torque working on angle theta:
        !   phi = 0      --> torqtheta =  torqyltot
        !   phi = pi/2   --> torqtheta = -torqxltot
        !   phi = pi     --> torqtheta = -torqyltot
        !   phi = 3*pi/2 --> torqtheta =  torqxltot
        ep(p)%torqtheta = (ep(p)%torqyltot*cos(ep(p)%phi)) - &
                          (ep(p)%torqxltot*sin(ep(p)%phi))
        ! phi is value of phic at time step n for which torqyltot and torqxltot are computed!
      endif
    enddo
    !$omp end parallel 
    !
    !call mpi_allreduce(maxavererror,maxavererror_all,1,mpi_real8,mpi_max,comm_cart,error) ! incorrect now!
    !call mpi_allreduce(maxerror,maxerror_all,1,mpi_real8,mpi_max,comm_cart,error)
    !maxavererror_all = 0._rp
    !maxerror_all     = 0._rp
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
      !        call MPI_IRECV(buf_recv,6,MPI_REAL_RP,neighbor(nbsend), &
      !                       tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
      !        anb(nbrecv,p)%intu=buf_recv(1)
      !        anb(nbrecv,p)%intv=buf_recv(2)
      !        anb(nbrecv,p)%intw=buf_recv(3)
      !        anb(nbrecv,p)%intomx=buf_recv(4)
      !        anb(nbrecv,p)%intomy=buf_recv(5)
      !        anb(nbrecv,p)%intomz=buf_recv(6)
              ! recv intu,intv,intw (...) -> 6 contiguous info
              ! (see definition of type pneighbor in the begining of the subroutine)
            endif
          endif
        endif
        if (ep(p)%mslv < 0) then
          ! myid is slave of particle -ep(p)%mslv
          if (ep(p)%nb(nbrecv) == 1) then
            ! neighbor(nbrecv) is rank of master of particle -ep(p)%mslv
            nrrequests = nrrequests + 1
    !        buf_send=[anb(0,p)%intu,anb(0,p)%intv,anb(0,p)%intw,anb(0,p)%intomx,anb(0,p)%intomy,anb(0,p)%intomz]
    !        call MPI_ISEND(buf_send,6,MPI_REAL_RP,neighbor(nbrecv), &
    !                       tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
            call MPI_ISEND(anb(0,p)%intu,6,MPI_REAL_RP,neighbor(nbrecv), &
                           tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
            ! send intu,intv,intw (...) -> 6 contiguous info
            ! (see definition of type pneighbor in the begining of the subroutine)
          endif
        endif
      enddo ! do nb=
      call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
    enddo
    !
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
        if (ibmit == 2) then
          fx_tot(p) = fx_tot(p) + ep(p)%fxltot*dt*dtitot
          fy_tot(p) = fy_tot(p) + ep(p)%fyltot*dt*dtitot
          fz_tot(p) = fz_tot(p) + ep(p)%fzltot*dt*dtitot
          tx_tot(p) = tx_tot(p) + ep(p)%torqxltot*dt*dtitot
          ty_tot(p) = ty_tot(p) + ep(p)%torqyltot*dt*dtitot
          tz_tot(p) = tz_tot(p) + ep(p)%torqzltot*dt*dtitot
          PRINT *, "fx_tot(p)", fx_tot(p)
          PRINT *, "fy_tot(p)", fy_tot(p)
          PRINT *, "fz_tot(p)", fz_tot(p)
          PRINT *, "tx_tot(p)", tx_tot(p)
          PRINT *, "ty_tot(p)", ty_tot(p)
          PRINT *, "tz_tot(p)", tz_tot(p)
        endif
      endif
    enddo
    !$omp end parallel
    !
    return
  end subroutine updtlagrforces
  !

#endif
#endif
end module prt_mod_forcing
