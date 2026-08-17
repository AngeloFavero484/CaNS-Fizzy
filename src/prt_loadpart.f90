module prt_mod_loadpart
#if defined(_PARTICLE)
  use mpi
  use mod_types
  use mod_param        , only: dims,l,datadir
  use mod_common_mpi   , only: ierr,myid,prt_comm_cart
  use prt_mod_common   , only: ep,tp,neighbor,pmax,npmax,npmstr,offset
  use prt_mod_param    , only: nqmax,radius,np,send_real,send_int
  !
  implicit none
  private
  public loadpart
  !
  contains
  !
  subroutine loadpart(io)
    implicit none
    character(len=1), intent(in) :: io
    type particle_restart
       real(rp) :: x,y,z,theta,phi, &
            u,v,w, &
            omx,omy,omz,omtheta, &
            intu,intv,intw, &
            intomx,intomy,intomz, &
            colfx,colfy,colfz, &
            coltx,colty,coltz ! 24
       real(rp), dimension(nqmax) :: dx,dy,dz, &
            dxt,dyt,dzt, &
            dut,dvt,dwt, &
            psi ! 10*nqmax
       ! total ammount of reals to be communicated: 24+10*nqmax
       real(rp) :: qmax,idp ! 2
       real(rp) , dimension(nqmax) :: firstc ! 1*nqmax
    end type particle_restart
    type(particle_restart), allocatable, dimension(:) :: glob
    integer, parameter :: skipr = 24+10*nqmax, &
                          skipi = 2+nqmax, & 
                          skip = skipr+skipi
    integer,dimension(0:dims(1)*dims(2)-1) :: npmstr_glob,npmstr_glob_all
    ! skipi is is the same as in the common flie + the global id of the particle
    integer i,p,idp   !,j
    integer :: fh
    integer :: proccoords(2),procrank
    real(rp) :: leftbound,rightbound,frontbound,backbound
    real(rp) :: xdum
    integer :: lenr
    real(rp) :: dist    !,angle
    !real(rp) :: xp,yp
    real(rp) :: ax
    real(rp) :: ay
    integer :: count_mstr_all
!    character(len=7) :: istepchar
    !character(len=4) rankpr
    !character(len=11) :: tempstr2
    integer :: count_slve_loc   !,counter
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    !logical :: found_mstr
    integer, dimension(send_int) :: itemp
    integer :: mydisp
    integer, dimension(np) :: rkmstr_glob,rkmstr_fake,locid_fake,locid_true
    integer, dimension(np) :: rkmstr_glob_all,rkmstr_fake_all,locid_fake_all,locid_true_all
    real(rp), dimension(np) :: xcglob,ycglob,zcglob
    real(rp), dimension(np) :: xcglob_all,ycglob_all,zcglob_all
    integer :: nrrequests
    integer :: arrayrequests(1:3)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    integer :: tag,nb,nbsend,nbrecv
    integer :: premain
    !
    !inquire (iolength=lenr) xdum
    !lenr = sizeof(xdum)
    lenr=storage_size(xdum)/8
    !
    select case(io)
    case('r')
      pmax    = np/product(dims)
      premain = np - pmax*product(dims)
      if( myid < premain ) pmax = pmax + 1
      allocate(glob(pmax))
      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(datadir)//'prt.bin', &
                         MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      npmstr_glob(:) = 0
      npmstr_glob(myid) = pmax
      call MPI_ALLREDUCE(npmstr_glob(0),npmstr_glob_all(0),product(dims),MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
      mydisp = 0
      if(myid /= 0) mydisp = sum(npmstr_glob_all(0:myid-1))
      disp = mydisp*skip*lenr
      call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', & 
           MPI_INFO_NULL, ierr)
      if (pmax > 0) then
        call MPI_FILE_READ(fh,glob(1)%x,skip*pmax,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      endif
      call MPI_FILE_CLOSE(fh,ierr)
      rkmstr_fake(:) = 0
      locid_fake(:) = 0
      locid_true(:) = 0
      rkmstr_glob(:) = 0
      xcglob(:) = 0.0_rp
      ycglob(:) = 0.0_rp
      zcglob(:) = 0.0_rp
      do p=1,pmax
        idp = nint(glob(p)%idp)
        rkmstr_fake(idp) = myid
        locid_fake(idp) = p
        if (glob(p)%x < 0._rp .or. glob(p)%x > l(1) .or. &
            glob(p)%y < 0._rp .or. glob(p)%y > l(2) .or. &
            glob(p)%z < 0._rp .or. glob(p)%z > l(3)) then
          if (myid == 0) then
            write(6,*) 'Fatal error in initialisation of particle positions - '
            write(6,*) 'particle outside the domain!'
            write(6,*) 'Program aborted...'
          endif
          call mpi_finalize(ierr)
          stop
        endif
        ax = 0.5_rp
        ay = 0.5_rp
        if (glob(p)%x == l(1))  ax = 0.51_rp
        if (glob(p)%x == 0._rp) ax = 0.49_rp
        if (glob(p)%y == l(2))  ay = 0.51_rp
        if (glob(p)%y == 0._rp) ay = 0.49_rp
        proccoords(1) = nint(glob(p)%x*dims(1)/l(1) - ax)
        proccoords(2) = nint(glob(p)%y*dims(2)/l(2) - ay)
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        rkmstr_glob(idp) = procrank
        xcglob(idp) = glob(p)%x
        ycglob(idp) = glob(p)%y
        zcglob(idp) = glob(p)%z
      enddo
      call MPI_ALLREDUCE(rkmstr_glob(1),rkmstr_glob_all(1),np,MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
      call MPI_ALLREDUCE(rkmstr_fake(1),rkmstr_fake_all(1),np,MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
      call MPI_ALLREDUCE(locid_fake(1) ,locid_fake_all(1) ,np,MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
      call MPI_ALLREDUCE(xcglob(1) ,xcglob_all(1) ,np,MPI_REAL_RP,MPI_SUM,prt_comm_cart,ierr)
      call MPI_ALLREDUCE(ycglob(1) ,ycglob_all(1) ,np,MPI_REAL_RP,MPI_SUM,prt_comm_cart,ierr)
      call MPI_ALLREDUCE(zcglob(1) ,zcglob_all(1) ,np,MPI_REAL_RP,MPI_SUM,prt_comm_cart,ierr)
      rkmstr_glob(:) = rkmstr_glob_all(:)
      rkmstr_fake(:) = rkmstr_fake_all(:)
      locid_fake(:)  = locid_fake_all(:)
      xcglob(:) = xcglob_all(:)
      ycglob(:) = ycglob_all(:)
      zcglob(:) = zcglob_all(:)
      i = 0
      do idp = 1,np
        nrrequests = 0
        if(rkmstr_glob(idp) == myid) then
          i = i + 1
          nrrequests = nrrequests + 1
          CALL MPI_IRECV(tp(i)%x,send_real,MPI_REAL_RP,rkmstr_fake(idp),idp, &
                         prt_comm_cart,arrayrequests((nrrequests-1)*2+1),ierr)
          CALL MPI_IRECV(tp(i)%qmax,send_int,MPI_INTEGER,rkmstr_fake(idp),idp+np, &
                         prt_comm_cart,arrayrequests((nrrequests-1)*2+2),ierr)
          tp(i)%mslv = idp
          locid_true(idp) = i
        endif
        if(rkmstr_fake(idp) == myid) then
          p = locid_fake(idp)
          nrrequests = nrrequests + 1
          CALL MPI_ISEND(glob(p)%x,send_real,MPI_REAL_RP,rkmstr_glob(idp),idp,prt_comm_cart,arrayrequests((nrrequests-1)*2+1),ierr)
          itemp(1) = nint(glob(p)%qmax)
          itemp(2:nqmax+1) =  nint(glob(p)%firstc(:))
          CALL MPI_ISEND(itemp(1),send_int,MPI_INTEGER,rkmstr_glob(idp),idp+np,prt_comm_cart,arrayrequests((nrrequests-1)*2+2),ierr)
        endif
        nrrequests = nrrequests*2
        call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
      enddo
      call MPI_ALLREDUCE(locid_true(1) ,locid_true_all(1) ,np,MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
      locid_true(:)  = locid_true_all(:)
      pmax = i
      ! 
      ! Determine master and slave processes for each particle.
      !
      ! initialisation
      !
      ep(1:npmax)%mslv = 0
      ep(1:npmax)%x = 0._rp
      ep(1:npmax)%y = 0._rp
      ep(1:npmax)%z = 0._rp
      ep(1:npmax)%theta = 0._rp
      ep(1:npmax)%phi = 0._rp
      ep(1:npmax)%u = 0._rp
      ep(1:npmax)%v = 0._rp
      ep(1:npmax)%w = 0._rp
      ep(1:npmax)%omx = 0._rp
      ep(1:npmax)%omy = 0._rp
      ep(1:npmax)%omz = 0._rp
      ep(1:npmax)%omtheta = 0._rp
      ep(1:npmax)%intu = 0._rp
      ep(1:npmax)%intv = 0._rp
      ep(1:npmax)%intw = 0._rp
      ep(1:npmax)%intomx = 0._rp
      ep(1:npmax)%intomy = 0._rp
      ep(1:npmax)%intomz = 0._rp
      ep(1:npmax)%colfx = 0._rp
      ep(1:npmax)%colfy = 0._rp
      ep(1:npmax)%colfz = 0._rp
      ep(1:npmax)%coltx = 0._rp
      ep(1:npmax)%colty = 0._rp
      ep(1:npmax)%coltz = 0._rp
      ep(1:npmax)%qmax = 0
      do p=1,npmax
        ep(p)%dx(1:nqmax) = 0._rp
        ep(p)%dy(1:nqmax) = 0._rp
        ep(p)%dz(1:nqmax) = 0._rp
        ep(p)%dxt(1:nqmax) = 0._rp
        ep(p)%dyt(1:nqmax) = 0._rp
        ep(p)%dzt(1:nqmax) = 0._rp
        ep(p)%firstc(1:nqmax) = 0
        ep(p)%nb(1:8) = 0
      enddo
      i = 0
      npmstr = 0
      do idp=1,np
        ax = 0.5_rp
        ay = 0.5_rp
        if (xcglob(idp) == l(1))  ax = 0.51_rp
        if (xcglob(idp) == 0._rp) ax = 0.49_rp
        if (ycglob(idp) == l(2))  ay = 0.51_rp
        if (ycglob(idp) == 0._rp) ay = 0.49_rp
        proccoords(1) = nint(xcglob(idp)*dims(1)/l(1) - ax)
        proccoords(2) = nint(ycglob(idp)*dims(2)/l(2) - ay)
        leftbound     = (proccoords(1)  )*l(1)/(1._rp*dims(1)) ! left  boundary of particle's master
        rightbound    = (proccoords(1)+1)*l(1)/(1._rp*dims(1)) ! right boundary of particle's master
        frontbound    = (proccoords(2)  )*l(2)/(1._rp*dims(2)) ! front boundary of particle's master
        backbound     = (proccoords(2)+1)*l(2)/(1._rp*dims(2)) ! back  boundary of particle's master
        if(rkmstr_glob(idp) == myid) then
          npmstr = npmstr + 1
          i = i + 1
          p = locid_true(idp)
          ep(i)%mslv      = tp(p)%mslv
          ep(i)%x         = tp(p)%x
          ep(i)%y         = tp(p)%y
          ep(i)%z         = tp(p)%z
          ep(i)%theta     = tp(p)%theta
          ep(i)%phi       = tp(p)%phi
          ep(i)%u         = tp(p)%u
          ep(i)%v         = tp(p)%v
          ep(i)%w         = tp(p)%w
          ep(i)%omx       = tp(p)%omx
          ep(i)%omy       = tp(p)%omy
          ep(i)%omz       = tp(p)%omz
          ep(i)%omtheta   = tp(p)%omtheta
          ep(i)%intu      = tp(p)%intu
          ep(i)%intv      = tp(p)%intv
          ep(i)%intw      = tp(p)%intw
          ep(i)%intomx    = tp(p)%intomx
          ep(i)%intomy    = tp(p)%intomy
          ep(i)%intomz    = tp(p)%intomz
          ep(i)%colfx     = tp(p)%colfx
          ep(i)%colfy     = tp(p)%colfy
          ep(i)%colfz     = tp(p)%colfz
          ep(i)%coltx     = tp(p)%coltx
          ep(i)%colty     = tp(p)%colty
          ep(i)%coltz     = tp(p)%coltz
          ep(i)%dx(1:nqmax) = tp(p)%dx(1:nqmax)
          ep(i)%dy(1:nqmax) = tp(p)%dy(1:nqmax)
          ep(i)%dz(1:nqmax) = tp(p)%dz(1:nqmax)
          ep(i)%dxt(1:nqmax) = tp(p)%dxt(1:nqmax)
          ep(i)%dyt(1:nqmax) = tp(p)%dyt(1:nqmax) 
          ep(i)%dzt(1:nqmax) = tp(p)%dzt(1:nqmax)
          ep(i)%qmax = tp(p)%qmax
          ep(i)%firstc(1:nqmax) = tp(p)%firstc(1:nqmax)
          ! neighbor 1
          if ( ep(i)%x > (rightbound-(radius+offset)) ) then 
            ep(i)%nb(1) = 1 ! neighbor 1 is slave of particle ep(p)%mslv 
          endif
          ! neighbor 2
          dist = sqrt( (rightbound-ep(i)%x)**2._rp + (frontbound-ep(i)%y)**2._rp )
          if ( abs(dist) < (radius+offset) ) then
            ep(i)%nb(2) = 1 ! neighbor 2 is slave of particle ep(p)%mslv
          endif
          ! neighbor 3
          if ( ep(i)%y < (frontbound+(radius+offset)) ) then
            ep(i)%nb(3) = 1 ! neighbor 3 is slave of particle ep(p)%mslv
          endif
          ! neighbor 4
          dist = sqrt( (leftbound-ep(i)%x)**2._rp + (frontbound-ep(i)%y)**2._rp ) 
          if ( abs(dist) < (radius+offset) ) then
            ep(i)%nb(4) = 1 ! neighbor 4 is slave of particle ep(p)%mslv
          endif
          ! neighbor 5
          if ( ep(i)%x < (leftbound+(radius+offset)) ) then
            ep(i)%nb(5) = 1 ! neighbor 5 is slave of particle ep(p)%mslv
          endif
          ! neighbor 6
          dist = sqrt( (leftbound-ep(i)%x)**2._rp + (backbound-ep(i)%y)**2._rp )
          if ( abs(dist) < (radius+offset) ) then
            ep(i)%nb(6) = 1 ! neighbor 6 is slave of particle ep(p)%mslv
          endif
          ! neighbor 7
          if ( ep(i)%y > (backbound-(radius+offset)) ) then
            ep(i)%nb(7) = 1 ! neighbor 7 is slave of particle ep(p)%mslv
          endif
          ! neighbor 8
          dist = sqrt( (rightbound-ep(i)%x)**2._rp + (backbound-ep(i)%y)**2._rp )
          if ( abs(dist) < (radius+offset) ) then
            ep(i)%nb(8) = 1 ! neighbor 8 is slave of particle ep(p)%mslv
          endif
        else
          count_slve_loc = 0
          ! neighbor 1 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(idp)/l(1) - ax ) + 1
          proccoords(2) = nint( dims(2)*ycglob(idp)/l(2) - ay ) 
          call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
          if (myid == procrank) then
            if ( xcglob(idp) > (rightbound-(radius+offset)) ) then
              if(count_slve_loc == 0 ) i = i+1
              count_slve_loc = count_slve_loc + 1
              ep(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
              ep(i)%nb(5) = 1      ! neighbor 5 of myid is particle's master
            endif
          endif
          ! neighbor 2 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(idp)/l(1) - ax ) + 1
          proccoords(2) = nint( dims(2)*ycglob(idp)/l(2) - ay ) - 1
          call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
          if (myid == procrank) then
            dist = sqrt( (rightbound-xcglob(idp))**2._rp + (frontbound-ycglob(idp))**2._rp ) 
            if ( abs(dist) < (radius+offset) ) then
              if(count_slve_loc == 0 ) i = i+1
              count_slve_loc = count_slve_loc + 1
              ep(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
              ep(i)%nb(6) = 1      ! neighbor 6 of myid is particle's master
            endif
          endif
          ! neighbor 3 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(idp)/l(1) - ax )
          proccoords(2) = nint( dims(2)*ycglob(idp)/l(2) - ay ) - 1
          call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
          if (myid == procrank) then
            if ( ycglob(idp) < (frontbound+(radius+offset)) ) then
              if(count_slve_loc == 0 ) i = i+1
              count_slve_loc = count_slve_loc + 1
              ep(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
              ep(i)%nb(7) = 1      ! neighbor 7 of myid is particle's master
            endif
          endif
          ! neighbor 4 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(idp)/l(1) - ax ) - 1
          proccoords(2) = nint( dims(2)*ycglob(idp)/l(2) - ay ) - 1
          call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
          if (myid == procrank) then
            dist = sqrt( (leftbound-xcglob(idp))**2._rp + (frontbound-ycglob(idp))**2._rp )
            if ( abs(dist) < (radius+offset) ) then
              if(count_slve_loc == 0 ) i = i+1
              count_slve_loc = count_slve_loc + 1
              ep(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
              ep(i)%nb(8) = 1      ! neighbor 8 of myid is particle's master
            endif
          endif
          ! neighbor 5 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(idp)/l(1) - ax ) - 1
          proccoords(2) = nint( dims(2)*ycglob(idp)/l(2) - ay )
          call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
          if (myid == procrank) then
            if ( xcglob(idp) < (leftbound+(radius+offset)) ) then
              if(count_slve_loc == 0 ) i = i+1
              count_slve_loc = count_slve_loc + 1
              ep(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
              ep(i)%nb(1) = 1      ! neighbor 1 of myid is particle's master
            endif
          endif
          ! neighbor 6 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(idp)/l(1) - ax ) - 1
          proccoords(2) = nint( dims(2)*ycglob(idp)/l(2) - ay ) + 1
          call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
          if (myid == procrank) then
            dist = sqrt( (leftbound-xcglob(idp))**2._rp + (backbound-ycglob(idp))**2._rp )
            if ( abs(dist) < (radius+offset) ) then
              if(count_slve_loc ==  0 ) i = i+1
              count_slve_loc = count_slve_loc + 1
              ep(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
              ep(i)%nb(2) = 1      ! neighbor 2 of myid is particle's master
            endif
          endif
          ! neighbor 7 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(idp)/l(1) - ax )
          proccoords(2) = nint( dims(2)*ycglob(idp)/l(2) - ay ) + 1
          call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
          if (myid == procrank) then
            if ( ycglob(idp) > (backbound-(radius+offset)) ) then
              if(count_slve_loc == 0 ) i = i+1
              count_slve_loc = count_slve_loc + 1
              ep(i)%mslv = -idp     ! myid is slave of particle p=ap(p)%mslv
              ep(i)%nb(3) = 1      ! neighbor 3 of myid is particle's master
            endif
          endif
          ! neighbor 8 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(idp)/l(1) - ax ) + 1
          proccoords(2) = nint( dims(2)*ycglob(idp)/l(2) - ay ) + 1
          call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
          if (myid == procrank) then
            dist = sqrt( (rightbound-xcglob(idp))**2._rp + (backbound-ycglob(idp))**2._rp )
            if ( abs(dist) < (radius+offset) ) then
              if(count_slve_loc == 0 ) i = i+1
              count_slve_loc = count_slve_loc + 1
              ep(i)%mslv = -idp     ! myid is slave of particle p=ap(p)%mslv
              ep(i)%nb(4) = 1      ! neighbor 4 of myid is particle's master
            endif
          endif
        endif
      enddo
      pmax = i 
      !
      write(6,'(A7,I5,A8,I7,A18,I7,A11,A8,I5)') 'Thread ', myid, ' masters ', npmstr, ' and is slave for ', &
                                                pmax-npmstr, ' particles. ', ' pmax = ', pmax
      !
      call MPI_ALLREDUCE(npmstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
      if (count_mstr_all /= np) then
        write(6,*) 'Fatal error in initialisation of particle positions!'
        write(6,*) 'Program aborted...'
        call mpi_abort(prt_comm_cart,ierr,ierr)
        stop
      elseif(pmax > npmax) then
        write(6,*) 'Size of local particle array of process ', myid, ' is too small!'
        write(6,*) 'Program aborted... (later I will simply write a warning and allocate more memory)'
        call mpi_abort(prt_comm_cart,ierr,ierr)
        stop
      else
        write(6,*) 'The particles were successfully initialized in thread ', myid, ' !'
      endif
!      print*,myid,'pre do p 1-pmax',pmax
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
          if (ep(p)%mslv > 0) then
            ! myid is master of particle ap(p)%mslv
            if (ep(p)%nb(nbsend) == 1) then
              ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
              if ( neighbor(nbsend) /= myid ) then
                nrrequests = nrrequests + 1
                call MPI_ISEND(ep(p)%x,11,MPI_REAL_RP,neighbor(nbsend), &
                              tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
                ! send x,y,z,theta,phi,u,v,w,omx,omy,omz
              endif
            endif
          endif
          if (ep(p)%mslv < 0) then
            ! myid is slave of particle -ap(p)%mslv
            if (ep(p)%nb(nbrecv) == 1) then
              ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
              nrrequests = nrrequests + 1
!              print*,'pre irecv'
              call MPI_IRECV(ep(p)%x,11,MPI_REAL_RP,neighbor(nbrecv), &
                             tag,prt_comm_cart,arrayrequests((nrrequests-1) + 1),ierr)
!              print*,'post irecv'
              ! recv x,y,z,theta,phi,u,v,w,omx,omy,omz
            endif
          endif
        enddo ! do nb=
!        print*,myid,'pre waitall'
        call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
!        print*,myid,'post waitall'
      enddo
      !
      ! initial particle positions written to file
      !
      !  write(rankpr,'(i4.4)') myid
      !  open(25,file=datadir//'mslv'//rankpr//'.txt')
      !  do p=1,pmax
      !    if (ap(p)%mslv .gt. 0) then
      !      counter = 0
      !      do i=1,8
      !        if (ap(p)%nb(i) .eq. 1) then
      !          counter = counter+1
      !          write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
      !                myid,' ',p,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
      !          write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
      !                myid,' ',p,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
      !        endif
      !      enddo
      !      ! in case of no overlap with any neighbor
      !      if (counter .eq. 0) then
      !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
      !                myid,' ',p,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
      !        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
      !                myid,' ',p,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
      !      endif
      !    endif
      !    if (ap(p)%mslv .lt. 0) then
      !      do i=1,8
      !        if (ap(p)%nb(i) .eq. 1) then
      !          write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
      !                myid,' ',p,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
      !          write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
      !                myid,' ',p,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
      !        endif
      !      enddo
      !    endif
      !  enddo
      !  close(25)
      !
      ! write tecplot output file with particle distribution
      !
      !  if (myid.eq.0) then
      !    open(25,file=datadir//'distr_particles.plt',form='formatted')
      !    write(25,*) 'TITLE="Tecplot Output"'
      !    write(25,*) 'VARIABLES= "X" "Y"'
      !    do p=1,pmax
      !      write(tempstr2,'(A,I2.2)') 'Particle ',p
      !      write(25,*) 'ZONE F=POINT T="', tempstr2, '" I=',250,' J=1'
      !      do i=1,250
      !        angle = 2.*pi*i/250.
      !        xp = glob(p)%x+radius*cos(angle)
      !        yp = glob(p)%y+radius*sin(angle)
      !          if (xp.lt.0) then
      !            xp = xp + lx
      !          elseif(xp.gt.lx) then
      !            xp = xp - lx
      !          endif
      !          if (yp.lt.0) then
      !            yp = yp + ly
      !          elseif(yp.gt.ly) then
      !            yp = yp - ly
      !          endif 
      !       write(25,*) xp/lx,yp/ly
      !      enddo
      !    enddo
      !    close(25) 
      !  endif
      !
      deallocate(glob)
    case('w')
      !
      ! write particle related data directly in parallel with MPI-IO
      !
      allocate(glob(npmstr))
      npmstr_glob(:) = 0
      npmstr_glob(myid) = npmstr
      call MPI_ALLREDUCE(npmstr_glob(0),npmstr_glob_all(0),product(dims),MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
      mydisp = 0
      if(myid /= 0) mydisp = sum(npmstr_glob_all(0:myid-1))
      i = 0
      do p=1,pmax
        if(ep(p)%mslv > 0) then
          idp = ep(p)%mslv
          i = i + 1
          glob(i)%x = ep(p)%x !WRITE THIS IN A MORE COMPACT WAY!
          glob(i)%y = ep(p)%y
          glob(i)%z = ep(p)%z
          glob(i)%theta = ep(p)%theta
          glob(i)%phi = ep(p)%phi
          glob(i)%u = ep(p)%u
          glob(i)%v = ep(p)%v
          glob(i)%w = ep(p)%w
          glob(i)%omx = ep(p)%omx
          glob(i)%omy = ep(p)%omy
          glob(i)%omz = ep(p)%omz
          glob(i)%omtheta = ep(p)%omtheta
          glob(i)%intu = ep(p)%intu
          glob(i)%intv = ep(p)%intv
          glob(i)%intw = ep(p)%intw
          glob(i)%intomx = ep(p)%intomx
          glob(i)%intomy = ep(p)%intomy
          glob(i)%intomz = ep(p)%intomz
          glob(i)%colfx = ep(p)%colfx
          glob(i)%colfy = ep(p)%colfy
          glob(i)%colfz = ep(p)%colfz
          glob(i)%coltx = ep(p)%coltx
          glob(i)%colty = ep(p)%colty
          glob(i)%coltz = ep(p)%coltz
          glob(i)%dx(1:nqmax) = ep(p)%dx(1:nqmax)
          glob(i)%dy(1:nqmax) = ep(p)%dy(1:nqmax)
          glob(i)%dz(1:nqmax) = ep(p)%dz(1:nqmax)
          glob(i)%dxt(1:nqmax) = ep(p)%dxt(1:nqmax)
          glob(i)%dyt(1:nqmax) = ep(p)%dyt(1:nqmax)
          glob(i)%dzt(1:nqmax) = ep(p)%dzt(1:nqmax)
          glob(i)%firstc(1:nqmax) = 1._rp*ep(p)%firstc(1:nqmax)
          glob(i)%qmax = 1._rp*ep(p)%qmax
          glob(i)%idp = 1._rp*idp
        endif
      enddo
      call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(datadir)//'prt.bin', &
                         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = mydisp*skip*lenr
      call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', & 
                             MPI_INFO_NULL, ierr)
      if (npmstr > 0) then
        call MPI_FILE_WRITE(fh,glob(1)%x,skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      endif
      call MPI_FILE_CLOSE(fh,ierr)
      deallocate(glob)
    end select
    !
    return
  end subroutine loadpart
  !
#endif
end module prt_mod_loadpart
