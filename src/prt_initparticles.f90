module prt_mod_initparticles
#if defined(_PARTICLE)
  use mpi
  use decomp_2d
  use mod_types
  use mod_param       , only: datadir,ng,dims,l
  use mod_common_mpi  , only: ierr,myid,status,prt_comm_cart
  use prt_mod_param   , only: mominert,np,radius,volp,ratiorho,rho_s
  use prt_mod_common  , only: ep,npmax,npmstr,pmax,offset,coords,neighbor
  !
  implicit none
  private
  public initparticles
  !
  contains
  !
  subroutine initparticles
    implicit none
    real(rp), dimension(np) :: xcglob,ycglob,zcglob,thetacglob,phicglob
    integer :: i,j,p,pp,rk  !,k
    integer :: proccoords(2),procrank  !proccoords(1:ndims)
    !integer, dimension(2) :: sbuf,rbuf
    real(rp) :: leftbound,rightbound,frontbound,backbound
    real(rp) :: dist,distx,disty,distz,distzw  !,angle
    !real(rp) :: xp,yp
    real(rp) :: ax   !,sx
    real(rp) :: ay   !,sy
    real(rp) :: rn
    integer :: counter,crys
    !integer :: idp
    !character(len=5) rankpr
    !character(len=7) :: tempstr
    !character(len=11) :: tempstr2
    integer :: count_mstr,count_mstr_all,count_slve_loc
    integer, allocatable, dimension(:) :: seed
    !integer, intent(in), dimension(3) :: lo
    !
    ! position of spheres: global initialisation by root (myid = 0)
    !
    allocate(seed(64))
    seed(:) = 18
    crys=0
    if (myid == 0) then
      call random_seed( put = seed )
      write(6,*) 'Particles position as crystal = 1 or random =  0. Selection: ',crys
      open(23,file=trim(datadir)//"position_spheres.txt")
      !
      if(crys == 1) then
          !    p=0
          !    ! pseudo-crystal
          !    do k=1,4
          !      do j=1,16
          !        do i=1,8
          !          p=p+1
          !          thetacglob(p) = 0.
          !          phicglob(p)   = 0.
          !          call random_number(rn)
          !          xcglob(p)     = (l(1)/4.0_rp)*((1.0_rp*i)-0.5_rp) + 0.25_rp*(rn-0.5_rp)*radius
          !          call random_number(rn)
          !          ycglob(p)     = (l(2)/8.0_rp)*((1.0_rp*j)-0.5_rp) + 0.25_rp*(rn-0.5)*radius
          !          call random_number(rn)
          !          zcglob(p)     = (l(3)/2.0_rp)*((1.0_rp*k)-0.5_rp) + (rn-0.5_rp)*radius
          !          write(6,*) 'Location of sphere #',p
          !          write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
          !          write(23,'(I3,5E32.16)') p,thetacglob(p),phicglob(p), &
          !                                  xcglob(p),ycglob(p),zcglob(p)
          !        enddo
          !      enddo
          !    enddo
      else
        ! pseudo-random
        counter = 0
        !
        do p=1,np
          !
          if(np == 1) then
            thetacglob(p) = 0._rp
            phicglob(p)   = 0._rp
            xcglob(p)     = l(1)*0.5_rp
            ycglob(p)     = l(2)*0.5_rp
            zcglob(p)     = l(3)*0.755_rp
            PRINT *, "z", zcglob(p)
          else      
111         continue
            thetacglob(p) = 0.0_rp
            phicglob(p)   = 0.0_rp
            call random_number(rn)
            xcglob(p)     = l(1)*rn
            call random_number(rn)
            ycglob(p)     = l(2)*rn
            call random_number(rn)
            zcglob(p)     = l(3)*rn
            distzw = min(abs(l(3)-zcglob(p)),abs(zcglob(p)))
            if(distzw.lt.1.05_rp*radius) goto 111
            !
            do pp=1,p
              if (pp == p) goto 444 ! could be changed by changing the loop limits!
              distz = abs(zcglob(p)- zcglob(pp))
              !
              do j=-1,1
                disty = abs(ycglob(p)- (ycglob(pp)+j*l(2)))
                if(disty > 2.05_rp*radius) goto 222
                !
                do i=-1,1
                  distx = abs(xcglob(p)- (xcglob(pp)+i*l(1)))
                  if(distx > 2.05_rp*radius) goto 333
                  if(distx > 2.05_rp*radius.or. &
                     disty > 2.05_rp*radius.or. &
                     distz > 2.05_rp*radius.or. &
                     p == pp) then
                     ! good particle
                  else
                     dist = distx**2+disty**2+distz**2
                     !
                     if((dist < (4.2_rp*radius**2))) then
                       !write(*,*)'RANDOM DEVIATION'
                       !write(*,*)dist,distw
                       counter=counter+1
                       goto 111
                     endif
                     !
                  endif
333               continue
                  !
                enddo
222             continue
                !
              enddo
444           continue
              !
            enddo
            !
          endif
          !
          write(6,*) 'Location of sphere #',p
          write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
          write(23,'(I6,5E32.16)') p,thetacglob(p),phicglob(p), &
                                   xcglob(p),ycglob(p),zcglob(p)
        enddo
        !
        write(*,*)'RANDOM DEVIATIONS: ',counter
        p=p-1
        !    p = np
      endif
      !
      close(23)
      if ( (p /= np) ) then
        print*,counter,np
        write(6,*) 'Fatal error in initialisation of particle positions!'
        write(6,*) 'Program aborted...'
        call mpi_finalize(ierr)
        stop
      endif
      !
      do rk=1,Nproc-1
        call MPI_SSEND(xcglob    ,np,MPI_REAL_RP,rk,rk+0*(Nproc-1),MPI_COMM_WORLD,ierr)
        call MPI_SSEND(ycglob    ,np,MPI_REAL_RP,rk,rk+1*(Nproc-1),MPI_COMM_WORLD,ierr)
        call MPI_SSEND(zcglob    ,np,MPI_REAL_RP,rk,rk+2*(Nproc-1),MPI_COMM_WORLD,ierr)
        call MPI_SSEND(thetacglob,np,MPI_REAL_RP,rk,rk+3*(Nproc-1),MPI_COMM_WORLD,ierr)
        call MPI_SSEND(phicglob  ,np,MPI_REAL_RP,rk,rk+4*(Nproc-1),MPI_COMM_WORLD,ierr)
      enddo
      !
    else ! if myid is not 0:
      call MPI_RECV(xcglob    ,np,MPI_REAL_RP,0,myid+0*(Nproc-1),MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(ycglob    ,np,MPI_REAL_RP,0,myid+1*(Nproc-1),MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(zcglob    ,np,MPI_REAL_RP,0,myid+2*(Nproc-1),MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(thetacglob,np,MPI_REAL_RP,0,myid+3*(Nproc-1),MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(phicglob  ,np,MPI_REAL_RP,0,myid+4*(Nproc-1),MPI_COMM_WORLD,status,ierr)
    endif
    !
    ! Determine master and slave processes for each particle.
    !
    !
    ! initialisation
    !
    ep(1:npmax)%x = 0.0_rp
    ep(1:npmax)%y = 0.0_rp
    ep(1:npmax)%z = 0.0_rp
    ep(1:npmax)%theta = 0.0_rp
    ep(1:npmax)%phi = 0.0_rp
    ep(1:npmax)%mslv = 0.0_rp
    do i=1,npmax
      !
      ep(i)%nb(1:8) = 0
      !
    enddo
    !
    count_mstr = 0
    i = 0
    !
    ax = 0.5_rp
    ay = 0.5_rp
    !
    pmax = 0
    do p=1,np
      !
      if (xcglob(p) < 0.0_rp .or. xcglob(p) > l(1) .or. &
          ycglob(p) < 0.0_rp .or. ycglob(p) > l(2) .or. &
          zcglob(p) < radius .or. zcglob(p) > l(3) - radius) then
        !
        if (myid == 0) then
          write(6,*) 'Fatal error in initialisation of particle positions - '
          write(6,*) 'particle outside the domain!'
          write(6,*) 'Program aborted...'
        endif
        !
        call mpi_finalize(ierr)
        stop
        !
      endif
      !
      if (xcglob(p) == l(1)  ) ax = 0.51_rp
      if (xcglob(p) == 0.0_rp) ax = 0.49_rp
      if (ycglob(p) == l(2)  ) ay = 0.51_rp
      if (ycglob(p) == 0.0_rp) ay = 0.49_rp
      !
      proccoords(1) = nint(xcglob(p)*dims(1)/l(1) - ax)
      proccoords(2) = nint(ycglob(p)*dims(2)/l(2) - ay)
      leftbound     = (proccoords(1)  )*l(1)/(1.0_rp*dims(1)) ! left  boundary of particle's master
      rightbound    = (proccoords(1)+1)*l(1)/(1.0_rp*dims(1)) ! right boundary of particle's master
      frontbound    = (proccoords(2)  )*l(2)/(1.0_rp*dims(2)) ! front boundary of particle's master
      backbound     = (proccoords(2)+1)*l(2)/(1.0_rp*dims(2)) ! back  boundary of particle's master
      call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)     
      !
      if (myid == procrank) then
        i = i + 1
        count_mstr = count_mstr + 1
        ep(i)%x = xcglob(p)
        ep(i)%y = ycglob(p)
        ep(i)%z = zcglob(p)
        ep(i)%theta = thetacglob(p)
        ep(i)%phi = phicglob(p)
        ep(i)%mslv = p
        !neighbor 1
        if   ( ep(i)%x > (rightbound-(radius+offset))) then 
           ep(i)%nb(1) = 1 !neighbor 1 is slave of particle ep(i)%mslv 
        endif
        !neighbor 2
        dist = sqrt( (rightbound-ep(i)%x)**2 + (frontbound-ep(i)%y)**2 )
        if ( abs(dist) < (radius+offset) ) then
           ep(i)%nb(2) = 1 !neighbor 2 is slave of particle ep(i)%mslv
        endif
        !neighbor 3
        if ( ep(i)%y < (frontbound+(radius+offset))) then
           ep(i)%nb(3) = 1 !neighbor 3 is slave of particle ep(i)%mslv
        endif
        !neighbor 4
        dist = sqrt( (leftbound-ep(i)%x)**2 + (frontbound-ep(i)%y)**2 ) 
        if ( abs(dist) < (radius+offset)) then
           ep(i)%nb(4) = 1 !neighbor 4 is slave of particle ep(i)%mslv
        endif
        !neighbor 5
        if ( ep(i)%x < (leftbound+(radius+offset)) ) then
           ep(i)%nb(5) = 1 !neighbor 5 is slave of particle ep(i)%mslv
        endif
        !neighbor 6
        dist = sqrt( (leftbound-ep(i)%x)**2 + (backbound-ep(i)%y)**2 )
        if ( abs(dist) < (radius+offset) ) then
           ep(i)%nb(6) = 1 !neighbor 6 is slave of particle ep(i)%mslv
        endif
        !neighbor 7
        if ( ep(i)%y > (backbound-(radius+offset)) ) then
           ep(i)%nb(7) = 1 !neighbor 7 is slave of particle ep(i)%mslv
        endif
        !neighbor 8
        dist = sqrt( (rightbound-ep(i)%x)**2 + (backbound-ep(i)%y)**2 )
        if ( abs(dist) < (radius+offset) ) then
           ep(i)%nb(8) = 1 !neighbor 8 is slave of particle ep(i)%mslv
        endif
      else
        count_slve_loc = 0
        !neighbor 1 of particle's master
        proccoords(1) = nint( dims(1)*xcglob(p)/l(1) - ax ) + 1
        proccoords(2) = nint( dims(2)*ycglob(p)/l(2) - ay ) 
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (myid == procrank) then
          if ( xcglob(p) > (rightbound-(radius+offset))) then
            if(count_slve_loc == 0) i = i+1
            ep(i)%mslv = -p      !myid is slave of particle abs(ep(i)%mslv)
            ep(i)%nb(5) = 1      !neighbor 5 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            ep(i)%x = xcglob(p)
            ep(i)%y = ycglob(p)
            ep(i)%z = zcglob(p)
            ep(i)%theta = thetacglob(p)
            ep(i)%phi = phicglob(p)
          endif
        endif
        !neighbor 2 of particle's master
        proccoords(1) = nint( dims(1)*xcglob(p)/l(1) - ax ) + 1
        proccoords(2) = nint( dims(2)*ycglob(p)/l(2) - ay ) - 1
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (myid == procrank) then
          dist = sqrt( (rightbound-xcglob(p))**2 + (frontbound-ycglob(p))**2 ) 
          if ( abs(dist) < (radius+offset) ) then
            if(count_slve_loc == 0) i = i+1
            ep(i)%mslv = -p      !myid is slave of particle abs(ep(i)%mslv)
            ep(i)%nb(6) = 1      !neighbor 6 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            ep(i)%x = xcglob(p)
            ep(i)%y = ycglob(p)
            ep(i)%z = zcglob(p)
            ep(i)%theta = thetacglob(p)
            ep(i)%phi = phicglob(p)
          endif
        endif
        !neighbor 3 of particle's master
        proccoords(1) = nint( dims(1)*xcglob(p)/l(1) - ax )
        proccoords(2) = nint( dims(2)*ycglob(p)/l(2) - ay ) - 1
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (myid == procrank) then
          if ( ycglob(p) < (frontbound+(radius+offset)) ) then
            if(count_slve_loc == 0) i = i+1
            ep(i)%mslv = -p      !myid is slave of particle abs(ep(i)%mslv)
            ep(i)%nb(7) = 1      !neighbor 7 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            ep(i)%x = xcglob(p)
            ep(i)%y = ycglob(p)
            ep(i)%z = zcglob(p)
            ep(i)%theta = thetacglob(p)
            ep(i)%phi = phicglob(p)
          endif
        endif
        !neighbor 4 of particle's master
        proccoords(1) = nint( dims(1)*xcglob(p)/l(1) - ax ) - 1
        proccoords(2) = nint( dims(2)*ycglob(p)/l(2) - ay ) - 1
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (myid == procrank) then
           dist = sqrt( (leftbound-xcglob(p))**2 + (frontbound-ycglob(p))**2 )
          if ( abs(dist) < (radius+offset) ) then
            if(count_slve_loc == 0) i = i+1
            ep(i)%mslv = -p      !myid is slave of particle abs(ep(i)%mslv)
            ep(i)%nb(8) = 1      !neighbor 8 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            ep(i)%x = xcglob(p)
            ep(i)%y = ycglob(p)
            ep(i)%z = zcglob(p)
            ep(i)%theta = thetacglob(p)
            ep(i)%phi = phicglob(p)
          endif
        endif
        !neighbor 5 of particle's master
        proccoords(1) = nint( dims(1)*xcglob(p)/l(1) - ax ) - 1
        proccoords(2) = nint( dims(2)*ycglob(p)/l(2) - ay )
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (myid == procrank) then
          if ( xcglob(p) < (leftbound+(radius+offset)) ) then
            if(count_slve_loc == 0) i = i+1
            ep(i)%mslv = -p      !myid is slave of particle abs(ep(i)%mslv)
            ep(i)%nb(1) = 1      !neighbor 1 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            ep(i)%x = xcglob(p)
            ep(i)%y = ycglob(p)
            ep(i)%z = zcglob(p)
            ep(i)%theta = thetacglob(p)
            ep(i)%phi = phicglob(p)
          endif
        endif
        !neighbor 6 of particle's master
        proccoords(1) = nint( dims(1)*xcglob(p)/l(1) - ax ) - 1
        proccoords(2) = nint( dims(2)*ycglob(p)/l(2) - ay ) + 1
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (myid == procrank) then
           dist = sqrt( (leftbound-xcglob(p))**2 + (backbound-ycglob(p))**2 )
          if ( abs(dist) < (radius+offset) ) then
            if(count_slve_loc == 0) i = i+1
            ep(i)%mslv = -p      !myid is slave of particle abs(ep(i)%mslv)
            ep(i)%nb(2) = 1      !neighbor 2 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            ep(i)%x = xcglob(p)
            ep(i)%y = ycglob(p)
            ep(i)%z = zcglob(p)
            ep(i)%theta = thetacglob(p)
            ep(i)%phi = phicglob(p)
          endif
        endif
        !neighbor 7 of particle's master
        proccoords(1) = nint( dims(1)*xcglob(p)/l(1) - ax )
        proccoords(2) = nint( dims(2)*ycglob(p)/l(2) - ay ) + 1
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (myid == procrank) then
          if ( ycglob(p) > (backbound-(radius+offset)) ) then
            if(count_slve_loc == 0) i = i+1
            ep(i)%mslv = -p      !myid is slave of particle p=ep(i)%mslv
            ep(i)%nb(3) = 1      !neighbor 3 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            ep(i)%x = xcglob(p)
            ep(i)%y = ycglob(p)
            ep(i)%z = zcglob(p)
            ep(i)%theta = thetacglob(p)
            ep(i)%phi = phicglob(p)
          endif
        endif
        !neighbor 8 of particle's master
        proccoords(1) = nint( dims(1)*xcglob(p)/l(1) - ax ) + 1
        proccoords(2) = nint( dims(2)*ycglob(p)/l(2) - ay ) + 1
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (myid == procrank) then
           dist = sqrt( (rightbound-xcglob(p))**2 + (backbound-ycglob(p))**2 )
          if ( abs(dist) < (radius+offset) ) then
            if(count_slve_loc == 0) i = i+1
            ep(i)%mslv = -p      !myid is slave of particle p=ep(i)%mslv
            ep(i)%nb(4) = 1      !neighbor 4 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            ep(i)%x = xcglob(p)
            ep(i)%y = ycglob(p)
            ep(i)%z = zcglob(p)
            ep(i)%theta = thetacglob(p)
            ep(i)%phi = phicglob(p)
          endif
        endif
      endif
    enddo
    !
    ! maximum number of particles in a thread is equal to the number of particles
    ! 'mastered' and 'slaved' by it
    !
    pmax = i 
    npmstr = count_mstr
    write(6,'(A7,I5,A8,I5,A18,I5,A11,A8,I5)') 'Thread ', myid, ' masters ', count_mstr, &
                                              ' and is slave for ', pmax-npmstr, ' particles. ', ' pmax = ', pmax
    !
    call MPI_ALLREDUCE(count_mstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
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
    !
    ! initial particle positions written to file
    !
    !write(rankpr,'(i5.5)') myid
    !open(25,file=datadir//'mslv'//rankpr//'.txt')
    !do p=1,pmax
    !  idp = abs(ap(p)%mslv)
    !  if (ap(p)%mslv .gt. 0) then
    !    counter = 0
    !    do i=1,8
    !      if (ap(p)%nb(i) .eq. 1) then
    !        counter = counter+1
    !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !      endif
    !    enddo
    !    ! in case of no overlap with any neighbor
    !    if (counter .eq. 0) then
    !      write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
    !      write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
    !    endif
    !  endif
    !  if (ap(p)%mslv .lt. 0) then
    !    do i=1,8
    !      if (ap(p)%nb(i) .eq. 1) then
    !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !      endif
    !    enddo
    !  endif
    !enddo
    !close(25)
    !!
    !! write topology of the domain decomposition in a tecplot file
    !!
    !if (myid==0) then
    !  open(25,file=datadir//'decomp.plt',form='formatted')
    !  write(25,*) 'TITLE="Tecplot Output"'
    !  write(25,*) 'VARIABLES= "X" "Y" "VAR"'
    !  write(25,*) 'ZONE F=POINT T="Rank 00" I=',imax,' J=',jmax
    !  do i=1,imax
    !    do j=1,jmax
    !      write(25,*) (i+zstart(1)-1)*dx/lx,(j+zstart(2)-1)*dy/ly,myid
    !    enddo
    !  enddo
    !  do rk=1,nproc-1
    !    CALL MPI_RECV(rbuf,2,MPI_INTEGER,rk,rk,MPI_COMM_WORLD,status,error)
    !    write(tempstr,'(A,I2.2)') 'Rank ',rk
    !  write(25,*) 'TITLE="Tecplot Output"'
    !  write(25,*) 'VARIABLES= "X" "Y" "VAR"'
    !  write(25,*) 'ZONE F=POINT T="', tempstr, '" I=',imax,' J=',jmax
    !    do i=1,imax
    !      do j=1,jmax
    !        write(25,*) (i+rbuf(1)-1)*dx/lx,(j+rbuf(2)-1)*dy/ly,rk
    !      enddo
    !    enddo
    !  enddo
    !  close(25)
    !else
    !  sbuf(1) = zstart(1)
    !  sbuf(2) = zstart(2)
    !  CALL MPI_SEND(sbuf,2,MPI_INTEGER,0,myid,MPI_COMM_WORLD,error)
    !endif
    !
    ! write initial particle distribution in a tecplot file
    !
    !if (myid.eq.0) then
    !  open(25,file=datadir//'distr_particles.plt',form='formatted')
    !  write(25,*) 'TITLE="Tecplot Output"'
    !  write(25,*) 'VARIABLES= "X" "Y"'
    !  do p=1,np
    !    write(tempstr2,'(A,I2.2)') 'Particle ',p
    !    write(25,*) 'ZONE F=POINT T="', tempstr2, '" I=',250,' J=1'
    !    do i=1,250
    !      angle = 2.*pi*i/250.
    !      xp = xcglob(p)+radius*cos(angle)
    !      yp = ycglob(p)+radius*sin(angle)
    !        if (xp.lt.0) then
    !          xp = xp + lx
    !        elseif(xp.gt.lx) then
    !          xp = xp - lx
    !        endif
    !        if (yp.lt.0) then
    !          yp = yp + ly
    !        elseif(yp.gt.ly) then
    !          yp = yp - ly
    !        endif 
    !      write(25,*) xp/lx,yp/ly
    !    enddo
    !  enddo
    !  close(25) 
    !endif
    !
    ep(1:npmax)%u = 0.0_rp
    ep(1:npmax)%v = 0.0_rp
    ep(1:npmax)%w = -28.78_rp
    ep(1:npmax)%omx = 0.0_rp
    ep(1:npmax)%omy = 0.0_rp
    ep(1:npmax)%omz = 0.0_rp
    ep(1:npmax)%omtheta = 0.0_rp
    ep(1:npmax)%vol = volp
    ep(1:npmax)%mominert = mominert
    ep(1:npmax)%ratiorho = ratiorho
    ep(1:npmax)%rho = rho_s
    ep(1:npmax)%intu = 0.0_rp
    ep(1:npmax)%intv = 0.0_rp
    ep(1:npmax)%intw = 0.0_rp
    ep(1:npmax)%intomx = 0.0_rp
    ep(1:npmax)%intomy = 0.0_rp
    ep(1:npmax)%intomz = 0.0_rp
    ep(1:npmax)%colfx = 0.0_rp
    ep(1:npmax)%colfy = 0.0_rp
    ep(1:npmax)%colfz = 0.0_rp
    ep(1:npmax)%coltx = 0.0_rp
    ep(1:npmax)%colty = 0.0_rp
    ep(1:npmax)%coltz = 0.0_rp
    ep(1:npmax)%qmax = 0
    do p=1,npmax
      ep(p)%dx(:) = 0.0_rp
      ep(p)%dy(:) = 0.0_rp
      ep(p)%dz(:) = 0.0_rp
      ep(p)%dxt(:) = 0.0_rp
      ep(p)%dyt(:) = 0.0_rp
      ep(p)%dzt(:) = 0.0_rp
      ep(p)%firstc(:) = 0
#if !defined(_EULER)
      ep(p)%xfp(:) = 0.0_rp
      ep(p)%yfp(:) = 0.0_rp
      ep(p)%zfp(:) = 0.0_rp
#endif
    enddo
    !
    return
  end subroutine initparticles
  !
#endif
end module prt_mod_initparticles
