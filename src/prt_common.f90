module prt_mod_common
#if defined(_PARTICLE)
  use mpi
  use mod_types
!  use mod_param       , only: pi,visc,dims,l,dli,is_forced,velf,ng
  use mod_param       , only: pi,dims,l,mu12,rho12,dli,ng,nh_wide
  use mod_common_mpi  , only: right,rightfront,front,leftfront,left,leftback,back,rightback, &
                              prt_comm_cart,mpi_comm_world,myid,ierr, &
                              boundleftmyid,boundfrontmyid,status
  use prt_mod_param   , only: nqmax,np,radius
  !
  implicit none
  public
  !
!  real ,dimension(0:i1,0:j1,0:k1) :: unew,vnew,wnew,pnew, &
!       dudtold,dvdtold,dwdtold
#if !defined(_EULER)
  real(rp), dimension(:,:,:), allocatable :: uf,vf,wf
#endif
  real(rp), target, dimension(:,:,:), allocatable :: alphac
  real(rp), target, dimension(:,:,:), allocatable :: norm_partx,norm_party,norm_partz
  real(rp), dimension(:,:,:), allocatable :: uphase,vphase,wphase
!  real, dimension(0:i1,0:j1,0:k1) :: dudt,dvdt,dwdt 
!  real(rp), dimension(:,:,:), allocatable :: forcex,forcey,forcez
!  real(mytype) :: time,dt
!  real(mytype) ::  wi(itot+15), wj(jtot+15)
!  real, dimension(imax,jmax) :: xyrt
!  real, dimension(kmax) :: a,b,c
!  real(rp) :: forcextot,forceytot,forceztot
!  real :: u_bulk,v_bulk,w_bulk
!  real wallshearold,wallshearnew
!  real :: dpdx_sumrk
!  integer :: rkiter
!  real :: rkparalpha
  !
  ! particles
  !
  ! for simplifying the communication between threads when there is a
  ! new master, the derived type 'particle' should be organized in this way:
  !
  !   (1st) real data that has to be communicated when there is a new master
  !         is defined contiguously
  !
  !   (2nd) integer data that has to be communicated when there is a new master
  !         is defined contiguously
  !
  !   (3rd) real data that does not have to be communicated
  !
  !   (4th) integer data that does not have to be communicated
  !
  type particle
     real(rp) :: x,y,z,theta,phi, &
                 u,v,w, &
                 omx,omy,omz,omtheta, &
                 intu,intv,intw, &
                 intomx,intomy,intomz, &
                 intrhox,intrhoy,intrhoz, &
                 colfx,colfy,colfz, &
                 coltx,colty,coltz ! 24
     real(rp), dimension(nqmax) :: dx,dy,dz, &
                                   dxt,dyt,dzt, & 
                                   dut,dvt,dwt, &
                                   psi ! 10*nqmax
     ! total ammount of reals to be communicated: 24+10*nqmax
     integer :: qmax ! 1
     integer, dimension(nqmax) :: firstc ! 1*nqmax
     ! total ammount of integers to be communicated: 1+nqmax
!     integer, dimension(nfriendsmax) :: friend ! p*nfriends
     integer :: nfriends
     real(rp) :: fxltot,fyltot,fzltot, &
                 torqxltot,torqyltot,torqzltot,torqtheta, &
                 fcapx,fcapy,fcapz
#if !defined(_EULER)
     real(rp), dimension(:), allocatable :: xfp,yfp,zfp, &
                                            ul,vl,wl, &
                                            dudtl,dvdtl,dwdtl, &
                                            fxl,fyl,fzl
#endif
     real(rp) :: vol,mominert,ratiorho,rho
     integer :: mslv
     !integer, dimension(8) :: nb !! dangerous accesses at nb(0) 
     integer, dimension(0:8) :: nb 
  end type particle
  !
  type(particle), dimension(:), allocatable :: ep ! 'an example particle' array
  ! this contains all the particle info
  !
  type particle_old
     real(rp) :: x,y,z,theta,phi, &
                 u,v,w, &
                 omx,omy,omz,omtheta, &
                 intu,intv,intw, &
                 intomx,intomy,intomz, &
                 colfx,colfy,colfz, &
                 coltx,colty,coltz, &
                 fxltot,fyltot,fzltot, &
                 torqxltot,torqyltot,torqzltot, &
                 fcapx,fcapy,fcapz ! 24
     real(rp), dimension(nqmax) :: dx,dy,dz, &
                                   dxt,dyt,dzt, &
                                   dut,dvt,dwt, &
                                   psi ! 10*nqmax
     integer :: qmax ! 1
     integer, dimension(nqmax) :: firstc ! 1*nqmax
     integer :: mslv
  end type particle_old
  !
  type(particle_old), dimension(:), allocatable :: op !old particle array (integration of N-E equations)
  type(particle_old), dimension(:), allocatable :: tp !send (transmitted) particle array (re-ordering of masters)
  !
  ! The structure of 'particle_sumrk' follows the same criterion as 'particle'
  ! for its organization
  !
  type particle_sumrk
     real(rp) :: fxltot,fyltot,fzltot, &
                 torqxltot,torqyltot,torqzltot,torqtheta, &
                 colfx,colfy,colfz, &
                 coltx,colty,coltz, &
                 dudt,dvdt,dwdt, &
                 domxdt,domydt,domzdt
     real(rp), dimension(:), allocatable :: fxl,fyl,fzl
  end type particle_sumrk
  !
  type(particle_sumrk), dimension(:), allocatable :: rkp ! this contains some particle
  ! data integrated over the rk3
  ! substeps
  type particle_interior
     real(rp) :: x,y,z
     real(rp), dimension(:), allocatable :: xfp,yfp,zfp,dvlagr
!     real, dimension(nl+nl2+nl3+nl4) :: xfp,yfp,zfp,dvlagr
  end type particle_interior
  !
  integer, dimension(1:2) :: coords
  integer, dimension(0:8) :: neighbor
  !
  real(rp) :: dVlagr,dVeul
  real(rp) :: radfp
  real(rp), dimension(:), allocatable :: thetarc,phirc ! angular position of the Lfps
  !
  integer :: pmax,npmstr

  integer, dimension(:), allocatable :: nla
  !
  integer :: npmax
  integer :: NL,NL2,NL3,NL4,NLtot
  real(rp) :: offset
  real(rp) :: lref,uref,tref
  real(rp) :: retrac
  real(rp) :: solidity
  real(rp) :: coeff_f,coeff_t
  !
  real(rp), dimension(:), allocatable :: fx_tot,fy_tot,fz_tot
  real(rp), dimension(:), allocatable :: tx_tot,ty_tot,tz_tot
  !
  contains
  !
  subroutine prt_InitMemo(n,lo)
    implicit none
    integer :: i
    integer, dimension(3), intent(in) :: n
    integer, dimension(1:3), intent(in) :: lo
    integer :: coordsleft(1:2) ,coordsright(1:2)
    integer :: coordsfront(1:2),coordsback(1:2)
    integer :: coordsneighbor(1:2)
    !
    npmax = nint(min(1._rp*np,max(1._rp,10._rp*np/(1._rp*dims(1)*dims(2)))))
#if !defined(_EULER)
    offset = ( sqrt(3._rp*(1.5_rp**2)) )/dli(1) + 0.01_rp/dli(1)
#else
    offset = (1._rp + 0.01_rp)/dli(1)
#endif
    !
!    lref = l(3)
!    uref = 1._rp   !bulk_v_sup, &
!    if(is_forced(1)) uref = velf(1)
!    tref = lref/uref
    !

    retrac = 0._rp/dli(1)
    !
#if !defined(_EULER)
    NL =nint((pi/3._rp)*(12._rp*(((radius-retrac)*dli(1))**2)+1._rp))
    !     nr lfp's of second shell
    NL2=nint((pi/3._rp)*(12._rp*(((radius-1.5_rp/dli(1))*dli(1))**2._rp) + 1._rp ))
    !     nr lfp's of third shell
    NL3=nint((pi/3._rp)*(12._rp*(((radius-2.5_rp/dli(1))*dli(1))**2._rp) + 1._rp ))
    !     nr lfp's of fourth shell
    NL4=nint((pi/3._rp)*(12._rp*(((radius-3.5_rp/dli(1))*dli(1))**2._rp) + 1._rp ))
    !     nr lfp's of 1st-4th shell
    NLtot=NL+NL2+NL3+NL4
#endif
    !
    solidity = (1._rp*np)*(4._rp/3._rp)*pi*(radius**3)/(l(1)*l(2)*l(3))
    !
    ! set parameters for the lubrication model
    !
    ! mu12(1) is a temporary solution to be used instead of visc
    coeff_f = 6._rp*pi*(mu12(1)/rho12(1))*radius
    coeff_t = 8._rp*pi*(mu12(1)/rho12(1))*radius**2._rp
    !
#if !defined(_EULER)
    allocate(uf(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide), &
             vf(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide), &
             wf(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide))
    !
    uf(:,:,:) = 0.0_rp
    vf(:,:,:) = 0.0_rp
    wf(:,:,:) = 0.0_rp
#endif
    allocate(alphac(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide))
    allocate(norm_partx(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide), &
             norm_party(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide), &
             norm_partz(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide))
    allocate(uphase(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide), &
             vphase(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide), &
             wphase(1-nh_wide:n(1)+nh_wide,1-nh_wide:n(2)+nh_wide,1-nh_wide:n(3)+nh_wide))

    !
    alphac(:,:,:) = 0.0_rp
    !
!    allocate(forcex(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2), &
!             forcey(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2), &
!             forcez(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2))
    !
    allocate(ep(1:npmax))
    do i = 1,npmax
#if !defined(_EULER)
      allocate(ep(i)%xfp(1:NL),   ep(i)%yfp(1:NL),   ep(i)%zfp(1:NL),   &
               ep(i)%ul(1:NL),    ep(i)%vl(1:NL),    ep(i)%wl(1:NL),    &
               ep(i)%dudtl(1:NL), ep(i)%dvdtl(1:NL), ep(i)%dwdtl(1:NL), &
               ep(i)%fxl(1:NL),   ep(i)%fyl(1:NL),   ep(i)%fzl(1:NL))
#endif
    enddo
    !
    !$acc enter data copyin(ep)
    do i = 1,npmax
      !$acc enter data &
      !$acc & copyin(ep(i)%xfp,ep(i)%yfp,ep(i)%zfp) &
      !$acc & copyin(ep(i)%ul,ep(i)%vl,ep(i)%wl) &
      !$acc & copyin(ep(i)%dudtl,ep(i)%dvdtl,ep(i)%dwdtl) &
      !$acc & copyin(ep(i)%fxl,ep(i)%fyl,ep(i)%fzl)
    enddo
    !
    allocate(op(1:npmax), &
             tp(1:npmax))
    !

    allocate(rkp(1:npmax))
#if !defined(_EULER)
    do i = 1,npmax
      allocate(rkp(i)%fxl(1:NL),   rkp(i)%fyl(1:NL),   rkp(i)%fzl(1:NL))
    enddo
    !
    allocate(thetarc(1:NL), phirc(1:NL))
    allocate(nla(1:npmax))
#endif
    !
    allocate(fx_tot(1:npmax),fy_tot(1:npmax),fz_tot(1:npmax))
    allocate(tx_tot(1:npmax),ty_tot(1:npmax),tz_tot(1:npmax))
    !
    ! Check for particle dimensions
    if ( (radius+offset) > l(1)/dims(1)) then
      if (myid == 0) then
        write(6,*) 'Radius spheres larger than x-dimension of processes'
        write(6,*) 'Change the value of dims(1) in "input.nml".'
        write(6,*) 'Program aborted...'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    
    !
    if ( (radius+offset) > l(2)/dims(2)) then
      if (myid == 0) then
        write(6,*) 'Radius spheres larger than y-dimension of processes'
        write(6,*) 'Change the value of dims(2) in "input.nml".'
        write(6,*) 'Program aborted...'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    !
    !
    ! Define a cartesian communicator for particles and other useful MPI stuff
    !
    coords(1) = (lo(1)-1)*dims(1)/ng(1)
    coords(2) = (lo(2)-1)*dims(2)/ng(2)
    boundleftmyid  = (lo(1)-1)/dli(1) ! left  boundary
    boundfrontmyid = (lo(2)-1)/dli(2) ! front boundary
    !
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims(1:2),(/.true., .true./),.false.,prt_comm_cart,ierr)
    !
    call MPI_CART_SHIFT(prt_comm_cart,0,1,left,right,ierr)
    call MPI_CART_SHIFT(prt_comm_cart,1,1,front,back,ierr)
    !
    call MPI_CART_COORDS(prt_comm_cart,right,2,coordsright,ierr)
    call MPI_CART_COORDS(prt_comm_cart,front,2,coordsfront,ierr)
    coordsneighbor(1) = coordsright(1)
    coordsneighbor(2) = coordsfront(2)
    call MPI_CART_RANK(prt_comm_cart,coordsneighbor,rightfront,ierr)
    !
    call MPI_CART_COORDS(prt_comm_cart,back,2,coordsback,ierr)
    coordsneighbor(1) = coordsright(1)
    coordsneighbor(2) = coordsback(2)
    call MPI_CART_RANK(prt_comm_cart,coordsneighbor,rightback,ierr)
    !
    call MPI_CART_COORDS(prt_comm_cart,left,2,coordsleft,ierr)
    call MPI_CART_COORDS(prt_comm_cart,front,2,coordsfront,ierr)
    coordsneighbor(1) = coordsleft(1)
    coordsneighbor(2) = coordsfront(2)
    call MPI_CART_RANK(prt_comm_cart,coordsneighbor,leftfront,ierr)
    !
    call MPI_CART_COORDS(prt_comm_cart,back,2,coordsback,ierr)
    coordsneighbor(1) = coordsleft(1)
    coordsneighbor(2) = coordsback(2)
    call MPI_CART_RANK(prt_comm_cart,coordsneighbor,leftback,ierr)
    !
    neighbor(0) = myid
    neighbor(1) = right
    neighbor(2) = rightfront
    neighbor(3) = front
    neighbor(4) = leftfront
    neighbor(5) = left
    neighbor(6) = leftback
    neighbor(7) = back
    neighbor(8) = rightback
    !
    return
  end subroutine prt_InitMemo
  !
#endif
end module prt_mod_common
