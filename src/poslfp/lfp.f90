program init_force_points
  !
  use mod_param
  use mod_loadd
  !
  !     Program computes the positions of the Lagrangian force points wrt center of sphere.
  !
  implicit none
  integer l,q,j,i,k
  real(rp) :: rn
  real(rp) :: absforce
  real(rp) :: thetanew,phinew
  real(rp) :: forcephi(1:NL),forcetheta(1:NL)
  real(rp) :: forcethetaold(1:NL),forcephiold(1:NL)
  real(rp) :: difvectorx,difvectory,difvectorz
  real(rp) :: lengthdifvector
  real(rp) :: ethetax,ethetay,ethetaz
  real(rp) :: ephix,ephiy,ephiz
  real(rp) :: erx,ery,erz
  real(rp) :: forcephimean,forcephivar
  real(rp) :: forcethetamean,forcethetavar
  real(rp) :: xfp(1:NL),yfp(1:NL),zfp(1:NL)
  integer :: begin
  real(rp),parameter :: dt=1.d-4 !1.e-3
  real(rp) :: dist1max,dist2max,dist
  integer :: lmin1,lmin2
  real(rp) :: phi(1:NL),theta(1:NL)
  integer,parameter :: restart=0
  real(rp),parameter :: radfp=radius-retraction !IMPORTANT, lfp's positioned at slightly lower radius!
  integer,parameter :: qmax=50
  real(rp) :: distr(1:qmax),lowbound,upbound
  integer :: sumq
  real(rp) ::  sumsquares
  integer :: rkiter
  !
  !     Initialization     
  !
  write(6,*) 'NL = ',NL
  write(6,*) 'Volume occupied by 1 lfp = ',(4._rp/3._rp)*pi*((radfp+0.5_rp/dxi)**3._rp-(radfp-0.5_rp/dxi)**3._rp )/(1._rp*NL)
  write(6,*) 'Volume of Eulerian grid cell = ', 1._rp/dxi/dyi/dzi
  !
  do q=1,qmax
    distr(q) = 0._rp
  enddo
  !
  sumq = 0
  rn = 10._rp
  !
  do l=1,NL
    call random_number(rn)
    call random_number(rn)
    call random_number(rn)
    q=0
    lowbound = -1._rp/qmax
    upbound  = 0._rp
777 q=q+1
    lowbound = lowbound + 1._rp/qmax
    upbound  = upbound + 1._rp/qmax
    !
    if ( (rn > lowbound) .and. (rn < upbound) ) then
      distr(q) = distr(q)+1._rp
      sumq     = sumq + 1
    else
      go to 777
    endif
    !
    phi(l) = rn*2._rp*pi !random value between 0 and 2*pi         
    call random_number(rn)
    call random_number(rn)
    call random_number(rn)
    q=0
    lowbound = -1._rp/qmax
    upbound  = 0._rp
888 q=q+1
    lowbound = lowbound + 1._rp/qmax
    upbound  = upbound + 1._rp/qmax
    !
    if ( (rn > lowbound) .and. (rn < upbound) ) then
      distr(q) = distr(q)+1._rp
      sumq     = sumq + 1
    else
      go to 888
    endif
    !
    theta(l)    = rn*pi    !random value between 0 and pi         
    xfp(l)      = radfp*sin(theta(l))*cos(phi(l))
    yfp(l)      = radfp*sin(theta(l))*sin(phi(l))
    zfp(l)      = radfp*cos(theta(l))
    forcetheta(l)    = 0._rp
    forcephi(l)      = 0._rp
  enddo
  !
  if (sumq /= 2*NL) then
    write(6,*) 'Error! Failure in binning of random number.'
    write(6,*) 'Program aborted...'
    stop
  endif
  !
  !open(18,file=datadir//'distribution.txt')
  !do q=1,qmax
  !  write(18,'(3E16.8)') (1.*q-0.5)/(1.*qmax),distr(q),sumq/(1.*qmax)
  !enddo
  !close(18)
!     data to file
  open(42,file=datadir//'lagrangianforcepoints_init')
  write(42,*) 'VARIABLES = "lfpx","lfpy","lfpz","theta","phi","radfp"'
  write(42,*) 'ZONE T="Zone1"',' I=',NL,', F=POINT'
  write(42,*) ''
  !
  do l=1,NL
    write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp
  enddo
  !
  close(42) 
!     data for sphere
  !open(42,file=datadir//'datasphere')
  !write(42,*) 'VARIABLES = "x","y","z","r"'
  !write(42,*) 'ZONE T="Zone1"',' I=',imax,' J=',jmax,' K=',kmax,', F=POINT'
  !write(42,*) ''
  !do k=1,kmax
  !  do j=1,jmax
  !    do i=1,imax
  !      write(42,'(4E16.8)') (i-imax/2)/dxi,(j-jmax/2)/dyi,(k-kmax/2)/dzi, &
  ! (sqrt( ((i-imax/2)/dxi)**2. + ((j-jmax/2)/dyi)**2. + ((k-kmax/2)/dzi)**2. ))/radfp
  !    enddo
  !  enddo
  !enddo
  !close(42)
!
!     Compute 'Coulomb force' acting on every charged particle
!
  begin = 0
  if (restart == 1) then 
    call loadd(0,begin,phi,theta)
    !
    do l=1,NL
      xfp(l) = radfp*sin(theta(l))*cos(phi(l))
      yfp(l) = radfp*sin(theta(l))*sin(phi(l))
      zfp(l) = radfp*cos(theta(l))
    enddo
    !
  endif
  !
  write(6,*) 'begin = ',begin
  !open(22,file=datadir//'variances',status='replace')
  !
  do j=begin+1,50000!100000
    if (mod(j,1000) == 0) then
      write(6,*) 'Iteration step = ',j
    endif
    rkiter = 0
999 rkiter = rkiter+1
    !
    !$omp parallel default(shared) &
    !$omp private(l,q,ethetax,ethetay,ethetaz,ephix,ephiy,ephiz) &
    !$omp private(difvectorx,difvectory,difvectorz,sumsquares,absforce,lengthdifvector)
    !$omp do
    do l=1,NL
!     unit vector in theta direction
      ethetax    =  cos( theta(l) )*cos( phi(l) )
      ethetay    =  cos( theta(l) )*sin( phi(l) )
      ethetaz    = -sin( theta(l) )
!     unit vector in phi direction
      ephix      = -sin( phi(l) )
      ephiy      =  cos( phi(l) )
      ephiz      =  0._rp
!     unit vector in radial direction
!     erx        = xfp(l)/radfp
!     ery        = yfp(l)/radfp
!     erz        = zfp(l)/radfp
      forcetheta(l) = 0._rp
      forcephi(l)   = 0._rp
      !
      do q=1,NL
        if (q /= l) then
          difvectorx  = xfp(l)-xfp(q) 
          difvectory  = yfp(l)-yfp(q)
          difvectorz  = zfp(l)-zfp(q)
          sumsquares  = difvectorx*difvectorx + difvectory*difvectory + difvectorz*difvectorz
          absforce   = 1._rp/( sumsquares ) !1/(distance**2)
!         normalized difference vector
          lengthdifvector = sqrt( sumsquares )
          difvectorx = difvectorx/lengthdifvector 
          difvectory = difvectory/lengthdifvector
          difvectorz = difvectorz/lengthdifvector
!         force component in theta direction              
          forcetheta(l) = forcetheta(l) + (ethetax*difvectorx + ethetay*difvectory + ethetaz*difvectorz)*absforce
!         force component in phi direction              
          forcephi(l) = forcephi(l) + (ephix*difvectorx + ephiy*difvectory + ephiz*difvectorz)*absforce
        endif
      enddo
      forcetheta(l) = forcetheta(l)/(1._rp*NL-1._rp) !averaged force 
      forcephi(l)   = forcephi(l)/(1._rp*NL-1._rp)   !averaged force
!     The forces are expected to scale with dxi**2. Since NL is proportional to dxi**2, 
!     they scale approximately with NL. By dividing the forces by NL, the forces become independent
!     of the resolution. This then implies that the time step is insensitive to the resolution.
    enddo
    !$omp end do
    !$omp end parallel
!   RK3 scheme
    if (rkiter == 1) then
      !
      !$omp parallel default(shared) &
      !$omp private(l)
      !$omp do
      do l=1,NL
        theta(l)         = theta(l) + dt*(32._rp/60._rp)*forcetheta(l)
        phi(l)           = phi(l) + dt*(32._rp/60._rp)*forcephi(l)
        xfp(l)           = radfp*sin(theta(l))*cos(phi(l))
        yfp(l)           = radfp*sin(theta(l))*sin(phi(l))
        zfp(l)           = radfp*cos(theta(l))
        forcethetaold(l) = forcetheta(l)
        forcephiold(l)   = forcephi(l)
      enddo
      !$omp end do
      !$omp end parallel
      go to 999
    endif
    if (rkiter == 2) then
      !$omp parallel default(shared) &
      !$omp private(l)
      !$omp do
      do l=1,NL
        theta(l)         = theta(l) + dt*( (25._rp/60._rp)*forcetheta(l) - (17._rp/60._rp)*forcethetaold(l) )
        phi(l)           = phi(l) + dt*( (25._rp/60._rp)*forcephi(l)-(17._rp/60._rp)*forcephiold(l) )
        xfp(l)           = radfp*sin(theta(l))*cos(phi(l))
        yfp(l)           = radfp*sin(theta(l))*sin(phi(l))
        zfp(l)           = radfp*cos(theta(l))
        forcethetaold(l) = forcetheta(l)
        forcephiold(l)   = forcephi(l)
      enddo
      !$omp end do
      !$omp end parallel
      go to 999
    endif
    if (rkiter == 3) then
      !$omp parallel default(shared) &
      !$omp private(l)
      !$omp do
      do l=1,NL
        theta(l) = theta(l) + dt*( (45._rp/60._rp)*forcetheta(l) - (25._rp/60._rp)*forcethetaold(l) )
        phi(l)   = phi(l) + dt*( (45._rp/60._rp)*forcephi(l)-(25._rp/60._rp)*forcephiold(l) )
        xfp(l)   = radfp*sin(theta(l))*cos(phi(l))
        yfp(l)   = radfp*sin(theta(l))*sin(phi(l))
        zfp(l)   = radfp*cos(theta(l))
      enddo
      !$omp end do
      !$omp end parallel
    endif
    !
!    if (mod(j,10) == 0) then
!      forcephimean = 0._rp
!      forcethetamean = 0._rp
!      !$omp parallel default(shared) &
!      !$omp private(l) &
!      !$omp reduction(+:forcephimean,forcethetamean,forcephivar,forcethetavar)
!      !$omp do
!      do l=1,NL
!        forcephimean = forcephimean + forcephi(l)
!        forcethetamean = forcethetamean + forcetheta(l)
!      enddo
!      !$omp end do
!      !$omp end parallel
!      forcephimean = forcephimean/(1._rp*NL)
!      forcethetamean = forcethetamean/(1._rp*NL)
!      forcephivar = 0._rp
!      forcethetavar = 0._rp
!      !$omp parallel default(shared) &
!      !$omp private(l) &
!      !$omp reduction(+:forcephivar,forcethetavar)
!      !$omp do
!      do l=1,NL
!        forcephivar = forcephivar + (forcephi(l)-forcephimean)**2
!        forcethetavar = forcethetavar + (forcetheta(l)-forcethetamean)**2
!      enddo
!      !$omp end do
!      !$omp end parallel
!      forcephivar = forcephivar/(1._rp*NL-1._rp)
!      forcethetavar = forcethetavar/(1._rp*NL-1._rp)
!      dist1max = 9999._rp
!      do l=2,NL
!        dist = sqrt( (xfp(l)-xfp(1))**2 + (yfp(l)-yfp(1))**2 + (zfp(l)-zfp(1))**2 )
!        if (dist < dist1max) then
!          lmin1    = l
!          dist1max = dist
!        endif
!      enddo
!      dist2max = 9999._rp
!      do l=2,NL
!        dist = sqrt( (xfp(l)-xfp(1))**2 + (yfp(l)-yfp(1))**2 + (zfp(l)-zfp(1))**2 )
!        if ( (dist < dist2max) .and. (l /= lmin1) ) then
!          lmin2    = l
!          dist2max = dist
!        endif
!      enddo
!      !open(22,file=datadir//'variances',position='append')
!      !thetanew = theta(NL)-theta(1) !position of l=NL wrt l=1
!      !phinew   = phi(NL)-phi(1)     !position of l=NL wrt l=1
!      !write(22,'(I5,7E16.8,I4,E16.8,I4,E16.8)') j,forcephimean,forcephivar,forcethetamean,forcethetavar, &
!      !                                          radfp*sin(thetanew)*cos(phinew), &
!      !                                          radfp*sin(thetanew)*sin(phinew), &
!      !                                          radfp*cos(thetanew),lmin1,dist1max,lmin2,dist2max
!      !close(22)
!    endif
    !
    if (mod(j,50) == 0) then
      call random_number(rn)
      open(42,file=datadir//'lagrangianforcepoints_shell')
      write(42,*) 'VARIABLES = "lfpx","lfpy","lfpz","theta","phi","radfp"'
      write(42,*) 'ZONE T="Zone1"',' I=',NL,', F=POINT'
      write(42,*) ''
      !
      do l=1,NL
        write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp
      enddo
      !
      close(42)
      call loadd(1,j,phi,theta)
    endif
  enddo
end program init_force_points
