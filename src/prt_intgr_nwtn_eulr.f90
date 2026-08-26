module prt_mod_intgr_nwtn_eulr
#if defined(_PARTICLE)
  use mpi
  use decomp_2d                  , only: nproc
  use mod_types
  use mod_param                  , only: dims,gacc,rho12,sigma,mu12,iout0d
  use mod_common_mpi             , only: prt_comm_cart,myid,ierr,status
  use prt_mod_common             , only: ep,op,tp,npmax, &
                                         boundfrontmyid,boundleftmyid,neighbor, &
                                         nl,npmstr,pmax,nla, &
                                         offset,radfp,coords, &
                                         phirc,thetarc
  use prt_mod_param              , only: colthr_pw,eps_cut_pp,eps_cut_pw, &
                                         eps_ini_pp,eps_ini_pw, &
                                         np,nqmax, &
                                         radius, &
                                         send_int,send_real, &
                                         r_dtcol,r_dtcoli,rho_s
  use prt_mod_intgr_over_sphere  , only: intgr_over_sphere
  use mod_collisions             , only: collisions,lubrication
  !
  implicit none
  !
  private
  !
  public intgr_nwtn_eulr
  !
  contains
  !
  subroutine intgr_nwtn_eulr(n,l,dl,dli,dt,rkpar,istep,psi,u,v,w,Fstot,Fstot_old,F_ibm,F_inertia,F_w,F_buoy,F_cap)
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in) :: dt
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in)               :: istep
    real(rp), intent(in), dimension(0:,0:,0:) :: psi
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(in), dimension(3) :: Fstot
    real(rp), intent(in), dimension(3) :: Fstot_old
    real(rp), intent(inout)            :: F_ibm,F_inertia,F_w,F_buoy,F_cap
    integer, parameter :: csv_unit = 5555
    real(rp) :: rkcoeffab
    real(rp) :: F_sup
    integer :: p,q,botw,topw,nb,nbsend,nbrec,iter,itermax,r
    real(rp) :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    real(rp), dimension(npmax) :: posx,posy,posz,postheta,posphi, &
                                  velx,vely,velz, &
                                  omgx,omgy,omgz
    integer, dimension(npmax) :: colrank
    integer :: sumcolrank,sumcolrank_all,rankmax
    real(rp), dimension(2,1) ::maxerr,maxerr_all
    real(rp) :: err,maxerror,coll_toll
    !
    integer, dimension(0:8) :: pmax_nb
    integer, dimension(1:npmax,0:8) :: mslv_nb,newmaster_nb
    ! newmaster:
    !   > 0 : nr of neighbor that has become the new master of particle p
    !     0 : a) myid didn't contain particle number p, or
    !         b) myid remains master of this particle, or
    !         c) myid remains slave of this particle
    real(rp), dimension(npmax) :: colflgx,colflgy,colflgz
    real(rp) :: deltax,deltay,deltaz,deltan,dist,nx,ny,nz
    !
    integer :: nprocs
    integer :: k,i
    type neighbour
       real(rp) :: x,y,z,theta,phi,u,v,w,omx,omy,omz
    end type neighbour
    type(neighbour), dimension(1:npmax,0:8) :: anb
    !type neighbour2
    !  real :: x,y,z
    !end type neighbour2
    !type(neighbour2), dimension(1:npmax,0:8) :: anb2
    integer :: tag
    integer :: nrrequests
    !integer :: arrayrequests(1:30)
    !integer :: arraystatuses(MPI_STATUS_SIZE,1:30)
    integer :: arrayrequests(1:2)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:2)
    real(rp) :: ax,ay
    integer :: lp
    real(rp) :: leftbound,rightbound,frontbound,backbound
    integer :: nbrec2
    integer, dimension(2) :: proccoords !dimension(size of dims)
    integer :: procrank
!    integer :: counter !delete after checking things
    real(rp) :: eps
    integer :: count_mstr,count_slve,count_mstr_all,count_slve_loc
    integer :: idp,idp_nb,idq
    logical :: found_mstr
    integer :: qlast,qq
!    character(len=3) :: rankpr
    real(rp), dimension(3) :: arrdeltax,arrdeltay
    integer :: ideltax,ideltay
!    real :: isperiodx,isperiody
    real(rp) :: xfploc,yfploc,zfploc
    real(rp) :: isperiodx,isperiody
    real(rp) :: coorxfp,cooryfp   !,coorzfp
    logical :: isout
    integer :: ll
    real(rp) :: dtp
    real(rp) :: We, Re
!    real(rp) :: fxcont, intucont
    !
    nprocs = dims(1)*dims(2)
    rkcoeffab=rkpar(1)+rkpar(2)
    !
    dtp = dt*r_dtcoli
    !
    ! Initialization: new --> old
    !
    !$omp parallel default(none) &
    !$omp shared(ap,op,pmax) &
    !$omp private(p)
    !$omp do 
    
    do p=1,pmax
      if (ep(p)%mslv > 0) then
        op(p)%x = ep(p)%x
        op(p)%y = ep(p)%y
        op(p)%z = ep(p)%z
        op(p)%theta = ep(p)%theta
        op(p)%phi = ep(p)%phi
        op(p)%u = ep(p)%u
        op(p)%v = ep(p)%v
        op(p)%w = ep(p)%w
        op(p)%omx = ep(p)%omx
        op(p)%omy = ep(p)%omy
        op(p)%omz = ep(p)%omz
        op(p)%omtheta = ep(p)%omtheta
        op(p)%intu = ep(p)%intu
        op(p)%intv = ep(p)%intv
        op(p)%intw = ep(p)%intw
        op(p)%intomx = ep(p)%intomx
        op(p)%intomy = ep(p)%intomy
        op(p)%intomz = ep(p)%intomz
        op(p)%colfx = ep(p)%colfx
        op(p)%colfy = ep(p)%colfy
        op(p)%colfz = ep(p)%colfz
        op(p)%coltx = ep(p)%coltx
        op(p)%colty = ep(p)%colty
        op(p)%coltz = ep(p)%coltz
        op(p)%fxltot = ep(p)%fxltot
        op(p)%fyltot = ep(p)%fyltot
        op(p)%fzltot = ep(p)%fzltot
        op(p)%torqxltot = ep(p)%torqxltot
        op(p)%torqyltot = ep(p)%torqyltot
        op(p)%torqzltot = ep(p)%torqzltot
        op(p)%fcapx = ep(p)%fcapx
        op(p)%fcapy = ep(p)%fcapy
        op(p)%fcapz = ep(p)%fcapz
        op(p)%dx(:) = ep(p)%dx(:) 
        op(p)%dy(:) = ep(p)%dy(:)
        op(p)%dz(:) = ep(p)%dz(:)
        op(p)%dxt(:) = ep(p)%dxt(:)
        op(p)%dyt(:) = ep(p)%dyt(:)
        op(p)%dzt(:) = ep(p)%dzt(:)
        op(p)%dut(:) = ep(p)%dut(:)
        op(p)%dvt(:) = ep(p)%dvt(:)
        op(p)%dwt(:) = ep(p)%dwt(:)
      endif
    enddo
    !$omp end parallel
    !
    ! compute integral of linear and angular momentum over sphere
    !
    !call intgr_over_sphere(ap(:)%intu,ap(:)%intv,ap(:)%intw,2)
    !call intgr_over_sphere(ap(:)%intomx,ap(:)%intomy,ap(:)%intomz,3)
!    call intgr_over_sphere(1,n,u,v,w)
    call intgr_over_sphere(2,n,psi,u,v,w)
    call intgr_over_sphere(3,n,psi,u,v,w)
    call intgr_over_sphere(4,n,psi,u,v,w)
    !
    ! exchange data with neighbors: pmasterslave
    !
    pmax_nb(:) = 0
    pmax_nb(0)=pmax
    !$omp workshare
    mslv_nb(:,:) = 0
    mslv_nb(1:pmax,0) = ep(1:pmax)%mslv
    !$omp end workshare
    do nb=1,8
      nbsend = 4+nb
      nbrec  = nb
      if (nbsend > 8) nbsend = nbsend - 8
      call MPI_SENDRECV(pmax_nb(0),1,MPI_INTEGER,neighbor(nbsend),1, &
                        pmax_nb(nbrec),1,MPI_INTEGER,neighbor(nbrec),1, &
                        prt_comm_cart,status,ierr)
      call MPI_SENDRECV(mslv_nb(1,0),pmax_nb(0),MPI_INTEGER,neighbor(nbsend),2, &
                        mslv_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),2, &
                        prt_comm_cart,status,ierr)
    enddo
    !
    ! begin sub-cycling loop: repeat the contact-convergence + Newton-Euler update
    ! r_dtcol times per macro step, each using the sub-step dtp = dt/r_dtcol
    !
    do r = 1,r_dtcol
    !
    ! begin iterative loop
    !
    iter = -1
    itermax = 10
    coll_toll = 1.0e-8_rp
    maxerror = 1
    sumcolrank_all = 1
    do while ((iter < itermax) .and. ((sumcolrank_all+npmax*Nproc) /= 0) .and. (maxerror*dli(1) > coll_toll))
      iter = iter + 1
      !$omp workshare
      colrank(1:npmax) = -1 ! -1 means that particle p is not involved in a collision
      !$omp end workshare

      !$omp parallel default(none) &
      !$omp shared(ap,anb,pmax) &
      !$omp shared(posx,posy,posz,postheta,posphi,velx,vely,velz,omgx,omgy,omgz) &
      !$omp private(p) &
      !$omp firstprivate(iter)  
      !$omp do 
      do p=1,pmax
        if (ep(p)%mslv > 0) then
          !myid is master of particle ep(p)%mslv
          posx(p)     = ep(p)%x
          posy(p)     = ep(p)%y
          posz(p)     = ep(p)%z
          postheta(p) = ep(p)%theta
          posphi(p)   = ep(p)%phi
          velx(p)     = ep(p)%u
          vely(p)     = ep(p)%v
          velz(p)     = ep(p)%w
          omgx(p)     = ep(p)%omx
          omgy(p)     = ep(p)%omy
          omgz(p)     = ep(p)%omz
        else
          posx(p)     = 0.0_rp
          posy(p)     = 0.0_rp
          posz(p)     = 0.0_rp
          postheta(p) = 0.0_rp 
          posphi(p)   = 0.0_rp 
          velx(p)     = 0.0_rp
          vely(p)     = 0.0_rp
          velz(p)     = 0.0_rp
          omgx(p)     = 0.0_rp
          omgy(p)     = 0.0_rp
          omgz(p)     = 0.0_rp ! I think I can simply put anb(p,0)%x here, check later
        endif
      enddo
      !
      ! exchange data with neighbors: posx,posy,posz
      !
      !$omp do  
      do p=1,pmax
        anb(p,0)%x = posx(p)
        anb(p,0)%y = posy(p)
        anb(p,0)%z = posz(p)
        anb(p,0)%theta = postheta(p) 
        anb(p,0)%phi   = posphi(p)
        anb(p,0)%u = velx(p)
        anb(p,0)%v = vely(p)
        anb(p,0)%w = velz(p)
        anb(p,0)%omx = omgx(p)
        anb(p,0)%omy = omgy(p)
        anb(p,0)%omz = omgz(p)
      enddo
      !$omp end parallel
      do nb=1,8
        nbsend = 4+nb
        nbrec  = nb
        if (nbsend > 8) nbsend = nbsend - 8
        call MPI_SENDRECV(anb(1,0)%x,pmax_nb(0)*11,MPI_REAL_RP,neighbor(nbsend),3, &
                          anb(1,nbrec)%x,pmax_nb(nbrec)*11,MPI_REAL_RP,neighbor(nbrec),3, &
                          prt_comm_cart,status,ierr)
        ! send x,y,z,u,v,w,omx,omy,omz -> 9*pmax contiguous info
        ! (see definition of type neighbor in the begining of the subroutine)
      enddo
      !
      ! recompute particle positions because of periodic b.c.'s in the x and y-direction
      !
      !$omp parallel default(none) &
      !$omp shared(ap,anb,pmax,pmax_nb,mslv_nb,coords) &
      !$omp private(p,nb,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)
      nb=1
      !$omp do
      do p=1,pmax_nb(nb)
        if (mslv_nb(p,nb) > 0) then
          boundleftnb  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! left  boundary of neighbor nb
          if (anb(p,nb)%x < boundleftnb) then
            anb(p,nb)%x = anb(p,nb)%x + l(1)
          endif
        endif
      enddo
      nb=2
      !$omp do
      do p=1,pmax_nb(nb)
        if (mslv_nb(p,nb) > 0) then
          boundleftnb  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! left  boundary of neighbor nb
          boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back  boundary of neighbor nb
          if (anb(p,nb)%x < boundleftnb) then
            anb(p,nb)%x = anb(p,nb)%x + l(1)
          endif
          if (anb(p,nb)%y > boundbacknb) then
            anb(p,nb)%y = anb(p,nb)%y - l(2)
          endif
        endif
      enddo
      nb=3
      !$omp do
      do p=1,pmax_nb(nb)
        if (mslv_nb(p,nb) > 0) then
          boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back  boundary of neighbor nb
          if (anb(p,nb)%y > boundbacknb) then
            anb(p,nb)%y = anb(p,nb)%y - l(2)
          endif
        endif
      enddo
      nb=4
      !$omp do
      do p=1,pmax_nb(nb)
        if (mslv_nb(p,nb) > 0) then
          boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
          boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back  boundary of neighbor nb
          if (anb(p,nb)%x > boundrightnb) then
            anb(p,nb)%x = anb(p,nb)%x - l(1)
          endif
          if (anb(p,nb)%y > boundbacknb) then
            anb(p,nb)%y = anb(p,nb)%y - l(2)
          endif
        endif
      enddo
      nb=5
      !$omp do
      do p=1,pmax_nb(nb)
        if (mslv_nb(p,nb) > 0) then
          boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
          if (anb(p,nb)%x > boundrightnb) then
            anb(p,nb)%x = anb(p,nb)%x - l(1)
          endif
        endif
      enddo
      nb=6
      !$omp do
      do p=1,pmax_nb(nb)
        if (mslv_nb(p,nb) > 0) then
          boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
          boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
          if (anb(p,nb)%x > boundrightnb) then
            anb(p,nb)%x = anb(p,nb)%x - l(1)
          endif
          if (anb(p,nb)%y < boundfrontnb) then
            anb(p,nb)%y = anb(p,nb)%y + l(2)
          endif
        endif
      enddo
      nb=7
      !$omp do
      do p=1,pmax_nb(nb)
        if (mslv_nb(p,nb) > 0) then
          boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
          if (anb(p,nb)%y < boundfrontnb) then
            anb(p,nb)%y = anb(p,nb)%y + l(2)
          endif
        endif
      enddo
      nb=8
      !$omp do
      do p=1,pmax_nb(nb)
        if (mslv_nb(p,nb) > 0) then
          boundleftnb  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! left  boundary of neighbor nb
          boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
          if (anb(p,nb)%x < boundleftnb) then
            anb(p,nb)%x = anb(p,nb)%x + l(1)
          endif
          if (anb(p,nb)%y < boundfrontnb) then
            anb(p,nb)%y = anb(p,nb)%y + l(2)
          endif
        endif
      enddo
      !$omp end parallel
      !
      ! collision model
      !
      botw = np + 1 ! id of the bottom wall
      topw = np + 2 ! id of the top wall
      !
      !$omp parallel default(none) &
      !$omp shared(ap,op,anb,pmax) &
      !$omp shared(pmax_nb,mslv_nb) &
      !$omp shared(colflgx,colflgy,colflgz,colrank) &
      !$omp shared(itermax,myid,dt,rkiter,topw,botw) &
      !$omp private(p,q,nb,dist,deltax,deltay,deltaz,deltan,eps) &
      !$omp private(nx,ny,nz) &
      !$omp private(idp,idq,qq,qlast) &
      !$omp private(arrdeltax,arrdeltay) &
      !$omp private(ideltax,ideltay) &
      !$omp firstprivate(iter)
      !$omp do
      do p=1,pmax
        qlast = ep(p)%qmax
        colflgx(p) = 0.0_rp
        colflgy(p) = 0.0_rp
        colflgz(p) = 0.0_rp
        if (ep(p)%mslv > 0) then
          idp = ep(p)%mslv
          qlast = ep(p)%qmax
          ep(p)%colfx = 0.0_rp
          ep(p)%colfy = 0.0_rp
          ep(p)%colfz = 0.0_rp
          ep(p)%coltx = 0.0_rp
          ep(p)%colty = 0.0_rp
          ep(p)%coltz = 0.0_rp
          !
          ! (i) collision with other particles
          !
          do nb=0,8
            do q=1,pmax_nb(nb)
              if((nb >= 0).and.(idp == mslv_nb(q,nb))) then
                dist = 0.0_rp ! dummy operation
              elseif(mslv_nb(q,nb) > 0) then
                idq = mslv_nb(q,nb)
                !          deltax = anb(q,nb)%x - ap(p)%x
                !          deltay = anb(q,nb)%y - ap(p)%y
                deltaz = anb(q,nb)%z - ep(p)%z
                arrdeltax(1) = anb(q,nb)%x - ep(p)%x
                arrdeltax(2) = anb(q,nb)%x - l(1) - ep(p)%x
                arrdeltax(3) = anb(q,nb)%x + l(1) - ep(p)%x
                ideltax = minloc(abs(arrdeltax(1:3)),1)
                deltax = arrdeltax(ideltax)
                !
                arrdeltay(1) = anb(q,nb)%y - ep(p)%y
                arrdeltay(2) = anb(q,nb)%y - l(2) - ep(p)%y
                arrdeltay(3) = anb(q,nb)%y + l(2) - ep(p)%y
                ideltay = minloc(abs(arrdeltay(1:3)),1)
                deltay = arrdeltay(ideltay)
                !
                dist = sqrt(deltax**2+deltay**2+deltaz**2)
                deltan = 2*radius-dist
                nx = deltax/dist
                ny = deltay/dist
                nz = deltaz/dist ! computed twice (here and in the subroutine collisions)
                eps = (dist-2.0_rp*radius)/radius
                if((eps < eps_ini_pp).and.(eps > eps_cut_pp)) then
                  colrank(p) = myid
                  call lubrication(p,idq,nx,ny,nz,eps, &
                                   anb(q,nb)%u,anb(q,nb)%v,anb(q,nb)%w, &
                                   anb(q,nb)%omx,anb(q,nb)%omy,anb(q,nb)%omz)
                endif
                qq = 0
                !
                ! Matching global id 'idq' with local id 'qq' of the particle in contact.
                ! This has to be done because of the oblique collison model, which has 
                ! memory: psi must be fixed at first contact and the tangential displacement
                ! has to be integrated from the imminence of contact
                !
                do i=1,ep(p)%qmax
                  if(idq == ep(p)%firstc(i)) then
                    !
                    ! ep(p)%firstc(i) == idq -> particles with ids 'idp' and 'idq' were in
                    ! contact in the previous substep
                    !
                    if(dist < (2.0_rp*radius)) then 
                      !
                      ! particles are still in contact
                      !
                      qq = i
                      call collisions(rkpar,p,qq,idq,dist,deltax,deltay,deltaz, &
                                      anb(q,nb)%u,anb(q,nb)%v,anb(q,nb)%w, &
                                      anb(q,nb)%omx,anb(q,nb)%omy,anb(q,nb)%omz,dtp)
                      !                if (abs(deltan*nx).gt.(colthr_pp)) colflgx(p) = 0. ! IBM force off
                      !                if (abs(deltan*ny).gt.(colthr_pp)) colflgy(p) = 0. ! IBM force off
                      !                if (abs(deltan*nz).gt.(colthr_pp)) colflgz(p) = 0. ! IBM force off
                      colrank(p) = myid
                    elseif(iter == itermax) then
                      !
                      ! particles ceased to be in contact, remove local id and put the last
                      ! element in position 'i' to avoid 'gaps' in the local contact arrays
                      !
                      if(qlast > 1) then
                        ep(p)%firstc(i) = ep(p)%firstc(qlast)
                        ep(p)%dxt(i)    = ep(p)%dxt(qlast)
                        ep(p)%dyt(i)    = ep(p)%dyt(qlast)
                        ep(p)%dzt(i)    = ep(p)%dzt(qlast)
                        ep(p)%dut(i)    = ep(p)%dut(qlast)
                        ep(p)%dvt(i)    = ep(p)%dvt(qlast)
                        ep(p)%dwt(i)    = ep(p)%dwt(qlast)
                        op(p)%firstc(i) = op(p)%firstc(qlast)
                        op(p)%dxt(i)    = op(p)%dxt(qlast)
                        op(p)%dyt(i)    = op(p)%dyt(qlast)
                        op(p)%dzt(i)    = op(p)%dzt(qlast)
                        op(p)%dut(i)    = op(p)%dut(qlast)
                        op(p)%dvt(i)    = op(p)%dvt(qlast)
                        op(p)%dwt(i)    = op(p)%dwt(qlast)
                      endif
                      ep(p)%firstc(qlast) = 0
                      ep(p)%dx(qlast) = 0.0_rp
                      ep(p)%dy(qlast) = 0.0_rp
                      ep(p)%dz(qlast) = 0.0_rp
                      ep(p)%dxt(qlast) = 0.0_rp
                      ep(p)%dyt(qlast) = 0.0_rp
                      ep(p)%dzt(qlast) = 0.0_rp
                      ep(p)%dut(qlast) = 0.0_rp
                      ep(p)%dvt(qlast) = 0.0_rp
                      ep(p)%dwt(qlast) = 0.0_rp
                      op(p)%firstc(qlast) = 0
                      op(p)%dx(qlast) = 0.0_rp
                      op(p)%dy(qlast) = 0.0_rp
                      op(p)%dz(qlast) = 0.0_rp
                      op(p)%dxt(qlast) = 0.0_rp
                      op(p)%dyt(qlast) = 0.0_rp
                      op(p)%dzt(qlast) = 0.0_rp
                      op(p)%dut(qlast) = 0.0_rp
                      op(p)%dvt(qlast) = 0.0_rp
                      op(p)%dwt(qlast) = 0.0_rp
                      qlast = qlast - 1
                    endif
                  endif
                enddo
                if(qq == 0 .and. dist < (2.0_rp*radius)) then 
                  !
                  ! particle 'idp' in contact with particle 'idq' for the first time:
                  ! increase extent of local contact array and initialize it
                  !
                  qlast = qlast + 1
                  qq = qlast
                  ep(p)%dx(qlast) = 0.0_rp
                  ep(p)%dy(qlast) = 0.0_rp
                  ep(p)%dz(qlast) = 0.0_rp
                  ep(p)%dxt(qlast) = 0.0_rp
                  ep(p)%dyt(qlast) = 0.0_rp
                  ep(p)%dzt(qlast) = 0.0_rp
                  ep(p)%dut(qlast) = 0.0_rp
                  ep(p)%dvt(qlast) = 0.0_rp
                  ep(p)%dwt(qlast) = 0.0_rp
                  op(p)%dx(qlast) = 0.0_rp
                  op(p)%dy(qlast) = 0.0_rp
                  op(p)%dz(qlast) = 0.0_rp
                  op(p)%dxt(qlast) = 0.0_rp
                  op(p)%dyt(qlast) = 0.0_rp
                  op(p)%dzt(qlast) = 0.0_rp
                  op(p)%dut(qlast) = 0.0_rp
                  op(p)%dvt(qlast) = 0.0_rp
                  op(p)%dwt(qlast) = 0.0_rp
                  call collisions(rkpar,p,qq,idq,dist,deltax,deltay,deltaz, &
                                  anb(q,nb)%u,anb(q,nb)%v,anb(q,nb)%w, &
                                  anb(q,nb)%omx,anb(q,nb)%omy,anb(q,nb)%omz,dtp)
                  !            if (abs(deltan*nx).gt.(colthr_pp)) colflgx(p) = 0. ! IBM force off
                  !            if (abs(deltan*ny).gt.(colthr_pp)) colflgy(p) = 0. ! IBM force off
                  !            if (abs(deltan*nz).gt.(colthr_pp)) colflgz(p) = 0. ! IBM force off
                  ep(p)%firstc(qlast) = idq
                  op(p)%firstc(qlast) = idq
                  colrank(p) = myid
                endif
                ep(p)%qmax = qlast ! size of local contact arrays updated
              endif
            enddo
          enddo
          !
          ! (ii) collision with walls
          ! change velocities of the wall if the wall has a prescribed velocity (e.g. Couette flow)
          !
          do q=botw,topw
            idq = q
            deltax = 0.0_rp
            deltay = 0.0_rp
            deltaz = (q-botw)*l(3) - ep(p)%z
            dist = sqrt(deltax**2+deltay**2+deltaz**2)
            deltan = radius-dist
            nx = deltax/dist
            ny = deltay/dist
            nz = deltaz/dist
            qq = 0
            eps = -deltan/radius 
            if((eps < eps_ini_pw) .and. (eps > eps_cut_pw)) then
              colrank(p) = myid
              call lubrication(p,idq,nx,ny,nz,eps, &
                               0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp)
            endif
            do i=1,ep(p)%qmax
              if(idq == ep(p)%firstc(i)) then
                !
                ! algorithm is analogous for particle-wall interactions, see comments
                ! for particle-particle interactions 
                !
                qq = i
                if(dist < radius) then
                  call collisions(rkpar,p,qq,idq,dist,deltax,deltay,deltaz, &
                                  0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,dtp)
                  if (dist < (radius-colthr_pw)) colflgz(p) = 1.0_rp !IBM force off
                  colrank(p) = myid
                elseif(iter == itermax) then
                  if(qlast > 1) then
                    ep(p)%firstc(i) = ep(p)%firstc(qlast)
                    ep(p)%dxt(i)    = ep(p)%dxt(qlast)
                    ep(p)%dyt(i)    = ep(p)%dyt(qlast)
                    ep(p)%dzt(i)    = ep(p)%dzt(qlast)
                    ep(p)%dut(i)    = ep(p)%dut(qlast)
                    ep(p)%dvt(i)    = ep(p)%dvt(qlast)
                    ep(p)%dwt(i)    = ep(p)%dwt(qlast)
                    op(p)%firstc(i) = op(p)%firstc(qlast)
                    op(p)%dxt(i)    = op(p)%dxt(qlast)
                    op(p)%dyt(i)    = op(p)%dyt(qlast)
                    op(p)%dzt(i)    = op(p)%dzt(qlast)
                    op(p)%dut(i)    = op(p)%dut(qlast)
                    op(p)%dvt(i)    = op(p)%dvt(qlast)
                    op(p)%dwt(i)    = op(p)%dwt(qlast)
                  endif
                  ep(p)%firstc(qlast) = 0
                  ep(p)%dx(qlast) = 0.0_rp
                  ep(p)%dy(qlast) = 0.0_rp
                  ep(p)%dz(qlast) = 0.0_rp
                  ep(p)%dxt(qlast) = 0.0_rp
                  ep(p)%dyt(qlast) = 0.0_rp
                  ep(p)%dzt(qlast) = 0.0_rp
                  ep(p)%dut(qlast) = 0.0_rp
                  ep(p)%dvt(qlast) = 0.0_rp
                  ep(p)%dwt(qlast) = 0.0_rp
                  op(p)%firstc(qlast) = 0
                  op(p)%dx(qlast) = 0.0_rp
                  op(p)%dy(qlast) = 0.0_rp
                  op(p)%dz(qlast) = 0.0_rp
                  op(p)%dxt(qlast) = 0.0_rp
                  op(p)%dyt(qlast) = 0.0_rp
                  op(p)%dzt(qlast) = 0.0_rp
                  op(p)%dut(qlast) = 0.0_rp
                  op(p)%dvt(qlast) = 0.0_rp
                  op(p)%dwt(qlast) = 0.0_rp
                  qlast = qlast - 1
                endif
              endif
            enddo
            if((qq == 0) .and. (dist < radius)) then
               qlast = qlast + 1
               qq = qlast
               ep(p)%dx(qlast) = 0.0_rp
               ep(p)%dy(qlast) = 0.0_rp
               ep(p)%dz(qlast) = 0.0_rp
               ep(p)%dxt(qlast) = 0.0_rp
               ep(p)%dyt(qlast) = 0.0_rp
               ep(p)%dzt(qlast) = 0.0_rp
               ep(p)%dut(qlast) = 0.0_rp
               ep(p)%dvt(qlast) = 0.0_rp
               ep(p)%dwt(qlast) = 0.0_rp
               op(p)%dx(qlast) = 0.0_rp
               op(p)%dy(qlast) = 0.0_rp
               op(p)%dz(qlast) = 0.0_rp
               op(p)%dxt(qlast) = 0.0_rp
               op(p)%dyt(qlast) = 0.0_rp
               op(p)%dzt(qlast) = 0.0_rp
               op(p)%dut(qlast) = 0.0_rp
               op(p)%dvt(qlast) = 0.0_rp
               op(p)%dwt(qlast) = 0.0_rp
               call collisions(rkpar,p,qq,idq,dist,deltax,deltay,deltaz, &
                               0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,dtp)
               if ((dist) < (radius-colthr_pw)) colflgz(p) = 1.0_rp !switch off IBM force
               ep(p)%firstc(qlast) = idq
               op(p)%firstc(qlast) = idq
               colrank(p) = myid
            endif
          enddo
          ep(p)%qmax = qlast
        endif
      enddo
      !
      !compute new particle positions
      !
!!$omp do schedule(dynamic)
      !$omp do 
      do p=1,pmax
        if (ep(p)%mslv > 0) then
          if (colrank(p) /= -1) then ! Particle with collision => dtp
            ep(p)%u = op(p)%u + &
                      (1.0_rp-colflgx(p))*( &
                      (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%fxltot+op(p)%fxltot)/(ep(p)%vol*rho_s) + &
                      r_dtcoli*(ep(p)%intu-op(p)%intu)/(ep(p)%vol*rho_s)) + &
                      rkcoeffab*dtp*gacc(1)*(1.0_rp-(ep(p)%intrhox/(ep(p)%vol*rho_s))) + &
!                      (-1._rp)*rkcoeffab*dtp*0.5_rp*(ep(p)%fcapx+op(p)%fcapx)/(ep(p)%vol*rho_s) + &
                      rkcoeffab*0.5_rp*dtp*(ep(p)%colfx+op(p)%colfx)/(ep(p)%vol*ep(p)%ratiorho) !+ &
!                      rkcoeffab*dtp*0.5_rp*(Fstot(1)+Fstot_old(1))/(ep(p)%vol*rho_s)
            ep(p)%x = op(p)%x + rkcoeffab*dtp*0.5_rp*(ep(p)%u+op(p)%u)
!            PRINT *, "ep(p)%x", ep(p)%x
            ep(p)%v = op(p)%v + &
                      (1.0_rp-colflgy(p))*( &
                      (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%fyltot+op(p)%fyltot)/(ep(p)%vol*rho_s) + &
                      r_dtcoli*(ep(p)%intv-op(p)%intv)/(ep(p)%vol*rho_s)) + &
                      rkcoeffab*dtp*gacc(2)*(1.0_rp-(ep(p)%intrhoy/(ep(p)%vol*rho_s))) + &
!                      (-1._rp)*rkcoeffab*dtp*0.5_rp*(ep(p)%fcapy+op(p)%fcapy)/(ep(p)%vol*rho_s) + &
                      rkcoeffab*0.5_rp*dtp*(ep(p)%colfy+op(p)%colfy)/(ep(p)%vol*ep(p)%ratiorho) !+ &
!                      rkcoeffab*dtp*0.5_rp*(Fstot(2)+Fstot_old(2))/(ep(p)%vol*rho_s)
            ep(p)%y = op(p)%y + rkcoeffab*dtp*0.5_rp*(ep(p)%v+op(p)%v)
            ep(p)%w = op(p)%w + &
                      (1.0_rp-colflgz(p))*( &
                      (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%fzltot+op(p)%fzltot)/(ep(p)%vol*rho_s) + &
                      r_dtcoli*(ep(p)%intw-op(p)%intw)/(ep(p)%vol*rho_s)) + &
                      rkcoeffab*dtp*gacc(3)*(1.0_rp-(ep(p)%intrhoz/(ep(p)%vol*rho_s))) + &
!                      (-1._rp)*rkcoeffab*dtp*0.5_rp*(ep(p)%fcapz+op(p)%fcapz)/(ep(p)%vol*rho_s) + &
                      rkcoeffab*0.5_rp*dtp*(ep(p)%colfz+op(p)%colfz)/(ep(p)%vol*ep(p)%ratiorho) !+ &
!                      rkcoeffab*dtp*0.5_rp*(Fstot(3)+Fstot_old(3))/(ep(p)%vol*rho_s)
            ep(p)%z = op(p)%z + 0.5_rp*rkcoeffab*dtp*(ep(p)%w+op(p)%w)
            ep(p)%omx = op(p)%omx + &
                        (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%torqxltot+op(p)%torqxltot)/(ep(p)%mominert*rho_s) + &
                        r_dtcoli*(ep(p)%intomx-op(p)%intomx)/(ep(p)%mominert*rho_s) !+ &
!                        rkcoeffab*0.5_rp*dtp*(ep(p)%coltx+op(p)%coltx)/(ep(p)%mominert*ep(p)%ratiorho)
            ep(p)%omy = op(p)%omy + &
                        (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%torqyltot+op(p)%torqyltot)/(ep(p)%mominert*rho_s) + &
                        r_dtcoli*(ep(p)%intomy-op(p)%intomy)/(ep(p)%mominert*rho_s) !+ &
!                        rkcoeffab*0.5_rp*dtp*(ep(p)%colty+op(p)%colty)/(ep(p)%mominert*ep(p)%ratiorho)
            ep(p)%omz = op(p)%omz + &
                        (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%torqzltot+op(p)%torqzltot)/(ep(p)%mominert*rho_s) + &
                        r_dtcoli*(ep(p)%intomz-op(p)%intomz)/(ep(p)%mominert*rho_s) !+ &
!                        rkcoeffab*0.5_rp*dtp*(ep(p)%coltz+op(p)%coltz)/(ep(p)%mominert*ep(p)%ratiorho)
            ep(p)%phi   = op(p)%phi + 0.5_rp*rkcoeffab*dtp*(ep(p)%omz+op(p)%omz)
            ep(p)%omtheta = (ep(p)%omy*cos(ep(p)%phi)) - &
                            (ep(p)%omx*sin(ep(p)%phi))
            ep(p)%theta = op(p)%theta + rkcoeffab*dtp*0.5_rp*(ep(p)%omtheta+op(p)%omtheta)
          else ! Particles without collision => dtp
            ep(p)%u = op(p)%u + &
                      (1.0_rp-colflgx(p))*( &
                      (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%fxltot+op(p)%fxltot)/(ep(p)%vol*rho_s) + &
                      r_dtcoli*(ep(p)%intu-op(p)%intu)/(ep(p)%vol*rho_s)) + &
                      rkcoeffab*dtp*gacc(1)*(1.0_rp-(ep(p)%intrhox/(ep(p)%vol*rho_s))) + &
!                      (-1._rp)*rkcoeffab*dtp*0.5_rp*(ep(p)%fcapx+op(p)%fcapx)/(ep(p)%vol*rho_s) + &
                      rkcoeffab*0.5_rp*dtp*(ep(p)%colfx+op(p)%colfx)/(ep(p)%vol*ep(p)%ratiorho) !+ &
!                      rkcoeffab*dtp*0.5_rp*(Fstot(1)+Fstot_old(1))/(ep(p)%vol*rho_s)
            ep(p)%x = op(p)%x + rkcoeffab*dtp*0.5_rp*(ep(p)%u+op(p)%u)
!            PRINT *, "ep(p)%x", ep(p)%x
            ep(p)%v = op(p)%v + &
                      (1.0_rp-colflgy(p))*( &
                      (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%fyltot+op(p)%fyltot)/(ep(p)%vol*rho_s) + &
                      r_dtcoli*(ep(p)%intv-op(p)%intv)/(ep(p)%vol*rho_s)) + &
                      rkcoeffab*dtp*gacc(2)*(1.0_rp-(ep(p)%intrhoy/(ep(p)%vol*rho_s))) + &
!                      (-1._rp)*rkcoeffab*dtp*0.5_rp*(ep(p)%fcapy+op(p)%fcapy)/(ep(p)%vol*rho_s) + &
                      rkcoeffab*0.5_rp*dtp*(ep(p)%colfy+op(p)%colfy)/(ep(p)%vol*ep(p)%ratiorho) !+ &
!                      rkcoeffab*dtp*0.5_rp*(Fstot(2)+Fstot_old(2))/(ep(p)%vol*rho_s)
            ep(p)%y = op(p)%y + rkcoeffab*dtp*0.5_rp*(ep(p)%v+op(p)%v)
            ep(p)%w = op(p)%w + &
                      (1.0_rp-colflgz(p))*( &
                      (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%fzltot+op(p)%fzltot)/(ep(p)%vol*rho_s) + &
                      r_dtcoli*(ep(p)%intw-op(p)%intw)/(ep(p)%vol*rho_s)) + &
                      rkcoeffab*dtp*gacc(3)*(1.0_rp-(ep(p)%intrhoz/(ep(p)%vol*rho_s))) + &
!                      (-1._rp)*rkcoeffab*dtp*0.5_rp*(ep(p)%fcapz+op(p)%fcapz)/(ep(p)%vol*rho_s) + &
                      rkcoeffab*0.5_rp*dtp*(ep(p)%colfz+op(p)%colfz)/(ep(p)%vol*ep(p)%ratiorho) !+ &
!                      rkcoeffab*dtp*0.5_rp*(Fstot(3)+Fstot_old(3))/(ep(p)%vol*rho_s)
            ep(p)%z = op(p)%z + 0.5_rp*rkcoeffab*dtp*(ep(p)%w+op(p)%w)
            !
            ep(p)%omx = op(p)%omx + &
                        (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%torqxltot+op(p)%torqxltot)/(ep(p)%mominert*rho_s) + &
                        r_dtcoli*(ep(p)%intomx-op(p)%intomx)/(ep(p)%mominert*rho_s) !+ &
!                        rkcoeffab*0.5_rp*dtp*(ep(p)%coltx+op(p)%coltx)/(ep(p)%mominert*ep(p)%ratiorho)
            ep(p)%omy = op(p)%omy + &
                        (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%torqyltot+op(p)%torqyltot)/(ep(p)%mominert*rho_s) + &
                        r_dtcoli*(ep(p)%intomy-op(p)%intomy)/(ep(p)%mominert*rho_s) !+ &
!                        rkcoeffab*0.5_rp*dtp*(ep(p)%colty+op(p)%colty)/(ep(p)%mominert*ep(p)%ratiorho)
            ep(p)%omz = op(p)%omz + &
                        (-1.0_rp)*rkcoeffab*0.5_rp*dtp*(ep(p)%torqzltot+op(p)%torqzltot)/(ep(p)%mominert*rho_s) + &
                        r_dtcoli*(ep(p)%intomz-op(p)%intomz)/(ep(p)%mominert*rho_s) !+ &
!                        rkcoeffab*0.5_rp*dtp*(ep(p)%coltz+op(p)%coltz)/(ep(p)%mominert*ep(p)%ratiorho)
            ep(p)%phi   = op(p)%phi + 0.5_rp*rkcoeffab*dtp*(ep(p)%omz+op(p)%omz)
            ep(p)%omtheta = (ep(p)%omy*cos(ep(p)%phi)) - &
                            (ep(p)%omx*sin(ep(p)%phi))
            ep(p)%theta = op(p)%theta + rkcoeffab*dtp*0.5_rp*(ep(p)%omtheta+op(p)%omtheta)
          endif
        endif
      enddo
      !$omp end parallel
      !
      !Check whether a collision occured: (sumcolrank_all + pmax*Nproc) .ne. 0
      !
      sumcolrank = sum(colrank(1:npmax))
      call MPI_ALLREDUCE(sumcolrank,sumcolrank_all,1,MPI_INTEGER,MPI_SUM, &
                         prt_comm_cart,ierr)
      if ((sumcolrank_all+npmax*Nproc) /= 0) then !criterion for occurence of collisions 
        !THIS HAS TO BE CORRECTED, WILL DO LATER!!!!
        maxerror = 0.0_rp
        do p=1,pmax
          if (ep(p)%mslv > 0) then
            err = sqrt( (posx(p)-ep(p)%x)**2 + (posy(p)-ep(p)%y)**2 + (posz(p)-ep(p)%z)**2 )
            if (err > maxerror) maxerror = err
          endif
        enddo
        maxerr(1,1) = maxerror
        maxerr(2,1) = myid*1.0_rp
        call MPI_ALLREDUCE(maxerr,maxerr_all,1,MPI_2REAL_RP,MPI_MAXLOC, &
                           prt_comm_cart,ierr)
        maxerror = maxerr_all(1,1)
        rankmax  = nint(maxerr_all(2,1))
        !if ( myid == rankmax .and. iter == 1) then
        if (myid == rankmax .and. (iter == 1 .or. (iter > 0 .and. mod(iter,5) == 0))) then
          write(6,*) 'Collision Iter = ',iter,'; (max error)/dx = ',maxerror*dli(1)
        endif
      endif
    enddo ! do while
    !
    ! log position, velocity and forces once per macro step (r == r_dtcol),
    ! after the collision-constraint iteration above has converged -- doing
    ! this inside the do-while would write one row per iteration instead of one.
    !
    if (r == r_dtcol) then
      do p=1,pmax
        if (ep(p)%mslv > 0) then
          F_sup = 0
          F_ibm = 0
          F_inertia = 0
          F_w = 0
          F_buoy = 0
          F_cap = 0
          F_sup = F_sup - rkcoeffab*0.5_rp*(ep(p)%fcapz+op(p)%fcapz)
          F_ibm = F_ibm - ep(p)%fzltot*rkcoeffab
          F_inertia = F_inertia + (ep(p)%intw-op(p)%intw)/dt
          F_w = F_w + (ep(p)%vol*rho_s)*gacc(3)*rkcoeffab
          F_buoy = F_buoy + gacc(3)*(-ep(p)%intrhoz)*rkcoeffab
          F_cap = F_cap + 0.5_rp*(Fstot(3)+Fstot_old(3))*rkcoeffab
          if (MOD(istep, iout0d) == 0) then
            write(csv_unit, '(7(E16.8, ","), E16.8)') &
                  F_sup, F_ibm, F_inertia, F_w, F_buoy, F_cap, ep(p)%z, ep(p)%w
            flush(csv_unit)
          endif
          PRINT *, "ep(p)%z", ep(p)%z
          PRINT *, "ep(p)%w", ep(p)%w
        endif
      enddo
    endif
    !
    ! new --> old, in preparation for the next collision sub-step.
    ! intu,intv,intw,intomx,intomy,intomz,fxltot,fyltot,fzltot,torqxltot,torqyltot,torqzltot,
    ! fcapx,fcapy,fcapz are intentionally left untouched: they are constant over the whole
    ! macro step, and the r_dtcoli-scaled terms above rely on them keeping their pre-loop value.
    !
    do p=1,pmax
      if (ep(p)%mslv > 0) then
        op(p)%x       = ep(p)%x
        op(p)%y       = ep(p)%y
        op(p)%z       = ep(p)%z
        op(p)%theta   = ep(p)%theta
        op(p)%phi     = ep(p)%phi
        op(p)%u       = ep(p)%u
        op(p)%v       = ep(p)%v
        op(p)%w       = ep(p)%w
        op(p)%omx     = ep(p)%omx
        op(p)%omy     = ep(p)%omy
        op(p)%omz     = ep(p)%omz
        op(p)%omtheta = ep(p)%omtheta
        op(p)%colfx   = ep(p)%colfx
        op(p)%colfy   = ep(p)%colfy
        op(p)%colfz   = ep(p)%colfz
        op(p)%coltx   = ep(p)%coltx
        op(p)%colty   = ep(p)%colty
        op(p)%coltz   = ep(p)%coltz
        op(p)%dx(:)   = ep(p)%dx(:)
        op(p)%dy(:)   = ep(p)%dy(:)
        op(p)%dz(:)   = ep(p)%dz(:)
        op(p)%dxt(:)  = ep(p)%dxt(:)
        op(p)%dyt(:)  = ep(p)%dyt(:)
        op(p)%dzt(:)  = ep(p)%dzt(:)
        op(p)%dut(:)  = ep(p)%dut(:)
        op(p)%dvt(:)  = ep(p)%dvt(:)
        op(p)%dwt(:)  = ep(p)%dwt(:)
      endif
    enddo
    enddo ! do r_dtcol
    !
    ! correction for periodic b.c.'s
    !
    !$omp parallel default(none) &
    !$omp shared(ap,pmax) &
    !$omp private(p)
    !$omp do
    do p=1,pmax
      if (ep(p)%mslv > 0) then
        if (ep(p)%x > l(1))   ep(p)%x = ep(p)%x-l(1)
        if (ep(p)%x < 0.0_rp) ep(p)%x = ep(p)%x+l(1)
        if (ep(p)%y > l(2))   ep(p)%y = ep(p)%y-l(2)
        if (ep(p)%y < 0.0_rp) ep(p)%y = ep(p)%y+l(2)
      endif
    enddo
    !$omp end parallel
    !
    ! particle positions were updated. Now check if there are new masters
    !
    do p=1,pmax
      newmaster_nb(p,0) = 0
      if (ep(p)%mslv > 0) then
        ! myid was master of particle ap(p)%mslv at previous time step n
        ax = 0.5_rp
        ay = 0.5_rp
        if (ep(p)%x == l(1))   ax = 0.51_rp
        if (ep(p)%x == 0.0_rp) ax = 0.49_rp
        if (ep(p)%y == l(2))   ay = 0.51_rp
        if (ep(p)%y == 0.0_rp) ay = 0.49_rp
        proccoords(1) = nint( dims(1)*ep(p)%x/l(1) - ax )
        proccoords(2) = nint( dims(2)*ep(p)%y/l(2) - ay ) 
        call MPI_CART_RANK(prt_comm_cart,proccoords,procrank,ierr)
        if (procrank /= myid) then
          ! particle ep(p)%mslv has a new master at time step n+1
          do nb=1,8
            if (procrank == neighbor(nb)) then
              newmaster_nb(p,0) = nb
            endif
          enddo
        endif
      endif
    enddo
    !
    ! exchange data
    !
    do nb=1,8
      nbsend = nb
      nbrec  = nb+4
      if (nbrec > 8) nbrec = nbrec-8
      call MPI_SENDRECV(newmaster_nb(1,0),pmax_nb(0),MPI_INTEGER,neighbor(nbsend),1, &
                        newmaster_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),1, &
                        prt_comm_cart,status,ierr)
    enddo
    !
    do p=1,pmax
      nrrequests = 0
      if (newmaster_nb(p,0) > 0) then
        nbsend = newmaster_nb(p,0)
        tag = ep(p)%mslv!*10+nbsend
        nrrequests = nrrequests + 1
        call MPI_ISEND(ep(p)%x,send_real,MPI_REAL_RP,neighbor(nbsend),tag,prt_comm_cart,arrayrequests((nrrequests-1)*2+1),ierr)
        call MPI_ISEND(ep(p)%qmax,send_int,MPI_INTEGER,neighbor(nbsend),tag+np,prt_comm_cart,arrayrequests((nrrequests-1)*2+2),ierr)
        ep(p)%mslv = -ep(p)%mslv ! master is a slave now
      endif
      if(mslv_nb(p,0) < 0) then
        idp = -ep(p)%mslv
        do nbrec = 1,8
          nbrec2 = nbrec + 4
          if(nbrec2 > 8) nbrec2 = nbrec2-8
          do i=1,pmax_nb(nbrec)
            idp_nb = mslv_nb(i,nbrec)
            if(newmaster_nb(i,nbrec) == nbrec2 .and. idp == idp_nb) then
              ep(p)%mslv = -ep(p)%mslv ! slave became a master
              nrrequests = nrrequests + 1
              tag = ep(p)%mslv!*10+nbsend
              call MPI_IRECV(ep(p)%x,send_real,MPI_REAL_RP,neighbor(nbrec),tag,prt_comm_cart,arrayrequests((nrrequests-1)*2+1),ierr)
              call MPI_IRECV(ep(p)%qmax,send_int,MPI_INTEGER,neighbor(nbrec),tag+np,prt_comm_cart, &
                             arrayrequests((nrrequests-1)*2+2),ierr)
            endif
          enddo
        enddo
      endif
      nrrequests = nrrequests*2
      call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
    enddo
    !
    ! Masters are known now: ap(p)%mslv > 0.  
    ! Next step: determine slaves and master/slave neighbors
    !
    do p=1,pmax
      do k=0,8
        mslv_nb(p,k) = 0
      enddo
    enddo
    !
    mslv_nb(1:pmax,0) = ep(1:pmax)%mslv  ! new masters, but new slaves not all determined yet!
    ! process might remain master, but slave might lose particle
    !$omp workshare
    anb(1:pmax,1:8)%x     = 0.0_rp
    anb(1:pmax,1:8)%y     = 0.0_rp
    anb(1:pmax,1:8)%z     = 0.0_rp
    anb(1:pmax,1:8)%theta = 0.0_rp
    anb(1:pmax,1:8)%phi   = 0.0_rp
    anb(1:pmax,1:8)%u     = 0.0_rp
    anb(1:pmax,1:8)%v     = 0.0_rp
    anb(1:pmax,1:8)%w     = 0.0_rp
    anb(1:pmax,1:8)%omx   = 0.0_rp
    anb(1:pmax,1:8)%omy   = 0.0_rp
    anb(1:pmax,1:8)%omz   = 0.0_rp
    anb(1:pmax,0)%x = ep(1:pmax)%x
    anb(1:pmax,0)%y = ep(1:pmax)%y
    anb(1:pmax,0)%z = ep(1:pmax)%z
    anb(1:pmax,0)%theta = ep(1:pmax)%theta
    anb(1:pmax,0)%phi = ep(1:pmax)%phi
    anb(1:pmax,0)%u = ep(1:pmax)%u
    anb(1:pmax,0)%v = ep(1:pmax)%v
    anb(1:pmax,0)%w = ep(1:pmax)%w
    anb(1:pmax,0)%omx = ep(1:pmax)%omx
    anb(1:pmax,0)%omy = ep(1:pmax)%omy
    anb(1:pmax,0)%omz = ep(1:pmax)%omz
    !$omp end workshare
    do nb=1,8
      nbsend = 4+nb
      if (nbsend > 8) nbsend = nbsend - 8
      nbrec  = nb
      !  nbsend = nb
      !  nbrec  = nb+4
      !  if (nbrec .gt. 8) nbrec = nbrec-8
      call MPI_SENDRECV(mslv_nb(1,0),pmax_nb(0),MPI_INTEGER,neighbor(nbsend),1, &
                        mslv_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),1, &
                        prt_comm_cart,status,ierr)
      call MPI_SENDRECV(anb(1,0)%x,pmax_nb(0)*11,MPI_REAL_RP,neighbor(nbsend),2, &
                        anb(1,nbrec)%x,pmax_nb(nbrec)*11,MPI_REAL_RP,neighbor(nbrec),2, &
                        prt_comm_cart,status,ierr)
      ! send x,y,z -> 3*pmax contiguous info
      ! (see definition of type neighbor2 in the begining of the subroutine)
    enddo
    !
    ! recompute particle positions because of periodic b.c.'s in the x and y-direction
    !
    !$omp parallel default(none) &
    !$omp shared(ap,anb,coords,pmax_nb,mslv_nb)  &
    !$omp private(p,nb,boundleftnb,boundbacknb,boundrightnb,boundfrontnb) 
    nb=1
    !$omp do
    do p=1,pmax_nb(nb)
      if (mslv_nb(p,nb) > 0) then
        boundleftnb  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! left  boundary of neighbor nb
        if (anb(p,nb)%x < boundleftnb) then
          anb(p,nb)%x = anb(p,nb)%x + l(1)
        endif
      endif
    enddo
    nb=2
    !$omp do
    do p=1,pmax_nb(nb)
      if (mslv_nb(p,nb) > 0) then
        boundleftnb  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! left  boundary of neighbor nb
        boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back  boundary of neighbor nb
        if (anb(p,nb)%x < boundleftnb) then
          anb(p,nb)%x = anb(p,nb)%x + l(1)
        endif
        if (anb(p,nb)%y > boundbacknb) then
          anb(p,nb)%y = anb(p,nb)%y - l(2)
        endif
      endif
    enddo
    nb=3
    !$omp do
    do p=1,pmax_nb(nb)
      if (mslv_nb(p,nb) > 0) then
        boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back  boundary of neighbor nb
        if (anb(p,nb)%y > boundbacknb) then
          anb(p,nb)%y = anb(p,nb)%y - l(2)
        endif
      endif
    enddo
    nb=4
    !$omp do
    do p=1,pmax_nb(nb)
      if (mslv_nb(p,nb) > 0) then
        boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
        boundbacknb  = (coords(2))*l(2)/(1.0_rp*dims(2)) ! back  boundary of neighbor nb
        if (anb(p,nb)%x > boundrightnb) then
          anb(p,nb)%x = anb(p,nb)%x - l(1)
        endif
        if (anb(p,nb)%y > boundbacknb) then
          anb(p,nb)%y = anb(p,nb)%y - l(2)
        endif
      endif
    enddo
    nb=5
    !$omp do
    do p=1,pmax_nb(nb)
      if (mslv_nb(p,nb) > 0) then
        boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
        if (anb(p,nb)%x > boundrightnb) then
          anb(p,nb)%x = anb(p,nb)%x - l(1)
        endif
      endif
    enddo
    nb=6
    !$omp do
    do p=1,pmax_nb(nb)
      if (mslv_nb(p,nb) > 0) then
        boundrightnb = (coords(1))*l(1)/(1.0_rp*dims(1)) ! right boundary of neighbor nb
        boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
        if (anb(p,nb)%x > boundrightnb) then
          anb(p,nb)%x = anb(p,nb)%x - l(1)
        endif
        if (anb(p,nb)%y < boundfrontnb) then
          anb(p,nb)%y = anb(p,nb)%y + l(2)
        endif
      endif
    enddo
    nb=7
    !$omp do
    do p=1,pmax_nb(nb)
      if (mslv_nb(p,nb) > 0) then
        boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
        if (anb(p,nb)%y < boundfrontnb) then
          anb(p,nb)%y = anb(p,nb)%y + l(2)
        endif
      endif
    enddo
    nb=8
    !$omp do
    do p=1,pmax_nb(nb)
      if (mslv_nb(p,nb) > 0) then
        boundleftnb  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! left  boundary of neighbor nb
        boundfrontnb = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! front boundary of neighbor nb
        if (anb(p,nb)%x < boundleftnb) then
          anb(p,nb)%x = anb(p,nb)%x + l(1)
        endif
        if (anb(p,nb)%y < boundfrontnb) then
          anb(p,nb)%y = anb(p,nb)%y + l(2)
        endif
      endif
    enddo
    !$omp end parallel
    !
    ! important info that should be kept:
    !
    !$omp workshare
    tp(1:pmax)%mslv    = ep(1:pmax)%mslv
    tp(1:pmax)%x       = ep(1:pmax)%x
    tp(1:pmax)%y       = ep(1:pmax)%y
    tp(1:pmax)%z       = ep(1:pmax)%z
    tp(1:pmax)%theta   = ep(1:pmax)%theta
    tp(1:pmax)%phi     = ep(1:pmax)%phi
    tp(1:pmax)%u       = ep(1:pmax)%u
    tp(1:pmax)%v       = ep(1:pmax)%v
    tp(1:pmax)%w       = ep(1:pmax)%w
    tp(1:pmax)%omx     = ep(1:pmax)%omx
    tp(1:pmax)%omy     = ep(1:pmax)%omy
    tp(1:pmax)%omz     = ep(1:pmax)%omz
    tp(1:pmax)%omtheta = ep(1:pmax)%omtheta
    tp(1:pmax)%intu    = ep(1:pmax)%intu
    tp(1:pmax)%intv    = ep(1:pmax)%intv
    tp(1:pmax)%intw    = ep(1:pmax)%intw
    tp(1:pmax)%intomx  = ep(1:pmax)%intomx
    tp(1:pmax)%intomy  = ep(1:pmax)%intomy
    tp(1:pmax)%intomz  = ep(1:pmax)%intomz
    tp(1:pmax)%colfx   = ep(1:pmax)%colfx
    tp(1:pmax)%colfy   = ep(1:pmax)%colfy
    tp(1:pmax)%colfz   = ep(1:pmax)%colfz
    tp(1:pmax)%coltx   = ep(1:pmax)%coltx
    tp(1:pmax)%colty   = ep(1:pmax)%colty
    tp(1:pmax)%coltz   = ep(1:pmax)%coltz
    tp(1:pmax)%qmax    = ep(1:pmax)%qmax
    do q=1,nqmax
      tp(1:pmax)%dx(q)     = ep(1:pmax)%dx(q)
      tp(1:pmax)%dy(q)     = ep(1:pmax)%dy(q)
      tp(1:pmax)%dz(q)     = ep(1:pmax)%dz(q)
      tp(1:pmax)%dxt(q)    = ep(1:pmax)%dxt(q)
      tp(1:pmax)%dyt(q)    = ep(1:pmax)%dyt(q)
      tp(1:pmax)%dzt(q)    = ep(1:pmax)%dzt(q)
      tp(1:pmax)%firstc(q) = ep(1:pmax)%firstc(q)
    enddo
    !
    ! clear structure ap for re-ordering:
    !
    ep(1:npmax)%mslv = 0
    ep(1:npmax)%x = 0.0_rp
    ep(1:npmax)%y = 0.0_rp
    ep(1:npmax)%z = 0.0_rp
    ep(1:npmax)%theta = 0.0_rp
    ep(1:npmax)%phi = 0.0_rp
    ep(1:npmax)%u = 0.0_rp
    ep(1:npmax)%v = 0.0_rp
    ep(1:npmax)%w = 0.0_rp
    ep(1:npmax)%omx = 0.0_rp
    ep(1:npmax)%omy = 0.0_rp
    ep(1:npmax)%omz = 0.0_rp
    ep(1:npmax)%omtheta = 0.0_rp
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
      ep(p)%nb(:) = 0
    enddo
    !$omp end workshare
    !
    leftbound   = (coords(1)  )*l(1)/(1.0_rp*dims(1)) ! left  boundary of process myid
    rightbound  = (coords(1)+1)*l(1)/(1.0_rp*dims(1)) ! right boundary of process myid
    frontbound  = (coords(2)  )*l(2)/(1.0_rp*dims(2)) ! front boundary of process myid
    backbound   = (coords(2)+1)*l(2)/(1.0_rp*dims(2)) ! back  boundary of process myid
    !
    i = 0
    count_mstr = 0
    count_slve = 0
    do idp = 1,np
      found_mstr = .false.
      call binsearch(idp,mslv_nb(1:pmax_nb(0),0),pmax_nb(0),found_mstr,p)
      !  do p=1,pmax
      !    if(sp(p)%mslv.eq.idp) then
      !      found_mstr = .true.
      if(found_mstr) then
        i = i + 1
        ep(i)%mslv    = tp(p)%mslv
        ep(i)%x       = tp(p)%x
        ep(i)%y       = tp(p)%y
        ep(i)%z       = tp(p)%z
        ep(i)%theta   = tp(p)%theta
        ep(i)%phi     = tp(p)%phi
        ep(i)%u       = tp(p)%u
        ep(i)%v       = tp(p)%v
        ep(i)%w       = tp(p)%w
        ep(i)%omx     = tp(p)%omx
        ep(i)%omy     = tp(p)%omy
        ep(i)%omz     = tp(p)%omz
        ep(i)%omtheta = tp(p)%omtheta
        ep(i)%intu    = tp(p)%intu
        ep(i)%intv    = tp(p)%intv
        ep(i)%intw    = tp(p)%intw
        ep(i)%intomx  = tp(p)%intomx
        ep(i)%intomy  = tp(p)%intomy
        ep(i)%intomz  = tp(p)%intomz
        ep(i)%colfx   = tp(p)%colfx
        ep(i)%colfy   = tp(p)%colfy
        ep(i)%colfz   = tp(p)%colfz
        ep(i)%coltx   = tp(p)%coltx
        ep(i)%colty   = tp(p)%colty
        ep(i)%coltz   = tp(p)%coltz
        ep(i)%qmax    = tp(p)%qmax
        ep(i)%dx(1:nqmax)     = tp(p)%dx(1:nqmax)
        ep(i)%dy(1:nqmax)     = tp(p)%dy(1:nqmax)
        ep(i)%dz(1:nqmax)     = tp(p)%dz(1:nqmax)
        ep(i)%dxt(1:nqmax)    = tp(p)%dxt(1:nqmax)
        ep(i)%dyt(1:nqmax)    = tp(p)%dyt(1:nqmax)
        ep(i)%dzt(1:nqmax)    = tp(p)%dzt(1:nqmax)
        ep(i)%firstc(1:nqmax) = tp(p)%firstc(1:nqmax)
        count_mstr = count_mstr + 1
        ! neighbor 1
        if ( tp(p)%x > (rightbound-(radius+offset)) ) then
          ep(i)%nb(1) = 1 ! neighbor 1 is slave of particle ap(i)%mslv
        endif
        ! neighbor 2
        dist = sqrt( (rightbound-tp(p)%x)**2 + (frontbound-tp(p)%y)**2 )
        if ( abs(dist) < (radius+offset) ) then
          ep(i)%nb(2) = 1 ! neighbor 2 is slave of particle ap(i)%mslv
        endif
        ! neighbor 3
        if ( tp(p)%y < (frontbound+(radius+offset)) ) then
          ep(i)%nb(3) = 1 ! neighbor 3 is slave of particle ap(i)%mslv
        endif
        ! neighbor 4
        dist = sqrt( (leftbound-tp(p)%x)**2 + (frontbound-tp(p)%y)**2 )
        if ( abs(dist) < (radius+offset) ) then
          ep(i)%nb(4) = 1 ! neighbor 4 is slave of particle ap(i)%mslv
        endif
        ! neighbor 5
        if ( tp(p)%x < (leftbound+(radius+offset)) ) then
          ep(i)%nb(5) = 1 ! neighbor 5 is slave of particle ap(i)%mslv
        endif
        ! neighbor 6
        dist = sqrt( (leftbound-tp(p)%x)**2 + (backbound-tp(p)%y)**2 )
        if ( abs(dist) < (radius+offset) ) then
          ep(i)%nb(6) = 1 ! neighbor 6 is slave of particle ap(i)%mslv
        endif
        ! neighbor 7
        if ( tp(p)%y > (backbound-(radius+offset)) ) then
          ep(i)%nb(7) = 1 ! neighbor 7 is slave of particle ap(i)%mslv
        endif
        ! neighbor 8
        dist = sqrt( (rightbound-tp(p)%x)**2 + (backbound-tp(p)%y)**2 )
        if ( abs(dist) < (radius+offset) ) then
          ep(i)%nb(8) = 1 ! neighbor 8 is slave of particle ap(i)%mslv
        endif
        !      call sumrk3(i,p)
        !    endif
        !  enddo
        ! if(.NOT.found_mstr) then
      else
        count_slve_loc = 0
        nb=5
        !    do k=1,pmax_nb(nb)
        !      if ( (mslv_nb(k,nb).eq.idp) ) then
        ! neighbor nb is master of particle mslv_nb(p,nb)
        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
          dist = sqrt( (leftbound-anb(k,nb)%x)**2 )
          if ( dist < (radius+offset) ) then
            if(count_slve_loc == 0 ) i = i+1
            ep(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
            ep(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            isperiodx = 0.0_rp
            isperiody = 0.0_rp
            if (anb(k,nb)%x < 0.0_rp) isperiodx =  1.0_rp
            if (anb(k,nb)%x > l(1))   isperiodx = -1.0_rp
            if (anb(k,nb)%y < 0.0_rp) isperiody =  1.0_rp
            if (anb(k,nb)%y > l(2))   isperiody = -1.0_rp
            ep(i)%x = anb(k,nb)%x+isperiodx*l(1)
            ep(i)%y = anb(k,nb)%y+isperiody*l(2)
            ep(i)%z = anb(k,nb)%z
            ep(i)%theta = anb(k,nb)%theta
            ep(i)%phi   = anb(k,nb)%phi
            ep(i)%u     = anb(k,nb)%u
            ep(i)%v     = anb(k,nb)%v
            ep(i)%w     = anb(k,nb)%w
            ep(i)%omx   = anb(k,nb)%omx
            ep(i)%omy   = anb(k,nb)%omy
            ep(i)%omz   = anb(k,nb)%omz
          endif
        endif
        !    enddo
        nb=6
        !    do k=1,pmax_nb(nb)
        !      if ( (mslv_nb(k,nb).eq.idp) ) then
        ! neighbor nb is master of particle mslv_nb(p,nb)
        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
          dist = sqrt( (leftbound-anb(k,nb)%x)**2 + (backbound-anb(k,nb)%y)**2 )
          if ( dist < (radius+offset) ) then
            if(count_slve_loc == 0 ) i = i+1
            ep(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
            ep(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            isperiodx = 0.0_rp
            isperiody = 0.0_rp
            if (anb(k,nb)%x < 0.0_rp) isperiodx =  1.0_rp
            if (anb(k,nb)%x > l(1))   isperiodx = -1.0_rp
            if (anb(k,nb)%y < 0.0_rp) isperiody =  1.0_rp
            if (anb(k,nb)%y > l(2))   isperiody = -1.0_rp
            ep(i)%x = anb(k,nb)%x+isperiodx*l(1)
            ep(i)%y = anb(k,nb)%y+isperiody*l(2)
            ep(i)%z = anb(k,nb)%z
            ep(i)%theta = anb(k,nb)%theta
            ep(i)%phi   = anb(k,nb)%phi
            ep(i)%u     = anb(k,nb)%u
            ep(i)%v     = anb(k,nb)%v
            ep(i)%w     = anb(k,nb)%w
            ep(i)%omx   = anb(k,nb)%omx
            ep(i)%omy   = anb(k,nb)%omy
            ep(i)%omz   = anb(k,nb)%omz
          endif
        endif
        !    enddo
        nb=7
        !    do k=1,pmax_nb(nb)
        !      if ( (mslv_nb(k,nb).eq.idp) ) then
        ! neighbor nb is master of particle mslv_nb(p,nb)
        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
          dist = sqrt( (backbound-anb(k,nb)%y)**2 )
          if ( dist < (radius+offset) ) then
            if(count_slve_loc == 0 ) i = i+1
            ep(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
            ep(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            isperiodx = 0.0_rp
            isperiody = 0.0_rp
            if (anb(k,nb)%x < 0.0_rp) isperiodx =  1.0_rp
            if (anb(k,nb)%x > l(1))   isperiodx = -1.0_rp
            if (anb(k,nb)%y < 0.0_rp) isperiody =  1.0_rp
            if (anb(k,nb)%y > l(2))   isperiody = -1.0_rp
            ep(i)%x = anb(k,nb)%x+isperiodx*l(1)
            ep(i)%y = anb(k,nb)%y+isperiody*l(2)
            ep(i)%z = anb(k,nb)%z
            ep(i)%theta = anb(k,nb)%theta
            ep(i)%phi   = anb(k,nb)%phi
            ep(i)%u     = anb(k,nb)%u
            ep(i)%v     = anb(k,nb)%v
            ep(i)%w     = anb(k,nb)%w
            ep(i)%omx   = anb(k,nb)%omx
            ep(i)%omy   = anb(k,nb)%omy
            ep(i)%omz   = anb(k,nb)%omz
          endif
        endif
        !    enddo
        nb=8
        !    do k=1,pmax_nb(nb)
        !      if ( (mslv_nb(k,nb).eq.idp) ) then
        !neighbor nb is master of particle mslv_nb(p,nb)
        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
          dist = sqrt( (rightbound-anb(k,nb)%x)**2 + (backbound-anb(k,nb)%y)**2 )
          if ( dist < (radius+offset) ) then
            if(count_slve_loc == 0 ) i = i+1
            ep(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
            ep(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            isperiodx = 0.0_rp
            isperiody = 0.0_rp
            if (anb(k,nb)%x < 0.0_rp) isperiodx =  1.0_rp
            if (anb(k,nb)%x > l(1))   isperiodx = -1.0_rp
            if (anb(k,nb)%y < 0.0_rp) isperiody =  1.0_rp
            if (anb(k,nb)%y > l(2))   isperiody = -1.0_rp
            ep(i)%x = anb(k,nb)%x+isperiodx*l(1)
            ep(i)%y = anb(k,nb)%y+isperiody*l(2)
            ep(i)%z = anb(k,nb)%z
            ep(i)%theta = anb(k,nb)%theta
            ep(i)%phi   = anb(k,nb)%phi
            ep(i)%u     = anb(k,nb)%u
            ep(i)%v     = anb(k,nb)%v
            ep(i)%w     = anb(k,nb)%w
            ep(i)%omx   = anb(k,nb)%omx
            ep(i)%omy   = anb(k,nb)%omy
            ep(i)%omz   = anb(k,nb)%omz
          endif
        endif
        !    enddo
        nb=1
        !    do k=1,pmax_nb(nb)
        !      if ( (mslv_nb(k,nb).eq.idp) ) then
        ! neighbor nb is master of particle mslv_nb(p,nb)
        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
          dist = sqrt( (rightbound-anb(k,nb)%x)**2 )
          if ( dist < (radius+offset) ) then
            if(count_slve_loc == 0 ) i = i+1
            ep(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
            ep(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            isperiodx = 0.0_rp
            isperiody = 0.0_rp
            if (anb(k,nb)%x < 0.0_rp) isperiodx =  1.0_rp
            if (anb(k,nb)%x > l(1))   isperiodx = -1.0_rp
            if (anb(k,nb)%y < 0.0_rp) isperiody =  1.0_rp
            if (anb(k,nb)%y > l(2))   isperiody = -1.0_rp
            ep(i)%x = anb(k,nb)%x+isperiodx*l(1)
            ep(i)%y = anb(k,nb)%y+isperiody*l(2)
            ep(i)%z = anb(k,nb)%z
            ep(i)%theta = anb(k,nb)%theta
            ep(i)%phi   = anb(k,nb)%phi
            ep(i)%u     = anb(k,nb)%u
            ep(i)%v     = anb(k,nb)%v
            ep(i)%w     = anb(k,nb)%w
            ep(i)%omx   = anb(k,nb)%omx
            ep(i)%omy   = anb(k,nb)%omy
            ep(i)%omz   = anb(k,nb)%omz
          endif
        endif
        !    enddo
        nb=2
        !    do k=1,pmax_nb(nb)
        !      if ( (mslv_nb(k,nb).eq.idp) ) then
        ! neighbor nb is master of particle mslv_nb(p,nb)
        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
          dist = sqrt( (rightbound-anb(k,nb)%x)**2 + (frontbound-anb(k,nb)%y)**2 )
          if ( dist < (radius+offset) ) then
            if(count_slve_loc == 0 ) i = i+1
            ep(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
            ep(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            isperiodx = 0.0_rp
            isperiody = 0.0_rp
            if (anb(k,nb)%x < 0.0_rp) isperiodx =  1.0_rp
            if (anb(k,nb)%x > l(1))   isperiodx = -1.0_rp
            if (anb(k,nb)%y < 0.0_rp) isperiody =  1.0_rp
            if (anb(k,nb)%y > l(2))   isperiody = -1.0_rp
            ep(i)%x = anb(k,nb)%x+isperiodx*l(1)
            ep(i)%y = anb(k,nb)%y+isperiody*l(2)
            ep(i)%z = anb(k,nb)%z
            ep(i)%theta = anb(k,nb)%theta
            ep(i)%phi   = anb(k,nb)%phi
            ep(i)%u     = anb(k,nb)%u
            ep(i)%v     = anb(k,nb)%v
            ep(i)%w     = anb(k,nb)%w
            ep(i)%omx   = anb(k,nb)%omx
            ep(i)%omy   = anb(k,nb)%omy
            ep(i)%omz   = anb(k,nb)%omz
          endif
        endif
        !    enddo
        nb=3
        !    do k=1,pmax_nb(nb)
        !      if ( (mslv_nb(k,nb).eq.idp) ) then
        ! neighbor nb is master of particle mslv_nb(p,nb)
        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
          dist = sqrt( (frontbound-anb(k,nb)%y)**2 )
          if ( dist < (radius+offset) ) then
            if(count_slve_loc == 0 ) i = i+1
            ep(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
            ep(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            isperiodx = 0.0_rp
            isperiody = 0.0_rp
            if (anb(k,nb)%x < 0.0_rp) isperiodx =  1.0_rp
            if (anb(k,nb)%x > l(1))   isperiodx = -1.0_rp
            if (anb(k,nb)%y < 0.0_rp) isperiody =  1.0_rp
            if (anb(k,nb)%y > l(2))   isperiody = -1.0_rp
            ep(i)%x = anb(k,nb)%x+isperiodx*l(1)
            ep(i)%y = anb(k,nb)%y+isperiody*l(2)
            ep(i)%z = anb(k,nb)%z
            ep(i)%theta = anb(k,nb)%theta
            ep(i)%phi   = anb(k,nb)%phi
            ep(i)%u     = anb(k,nb)%u
            ep(i)%v     = anb(k,nb)%v
            ep(i)%w     = anb(k,nb)%w
            ep(i)%omx   = anb(k,nb)%omx
            ep(i)%omy   = anb(k,nb)%omy
            ep(i)%omz   = anb(k,nb)%omz
          endif
        endif
        !    enddo
        nb=4
        !    do k=1,pmax_nb(nb)
        !      if ( (mslv_nb(k,nb).eq.idp) ) then
        ! neighbor nb is master of particle mslv_nb(p,nb)
        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
          dist = sqrt( (leftbound-anb(k,nb)%x)**2 + (frontbound-anb(k,nb)%y)**2 )
          if ( dist < (radius+offset) ) then
            if(count_slve_loc == 0 ) i = i+1
            ep(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
            ep(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            isperiodx = 0.0_rp
            isperiody = 0.0_rp
            if (anb(k,nb)%x < 0.0_rp) isperiodx =  1.0_rp
            if (anb(k,nb)%x > l(1))   isperiodx = -1.0_rp
            if (anb(k,nb)%y < 0.0_rp) isperiody =  1.0_rp
            if (anb(k,nb)%y > l(2))   isperiody = -1.0_rp
            ep(i)%x = anb(k,nb)%x+isperiodx*l(1)
            ep(i)%y = anb(k,nb)%y+isperiody*l(2)
            ep(i)%z = anb(k,nb)%z
            ep(i)%theta = anb(k,nb)%theta
            ep(i)%phi   = anb(k,nb)%phi
            ep(i)%u     = anb(k,nb)%u
            ep(i)%v     = anb(k,nb)%v
            ep(i)%w     = anb(k,nb)%w
            ep(i)%omx   = anb(k,nb)%omx
            ep(i)%omy   = anb(k,nb)%omy
            ep(i)%omz   = anb(k,nb)%omz
          endif
        endif
        !    enddo
        if(count_slve_loc.ne.0) count_slve = count_slve + 1
      endif
    enddo
    !
    ! the new value of pmax is:
    !
    pmax = count_mstr + count_slve

    npmstr = count_mstr
    !
    !check if number of masters yield np
    !
    call MPI_ALLREDUCE(count_mstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
!!!!    if(myid == 0) write(6,*) 'num of particles = num of masters?',np,' = ', count_mstr_all
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
    endif
    !
    !write new master/slave configuration to a file for debugging purposes
    !
    !write(rankpr,'(i3.3)') myid
    !open(25,file=datadir//'mslv'//rankpr//'.txt')
    !do p=1,pmax
    !  idp = abs(ap(p)%mslv)
    !  if (ap(p)%mslv .gt. 0) then
    !    counter = 0
    !    do i=1,8
    !      if (ap(p)%nb(i) .eq. 1) then
    !        counter = counter+1
    !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !!        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !!              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !      endif
    !    enddo
    !    !in case of no overlap with any neighbor
    !    if (counter .eq. 0) then
    !      write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
    !            myid,' ',idp,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
    !!      write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !!              myid,' ',idp,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
    !    endif
    !  endif
    !  if (ap(p)%mslv .lt. 0) then
    !    do i=1,8
    !      if (ap(p)%nb(i) .eq. 1) then
    !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !!        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !!              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !      endif
    !    enddo
    !  endif
    !enddo
    !close(25)
    !
    ! masters: new positions and velocities of Lagrangian forcing points
    !
#if !defined(_EULER)
    !$omp workshare
    nla(:) = 0 ! set to zero from 1 to npmax
    !$omp end workshare
!
    !$omp parallel default(none)               &
    !$omp shared(ap,nla,pmax)                  &
    !$omp shared(radfp,phirc,thetarc)          &
    !$omp shared(boundleftmyid,boundfrontmyid) &
    !$omp private(p,l,ll,coorxfp,cooryfp,xfploc,yfploc,zfploc)    &
    !$omp private(isperiodx,isperiody,isout)
    !$omp do 
    do p=1,pmax
      if (ep(p)%mslv /= 0) then
        ! myid is master of particle ep(p)%mslv
        ll = 0
        do lp=1,NL
          xfploc = ep(p)%x + radfp*sin(ep(p)%theta*0.0_rp+thetarc(lp))*cos(ep(p)%phi*0.0_rp+phirc(lp))
          yfploc = ep(p)%y + radfp*sin(ep(p)%theta*0.0_rp+thetarc(lp))*sin(ep(p)%phi*0.0_rp+phirc(lp))
          zfploc = ep(p)%z + radfp*cos(ep(p)%theta*0.0_rp+thetarc(lp))
          isperiodx = 0.0_rp
          isperiody = 0.0_rp
          if (xfploc < 0.0_rp+0.5_rp*dl(1)) isperiodx =  1.0_rp
          if (xfploc > l(1)+0.5_rp*dl(1))   isperiodx = -1.0_rp
          if (yfploc < 0.0_rp+0.5_rp*dl(2)) isperiody =  1.0_rp
          if (yfploc > l(2)+0.5_rp*dl(2))   isperiody = -1.0_rp
          isout = .false.
          coorxfp = (xfploc+isperiodx*l(1)-boundleftmyid )*dli(1)
          if( nint(coorxfp) < 1 .or. nint(coorxfp)  > n(1) ) isout = .true.
          cooryfp = (yfploc+isperiody*l(2)-boundfrontmyid)*dli(2)
          if( nint(cooryfp) < 1 .or. nint(cooryfp) > n(2) ) isout = .true.
!
          if(.not.isout) then
            ll = ll + 1
            ep(p)%xfp(ll) = xfploc
            ep(p)%yfp(ll) = yfploc
            ep(p)%zfp(ll) = zfploc
            ep(p)%ul(ll) = ep(p)%u + ep(p)%omy*(ep(p)%zfp(ll)-ep(p)%z) &
                                   - ep(p)%omz*(ep(p)%yfp(ll)-ep(p)%y)
            ep(p)%vl(ll) = ep(p)%v + ep(p)%omz*(ep(p)%xfp(ll)-ep(p)%x) &
                                   - ep(p)%omx*(ep(p)%zfp(ll)-ep(p)%z)
            ep(p)%wl(ll) = ep(p)%w + ep(p)%omx*(ep(p)%yfp(ll)-ep(p)%y) &
                                   - ep(p)%omy*(ep(p)%xfp(ll)-ep(p)%x)
          endif
        enddo
        nla(p) = ll
      endif
    enddo
    !$omp end parallel
#endif
    !
    return
  end subroutine intgr_nwtn_eulr
  !


  !
  subroutine binsearch(key,array,idmax,found,indx) 
    !! this is not a real binary search !!
    ! it searches key with absolute value, then checks its positivity
    integer, intent(in), dimension(1:) :: array
    integer, intent(in) :: key,idmax
    logical, intent(out) :: found
    integer, intent(out) :: indx
    integer imin, imax, imid

    found = .false.
    indx = 0
    imin = 1
    imax = idmax

    ! continue searching while [imin,imax] is not empty
    do while (imin <= imax)
      !calculate the midpoint for roughly equal partition
      imid = (imin + imax)/2;
      if(abs(array(imid)) == key) then
        !key found at index imid 
        if (array(imid) == key) then
          !! this is what makes it not a real binary search !!
          indx = imid
          found = .true.
        endif
        return
        ! determine which subarray to search
      else if (abs(array(imid)) < key) then
        !change min index to search upper subarray
        imin = imid + 1;
      else          
        !change max index to search lower subarray
        imax = imid - 1;
      end if
    end do
    return
  end subroutine binsearch
  !

  !
#endif
end module prt_mod_intgr_nwtn_eulr
