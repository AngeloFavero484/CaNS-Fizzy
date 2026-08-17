module mod_forces
  !
  use mod_param     , only: nl,rp,visci,datadir,cx,cy,cz
  use mod_common    , only: lpx,lpy,lpz,lptheta,lpphi,nx,ny,nz,lpA,rad
  use mod_interp    , only: InterpFastVel,InterpFastSca,InterpFastGRP
  !
  implicit none
  !
  real(rp) :: Lprb
  !
!  integer :: N_FCS
!  !
!  type stl
!    real(rp) :: c0_x
!    real(rp) :: c0_y
!    real(rp) :: c0_z
!    real(rp) :: c1_x
!    real(rp) :: c1_y
!    real(rp) :: c1_z
!    real(rp) :: c2_x
!    real(rp) :: c2_y
!    real(rp) :: c2_z
!    real(rp) :: c3_x
!    real(rp) :: c3_y
!    real(rp) :: c3_z
!    real(rp) :: A
!    real(rp) :: f_x
!    real(rp) :: f_y
!    real(rp) :: f_z
!    real(rp) :: n_x
!    real(rp) :: n_y
!    real(rp) :: n_z
!  end type stl
!  !
!  type(stl),dimension(:),allocatable :: geo
  !
  private
  !
!  public :: InitForces
  public :: CmptDNSForces
!  public :: CmptWMLESForces
!  public :: ForcesSetCyl
!  public :: ForcesSetChn
!  public :: ForcesSetPipe
  !
  contains
  !
  !===================================================================================================================
  !
!  subroutine InitForces
!    implicit none
!    integer :: fid,err,dum
!    integer :: f
!    real(rp),dimension(1:3) :: vecA,vecB,vecC
!    real(rp),dimension(1:3) :: norm
!    !
!    open(newunit=fid,file='./data/geom.dat',status='old',action='read',iostat=err)
!      !
!      if(err==0) then
!        !
!        read(fid,*,iostat=err) N_FCS
!        !
!        allocate(geo(1:N_FCS))
!        !
!        do f=1,N_FCS
!          !
!          read(fid,*,iostat=err) dum,geo(f)%c1_x,geo(f)%c1_y,geo(f)%c1_z,&
!                                     geo(f)%c2_x,geo(f)%c2_y,geo(f)%c2_z,&
!                                     geo(f)%c3_x,geo(f)%c3_y,geo(f)%c3_z
!          !
!          geo(f)%c0_x=1.0_rp/3.0_rp*(geo(f)%c1_x+geo(f)%c2_x+geo(f)%c3_x)
!          geo(f)%c0_y=1.0_rp/3.0_rp*(geo(f)%c1_y+geo(f)%c2_y+geo(f)%c3_y)
!          geo(f)%c0_z=1.0_rp/3.0_rp*(geo(f)%c1_z+geo(f)%c2_z+geo(f)%c3_z)
!          !
!          vecA(1)=geo(f)%c2_x-geo(f)%c1_x
!          vecA(2)=geo(f)%c2_y-geo(f)%c1_y
!          vecA(3)=geo(f)%c2_z-geo(f)%c1_z
!          !
!          vecB(1)=geo(f)%c3_x-geo(f)%c1_x
!          vecB(2)=geo(f)%c3_y-geo(f)%c1_y
!          vecB(3)=geo(f)%c3_z-geo(f)%c1_z
!          !
!          vecC=cross(vecA,vecB)
!          !
!          geo(f)%A=0.5_rp*norm_2(vecC)
!          !
!        enddo
!        !
!      else
!        !
!        if(myid==0) print*, 'Error reading the geometry input file.'
!        if(myid==0) print*, 'Aborting ...'
!        !
!        call MPI_FINALIZE(err)
!        !
!        error stop
!        !
!      end if
!      !
!    close(fid)
!    !
!    do f=1,N_FCS
!      !
!      vecA(1)=geo(f)%c2_x-geo(f)%c1_x
!      vecA(2)=geo(f)%c2_y-geo(f)%c1_y
!      vecA(3)=geo(f)%c2_z-geo(f)%c1_z
!      !
!      vecB(1)=geo(f)%c3_x-geo(f)%c1_x
!      vecB(2)=geo(f)%c3_y-geo(f)%c1_y
!      vecB(3)=geo(f)%c3_z-geo(f)%c1_z
!      !
!      vecC=cross(vecA,vecB)
!      !
!      norm=vecC/norm_2(vecC)
!      !
!      geo(f)%n_x=norm(1)
!      geo(f)%n_y=norm(2)
!      geo(f)%n_z=norm(3)
!      !
!    enddo
!    !
!    return
!  end subroutine InitForces
  !
  !===================================================================================================================
  !
  subroutine CmptDNSForces(ng,dl,u,v,w,p,f_x,f_y,f_z)
    implicit none
    integer,dimension(1:3),intent(in) :: ng
    real(rp),dimension(1:3),intent(in) :: dl
!    real(rp),dimension(-3:n(3)+4),intent(in) :: zc,dzc
!    real(rp),dimension(-3:n(3)+4),intent(in) :: zf,dzf
!    real(rp),intent(in) :: time
    real(rp),dimension(1:ng(1),1:ng(2),1:ng(3)),intent(in) :: u,v,w,p
    real(rp),intent(out) :: f_x,f_y,f_z
    !
    integer :: idx
    integer :: err
    integer :: i,j,k,ll
    real(rp) :: stn_c_x,stn_c_y,stn_c_z
    real(rp) :: stn_t_x,stn_t_y,stn_t_z
    real(rp) :: stn_s_x,stn_s_y,stn_s_z
    real(rp) :: stn_n_x,stn_n_y,stn_n_z
    real(rp) :: stn_e_x,stn_e_y,stn_e_z
    real(rp) :: stn_w_x,stn_w_y,stn_w_z
    real(rp) :: vel_c_x,vel_c_y,vel_c_z,vel_c_t,vel_c_b
    real(rp) :: vel_t_x,vel_t_y,vel_t_z,vel_t_t,vel_t_b
    real(rp) :: vel_s_x,vel_s_y,vel_s_z,vel_s_n
    real(rp) :: vel_n_x,vel_n_y,vel_n_z,vel_n_n
    real(rp) :: vel_e_x,vel_e_y,vel_e_z,vel_e_n
    real(rp) :: vel_w_x,vel_w_y,vel_w_z,vel_w_n
    real(rp) :: vel_r_t,vel_r_b,vel_c_n
    real(rp) :: pre_t
    real(rp) :: dp_dx,dp_dy,dp_dz
    real(rp) :: dun_dt,dut_dn
    real(rp) :: dun_db,dub_dn
    real(rp) :: dp_dn_m,dp_dn_p
    real(rp) :: tau_t,tau_n,tau_b
    real(rp) :: tau_x,tau_y,tau_z
    real(rp) :: delta,fct
    real(rp),dimension(1:3) :: norm,tang,binr
    real(rp),dimension(1:ng(1),1:ng(2),1:ng(3)) :: vel_x,vel_y,vel_z,sca,gp_x,gp_y,gp_z
    !
    fct=1.99_rp
    !
    vel_x=0.0_rp
    vel_y=0.0_rp
    vel_z=0.0_rp
    !
    sca=0.0_rp
    !
    gp_x=0.0_rp
    gp_y=0.0_rp
    gp_z=0.0_rp
    !
    vel_x(1:ng(1),1:ng(2),1:ng(3))=u(1:ng(1),1:ng(2),1:ng(3))
    vel_y(1:ng(1),1:ng(2),1:ng(3))=v(1:ng(1),1:ng(2),1:ng(3))
    vel_z(1:ng(1),1:ng(2),1:ng(3))=w(1:ng(1),1:ng(2),1:ng(3))
    !
    sca(1:ng(1),1:ng(2),1:ng(3))=p(1:ng(1),1:ng(2),1:ng(3))
    !
!    call MPI_UpdtHalosVel(n,4,nb,MPI_halo4,vel_x,vel_y,vel_z)
!    call MPI_UpdtHalosSca(n,4,nb,MPI_halo4,sca)
    !
    do k=1,ng(3)-1
      do j=1,ng(2)-1
        do i=1,ng(1)-1
          !
          gp_x(i,j,k)=(sca(i+1,j  ,k  )-sca(i,j,k))/dl(1)
          gp_y(i,j,k)=(sca(i  ,j+1,k  )-sca(i,j,k))/dl(2)
          gp_z(i,j,k)=(sca(i  ,j  ,k+1)-sca(i,j,k))/dl(3)
          !
        enddo
      enddo
    enddo
    !
    f_x=0.0_rp
    f_y=0.0_rp
    f_z=0.0_rp
    !
    open(42,file=trim(datadir)//'pressure.dat')
    !
    do ll=1,nl
      !
      norm(1)=nx(ll)
      norm(2)=ny(ll)
      norm(3)=nz(ll)
      !
      idx=minloc(abs(norm),1)
      !  
      if(idx==1) tang=cross(norm,[1.0_rp,0.0_rp,0.0_rp])
      if(idx==2) tang=cross(norm,[0.0_rp,1.0_rp,0.0_rp])
      if(idx==3) tang=cross(norm,[0.0_rp,0.0_rp,1.0_rp])
      !
      tang=tang/norm_2(tang)
      !
      binr=cross(tang,norm)
      binr=binr/norm_2(binr)
      !
      delta=abs(norm(1)*dl(1))+abs(norm(2)*dl(2))+abs(norm(3)*dl(3))
      !
      Lprb=fct*delta
      !       
      stn_t_x=lpx(ll)+norm(1)*Lprb
      stn_t_y=lpy(ll)+norm(2)*Lprb
      stn_t_z=lpz(ll)+norm(3)*Lprb
      !
      stn_c_x=lpx(ll)+norm(1)*0.5_rp*Lprb
      stn_c_y=lpy(ll)+norm(2)*0.5_rp*Lprb
      stn_c_z=lpz(ll)+norm(3)*0.5_rp*Lprb
      !
      stn_w_x=lpx(ll)+(norm(1)-tang(1))*0.5_rp*Lprb
      stn_w_y=lpy(ll)+(norm(2)-tang(2))*0.5_rp*Lprb
      stn_w_z=lpz(ll)+(norm(3)-tang(3))*0.5_rp*Lprb
      !
      stn_e_x=lpx(ll)+(norm(1)+tang(1))*0.5_rp*Lprb
      stn_e_y=lpy(ll)+(norm(2)+tang(2))*0.5_rp*Lprb
      stn_e_z=lpz(ll)+(norm(3)+tang(3))*0.5_rp*Lprb
      !
      stn_s_x=lpx(ll)+(norm(1)-binr(1))*0.5_rp*Lprb
      stn_s_y=lpy(ll)+(norm(2)-binr(2))*0.5_rp*Lprb
      stn_s_z=lpz(ll)+(norm(3)-binr(3))*0.5_rp*Lprb
      !
      stn_n_x=lpx(ll)+(norm(1)+binr(1))*0.5_rp*Lprb
      stn_n_y=lpy(ll)+(norm(2)+binr(2))*0.5_rp*Lprb
      stn_n_z=lpz(ll)+(norm(3)+binr(3))*0.5_rp*Lprb
      ! 
      call InterpFastVel(ng,dl,(/stn_t_x,stn_t_y,stn_t_z/),vel_x,vel_y,vel_z,vel_t_x,vel_t_y,vel_t_z)
      call InterpFastVel(ng,dl,(/stn_c_x,stn_c_y,stn_c_z/),vel_x,vel_y,vel_z,vel_c_x,vel_c_y,vel_c_z)
      call InterpFastVel(ng,dl,(/stn_w_x,stn_w_y,stn_w_z/),vel_x,vel_y,vel_z,vel_w_x,vel_w_y,vel_w_z)
      call InterpFastVel(ng,dl,(/stn_e_x,stn_e_y,stn_e_z/),vel_x,vel_y,vel_z,vel_e_x,vel_e_y,vel_e_z)
      call InterpFastVel(ng,dl,(/stn_s_x,stn_s_y,stn_s_z/),vel_x,vel_y,vel_z,vel_s_x,vel_s_y,vel_s_z)
      call InterpFastVel(ng,dl,(/stn_n_x,stn_n_y,stn_n_z/),vel_x,vel_y,vel_z,vel_n_x,vel_n_y,vel_n_z)
      !
      call InterpFastSca(ng,dl,(/stn_t_x,stn_t_y,stn_t_z/),sca,pre_t)
      call InterpFastGRP(ng,dl,(/stn_t_x,stn_t_y,stn_t_z/),gp_x,gp_y,gp_z,dp_dx,dp_dy,dp_dz)
      !
      vel_c_n=vel_c_x*norm(1)+vel_c_y*norm(2)+vel_c_z*norm(3)
      !
      vel_e_n=vel_e_x*norm(1)+vel_e_y*norm(2)+vel_e_z*norm(3)
      vel_w_n=vel_w_x*norm(1)+vel_w_y*norm(2)+vel_w_z*norm(3)
      !
      vel_n_n=vel_n_x*norm(1)+vel_n_y*norm(2)+vel_n_z*norm(3)
      vel_s_n=vel_s_x*norm(1)+vel_s_y*norm(2)+vel_s_z*norm(3)
      !
      vel_t_t=vel_t_x*tang(1)+vel_t_y*tang(2)+vel_t_z*tang(3)
      vel_c_t=vel_c_x*tang(1)+vel_c_y*tang(2)+vel_c_z*tang(3)
      !
      vel_t_b=vel_t_x*binr(1)+vel_t_y*binr(2)+vel_t_z*binr(3)
      vel_c_b=vel_c_x*binr(1)+vel_c_y*binr(2)+vel_c_z*binr(3)
      !
      vel_r_t=0.0_rp
      vel_r_b=0.0_rp
      !
      dun_dt=0.0_rp!(vel_e_n-vel_w_n)/Lprb
      !
      dun_db=0.0_rp!(vel_n_n-vel_s_n)/Lprb
      !
      dut_dn=(-3.0_rp*vel_r_t+4.0_rp*vel_c_t-vel_t_t)/Lprb
      !
      dub_dn=(-3.0_rp*vel_r_b+4.0_rp*vel_c_b-vel_t_b)/Lprb
      !
      dp_dn_m=0.0_rp
      !
      dp_dn_p=dp_dx*norm(1)+dp_dy*norm(2)+dp_dz*norm(3)
      !
      tau_t=1.0_rp/visci*(dut_dn+dun_dt)
      !
      tau_b=1.0_rp/visci*(dub_dn+dun_db)
      !
      tau_n=-(pre_t-0.5_rp*(dp_dn_m+dp_dn_p)*Lprb)
      !
      tau_x=tang(1)*tau_t+norm(1)*tau_n+binr(1)*tau_b
      tau_y=tang(2)*tau_t+norm(2)*tau_n+binr(2)*tau_b
      tau_z=tang(3)*tau_t+norm(3)*tau_n+binr(3)*tau_b
      !
      f_x=f_x+tau_x*lpA
      f_y=f_y+tau_y*lpA
      f_z=f_z+tau_z*lpA
      !
      write(42,'(7E16.8)') lpx(ll)-cx,lpy(ll)-cy,lpz(ll)-cz,lptheta(ll),lpphi(ll),rad,tau_n
      !
    enddo
    close(42)
    print*,f_x,f_y,f_z
    !
    return
  end subroutine CmptDNSForces
  !
  !===================================================================================================================
  !
!  subroutine CmptWMLESForces(n,dl,zc,dzc,zf,dzf,time,u,v,w,p,phi_x,phi_y,phi_z,f_x,f_y,f_z)
!    implicit none
!    integer,dimension(1:3),intent(in) :: n
!    real(rp),dimension(1:3),intent(in) :: dl
!    real(rp),dimension(-3:n(3)+4),intent(in) :: zc,dzc
!    real(rp),dimension(-3:n(3)+4),intent(in) :: zf,dzf
!    real(rp),intent(in) :: time
!    real(rp),dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1),intent(in) :: u,v,w,p
!    real(rp),dimension(-3:n(1)+4,-3:n(2)+4,-3:n(3)+4),intent(in) :: phi_x,phi_y,phi_z
!    real(rp),intent(out) :: f_x,f_y,f_z
!    !
!    integer :: idx,err,chk
!    integer :: i,j,k,f
!    real(rp) :: pre_t
!    real(rp) :: dp_dx,dp_dy,dp_dz
!    real(rp) :: dp_dn_m,dp_dn_p
!    real(rp) :: tau_t,tau_n,tau_b
!    real(rp) :: tau_x,tau_y,tau_z
!    real(rp) :: delta,fct
!    real(rp) :: tau_w,tau_w_t,tau_w_b,u_tau
!    real(rp) :: U_tang,U_binr,U_mod
!    real(rp) :: v_x,v_y,v_z
!    real(rp),dimension(1:3) :: point,probe,norm,tang,binr
!    real(rp),dimension(-3:n(1)+4,-3:n(2)+4,-3:n(3)+4) :: vel_x,vel_y,vel_z,sca,gp_x,gp_y,gp_z
!    !
!    fct=1.99_rp
!    !
!    vel_x=0.0_rp
!    vel_y=0.0_rp
!    vel_z=0.0_rp
!    !
!    sca=0.0_rp
!    !
!    gp_x=0.0_rp
!    gp_y=0.0_rp
!    gp_z=0.0_rp
!    !
!    vel_x(0:n(1)+1,0:n(2)+1,0:n(3)+1)=u(0:n(1)+1,0:n(2)+1,0:n(3)+1)
!    vel_y(0:n(1)+1,0:n(2)+1,0:n(3)+1)=v(0:n(1)+1,0:n(2)+1,0:n(3)+1)
!    vel_z(0:n(1)+1,0:n(2)+1,0:n(3)+1)=w(0:n(1)+1,0:n(2)+1,0:n(3)+1)
!    !
!    sca(0:n(1)+1,0:n(2)+1,0:n(3)+1)=p(0:n(1)+1,0:n(2)+1,0:n(3)+1)
!    !
!    call MPI_UpdtHalosVel(n,4,nb,MPI_halo4,vel_x,vel_y,vel_z)
!    call MPI_UpdtHalosSca(n,4,nb,MPI_halo4,sca)
!    !
!    do k=-3,n(3)+3
!      do j=-3,n(2)+3
!        do i=-3,n(1)+3
!          !
!          gp_x(i,j,k)=(sca(i+1,j  ,k  )-sca(i,j,k))/dl(1) 
!          gp_y(i,j,k)=(sca(i  ,j+1,k  )-sca(i,j,k))/dl(2) 
!          gp_z(i,j,k)=(sca(i  ,j  ,k+1)-sca(i,j,k))/dzc(k) 
!          !
!        enddo
!      enddo
!    enddo 
!    !
!    f_x=0.0_rp
!    f_y=0.0_rp
!    f_z=0.0_rp
!    !
!    do f=1,N_FCS
!      !
!      norm(1)=geo(f)%n_x
!      norm(2)=geo(f)%n_y
!      norm(3)=geo(f)%n_z
!      !
!      idx=minloc(abs(norm),1)
!      !  
!      if(idx==1) tang=cross(norm,[1.0_rp,0.0_rp,0.0_rp])
!      if(idx==2) tang=cross(norm,[0.0_rp,1.0_rp,0.0_rp])
!      if(idx==3) tang=cross(norm,[0.0_rp,0.0_rp,1.0_rp])
!      !
!      tang=tang/norm_2(tang)
!      !
!      binr=cross(tang,norm)
!      binr=binr/norm_2(binr)
!      !
!      delta=abs(norm(1)*dl(1))+abs(norm(2)*dl(2))+abs(norm(3)*dl(3))
!      !
!      Lprb=fct*delta
!      !
!      point=(/geo(f)%c0_x,geo(f)%c0_y,geo(f)%c0_z/)
!      !
!      probe=norm*Lprb
!      !
!      call InterpIBMVelX(n,dl,point+probe,zc,zf,dzc,dzf,vel_x,phi_x,v_x)
!      call InterpIBMVelY(n,dl,point+probe,zc,zf,dzc,dzf,vel_y,phi_y,v_y)
!      call InterpIBMVelZ(n,dl,point+probe,zc,zf,dzc,dzf,vel_z,phi_z,v_z)
!      !
!      U_tang=v_x*tang(1)+v_y*tang(2)+v_z*tang(3)
!      U_binr=v_x*binr(1)+v_y*binr(2)+v_z*binr(3)
!      !
!      U_mod=sqrt(U_tang**2+U_binr**2)
!      !
!      call NewtonRaphson(U_mod,Lprb,u_tau,chk)
!      !
!      tau_w=u_tau**2
!      !
!      tau_w_t=tau_w*U_tang/U_mod
!      tau_w_b=tau_w*U_binr/U_mod
!      !
!      call InterpFastSca(n,dl,zc,dzc,point+probe,sca,pre_t)
!      call InterpFastGRP(n,dl,zc,dzc,point+probe,gp_x,gp_y,gp_z,dp_dx,dp_dy,dp_dz)
!      !
!      dp_dn_m=0.0_rp
!      !
!      dp_dn_p=dp_dx*norm(1)+dp_dy*norm(2)+dp_dz*norm(3)
!      !
!      tau_t=tau_w_t
!      !
!      tau_b=tau_w_b
!      !
!      tau_n=-(pre_t-0.5_rp*(dp_dn_m+dp_dn_p)*Lprb)
!      !
!      tau_x=tang(1)*tau_t+norm(1)*tau_n+binr(1)*tau_b
!      tau_y=tang(2)*tau_t+norm(2)*tau_n+binr(2)*tau_b
!      tau_z=tang(3)*tau_t+norm(3)*tau_n+binr(3)*tau_b
!      !
!      geo(f)%f_x=tau_x*geo(f)%A
!      geo(f)%f_y=tau_y*geo(f)%A
!      geo(f)%f_z=tau_z*geo(f)%A
!      !
!      f_x=f_x+geo(f)%f_x
!      f_y=f_y+geo(f)%f_y
!      f_z=f_z+geo(f)%f_z
!      !
!    enddo
!    !
!    call MPI_ALLREDUCE(MPI_IN_PLACE,f_x,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,err)
!    call MPI_ALLREDUCE(MPI_IN_PLACE,f_y,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,err)
!    call MPI_ALLREDUCE(MPI_IN_PLACE,f_z,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,err)
!    !
!    return
!  end subroutine CmptWMLESForces
!  !
!  !===================================================================================================================
!  !
!  subroutine ForcesSetCyl(dl)
!    implicit none
!    real(rp),dimension(1:3),intent(in) :: dl
!    integer :: i,f
!    real(rp) :: cr
!    real(rp) cx,cy,cz
!    real(rp) :: t,dt
!    real(rp) :: x,y,z
!    real(rp) :: norm
!    !
!    cr=0.25_rp
!    !
!    cx=0.5_rp
!    cy=0.5_rp
!    cz=0.5_rp
!    !
!    N_FCS=100
!    !
!    allocate(geo(1:N_FCS))
!    !
!    geo(1:N_FCS)%c0_x=0.0_rp
!    geo(1:N_FCS)%c0_y=0.0_rp
!    geo(1:N_FCS)%c0_z=0.0_rp
!    geo(1:N_FCS)%c1_x=0.0_rp
!    geo(1:N_FCS)%c1_y=0.0_rp
!    geo(1:N_FCS)%c1_z=0.0_rp
!    geo(1:N_FCS)%c2_x=0.0_rp
!    geo(1:N_FCS)%c2_y=0.0_rp
!    geo(1:N_FCS)%c2_z=0.0_rp
!    geo(1:N_FCS)%c3_x=0.0_rp
!    geo(1:N_FCS)%c3_y=0.0_rp
!    geo(1:N_FCS)%c3_z=0.0_rp
!    geo(1:N_FCS)%A=0.0_rp
!    geo(1:N_FCS)%f_x=0.0_rp
!    geo(1:N_FCS)%f_y=0.0_rp
!    geo(1:N_FCS)%f_z=0.0_rp
!    geo(1:N_FCS)%n_x=0.0_rp
!    geo(1:N_FCS)%n_y=0.0_rp
!    geo(1:N_FCS)%n_z=0.0_rp
!    !
!    dt=2.0_rp*acos(-1.0_rp)/real(N_FCS,rp)
!    !
!    f=0
!    !
!    do i=1,N_FCS
!      !
!      t=dt*(real(i,rp)-0.5_rp)
!      !
!      x=cx+cr*cos(t)
!      y=cy
!      z=cz+cr*sin(t)
!      !
!      if((x<x_lmin).or.(x>=x_lmax)) cycle
!      if((y<y_lmin).or.(y>=y_lmax)) cycle
!      if((z<z_lmin).or.(z>=z_lmax)) cycle
!      !
!      f=f+1
!      !
!      geo(f)%c0_x=x
!      geo(f)%c0_y=y
!      geo(f)%c0_z=z
!      !
!      norm=sqrt((x-cx)**2+(z-cz)**2)
!      !
!      geo(f)%n_x=(x-cx)/norm
!      geo(f)%n_y=0.0_rp
!      geo(f)%n_z=(z-cz)/norm
!      !
!      geo(f)%A=dt*cr*dl(2)
!      !
!    enddo
!    !
!    N_FCS=f
!    !
!    return
!  end subroutine ForcesSetCyl
!  !
!  !===================================================================================================================v
!  !
!  subroutine ForcesSetChn(n,l,dl)
!    implicit none
!    integer,dimension(1:3),intent(in) :: n
!    real(rp),dimension(1:3),intent(in) :: l,dl
!    integer :: i,j,f
!    real(rp) :: H_w
!    real(rp) :: x,y,z
!    !
!    H_w=1.0_rp/64.0_rp*5.75_rp
!    !
!    N_FCS=n(1)*n(2)*2
!    !
!    allocate(geo(1:N_FCS))
!    !
!    geo(1:N_FCS)%c0_x=0.0_rp
!    geo(1:N_FCS)%c0_y=0.0_rp
!    geo(1:N_FCS)%c0_z=0.0_rp
!    geo(1:N_FCS)%c1_x=0.0_rp
!    geo(1:N_FCS)%c1_y=0.0_rp
!    geo(1:N_FCS)%c1_z=0.0_rp
!    geo(1:N_FCS)%c2_x=0.0_rp
!    geo(1:N_FCS)%c2_y=0.0_rp
!    geo(1:N_FCS)%c2_z=0.0_rp
!    geo(1:N_FCS)%c3_x=0.0_rp
!    geo(1:N_FCS)%c3_y=0.0_rp
!    geo(1:N_FCS)%c3_z=0.0_rp
!    geo(1:N_FCS)%A=0.0_rp
!    geo(1:N_FCS)%f_x=0.0_rp
!    geo(1:N_FCS)%f_y=0.0_rp
!    geo(1:N_FCS)%f_z=0.0_rp
!    geo(1:N_FCS)%n_x=0.0_rp
!    geo(1:N_FCS)%n_y=0.0_rp
!    geo(1:N_FCS)%n_z=0.0_rp
!    !
!    f=0
!    !
!    do j=1,n(2)
!      do i=1,n(1)
!        !
!        x=x_lmin+(real(i,rp)-0.5_rp)*dl(1)
!        y=y_lmin+(real(j,rp)-0.5_rp)*dl(2)
!        z=z_lmin+H_w
!        !
!        if((x<x_lmin).or.(x>=x_lmax)) cycle
!        if((y<y_lmin).or.(y>=y_lmax)) cycle
!        if((z<z_lmin).or.(z>=z_lmax)) cycle
!        !
!        f=f+1
!        !
!        geo(f)%c0_x=x
!        geo(f)%c0_y=y
!        geo(f)%c0_z=z
!        !
!        geo(f)%n_x=0.0_rp
!        geo(f)%n_y=0.0_rp
!        geo(f)%n_z=+1.0_rp
!        !
!        geo(f)%A=dl(1)*dl(2)
!        !
!      enddo
!    enddo
!    !
!    do j=1,n(2)
!      do i=1,n(1)
!        !
!        x=x_lmin+(real(i,rp)-0.5_rp)*dl(1)
!        y=y_lmin+(real(j,rp)-0.5_rp)*dl(2)
!        z=z_lmin+l(3)-H_w
!        !
!        if((x<x_lmin).or.(x>=x_lmax)) cycle
!        if((y<y_lmin).or.(y>=y_lmax)) cycle
!        if((z<z_lmin).or.(z>=z_lmax)) cycle
!        !
!        f=f+1
!        !
!        geo(f)%c0_x=x
!        geo(f)%c0_y=y
!        geo(f)%c0_z=z
!        !
!        geo(f)%n_x=0.0_rp
!        geo(f)%n_y=0.0_rp
!        geo(f)%n_z=-1.0_rp
!        !
!        geo(f)%A=dl(1)*dl(2)
!        !
!      enddo
!    enddo
!    !
!    N_FCS=f
!    !
!    return
!  end subroutine ForcesSetChn
!  !
!  !===================================================================================================================
!  !
!  subroutine ForcesSetPipe(n,l,dl)
!    implicit none
!    integer,dimension(1:3),intent(in) :: n
!    real(rp),dimension(1:3),intent(in) :: l,dl
!    integer :: i,j,f
!    integer :: N_t,N_l
!    real(rp) :: R_c,dt
!    real(rp) :: y_c,z_c
!    real(rp) :: x,y,z,a,t
!    !
!    R_c=0.6_rp
!    y_c=0.75_rp
!    z_c=0.75_rp
!    !
!    N_t=2*nint(pi*R_c/min(dl(1),dl(2),dl(3)))
!    N_l=n(1)
!    !
!    dt=2.0_rp*pi/N_t!real(N_t,rp)
!    !
!    N_FCS=N_t*N_l*2
!    !
!    allocate(geo(1:N_FCS))
!    !
!    geo(1:N_FCS)%c0_x=0.0_rp
!    geo(1:N_FCS)%c0_y=0.0_rp
!    geo(1:N_FCS)%c0_z=0.0_rp
!    geo(1:N_FCS)%c1_x=0.0_rp
!    geo(1:N_FCS)%c1_y=0.0_rp
!    geo(1:N_FCS)%c1_z=0.0_rp
!    geo(1:N_FCS)%c2_x=0.0_rp
!    geo(1:N_FCS)%c2_y=0.0_rp
!    geo(1:N_FCS)%c2_z=0.0_rp
!    geo(1:N_FCS)%c3_x=0.0_rp
!    geo(1:N_FCS)%c3_y=0.0_rp
!    geo(1:N_FCS)%c3_z=0.0_rp
!    geo(1:N_FCS)%A=0.0_rp
!    geo(1:N_FCS)%f_x=0.0_rp
!    geo(1:N_FCS)%f_y=0.0_rp
!    geo(1:N_FCS)%f_z=0.0_rp
!    geo(1:N_FCS)%n_x=0.0_rp
!    geo(1:N_FCS)%n_y=0.0_rp
!    geo(1:N_FCS)%n_z=0.0_rp
!    !
!    f=0
!    !
!    do i=1,N_l
!      do j=1,N_t
!        !
!        a=(real(i,rp)-0.5_rp)*dl(1)
!        t=(real(j,rp)-0.5_rp)*dt
!        !
!        x=a
!        y=y_c+R_c*cos(t)
!        z=z_c+R_c*sin(t)
!        !
!        if((x<x_lmin).or.(x>=x_lmax)) cycle
!        if((y<y_lmin).or.(y>=y_lmax)) cycle
!        if((z<z_lmin).or.(z>=z_lmax)) cycle
!        !
!        f=f+1
!        !
!        geo(f)%c0_x=x
!        geo(f)%c0_y=y
!        geo(f)%c0_z=z
!        !
!        geo(f)%n_x=0.0_rp
!        geo(f)%n_y=-cos(t)
!        geo(f)%n_z=-sin(t)
!        !
!        geo(f)%A=R_c*dt*dl(1)
!        !
!      enddo
!    enddo
!    !
!    N_FCS=f
!    !
!    return
!  end subroutine ForcesSetPipe
  !
  !===================================================================================================================
  !
  function cross(vec_1,vec_2) result(res)
    implicit none
    real(rp),dimension(1:3) :: res
    real(rp),dimension(1:3),intent(in) :: vec_1,vec_2
    !
    res(1)=vec_1(2)*vec_2(3)-vec_1(3)*vec_2(2)
    res(2)=vec_1(3)*vec_2(1)-vec_1(1)*vec_2(3)
    res(3)=vec_1(1)*vec_2(2)-vec_1(2)*vec_2(1)
    !
    return
  end function cross
  !
  !===================================================================================================================
  !
  function norm_2(vec) result(res)
    implicit none
    real(rp) :: res
    real(rp),dimension(1:3),intent(in) :: vec
    !
    res=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
    !
    return
  end function norm_2
  !
  !===================================================================================================================
  !
end module mod_forces
