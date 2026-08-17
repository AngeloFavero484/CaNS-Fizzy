module mod_interp
  !
!  use mpi
!  use mod_mpi, only: x_lmin,y_lmin,z_lmin
!  use mod_types, only: rp
!  use mod_common_mpi, only: myid,ierr
!  use mod_misc, only: Cmpt3DNorm
  use mod_param     , only: rp
  !
  !===================================================================================================================
  !
  implicit none
  !
!  integer,dimension(:),allocatable :: map_c_z
!  integer :: map_c_n
!  real(rp) :: map_c_dz
!  real(rp) :: map_c_z_min
!  !
!  integer,dimension(:),allocatable :: map_f_z
!  integer :: map_f_n
!  real(rp) :: map_f_dz
!  real(rp) :: map_f_z_min
  !
  !===================================================================================================================
  !
  ! public/private access specifications
  !
  private
  !
!  public :: x_lmin,y_lmin,z_lmin
!  !
!  public :: InterpInit
!  public :: InterpIBMVelX,InterpIBMVelY,InterpIBMVelZ
!  public :: InterpScaDir
!  public :: InterpGradCCC,InterpGradCFF
!  public :: InterpGradFCF,InterpGradFFC
  public :: InterpFastVel
  public :: InterpFastSca
  public :: InterpFastGRP
  !
  !===================================================================================================================
  !
  contains
  !
!  subroutine InterpInit(n,lo,hi,l,dl,zc,zf,dzc,dzf)
!    implicit none
!    integer,dimension(1:3),intent(in) :: n
!    integer,dimension(1:3),intent(in) :: lo,hi
!    real(rp),dimension(1:3),intent(in) :: l,dl
!    real(rp),dimension(-3:n(3)+4),intent(in) :: zc,zf,dzc,dzf
!    !
!    integer :: k,k_map
!    real(rp) :: dum
!    !
!    map_c_dz=minval(dzc)
!    !
!    map_c_z_min=minval(zc)
!    !
!    map_c_n=int((maxval(zc)-minval(zc))/map_c_dz)
!    ! 
!    map_c_dz=(maxval(zc)-minval(zc))/real(map_c_n,rp)
!    !
!    allocate(map_c_z(1:map_c_n))
!    !
!    map_c_z=0
!    !
!    do k_map=1,map_c_n
!      !
!      dum=minval(zc)+(real(k_map,rp)-0.5_rp)*map_c_dz
!      !
!      do k=-3,n(3)+4
!        !
!        if(k>-3) then
!          !
!          if((dum.ge.zc(k)-0.5_rp*dzc(k-1)).and.(dum.lt.zc(k)+0.5_rp*dzc(k))) map_c_z(k_map)=k
!          !
!        elseif(k==-3) then
!          !
!          if((dum.ge.zc(-3)-0.5_rp*2.0_rp*dzc(-3)).and.(dum.lt.zc(-3)+0.5_rp*dzc(-3))) map_c_z(k_map)=k
!          !
!        endif
!        !
!      enddo
!      !
!    enddo
!    !
!    map_f_dz=minval(dzf)
!    !
!    map_f_z_min=minval(zf)
!    !
!    map_f_n=int((maxval(zf)-minval(zf))/map_f_dz)
!    ! 
!    map_f_dz=(maxval(zf)-minval(zf))/real(map_f_n,rp)
!    !
!    allocate(map_f_z(1:map_f_n))
!    !
!    map_f_z=0
!    !
!    do k_map=1,map_f_n
!      !
!      dum=minval(zf)+(real(k_map,rp)-0.5_rp)*map_f_dz
!      !
!      do k=-3,n(3)+4
!        !
!        if(k>-3) then
!          !
!          if((dum.ge.zf(k)-0.5_rp*dzf(k-1)).and.(dum.lt.zf(k)+0.5_rp*dzf(k))) map_f_z(k_map)=k
!          !
!        elseif(k==-3) then
!          !
!          if((dum.ge.zf(-3)-0.5_rp*2.0_rp*dzf(-3)).and.(dum.lt.zf(-3)+0.5_rp*dzf(-3))) map_f_z(k_map)=k
!          !
!        endif
!        !
!      enddo
!      !
!    enddo
!    !
!    return
!  end subroutine InterpInit
  !
  !===================================================================================================================
  !
  subroutine InterpFastVel(n,dl,point,u,v,w,v_x,v_y,v_z)
    implicit none
    integer,dimension(1:3),intent(in) :: n
    real(rp),dimension(1:3),intent(in) :: dl
    real(rp),dimension(1:3),intent(in) :: point
    real(rp),dimension(1:n(1),1:n(2),1:n(3)),intent(in) :: u,v,w
    real(rp),intent(out) :: v_x,v_y,v_z
    integer :: i,j,k
    integer :: ii,jj,kk
    !
    real(rp),dimension(-1:1,-1:1) :: C1_x,C1_y,C1_z
    real(rp),dimension(-1:0) :: C2_x,C2_y,C2_z
    !
    real(rp) :: xc_0,xc_p1,xc_m1
    real(rp) :: yc_0,yc_p1,yc_m1
    real(rp) :: zc_0,zc_p1,zc_m1
    !
    real(rp) :: x,y,z
    real(rp) :: x_a,y_a,z_a
    !
    x=point(1)
    y=point(2)
    z=point(3)
    !
    ! setting indexes
    !
    i=floor(x/dl(1))+1
    j=floor(y/dl(2))+1
    k=floor(z/dl(3))+1
    !
    i=max(i,1)
    j=max(j,1)
    k=max(k,1)
    !
    i=min(i,n(1))
    j=min(j,n(2))
    k=min(k,n(3))
    !
    xc_0=(real(i,rp)-0.5_rp)*dl(1)
    yc_0=(real(j,rp)-0.5_rp)*dl(2)
    zc_0=(real(k,rp)-0.5_rp)*dl(3)
    !
    xc_p1=xc_0+dl(1)
    yc_p1=yc_0+dl(2)
    zc_p1=zc_0+dl(3)
    !
    xc_m1=xc_0-dl(1)
    yc_m1=yc_0-dl(2)
    zc_m1=zc_0-dl(3)
    !
    x_a=2.0_rp*(x-xc_0)/(xc_p1-xc_m1)
    y_a=2.0_rp*(y-yc_0)/(yc_p1-yc_m1)
    z_a=2.0_rp*(z-zc_0)/(zc_p1-zc_m1)
    !
    C1_x(-1,-1)=0.25_rp*y_a*(y_a-1.0_rp)*z_a*(z_a-1.0_rp)
    C1_x(-1, 0)=0.50_rp*y_a*(y_a-1.0_rp)*(1.0_rp-z_a**2)
    C1_x(-1, 1)=0.25_rp*y_a*(y_a-1.0_rp)*z_a*(z_a+1.0_rp)
    C1_x( 0,-1)=0.50_rp*(1.0_rp-y_a**2)*z_a*(z_a-1.0_rp)
    C1_x( 0, 0)=1.00_rp*(1.0_rp-y_a**2)*(1.0_rp-z_a**2)
    C1_x( 0, 1)=0.50_rp*(1.0_rp-y_a**2)*z_a*(z_a+1.0_rp)
    C1_x( 1,-1)=0.25_rp*y_a*(y_a+1.0_rp)*z_a*(z_a-1.0_rp)
    C1_x( 1, 0)=0.50_rp*y_a*(y_a+1.0_rp)*(1.0_rp-z_a**2)
    C1_x( 1, 1)=0.25_rp*y_a*(y_a+1.0_rp)*z_a*(z_a+1.0_rp)
    !
    C2_x( 0)=+2.0_rp*(x-0.5_rp*(xc_0+xc_m1))/(xc_p1-xc_m1)
    C2_x(-1)=-2.0_rp*(x-0.5_rp*(xc_0+xc_p1))/(xc_p1-xc_m1)
    !
    C1_y(-1,-1)=0.25_rp*x_a*(x_a-1.0_rp)*z_a*(z_a-1.0_rp)
    C1_y(-1, 0)=0.50_rp*x_a*(x_a-1.0_rp)*(1.0_rp-z_a**2)
    C1_y(-1, 1)=0.25_rp*x_a*(x_a-1.0_rp)*z_a*(z_a+1.0_rp)
    C1_y( 0,-1)=0.50_rp*(1.0_rp-x_a**2)*z_a*(z_a-1.0_rp)
    C1_y( 0, 0)=1.00_rp*(1.0_rp-x_a**2)*(1.0_rp-z_a**2)
    C1_y( 0, 1)=0.50_rp*(1.0_rp-x_a**2)*z_a*(z_a+1.0_rp)
    C1_y( 1,-1)=0.25_rp*x_a*(x_a+1.0_rp)*z_a*(z_a-1.0_rp)
    C1_y( 1, 0)=0.50_rp*x_a*(x_a+1.0_rp)*(1.0_rp-z_a**2)
    C1_y( 1, 1)=0.25_rp*x_a*(x_a+1.0_rp)*z_a*(z_a+1.0_rp)
    !
    C2_y( 0)=+2.0_rp*(y-0.5_rp*(yc_0+yc_m1))/(yc_p1-yc_m1)
    C2_y(-1)=-2.0_rp*(y-0.5_rp*(yc_0+yc_p1))/(yc_p1-yc_m1)
    !
    C1_z(-1,-1)=0.25_rp*x_a*(x_a-1.0_rp)*y_a*(y_a-1.0_rp)
    C1_z(-1, 0)=0.50_rp*x_a*(x_a-1.0_rp)*(1.0_rp-y_a**2)
    C1_z(-1, 1)=0.25_rp*x_a*(x_a-1.0_rp)*y_a*(y_a+1.0_rp)
    C1_z( 0,-1)=0.50_rp*(1.0_rp-x_a**2)*y_a*(y_a-1.0_rp)
    C1_z( 0, 0)=1.00_rp*(1.0_rp-x_a**2)*(1.0_rp-y_a**2)
    C1_z( 0, 1)=0.50_rp*(1.0_rp-x_a**2)*y_a*(y_a+1.0_rp)
    C1_z( 1,-1)=0.25_rp*x_a*(x_a+1.0_rp)*y_a*(y_a-1.0_rp)
    C1_z( 1, 0)=0.50_rp*x_a*(x_a+1.0_rp)*(1.0_rp-y_a**2)
    C1_z( 1, 1)=0.25_rp*x_a*(x_a+1.0_rp)*y_a*(y_a+1.0_rp)
    !
    C2_z( 0)=+2.0_rp*(z-0.5_rp*(zc_0+zc_m1))/(zc_p1-zc_m1)
    C2_z(-1)=-2.0_rp*(z-0.5_rp*(zc_0+zc_p1))/(zc_p1-zc_m1)
    !
    v_x=0.0_rp
    v_y=0.0_rp
    v_z=0.0_rp
    !
    do kk=-1,1
      do jj=-1,1
        do ii=-1,0
          !
          v_x=v_x+C1_x(jj,kk)*C2_x(ii)*u(i+ii,j+jj,k+kk)
          !
        enddo
      enddo
    enddo
    !
    do kk=-1,1
      do jj=-1,0
        do ii=-1,1
          !
          v_y=v_y+C1_y(ii,kk)*C2_y(jj)*v(i+ii,j+jj,k+kk)
          !
        enddo
      enddo
    enddo
    !
    do kk=-1,0
      do jj=-1,1
        do ii=-1,1
          !
          v_z=v_z+C1_z(ii,jj)*C2_z(kk)*w(i+ii,j+jj,k+kk)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine InterpFastVel
  !
  !===================================================================================================================
  !
  subroutine InterpFastSca(n,dl,point,s,q)
    implicit none
    integer,dimension(1:3),intent(in) :: n
    real(rp),dimension(1:3),intent(in) :: dl
    real(rp),dimension(1:3),intent(in) :: point
    real(rp),dimension(1:n(1),1:n(2),1:n(3)),intent(in) :: s
    real(rp),intent(out) :: q
    !
    integer :: i,j,k
    integer :: ii,jj,kk
    !
    real(rp),dimension(-1:1,-1:1,-1:1) :: C
    !
    real(rp) :: xc_0,xc_p1,xc_m1
    real(rp) :: yc_0,yc_p1,yc_m1
    real(rp) :: zc_0,zc_p1,zc_m1
    !
    real(rp) :: x,y,z
    real(rp) :: x_a,y_a,z_a
    !
    x=point(1)
    y=point(2)
    z=point(3)
    !
    ! setting indexes
    !
    i=floor(x/dl(1))+1 
    j=floor(y/dl(2))+1
    k=floor(z/dl(3))+1
    !
    i=max(i,1)
    j=max(j,1)
    k=max(k,1)
    !
    i=min(i,n(1))
    j=min(j,n(2))
    k=min(k,n(3))
    !
    xc_0=(real(i,rp)-0.5_rp)*dl(1)
    yc_0=(real(j,rp)-0.5_rp)*dl(2)
    zc_0=(real(k,rp)-0.5_rp)*dl(3)
    !
    xc_p1=xc_0+dl(1)
    yc_p1=yc_0+dl(2)
    zc_p1=zc_0+dl(3)
    !
    xc_m1=xc_0-dl(1)
    yc_m1=yc_0-dl(2)
    zc_m1=zc_0-dl(3)
    !
    x_a=2.0_rp*(x-xc_0)/(xc_p1-xc_m1)
    y_a=2.0_rp*(y-yc_0)/(yc_p1-yc_m1)
    z_a=2.0_rp*(z-zc_0)/(zc_p1-zc_m1)
    !
    C(-1,-1,-1)=0.125_rp*x_a*(x_a-1.0_rp)*y_a*(y_a-1.0_rp   )*z_a*(z_a-1.0_rp)
    C(-1,-1, 0)=0.250_rp*x_a*(x_a-1.0_rp)*y_a*(y_a-1.0_rp   )*(1.0_rp-z_a**2 )
    C(-1,-1, 1)=0.125_rp*x_a*(x_a-1.0_rp)*y_a*(y_a-1.0_rp   )*z_a*(z_a+1.0_rp)
    C(-1, 0,-1)=0.250_rp*x_a*(x_a-1.0_rp)*    (1.0_rp-y_a**2)*z_a*(z_a-1.0_rp)
    C(-1, 0, 0)=0.500_rp*x_a*(x_a-1.0_rp)*    (1.0_rp-y_a**2)*(1.0_rp-z_a**2 )
    C(-1, 0, 1)=0.250_rp*x_a*(x_a-1.0_rp)*    (1.0_rp-y_a**2)*z_a*(z_a+1.0_rp)
    C(-1, 1,-1)=0.125_rp*x_a*(x_a-1.0_rp)*y_a*(y_a+1.0_rp   )*z_a*(z_a-1.0_rp)
    C(-1, 1, 0)=0.250_rp*x_a*(x_a-1.0_rp)*y_a*(y_a+1.0_rp   )*(1.0_rp-z_a**2 )
    C(-1, 1, 1)=0.125_rp*x_a*(x_a-1.0_rp)*y_a*(y_a+1.0_rp   )*z_a*(z_a+1.0_rp)
    C( 0,-1,-1)=0.250_rp*(1.0_rp-x_a**2 )*y_a*(y_a-1.0_rp   )*z_a*(z_a-1.0_rp)
    C( 0,-1, 0)=0.500_rp*(1.0_rp-x_a**2 )*y_a*(y_a-1.0_rp   )*(1.0_rp-z_a**2 )
    C( 0,-1, 1)=0.250_rp*(1.0_rp-x_a**2 )*y_a*(y_a-1.0_rp   )*z_a*(z_a+1.0_rp)
    C( 0, 0,-1)=0.500_rp*(1.0_rp-x_a**2 )*    (1.0_rp-y_a**2)*z_a*(z_a-1.0_rp)
    C( 0, 0, 0)=         (1.0_rp-x_a**2 )*    (1.0_rp-y_a**2)*(1.0_rp-z_a**2 )
    C( 0, 0, 1)=0.500_rp*(1.0_rp-x_a**2 )*    (1.0_rp-y_a**2)*z_a*(z_a+1.0_rp)
    C( 0, 1,-1)=0.250_rp*(1.0_rp-x_a**2 )*y_a*(y_a+1.0_rp   )*z_a*(z_a-1.0_rp)
    C( 0, 1, 0)=0.500_rp*(1.0_rp-x_a**2 )*y_a*(y_a+1.0_rp   )*(1.0_rp-z_a**2 )
    C( 0, 1, 1)=0.250_rp*(1.0_rp-x_a**2 )*y_a*(y_a+1.0_rp   )*z_a*(z_a+1.0_rp)
    C( 1,-1,-1)=0.125_rp*x_a*(x_a+1.0_rp)*y_a*(y_a-1.0_rp   )*z_a*(z_a-1.0_rp)
    C( 1,-1, 0)=0.250_rp*x_a*(x_a+1.0_rp)*y_a*(y_a-1.0_rp   )*(1.0_rp-z_a**2 )
    C( 1,-1, 1)=0.125_rp*x_a*(x_a+1.0_rp)*y_a*(y_a-1.0_rp   )*z_a*(z_a+1.0_rp)
    C( 1, 0,-1)=0.250_rp*x_a*(x_a+1.0_rp)*    (1.0_rp-y_a**2)*z_a*(z_a-1.0_rp)
    C( 1, 0, 0)=0.500_rp*x_a*(x_a+1.0_rp)*    (1.0_rp-y_a**2)*(1.0_rp-z_a**2 )
    C( 1, 0, 1)=0.250_rp*x_a*(x_a+1.0_rp)*    (1.0_rp-y_a**2)*z_a*(z_a+1.0_rp)
    C( 1, 1,-1)=0.125_rp*x_a*(x_a+1.0_rp)*y_a*(y_a+1.0_rp   )*z_a*(z_a-1.0_rp)
    C( 1, 1, 0)=0.250_rp*x_a*(x_a+1.0_rp)*y_a*(y_a+1.0_rp   )*(1.0_rp-z_a**2 )
    C( 1, 1, 1)=0.125_rp*x_a*(x_a+1.0_rp)*y_a*(y_a+1.0_rp   )*z_a*(z_a+1.0_rp)
    !
    q=0.0_rp
    !
    do kk=-1,1
      do jj=-1,1
        do ii=-1,1
          !
          q=q+C(ii,jj,kk)*s(i+ii,j+jj,k+kk)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine InterpFastSca
  !
  !===================================================================================================================
  !
  subroutine InterpFastGRP(n,dl,point,gp_x,gp_y,gp_z,gq_x,gq_y,gq_z)
    implicit none
    integer,dimension(1:3),intent(in) :: n
    real(rp),dimension(1:3),intent(in) :: dl
    real(rp),dimension(1:3),intent(in) :: point
    real(rp),dimension(1:n(1),1:n(2),1:n(3)),intent(in) :: gp_x,gp_y,gp_z
    real(rp),intent(out) :: gq_x,gq_y,gq_z
    !
    integer :: i,j,k
    integer :: ii,jj,kk
    !
    real(rp),dimension(-1:1,-1:1) :: C1_x,C1_y,C1_z
    real(rp),dimension(-1:0) :: C2_x,C2_y,C2_z
    !
    real(rp) :: xc_0,xc_p1,xc_m1
    real(rp) :: yc_0,yc_p1,yc_m1
    real(rp) :: zc_0,zc_p1,zc_m1
    !
    real(rp) :: x,y,z
    real(rp) :: x_a,y_a,z_a
    !
    x=point(1)
    y=point(2)
    z=point(3)
    !
    i=floor(x/dl(1))+1 
    j=floor(y/dl(2))+1
    k=floor(z/dl(3))+1
    !
    i=max(i,1)
    j=max(j,1)
    k=max(k,1)
    !
    i=min(i,n(1))
    j=min(j,n(2))
    k=min(k,n(3))
    !
    xc_0=(real(i,rp)-0.5_rp)*dl(1)
    yc_0=(real(j,rp)-0.5_rp)*dl(2)
    zc_0=(real(k,rp)-0.5_rp)*dl(3)
    !
    xc_p1=xc_0+dl(1)
    yc_p1=yc_0+dl(2)
    zc_p1=zc_0+dl(3)
    !
    xc_m1=xc_0-dl(1)
    yc_m1=yc_0-dl(2)
    zc_m1=zc_0-dl(3)
    !
    x_a=2.0_rp*(x-xc_0)/(xc_p1-xc_m1)
    y_a=2.0_rp*(y-yc_0)/(yc_p1-yc_m1)
    z_a=2.0_rp*(z-zc_0)/(zc_p1-zc_m1)
    !
    C1_x(-1,-1)=0.25_rp*y_a*(y_a-1.0_rp)*z_a*(z_a-1.0_rp)
    C1_x(-1, 0)=0.50_rp*y_a*(y_a-1.0_rp)*(1.0_rp-z_a**2)
    C1_x(-1, 1)=0.25_rp*y_a*(y_a-1.0_rp)*z_a*(z_a+1.0_rp)
    C1_x( 0,-1)=0.50_rp*(1.0_rp-y_a**2)*z_a*(z_a-1.0_rp)
    C1_x( 0, 0)=1.00_rp*(1.0_rp-y_a**2)*(1.0_rp-z_a**2)
    C1_x( 0, 1)=0.50_rp*(1.0_rp-y_a**2)*z_a*(z_a+1.0_rp)
    C1_x( 1,-1)=0.25_rp*y_a*(y_a+1.0_rp)*z_a*(z_a-1.0_rp)
    C1_x( 1, 0)=0.50_rp*y_a*(y_a+1.0_rp)*(1.0_rp-z_a**2)
    C1_x( 1, 1)=0.25_rp*y_a*(y_a+1.0_rp)*z_a*(z_a+1.0_rp)
    !
    C2_x( 0)=+2.0_rp*(x-0.5_rp*(xc_0+xc_m1))/(xc_p1-xc_m1)
    C2_x(-1)=-2.0_rp*(x-0.5_rp*(xc_0+xc_p1))/(xc_p1-xc_m1)
    !
    C1_y(-1,-1)=0.25_rp*x_a*(x_a-1.0_rp)*z_a*(z_a-1.0_rp)
    C1_y(-1, 0)=0.50_rp*x_a*(x_a-1.0_rp)*(1.0_rp-z_a**2)
    C1_y(-1, 1)=0.25_rp*x_a*(x_a-1.0_rp)*z_a*(z_a+1.0_rp)
    C1_y( 0,-1)=0.50_rp*(1.0_rp-x_a**2)*z_a*(z_a-1.0_rp)
    C1_y( 0, 0)=1.00_rp*(1.0_rp-x_a**2)*(1.0_rp-z_a**2)
    C1_y( 0, 1)=0.50_rp*(1.0_rp-x_a**2)*z_a*(z_a+1.0_rp)
    C1_y( 1,-1)=0.25_rp*x_a*(x_a+1.0_rp)*z_a*(z_a-1.0_rp)
    C1_y( 1, 0)=0.50_rp*x_a*(x_a+1.0_rp)*(1.0_rp-z_a**2)
    C1_y( 1, 1)=0.25_rp*x_a*(x_a+1.0_rp)*z_a*(z_a+1.0_rp)
    !
    C2_y( 0)=+2.0_rp*(y-0.5_rp*(yc_0+yc_m1))/(yc_p1-yc_m1)
    C2_y(-1)=-2.0_rp*(y-0.5_rp*(yc_0+yc_p1))/(yc_p1-yc_m1)
    !
    C1_z(-1,-1)=0.25_rp*x_a*(x_a-1.0_rp)*y_a*(y_a-1.0_rp)
    C1_z(-1, 0)=0.50_rp*x_a*(x_a-1.0_rp)*(1.0_rp-y_a**2)
    C1_z(-1, 1)=0.25_rp*x_a*(x_a-1.0_rp)*y_a*(y_a+1.0_rp)
    C1_z( 0,-1)=0.50_rp*(1.0_rp-x_a**2)*y_a*(y_a-1.0_rp)
    C1_z( 0, 0)=1.00_rp*(1.0_rp-x_a**2)*(1.0_rp-y_a**2)
    C1_z( 0, 1)=0.50_rp*(1.0_rp-x_a**2)*y_a*(y_a+1.0_rp)
    C1_z( 1,-1)=0.25_rp*x_a*(x_a+1.0_rp)*y_a*(y_a-1.0_rp)
    C1_z( 1, 0)=0.50_rp*x_a*(x_a+1.0_rp)*(1.0_rp-y_a**2)
    C1_z( 1, 1)=0.25_rp*x_a*(x_a+1.0_rp)*y_a*(y_a+1.0_rp)
    !
    C2_z( 0)=+2.0_rp*(z-0.5_rp*(zc_0+zc_m1))/(zc_p1-zc_m1)
    C2_z(-1)=-2.0_rp*(z-0.5_rp*(zc_0+zc_p1))/(zc_p1-zc_m1)
    !
    gq_x=0.0_rp
    !
    do kk=-1,1
      do jj=-1,1
        do ii=-1,0
          !
          gq_x=gq_x+C1_x(jj,kk)*C2_x(ii)*gp_x(i+ii,j+jj,k+kk)
          !
        enddo
      enddo
    enddo
    !
    gq_y=0.0_rp
    !
    do kk=-1,1
      do jj=-1,0
        do ii=-1,1
          !
          gq_y=gq_y+C1_y(ii,kk)*C2_y(jj)*gp_y(i+ii,j+jj,k+kk)
          !
        enddo
      enddo
    enddo
    !
    gq_z=0.0_rp
    !
    do kk=-1,0
      do jj=-1,1
        do ii=-1,1
          !
          gq_z=gq_z+C1_z(ii,jj)*C2_z(kk)*gp_z(i+ii,j+jj,k+kk)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine InterpFastGRP
  !
  !===================================================================================================================
  !
end module mod_interp
