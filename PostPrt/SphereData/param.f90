module mod_param
  !
  implicit none
  !
  integer,parameter :: rp = selected_real_kind(15,307)
  !
  integer,parameter :: ng(3)=(/128,128,128/) !itot=128, jtot=128, kmax=128
  integer,parameter :: iter=16
  !
  real(rp),parameter :: visci=100
  real(rp),parameter :: l(3)=(/2.0,2.0,2.0/) !lx=2.0_rp, ly=2.0_rp, lz=2.0_rp
  real(rp),parameter :: dl(3)=l/ng !dxi=itot/lx, dyi=jtot/ly, dzi=kmax/lz
  real(rp),parameter :: dli(3)=1.0_rp/dl
  !
  real(rp),parameter :: radius=0.1_rp
  real(rp),parameter :: cx = 1.0_rp
  real(rp),parameter :: cy = 1.0_rp
  real(rp),parameter :: cz = 1.0_rp
  !
  real(rp),parameter :: pi=acos(-1.0_rp)
  !
  real(rp),parameter :: retraction=0.3_rp*dl(1) !lfps retracted from interface into interior of obstacle
  real(rp),parameter :: radfp=radius-retraction
  !
  integer,parameter :: NL=nint((pi/3.0_rp)*(12.0_rp*(((radius-retraction)*dli(1))**2)+1.0_rp))
  !
  character(len=5),parameter :: datadir = 'data/'
  !
end module mod_param
