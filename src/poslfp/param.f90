module mod_param
  !
  implicit none
  !
  integer,parameter :: rp = selected_real_kind(15,307)
  !
!  integer,parameter :: ndims = 2
!  integer,dimension(ndims),parameter :: dims = (/2,2/)
  !
  integer,parameter :: itot=128, jtot=128, kmax=128
!  integer,parameter :: imax=itot/dims(1), jmax=jtot/dims(2)
  !
  real(rp),parameter :: lx=2.0_rp, ly=2.0_rp, lz=2.0_rp
  real(rp),parameter :: dxi=itot/lx, dyi=jtot/ly, dzi=kmax/lz
  !
  real(rp),parameter :: radius=0.1_rp
  real(rp),parameter :: pi=acos(-1.0_rp)
  !
  real(rp),parameter :: retraction=0.3_rp/dxi !lfps retracted from interface into interior of obstacle
  !
  integer,parameter :: NL=nint((pi/3.0_rp)*(12.0_rp*(((radius-retraction)*dxi)**2)+1.0_rp))
  !
  character(len=5),parameter :: datadir = 'data/'
end module mod_param
