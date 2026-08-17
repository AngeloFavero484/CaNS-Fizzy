module mod_common
  !
  implicit none
  !
  ! ===================================================================================================================
  !
  ! kinds
  !
  integer,parameter :: rp=selected_real_kind(15,307)
  !
  ! others
  !
  integer :: Np
  real(rp) :: radius
  integer :: it,it_min,it_max,it_out
  real(rp) :: dt,t0,time
  !
  ! particle data
  !
  real(rp), allocatable :: prt(:),rad(:)
  real(rp), allocatable :: prt_x(:),prt_y(:),prt_z(:),prt_p(:)
  !
  !
  ! ===================================================================================================================
  !
end module mod_common
