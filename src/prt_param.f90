module prt_mod_param
#if defined(_PARTICLE)
!  use decomp_2d
  use mod_types
  use mod_param  , only: pi,l
  !
  implicit none
  !
!  integer, parameter :: ndims = 2
!  integer, dimension(ndims), parameter :: dims = (/ 100, 200/)
!  integer, parameter :: itot = 200*16, jtot = 300*16, ktot = 100*16
!  integer, parameter :: it1 = itot+1, jt1 = jtot+1, kt1 = ktot+1
!  integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
!  integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
!  real, parameter :: lx = 200.,ly = 300.,lz = 100.
!  real, parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
!  real, parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
!  real, parameter :: pi = acos(-1.)
! Width hyperbola
#if defined(_EULER)
    real(rp), protected :: eps_sol
#endif
  !
  integer, protected :: np
  ! fixed compile-time cap on simultaneous collision partners+walls per particle
  ! (sizes fixed-shape array components of type particle/particle_old in prt_common.f90,
  ! so it cannot be derived from the runtime np; large enough for any np of practical interest)
  integer, parameter :: nqmax = 15
!  integer, parameter :: nfriendsmax = 80
  !
  integer, parameter :: send_real = 24+10*nqmax, &
                        send_int = 1+nqmax
!  real, parameter ::Reb=5600.
!  ! desired bulk Reynolds number based on superficial bulk velocity and channel height
!  real, parameter :: visc=lz/Reb
!  real, parameter :: bulk_v_sup = Reb*visc/lz !equal to 1. for chosen value of visc
!  real, parameter :: gacc = -9.81  !*lz/bulk_v_sup**2.*0!diam/vscale**2.
!  real(rp), parameter :: gaccx = 0.0_rp, &  !gacc*cos(pi/2.), &
!                         gaccy = 0.0_rp, &     !gacc*cos(pi/2.), &
!                         gaccz = -9.81_rp     !gacc*sin(pi/2.)
  !
  ! particle physical inputs (read from input.nml, see read_particle_input below)
  !
  real(rp), protected :: ratiorho
  real(rp), protected :: rho_s
  real(rp), protected :: radius
  real(rp), protected :: volp
  real(rp), protected :: mominert
  real(rp), protected :: u_ini,v_ini,w_ini ! initial particle velocity (same for all particles)
  ! initial particle position, only used when np=1 (np>1 is placed pseudo-randomly, see initparticles);
  ! left at pos_ini_unset falls back to the legacy default placement (domain-centered x,y, z=0.755*l(3))
  real(rp), protected :: x_ini,y_ini,z_ini
  real(rp), parameter :: pos_ini_unset = -1.0e30_rp
!  character(len=5), parameter :: datadir = 'data/'
!  !
!  real, parameter, dimension(3,2) :: rkcoeff = reshape((/ 32./60., 25./60., 45./60., 0., -17./60., -25./60. /), shape(rkcoeff))
!  real, parameter, dimension(3) :: rkcoeffab = rkcoeff(:,1)+rkcoeff(:,2)
!  !
!  !Output parameters
!  integer,parameter :: ioutchk = 10/1, iout1d = 50/1,iout2d = 3000 ,ioutfld = iout2d 
  !
  ! set collision parameters (read from input.nml, see read_particle_input below)
  !
  real(rp), protected :: Nstretch,dt_estim,r_dtcol,r_dtcoli
  real(rp), protected :: en,et,muc
  ! sphere/sphere (derived from radius/ratiorho, see read_particle_input below)
  real(rp), protected :: colthr_pp
  real(rp), protected :: meffn_ss,mefft_ss,kn_ss,kt_ss,etan_ss,etat_ss,muc_ss,psi_crit_ss
  ! sphere/wall (derived from radius/ratiorho, see read_particle_input below)
  real(rp), protected :: colthr_pw
  real(rp), protected :: meffn_sw,mefft_sw,kn_sw,kt_sw,etan_sw,etat_sw,muc_sw,psi_crit_sw
  !
  ! set parameters for the lubrication model
  !
  ! sphere/sphere
  real(rp), parameter :: eps_ini_pp = 0.025_rp, &
                         eps_sat_pp = 0.001_rp, &
                         eps_cut_pp = 0.0_rp
  ! Correction with Stokes amplification factor is added for values of epsilon smaller than
  ! eps_ini_pp
  real(rp), parameter :: &
       a11_ini_pp = -1._rp/4._rp*eps_ini_pp**(-1._rp)+9._rp/40._rp*log(eps_ini_pp)+3._rp/112._rp*eps_ini_pp*log(eps_ini_pp), &
       a11a_ini_pp = a11_ini_pp-0.995_rp, &
       a11b_ini_pp = -a11_ini_pp+0.350_rp, &
       a22_ini_pp = 1._rp/6._rp*log(eps_ini_pp), &
       a22a_ini_pp = a22_ini_pp-0.998_rp, &
       a22b_ini_pp = -a22_ini_pp+0.274_rp, &
       a33a_ini_pp = a22a_ini_pp, &
       a33b_ini_pp = a22b_ini_pp, &
       b23_ini_pp = -1._rp/6._rp*log(eps_ini_pp)-1._rp/12._rp*eps_ini_pp*log(eps_ini_pp), &
       b23a_ini_pp = b23_ini_pp-0.159_rp, &
       b23b_ini_pp = -b23_ini_pp+0.001_rp, &
       b32a_ini_pp = -b23a_ini_pp, &
       b32b_ini_pp = -b23b_ini_pp, &
       c23a_ini_pp = b32a_ini_pp, &
       c23b_ini_pp = b32b_ini_pp, &
       c32a_ini_pp = b23a_ini_pp, &
       c32b_ini_pp = b23b_ini_pp, &
       d11_ini_pp = 1._rp/8._rp*eps_ini_pp*log(eps_ini_pp), &
       d11a_ini_pp = 1._rp/8._rp*eps_ini_pp*log(eps_ini_pp), &
       d11b_ini_pp = -1._rp/8._rp*eps_ini_pp*log(eps_ini_pp), &
       d22a_ini_pp = 1._rp/5._rp*log(eps_ini_pp)+47._rp/250._rp*eps_ini_pp*log(eps_ini_pp)-0.703_rp, &
       d22b_ini_pp = -1._rp/20._rp*log(eps_ini_pp)+31._rp/500._rp*eps_ini_pp*log(eps_ini_pp)-0.027_rp, &
       d33a_ini_pp = 1._rp/5._rp*log(eps_ini_pp)+47._rp/250._rp*eps_ini_pp*log(eps_ini_pp)-0.703_rp, &
       d33b_ini_pp = -1._rp/20._rp*log(eps_ini_pp)+31._rp/500._rp*eps_ini_pp*log(eps_ini_pp)-0.027_rp
  ! Stokes amplification factor saturated for values of epsilon smaller than eps_sat_pp
  real(rp), parameter :: &
       a11_sat_pp = -1._rp/4._rp*eps_sat_pp**(-1._rp)+9._rp/40._rp*log(eps_sat_pp)+3._rp/112._rp*eps_sat_pp*log(eps_sat_pp), &
       a11a_sat_pp = a11_sat_pp-0.995_rp, &
       a11b_sat_pp = -a11_sat_pp+0.350_rp, &
       a22_sat_pp = 1._rp/6._rp*log(eps_sat_pp), &
       a22a_sat_pp = a22_sat_pp-0.998_rp, &  
       a22b_sat_pp = -a22_sat_pp+0.274_rp, & 
       a33a_sat_pp = a22a_sat_pp, &
       a33b_sat_pp = a22b_sat_pp, &
       b23_sat_pp = -1._rp/6._rp*log(eps_sat_pp)-1._rp/12._rp*eps_sat_pp*log(eps_sat_pp), &
       b23a_sat_pp = b23_sat_pp-0.159_rp, &
       b23b_sat_pp = -b23_sat_pp+0.001_rp, &
       b32a_sat_pp = -b23a_sat_pp, &
       b32b_sat_pp = -b23b_sat_pp, &
       c23a_sat_pp = b32a_sat_pp, &
       c23b_sat_pp = b32b_sat_pp, &
       c32a_sat_pp = b23a_sat_pp, &
       c32b_sat_pp = b23b_sat_pp, &
       d11_sat_pp = 1._rp/8._rp*eps_sat_pp*log(eps_sat_pp), &
       d11a_sat_pp = 1._rp/8._rp*eps_sat_pp*log(eps_sat_pp), &
       d11b_sat_pp = -1._rp/8._rp*eps_sat_pp*log(eps_sat_pp), &
       d22a_sat_pp = 1._rp/5._rp*log(eps_sat_pp)+47._rp/250._rp*eps_sat_pp*log(eps_sat_pp)-0.703_rp, &
       d22b_sat_pp = -1._rp/20._rp*log(eps_sat_pp)+31._rp/500._rp*eps_sat_pp*log(eps_sat_pp)-0.027_rp, &
       d33a_sat_pp = 1._rp/5._rp*log(eps_sat_pp)+47._rp/250._rp*eps_sat_pp*log(eps_sat_pp)-0.703_rp, &
       d33b_sat_pp = -1._rp/20._rp*log(eps_sat_pp)+31._rp/500._rp*eps_sat_pp*log(eps_sat_pp)-0.027_rp
  ! sphere/wall
  real(rp), parameter :: eps_ini_pw = 0.075_rp, &
                         eps_sat_pw = 0.001_rp, &
                         eps_cut_pw = 0._rp
  ! Correction with Stokes amplification factor is added for values of epsilon smaller than
  ! eps_ini_pw
  real(rp), parameter :: &
       a11_ini_pw = -1._rp/eps_ini_pw+1._rp/5._rp*log(eps_ini_pw)+1._rp/21._rp*eps_ini_pw*log(eps_ini_pw)-0.9713_rp, &
       a22_ini_pw = 8._rp/15._rp*log(eps_ini_pw)+64._rp/375._rp*eps_ini_pw*log(eps_ini_pw)-0.952_Rp, &
       a33_ini_pw = a22_ini_pw, &
       b23_ini_pw = -2._rp/15._rp*log(eps_ini_pw)-86._rp/375._rp*eps_ini_pw*log(eps_ini_pw)-0.257_rp, &
       b32_ini_pw = -b23_ini_pw, &
       c23_ini_pw = b32_ini_pw, &
       c32_ini_pw = b23_ini_pw, &
       d11_ini_pw = 1._rp/2._rp*eps_ini_pw*log(eps_ini_pw)-1.202_rp, &
       d22_ini_pw = 2._rp/5._rp*log(eps_ini_pw)+66._rp/125._rp*eps_ini_pw*log(eps_ini_pw)-0.371_rp, &
       d33_ini_pw = 2._rp/5._rp*log(eps_ini_pw)+66._rp/125._rp*eps_ini_pw*log(eps_ini_pw)-0.371_rp
  ! Stokes amplification factor saturated for values of epsilon smaller than eps_sat_pw
  real(rp), parameter :: &
       a11_sat_pw = -1._rp/eps_sat_pw+1._rp/5._rp*log(eps_sat_pw)+1._rp/21._rp*eps_sat_pw*log(eps_sat_pw)-0.9713_rp, &
       a22_sat_pw = 8._rp/15._rp*log(eps_sat_pw)+64._rp/375._rp*eps_sat_pw*log(eps_sat_pw)-0.952_rp, &
       a33_sat_pw = a22_sat_pw, &
       b23_sat_pw = -2._rp/15._rp*log(eps_sat_pw)-86._rp/375._rp*eps_sat_pw*log(eps_sat_pw)-0.257_rp, &
       b32_sat_pw = -b23_sat_pw, &
       c23_sat_pw = b32_sat_pw, &
       c32_sat_pw = b23_sat_pw, &
       d11_sat_pw = 1._rp/2._rp*eps_sat_pw*log(eps_sat_pw)-1.202_rp, &
       d22_sat_pw = 2._rp/5._rp*log(eps_sat_pw)+66._rp/125._rp*eps_sat_pw*log(eps_sat_pw)-0.371_rp, &
       d33_sat_pw = 2._rp/5._rp*log(eps_sat_pw)+66._rp/125._rp*eps_sat_pw*log(eps_sat_pw)-0.371_rp
  !
  contains
  subroutine read_particle_input(myid)
    use mpi
    implicit none
    integer, intent(in) :: myid
    integer :: iunit,ierr
    namelist /particle/ np,radius,rho_s,ratiorho,u_ini,v_ini,w_ini,x_ini,y_ini,z_ini
    namelist /collision_parameters/ Nstretch,dt_estim,r_dtcol,en,et,muc
#if defined(_EULER)
    namelist /particle_euler/ eps_sol
#endif
    !
    ! set-up default parameters
    !
    np = 1
    radius = 1._rp
    rho_s = 1320._rp
    ratiorho = 5._rp
    u_ini = 0.0_rp
    v_ini = 0.0_rp
    w_ini = -28.78_rp
    x_ini = pos_ini_unset
    y_ini = pos_ini_unset
    z_ini = pos_ini_unset
    Nstretch = 8.0_rp
    dt_estim = 0.003_rp !0.05_rp !0.003
    r_dtcol  = 50.0_rp   !=dt/dtp
    en = 0.97_rp
    et = 0.10_rp
    muc = 0.0_rp
#if defined(_EULER)
    eps_sol = 1.5_rp
#endif
    !
    ! read input file
    !
    open(newunit=iunit,file='input.nml',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,nml=particle,iostat=ierr)
        rewind(iunit)
        read(iunit,nml=collision_parameters,iostat=ierr)
#if defined(_EULER)
        rewind(iunit)
        read(iunit,nml=particle_euler,iostat=ierr)
#endif
      else
        if(myid == 0) print*, 'Error reading the input file'
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        error stop
      end if
    close(iunit)
    !
    ! derived quantities
    !
    if(np == 1) then
      if(x_ini == pos_ini_unset) x_ini = l(1)*0.5_rp
      if(y_ini == pos_ini_unset) y_ini = l(2)*0.5_rp
      if(z_ini == pos_ini_unset) z_ini = l(3)*0.755_rp
    endif
    r_dtcoli = 1.0_rp/r_dtcol
    volp = (4._rp/3._rp)*pi*radius**3._rp
    mominert = (2._rp/5._rp)*volp*radius**2._rp
    !
    ! sphere/sphere
    !
    colthr_pp = 0._rp*radius
    meffn_ss = ratiorho*volp/2._rp
    mefft_ss = 2._rp/7._rp*meffn_ss
    kn_ss = (pi**2._rp + abs(log(en))**2._rp)*meffn_ss/((Nstretch*dt_estim)**2._rp)
    kt_ss = (pi**2._rp + abs(log(et))**2._rp)*mefft_ss/((Nstretch*dt_estim)**2._rp)
    etan_ss = -2._rp*(log(en))*meffn_ss/(Nstretch*dt_estim)
    etat_ss = -2._rp*(log(et))*mefft_ss/(Nstretch*dt_estim)
    muc_ss = muc
    psi_crit_ss = 7._rp/2._rp*(1._rp+en)/(1._rp+et)*muc_ss
    !
    ! sphere/wall
    !
    colthr_pw = 0.001_rp*radius
    meffn_sw = ratiorho*volp
    mefft_sw = 2._rp/7._rp*meffn_sw
    kn_sw = (pi**2._rp + abs(log(en))**2._rp)*meffn_sw/((Nstretch*dt_estim)**2._rp)
    kt_sw = (pi**2._rp + abs(log(et))**2._rp)*mefft_sw/((Nstretch*dt_estim)**2._rp)
    etan_sw = -2._rp*(log(en))*meffn_sw/(Nstretch*dt_estim)
    etat_sw = -2._rp*(log(et))*mefft_sw/(Nstretch*dt_estim)
    muc_sw = muc
    psi_crit_sw = 7._rp/2._rp*(1._rp+en)/(1._rp+et)*muc_sw
  end subroutine read_particle_input
#endif
end module prt_mod_param
