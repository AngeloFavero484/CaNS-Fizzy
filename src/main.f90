! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
!
!        CCCCCCCCCCCCC                    NNNNNNNN        NNNNNNNN    SSSSSSSSSSSSSSS
!     CCC::::::::::::C                    N:::::::N       N::::::N  SS:::::::::::::::S
!   CC:::::::::::::::C                    N::::::::N      N::::::N S:::::SSSSSS::::::S
!  C:::::CCCCCCCC::::C                    N:::::::::N     N::::::N S:::::S     SSSSSSS
! C:::::C       CCCCCC   aaaaaaaaaaaaa    N::::::::::N    N::::::N S:::::S
!C:::::C                 a::::::::::::a   N:::::::::::N   N::::::N S:::::S
!C:::::C                 aaaaaaaaa:::::a  N:::::::N::::N  N::::::N  S::::SSSS
!C:::::C                          a::::a  N::::::N N::::N N::::::N   SS::::::SSSSS
!C:::::C                   aaaaaaa:::::a  N::::::N  N::::N:::::::N     SSS::::::::SS
!C:::::C                 aa::::::::::::a  N::::::N   N:::::::::::N        SSSSSS::::S
!C:::::C                a::::aaaa::::::a  N::::::N    N::::::::::N             S:::::S
! C:::::C       CCCCCC a::::a    a:::::a  N::::::N     N:::::::::N             S:::::S
!  C:::::CCCCCCCC::::C a::::a    a:::::a  N::::::N      N::::::::N SSSSSSS     S:::::S
!   CC:::::::::::::::C a:::::aaaa::::::a  N::::::N       N:::::::N S::::::SSSSSS:::::S
!     CCC::::::::::::C  a::::::::::aa:::a N::::::N        N::::::N S:::::::::::::::SS
!        CCCCCCCCCCCCC   aaaaaaaaaa  aaaa NNNNNNNN         NNNNNNN  SSSSSSSSSSSSSSS
!-------------------------------------------------------------------------------------
! CaNS -- Canonical Navier-Stokes Solver
!-------------------------------------------------------------------------------------
program cans
#if defined(_DEBUG)
  use, intrinsic :: iso_fortran_env, only: compiler_version,compiler_options
#endif
  use, intrinsic :: iso_c_binding  , only: C_PTR
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan
  use mpi
  use decomp_2d
  use mod_bound          , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,ierr
  use mod_correc         , only: correc
  use mod_fft            , only: fftini,fftend
  use mod_fillps         , only: fillps
  use mod_forcing        , only: lscale_forcing
  use mod_initflow       , only: initflow,initscal
  use mod_initgrid       , only: initgrid
  use mod_initmpi        , only: initmpi
  use mod_initsolver     , only: initsolver
  use mod_load           , only: load_one
  use mod_rk             , only: tm => rk,tm_scal => rk_scal,tm_2fl => rk_2fl
  use mod_output         , only: out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  use mod_param          , only: small, &
                                 nb,is_bound,cbcvel,bcvel,cbcpre,bcpre,cbcsca,bcsca,cbcpsi,bcpsi,cbcnor,bcnor,cbccur,bccur, &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 nstep,time_max,tw_max,stop_type,restart,is_overwrite_save,nsaves_max, &
                                 datadir,   &
                                 is_solve_ns,is_track_interface, &
                                 cfl,dtmax,dt_f, &
                                 inivel,inisca,inipsi, &
                                 is_wallturb,is_forced_hit, &
                                 dims, &
                                 gtype,gr, &
                                 bforce,ssource, &
                                 ng,l,dl,dli,pi, &
                                 read_input, &
                                 rho0,rho12,mu12,theta,sigma,gacc,ka12,cp12,beta12, &
                                 psi_thickness_factor, &
                                 acdi_gam_factor,acdi_gam_min, &
                                 vof_thinc_beta, &
                                 max_pseudo_iter,dtau_cfl
  use mod_rotnorm        , only: rot_norm
  use mod_extend         , only: compute_uextend, advect_vof_upwind
#if 1
  use mod_sanity         , only: test_sanity_input
#endif
#if !defined(_INTERFACE_CAPTURING_VOF)
  use mod_acdi           , only: acdi_set_epsilon,acdi_set_gamma,acdi_cmpt_phi
#elif defined(_SDF_NORMALS)
  use mod_vof_thinc_qq   , only: vof_thinc_cmpt_phi
#endif
  use mod_two_fluid      , only: init2fl,cmpt_norm_curv => cmpt_norm_curv_youngs
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mod_solver_vc      , only: solver_vc
#endif
#if !defined(_OPENACC)
  use mod_solver         , only: solver
#else
  use mod_solver_gpu     , only: solver => solver_gpu
  use mod_workspaces     , only: init_wspace_arrays,set_cufft_wspace
  use mod_common_cudecomp, only: istream_acc_queue_1
#endif
  use mod_timer          , only: timer_tic,timer_toc,timer_print
  use mod_updatep        , only: updatep
  use mod_utils          , only: bulk_mean_12_stag
  !@acc use mod_utils    , only: device_memory_footprint
  use mod_types
#if defined(_PARTICLE)
  use mod_param                    , only: is_ibm !,nh_wide
  use prt_mod_kernel               , only: kerneltest
#if !defined(_EULER)
  use prt_mod_common               , only: prt_InitMemo,uf,vf,wf,solidity
#else
  use prt_mod_common               , only: prt_InitMemo,uphase,vphase,wphase,solidity,alphac
#endif
  use prt_mod_initparticles        , only: initparticles
  use prt_mod_initeul              , only: initeul
  use prt_mod_initvof              , only: initvof
  use prt_mod_param                , only: radius,rho_s,read_particle_input
  use prt_mod_coordsfp             , only: coordsfp
  use prt_mod_intgr_over_sphere    , only: intgr_over_sphere
#if !defined(_EULER)
  use prt_mod_interp_spread        , only: eulr2lagr,lagr2eulr
  use prt_mod_forcing              , only: complagrforces,updtlagrforces
#else
  use prt_mod_eulint               , only: eulint
#endif
  use prt_mod_phase_indicator      , only: phase_indicator
  use prt_mod_output               , only: outpart
  use prt_mod_loadpart             , only: loadpart
  use prt_mod_intgr_nwtn_eulr      , only: intgr_nwtn_eulr
#endif
  implicit none
  integer , dimension(3) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,pp,pn,po
  real(rp), allocatable, dimension(:,:,:) :: u_ext,v_ext,w_ext
  real(rp), dimension(3) :: rho_av
  logical , dimension(3) :: is_cmpt_rho_av
  real(rp) :: rho_ratio,mu_ratio,Ga,Ar,Mo,Eo,Oh,radius_sphere,cap_length, dt_cap
  real(rp), dimension(3) :: Fs,Fstot,Fstot_old
  real(rp)  :: F_ibm,F_inertia,F_w,F_buoy,F_cap
#if !defined(_OPENACC)
  type(C_PTR), dimension(2,2) :: arrplanp
#else
  integer    , dimension(2,2) :: arrplanp
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  real(rp) :: normfftp
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  type(rhs_bound) :: rhsbp
  real(rp) :: dt,dto,dt_r,dti,dt_cfl,dtrk,dtrki,time,divtot,divmax
  real(rp) :: gam,seps
  real(rp) :: mass,mass_tot
  integer :: irk,istep
  real(rp) :: dtau
  integer :: iter
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g, &
                                         grid_vol_ratio_c,grid_vol_ratio_f
  real(rp), dimension(3) :: dpdl
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(42) :: var
#if defined(_TIMING)
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  real(rp) :: twi,tw
  !
  integer  :: savecounter
  character(len=7  ) :: fldnum
  character(len=4  ) :: chkptnum
  character(len=100) :: filename,fexts(10)
  integer :: k,kk
  integer :: i,j
  logical :: is_done,kill
  real(rp), dimension(2) :: tm_coeff
  !
  integer, parameter :: csv_unit = 5555
  !
  real(rp), allocatable, dimension(:,:,:) :: s
  !
  ! two-fluid solver specific
  !
  real(rp), allocatable, dimension(:,:,:) :: psi,psio,phi,kappa,normx,normy,normz, &
                                             psiflx_x,psiflx_y,psiflx_z,fx_old,fy_old,fz_old
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
  real(rp), allocatable, dimension(:,:,:) :: surfx_n,surfy_n,surfz_n, &
                                             surfx_o,surfy_o,surfz_o
#endif
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  !
  ! forces_data.csv is opened, written and closed by rank 0 only.
  ! The rows come from whichever rank currently masters the particle, and that
  ! rank changes as the particle centre crosses a pencil boundary; the data is
  ! therefore reduced onto rank 0 in intgr_nwtn_eulr before being written here.
  ! (Opening on every rank, as was done before, gave each rank its own handle
  ! and its own file offset on the same path, so at every change of master the
  ! new writer resumed just after the header and overwrote the file from the
  ! top. The open also has to come after MPI_INIT for myid to be defined.)
  !
  if(myid == 0) then
    open(unit=csv_unit, file='forces_data.csv', status='replace', action='write')
    write(csv_unit, '(A)') "time,F_cap_ibm,F_ibm,F_inertia,F_w,F_bouy,F_cap,ep_z,ep_w"
    flush(csv_unit)
  endif
  !
  ! read parameter file
  !
  call read_input(myid)
#if defined(_PARTICLE)
  call read_particle_input(myid)
#endif
  !
  ! a fresh (non-restart) run must not append to time.out/log_visu_*.out left
  ! over from a previous run in the same datadir: out0d and write_log_output
  ! always open with position='append', so stale rows/entries from before
  ! would otherwise survive into the new run's log.
  !
  if(.not.restart) then
    if(myid == 0) call execute_command_line('rm -f '//trim(datadir)//'time.out '//trim(datadir)//'log_visu_*.out')
  end if
  !
  ! initialize MPI/OpenMP
  !
  call initmpi(ng,dims,cbcpre,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,nb,is_bound)
  twi = MPI_WTIME()
  savecounter = 0
  !
  ! allocate variables
  !
#if defined(_PARTICLE)
  is_ibm = .false.
  !
  call prt_InitMemo(n,lo)
#endif
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           u_ext( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v_ext( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w_ext( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
#if !defined(_CONSTANT_COEFFS_POISSON)
  allocate(po,mold=pp)
  po(:,:,:) = 0._rp
#else
  pp(:,:,:) = 0._rp
  allocate(pn,mold=pp)
  allocate(po,mold=pp)
  pn(:,:,:) = pp(:,:,:)
  po(:,:,:) = pn(:,:,:)
#endif
#if defined(_SCALAR)
  allocate(s,mold=pp)
#endif
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  allocate(dzc( 0:n(3)+1), &
           dzf( 0:n(3)+1), &
           zc(  0:n(3)+1), &
           zf(  0:n(3)+1), &
           dzci(0:n(3)+1), &
           dzfi(0:n(3)+1))
  allocate(dzc_g( 0:ng(3)+1), &
           dzf_g( 0:ng(3)+1), &
           zc_g(  0:ng(3)+1), &
           zf_g(  0:ng(3)+1), &
           dzci_g(0:ng(3)+1), &
           dzfi_g(0:ng(3)+1))
  allocate(grid_vol_ratio_c,mold=dzc)
  allocate(grid_vol_ratio_f,mold=dzf)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
  allocate(psi,kappa,normx,normy,normz,mold=pp)
  allocate(psio,mold=pp)
  allocate(fx_old,fy_old,fz_old,mold=pp)
  Fs(:) = 0._rp
  Fstot(:) = 0._rp
  Fstot_old(:) = 0._rp
#if defined(_SDF_NORMALS)
  allocate(phi,mold=pp)
#endif
  allocate(psiflx_x,psiflx_y,psiflx_z,mold=pp)
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
  allocate(surfx_n(n(1),n(2),n(3)),surfy_n(n(1),n(2),n(3)),surfz_n(n(1),n(2),n(3)), &
           surfx_o(n(1),n(2),n(3)),surfy_o(n(1),n(2),n(3)),surfz_o(n(1),n(2),n(3)))
  surfx_n(:,:,:) = 0._rp
  surfy_n(:,:,:) = 0._rp
  surfz_n(:,:,:) = 0._rp
  surfx_o(:,:,:) = 0._rp
  surfy_o(:,:,:) = 0._rp
  surfz_o(:,:,:) = 0._rp
#endif
#if defined(_DEBUG)
  if(myid == 0) print*, 'This executable of CaNS was built with compiler: ', compiler_version()
  if(myid == 0) print*, 'Using the options: ', compiler_options()
  block
    character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version
    integer :: ilen
    call MPI_GET_LIBRARY_VERSION(mpi_version,ilen,ierr)
    if(myid == 0) print*, 'MPI Version: ', trim(mpi_version)
  end block
  if(myid == 0) print*, ''
#endif
!Dimensionless parameters
if (myid == 0) then
!radius_sphere=0.15
radius_sphere=radius
rho_ratio=rho12(2)/rho12(1)
mu_ratio=mu12(2)/mu12(1)
Eo=(abs(rho12(1)-rho12(2))*max(abs(gacc(1)),abs(gacc(2)),abs(gacc(3)))*((2*radius_sphere)**2))/sigma
!Goccia
Ga=(rho12(1)*abs(rho12(1)-rho12(2))*radius_sphere*max(abs(gacc(1)),abs(gacc(2)),abs(gacc(3))))/(mu12(1)**2)
!Bolla
!Ga=(rho12(2)*abs(rho12(1)-rho12(2))*radius_sphere*max(abs(gacc(1)),abs(gacc(2)),abs(gacc(3))))/(mu12(2)**2)
Ar=sqrt(Ga)
Oh = mu12(1)/sqrt(2*radius_sphere*sigma*rho12(1))
Mo=(Eo**3)/(Ar**2)
cap_length = sqrt(sigma/(rho12(1)*max(abs(gacc(1)),abs(gacc(2)),abs(gacc(3)))))

dt_cap = sqrt(((rho12(1)+rho12(2))*dl(1)*dl(2)*dl(3))/(4*pi*sigma))

PRINT *, "Dimensionless parameters"
PRINT *, "Bubble initial radius", radius_sphere
PRINT *, "rho_ratio", rho_ratio 
PRINT *, "mu_ratio", mu_ratio 
PRINT *, "Ar", Ar
PRINT *, "Eo", Eo 
PRINT *, "Oh", Oh
PRINT *, "Mo", Mo
PRINT *, "Capillary length", cap_length
PRINT *, "Capillary max Dt", dt_cap

endif

  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, '*** Beginning of simulation ***'
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, ''
  call initgrid(gtype,ng(3),gr,l(3),dzc_g,dzf_g,zc_g,zf_g)
  if(myid == 0) then
    open(99,file=trim(datadir)//'grid.bin',action='write',form='unformatted',access='stream',status='replace')
    write(99) dzc_g(1:ng(3)),dzf_g(1:ng(3)),zc_g(1:ng(3)),zf_g(1:ng(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=0,ng(3)+1
      write(99,*) 0.,zf_g(kk),zc_g(kk),dzf_g(kk),dzc_g(kk)
    end do
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3)
      write(99,*) l(1),l(2),l(3)
    close(99)
  end if
  !$acc enter data copyin(lo,hi,n) async
  !$acc enter data copyin(bforce,gacc,dl,dli,l) async
  !$acc enter data copyin(zc_g,zf_g,dzc_g,dzf_g) async
  !$acc enter data create(zc,zf,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g) async
  !
  !$acc enter data copyin(rho12,mu12,ka12,cp12,beta12) async
  !
  !$acc parallel loop default(present) private(k) async
  do kk=lo(3)-1,hi(3)+1
    k = kk-(lo(3)-1)
    zc( k) = zc_g(kk)
    zf( k) = zf_g(kk)
    dzc(k) = dzc_g(kk)
    dzf(k) = dzf_g(kk)
    dzci(k) = dzc(k)**(-1)
    dzfi(k) = dzf(k)**(-1)
  end do
  !$acc kernels default(present) async
  dzci_g(:) = dzc_g(:)**(-1)
  dzfi_g(:) = dzf_g(:)**(-1)
  !$acc end kernels
  !$acc enter data create(grid_vol_ratio_c,grid_vol_ratio_f) async
  !$acc kernels default(present) async
  grid_vol_ratio_c(:) = dl(1)*dl(2)*dzc(:)/(l(1)*l(2)*l(3))
  grid_vol_ratio_f(:) = dl(1)*dl(2)*dzf(:)/(l(1)*l(2)*l(3))
  !$acc end kernels
  !$acc update self(zc,zf,dzc,dzf,dzci,dzfi) async
  !$acc exit data copyout(zc_g,zf_g,dzc_g,dzf_g,dzci_g,dzfi_g) async ! not needed on the device
  !$acc wait
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre)
  !
  ! initialize Poisson solver
  !
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
                  lambdaxyp,['c','c','c'],ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
  !$acc enter data copyin(lambdaxyp,ap,bp,cp) async
  !$acc enter data copyin(rhsbp,rhsbp%x,rhsbp%y,rhsbp%z) async
  !$acc wait
#if defined(_OPENACC)
  !
  ! determine workspace sizes and allocate the memory
  !
  call init_wspace_arrays()
  call set_cufft_wspace(pack(arrplanp,.true.),istream_acc_queue_1)
  if(myid == 0) print*,'*** Device memory footprint (Gb): ', &
                  device_memory_footprint(n,n_z)/(1._sp*1024**3), ' ***'
#endif
#if defined(_DEBUG_SOLVER)
  call test_sanity_solver(ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,dli,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g, &
                          nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
#endif
  !
#if !defined(_INTERFACE_CAPTURING_VOF)
  call acdi_set_epsilon(dl,dzfi,psi_thickness_factor,seps)
#endif
  !
  fexts(1) = 'u'
  fexts(2) = 'v'
  fexts(3) = 'w'
  fexts(4) = 'p'
  fexts(5) = 'psi'
  fexts(6) = 'fx_old'
  fexts(7) = 'fy_old'
  fexts(8) = 'fz_old'
#if defined(_PARTICLE)
  fexts(9) = 'alphac'
#endif
!  fexts(6) = 'normx'
!  fexts(7) = 'normy'
!  fexts(8) = 'normz'
!  fexts(9) = 'kappa'
#if defined(_SCALAR)
  fexts(10) = 's'
#endif
  if(.not.restart) then
    istep = 0
    time = 0.
    fx_old(:,:,:) = 0._rp
    fy_old(:,:,:) = 0._rp
    fz_old(:,:,:) = 0._rp
    call initflow(inivel,bcvel,ng,lo,l,dl,zc,zf,dzc,dzf,rho12(2),mu12(2),bforce,is_wallturb,time,u,v,w,p)
#if defined(_SCALAR)
    call initscal(inisca,bcsca,ng,lo,l,dl,dzf,zc,s)
#endif
#if defined(_PARTICLE)
    if(myid == 0) print*, '*** Particle initialization  ***'
    !
#if !defined(_EULER)
    if (myid == 0) then
      call kerneltest(sumk)
      write(6,*) 'Integral over kernel = ', sumk, ' (~ 1)'
    endif
#endif
    call initparticles
    call initeul(n)
    if(myid == 0) print*, '*** Particle initial condition succesfully set  ***'
    if(myid == 0) print*, ''
#endif
    !
    call init2fl(inipsi,cbcpsi,seps,lo,hi,l,dl,zc_g,psi)
    if(myid == 0) print*, '*** Initial condition succesfully set ***'
  else
#if defined(_PARTICLE)
    call loadpart('r')
#endif
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(1))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,u,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(2))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,v,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(3))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,w,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(4))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,p,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(5))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,psi,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(6))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,fx_old,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(7))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,fy_old,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(8))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,fz_old,time,istep)
!    call load_one('r',trim(datadir)//'fld_'//trim(fexts(9))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,kappa,time,istep)
#if defined(_PARTICLE)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(9))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,alphac,time,istep)
#endif
#if defined(_SCALAR)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(10))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,s,time,istep)
#endif
    if(myid == 0) print*, '*** Checkpoints loaded at time = ', time, 'time step = ', istep, '. ***'
  end if
  !$acc enter data copyin(u,v,w,p,pp,pn,po) async
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
  !$acc enter data copyin(surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o) async
#endif
  !$acc wait
  call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
#if defined(_PARTICLE)
#if defined(_EULER)
  call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,alphac)
#endif
#endif
#if defined(_CONSTANT_COEFFS_POISSON)
  !$acc kernels default(present) async(1)
  pn(:,:,:) =  p(:,:,:)
  po(:,:,:) = pn(:,:,:)
  !$acc end kernels
#endif
#if defined(_SCALAR)
  !$acc enter data copyin(s) async(1)
  call boundp(cbcsca,n,bcsca,nb,is_bound,dl,dzc,s)
  !$acc wait
#endif
  !$acc enter data copyin(psi) create(kappa,normx,normy,normz) async(1)
#if defined(_SDF_NORMALS)
  !$acc enter data create(phi) async(1)
#endif
  !$acc enter data create(psio,psiflx_x,psiflx_y,psiflx_z) async(1)
  call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,psi)
  !$acc wait
  !
#if defined(_SDF_NORMALS)
#if !defined(_INTERFACE_CAPTURING_VOF)
  call acdi_cmpt_phi(n,seps,psi,phi)
#else
  call vof_thinc_cmpt_phi(n,vof_thinc_beta,psi,phi)
#endif
#endif
#if defined(_SDF_NORMALS)
  call cmpt_norm_curv(n,dli,dzci,dzfi,phi,normx,normy,normz,kappa)
#else
  call cmpt_norm_curv(n,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
#endif
  call boundp(cbccur,n,bccur,nb,is_bound,dl,dzc,kappa)
  call boundp(cbcnor(:,:,1),n,bcnor(:,:,1),nb,is_bound,dl,dzc,normx)
  call boundp(cbcnor(:,:,2),n,bcnor(:,:,2),nb,is_bound,dl,dzc,normy)
  call boundp(cbcnor(:,:,3),n,bcnor(:,:,3),nb,is_bound,dl,dzc,normz)

  ! pseudo-time step for the contact-line relaxation below; u_ext is a unit
  ! vector, so dtau_cfl is a CFL number on the smallest cell size
  dtau = dtau_cfl / maxval(dli(1:3))
  do iter = 1, max_pseudo_iter
    u_ext=0
    v_ext=0
    w_ext=0
    call compute_uextend(n, theta, normx, normy, normz, u_ext, v_ext, w_ext)
    call advect_vof_upwind(n, dli, dtau, u_ext, v_ext, w_ext, psi)
    call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,psi)
    call cmpt_norm_curv(n,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
    call boundp(cbcnor(:,:,1),n,bcnor(:,:,1),nb,is_bound,dl,dzc,normx)
    call boundp(cbcnor(:,:,2),n,bcnor(:,:,2),nb,is_bound,dl,dzc,normy)
    call boundp(cbcnor(:,:,3),n,bcnor(:,:,3),nb,is_bound,dl,dzc,normz)
    call boundp(cbcpsi,n,bcpre,nb,is_bound,dl,dzc,kappa)
  end do
  call rot_norm(n,dli,dzci,psi,theta,is_bound,normx,normy,normz,kappa,Fs)
  Fstot_old=Fstot
  call MPI_ALLREDUCE(Fs, Fstot, 3, MPI_REAL_RP, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (myid==0) then
    PRINT *, "Fstot", Fstot
  end if
  !
#if !defined(_INTERFACE_CAPTURING_VOF)
  call acdi_set_gamma(n,acdi_gam_factor,u,v,w,gam)
  gam = max(gam,acdi_gam_min)
  if(myid == 0) print*, 'ACDI parameters. Gamma: ', gam, 'Epsilon: ', seps
#endif
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !$acc wait
  !$acc update self(u,v,w,p,psi,kappa,s)
  if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
#include "out1d.h90"
  end if
  if(iout2d > 0.and.mod(istep,max(iout2d,1)) == 0) then
#include "out2d.h90"
  end if
  if(iout3d > 0.and.mod(istep,max(iout3d,1)) == 0) then
#include "out3d.h90"
#if defined(_PARTICLE)
    include 'prt_out.h90'
#endif
  end if
#if defined(_PARTICLE)
  !$acc enter data create(uf,vf,wf)
  if (myid .eq. 0) write(6,*) 'Solid volume fraction of particles = ',solidity
  !
  call coordsfp(n)
  !
  if(.not.restart) then
    call intgr_over_sphere(1,n,psi,u,v,w)
    call intgr_over_sphere(2,n,psi,u,v,w)
    call intgr_over_sphere(3,n,psi,u,v,w)
  endif
#endif
  !
  call chkdt(n,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dt_cfl,gam,seps,ka12,cp12)
  dt = min(cfl*dt_cfl,dtmax)
  if(dt_f > 0.) then
    if(dt_f > dt) then
      if(myid == 0) print*, 'WARNING: fixed time step exceeds estimated stability limit.'
      if(myid == 0) print*, 'dt_f = ', dt_f, 'Estimated stable dt = ', dt
    end if
    dt = dt_f
  end if
  if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
  dto = dt
  dti = 1./dt
  kill = .false.
  !
!  mass=0
!  do k=1,n(3)
!    do j=1,n(2)
!      do i=1,n(1)
!        mass = mass + psi(i,j,k) !(rho12(1)*psi(i,j,k)+rho12(2)*(1-psi(i,j,k)))*dl(1)*dl(2)*dl(3)
!      end do
!    end do
!  end do
!  PRINT *, "Mass", mass
  ! main loop
  !
  if(myid == 0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)
#if defined(_TIMING)
    !$acc wait(1)
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time = time + dt
    if(myid == 0) print*, 'Time step #', istep, 'Time = ', time
! Runge-Kutta
!    do irk=1,3
!      tm_coeff(:) = rkcoeff(:,irk)
!      dtrk = sum(rkcoeff(:,irk))*dt
!Adams-Bashfort
    do irk=1,1
      tm_coeff(:) = [2._rp+dt/dto,-dt/dto]/2._rp
      dtrk = dt
      !
      dtrki = dtrk**(-1)
      dt_r = dtrk/dto
      !
      ! phase field update
      !
      !$acc kernels default(present) async(1)
      psio(:,:,:)   = psi(:,:,:)
      !$acc end kernels
      if(is_track_interface) then
        uphase(1:n(1),1:n(2),1:n(3)) = u(1:n(1),1:n(2),1:n(3))
        vphase(1:n(1),1:n(2),1:n(3)) = v(1:n(1),1:n(2),1:n(3))
        wphase(1:n(1),1:n(2),1:n(3)) = w(1:n(1),1:n(2),1:n(3))
        call initvof(n,u,v,w,uphase,vphase,wphase)
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,uphase,vphase,wphase)
        call tm_2fl(tm_coeff,n,dli,dzci,dzfi,dt,gam,seps,vof_thinc_beta,uphase,vphase,wphase, &
        normx,normy,normz,phi,psi,psiflx_x,psiflx_y,psiflx_z)
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
        call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,psi)
#if !defined(_INTERFACE_CAPTURING_VOF)
        call acdi_cmpt_phi(n,seps,psi,phi)
#else
        call vof_thinc_cmpt_phi(n,vof_thinc_beta,psi,phi)
#endif
#if defined(_SDF_NORMALS)
        call cmpt_norm_curv(n,dli,dzci,dzfi,phi,normx,normy,normz,kappa)
#else
        call cmpt_norm_curv(n,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
#endif
        call boundp(cbccur,n,bccur,nb,is_bound,dl,dzc,kappa)
        call boundp(cbcnor(:,:,1),n,bcnor(:,:,1),nb,is_bound,dl,dzc,normx)
        call boundp(cbcnor(:,:,2),n,bcnor(:,:,2),nb,is_bound,dl,dzc,normy)
        call boundp(cbcnor(:,:,3),n,bcnor(:,:,3),nb,is_bound,dl,dzc,normz)
        !
        ! pseudo-time step for the contact-line relaxation below; u_ext is a
        ! unit vector, so dtau_cfl is a CFL number on the smallest cell size
        dtau = dtau_cfl / maxval(dli(1:3))
        !
        call initeul(n)
        call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,alphac)
        do iter = 1, max_pseudo_iter
          u_ext=0
          v_ext=0
          w_ext=0
          call compute_uextend(n, theta, normx, normy, normz, u_ext, v_ext, w_ext)
          call advect_vof_upwind(n, dli, dtau, u_ext, v_ext, w_ext, psi)
          call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,psi)
          call cmpt_norm_curv(n,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
          call boundp(cbcnor(:,:,1),n,bcnor(:,:,1),nb,is_bound,dl,dzc,normx)
          call boundp(cbcnor(:,:,2),n,bcnor(:,:,2),nb,is_bound,dl,dzc,normy)
          call boundp(cbcnor(:,:,3),n,bcnor(:,:,3),nb,is_bound,dl,dzc,normz)
          call boundp(cbcpsi,n,bcpre,nb,is_bound,dl,dzc,kappa)
        end do
        call rot_norm(n,dli,dzci,psi,theta,is_bound,normx,normy,normz,kappa,Fs)
        Fstot_old=Fstot
        call MPI_ALLREDUCE(Fs, Fstot, 3, MPI_REAL_RP, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (myid==0) then
          PRINT *, "Fstot", Fstot
        end if
      else
        call initeul(n)
        call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,alphac)
        !$acc kernels default(present) async(1)
        psiflx_x(:,:,:) = 0.
        psiflx_y(:,:,:) = 0.
        psiflx_z(:,:,:) = 0.
        !$acc end kernels
      end if
#if defined(_SCALAR)
    call tm_scal(tm_coeff,n,dli,dzci,dzfi,dt,ssource,rho12,ka12,cp12,psi,u,v,w,psio,psiflx_x,psiflx_y,psiflx_z,s)
    call boundp(cbcsca,n,bcsca,nb,is_bound,dl,dzc,s)
#endif
      if(.not.is_solve_ns) then
        call initflow(inivel,bcvel,ng,lo,l,dl,zc,zf,dzc,dzf,rho12(2),mu12(2),bforce,is_wallturb,time,u,v,w,p)
        !$acc wait(1)
        !$acc update device(u,v,w,p) async(1)
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
      else
        rho_av(:) = 0.
        is_cmpt_rho_av(:) = (abs(gacc(:)) > 0.) .and. cbcpre(0,:)//cbcpre(1,:) == 'PP'
        if(is_cmpt_rho_av(1)) &
          call bulk_mean_12_stag(n,1,grid_vol_ratio_c,psi,rho12,rho_av(1))
        if(is_cmpt_rho_av(2)) &
          call bulk_mean_12_stag(n,2,grid_vol_ratio_c,psi,rho12,rho_av(2))
        if(is_cmpt_rho_av(3)) &
          call bulk_mean_12_stag(n,3,grid_vol_ratio_f,psi,rho12,rho_av(3))
        !
        call tm(tm_coeff,n,dli,dzci,dzfi,dt,dt_r, &
                bforce,gacc,sigma,rho_av,rho12,mu12,beta12,rho0,psi,kappa,p,pn,po,s, &
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
                surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o, &
#endif
                psio,psiflx_x,psiflx_y,psiflx_z,u,v,w)
        !
        if(is_forced_hit) then
          call lscale_forcing(2,lo,hi,0.5_rp,dtrk,l,dl,zc,zf,u,v,w)
        end if
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
#if defined(_PARTICLE)
      !
      !$acc wait
      !$acc update self(u,v,w)
      is_ibm = .true.
      !
#if !defined(_EULER)
      ibmiter = 0
      !
      !$acc kernels default(present) async(1)
      uf(1:n(1),1:n(2),1:n(3)) = u(1:n(1),1:n(2),1:n(3))
      vf(1:n(1),1:n(2),1:n(3)) = v(1:n(1),1:n(2),1:n(3))
      wf(1:n(1),1:n(2),1:n(3)) = w(1:n(1),1:n(2),1:n(3))
#endif
      !$acc end kernels
      !
#if !defined(_EULER)
!      do while(ibmiter <= 2)
!      !  call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,uf,vf,wf)
!        !$acc wait
!        !$acc update self(uf,vf,wf)
!        call eulr2lagr(n,lo,uf,vf,wf)
!        !
!        if (ibmiter == 0) then
!          call complagrforces(dtrki,irk)
!        else
!          call updtlagrforces(dtrki,dti,ibmiter)
!        endif
!        !
!        call lagr2eulr(dtrk,n,lo,u,v,w,uf,vf,wf)
!        !
!        !$acc update device(uf,vf,wf)
!        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,uf,vf,wf)
!        ibmiter = ibmiter+1
!      enddo
#else
      call eulint(tm_coeff,dt,irk,n,psi,psio,kappa,fx_old,fy_old,fz_old,u,v,w)
      call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,alphac)
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
#endif
      !
#if !defined(_EULER)
      !$acc kernels default(present) async(1)
      u(1:n(1),1:n(2),1:n(3)) = uf(1:n(1),1:n(2),1:n(3))
      v(1:n(1),1:n(2),1:n(3)) = vf(1:n(1),1:n(2),1:n(3))
      w(1:n(1),1:n(2),1:n(3)) = wf(1:n(1),1:n(2),1:n(3))
      !$acc end kernels
#endif
      !
      is_ibm = .false.
      !
#endif
      !$acc kernels default(present) async(1)
      pp(:,:,:) = p(:,:,:)
      !$acc end kernels
      call fillps(n,dli,dzfi,dtrki,rho0,u,v,w,p)
#if defined(_CONSTANT_COEFFS_POISSON)
      call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,p)
      call solver(n,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],p)
#else
      call solver_vc(ng,lo,hi,cbcpre,bcpre,dli,dzci,dzfi,is_bound,rho12,psi,p,po)
#endif
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
      call correc(n,dli,dzci,rho0,rho12,dtrk,p,psi,u,v,w)
      call bounduvw(cbcvel,n,bcvel,nb,is_bound,.true.,dl,dzc,dzf,u,v,w)
      call updatep(pp,p)
      call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
#if defined(_PARTICLE)
      call intgr_nwtn_eulr(n,l,dl,dli,dt,time,tm_coeff,istep,psi,u,v,w,Fstot,Fstot_old,F_ibm,F_inertia,F_w,F_buoy,F_cap)
#endif
      end if
    end do
#if defined(_CONSTANT_COEFFS_POISSON)
    !$acc kernels default(present) async(1)
    po(:,:,:) = pn(:,:,:)
    pn(:,:,:) =  p(:,:,:)
    !$acc end kernels
#endif
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep >= nstep   ) is_done = is_done.or..true.
    end if
    if(stop_type(2)) then ! maximum simulation time reached
      if(time  >= time_max) is_done = is_done.or..true.
    end if
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600.
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    end if
    dto = dt
!  mass=0
!  do k=1,n(3)
!    do j=1,n(2)
!      do i=1,n(1)
!        mass = mass + rho12(1)*psi(i,j,k)*(1-alphac(i,j,k)) + rho12(2)*(1-psi(i,j,k))*(1-alphac(i,j,k)) + &
!               rho_s*alphac(i,j,k) !(rho12(1)*psi(i,j,k)+rho12(2)*(1-psi(i,j,k)))*dl(1)*dl(2)*dl(3)
!      end do
!    end do
!  end do
!  call MPI_ALLREDUCE(mass, mass_tot, 1, MPI_REAL_RP, MPI_SUM, MPI_COMM_WORLD, ierr)
!  if (myid==0) then 
!    PRINT *, "Mass", mass_tot
!  end if
#if !defined(_INTERFACE_CAPTURING_VOF)
    if(mod(istep,1) == 0) then
      call acdi_set_gamma(n,acdi_gam_factor,u,v,w,gam)
      gam = max(gam,acdi_gam_min)
    end if
#endif
    if(mod(istep,icheck) == 0) then
      if(myid == 0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dt_cfl,gam,seps,ka12,cp12)
      dt = min(cfl*dt_cfl,dtmax)
      if(dt_f > 0.) then
        if(dt_f > dt) then
          if(myid == 0) print*, 'WARNING: fixed time step exceeds estimated stability limit.'
          if(myid == 0) print*, 'dt_f = ', dt_f, 'Estimated stable dt = ', dt
          is_done = .true.
          kill = .true.
        else
          dt = dt_f
        end if
      end if
      if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
      if(dt_cfl < small) then
        if(myid == 0) print*, 'ERROR: time step is too small.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
      dti = 1./dt
      call chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
      if(myid == 0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
#if !defined(_MASK_DIVERGENCE_CHECK)
      if(divmax > small.or.is_nan(divtot)) then
        if(myid == 0) print*, 'ERROR: maximum divergence is too large.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
#endif
    end if
    !
    ! output routines below
    !
    if(iout0d > 0.and.mod(istep,max(iout0d,1)) == 0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
#if !defined(_INTERFACE_CAPTURING_VOF)
      var(1) = 1.*istep
      var(2) = time
      var(3) = gam
      var(4) = seps
      call out0d(trim(datadir)//'log_acdi.out',4,var)
#endif
    end if
    write(fldnum,'(i7.7)') istep
    if(iout1d > 0.and.mod(istep,max(iout1d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa,s)
#include "out1d.h90"
    end if
    if(iout2d > 0.and.mod(istep,max(iout2d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa,s)
#include "out2d.h90"
    end if
    if(iout3d > 0.and.mod(istep,max(iout3d,1)) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa,s)
#include "out3d.h90"
#if defined(_PARTICLE)
      include 'prt_out.h90'
#endif
    end if
    if(isave > 0.and.((mod(istep,max(isave,1)) == 0).or.(is_done.and..not.kill))) then
      if(is_overwrite_save) then
        filename = 'fld'
      else
        filename = 'fld_'//fldnum
        if(nsaves_max > 0) then
          if(savecounter >= nsaves_max) savecounter = 0
          savecounter = savecounter + 1
          write(chkptnum,'(i4.4)') savecounter
          filename = 'fld_'//chkptnum
          var(1) = 1.*istep
          var(2) = time
          var(3) = 1.*savecounter
          call out0d(trim(datadir)//'log_checkpoints.out',3,var)
        end if
      end if
      !$acc wait
      !$acc update self(u,v,w,p,psi,s)
   call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(1))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,u,time,istep)
   call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(2))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,v,time,istep)
   call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(3))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,w,time,istep)
   call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(4))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,p,time,istep)
   call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(5))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,psi,time,istep)
   call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(6))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,fx_old,time,istep)
   call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(7))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,fy_old,time,istep)
   call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(8))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,fz_old,time,istep)
#if defined(_PARTICLE)
    call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(9))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,alphac,time,istep)
#endif
!      call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(9))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,kappa,time,istep)
#if defined(_SCALAR)
      call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(10))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,s,time,istep)
#endif
#if defined(_PARTICLE)
      call loadpart('w')
#endif
      if(.not.is_overwrite_save) then
        !
        ! fld_?.bin -> last checkpoint file (symbolic link)
        !
        do k = 1,9
          call gen_alias(myid,trim(datadir),trim(filename)//'_'//trim(fexts(k))//'.bin','fld_'//trim(fexts(k))//'.bin')
        end do
#if defined(_SCALAR)
        call gen_alias(myid,trim(datadir),trim(filename)//'_'//trim(fexts(k))//'.bin','fld_'//trim(fexts(k))//'.bin') ! k = 6 now
#endif
      end if
      if(myid == 0) print*, '*** Checkpoints saved at time = ', time, 'time step = ', istep, '. ***'
    end if
#if defined(_TIMING)
      !$acc wait(1)
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) print*, 'Avrg, min & max elapsed time: '
      if(myid == 0) print*, dt12av/(1.*product(dims)),dt12min,dt12max
#endif
  end do
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  if(myid == 0) close(csv_unit)
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
end program cans
