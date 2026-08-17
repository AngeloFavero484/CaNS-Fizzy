module prt_mod_output
#if defined(_PARTICLE)
  use mpi
  use mod_types
  use mod_param        , only: datadir
  use mod_common_mpi   , only: myid,ierr,prt_comm_cart
  use mod_param        , only: dims
  use prt_mod_param    , only: np,radius
  use prt_mod_common   , only: ep,rkp,npmstr,pmax, &
                               fx_tot,fy_tot,fz_tot
  !
  implicit none
  private
  public outpart
  !
  contains
  !
  subroutine outpart(nr)
    implicit none
    integer :: skip,skipacc
    integer,dimension(0:dims(1)*dims(2)-1) :: npmstr_glob,npmstr_glob_all
    real(rp), allocatable, dimension(:,:) :: posp,velp,ridp,angpos,angvel,accel,angaccel,extravec,extrascal,force
    real(rp) :: aux_x,aux_y,aux_z
    integer i,p,idp
    integer :: fh
    integer, intent(in) :: nr
    real(rp) :: xdum
    integer :: lenr
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    integer :: mydisp
    character(len=7) :: istepchar
    !
    !inquire (iolength=lenr) xdum
    !lenr = sizeof(xdum)
    lenr=storage_size(xdum)/8
    write(istepchar,'(i7.7)') nr
    !
    ! write particle related data directly in parallel with MPI-IO
    !
    npmstr_glob(:) = 0
    npmstr_glob(myid) = npmstr
    call MPI_ALLREDUCE(npmstr_glob(0),npmstr_glob_all(0),product(dims),MPI_INTEGER,MPI_SUM,prt_comm_cart,ierr)
    mydisp = 0
    if(myid /= 0) mydisp = sum(npmstr_glob_all(0:myid-1))
    allocate(posp(3,npmstr))
    allocate(velp(3,npmstr))
    allocate(ridp(1,npmstr))
    allocate(angpos(3,npmstr))
    allocate(angvel(3,npmstr))
    allocate(accel(3,npmstr))
    allocate(angaccel(3,npmstr))
    allocate(force(3,npmstr))
    allocate(extravec(3,npmstr))
    allocate(extrascal(1,npmstr))
    i = 0
    do p=1,pmax
       if(ep(p)%mslv > 0) then
          idp = ep(p)%mslv
          i = i + 1
          posp(1,i)      = ep(p)%x
          posp(2,i)      = ep(p)%y
          posp(3,i)      = ep(p)%z
          ridp(1,i)      = 1.0_rp*ep(p)%mslv
          velp(1,i)      = ep(p)%u
          velp(2,i)      = ep(p)%v
          velp(3,i)      = ep(p)%w
          aux_x          = sin(ep(p)%theta)*cos(ep(p)%phi)
          aux_y          = sin(ep(p)%theta)*sin(ep(p)%phi)
          aux_z          = cos(ep(p)%phi)
          angpos(1,i)    = acos(aux_y/sqrt(aux_y**2.0_rp+aux_z**2.0_rp))
          angpos(2,i)    = acos(aux_z/sqrt(aux_x**2.0_rp+aux_z**2.0_rp))
          angpos(3,i)    = ep(p)%phi !acos(aux_x/sqrt(aux_x**2.0_rp+aux_y**2.0_rp)) or ep(p)%phi
          angvel(1,i)    = ep(p)%omx
          angvel(2,i)    = ep(p)%omy
          angvel(3,i)    = ep(p)%omz
          accel(1,i)     = rkp(p)%dudt
          accel(2,i)     = rkp(p)%dvdt
          accel(3,i)     = rkp(p)%dwdt
          angaccel(1,i)  = rkp(p)%domxdt
          angaccel(2,i)  = rkp(p)%domydt
          angaccel(3,i)  = rkp(p)%domzdt
          force(1,i)     = fx_tot(p) !ep(p)%fxltot
          force(2,i)     = fy_tot(p) !ep(p)%fyltot
          force(3,i)     = fz_tot(p) !ep(p)%fzltot
          extravec(1,i)  = radius*aux_x
          extravec(2,i)  = radius*aux_y
          extravec(3,i)  = radius*aux_z
          extrascal(1,i) = 0.0_rp
       endif
    enddo
    !
    skipacc = 0
    call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(datadir)//'partvis_'//istepchar//'.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    !
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,posp(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif 
    !
    skipacc = skipacc + skip
    skip = 1 ! scalar
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,ridp(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    !
    skipacc = skipacc + skip
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,velp(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    !
    skipacc = skipacc + skip
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,angpos(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    !
    skipacc = skipacc + skip
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,angvel(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    !
    skipacc = skipacc + skip
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,accel(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    !
    skipacc = skipacc + skip
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,angaccel(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    !
    skipacc = skipacc + skip
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,force(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    !
    skipacc = skipacc + skip
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL, ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,extravec(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    !
    skipacc = skipacc + skip
    skip = 1 ! scalar
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL_RP,MPI_REAL_RP, 'native', &
         MPI_INFO_NULL,ierr)
    if (npmstr > 0) then
      call MPI_FILE_WRITE(fh,extrascal(1,1),skip*npmstr,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
    endif
    call MPI_FILE_CLOSE(fh,ierr)
    !
    deallocate(posp)
    deallocate(velp)
    deallocate(ridp)
    deallocate(angpos)
    deallocate(angvel)
    deallocate(accel)
    deallocate(angaccel)
    deallocate(extravec)
    deallocate(extrascal)
    !
    return
  end subroutine outpart
#endif
end module prt_mod_output
