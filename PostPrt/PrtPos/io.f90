module mod_io
  ! 
  use mod_common
  !
  ! ===================================================================================================================
  !
  implicit none
  !
  integer :: fid,err
  !
  logical :: open_check
  character(256) :: open_path
  character(256) :: open_file
  character(32) :: fileres
  !
  character(32),parameter :: fileres_01='partvis_1234567.bin'
  !
  ! ===================================================================================================================
  ! private/public access specifications
  !
  private
  !
  public :: PRT_IO_LoadOpt
  public :: PRT_IO_LoadPrt
  public :: PRT_IO_SaveVTK
  public :: PRT_IO_PrintPVD
  !
  ! ===================================================================================================================
  !
  contains
  !
  subroutine PRT_IO_LoadOpt
    implicit none
    !
    print'(A)','Reading ./param.dat ...'
    !
    open(newunit=fid,file='./param.dat',form='formatted',action='read',iostat=err)
    !
    if(err==0) then
      !
      read(fid,*) open_path
      read(fid,*) Np
      read(fid,*) radius
      read(fid,*) it_min
      read(fid,*) it_max
      read(fid,*) it_out
      read(fid,*) t0
      read(fid,*) dt
      close(fid)
      !
      print'(A)','Reading completed.'
      !
    else
      !
      print'(A)','ERROR: File ./opt.dat not found.'
      print'(A)','ABORTED!'
      stop
      !
    endif
    !
    return   
  end subroutine PRT_IO_LoadOpt
  !
  ! ===================================================================================================================
  !
  subroutine PRT_IO_LoadPrt(iostep)
    implicit none
    integer,intent(in) :: iostep
    !
    integer :: i
    !
    fileres=fileres_01
    call PRT_IO_PrintStep(iostep)
    call PRT_IO_PrintPath
    !
    print'(A)','Reading '//trim(open_file)
    !
    inquire(file=open_file,exist=open_check) 
    !
    if(open_check) then
      !
      open(newunit=fid,file=open_file,form='unformatted',access='stream',status='old',action='read')
      !
      read(fid) prt
      !
      do i=1,Np
        prt_x(i) = prt(3*i - 2)
        prt_y(i) = prt(3*i - 1)
        prt_z(i) = prt(3*i - 0)
        rad(i) = radius
      enddo
      !
      prt_p(1:3*Np) = prt(1:3*Np)
      !
      close(fid)
    else
      print'(A)','ERROR: File '//trim(open_file)//' not found.'
      print'(A)','ABORTED!'
      stop
    endif
    !
    print'(A)', 'Reading completed.'
    !
    return
  end subroutine PRT_IO_LoadPrt
  !
  ! ===================================================================================================================
  !
  subroutine PRT_IO_SaveVTK(iostep)
    implicit none
    integer,intent(in) :: iostep
    !
    character(len=1),parameter :: newline = char(10)
    character(len=100) :: buffer
    integer :: i
    integer :: iostatus
    !
    fileres='prt_1234567.vtk'
    write(fileres(5:11),'(I7.7)') iostep
    call PRT_IO_PrintPath
    !
    buffer=' '
    !
    print'(A)','Writing '//trim(open_file)
    !
    open(newunit=fid,file=open_file,status='replace',form='unformatted',access='stream',&
         action='write',convert='BIG_ENDIAN',iostat=iostatus)
    !
    write(fid) '# vtk DataFile Version 3.0'//char(10)
    write(fid) 'Binary particle data'//char(10)
    write(fid) 'BINARY'//char(10)
    write(fid) char(10)
    !
    write(fid) 'DATASET POLYDATA'//char(10)
    write(buffer,fmt='(A,1I12,A)') 'POINTS ',Np,' double'
    write(fid) trim(buffer)//char(10)
    !
    do i=1,Np
      write(fid,iostat=iostatus) prt_x(i),prt_y(i),prt_z(i)
    enddo
    !
    write(fid) char(10)
    write(buffer,fmt='(A11,I12)') 'POINT_DATA ',Np
    write(fid) trim(buffer)//char(10)
    write(fid) 'SCALARS radius double'//char(10)
    write(fid) 'LOOKUP_TABLE default'//char(10)
    !
    do i = 1, Np
      write(fid,iostat=iostatus) rad(i)
    end do
    !
    close(fid)
    !
    print'(A)', 'Writing completed.'
    !
    return
  end subroutine PRT_IO_SaveVTK
  !
  ! ===================================================================================================================
  !
  subroutine PRT_IO_PrintPVD
    implicit none
    character(len=1),parameter :: newline = char(10)
    !
    character(len=100) :: buffer
    character(len=7) :: itchar
    character(len=15) :: filevtk
    integer :: i
    integer :: indent
    !
    fileres='prt.pvd'
    call PRT_IO_PrintPath
    !
    print'(A)','---------------------------------------------'
    print'(A)','Writing '//trim(open_file)
    !
    indent=0
    !
    open(newunit=fid,file=open_file,status='replace',form='formatted',action='write')
    !
    write(fid,fmt='(A)') repeat(' ',indent)//'<?xml version="1.0"?>'
    write(fid,fmt='(A)') repeat(' ',indent)//'<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">'
    indent = indent + 4
      write(fid,fmt='(A)') repeat(' ',indent)//'<Collection>'
      !
      print*,it_min,it_max,it_out
      !
      do it=it_min,it_max,it_out
        write(itchar,fmt='(i7.7)') it
        indent = indent + 4
          filevtk='prt_1234567.vtk'
          write(filevtk(5:11),'(I7.7)') it
          write(fid,fmt='(A,F6.3,3A)') repeat(' ',indent)//'<DataSet timestep="',t0 + 1.0_rp*(it-it_out)*dt, &
                                       '" group="" part="0" file="',filevtk,'"/>'
          indent = indent - 4
      enddo
      !
      write(fid,fmt='(A)') repeat(' ',indent)//'</Collection>'
      indent = indent - 4
    write(fid,fmt='(A)') repeat(' ',indent)//'</VTKFile>'
    !
    close(fid)
    !
    print'(A)', 'Writing completed.'
    !
  end subroutine PRT_IO_PrintPVD
  !
  ! ===================================================================================================================
  !
  subroutine PRT_IO_PrintStep(iostep)
    implicit none
    integer,intent(in) :: iostep
    !
    write(fileres(9:15),1) iostep
    !
  1 format(1I7.7)
    !
    return
  end subroutine PRT_IO_PrintStep
  !
  ! ===================================================================================================================
  !
  subroutine PRT_IO_PrintPath
    implicit none
    !
    open_file=trim(open_path)//trim(fileres) 
    !
    return
  end subroutine PRT_IO_PrintPath
  !
  ! ===================================================================================================================
  !
end module mod_io
