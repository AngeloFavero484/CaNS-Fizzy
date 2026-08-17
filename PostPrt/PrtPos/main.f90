program prt_post
  !
  use mod_common
  use mod_io
  !
  implicit none
  !
  call PRT_IO_LoadOpt
  !
  allocate(prt(26*Np))
  allocate(prt_x(Np))
  allocate(prt_y(Np))
  allocate(prt_z(Np))
  allocate(prt_p(3*Np))
  allocate(rad(Np))
  !
  do it=it_min,it_max,it_out
    !
    print'(A)','---------------------------------------------'
    !
    write(*,'(A,i9)') 'It = ',it
    !
    call PRT_IO_LoadPrt(it)
    !
    call PRT_IO_SaveVTK(it)
    !
  enddo
  !
!  call PRT_IO_PrintPVD
  !
  print'(A)','---------------------------------------------'
  print'(A)','Normal end of execution.'
  !
end program prt_post
