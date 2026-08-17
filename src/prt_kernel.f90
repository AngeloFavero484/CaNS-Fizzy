module prt_mod_kernel
#if defined(_PARTICLE)
  use mod_types
  !
  implicit none
  private
  public kernel,kernelchk,kerneltest
  !
  contains
  !
  real(rp) function kernel(r)
    !
    ! Discretized delta function. (N.B.: multiplied by the grid spacing!)
    !
    implicit none
    real(rp) :: r
    !
    if (abs(r) > 1.5_rp) then
      kernel=0.0_rp
    else
      if (abs(r) > 0.5_rp) then
        kernel = (1.0_rp/6.0_rp)*(5.0_rp-3.0_rp*abs(r)-sqrt(-3.0_rp*((1.0_rp-abs(r))**2)+1.0_rp))
      else
        kernel = (1.0_rp/3.0_rp)*(1.0_rp+sqrt(-3.0_rp*(r**2)+1.0_rp))
      endif
    endif
    !
    return
  end function kernel
  !
  integer function kernelchk(r)
    !
    ! checks if r is within the kernel stencil. if so kernelchk(r) = 1
    !
    implicit none
    real(rp) :: r
    !
    if (abs(r) > 1.5) then
      kernelchk = 0
    else
      kernelchk = 1
    endif
    !
    return
  end function kernelchk
  !
  subroutine kerneltest(sum)
    integer :: i,j,k
    real(rp) :: coorx,coory,coorz
    real(rp), intent(out) :: sum
    !
    ! tests if function kernel is working properly (sum ~ 1)
    !
    sum = 0.0_rp
    do k=-25,25
      coorz=k/12.0_rp
      do j=-25,25
        coory=j/12.0_rp
        do i=-25,25
          coorx=i/12.0_rp
          sum = sum + (1.0_rp/12._rp)*(1.0_rp/12.0_rp)*(1.0_rp/12.0_rp)*kernel(coorx)*kernel(coory)*kernel(coorz)
        enddo
      enddo
    enddo
    !
    return
  end subroutine kerneltest
  !
#endif
end module prt_mod_kernel
