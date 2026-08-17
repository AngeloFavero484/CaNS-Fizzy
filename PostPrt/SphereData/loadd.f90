module mod_loadd
  !
  use mod_param      , only: nl,ng,iter,rp,cx,cy,cz,pi
  use mod_common     , only: lpx,lpy,lpz,nx,ny,nz,lpA
  !
  implicit none
  private
  public flpcoord,loadfld
  !
  contains
  !
  subroutine flpcoord(thetarc,phirc,rad)
    !
    implicit none
    !
    real(rp),dimension(nl),intent(out) :: thetarc,phirc
    real(rp),intent(in) :: rad
    integer :: ll
    real(rp) :: dummy
    real(rp) :: normx,normy,normz
    !
    open(25,file='../../src/poslfp/data/lagrangianforcepoints_shell')
    read(25,*)
    read(25,*)
    read(25,*)
    do ll=1,nl
      read(25,'(6E16.8)') dummy,dummy,dummy,thetarc(ll),phirc(ll),dummy
    enddo
    close(25)
    !
    do ll=1,nl
      lpx(ll) = cx + rad*sin(thetarc(ll))*cos(phirc(ll))
      lpy(ll) = cy + rad*sin(thetarc(ll))*sin(phirc(ll))
      lpz(ll) = cz + rad*cos(thetarc(ll))
      !
      normx=lpx(ll)-cx
      normy=lpy(ll)-cy
      normz=lpz(ll)-cz
      !
      nx(ll)=normx/sqrt(normx**2+normy**2+normz**2)
      ny(ll)=normy/sqrt(normx**2+normy**2+normz**2)
      nz(ll)=normz/sqrt(normx**2+normy**2+normz**2)
    enddo
    !
    lpA=4.0_rp*pi*rad**2/nl
    !
    return
  end subroutine flpcoord
  !
  subroutine loadfld(vel_x,vel_y,vel_z,pre)
    !
    implicit none
    !
    real(rp),dimension(1:ng(1),1:ng(2),1:ng(3)),intent(inout) :: vel_x,vel_y,vel_z,pre
    !
    integer :: iunit 
    character(len=15) :: path
    character(len=7) :: iterchar
    !
    path='../../run/data/'
    write(iterchar,'(i7.7)') iter
    !
    ! read u
    open(newunit=iunit,file=trim(path)//'vex_fld_'//trim(iterchar)//'.bin', &
         action='read',access='stream',form='unformatted',status='old')
    read(iunit) vel_x(1:ng(1),1:ng(2),1:ng(3))
    close(iunit)
    ! read v
    open(newunit=iunit,file=trim(path)//'vey_fld_'//trim(iterchar)//'.bin', &
         action='read',access='stream',form='unformatted',status='old')
    read(iunit) vel_y(1:ng(1),1:ng(2),1:ng(3))
    close(iunit)
    ! read w
    open(newunit=iunit,file=trim(path)//'vez_fld_'//trim(iterchar)//'.bin', &
         action='read',access='stream',form='unformatted',status='old')
    read(iunit) vel_z(1:ng(1),1:ng(2),1:ng(3))
    close(iunit)
    ! read p
    open(newunit=iunit,file=trim(path)//'pre_fld_'//trim(iterchar)//'.bin', &
         action='read',access='stream',form='unformatted',status='old')
    read(iunit) pre(1:ng(1),1:ng(2),1:ng(3))
    close(iunit)
    !
    return
  end subroutine loadfld
  !
!  subroutine loadd(rw,nr,phi,theta)
!    !
!    implicit none
!    !
!    integer :: l
!    integer :: rw,reclengte,nr
!    real(rp) :: phi(:),theta(:)
!    real(rp) :: xdum
!    !
!    inquire (iolength=reclengte) xdum
!    !
!    if (rw == 0) then
!      open(15,file=datadir//'stopiter',access='direct',recl=(2*NL+1)*reclengte)
!      read(15,rec=1) (theta(l),l=1,NL),(phi(l),l=1,NL),nr
!      close(15)
!    endif
!    !
!    if (rw == 1) then
!      open(25,file=datadir//'stopiter',access='direct',recl=(2*NL+1)*reclengte)
!      write(25,rec=1) (theta(l),l=1,NL),(phi(l),l=1,NL),nr
!      close(25)
!    endif
!    !
!    return
!  end subroutine loadd
end module mod_loadd
