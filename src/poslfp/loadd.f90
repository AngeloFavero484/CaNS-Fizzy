module mod_loadd
  !
  use mod_param
  !
  implicit none
  !
  contains
    !
    subroutine loadd(rw,nr,phi,theta)
      !
      implicit none
      !
      integer :: l
      integer :: rw,reclengte,nr
      real(rp) :: phi(:),theta(:)
      real(rp) :: xdum
      !
      inquire (iolength=reclengte) xdum
      !
      if (rw == 0) then
        open(15,file=datadir//'stopiter',access='direct',recl=(2*NL+1)*reclengte)
        read(15,rec=1) (theta(l),l=1,NL),(phi(l),l=1,NL),nr
        close(15)
      endif
      !
      if (rw == 1) then
        open(25,file=datadir//'stopiter',access='direct',recl=(2*NL+1)*reclengte)
        write(25,rec=1) (theta(l),l=1,NL),(phi(l),l=1,NL),nr
        close(25)
      endif
      !
      return
    end subroutine loadd
end module mod_loadd
