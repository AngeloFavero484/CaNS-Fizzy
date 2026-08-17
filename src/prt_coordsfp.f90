module prt_mod_coordsfp
#if defined(_PARTICLE)
  use mod_types
  use mod_common_mpi  , only: myid,boundfrontmyid,boundleftmyid
  use mod_param       , only: pi,l,dl,dli
  use prt_mod_common  , only: ep,nl,dVlagr,dVeul,npmax,pmax, &
                              thetarc,phirc,radfp,nla
  use prt_mod_param   , only: volp,mominert,ratiorho
  implicit none
  private
  public coordsfp !,coordsfp_interior
contains
  !
  subroutine coordsfp(n)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer :: p
#if !defined(_EULER)
    integer :: lp,ll
    real(rp) :: dummy
    real(rp) :: xfploc,yfploc,zfploc
    real(rp) :: isperiodx,isperiody
    real(rp) :: coorxfp,cooryfp   !,coorzfp
    logical :: isout
    !
    ! position of Lfps wrt the center of sphere
    !
    open(25,file='../src/poslfp/data/lagrangianforcepoints_shell')
    read(25,*)
    read(25,*)
    read(25,*)
    do lp=1,NL
      read(25,'(6E16.8)') dummy,dummy,dummy,thetarc(lp),phirc(lp),radfp
    enddo
    close(25)
    !
    ! position of Lfps wrt to the center of computational domain
    !
    !$omp workshare
    nla(:) = 0 ! set to zero from 1 to npmax
    !$omp end workshare
    !
    !$omp parallel default(none) &
    !$omp shared(ep,nla,pmax,nl,l,n,dli,dl,npmax) &
    !$omp shared(radfp,phirc,thetarc) &
    !$omp shared(boundleftmyid,boundfrontmyid) &
    !$omp private(p,lp,ll,coorxfp,cooryfp,xfploc,yfploc,zfploc) &
    !$omp private(isperiodx,isperiody,isout)
    !$omp do 
    do p=1,pmax
       if (ep(p)%mslv /= 0) then
          ! myid is master of particle ep(p)%mslv
          ll = 0
          do lp=1,NL
             xfploc = ep(p)%x + radfp*sin(ep(p)%theta*0._rp+thetarc(lp))*cos(ep(p)%phi*0._rp+phirc(lp))
             yfploc = ep(p)%y + radfp*sin(ep(p)%theta*0._rp+thetarc(lp))*sin(ep(p)%phi*0._rp+phirc(lp))
             zfploc = ep(p)%z + radfp*cos(ep(p)%theta*0._rp+thetarc(lp))
             !
             isperiodx = 0._rp
             isperiody = 0._rp
             !
             if (xfploc < 0._rp+0.5_rp*dl(1)) isperiodx =  1._rp
             if (xfploc > l(1)+0.5_rp*dl(1))  isperiodx = -1._rp
             if (yfploc < 0._rp+0.5_rp*dl(2)) isperiody =  1._rp
             if (yfploc > l(2)+0.5_rp*dl(2))  isperiody = -1._rp
             !
             isout = .false.
             coorxfp = (xfploc+isperiodx*l(1)-boundleftmyid )*dli(1)
             if( nint(coorxfp) < 1 .or. nint(coorxfp) > n(1) ) isout = .true.
             cooryfp = (yfploc+isperiody*l(2)-boundfrontmyid)*dli(2)
             if( nint(cooryfp) < 1 .or. nint(cooryfp) > n(2) ) isout = .true.
             !
             if(.not.isout) then
               ll = ll + 1
               ep(p)%xfp(ll) = xfploc
               ep(p)%yfp(ll) = yfploc
               ep(p)%zfp(ll) = zfploc
               ep(p)%ul(ll) = ep(p)%u + ep(p)%omy*(ep(p)%zfp(ll)-ep(p)%z) &
                                      - ep(p)%omz*(ep(p)%yfp(ll)-ep(p)%y)
               ep(p)%vl(ll) = ep(p)%v + ep(p)%omz*(ep(p)%xfp(ll)-ep(p)%x) &
                                      - ep(p)%omx*(ep(p)%zfp(ll)-ep(p)%z)
               ep(p)%wl(ll) = ep(p)%w + ep(p)%omx*(ep(p)%yfp(ll)-ep(p)%y) &
                                      - ep(p)%omy*(ep(p)%xfp(ll)-ep(p)%x)
             endif
          enddo
          nla(p) = ll
       endif
    enddo
    !$omp end parallel
    !
    ! dVlagr of outer shell
    !
    dVlagr = ( (4._rp/3._rp)*pi*(radfp+0.5_rp/dli(1))**3 - &
               (4._rp/3._rp)*pi*(radfp-0.5_rp/dli(1))**3 )/(1._rp*nl)
    !
    ! Volume of an Eulerian grid cell
    !
#endif
    dVeul  = 1._rp/(dli(1)*dli(2)*dli(3))
    !
#if !defined(_EULER)
    if (myid .eq. 0) then
       write(6,*) 'dVlagr, dVeul = ',dVlagr, dVeul
    endif
#endif
    !
    ! Define geometrical properties of the particles
    !
    !$omp parallel default(none)  &
    !$omp shared(ep,nla,npmax)     &
    !$omp private(p)
    !$omp do 
    do p=1,npmax
      ep(p)%vol = volp
      ep(p)%ratiorho = ratiorho
      ep(p)%mominert = mominert
    enddo
    !$omp end parallel
    !
    ! data to file
    !
    !if (ap(1)%mslv .eq. 1) then
    !  open(25,file=datadir//'lagrangianforcepoints')
    !  write(25,*) 'VARIABLES = "lfpx","lfpy","lfpz","radius"'
    !  write(25,*) 'ZONE T="Zone1"',' I=',NL,', F=POINT'
    !  write(25,*) ''
    !  do l=1,nl
    !    write(25,'(4E16.8)') ap(1)%xfp(l),ap(1)%yfp(l),ap(1)%zfp(l), &
    !                         sqrt( (ap(1)%xfp(l)-ap(1)%x)**2. + (ap(1)%yfp(l)-ap(1)%y)**2. + &
    !                               (ap(1)%zfp(l)-ap(1)%z)**2. )
    !  enddo
    !  close(25)
    !endif
    !
    return
  end subroutine coordsfp
  !
!  subroutine coordsfp_interior(api)
!    implicit none
!    integer :: l,p,q
!    type(particle_interior), dimension(npmax), intent(out) :: api
!    real(rp) ::  dummy,thetarcint(1:NLtot),phircint(1:NLtot),radfp2,radfp3,radfp4
!    !
!    ! position of Lagrangian force points wrt center of sphere
!    !
!    open(42,file='poslfp/4ringstogether/data/lagrangianforcepoints2')
!    read(42,*)
!    read(42,*)
!    read(42,*)
!    do l=1,NL
!       read(42,'(6E16.8)') dummy,dummy,dummy,thetarc(l),phirc(l),radfp
!       thetarcint(l) = thetarc(l)
!       phircint(l)   = phirc(l)
!    enddo
!    do q=1,NL2
!       l=NL+q
!       read(42,'(6E16.8)') dummy,dummy,dummy,thetarcint(l),phircint(l),radfp2
!    enddo
!    do q=1,NL3
!       l=NL+NL2+q
!       read(42,'(6E16.8)') dummy,dummy,dummy,thetarcint(l),phircint(l),radfp3
!    enddo
!    do q=1,NL4
!       l=NL+NL2+NL3+q
!       read(42,'(6E16.8)') dummy,dummy,dummy,thetarcint(l),phircint(l),radfp4
!    enddo
!    close(42)
!    !
!    ! position of lfp's wrt to center of computational domain
!    !
!    do p=1,pmax
!       if (ap(p)%mslv .gt. 0) then
!          ! myid is master of particle ap(p)%mslv
!          api(p)%x = ap(p)%x
!          api(p)%y = ap(p)%y
!          api(p)%z = ap(p)%z
!          do l=1,NL
!             api(p)%xfp(l)      = ap(p)%x + radfp*sin(ap(p)%theta+thetarc(l))*cos(ap(p)%phi+phirc(l))
!             api(p)%yfp(l)      = ap(p)%y + radfp*sin(ap(p)%theta+thetarc(l))*sin(ap(p)%phi+phirc(l))
!             api(p)%zfp(l)      = ap(p)%z + radfp*cos(ap(p)%theta+thetarc(l))
!          enddo
!          do q=1,NL2
!             l=NL+q
!             api(p)%xfp(l) = ap(p)%x + radfp2*sin(ap(p)%theta+thetarcint(l))*cos(ap(p)%phi+phircint(l))
!             api(p)%yfp(l) = ap(p)%y + radfp2*sin(ap(p)%theta+thetarcint(l))*sin(ap(p)%phi+phircint(l))
!             api(p)%zfp(l) = ap(p)%z + radfp2*cos(ap(p)%theta+thetarcint(l))
!          enddo
!          do q=1,NL3
!             l=NL+NL2+q
!             api(p)%xfp(l) = ap(p)%x + radfp3*sin(ap(p)%theta+thetarcint(l))*cos(ap(p)%phi+phircint(l))
!             api(p)%yfp(l) = ap(p)%y + radfp3*sin(ap(p)%theta+thetarcint(l))*sin(ap(p)%phi+phircint(l))
!             api(p)%zfp(l) = ap(p)%z + radfp3*cos(ap(p)%theta+thetarcint(l))
!          enddo
!          do q=1,NL4
!             l=NL+NL2+NL3+q
!             api(p)%xfp(l) = ap(p)%x + radfp4*sin(ap(p)%theta+thetarcint(l))*cos(ap(p)%phi+phircint(l))
!             api(p)%yfp(l) = ap(p)%y + radfp4*sin(ap(p)%theta+thetarcint(l))*sin(ap(p)%phi+phircint(l))
!             api(p)%zfp(l) = ap(p)%z + radfp4*cos(ap(p)%theta+thetarcint(l))
!          enddo
!       else
!          ! myid is either slave of particle ap(p)%mslv or does not contain this particle
!          do l=1,NL+nl2+nl3+nl4
!             api(p)%xfp(l) = 0.
!             api(p)%yfp(l) = 0.
!             api(p)%zfp(l) = 0.
!          enddo
!       endif
!    enddo
!    !
!    ! Volume of a Lagr. force point cell (lfp-cell)
!    !
!    do p=1,pmax
!       do l=1,NL
!          api(p)%dvlagr(l)= ( (4./3.)*pi*(radfp+0.5/dxi)**3. - &
!               (4./3.)*pi*(radfp-0.5/dxi)**3. )/(1.*NL)
!       enddo
!       do q=1,NL2
!          l=NL+q
!          api(p)%dvlagr(l) = ( (4./3.)*pi*(radfp-0.5/dxi)**3. - &
!               (4./3.)*pi*(radfp-1.5/dxi)**3. )/(1.*NL2)
!       enddo
!       do q=1,NL3
!          l=NL+NL2+q
!          api(p)%dvlagr(l) = ( (4./3.)*pi*(radfp-1.5/dxi)**3. - &
!               (4./3.)*pi*(radfp-2.5/dxi)**3. )/(1.*NL3)
!       enddo
!       do q=1,NL4
!          l=NL+NL2+NL3+q
!          api(p)%dvlagr(l) = ( (4./3.)*pi*(radfp-2.5/dxi)**3. - &
!               (4./3.)*pi*(radfp-3.5/dxi)**3. )/(1.*NL4)
!       enddo
!    enddo
!    !
!    return
!  end subroutine coordsfp_interior
  !
#endif
end module prt_mod_coordsfp
