! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_massbal
  !
  ! volume budget of the fluid-1 phase, split by the solid indicator alphac.
  !
  ! the contact-line relaxation in mod_extend advects psi in non-conservative
  ! (advective) form over the band alpha_min < alphac < 1, so unlike the THINC
  ! flux update in mod_two_fluid it is a genuine source/sink of psi. Part of
  ! what it injects lands inside the particle, where it is not physical fluid.
  ! Splitting the integral separates the two:
  !
  !   vol(1) = sum psi             *dV   total fluid-1 volume in the domain
  !   vol(2) = sum psi*(1-alphac)  *dV   fluid-1 volume outside the solid
  !   vol(3) = sum psi*alphac      *dV   fluid-1 volume buried in the solid
  !
  ! vol(1) = vol(2) + vol(3) identically, so vol(2) is the quantity that should
  ! stay constant if the extension only writes psi into the solid interior.
  !
  ! vol(4) and vol(5) are the out-of-range content of psi,
  !
  !   vol(4) = sum max(psi-1,0)*dV   vol(5) = sum max(-psi,0)*dV
  !
  ! which mod_extend can produce because it advects psi with no bound: nothing
  ! clips between the relaxation loop and the end of the next rk_2fl, where
  ! clip_field (rk.f90:252) truncates it. That truncation is a volume sink the
  ! advection gets blamed for, so it has to be measured separately.
  !
  ! called twice per step from main.f90 -- once after the interface advection
  ! and once after the relaxation loop -- so that the two contributions to the
  ! drift can be attributed separately.
  !
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr
#if defined(_PARTICLE)
  use prt_mod_common, only: alphac
#endif
  implicit none
  private
  public cmpt_massbal
  contains
  subroutine cmpt_massbal(n,dl,dzf,psi,vol)
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dl
    real(rp), intent(in ), dimension(0:)       :: dzf
    real(rp), intent(in ), dimension(0:,0:,0:) :: psi
    real(rp), intent(out), dimension(5)        :: vol
    real(rp), dimension(5) :: vol_l
    real(rp) :: vtot,vout,vin,vover,vundr
    real(rp) :: dv,pf,af
    integer  :: i,j,k
    !
    vtot  = 0._rp
    vout  = 0._rp
    vin   = 0._rp
    vover = 0._rp
    vundr = 0._rp
    !$acc parallel loop collapse(3) default(present) private(dv,pf,af) &
    !$acc reduction(+:vtot,vout,vin,vover,vundr) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dv = dl(1)*dl(2)*dzf(k)
          pf = psi(i,j,k)
#if defined(_PARTICLE)
          af = alphac(i,j,k)
#else
          af = 0._rp
#endif
          vtot = vtot + pf*dv
          vout = vout + pf*(1._rp-af)*dv
          vin  = vin  + pf*af*dv
          vover = vover + max(pf-1._rp,0._rp)*dv
          vundr = vundr + max(-pf,0._rp)*dv
        end do
      end do
    end do
    !$acc wait(1)
    vol_l(1:5) = [vtot,vout,vin,vover,vundr]
    call MPI_ALLREDUCE(vol_l,vol,5,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_massbal
end module mod_massbal
