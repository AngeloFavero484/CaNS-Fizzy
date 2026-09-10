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
  public cmpt_massbal, crrct_vout
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
  !
  subroutine crrct_vout(n,dl,dzf,dvol,psi,dvol_res)
    !
    ! removes dvol of fluid-1 volume from *outside* the solid, so that
    ! vol(2) = sum psi*(1-alphac)*dV is restored to what it was before the
    ! contact-line relaxation ran. dvol is that relaxation's injection,
    ! vol_ext(2)-vol_adv(2) as measured by cmpt_massbal above.
    !
    ! the relaxation has no discrete conservation property (see the module
    ! header), so the volume it injects cannot be traced back to donor cells.
    ! It is taken out instead as a constrained projection: find the smallest
    ! correction, in the weighted sense below, satisfying
    !
    !   sum dpsi*(1-alphac)*dV = -dvol
    !
    ! with dpsi = -c*g and the weight
    !
    !   g = psi*(1-psi)*(1-alphac)
    !
    ! so that c = dvol/sum g*(1-alphac)*dV. g is the Lagrange-multiplier
    ! weight and it is what makes the correction land where the error came
    ! from rather than smeared over the drop:
    !
    !   psi*(1-psi) vanishes in both bulk phases -> only interface cells pay
    !   (1-alphac)  vanishes in the solid interior -> nothing is taken from
    !               volume that is not physical fluid anyway, and the same
    !               factor appears in the constraint, so a cell is debited in
    !               proportion to how much it actually contributes to vol(2)
    !
    ! the sweep is restricted to the diffuse solid shell 0 < alphac < 1: that
    ! is the relaxation band plus the donor cells its upwind stencil reads,
    ! i.e. the support of the error. The drop's free surface away from the
    ! particle is deliberately left alone.
    !
    ! psi*(1-psi) also makes the correction bound-preserving on its own. With
    ! psi in [0,1] and alphac in (0,1), driving a cell below 0 needs
    ! c*(1-psi)*(1-alphac) > 1 and above 1 needs |c|*psi*(1-alphac) > 1, so any
    ! |c| < 1 is safe for every cell at once. c_max sits just under that.
    !
    ! one pass is exact whenever |c| <= c_max, which is the normal case. When
    ! the injection is too large for a single capped pass -- the seeding
    ! relaxation before the time loop, mainly -- the pass is repeated: psi has
    ! moved, so wsum and the remaining debt are both recomputed and the debt
    ! falls geometrically. Whatever n_iter passes still leave over comes back
    ! in dvol_res for the caller to defer rather than silently drop.
    !
    implicit none
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(in   ), dimension(3)        :: dl
    real(rp), intent(in   ), dimension(0:)       :: dzf
    real(rp), intent(in   )                      :: dvol
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    real(rp), intent(out  )                      :: dvol_res
    real(rp), parameter :: c_max = 0.9_rp
    integer , parameter :: n_iter = 20
    real(rp) :: wsum_l,wsum,c
    real(rp) :: dv,pf,af,g
    integer  :: i,j,k,it
    !
    dvol_res = dvol
    !
    do it = 1,n_iter
      if(dvol_res == 0._rp) exit
      wsum_l = 0._rp
      !$acc parallel loop collapse(3) default(present) private(dv,pf,af,g) &
      !$acc reduction(+:wsum_l) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
#if defined(_PARTICLE)
            af = alphac(i,j,k)
#else
            af = 0._rp
#endif
            if(af > 0._rp .and. af < 1._rp) then
              dv = dl(1)*dl(2)*dzf(k)
              pf = psi(i,j,k)
              g  = pf*(1._rp-pf)*(1._rp-af)
              wsum_l = wsum_l + g*(1._rp-af)*dv
            end if
          end do
        end do
      end do
      !$acc wait(1)
      call MPI_ALLREDUCE(wsum_l,wsum,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      !
      ! no interface anywhere in the shell: nothing to correct against, the whole
      ! injection is left on the table and reported
      !
      if(wsum <= 0._rp) return
      !
      c = dvol_res/wsum
      if(abs(c) > c_max) then
        c = sign(c_max,c)
        dvol_res = dvol_res - c*wsum
      else
        dvol_res = 0._rp
      end if
      !
      !$acc parallel loop collapse(3) default(present) private(pf,af,g) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
#if defined(_PARTICLE)
            af = alphac(i,j,k)
#else
            af = 0._rp
#endif
            if(af > 0._rp .and. af < 1._rp) then
              pf = psi(i,j,k)
              g  = pf*(1._rp-pf)*(1._rp-af)
              psi(i,j,k) = pf - c*g
            end if
          end do
        end do
      end do
      !$acc wait(1)
    end do
  end subroutine crrct_vout
end module mod_massbal
