module mod_extend
  !
  use mpi
  use mod_types
  use mod_param         , only: pi,sigma,alpha_min,alpha_ramp
#if defined(_PARTICLE)
    use prt_mod_common    , only: alphac,norm_partx,norm_party,norm_partz
#endif
  implicit none
  private
  public compute_uextend, advect_vof_upwind
  contains

  subroutine compute_uextend(n, theta, normx, normy, normz, u_ext, v_ext, w_ext)
    implicit none
    ! Input
    integer,  intent(in), dimension(3) :: n
    real(rp), intent(in)  :: theta
    real(rp), intent(in), dimension(0:,0:,0:) :: normx,normy,normz
    ! Output
    real(rp), intent(inout), dimension(0:,0:,0:) :: u_ext
    real(rp), intent(inout), dimension(0:,0:,0:) :: v_ext
    real(rp), intent(inout), dimension(0:,0:,0:) :: w_ext
    
    ! Variabili locali
    real(rp) :: n_int(3), n_wall(3), n1(3), n2(3)
    real(rp) :: norm_n1, norm_n2, norm_uext
    real(rp) :: c, theta_rad, cot_theta, sin_theta
    real(rp), parameter :: eps = epsilon(1._rp)
    !
    ! floor on sin(theta) in the cotangent below. theta -> 0 or 180 makes
    ! cot(pi-theta) singular; u_ext is normalised a few lines later, so once
    ! |cot_theta| is this large the result is already indistinguishable from
    ! the purely tangential +/-n2 and the floor changes nothing physical --
    ! it only keeps the two end points from producing Inf/NaN.
    !
    real(rp), parameter :: sin_min = 1.e-6_rp
    integer  :: i, j, k
    !
    ! theta is a run constant: the cotangent is loop-invariant
    !
    theta_rad = theta * pi / 180.0_rp
    sin_theta = max(sin(pi-theta_rad),sin_min)
    cot_theta = cos(pi-theta_rad) / sin_theta
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if (alphac(i,j,k) > alpha_min .and. alphac(i,j,k) < 1._rp ) then
            n_wall(1)=-norm_partx(i,j,k)
            n_wall(2)=-norm_party(i,j,k)
            n_wall(3)=-norm_partz(i,j,k)
            n_int(1)=normx(i,j,k)
            n_int(2)=normy(i,j,k)
            n_int(3)=normz(i,j,k)
            ! Vettore parallelo alla linea di contatto
            n1(1) = n_int(2)*n_wall(3) - n_int(3)*n_wall(2)
            n1(2) = n_int(3)*n_wall(1) - n_int(1)*n_wall(3)
            n1(3) = n_int(1)*n_wall(2) - n_int(2)*n_wall(1)
            !
            norm_n1 = max(sqrt(n1(1)**2 + n1(2)**2 + n1(3)**2),eps)
            n1(1) = -n1(1) / norm_n1
            n1(2) = -n1(2) / norm_n1
            n1(3) = -n1(3) / norm_n1
            ! Vettore tangente al solido, ortogonale alla linea di contatto
            n2(1) = n1(2)*n_wall(3) - n1(3)*n_wall(2)
            n2(2) = n1(3)*n_wall(1) - n1(1)*n_wall(3)
            n2(3) = n1(1)*n_wall(2) - n1(2)*n_wall(1)
            !
            norm_n2 = max(sqrt(n2(1)**2 + n2(2)**2 + n2(3)**2),eps)
            n2(1) = -n2(1) / norm_n2
            n2(2) = -n2(2) / norm_n2
            n2(3) = -n2(3) / norm_n2
            c = n_int(1)*n2(1) + n_int(2)*n2(2) + n_int(3)*n2(3)
            if (abs(c) < eps) then
              u_ext(i,j,k) = n_wall(1)
              v_ext(i,j,k) = n_wall(2)
              w_ext(i,j,k) = n_wall(3)
            else if (c < 0.0_rp) then
              u_ext(i,j,k) = n_wall(1) - cot_theta * n2(1)
              v_ext(i,j,k) = n_wall(2) - cot_theta * n2(2)
              w_ext(i,j,k) = n_wall(3) - cot_theta * n2(3)
            else
              u_ext(i,j,k) = n_wall(1) + cot_theta * n2(1)
              v_ext(i,j,k) = n_wall(2) + cot_theta * n2(2)
              w_ext(i,j,k) = n_wall(3) + cot_theta * n2(3)
            end if
            norm_uext = max(sqrt(u_ext(i,j,k)**2 + v_ext(i,j,k)**2 + w_ext(i,j,k)**2),eps)
            u_ext(i,j,k) = u_ext(i,j,k) / norm_uext
            v_ext(i,j,k) = v_ext(i,j,k) / norm_uext
            w_ext(i,j,k) = w_ext(i,j,k) / norm_uext
          end if
        end do
      end do
    end do
  end subroutine compute_uextend

  subroutine advect_vof_upwind(n, dli, dtau, u_ext, v_ext, w_ext, psi)
    !
    ! the update is weighted by wgt(alphac), which rises smoothly from 0 at the
    ! outer band edge alphac = alpha_min to 1 at
    ! alphac = alpha_min + alpha_ramp*(1-alpha_min).
    !
    ! this is not a strength knob, it is a smoothness one. With the hard on/off
    ! mask the band used to have, cells just inside the edge were relaxed and
    ! cells just outside were not, so the field left behind had a kink at
    ! alphac = alpha_min -- and cmpt_norm_curv, which runs on that field, takes
    ! two derivatives across it. The resulting spurious kappa feeds
    ! sigma*kappa*grad psi and pulls fluid into the band, pitting the cells just
    ! outside it. That is the near-wall nucleation, and it is why alpha_min
    ! (which moves the edge) shifts the voids while max_pseudo_iter and dtau_cfl
    ! (which only change the strength) barely touch them. See
    ! .claude/references/contact-line-model.md.
    !
    ! the ramp is the quintic smootherstep 6x^5-15x^4+10x^3, whose first *and*
    ! second derivatives vanish at both ends. A linear or cubic ramp would only
    ! remove the jump in psi or in grad psi; kappa needs the second derivative
    ! continuous too, so the quintic is the lowest order that actually helps.
    !
    ! alpha_ramp = 0 restores the old hard edge exactly.
    !
    implicit none
    ! Input
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3):: dli
    real(rp), intent(in)              :: dtau
    real(rp), intent(in), dimension(0:,0:,0:) :: u_ext
    real(rp), intent(in), dimension(0:,0:,0:) :: v_ext
    real(rp), intent(in), dimension(0:,0:,0:) :: w_ext
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    integer  :: i, j, k
    real(rp) :: u, v, w
    real(rp) :: dpsidx, dpsidy, dpsidz
    real(rp) :: x, wgt, rampi
    logical  :: is_ramp
    !
    is_ramp = alpha_ramp > 0._rp
    rampi   = 0._rp
    if (is_ramp) rampi = 1._rp/(alpha_ramp*(1._rp-alpha_min))
    !
    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          if (alphac(i,j,k) > alpha_min .and. alphac(i,j,k) < 1._rp) then
            u = u_ext(i,j,k)
            v = v_ext(i,j,k)
            w = w_ext(i,j,k)
            if (is_ramp) then
              x   = min((alphac(i,j,k)-alpha_min)*rampi,1._rp)
              wgt = x*x*x*(x*(6._rp*x - 15._rp) + 10._rp)
            else
              wgt = 1._rp
            end if
          if (u > 0.0_rp) then
            dpsidx = (psi(i,j,k) - psi(i-1,j,k)) * dli(1)
          else
            dpsidx = (psi(i+1,j,k) - psi(i,j,k)) * dli(1)
          end if
          if (v > 0.0_rp) then
            dpsidy = (psi(i,j,k) - psi(i,j-1,k)) * dli(2)
          else
            dpsidy = (psi(i,j+1,k) - psi(i,j,k)) * dli(2)
          end if
          if (w > 0.0_rp) then
            dpsidz = (psi(i,j,k) - psi(i,j,k-1)) * dli(3)
          else
            dpsidz = (psi(i,j,k+1) - psi(i,j,k)) * dli(3)
          end if
          !
          psi(i,j,k) = psi(i,j,k) - wgt * dtau * (u*dpsidx + v*dpsidy + w*dpsidz)
          !
          end if
          
        end do
      end do
    end do
    
  end subroutine advect_vof_upwind

end module mod_extend
