module mod_extend
  !
  use mpi
  use mod_types
  use mod_param         , only: pi,sigma,alpha_min
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
    real(rp) :: c, theta_rad, cot_theta
    real(rp), parameter :: eps = epsilon(1._rp)
    integer  :: i, j, k
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
            norm_n1 = sqrt(n1(1)**2 + n1(2)**2 + n1(3)**2) + eps
            n1(1) = -n1(1) / norm_n1
            n1(2) = -n1(2) / norm_n1
            n1(3) = -n1(3) / norm_n1
            ! Vettore tangente al solido, ortogonale alla linea di contatto
            n2(1) = n1(2)*n_wall(3) - n1(3)*n_wall(2)
            n2(2) = n1(3)*n_wall(1) - n1(1)*n_wall(3)
            n2(3) = n1(1)*n_wall(2) - n1(2)*n_wall(1)
            !
            norm_n2 = sqrt(n2(1)**2 + n2(2)**2 + n2(3)**2) + eps
            n2(1) = -n2(1) / norm_n2
            n2(2) = -n2(2) / norm_n2
            n2(3) = -n2(3) / norm_n2
            c = n_int(1)*n2(1) + n_int(2)*n2(2) + n_int(3)*n2(3)
            theta_rad = theta * pi / 180.0_rp
            if (abs(c) < eps) then
              u_ext(i,j,k) = n_wall(1)
              v_ext(i,j,k) = n_wall(2)
              w_ext(i,j,k) = n_wall(3)
            else if (c < 0.0_rp) then
              cot_theta = cos(pi-theta_rad) / sin(pi-theta_rad)
              u_ext(i,j,k) = n_wall(1) - cot_theta * n2(1)
              v_ext(i,j,k) = n_wall(2) - cot_theta * n2(2)
              w_ext(i,j,k) = n_wall(3) - cot_theta * n2(3)
            else
              cot_theta = cos(pi-theta_rad) / sin(pi-theta_rad)
              u_ext(i,j,k) = n_wall(1) + cot_theta * n2(1)
              v_ext(i,j,k) = n_wall(2) + cot_theta * n2(2)
              w_ext(i,j,k) = n_wall(3) + cot_theta * n2(3)
            end if
            norm_uext = sqrt(u_ext(i,j,k)**2 + v_ext(i,j,k)**2 + w_ext(i,j,k)**2) + eps
            u_ext(i,j,k) = u_ext(i,j,k) / norm_uext
            v_ext(i,j,k) = v_ext(i,j,k) / norm_uext
            w_ext(i,j,k) = w_ext(i,j,k) / norm_uext
          end if
        end do
      end do
    end do
  end subroutine compute_uextend

  subroutine advect_vof_upwind(n, dli, dtau, u_ext, v_ext, w_ext, psi)
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

    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          if (alphac(i,j,k) > alpha_min .and. alphac(i,j,k) < 1._rp) then
            u = u_ext(i,j,k)
            v = v_ext(i,j,k)
            w = w_ext(i,j,k)
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
          psi(i,j,k) = psi(i,j,k) - dtau * (u*dpsidx + v*dpsidy + w*dpsidz)
          !
          end if
          
        end do
      end do
    end do
    
  end subroutine advect_vof_upwind

end module mod_extend
