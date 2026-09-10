! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_rotnorm
  !
  use mpi
  use mod_types
  use mod_param         , only: pi,sigma
#if defined(_PARTICLE)
    use prt_mod_common    , only: alphac,norm_partx,norm_party,norm_partz
#endif
  implicit none
  private
  public rot_norm
  contains
  subroutine rot_norm(n,dli,dzci,psi,theta,is_bound,normx,normy,normz,kappa,Fs)
    implicit none
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(in   ), dimension(3)        :: dli
    real(rp), intent(in   ), dimension(0:)       :: dzci
    real(rp), intent(in   )                      :: theta
    logical , intent(in   ), dimension(0:1,3  )  :: is_bound
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: normx,normy,normz,kappa
    real(rp), intent(out)  , dimension(3)        :: Fs
    integer  :: i,j,k
    real(rp) :: normt
    real(rp) :: dot
    real(rp) :: ntx,nty,ntz
    real(rp) :: theta_rad
    real(rp) :: psixp,psixm,psiyp,psiym,psizp,psizm,dpsidx,dpsidy,dpsidz,norm
    real(rp) :: alphaxp,alphaxm,alphayp,alphaym,alphazp,alphazm,dalphadx,dalphady,dalphadz
    real(rp) :: vecti,vectj,vectk,prod
    real(rp) :: tx,ty,tz
    real(rp) :: norm_t
    real(rp) :: norm_alpha
    real(rp) :: dA,dV
    real(rp) :: perimetro_numerico
    dV=(1/(dli(1)))*(1/(dli(2)))*(1/(dli(3)))
    !
    theta_rad = theta*(pi/180)
      !
!      do k=0,n(3)
!        do j=0,n(2)
!          do i=0,n(1)
!            if (alphac(i,j,k)>0._rp .and. alphac(i,j,k)<1._rp) then
!              dot = normx(i,j,k)*norm_partx(i,j,k) + normy(i,j,k)*norm_party(i,j,k) + normz(i,j,k)*norm_partz(i,j,k)
!              ntx = normx(i,j,k) - dot*norm_partx(i,j,k)
!              nty = normy(i,j,k) - dot*norm_party(i,j,k)
!              ntz = normz(i,j,k) - dot*norm_partz(i,j,k)
!              !
!              normt = sqrt(ntx**2 + nty**2 + ntz**2) + epsilon(1._rp)
!              !
!              ntx = ntx/normt
!              nty = nty/normt
!              ntz = ntz/normt
!              !
!              normx(i,j,k) = -cos(theta_rad)*norm_partx(i,j,k) + sin(theta_rad)*ntx
!              normy(i,j,k) = -cos(theta_rad)*norm_party(i,j,k) + sin(theta_rad)*nty
!              normz(i,j,k) = -cos(theta_rad)*norm_partz(i,j,k) + sin(theta_rad)*ntz
!            end if
!          end do
!        end do
!      end do
!      !
!      do k=1,n(3)
!        do j=1,n(2)
!          do i=1,n(1)
!            if (alphac(i,j,k)>0._rp .and. alphac(i,j,k)<1._rp) then
!              kappa(i,j,k) = - ( (normx(i+1,j,k)-normx(i-1,j,k))*0.5*dli(1) &
!                             +   (normy(i,j+1,k)-normy(i,j-1,k))*0.5*dli(2) &
!                             +   (normz(i,j,k+1)-normz(i,j,k-1))*0.5*dzci(k) )
!           end if
!         end do
!        end do
!      end do
      !
      Fs(:)=0
      perimetro_numerico=0
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if (alphac(i,j,k)>0._rp .and. alphac(i,j,k)<1._rp ) then
              psixp = 0.5*(psi(i+1,j,k)+psi(i  ,j,k))
              psixm = 0.5*(psi(i  ,j,k)+psi(i-1,j,k))
              psiyp = 0.5*(psi(i,j+1,k)+psi(i,j  ,k))
              psiym = 0.5*(psi(i,j  ,k)+psi(i,j-1,k))
              psizp = 0.5*(psi(i,j,k+1)+psi(i,j,k  ))
              psizm = 0.5*(psi(i,j,k  )+psi(i,j,k-1))
              !
              alphaxp = 0.5*(alphac(i+1,j,k)+alphac(i  ,j,k))
              alphaxm = 0.5*(alphac(i  ,j,k)+alphac(i-1,j,k))
              alphayp = 0.5*(alphac(i,j+1,k)+alphac(i,j  ,k))
              alphaym = 0.5*(alphac(i,j  ,k)+alphac(i,j-1,k))
              alphazp = 0.5*(alphac(i,j,k+1)+alphac(i,j,k  ))
              alphazm = 0.5*(alphac(i,j,k  )+alphac(i,j,k-1))
              !
              dpsidx = (psixp-psixm)*dli(1)
              dpsidy = (psiyp-psiym)*dli(2)
              dpsidz = (psizp-psizm)*dli(3)
              !
              dalphadx = (alphaxp-alphaxm)*dli(1)
              dalphady = (alphayp-alphaym)*dli(2)
              dalphadz = (alphazp-alphazm)*dli(3)
              !
! This formulation uses the actual normals to the sphere
              norm_alpha=sqrt(dalphadx**2+dalphady**2+dalphadz**2)
              dalphadx = norm_partx(i,j,k)*norm_alpha
              dalphady = norm_party(i,j,k)*norm_alpha
              dalphadz = norm_partz(i,j,k)*norm_alpha
              !
              vecti = dpsidy*dalphadz-dpsidz*dalphady
              vectj = dpsidz*dalphadx-dpsidx*dalphadz
              vectk = dpsidx*dalphady-dpsidy*dalphadx
              !
              prod = sqrt(vecti**2+vectj**2+vectk**2)
              !
              tx = normy(i,j,k)*vectk - normz(i,j,k)*vectj
              ty = normz(i,j,k)*vecti - normx(i,j,k)*vectk
              tz = normx(i,j,k)*vectj - normy(i,j,k)*vecti
              !
              norm_t = sqrt(tx**2 + ty**2 + tz**2) + epsilon(1._rp)
              !
              tx = tx / norm_t
              ty = ty / norm_t
              tz = tz / norm_t
              !
              Fs(1) = Fs(1) - sigma * prod * tx * dV
              Fs(2) = Fs(2) - sigma * prod * ty * dV
              Fs(3) = Fs(3) - sigma * prod * tz * dV
              perimetro_numerico = perimetro_numerico + prod * dV
            end if
          end do
        end do
      end do
      !
!      Fs=0
!      PRINT *, "Perimetro", perimetro_numerico
!      PRINT *, "Forcex", Fs(1)
!      PRINT *, "Forcey", Fs(2)
!      PRINT *, "Forcez", Fs(3)
  end subroutine rot_norm
  !
end module mod_rotnorm
