module prt_mod_interp_spread
#if defined(_PARTICLE)
#if !defined(_EULER)
  !
  use mpi
  use mod_types
  use mod_param       , only: nh_wide,ng,dl,dli,l
  use mod_common_mpi  , only: prt_comm_cart
  use prt_mod_common  , only: ep,dVlagr,dVeul,pmax,nl,nla
  use prt_mod_kernel  , only: kernel
  !
  implicit none
  !
  private
  !
  public eulr2lagr,lagr2eulr
  !
  contains
  !
  subroutine eulr2lagr(n,lo,uf,vf,wf)
    implicit none
    integer , intent(in), dimension(3) :: n,lo
    real(rp), intent(inout), dimension(1-nh_wide:,1-nh_wide:,1-nh_wide:) :: uf,vf,wf
    integer :: i,j,k,ll,p
    integer :: ilow,jlow,klow, ilows,jlows,klows    !,ihigh,jhigh,khigh
    real(rp) :: coorxfp,cooryfp,coorzfp
    real(rp) :: phix(1:3), phiy(1:3), phiz(1:3), phix_half(1:3), phiy_half(1:3), phiz_half(1:3)
    integer :: val_k_min, val_k_max
    !
    val_k_min = 1
    val_k_max = n(3)
    !
    ! third step: perform partial integration.
    !    
    do ll=1,NL
      do p=1,pmax
        ep(p)%dudtl(ll) = 0._rp
        ep(p)%dvdtl(ll) = 0._rp
        ep(p)%dwdtl(ll) = 0._rp
        !    sumu(:,p)      = 0._rp
      enddo
    enddo
    !
!    call updthalos_ibm(dudtf,1,1)
!    call updthalos_ibm(dvdtf,1,1)
!    call updthalos_ibm(dwdtf,1,1)
!    call updthalos_ibm(dudtf,2,1)
!    call updthalos_ibm(dvdtf,2,1)
!    call updthalos_ibm(dwdtf,2,1)
!    !
    !$omp parallel default(none)             &
    !$omp shared(ep,nla,pmax)                &
    !$omp shared(uf,vf,wf)          &
    !$omp shared(lo,l,ng,dl,dli,val_k_min,val_k_max) &
    !$omp private(p,nbrecv,ll,coorxfp,cooryfp,coorzfp,ilow,ihigh,jlow,jhigh,klow,khigh) &
    !$omp private(i,j,k,nbsend,nb) &
    !$omp private(coorzs,coorz,kernelzs,kernelz,coorys,coory,kernelys,kernely,coorxs,coorx,kernelxs,kernelx) &
    !$omp private(phix,phiy,phiz,phix_half,phiy_half,phiz_half) &
    !$omp private(ilows,jlows,klows)
    !$omp do 
    do p=1,pmax
      if (ep(p)%mslv /= 0) then
        do ll=1,nla(p)
          ! description of the coordinates on a local frame of reference with normalized space steps, 
          ! taking into account the periodicities in the axis x and y 
          coorxfp = ep(p)%xfp(ll)*dli(1) - floor((ep(p)%xfp(ll)-0.5_rp*dl(1))/l(1))*(ng(1)) - (lo(1)-1)
          cooryfp = ep(p)%yfp(ll)*dli(2) - floor((ep(p)%yfp(ll)-0.5_rp*dl(2))/l(2))*(ng(2)) - (lo(2)-1)
          coorzfp = ep(p)%zfp(ll)*dli(3)

          !          
!!! DEBUGGING SECTION
!          isperiodx = 0._rp
!          isperiody = 0._rp 
!          if (ep(p)%xfp(ll) < 0._rp+0.5_rp*dl(1)) isperiodx =  1._rp
!          if (ep(p)%xfp(ll) > l(1) +0.5_rp*dl(1)) isperiodx = -1._rp
!          if (ep(p)%yfp(ll) < 0._rp+0.5_rp*dl(2)) isperiody =  1._rp
!          if (ep(p)%yfp(ll) > l(2) +0.5_rp*dl(2)) isperiody = -1._rp
!          if (abs(coorxfp - (ep(p)%xfp(ll)+isperiodx*l(1)-boundleftmyid )*dli(1)) > (1d-12)) then
!            print*,'CALCULATION MISTAKE ON COORXFP ',coorxfp,(ep(p)%xfp(ll)+isperiodx*l(1)-boundleftmyid )*dli(1)
!          endif
!          !   
!          if (abs(cooryfp - (ep(p)%yfp(ll)+isperiody*l(2)-boundfrontmyid)*dli(2)) > (1d-12)) then
!            print*,'CALCULATION MISTAKE ON COORYFP ',cooryfp,(ep(p)%yfp(ll)+isperiody*l(2)-boundfrontmyid)*dli(2)
!          endif
!!! END DEBUGGING SECTION
          !
          ! Indexes for the beginning of the stencil for coorxs, coorys,coorzs
          ! (the only coordinates for which the kernel calculation is non null) 
          ilows =  nint(coorxfp) -1
          jlows =  nint(cooryfp) -1
          klows =  max( nint(coorzfp) -1 , val_k_min )
          ! indexes for the beginning of the stencil for coorx, coory, coorz
          ilow = floor(coorxfp)
          jlow = floor(cooryfp)
          klow = max( min( floor(coorzfp) , val_k_max-2 ) , val_k_min )
          !             
!!!!THIRD VERSION :   
          phix(1) = kernel(  1._rp*ilows        - coorxfp)
          phix(2) = kernel(  1._rp*ilows +1._rp - coorxfp)
          phix(3) = kernel(  1._rp*ilows +2._rp - coorxfp)
          !
          phiy(1) = kernel(  1._rp*jlows        - cooryfp)
          phiy(2) = kernel(  1._rp*jlows +1._rp - cooryfp)
          phiy(3) = kernel(  1._rp*jlows +2._rp - cooryfp)
          !
          phiz(1) = kernel(  1._rp*klows        - coorzfp)
          phiz(2) = kernel(  1._rp*klows +1._rp - coorzfp)
          phiz(3) = kernel(  1._rp*klows +2._rp - coorzfp)
          !
          phix_half(1) = kernel(  1._rp*ilow        -0.5 - coorxfp)
          phix_half(2) = kernel(  1._rp*ilow +1._rp -0.5 - coorxfp)
          phix_half(3) = kernel(  1._rp*ilow +2._rp -0.5 - coorxfp)
          !
          phiy_half(1) = kernel(  1._rp*jlow        -0.5_rp - cooryfp)
          phiy_half(2) = kernel(  1._rp*jlow +1._rp -0.5_rp - cooryfp)
          phiy_half(3) = kernel(  1._rp*jlow +2._rp -0.5_rp - cooryfp)
          !
          phiz_half(1) = kernel(  1._rp*klow        -0.5_rp - coorzfp)
          phiz_half(2) = kernel(  1._rp*klow +1._rp -0.5_rp - coorzfp)
          phiz_half(3) = kernel(  1._rp*klow +2._rp -0.5_rp - coorzfp)
          !
          do k=0,2
            do j=0,2
              do i=0,2
                ep(p)%dudtl(ll) = ep(p)%dudtl(ll) + uf(ilows+i,jlow+j,klow+k) * phix(i+1) * phiy_half(j+1) * phiz_half(k+1)
              enddo
            enddo
          enddo
          !
          do k=0,2
            do j=0,2
              do i=0,2
                ep(p)%dvdtl(ll) = ep(p)%dvdtl(ll) + vf(ilow+i,jlows+j,klow+k) * phix_half(i+1) * phiy(j+1) * phiz_half(k+1)
              enddo
            enddo
          enddo
          !
          do k=0,2
            do j=0,2
              do i=0,2
                ep(p)%dwdtl(ll) = ep(p)%dwdtl(ll) + wf(ilow+i,jlow+j,klows+k) * phix_half(i+1) * phiy_half(j+1) * phiz(k+1)
              enddo
            enddo
          enddo
!!!END THIRD VERSION

!!! FIRST VERSION
             !isperiodx = 0.                                             
             !isperiody = 0.                                             
             !if (ap(p)%xfp(l).lt.0.+0.5*dx) isperiodx =  1.             
             !if (ap(p)%xfp(l).ge.lx+0.5*dx) isperiodx = -1.             
             !if (ap(p)%yfp(l).lt.0.+0.5*dy) isperiody =  1.             
             !if (ap(p)%yfp(l).ge.ly+0.5*dy) isperiody = -1.             
             !      isout = .false.                                      
             !coorxfp = (ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi   
             !      if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
             !cooryfp = (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi 
             !      if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
             !      if (.not.isout) then
             !ilow  = nint( coorxfp - 1.5)
             !if( ((1.*ilow     ) - coorxfp) .lt. -1.5) ilow  = ilow + 1
             !ihigh = nint( coorxfp + 1.5)
             !if( ((1.*ihigh-0.5) - coorxfp) .gt.  1.5) ihigh = ihigh - 1
             !jlow  = nint( cooryfp - 1.5)
             !if( ((1.*jlow     ) - cooryfp) .lt. -1.5) jlow  = jlow + 1
             !jhigh = nint( cooryfp + 1.5)
             !if( ((1.*jhigh-0.5) - cooryfp) .gt.  1.5) jhigh = jhigh - 1
             !klow  = nint( coorzfp - 1.5)
             !if( ((1.*klow     ) - coorzfp) .lt. -1.5) klow  = klow + 1
             !khigh = nint( coorzfp + 1.5)
             !if( ((1.*khigh-0.5) - coorzfp) .gt.  1.5) khigh = khigh - 1
             !        if (ilow .lt. -1) ilow = -1
             !        if (ihigh .gt. i1+1) ihigh = i1+1
             !        if (jlow .lt. -1) jlow = -1
             !        if (jhigh .gt. j1+1) jhigh = j1+1
             !if (klow .lt. 1) klow = 1 
             !if (khigh .gt. kmax) khigh = kmax
 
             !do k=klow,khigh
             !   coorzs   = (1.*k)-coorzfp
             !   coorz    = coorzs-0.5           !vert.    distance in grid points
             !   kernelzs = kernel(coorzs)
             !   kernelz  = kernel(coorz)
             !   coeffx1  = kernelz
             !   coeffy1  = kernelz
             !   coeffz1  = kernelzs
             !   do j=jlow,jhigh
             !      coorys   = (1.*j)-cooryfp
             !      coory    = coorys-0.5         !spanw.   distance in grid points
             !      kernelys = kernel(coorys)
             !      kernely  = kernel(coory)
             !      coeffx2  = coeffx1*kernely
             !      coeffy2  = coeffy1*kernelys
             !      coeffz2  = coeffz1*kernely
             !      do i=ilow,ihigh
             !         coorxs   = (1.*i)-coorxfp
             !         coorx    = coorxs-0.5       !streamw. distance in grid points
             !         kernelxs = kernel(coorxs)
             !         kernelx  = kernel(coorx)
             !         coeffx3  = coeffx2*kernelxs
             !         coeffy3  = coeffy2*kernelx
             !         coeffz3  = coeffz2*kernelx
             !         ap(p)%dudtl(l) = ap(p)%dudtl(l) + dudtf(i,j,k)*coeffx3
             !         ap(p)%dvdtl(l) = ap(p)%dvdtl(l) + dvdtf(i,j,k)*coeffy3
             !         ap(p)%dwdtl(l) = ap(p)%dwdtl(l) + dwdtf(i,j,k)*coeffz3
             !         !              sumu(l,p) = sumu(l,p) + kernelxs*kernely*kernelz
             !      enddo
             !   enddo
             !enddo
!!! END FIRST VERSION                                                                                                                                                   
             !!      endif
        enddo ! do l=
      endif
    enddo
    !$omp end parallel
    !
    return
  end subroutine eulr2lagr
  !
  subroutine lagr2eulr(dt,n,lo,dudtold,dvdtold,dwdtold,uf,vf,wf)
    implicit none
    real(rp), intent(in) :: dt
    integer, intent(in), dimension(3) :: n,lo
    real(rp), intent(in), dimension(0:,0:,0:) :: dudtold,dvdtold,dwdtold
    real(rp), intent(inout), dimension(1-nh_wide:,1-nh_wide:,1-nh_wide:) :: uf,vf,wf
    integer :: i,j,k,ll,p
    integer :: ilow,jlow,klow, ilows,jlows,klows    !, ihigh,jhigh,khigh
    real(rp) :: coorxfp,cooryfp,coorzfp
    real(rp) :: dVlagrdVeuli
    real(rp) :: coeffx0,coeffy0,coeffz0    !, &
    real(rp) :: phix(1:3), phiy(1:3), phiz(1:3), phix_half(1:3), phiy_half(1:3), phiz_half(1:3)
    integer :: val_k_min, val_k_max
    !
    val_k_min = 1
    val_k_max = n(3)
    !
    dVlagrdVeuli = dVlagr/dVeul
    !
    uf(:,:,:) = 0._rp
    vf(:,:,:) = 0._rp
    wf(:,:,:) = 0._rp
    !
    coeffx0 = dVlagrdVeuli*dt
    coeffy0 = dVlagrdVeuli*dt
    coeffz0 = dVlagrdVeuli*dt
!!!!!!!!!    forceytot = 0._rp
    !
    !$omp parallel default(none)             &
    !$omp shared(ep,nla,pmax)                &
    !$omp shared(uf,vf,wf)          &     
    !$omp shared(coeffx0,coeffy0,coeffz0)    &
    !$omp shared(val_k_min,val_k_max,dl,l,dli,ng) &
    !$omp private(p,nbrecv,ll,lo,coorxfp,cooryfp,coorzfp,ilow,ihigh,jlow,jhigh,klow,khigh) &
    !$omp private(i,j,k,nbsend,nb) &
    !$omp private(coorzs,coorz,kernelzs,kernelz,coorys,coory,kernelys,kernely) &
    !$omp private(coorxs,coorx,kernelxs,kernelx,forcex_sc,forcey_sc,forcez_sc) &
    !$omp private(phix,phiy,phiz,phix_half,phiy_half,phiz_half) &
    !$omp private(ilows,jlows,klows) &
    !$omp reduction(+:forceytot)
    !$omp do 
    do p=1,pmax
      if (ep(p)%mslv /= 0) then
        do ll=1,nla(p)
          !!!THIRD VERSION
          coorxfp = ep(p)%xfp(ll)*dli(1) - floor((ep(p)%xfp(ll)-0.5_rp*dl(1))/l(1))*(ng(1)) - (lo(1)-1)
          cooryfp = ep(p)%yfp(ll)*dli(2) - floor((ep(p)%yfp(ll)-0.5_rp*dl(2))/l(2))*(ng(2)) - (lo(2)-1)
          coorzfp = ep(p)%zfp(ll)*dli(3)
          !             
          ! Indexes for the beginning of the stencil for coorxs, coorys,coorzs 
          ! (the only coordinates for which the kernel calculation is non null)
          ilows =  nint(coorxfp) -1
          jlows =  nint(cooryfp) -1
          klows =  max( nint(coorzfp) -1 , val_k_min )
          ! indexes for the beginning of the stencil for coorx, coory, coorz 
          ilow = floor(coorxfp)
          jlow = floor(cooryfp)
          klow = max( min( floor(coorzfp) , val_k_max-2 ) , val_k_min )
          !
          phix(1) = kernel(  1._rp*ilows    - coorxfp)
          phix(2) = kernel(  1._rp*ilows +1._rp - coorxfp)
          phix(3) = kernel(  1._rp*ilows +2._rp - coorxfp)
          !        
          phiy(1) = kernel(  1._rp*jlows    - cooryfp)
          phiy(2) = kernel(  1._rp*jlows +1._rp - cooryfp)
          phiy(3) = kernel(  1._rp*jlows +2._rp - cooryfp)
          !
          phiz(1) = kernel(  1._rp*klows    - coorzfp)
          phiz(2) = kernel(  1._rp*klows +1._rp - coorzfp)
          phiz(3) = kernel(  1._rp*klows +2._rp - coorzfp)
          !
          phix_half(1) = kernel(  1._rp*ilow    -.5_rp - coorxfp)
          phix_half(2) = kernel(  1._rp*ilow +1._rp -.5_rp - coorxfp)
          phix_half(3) = kernel(  1._rp*ilow +2._rp -.5_rp - coorxfp)
          !
          phiy_half(1) = kernel(  1._rp*jlow    -.5_rp - cooryfp)
          phiy_half(2) = kernel(  1._rp*jlow +1._rp -.5_rp - cooryfp)
          phiy_half(3) = kernel(  1._rp*jlow +2._rp -.5_rp - cooryfp)
          !
          phiz_half(1) = kernel(  1._rp*klow    -.5_rp - coorzfp)
          phiz_half(2) = kernel(  1._rp*klow +1._rp -.5_rp - coorzfp)
          phiz_half(3) = kernel(  1._rp*klow +2._rp -.5_rp - coorzfp)
          !
          do k=0,2
            do j=0,2
              do i=0,2
                !$omp atomic
                uf(ilows+i,jlow+j,klow+k) = uf(ilows+i,jlow+j,klow+k) + &
                                            ep(p)%fxl(ll) * coeffx0 * phix(i+1) * phiy_half(j+1) * phiz_half(k+1)
              enddo
            enddo
          enddo
          !
          do k=0,2
            do j=0,2
              do i=0,2
                !$omp atomic
                vf(ilow+i,jlows+j,klow+k) = vf(ilow+i,jlows+j,klow+k) + &
                                            ep(p)%fyl(ll) * coeffy0 * phix_half(i+1) * phiy(j+1) * phiz_half(k+1)
!!!!!!!!!!!!                forceytot = forceytot + ep(p)%fyl(ll) * coeffy0 * phix_half(i+1) * phiy(j+1) * phiz_half(k+1)
              enddo
            enddo
          enddo
          !
          do k=0,2
            do j=0,2
              do i=0,2
                !$omp atomic
                wf(ilow+i,jlow+j,klows+k) = wf(ilow+i,jlow+j,klows+k) + &
                                            ep(p)%fzl(ll) * coeffz0 * phix_half(i+1) * phiy_half(j+1) * phiz(k+1)
              enddo
            enddo
          enddo
          !
             !!!FIRST VERSION
             !isperiodx = 0.
             !isperiody = 0.
             !if (ap(p)%xfp(l).lt.0.+0.5*dx) isperiodx =  1.
             !if (ap(p)%xfp(l).ge.lx+0.5*dx) isperiodx = -1.
             !if (ap(p)%yfp(l).lt.0.+0.5*dy) isperiody =  1.
             !if (ap(p)%yfp(l).ge.ly+0.5*dy) isperiody = -1.
             !!      isout = .false.
             !coorxfp = (ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi
             !!      if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
             !cooryfp = (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi
             !!      if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
             !!      if (.not.isout) then
             !coorzfp = ap(p)%zfp(l)*dzi
             !ilow  = nint( coorxfp - 1.5)
             !if( ((1.*ilow     ) - coorxfp) .lt. -1.5) ilow  = ilow + 1
             !ihigh = nint( coorxfp + 1.5)
             !if( ((1.*ihigh-0.5) - coorxfp) .gt.  1.5) ihigh = ihigh - 1
             !jlow  = nint( cooryfp - 1.5)
             !if( ((1.*jlow     ) - cooryfp) .lt. -1.5) jlow  = jlow + 1
             !jhigh = nint( cooryfp + 1.5)
             !if( ((1.*jhigh-0.5) - cooryfp) .gt.  1.5) jhigh = jhigh - 1
             !klow  = nint( coorzfp - 1.5)
             !if( ((1.*klow     ) - coorzfp) .lt. -1.5) klow  = klow + 1
             !khigh = nint( coorzfp + 1.5)
             !if( ((1.*khigh-0.5) - coorzfp) .gt.  1.5) khigh = khigh - 1
             !!        if (ilow .lt. -1) ilow = -1
             !!        if (ihigh .gt. i1+1) ihigh = i1+1
             !!        if (jlow .lt. -1) jlow = -1 
             !!        if (jhigh .gt. j1+1) jhigh = j1+1 
             !if (klow .lt. 1) klow = 1
             !if (khigh .gt. kmax) khigh = kmax
            ! 
            ! do k=klow,khigh
            !    coorzs   = (1.*k)-coorzfp
            !    coorz    = coorzs-0.5           !vert.    distance in grid points
            !    kernelzs = kernel(coorzs)
            !    kernelz  = kernel(coorz)
            !    coeffx1   = coeffx0*kernelz
            !    coeffy1   = coeffy0*kernelz
            !    coeffz1   = coeffz0*kernelzs
            !    do j=jlow,jhigh
            !       coorys   = (1.*j)-cooryfp
            !       coory    = coorys-0.5         !spanw.   distance in grid points
            !       kernelys = kernel(coorys)
            !       kernely  = kernel(coory)
            !       coeffx2   = coeffx1*kernely
            !       coeffy2   = coeffy1*kernelys
            !       coeffz2   = coeffz1*kernely
            !       do i=ilow,ihigh
            !          coorxs   = (1.*i)-coorxfp
            !          coorx    = coorxs-0.5       !streamw. distance in grid points
            !          kernelxs = kernel(coorxs)
            !          kernelx  = kernel(coorx)
            !          coeffx3  = coeffx2*kernelxs
            !          coeffy3  = coeffy2*kernelx
            !          coeffz3  = coeffz2*kernelx
            !          !  !$omp atomic
            !          !              forcex(i,j,k) = forcex(i,j,k) + ap(p)%fxl(l)*kernelxs*kernely*kernelz*dVlagrdVeuli
            !          !  !$omp atomic
            !          !              forcey(i,j,k) = forcey(i,j,k) + ap(p)%fyl(l)*kernelx*kernelys*kernelz*dVlagrdVeuli
            !          !  !$omp atomic
            !          !              forcez(i,j,k) = forcez(i,j,k) + ap(p)%fzl(l)*kernelx*kernely*kernelzs*dVlagrdVeuli
            !          !$omp atomic
            !             !              forcex(i,j,k) = forcex(i,j,k) + ap(p)%fxl(l)*kernelxs*kernely*kernelz*dVlagrdVeuli
            !          !  !$omp atomic
            !          !              forcey(i,j,k) = forcey(i,j,k) + ap(p)%fyl(l)*kernelx*kernelys*kernelz*dVlagrdVeuli
            !          !  !$omp atomic
            !          !              forcez(i,j,k) = forcez(i,j,k) + ap(p)%fzl(l)*kernelx*kernely*kernelzs*dVlagrdVeuli
            !          !$omp atomic
            !          dudtf(i,j,k) = dudtf(i,j,k) + ap(p)%fxl(l)*coeffx3
            !          !$omp atomic
            !          dvdtf(i,j,k) = dvdtf(i,j,k) + ap(p)%fyl(l)*coeffy3
            !          !$omp atomic
            !          dwdtf(i,j,k) = dwdtf(i,j,k) + ap(p)%fzl(l)*coeffz3
            !       enddo
            !    enddo
            ! enddo
            ! !       endif
        enddo !do l=
      endif
    enddo
    !$omp end parallel
    !
!    call updthalos_ibm(dudtf ,1,2)
!    call updthalos_ibm(dvdtf ,1,2)
!    call updthalos_ibm(dwdtf ,1,2)
!    call updthalos_ibm(dudtf ,2,2)
!    call updthalos_ibm(dvdtf ,2,2)
!    call updthalos_ibm(dwdtf ,2,2)
    !
!!!!!!!!!!    call MPI_ALLREDUCE(forceytot,forceytot_all,1,MPI_REAL_RP,mpi_sum,prt_comm_cart,ierr)
!!!!!!!!!!    forceytot = forceytot_all/dt/(1._rp*ng(1)*ng(2)*ng(3)) !+ wallshearnew
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          uf(i,j,k) = uf(i,j,k) + dudtold(i,j,k)
          vf(i,j,k) = vf(i,j,k) + dvdtold(i,j,k) !+ (-forceytot)*dt
          wf(i,j,k) = wf(i,j,k) + dwdtold(i,j,k)
        enddo
      enddo
    enddo
    !call updthalos_ibm(forcex,1,2)
    !call updthalos_ibm(forcey,1,2)
    !call updthalos_ibm(forcez,1,2)
    !call updthalos_ibm(forcex,2,2)
    !call updthalos_ibm(forcey,2,2)
    !call updthalos_ibm(forcez,2,2)
    !
!    !$omp workshare
!    dvdtf(1:imax,1:jmax,1:kmax) = dvdtf(1:imax,1:jmax,1:kmax) + 1.0*(bulk_v_sup-v_bulk)
!    !$omp end workshare
    !
    return
  end subroutine lagr2eulr
  !
#endif
#endif
end module prt_mod_interp_spread
