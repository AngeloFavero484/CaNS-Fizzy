! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_common_mpi
#if defined(_PARTICLE)
  use mpi
  use mod_types
  use decomp_2d, only: decomp_info
#endif
  implicit none
  public
  integer :: myid,ierr
  integer :: halo(3)
  integer :: ipencil_axis
#if defined(_PARTICLE)
  integer :: comm_cart
  integer :: status(MPI_STATUS_SIZE)
  integer :: prt_comm_cart
  integer :: halo_wide(3)
  integer :: right,rightfront,front,leftfront,left,leftback,back,rightback
  real(rp) :: boundleftmyid,boundfrontmyid
  type(decomp_info) :: dinfo_ptdma
#endif
end module mod_common_mpi
