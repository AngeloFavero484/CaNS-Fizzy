# Reference: coding conventions in this repository

Match the surrounding code. This file records what "surrounding code" means
here, since the repo mixes two lineages: upstream CaNS/Fizzy (clean, modern,
OpenACC-annotated) and the particle module (older, ported from a separate IBM
code, heavier and more repetitive).

---

## Layout and formatting

- **Fortran 90 free form**, `.f90` extension. Include files use `.h90`.
- **2-space indentation**. Continuation lines align under the opening construct.
- Line continuation is a trailing `&`; the continued line does **not** carry a
  leading `&` in the core solver, but **does** in some `prt_*` files. Follow the
  file you are in.
- Lowercase keywords (`subroutine`, `implicit none`, `end do`).
- `end do` / `end if` in core files; `enddo` / `endif` appear in `prt_*` files.
  Again: follow the file.
- Comment banners are a bare `!` line above and below:
  ```fortran
  !
  ! compute right-hand side of the phase field transport equation
  !
  ```
- A lone `!` is used liberally as a paragraph separator inside routine bodies.

## Module structure

Core solver:
```fortran
module mod_name
  use mpi
  use mod_types
  use mod_param, only: pi,sigma
  implicit none
  private
  public routine_name
  contains
  subroutine routine_name(...)
    implicit none          ! repeated, even though the module already has it
    ...
  end subroutine routine_name
end module mod_name
```

Particle modules:
```fortran
module prt_mod_name
#if defined(_PARTICLE)
#if defined(_EULER)          ! when the routine is coupling-specific
  use mod_types
  ...
  implicit none
  private
  public :: routine_name
  contains
  ...
#endif
#endif
end module prt_mod_name
```

Note the guards wrap the **entire module body**, leaving an empty module when the
feature is off. Keep that pattern — an empty module still satisfies `use`
statements elsewhere that are themselves guarded.

## Naming

| thing | convention | example |
|---|---|---|
| core module | `mod_<name>` in `<name>.f90` | `mod_bound` in `bound.f90` |
| particle module | `prt_mod_<name>` in `prt_<name>.f90` | `prt_mod_eulint` in `prt_eulint.f90` |
| preprocessor define | `_UPPERCASE` (leading underscore) | `_PARTICLE`, `_SDF_NORMALS` |
| logical flags | `is_*` | `is_bound`, `is_track_interface`, `is_ibm` |
| inverse quantities | `*i` suffix | `dli`, `dzci`, `dzfi`, `dti`, `r_dtcoli` |
| old/previous values | `*o` or `*_old` suffix | `psio`, `dto`, `fx_old`, `Fstot_old` |
| the two phases | `*12` suffix, index 1 = `psi=1` | `rho12(1:2)`, `mu12`, `ka12` |
| loop indices | `i,j,k` spatial, `p` particle, `l`/`lp` Lagrangian point, `nb` neighbour | |

## Numeric literals

**Always kind-suffixed:**
```fortran
0.5_rp, 1._rp, 0._rp, 2._rp/3._rp
```
Never `0.5d0`, never a bare `0.5` in an expression assigned to `real(rp)`.
(There are legacy violations in `prt_*` files — e.g. `alphac(i,j,k)=0.5*(...)`.
Do not propagate them, but do not churn the file to fix them either.)

Guard divisions with machine epsilon rather than a magic constant:
```fortran
norm = sqrt(x**2 + y**2 + z**2) + epsilon(1._rp)
```
`mod_param` also exports `eps` (which is `0._rp` in double precision, non-zero in
single) and `small = epsilon(1._rp)*10**(precision(1._rp)/2)`.

## OpenACC

Directives are interleaved throughout, even though the user currently builds with
`GPU=0`. **Keep them consistent when editing loops** — a mismatched
`enter data` / `exit data` pair is a latent GPU bug that a CPU build will not catch.

Common idioms:
```fortran
!$acc enter data copyin(psi) create(kappa,normx,normy,normz) async(1)
!$acc parallel loop collapse(3) default(present) private(uc,vc,wc) async(1)
!$acc kernels default(present) async(1)
  psio(:,:,:) = psi(:,:,:)
!$acc end kernels
!$acc wait(1)
!$acc update self(u,v,w,p,psi,kappa,s)     ! device -> host
!$acc update device(u,v,w,p) async(1)      ! host -> device
```

Rules observed here:
- `async(1)` on nearly everything, with explicit `!$acc wait` before any host
  access or MPI call.
- `!$acc update self(...)` before every output block and before the particle
  block; `!$acc update device(...)` after host-side modification.
- Reductions declared explicitly: `reduction(max:velmax)`, `reduction(+:intu,intv,intw)`.

Some `prt_*` files also carry `!$omp parallel` blocks with `shared`/`private`
clauses. Several of these list variables that no longer exist in the routine
(e.g. `cas`, `coorxmin` in `prt_initvof.f90`) — they are inert because `OPENMP`
is not enabled. **Do not enable `OPENMP=1` without auditing those clauses first.**

## MPI

```fortran
call MPI_ALLREDUCE(Fs, Fstot, 3, MPI_REAL_RP, MPI_SUM, MPI_COMM_WORLD, ierr)
```
- `MPI_REAL_RP` always, never a hard-coded datatype.
- `ierr` comes from `mod_common_mpi`; it is a module variable, not declared locally.
- Particle communication uses `prt_comm_cart`, not `MPI_COMM_WORLD`.
- Non-blocking sends/receives are collected into `arrayrequests(1:3)` and closed
  with a single `MPI_WAITALL`. The `3` is "a master may have at most 3 slaves".

## Rank-0 output

```fortran
if(myid == 0) print*, 'Time step #', istep, 'Time = ', time
```
Uppercase `PRINT *,` appears throughout the fork's additions (`PRINT *, "Fstot", Fstot`)
— that is the user's own debugging style. It is **not** guarded by `myid == 0` in
several places in `prt_intgr_nwtn_eulr.f90`, which is noisy but harmless on the
single-master path. Match the local style; do not silently remove the prints.

## Commented-out code

There is a lot of it, and it is **deliberate**:

- the RK3 driver in `main.f90:596-598`
- the entire Lagrangian IBM block in `main.f90:715-732`
- the normal-rotation contact-angle formulation in `rotnorm.f90:47-78`
- mass-conservation diagnostics in `main.f90:574-582` and `790-802`
- alternative outputs in `out2d.h90` / `out3d.h90` / `prt_out.h90`
- the capillary feedback terms in `prt_intgr_nwtn_eulr.f90`

**Do not delete any of it unless asked.** It documents alternatives the user has
tried and may return to. When adding a new alternative, comment out the old one
in the same style rather than removing it.

## Adding a new source file

Drop it in `src/`. The Makefile globs `src/*.f90` and regenerates
`src/.depend.mk` via `src/.gen-deps.awk`, so no build-system edit is needed.
Run `make clean && make` after adding a file so dependencies are re-derived.
