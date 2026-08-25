# Reference: preprocessor flags and how they gate the code

Flow: `build.conf` (plain `NAME=value`) → `configs/flags.mk` (maps to `-D_NAME`)
→ `#if defined(_NAME)` guards in `src/*.f90`.

**Any change to `build.conf` requires `make clean && make`** — the Makefile does
not track `-D` changes as dependencies.

---

## Currently set in `build.conf`

```
FCOMP=GNU
FFLAGS_OPT=1
DEBUG=1
TIMING=1
PENCIL_AXIS=3
SINGLE_PRECISION=0
CONSTANT_COEFFS_POISSON=1
INTERFACE_CAPTURING_VOF=1
SDF_NORMALS=1
GPU=0
PARTICLE=1
EULER=1
```

---

## Flags that materially change the physics

### `INTERFACE_CAPTURING_VOF` → `-D_INTERFACE_CAPTURING_VOF`

| value | active interface scheme | active file |
|---|---|---|
| `1` (**current**) | THINC/QQ volume-of-fluid | `src/vof_thinc_qq.f90` |
| `0` | ACDI diffuse interface | `src/acdi.f90` |

Consequences of `=1`:
- `mod_acdi` is **never used**. `acdi_set_epsilon`, `acdi_set_gamma`,
  `acdi_cmpt_phi`, `acdi_transport_pf` are all inside `#if !defined(...)` guards.
- `seps` and `gam` are left uninitialised/unused; `acdi_gam_factor` and
  `acdi_gam_min` have no effect.
- `log_acdi.out` is not written.
- **`psi_thickness_factor` still matters** — it is *not* ACDI-only. `param.f90:183`
  derives `vof_thinc_beta = 1/(2*psi_thickness_factor)`, so it sets the THINC
  sharpness. Its default also changes to `0.50` (from `0.51`) under this flag.

> **This is the single most common source of a wrong diagnosis in this repo.**
> Check `build.conf` before reasoning about interface behaviour.

### `SDF_NORMALS` → `-D_SDF_NORMALS` (defaults to `1` even if absent)

Selects the field used to compute normals and curvature:

| value | `cmpt_norm_curv` argument |
|---|---|
| `1` (**current**) | `phi` — the reconstructed signed-distance field |
| `0` | `psi` — the raw volume fraction |

With `=1`, `phi` is allocated and filled by `vof_thinc_cmpt_phi` (or
`acdi_cmpt_phi`) every step before the normals are taken. Note the
contact-line pseudo-loop in `main.f90:504–516` and `642–654` calls
`cmpt_norm_curv(..., psi, ...)` **directly on `psi`**, bypassing `phi` — an
inconsistency worth remembering when the contact-line normals look noisy.

### `PARTICLE` → `-D_PARTICLE`

Gates the entire `prt_*` module set, plus in the core:
- `param.f90`: `nh_wide` and `is_ibm` declarations
- `main.f90`: allocation of `alphac`, the particle init block, `eulint`,
  `intgr_nwtn_eulr`, the `alphac` checkpoint slot (`fexts(9)`)
- `out2d.h90` / `out3d.h90`: the `Alpha_C` output
- `rotnorm.f90` / `extend.f90`: `use prt_mod_common, only: alphac, norm_part*`

**`rotnorm.f90` and `extend.f90` cannot compile without `_PARTICLE`** — their
`use` statements are guarded but the loop bodies reference `alphac`
unconditionally. Turning `PARTICLE=0` in this fork will not build.

### `EULER` → `-D_EULER`

Selects the particle↔fluid coupling. Only meaningful with `_PARTICLE`.

| | `EULER=1` (**current**) | `EULER=0` |
|---|---|---|
| coupling | Eulerian IBM on `alphac` | Lagrangian forcing points on the sphere surface |
| live files | `prt_initeul`, `prt_initvof`, `prt_eulint`, `prt_digitiser` | `prt_interp_spread`, `prt_forcing`, `prt_kernel`, `prt_coordsfp` |
| velocity buffers | `uphase,vphase,wphase` | `uf,vf,wf` |
| `offset` (`prt_common`) | `(1+0.01)/dli(1)` | `sqrt(3*1.5²)/dli(1) + 0.01/dli(1)` |
| shells of markers `NL..NL4` | not allocated | allocated |

The `EULER=0` path in `main.f90` (lines ~714–732) is **additionally commented
out**, so it is doubly dead. Do not attempt to revive it without the user asking.

### `CONSTANT_COEFFS_POISSON` → `-D_CONSTANT_COEFFS_POISSON`

| value | Poisson solve |
|---|---|
| `1` (**current**) | FFT-based constant-coefficient (`solver.f90`) + pressure-splitting; allocates `pn`, `po` |
| `0` | variable-coefficient via HYPRE (`solver_vc.f90`), CPU only, much slower |

Sub-flag `BALANCED_CAPILLARY_PRESSURE_SPLIT` (not set) would allocate
`surfx_n/…/surfz_o` and pass them through `tm`.

---

## Flags that change numerics/precision

| flag | define | effect |
|---|---|---|
| `SINGLE_PRECISION` | `_SINGLE_PRECISION` | `rp = sp`; `MPI_REAL_RP = MPI_REAL`. **Note `eps` in `param.f90` changes with it**, and the machine-epsilon `psi` noise would move from `1e-16` to `~1e-7` — which would no longer be negligible. Keep at `0`. |
| `PENCIL_AXIS` | `_DECOMP_X/Y/Z` | `3` → `_DECOMP_Z`. Sets the default pencil orientation for the transposes. |
| `SPLIT_VISCOUS_DIFFUSION` | `_SPLIT_VISCOUS_DIFFUSION` | implicit/split treatment of viscous terms |
| `CAPILLARY_BRACKBILL_NORMALIZATION` | `_CAPILLARY_BRACKBILL_NORMALIZATION` | **defaults to 1** even when absent from `build.conf`. Brackbill-style CSF normalisation of the surface-tension force. |
| `MASK_DIVERGENCE_CHECK` | `_MASK_DIVERGENCE_CHECK` | suppresses the abort on `divmax > small`. Useful to *observe* a blowing-up run; never for production. |

---

## Diagnostics / build-profile flags

| flag | define | effect |
|---|---|---|
| `DEBUG` | `_DEBUG` | prints compiler + MPI version banners; enables `mod_debug` helpers. Cheap — leave on. |
| `DEBUG_SOLVER` | `_DEBUG_SOLVER` | runs `test_sanity_solver` at startup |
| `TIMING` | `_TIMING` | per-step wall-clock via `MPI_WTIME`, printed as min/max/avg |
| `SCALAR` | `_SCALAR` | scalar transport (`scal.f90`, `rk_scal`, `s` field, `fexts(10)`). **Off** — so `&scalar` in `input.nml` is inert. |
| `BOUSSINESQ_BUOYANCY` | `_BOUSSINESQ_BUOYANCY` | buoyancy from the scalar field |
| `USE_NVTX` | `_USE_NVTX` | NVIDIA profiler ranges |
| `OPENMP` | *(flag only)* | adds `-fopenmp`. Several `prt_*` files carry `!$omp` directives that are inert without it. |

---

## Compiler profiles (`configs/flags.mk`)

Exactly one of these should be `1`:

| `build.conf` | GNU flags |
|---|---|
| `FFLAGS_OPT=1` (**current**) | `-O3` |
| `FFLAGS_OPT_MAX=1` | `-Ofast -march=native` |
| `FFLAGS_DEBUG=1` | `-O0 -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -finit-real=snan -ffpe-trap=invalid -std=f2018` |
| `FFLAGS_DEBUG_MAX=1` | as above plus `-ffpe-trap=invalid,zero,overflow -finit-integer=-99999999` |

**Caution when switching to a debug profile:** `-ffpe-trap=invalid` will trap on
the `0/0` that the `+ epsilon(1._rp)` guards in `extend.f90` / `rotnorm.f90` are
*designed* to avoid but do not eliminate everywhere. A run that is fine at `-O3`
may trap under `FFLAGS_DEBUG=1`. That is usually the debug build being correct,
not a new bug.

Also: `-std=f2018` rejects some of the older constructs in the `prt_*` files that
`-O3` accepts silently. Expect warnings, and possibly errors, in a debug build.

---

## GPU

`GPU=1` requires `FCOMP=NVIDIA`, which adds
`-acc -cuda -Minfo=accel -gpu=cc60,cc70,cc80` and switches the Poisson solver to
`solver_gpu.f90` + cuFFT, and the decomposition to cuDecomp.

The `prt_*` files carry `!$acc` directives but the particle module has **not been
validated on GPU** in this fork. `main.f90` does `!$acc update self(u,v,w)`
before the particle block and `!$acc update device(...)` after — i.e. the
particle work is host-side even in a GPU build, so it would run but slowly.
