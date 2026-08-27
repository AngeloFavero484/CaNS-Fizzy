# Architecture

How this repository is laid out, and how the pieces fit together at runtime.

For the *why* and the session-level context, see [`CLAUDE.md`](CLAUDE.md).
For per-topic depth, see [`.claude/references/`](.claude/references/).

---

## 1. Top-level layout

```
.
├── build.conf              # THE knob file: compiler, precision, feature defines
├── Makefile                # thin driver; globs src/*.f90, auto-generates deps
├── configs/
│   ├── compilers.mk        # per-compiler flag sets (GNU / NVIDIA / INTEL)
│   ├── flags.mk            # maps build.conf switches -> -D_DEFINES
│   └── libs.mk             # FFTW / cuFFT / HYPRE / NCCL linking
├── dependencies/           # GIT SUBMODULES — must be init'd before building
│   ├── 2decomp-fft/        # CPU pencil decomposition + transposes
│   └── cuDecomp/           # GPU pencil decomposition (NVIDIA)
├── src/                    # all solver sources (flat, no subdirs)
│   ├── main.f90            # the entire driver + timeloop
│   ├── mod_*.f90 (as *.f90)# core two-phase Navier–Stokes solver (from CaNS/Fizzy)
│   ├── prt_*.f90           # THE PARTICLE MODULE — this fork's main addition
│   ├── extend.f90          # contact-line relaxation  (fork addition)
│   ├── rotnorm.f90         # capillary force at contact line (fork addition)
│   ├── out{1,2,3}d.h90     # include files: choose what gets written
│   ├── prt_out.h90         # include file: particle output hook
│   └── poslfp/             # standalone Lagrangian-point post-processor
├── examples/
│   ├── Two_Phase/          # upstream Fizzy validation cases
│   ├── Three_Phase/        # THIS FORK's cases (particle + 2 fluids)
│   └── _CaNS-example-files/# single-phase CaNS heritage cases
├── PostPrt/                # standalone post-processing programs
│   ├── PrtPos/             # particle trajectory reader
│   └── SphereData/         # force/interpolation analysis on the sphere
├── run/                    # where the executable + input.nml + data/ live
├── tests/                  # reference regression cases
├── utils/                  # visualisation + binary-reading helpers
└── docs/                   # INFO_COMPILING / INFO_INPUT / INFO_VISU
```

**`src/` is flat and the Makefile globs it.** Adding a `.f90` file there is
enough to get it compiled — no Makefile edit needed. Module dependency order is
resolved automatically by `src/.gen-deps.awk` into `src/.depend.mk`.

---

## 2. The two layers of the solver

The code is a **one-fluid formulation**: a single velocity/pressure field on a
staggered Cartesian grid, with material properties blended by a phase field.
This fork adds a *third* material — a rigid particle — via a second indicator.

```
      psi   ∈ [0,1]      fluid-1 / fluid-2 volume fraction   (VOF or ACDI)
    alphac  ∈ [0,1]      solid / not-solid indicator          (particle)
```

Density and viscosity anywhere in the domain follow from `psi`; the particle is
imposed on top by forcing the velocity toward rigid-body motion wherever
`alphac > 0`.

### Layer A — two-phase Navier–Stokes (inherited from CaNS-Fizzy)

| file | role |
|---|---|
| `initmpi.f90`, `initgrid.f90`, `initsolver.f90` | domain decomposition, grid, Poisson setup |
| `initflow.f90`, `two_fluid.f90` | initial conditions for `u,v,w,p` and `psi` |
| `rk.f90` | momentum / scalar / phase-field time advance |
| `mom.f90` | advective + diffusive momentum terms |
| `vof_thinc_qq.f90` | **THINC/QQ VOF interface capturing (ACTIVE)** |
| `acdi.f90` | ACDI diffuse-interface alternative (**inactive** in this build) |
| `two_fluid.f90::cmpt_norm_curv_youngs` | interface normals + curvature |
| `fillps.f90`, `solver.f90`, `solver_vc.f90`, `correc.f90`, `updatep.f90` | pressure projection |
| `bound.f90` | all boundary conditions + halo exchange |
| `fft.f90`, `fftw.f90`, `solver_gpu.f90`, `workspaces.f90` | FFT machinery |
| `chkdt.f90`, `chkdiv.f90`, `sanity.f90` | stability / divergence / input checks |
| `output.f90`, `load.f90`, `post.f90` | I/O, checkpointing, statistics |

### Layer B — particle + contact line (this fork)

Gated behind `-D_PARTICLE`, with a second switch `-D_EULER` selecting the
coupling strategy.

| file | role |
|---|---|
| `prt_param.f90` | particle & collision namelists; derived spring/damper constants |
| `prt_common.f90` | `type(particle)` definitions, global arrays, `prt_InitMemo` |
| `prt_initparticles.f90` | initial placement (single explicit, or pseudo-random for `np>1`) |
| `prt_digitiser.f90` | **`alphac` from signed distance** — the diffuse solid indicator |
| `prt_initeul.f90` | rebuild `alphac` + particle normals each step (Eulerian path) |
| `prt_initvof.f90` | blend rigid-body velocity into `uphase,vphase,wphase` |
| `prt_eulint.f90` | **Eulerian IBM forcing** — the live fluid→particle coupling |
| `prt_intgr_over_sphere.f90` | volume integrals over the sphere (momentum, density) |
| `prt_intgr_nwtn_eulr.f90` | **Newton–Euler update** of particle position/velocity/spin |
| `prt_collisions.f90` | soft-sphere DEM collisions + lubrication correction |
| `prt_output.f90`, `prt_loadpart.f90` | particle MPI-IO output and restart |
| `prt_phase_indicator.f90` | sharp phase indicator (diagnostics / legacy) |
| `prt_interp_spread.f90`, `prt_forcing.f90`, `prt_kernel.f90`, `prt_coordsfp.f90` | **Lagrangian** marker path — `#if !defined(_EULER)`, *not compiled now* |

### Layer C — the contact-line extension (this fork's core novelty)

Only two files, both called directly from `main.f90`:

| file | role |
|---|---|
| `extend.f90` | `compute_uextend` builds an extension velocity from the prescribed contact angle `theta`; `advect_vof_upwind` relaxes `psi` along it in a pseudo-time loop, forcing the interface to meet the sphere at `theta` |
| `rotnorm.f90` | `rot_norm` integrates the capillary force on the particle along the numerical contact line: `Fs = -sigma * |∇psi × ∇alphac| * t̂ * dV` |

---

## 3. Runtime flow

### Initialisation (`main.f90:185–580`)

```
MPI_INIT
read_input                    (param.f90)        <- &dns &scalar &two_fluid &cudecomp
read_particle_input           (prt_param.f90)    <- &particle &collision_parameters &particle_euler
initmpi -> initgrid -> test_sanity_input -> initsolver
prt_InitMemo                  allocate alphac, norm_part*, uphase/vphase/wphase, ep(:)
  ├─ restart=F : initflow, initparticles, initeul, init2fl
  └─ restart=T : loadpart + load_one x8  (u,v,w,p,psi,fx_old,fy_old,fz_old,alphac)
vof_thinc_cmpt_phi -> cmpt_norm_curv          normals/curvature from SDF phi
[contact-line pseudo-time loop, 5 iterations]  <- same block as in the timeloop
rot_norm -> MPI_ALLREDUCE -> Fstot
chkdt -> dt
```

### Timeloop (`main.f90:587–950`)

Time advance is **Adams–Bashforth** — `do irk=1,1` with
`tm_coeff = [2+dt/dto, -dt/dto]/2`. The RK3 loop is commented out directly above.

```
istep++, time += dt
│
├─ PHASE FIELD  (if is_track_interface)
│    psio <- psi
│    initvof            blend rigid-body velocity -> uphase,vphase,wphase
│    tm_2fl             THINC/QQ transport of psi with the blended velocity
│    vof_thinc_cmpt_phi -> cmpt_norm_curv -> boundp(kappa, normx/y/z)
│    initeul            rebuild alphac + sphere normals at the new particle position
│    ┌ 5x pseudo-time contact-line relaxation ─────────────┐
│    │  compute_uextend    (theta -> extension velocity)   │   <- extend.f90
│    │  advect_vof_upwind  (relax psi toward angle theta)  │
│    │  cmpt_norm_curv + halo updates                      │
│    └──────────────────────────────────────────────────────┘
│    rot_norm -> Fs -> MPI_ALLREDUCE -> Fstot          <- rotnorm.f90 (capillary force)
│
├─ MOMENTUM
│    bulk_mean_12_stag  -> rho_av        (only for periodic + gravity directions)
│    tm                 -> u,v,w predictor   (mom.f90 via rk.f90)
│    lscale_forcing     (only if is_forced_hit)
│
├─ PARTICLE COUPLING          (_PARTICLE + _EULER)
│    eulint             Eulerian IBM force -> corrects u,v,w toward rigid body
│
├─ PRESSURE PROJECTION
│    fillps -> solver (FFT, constant-coeff) -> correc -> updatep
│
└─ PARTICLE MOTION
     intgr_nwtn_eulr    Newton–Euler: integrates ep(p)%{x,y,z,u,v,w,om*}
                        includes collisions + lubrication; writes forces_data.csv
```

Then: stopping criteria → `chkdt`/`chkdiv` every `icheck` → output every
`iout0d/1d/2d/3d` → checkpoint every `isave`.

---

## 4. Data layout and MPI

- **Pencil decomposition** in x–y (`dims(1:2)`), z is never split
  (`PENCIL_AXIS=3` here). `dims(1:2) = 0,0` lets MPI choose.
- Core fields are `(0:n(1)+1, 0:n(2)+1, 0:n(3)+1)` — **one halo cell**.
- Particle fields (`alphac`, `norm_part*`, `uphase/vphase/wphase`) use
  `(1-nh_wide : n+nh_wide)` with `nh_wide = 1` (`param.f90:24`).
- The particle module keeps its **own Cartesian communicator** `prt_comm_cart`
  (periodic in both split directions) and an 8-neighbour map `neighbor(0:8)`.
- Each particle has exactly one **master** rank (`ep(p)%mslv > 0`) and zero or
  more **slaves** (`mslv < 0`) — the ranks whose subdomain the sphere overlaps.
  Nearly every `prt_*` routine begins with the same master→slave broadcast of
  particle position/velocity, then a periodic-image correction. This boilerplate
  is repeated (not factored out) across `prt_initeul`, `prt_initvof`,
  `prt_intgr_over_sphere`, `prt_eulint`, `prt_phase_indicator`.

---

## 5. Where the boundaries between components are

| boundary | mechanism |
|---|---|
| solver ↔ particle | `alphac` (solid indicator) and `uphase/vphase/wphase` (rigid-body velocity), both in `prt_mod_common` |
| solver ↔ contact line | `psi` mutated in place by `advect_vof_upwind`; `Fstot` returned by `rot_norm` |
| particle ↔ contact line | `norm_partx/y/z` (sphere normals) consumed by both `extend.f90` and `rotnorm.f90` |
| compile-time feature gating | `build.conf` → `configs/flags.mk` → `-D_*` → `#if defined(...)` |
| runtime configuration | `input.nml`, read by `param.f90::read_input` + `prt_param.f90::read_particle_input` |

There is no abstraction layer or plugin system — coupling is by **shared module
variables**. `use prt_mod_common, only: alphac` in a solver file is the normal
and expected way to reach particle state.

---

## 6. Output artefacts

Written into `data/` relative to the executable (i.e. `run/data/`):

| file | content |
|---|---|
| `fld_{u,v,w,p,psi,alphac,...}.bin` | checkpoint fields (`isave`) |
| `*_fld_NNNNNNN.bin` + `log_visu_3d.out` | 3-D visualisation series (`iout3d`) |
| `*_slice_fld*.bin` + `log_visu_2d_*.out` | 2-D slices (`iout2d`) |
| `time.out` | istep, dt, time (`iout0d`) |
| `grid.bin`, `grid.out`, `geometry.out` | grid metadata, written once |
| particle files via `outpart` | per-particle state, MPI-IO |
| `forces_data.csv` (in `run/`, **not** `data/`) | `F_cap_ibm,F_ibm,F_inertia,F_w,F_bouy,F_cap,ep_z,ep_w` |

`out1d.h90`, `out2d.h90`, `out3d.h90`, `prt_out.h90` are **include files
compiled into `main.f90`** — editing them requires a recompile.
