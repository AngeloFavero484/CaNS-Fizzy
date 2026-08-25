# Reference: what this fork adds on top of CaNS-Fizzy

Baseline: `origin/main` = https://github.com/CaNS-World/CaNS-Fizzy
Fork: `myfork/main` = https://github.com/AngeloFavero484/CaNS-Fizzy

```
git diff --stat origin/main..main
  97 files changed, 12044 insertions(+), 81 deletions(-)
```

The fork is **additive**: only ~81 lines of upstream code were removed. Almost
everything is new files or new blocks guarded by `#if defined(_PARTICLE)`.

---

## Commit history of the fork (oldest → newest)

| commit | what |
|---|---|
| `ae2929a` | **Extend CaNS-Fizzy with finite-size particle module and contact-line model** — the big one; imports the whole `prt_*` set, `extend.f90`, `rotnorm.f90`, `PostPrt/`, `src/poslfp/` |
| `c937604` | Merge 17 upstream commits, incl. **Adams–Bashforth time integration** |
| `70539c5` | Name the VOF pseudo-iteration threshold in `extend.f90` |
| `183391a` | Make particle physical and collision parameters runtime inputs (`&particle`, `&collision_parameters`) |
| `9073074` | Reorganise `examples/` by phase count; add `Sinking_Sphere` |
| `7fd5ccb` | Make initial particle velocity a runtime input (`u_ini,v_ini,w_ini`) |
| `f712bd1` | Corrected velocity of the sinking sphere |
| `3c63107` | Add `Head_On`, `Particle_Capture`, `Sessile_Drop` example cases |
| `fdcd83b` | **Fix lubrication coefficient and duplicate torque assignment typos** |
| `7e46c21` | Allow the initial position of a single particle (`np=1`) from `input.nml` (`x_ini,y_ini,z_ini`) |
| `5044311` | Explicit initial position for `Head_On` |
| `344a400` | Explicit initial position for `Sessile_Drop` |
| `26b6dfb` | Add `Wall_Collision` example case |

Direction of travel is clear: **move hard-coded constants into `input.nml`**, and
**add validation cases**. Follow that pattern for new work.

---

## New files (not in upstream at all)

### Particle module — `src/prt_*.f90` (~7 000 lines)
`prt_param`, `prt_common`, `prt_initparticles`, `prt_digitiser`, `prt_initeul`,
`prt_initvof`, `prt_eulint`, `prt_intgr_over_sphere`, `prt_intgr_nwtn_eulr`,
`prt_collisions`, `prt_output`, `prt_loadpart`, `prt_phase_indicator`,
`prt_coordsfp`, `prt_interp_spread`, `prt_forcing`, `prt_kernel`, `prt_out.h90`

### Contact-line model
- `src/extend.f90` (135 lines) — extension velocity + upwind relaxation
- `src/rotnorm.f90` (147 lines) — capillary force integration

### Post-processing
- `PostPrt/PrtPos/` — particle trajectory reader
- `PostPrt/SphereData/` — surface force decomposition & interpolation
- `src/poslfp/` — Lagrangian forcing-point post-processor

### Examples — `examples/Three_Phase/`
`Bouncing_Sphere`, `Sinking_Sphere`, `Head_On`, `Particle_Capture`,
`Sessile_Drop`, `Wall_Collision`

---

## Modifications to upstream files

| file | change |
|---|---|
| `main.f90` | +433 lines. Particle init/coupling blocks, the contact-line pseudo-loop (twice), Adams–Bashforth driver, dimensionless-number banner, `forces_data.csv`, `Fstot` plumbing |
| `param.f90` | +`theta` in `&two_fluid`; +`nh_wide`, `is_ibm` under `_PARTICLE` |
| `types.f90` | +9 lines of kind aliases |
| `bound.f90` | +77 lines — BC handling extended for the particle fields |
| `initmpi.f90` | +36 lines — extra neighbour/topology info for `prt_comm_cart` |
| `rk.f90` | +`vofmin` and a commented-out Gauss-quadrature THINC variant |
| `two_fluid.f90` | small changes to `init2fl` |
| `vof_thinc_qq.f90`, `acdi.f90` | minor |
| `out2d.h90`, `out3d.h90` | `Alpha_C` output added; `Kappa` 3-D output commented out |
| `build.conf` | `PARTICLE`, `EULER` switches; VOF+SDF defaults |
| `configs/flags.mk` | `-D_PARTICLE`, `-D_EULER` mapping |
| `examples/` | upstream cases moved under `Two_Phase/` |

---

## Merging future upstream changes

`git remote` already has `origin` pointing at upstream. The last sync was
`c937604`. To pull new upstream work:

```bash
git fetch origin
git merge origin/main          # expect conflicts concentrated in main.f90
```

**Expect conflicts almost exclusively in `src/main.f90`** — it is the only file
with substantial interleaved changes from both sides. The `prt_*` files,
`extend.f90` and `rotnorm.f90` are fork-only and will never conflict.

Conflict-prone regions of `main.f90`:
- the `use` block (fork adds ~20 `use prt_mod_*` lines)
- the allocation block
- the initialisation sequence around `init2fl`
- the timeloop body, especially around `tm` and the projection

When resolving, the invariants to preserve are:
1. `initeul` must run **before** the contact-line pseudo-loop (it supplies `alphac`
   and `norm_part*`).
2. `eulint` must run **after** the momentum predictor `tm` and **before**
   `fillps`/`solver` (it modifies `u,v,w` which then get projected).
3. `intgr_nwtn_eulr` must run **after** the projection (it needs the corrected
   velocity for the fluid-momentum integrals).
4. `rot_norm` must run **after** the pseudo-loop (it integrates the relaxed `psi`).

---

## Known latent issues in fork-specific code

Recorded here so they are not re-derived each session. **None of these have been
fixed** — check the file before assuming otherwise.

1. **`extend.f90:110-122` — gradients divided by `dli` instead of multiplied.**
   `dli` is the inverse spacing, so the upwind derivative is off by `dx²`.
   Partially compensated by `dtau = 0.3*minval(dli)` also being inverse-scaled.
   See `contact-line-model.md`. Fixing one without the other will change results.

2. **Capillary force is computed but not applied.** `Fstot` → `F_cap` is logged
   to `forces_data.csv`; the momentum terms are commented out at
   `prt_intgr_nwtn_eulr.f90:675,685,694,725,735,744`.

3. **Band mismatch** between `extend.f90` (`alphac > 0.5`) and `rotnorm.f90`
   (`alphac > 0`).

4. **`rkcoeffab` is identically 1** under Adams–Bashforth, so every
   `*rkcoeffab` in `prt_intgr_nwtn_eulr.f90` is a no-op left from RK3.

5. **`forces_data.csv` is opened before `MPI_INIT`** (`main.f90:181`) with
   `status='replace'` — every rank truncates the same file.

6. **`ratiorho` and `rho_s` are independent inputs** that both describe the
   particle density. `rho_s` drives the Newton–Euler balance; `ratiorho` drives
   the collision effective masses. Nothing checks they agree.

7. **Stale `!$omp` clauses** in `prt_initvof.f90` and others list variables that
   no longer exist. Harmless while `OPENMP=0`; would fail to compile at `=1`.
