---
name: debug-blowup
description: Diagnose a simulation that crashed, diverged, produced NaN, aborted on divergence or timestep, or shows unphysical fields. Use when the user reports a run that stopped, blew up, gives weird values, or asks "why is X appearing in my simulation".
---

# Diagnosing a failed or suspicious run

## First: read `build.conf`

Before reasoning about any interface or particle behaviour, confirm which code
paths are actually compiled.

```bash
grep -E "INTERFACE_CAPTURING_VOF|SDF_NORMALS|PARTICLE|EULER|CONSTANT_COEFFS|SINGLE_PRECISION" build.conf
```

With `INTERFACE_CAPTURING_VOF=1` (the current setting) **`src/acdi.f90` is dead
code** — the interface is handled by `src/vof_thinc_qq.f90`. With `EULER=1` the
Lagrangian particle files are not compiled. Diagnosing from the wrong file has
already produced one wrong answer in this project's history.

## The three abort paths in `main.f90`

All fire inside the `mod(istep,icheck) == 0` block (`main.f90:809-841`).

### 1. `ERROR: maximum divergence is too large`
```
Total divergence = ...  | Maximum divergence = ...
```
The pressure projection failed to make the velocity solenoidal, or a NaN has
appeared. Usual causes, in order of likelihood:

- **`dt` too large for surface tension.** Compare the running `dt` against
  `dt_cap` printed in the startup banner. If `dt > dt_cap` the capillary wave
  is unresolved. Lower `cfl`, or set `dt_f` below `dt_cap`.
- **Collision instability.** `dt_estim` in `&collision_parameters` far from the
  actual `dt` makes the DEM springs numerically rigid. Read the reported `dt`
  and set `dt_estim` to match.
- **Density/viscosity ratio too extreme** for the grid. `rho12 = 935., 1.` is
  already ~1000:1.
- **BC inconsistency** that `test_sanity_input` did not catch — check `cbcvel`
  against `cbcpre` (wall = `'D'`/`'N'`, periodic = `'P'`/`'P'`).

To *watch* the blow-up develop instead of aborting, set
`MASK_DIVERGENCE_CHECK=1` in `build.conf` (then `make clean && make`) and write
frequent 3-D output. Never leave it on for production.

### 2. `ERROR: time step is too small`
`dt_cfl < small`. Something has already diverged — the velocity field contains
huge or NaN values. Treat as case 1; the divergence check usually fires first.

### 3. `WARNING: fixed time step exceeds estimated stability limit`
Only when `dt_f > 0`. The run **aborts**. Either lower `dt_f` or switch to
adaptive (`dt_f = -1`).

## Startup aborts

| message | cause | fix |
|---|---|---|
| `Radius spheres larger than x-dimension of processes` | `radius + offset > l(1)/dims(1)` | **use fewer MPI ranks**, or enlarge the domain. Do not edit the check in `prt_common.f90`. |
| `Error reading the input file` | `input.nml` missing or malformed in the run directory | check you are running from `run/` and the file is there |
| sanity-check failures | BC/decomposition inconsistency | read the message from `sanity.f90` |

## "Another phase is appearing where there should be none"

**Already diagnosed — do not re-derive.** If the spurious volume fraction is
`~1e-16`, it is machine-epsilon round-off, not physics:

- Origin: the contact-line relaxation loop (`src/extend.f90`) runs
  **unconditionally every timestep** in every cell with `0 < alphac < 1`,
  whenever `is_track_interface = T`. It performs several vector normalisations
  per cell per iteration, 5 iterations per step, each guarded only by
  `+ epsilon(1._rp)`.
- It appears as faint concentric "levels" hugging the particle, following the
  `alphac` shells.
- `1e-16` is ~15 orders below any meaningful volume fraction. Density blending
  `rho12(1)*psi + rho12(2)*(1-psi)` is unaffected.
- **Fix only if it pollutes a diagnostic:** clip after the pseudo-loop with
  `where(abs(psi) < 1e-12) psi = 0._rp`.
- If it is *growing* or spreading into the bulk, that is different — then check
  `dtau`/`max_pseudo_iter` (`main.f90:637-638`) and the THINC `beta`.
- **This would stop being negligible under `SINGLE_PRECISION=1`** (`~1e-7`).

## Particle behaving wrongly

| symptom | check |
|---|---|
| particle passes through a wall | `dt_estim` vs actual `dt`; `en`/`et`; `colthr_pw` |
| particle does not feel the interface | the capillary feedback is **commented out** in `prt_intgr_nwtn_eulr.f90:675,685,694,725,735,744` — `F_cap` is logged but not applied |
| particle drifts / wrong terminal velocity | `rho_s` vs `ratiorho` consistency; `intrho*` buoyancy integrals |
| contact angle not respected | `theta` is in **degrees**; check the sign conventions in `extend.f90` (`n_wall = -norm_part`, `pi - theta_rad`) |
| jumps at a restart | `Fstot`/`Fstot_old` are **not** checkpointed — one-step `F_cap` transient is expected |

## Instruments available

1. **`forces_data.csv`** in `run/` — written every `iout0d` steps:
   `F_drag,F_ibm,F_inertia,F_w,F_bouy,F_cap,ep_z,ep_w`.
   Plot `ep_z`/`ep_w` against time first; a force term exploding before the
   divergence abort localises the problem.

2. **stdout**. `Fstot` is printed every step from `main.f90:659`, and
   `F_sup/F_ibm/F_inertia/F_w/F_bouy/F_cap/ep_z/ep_w` from
   `prt_intgr_nwtn_eulr.f90:770-777`. Capture it:
   `mpirun -n N ./cans 2>&1 | tee run.log`

3. **`data/time.out`** — `istep, dt, time`. A collapsing `dt` is the earliest
   warning sign.

4. **`icheck`** — lower it to `1` to get divergence reported every step.

5. **Debug build** — `FFLAGS_DEBUG=1`, `FFLAGS_OPT=0`, `make clean && make`.
   Gives `-fcheck=all -ffpe-trap=invalid -finit-real=snan`, which turns silent
   NaN propagation into an immediate backtrace. Expect it to also trap on the
   guarded `0/0` cases in `extend.f90`/`rotnorm.f90` — that is the debug build
   being correct, not a new bug.

## Method

Work backwards from the last good output:

1. `data/time.out` → when did `dt` start collapsing?
2. `forces_data.csv` → which force term diverged first, and at what `ep_z`?
3. last written 3-D/2-D field → **where** in space did it start?
   (`alx_fld_*` = `alphac`, `psi_fld_*` = `psi`)
4. Only then read the routine that owns that region.

Do not start by reading code. Start by localising in time and space.
