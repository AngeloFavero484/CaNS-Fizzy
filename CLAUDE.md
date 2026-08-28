# CLAUDE.md — orientation for future sessions

This file is written for Claude. It captures the context that is *not*
derivable by skimming the code, so a new session can be productive immediately.

## What this repository actually is

This is **a personal research fork of [CaNS-Fizzy](https://github.com/CaNS-World/CaNS-Fizzy)**,
extended by Angelo Raimondo Favero (PhD work) into a **three-phase solver**:
two immiscible fluids **plus a finite-size rigid spherical particle**, with an
**extended contact-line model** that imposes a prescribed contact angle where the
fluid–fluid interface meets the particle surface.

Upstream Fizzy is two-phase only. Everything particle-related, the contact-line
model, and the Adams–Bashforth driver in `main.f90` are local additions.
See `.claude/references/upstream-divergence.md` for the exact delta.

The directory path itself encodes the research lineage:
`Multiphase/Three_Phase/Extended_CL/Adams_Bashfort/Bouncing_Sphere` —
three-phase, extended contact line, Adams–Bashforth time integration,
bouncing-sphere case. **Do not treat the path as throwaway**; sibling directories
elsewhere on the user's disk are other variants of the same code.

## Git remotes — important

| remote | URL | role |
|---|---|---|
| `origin` | `https://github.com/CaNS-World/CaNS-Fizzy` | **upstream**, read-only |
| `myfork` | `git@github.com:AngeloFavero484/CaNS-Fizzy.git` | **the user's fork — push here** |

`git push` with no remote goes to `origin`, which is upstream and wrong.
**Always `git push myfork main`.** The local `main` is many commits ahead of
`origin/main`; that is expected and permanent, not a state to "fix".

## Where the work happens

The user runs production simulations on **CINECA Galileo100** over SSH. The
cluster is split deliberately: the git clone and the build live in
`$WORK/CaNS-Fizzy` (`/g100_work/IscrC_TP-PBR`, 1 TB, survives the 40-day
`$SCRATCH` cleanup), and each run is staged as a self-contained directory under
`$SCRATCH/runs/<case>/` (no quota, right filesystem for parallel I/O) by
`$WORK/bin/stage-run`. See `.claude/skills/run-on-galileo/`.

`$WORK` is bound to the project grant and dies with it (`IscrC_TP-PBR` ends
2026-12-02), so it is staging, not an archive.

An older clone at `/g100_scratch/userexternal/afavero0/New_Version/CaNS-Fizzy`
is superseded but may still hold live output — check `squeue` before touching it.
The local machine is used for editing and analysis. Consequence:

- Changes must round-trip through `myfork`. Commit + push locally, pull on the cluster.
- A fresh clone on the cluster **needs `git submodule update --init --recursive`** —
  `dependencies/2decomp-fft` and `dependencies/cuDecomp` are submodules, and a plain
  clone leaves them empty (the symptom is `make: Nessuna regola per generare
  l'obiettivo "clean"` / `No rule to make target 'clean'`).
- HPC login nodes often block outbound port 22; GitHub over `ssh.github.com:443`
  is the workaround.

## Build configuration in use

`build.conf` is the single knob file. The current, meaningful settings:

```
FCOMP=GNU
INTERFACE_CAPTURING_VOF=1   # THINC/QQ volume-of-fluid, NOT ACDI
SDF_NORMALS=1               # normals/curvature from the SDF phi, not psi
CONSTANT_COEFFS_POISSON=1   # FFT-based constant-coefficient Poisson
PARTICLE=1                  # particle module on
EULER=1                     # Eulerian (not Lagrangian-marker) particle coupling
GPU=0
```

**`INTERFACE_CAPTURING_VOF=1` means `src/acdi.f90` is dead code in this
configuration.** Do not diagnose interface behaviour by reading `acdi.f90` —
the live path is `src/vof_thinc_qq.f90`. This has already caused one wrong
diagnosis in a past session; the user corrected it.

Likewise `EULER=1` means the `#if !defined(_EULER)` branches (Lagrangian forcing
points: `prt_interp_spread.f90`, `prt_forcing.f90`, `prt_kernel.f90`) are **not
compiled**. The live particle coupling is `prt_eulint.f90` + `prt_initeul.f90`.

## Things already established with the user (do not re-derive)

- **`psi ~ 1e-16` noise hugging the particle is not a bug.** It is machine-epsilon
  round-off injected by the repeated vector normalisations in the contact-line
  relaxation loop (`src/extend.f90`), which runs unconditionally every timestep in
  cells with `0 < alphac < 1`. It is ~15 orders below anything physical. Clip with
  `where(abs(psi) < 1e-12) psi = 0._rp` only if it pollutes a diagnostic.
- The user works in **Italian locale** — compiler/make errors come back in Italian.
- The user is a PhD researcher in multiphase CFD. Answer at that level: name the
  routine and line, state the mechanism, skip the tutorial.

## Conventions to respect when editing

- **Fortran 90 free form**, 2-space indent, `!` comment banners, lowercase keywords.
- Real kind is always `rp` from `mod_types`; write literals as `0.5_rp`, never `0.5d0`.
- Preprocessor guards are `#if defined(_PARTICLE)` etc. — leading underscore, set
  from `build.conf` via `configs/flags.mk`.
- Particle modules are named `prt_mod_*` in files `prt_*.f90`. Core solver modules
  are `mod_*` in files without the prefix.
- OpenACC directives (`!$acc`) are interleaved everywhere. Keep them consistent
  when editing loops even though the user currently builds with `GPU=0` —
  they must not rot.
- There is a lot of **commented-out code** (old Lagrangian IBM path, diagnostics,
  Runge–Kutta driver). It is deliberate history, not clutter. Do not delete it
  unless asked.

## Current state of the physics (as of 2026-08)

- Time integration is **Adams–Bashforth**, not RK3: `main.f90:600` is
  `do irk=1,1` with `tm_coeff = [2+dt/dto, -dt/dto]/2`. The RK3 loop above it is
  commented out. Anything assuming three substeps is stale.
- The capillary force on the particle, `Fstot`, is computed in `src/rotnorm.f90`
  by integrating `sigma * |grad psi x grad alphac| * t` over contact-line cells,
  then `MPI_ALLREDUCE`d.
- In `prt_intgr_nwtn_eulr.f90` the `Fstot` contribution to the particle momentum
  update is **currently commented out**, and that is deliberate, not an oversight:
  the capillary force is *already* inside `F_ibm`. `eulint` runs after the
  momentum step, so `fzltot = alpha_eulz*rhoz*(wl-wnew)*dti` is a reaction against
  a velocity that has already felt `momz_sigma`'s `sigma*kappa*grad psi`. The
  commented block is a **substitution** — twelve paired lines, `-fcapz` *and*
  `+Fstot` — not six lines that "switch the capillary force on". Enabling only the
  `+Fstot` half double-counts. See `.claude/references/contact-line-model.md`.
- `forces_data.csv` is opened by **rank 0 only**, after `MPI_INIT`, and written per
  `iout0d` steps from `prt_intgr_nwtn_eulr.f90`, where the rows are `MPI_REDUCE`d
  onto rank 0 (the particle's master rank changes as it crosses pencil
  boundaries). Columns:
  `Time,F_cap_ibm,F_ibm,F_inertia,F_w,F_bouy,F_cap,ep_z,ep_w`. Note `F_cap_ibm` is the
  CSF capillary force the IBM absorbed and `F_cap` is rotnorm's contact-line
  integral: they are two estimates of the same force, so `F_cap_ibm + F_cap ~ 0`
  is the check on whether the substitution above would change anything.
- A fresh (non-restart) run deletes any leftover `time.out`/`log_visu_*.out` in
  `datadir` right after `read_input` in `main.f90`, before either file is first
  written. Both are opened with `position='append'` by `out0d`/`write_log_output`,
  so without this cleanup a new run starting in a directory with old output would
  silently append after stale rows instead of starting clean. `forces_data.csv`
  is unconditionally reset with `status='replace'` regardless of `restart`, so it
  needed no equivalent change.

## Reading order for a cold start

1. `ARCHITECTURE.md` — repo structure and the timeloop.
2. `.claude/references/timeloop.md` — what happens each step, in call order.
3. `.claude/references/contact-line-model.md` — the physics that is unique to this fork.
4. `.claude/references/input-namelists.md` — every runtime knob.

Skills in `.claude/skills/` cover the recurring operational tasks
(building, adding a case, debugging a blown-up run, syncing to the cluster).
