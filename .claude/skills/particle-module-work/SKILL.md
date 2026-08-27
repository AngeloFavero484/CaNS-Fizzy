---
name: particle-module-work
description: Modify the particle module (prt_*.f90) or the contact-line model (extend.f90, rotnorm.f90). Use when adding a particle parameter, changing IBM coupling, collisions, the contact angle treatment, capillary force, or touching anything under _PARTICLE / _EULER.
---

# Working on the particle and contact-line code

Background reading before editing:
`.claude/references/data-structures.md` (types, MPI protocol),
`.claude/references/contact-line-model.md` (the physics),
`.claude/references/timeloop.md` (call order).

## Which files are live

Build is `PARTICLE=1`, `EULER=1`. That means:

**Live:** `prt_param`, `prt_common`, `prt_initparticles`, `prt_digitiser`,
`prt_initeul`, `prt_initvof`, `prt_eulint`, `prt_intgr_over_sphere`,
`prt_intgr_nwtn_eulr`, `prt_collisions`, `prt_output`, `prt_loadpart`,
`extend.f90`, `rotnorm.f90`.

**Not compiled** (`#if !defined(_EULER)`): `prt_interp_spread`, `prt_forcing`,
`prt_kernel`. The `EULER=0` block in `main.f90:715-732` is *additionally*
commented out. Do not revive that path unless asked.

## Adding a runtime parameter

This is the most common request, and there is an established pattern — see
commits `183391a`, `7fd5ccb`, `7e46c21`.

In `src/prt_param.f90`:

1. Declare it: `real(rp), protected :: my_param` (keep `protected`).
2. Add to the right `namelist` — `/particle/`, `/collision_parameters/`, or
   `/particle_euler/`.
3. **Set a default** in the "set-up default parameters" block. Every parameter
   has one; a missing namelist group must not break existing cases.
4. If it is derived from others, compute it in the "derived quantities" block
   after the `read` calls.
5. `use prt_mod_param, only: my_param` wherever it is consumed.
6. Add it to the relevant `examples/Three_Phase/*/input.nml` files.
7. Document it in `.claude/references/input-namelists.md`.

For a sentinel-style optional parameter (only meaningful in some configurations),
copy the `x_ini`/`pos_ini_unset` pattern: a `parameter` sentinel value, then a
conditional fallback in the derived-quantities block.

## Adding a field to `type particle`

**The field order in `type particle` is load-bearing.** `prt_common.f90` documents
the rule; `prt_param.f90` encodes the counts:

```fortran
integer, parameter :: send_real = 24 + 10*nqmax
integer, parameter :: send_int  = 1  + nqmax
```

- A field that must survive a **master change** goes inside the first contiguous
  real block (currently 24 reals, `x` through `coltz`) — and you must bump
  `send_real`.
- A field that does not need communicating goes **after** those blocks.
- Getting this wrong silently corrupts MPI transfers when a particle crosses a
  subdomain boundary. It will look like a physics bug.

Mirror the addition in `particle_old` if the previous-step value is needed.

## The master/slave boilerplate

Every routine that touches the sphere opens with ~150 lines of:
pack into `anb(0:8,1:npmax)` → `MPI_ISEND`/`MPI_IRECV` between master and slaves
→ `MPI_WAITALL` → periodic-image correction (`±l(1)`, `±l(2)`) → compute and
clamp `ilow..khigh` → loop.

It is duplicated verbatim in `prt_initeul`, `prt_initvof`,
`prt_intgr_over_sphere`, `prt_eulint`, `prt_phase_indicator`,
`prt_intgr_nwtn_eulr`. **A change to the protocol must be applied to all of
them.** Resist the urge to factor it out mid-task — that is a large, risky
refactor and the user has not asked for it.

Note `prt_comm_cart` is created **periodic in both split directions regardless of
the physical BCs**; non-periodic walls are handled separately via `is_bound`.

## Grid conventions that bite

- **`dli`, `dzci`, `dzfi` are INVERSE spacings.** A derivative is `* dli(1)`,
  never `/ dli(1)`. `extend.f90` got this wrong until 2026-08-25 — see below.
- `alphac` and friends are declared `dimension(1-nh_wide:,1-nh_wide:,1-nh_wide:)`
  while core fields are `(0:,0:,0:)`. `nh_wide = 1` so the bounds coincide, but
  do not mix the declaration styles across a call boundary.
- Staggering: digitise `alphac` **separately** at cell centres and at each face
  when it multiplies a velocity component. `prt_initvof.f90` shows the pattern
  (`alpha_eulx`, `alpha_euly`, `alpha_eulz`, plus `alpha_eulc`).

## Known latent issues — check before "fixing"

1. ~~`extend.f90` divides by `dli` instead of multiplying.~~ **Fixed 2026-08-25**,
   together with the paired `dtau` scaling in `main.f90`. *Still latent:*
   `advect_vof_upwind` uses the uniform `dli(3)` for z instead of `dzci`/`dzfi`,
   so a clustered z-grid (`gtype`/`gr` ≠ uniform) is mishandled there. All
   current cases are uniform, so this is inactive.

2. **The capillary force is computed but never applied.** `Fstot` from
   `rot_norm` is reduced, printed, and logged as `F_cap`, but the momentum terms
   are commented out at `prt_intgr_nwtn_eulr.f90:675,685,694,725,735,744`.
   Uncommenting those six lines is the switch to enable feedback.

3. **`rkcoeffab` is identically 1** under Adams–Bashforth. Every `*rkcoeffab`
   in `prt_intgr_nwtn_eulr.f90` is a no-op left over from RK3.

4. **Band mismatch:** `extend.f90` relaxes on `alphac > 0.5`, `rotnorm.f90`
   integrates over `alphac > 0`.

5. **`rho_s` vs `ratiorho`** are independent inputs both describing particle
   density — `rho_s` drives Newton–Euler, `ratiorho` drives collision effective
   masses. Nothing checks consistency.

6. **Stale `!$omp` clauses** list variables that no longer exist (e.g. `cas`,
   `coorxmin` in `prt_initvof.f90`). Inert at `OPENMP=0`; would fail to compile
   at `=1`.

## Preserving invariants in `main.f90`

If you reorder anything in the timeloop:

1. `initeul` before the contact-line loop — it supplies `alphac` and `norm_part*`.
2. `eulint` after the momentum predictor `tm`, before `fillps`/`solver`.
3. `intgr_nwtn_eulr` after the pressure projection.
4. `rot_norm` after the pseudo-loop.

Also: the contact-line block exists **twice** — once for the initial condition
(`main.f90:516-532`) and once in the loop (`652-671`). `dtau_cfl` and
`max_pseudo_iter` now come from `&contact_line` in `input.nml`, so retuning them
needs no code edit -- but changes to the **call sequence** must still be made in
**both** copies.

## OpenACC

Keep `!$acc` directives consistent even though `GPU=0`. Match the surrounding
idiom: `async(1)` everywhere, explicit `!$acc wait` before host access or MPI,
`update self` / `update device` around host-side work. A mismatched
`enter data`/`exit data` pair is a latent bug a CPU build will not catch.

## After editing

```bash
make clean && make          # required if you touched build.conf or a .h90
cd run && mpirun -n 1 ./cans
```

Verify against a known case before trusting a change — `Bouncing_Sphere` and
`Sinking_Sphere` differ only in surface tension and should bounce vs. penetrate
respectively. That contrast is the cheapest regression test available.
