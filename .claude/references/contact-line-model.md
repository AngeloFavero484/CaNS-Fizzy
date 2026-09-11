# Reference: the extended contact-line model

This is the physics that makes this fork different from upstream CaNS-Fizzy.
Two files, ~280 lines total, called directly from `main.f90`.

---

## The problem being solved

A fluid–fluid interface (`psi`) meets a rigid sphere (`alphac`) along a curve —
the **contact line**. Physically the interface must meet the solid at a
prescribed **contact angle** `theta`. Nothing in the Navier–Stokes / VOF
discretisation enforces that, so it is imposed explicitly, each timestep, by
relaxing `psi` in a narrow band around the sphere.

Three fields participate:

| field | source | meaning |
|---|---|---|
| `psi` | `vof_thinc_qq.f90` | fluid-1 volume fraction, `[0,1]` |
| `alphac` | `prt_digitiser.f90` via `prt_initeul.f90` | solid indicator, `[0,1]`, diffuse over `eps_sol` cells |
| `norm_partx/y/z` | `prt_initeul.f90` | **exact** analytic outward normal of the sphere at each cell |

The contact-line band is defined as the cells where `alphac` is strictly between
0 and 1 — i.e. the diffuse solid shell.

---

## Step 1 — the extension velocity (`src/extend.f90::compute_uextend`)

Runs on cells with `alphac > alpha_min .and. alphac < 1` (`alpha_min` defaults
to **0.5**, a narrower band than `rot_norm` uses; `&contact_line` in
`input.nml`).

Given the solid normal `n_wall = -norm_part` (pointing *into* the solid) and the
interface normal `n_int = (normx,normy,normz)`:

```fortran
n1 = -normalize( n_int × n_wall )     ! tangent to the contact line
n2 = -normalize( n1  × n_wall )       ! tangent to the solid, ⟂ to the contact line
c  = n_int · n2                       ! which side of the solid the interface leans to
theta_rad = theta * pi / 180
cot_theta = cos(pi - theta_rad) / sin(pi - theta_rad)

if (|c| < eps)     u_ext = n_wall
else if (c < 0)    u_ext = n_wall - cot_theta * n2
else               u_ext = n_wall + cot_theta * n2

u_ext = normalize(u_ext)              ! unit vector
```

`u_ext` is the direction in which `psi` must be transported so that the
interface, once relaxed, meets the wall at `theta`. It is a **unit vector by
construction**, so its magnitude carries no information — only `dtau` and the
iteration count set how far the relaxation goes.

Sign conventions to be careful with: `n_wall` is the **negated** particle normal,
both `n1` and `n2` are negated after normalisation, and the angle used is
`pi - theta_rad`, not `theta_rad`. Changing any one of these flips wetting to
non-wetting.

## Step 2 — the relaxation (`src/extend.f90::advect_vof_upwind`)

First-order upwind advection of `psi` along `u_ext`, on the same band as
step 1, `alphac > alpha_min .and. alphac < 1`:

```fortran
psi(i,j,k) = psi(i,j,k) - dtau * (u*dpsidx + v*dpsidy + w*dpsidz)
```

Driven from `main.f90` as `max_pseudo_iter` iterations with
`dtau = dtau_cfl / maxval(dli(1:3))` (i.e. `dtau_cfl` is a CFL number on the
smallest cell). Both come from `&contact_line` in `input.nml` and default to
`5` and `0.3`, the values they were hard-coded to before promotion. The loop
appears twice in `main.f90`: once in the timeloop and once for the IC.

> **Fixed (2026-08).** This routine previously divided by `dli` instead of
> multiplying — `dli` is the *inverse* spacing, so the upwind derivative was
> computing `dpsi * dx` instead of `dpsi / dx`. It was paired with
> `dtau = 0.3 * minval(dli)`, which is also inverse-scaled, and on an
> **isotropic uniform grid the two errors cancel exactly**:
> `(0.3/h)·u·(Δpsi·h) == (0.3·h)·u·(Δpsi/h)`. Every example case uses such a
> grid, which is why the model behaved correctly and the defect stayed hidden.
> Both halves were corrected together, so results on uniform grids are
> unchanged (to within floating-point ordering); on an **anisotropic or
> z-stretched grid the old code was wrong** and the new code is right.
>
> Remaining limitation: `advect_vof_upwind` still uses the uniform `dli(3)`
> for z rather than the local `dzci`/`dzfi`, so a clustered z-grid
> (`gtype`/`gr` ≠ uniform) is still not handled correctly here. All current
> cases use `gtype = 1, gr = 0.`, so this is latent, not active. Fixing it
> requires passing `dzci`/`dzfi` into the routine.

## Step 3 — the capillary force (`src/rotnorm.f90::rot_norm`)

Runs on the **wider** band `alphac > 0 .and. alphac < 1`.

The contact line is located as the intersection of the two interfaces, and its
local direction is the cross product of the two gradients:

```fortran
∇psi    = centred difference of psi            (correct: * dli)
∇alphac = |∇alphac| * norm_part                 ! magnitude from the field,
                                                ! DIRECTION from the exact sphere normal
vect = ∇psi × ∇alphac
prod = |vect|                                   ! ~ contact-line length density
```

Using the analytic `norm_part` for the direction instead of the digitised
gradient is deliberate — it removes the staircase noise of a Cartesian-digitised
sphere from the force integral.

The force direction is the tangent that pulls along the interface, away from the
solid:

```fortran
t = normalize( n_int × vect )
Fs = Fs - sigma * prod * t * dV
```

then `MPI_ALLREDUCE(Fs → Fstot, SUM)` in `main.f90:657`.

`perimetro_numerico = Σ prod·dV` is accumulated as a diagnostic (the numerical
contact-line perimeter) but only printed from commented-out lines at the end of
the routine — useful to re-enable when validating.

### Where `Fstot` goes

- `main.f90` prints it every step: `PRINT *, "Fstot", Fstot`
- passed into `intgr_nwtn_eulr`, logged to `forces_data.csv` as `F_cap`
- **NOT applied to the particle** — see below, this is deliberate.

#### Why it is disabled: the capillary force is already in `F_ibm`

The call order decides it. `main.f90:703` runs the momentum step, which adds the
CSF term `sigma*kappa*grad psi` (`mom.f90::momz_sigma`). Only then, at
`main.f90:750`, does `eulint` compute the IBM reaction

```fortran
fz = alpha_eulz*rhoz*(wl - wnew(i,j,k))*dti          ! prt_eulint.f90:462
```

`wnew` there has *already* felt the surface-tension force. So `fzltot` — and
therefore the `F_ibm` column — is the force needed to restore rigid-body motion
against a velocity field that already carries the capillary contribution.
Adding `Fstot` on top of it double-counts.

This is confirmed by the shape of the commented-out block itself. It is **not**
six lines that switch the capillary force on; it is **twelve, in matched pairs**:

| lines | term |
|---|---|
| `679, 689, 698, 724, 734, 743` | `-(fcap+fcap_old)/2 · dtp/(vol·rho_s)` |
| `681, 691, 700, 726, 736, 745` | `+(Fstot+Fstot_old)/2 · dtp/(vol·rho_s)` |

A **substitution**: remove the CSF capillary force the IBM absorbed, then put
`rotnorm.f90`'s explicit contact-line integral in its place. Enabling only the
`+Fstot` half is a bug, not a switch.

`ep(p)%fcapz` (`prt_eulint.f90:492`) is the `alpha_eul`-weighted `sigma*kappa*grad psi`
— a direct measurement of how much capillary force the IBM picked up. Both
estimates sit side by side in `forces_data.csv`:

- `F_cap_ibm` = `-½(fcapz+fcapz_old)·rkcoeffab`
- `F_cap`     = `+½(Fstot+Fstot_old)·rkcoeffab`

If the two routes measure the same physics then `F_cap_ibm ~ -F_cap`, i.e.
**their sum is ~0 and the substitution is a no-op**. How far from zero it runs is
a direct measure of how much re-enabling the block would actually change the
dynamics. They are different discretisations though — `fcapz` over the
`alpha_eul` band truncated at `2*eps_sol`, `Fstot` over `alphac > 0` (the band
mismatch noted above) — so expect agreement in trend and magnitude, not to
round-off.

Empirically the current arrangement is the validated one: Bouncing_Sphere
(`sigma = 79112.9`) and Sinking_Sphere (`sigma = 96316.4`) both reproduce their
reference results with `Fstot` disabled and surface tension fully active.

So the contact-line model **shapes the interface** (via `extend.f90`) and
**measures** the capillary force independently, while the force actually felt by
the particle arrives through the IBM reaction.

---

## Dead code kept in `rotnorm.f90`

Lines ~47–78 hold an earlier formulation, commented out: it rotated `normx/y/z`
directly at the contact line,
```
n_new = -cos(theta)·n_wall + sin(theta)·n_tangential
```
and then recomputed `kappa` as `-div(n)`. This is the more standard
"normal rotation" contact-angle approach (and explains the module name
`mod_rotnorm` / `rot_norm`, which no longer describes what the routine does).
It was superseded by the extension-velocity approach in `extend.f90`.
Keep it — it is the documented alternative if the current model misbehaves.

---

## The extension as a boundary condition (`psi_cl`, 2026-09-10)

`main.f90` allocates `psi_cl` alongside `psi`. Three things define it:

- **Seeded once**, in the startup block before the time loop, from
  `psi_cl = psi` followed by
  `n_seed = max(max_pseudo_iter, ceiling(3*eps_sol/dtau_cfl))` pseudo-iterations
  — enough pseudo-time to cross the whole band from scratch (the extension
  advances `dtau_cfl` cells per iteration; the band is ~`eps_sol` cells wide).
  This also covers restarts: `psi` comes from the checkpoint and `psi_cl` is
  rebuilt from it, so `psi_cl` is **not** in the checkpoint.
- **Persistent.** Each step only the part *outside* the band is re-seeded,
  `where alphac <= alpha_min: psi_cl = psi`. Inside the band the previous
  step's extension is kept, so `max_pseudo_iter` iterations per step only have
  to follow the interface's motion — same cost and same convergence as before.
- **Host-only**, like the rest of the contact-line path: `mod_extend` and
  `mod_rotnorm` carry no `!$acc` and `alphac` is in no device data region, so
  `psi_cl` deliberately carries none either and is not in an `enter data`.

`advect_vof_upwind`, `cmpt_norm_curv` inside the loop, and `rot_norm` all take
`psi_cl`. `psi` is untouched between the two `cmpt_massbal` probes.

### What changed physically

`rho`, `mu` and the transported `psi` are now the *un-extended* field — the
extension no longer bleeds into the fluid properties. The angle is enforced
*dynamically* (unbalanced `sigma*kappa*grad psi` drives the interface until the
extended field is a smooth continuation, which happens at theta) rather than
*kinematically* (psi dragged to the angle). The fixed point should be the same;
**the equilibrium apparent angle has not yet been re-measured** — that is the
outstanding validation, a theta sweep on Sessile_Drop out to `t ~ 10`.

### What was verified (2026-09-10, local, Sessile_Drop 64x64x48, sigma = 1000)

- `V_out_ext - V_out_adv = 0` and `V_tot_ext - V_tot_adv = 0` **bit-for-bit**,
  every step. The relaxation is exactly volume-preserving by construction.
- theta = 150, 30 steps: whole-step `dV_out` `-0.146 %` -> `-0.064 %`,
  `dV_tot` `-0.641 %` -> `-0.063 %`, `V_in` frozen (`17.82 -> 17.25` becomes
  `19.75 -> 19.74`).
- theta = 30, matched step count: the relaxation's injection into `V_out` goes
  from `+1.5e-2` to exactly `0`. The baseline aborted on divergence at step 30
  in this example configuration; the new version ran 300 steps. Not claimed as
  a stability fix — the flow differs.
- **A residual `V_out` drift remains** (`-0.98 %` by `t = 0.26` at theta = 30,
  in the violent initial transient). It is *not* the extension: it is THINC
  transport plus `clip_field` losing psi into the solid. That is now the next
  thing to measure, on a settled drop rather than a startup transient.

The `V_out_adv` / `V_out_ext` column pair in `mass_data.csv` is kept as the
regression check on the bit-for-bit invariance. The whole-step drift is the
difference of `V_out_adv` between consecutive rows.

---

## Known numerical behaviours

1. **`psi ≈ 1e-16` noise around the particle.** Machine epsilon injected by the
   repeated normalisations (`/ (norm + epsilon(1._rp))`) across the
   `max_pseudo_iter` (default 5) iterations per step. Appears as faint concentric "levels" following the `alphac` shells.
   Physically meaningless — 15 orders below any real volume fraction. Already
   diagnosed and dismissed with the user. Clip with
   `where(abs(psi) < 1e-12) psi = 0._rp` only if it pollutes a diagnostic.
   *Would stop being negligible under `SINGLE_PRECISION=1`* (`~1e-7`).

2. **Band mismatch — open.** `rotnorm.f90` uses a hard-coded `alphac > 0`
   while `extend.f90` uses `alphac > alpha_min`, so the capillary force is
   integrated over a wider shell than the one the interface is relaxed on.
   Switching `rotnorm.f90` to `alpha_min` was committed in `f45d41a` and
   reverted on 2026-09-10 with the rest of that day's code. **Builds from
   `f45d41a` through `5eabb52` carry the `alpha_min` band**, so their `F_cap`
   column in `forces_data.csv` is not comparable with any other build's.

3. **Normals inconsistency — worse than it reads.** The main phase-field step
   computes normals and curvature from `phi` (the SDF) under `_SDF_NORMALS`,
   and then the pseudo-loop recomputes them from `psi` directly, with no guard.
   `cmpt_norm_curv_youngs` (`two_fluid.f90:81`) loops the **whole domain**, so
   this is not a local effect: whenever the contact-line relaxation runs,
   `_SDF_NORMALS` is discarded everywhere and every `kappa` in the momentum
   equation is a Youngs difference of the raw VOF field. `build.conf` sets
   `SDF_NORMALS=1`; the loop throws it away. Leading suspect for the
   nucleation — see `planned-changes.md` item 2.

4. **No `psi` clipping inside the loop** — but it does not bite at the default
   settings. `clip_field` (in `two_fluid.f90`) is only applied at the end of
   `rk_2fl` (`rk.f90:252`), so anything `advect_vof_upwind` writes survives
   until the *next* step's advection. Measured over 12 Sessile_Drop runs
   (2026-09-09, see the mass-budget section below), the out-of-range content of
   `psi` after the relaxation loop is **exactly zero** in every step of every
   run: first-order upwind at `dtau_cfl = 0.3` is monotone, so it cannot create
   a new extremum. This would stop holding if `dtau_cfl` were pushed past 1.

5. **The relaxation is a volume source, and it pits the cells next to its own
   band.** The two sections below are the study of both. **Neither is fixed in
   `HEAD`** — see the reverted-attempts sections below and
   `planned-changes.md`.

---

## Mass conservation: what the relaxation does to the fluid volume

Measured 2026-09-09 with `mod_massbal` (`src/massbal.f90`), which writes
`mass_data.csv` — the fluid-1 volume probed twice per step, after `tm_2fl` and
after the relaxation loop, split by the solid indicator:

```
V_tot = sum psi *dV     V_out = sum psi*(1-alphac) *dV     V_in = sum psi*alphac *dV
```

Sessile_Drop at 16 points per diameter (`ng = 64,64,48`, `l = 16,16,12`,
`radius = 2`), particle held fixed (`is_solve_nwtn_eulr = F`), `gacc = 0`,
`alpha_min = 0.5`, over `0 < t < 10`:

| theta | dV_tot | dV_out | V_in(T)/V_particle |
|---|---|---|---|
| 30 | **+14.04 %** | +0.348 % | 14.4 % |
| 60 | +10.19 % | +0.180 % | 10.4 % |
| 90 | +4.70 % | +0.042 % | 4.9 % |
| 120 | +1.50 % | -0.022 % | 1.7 % |
| 150 | +0.14 % | -0.061 % | 0.30 % |
| 150, `max_pseudo_iter = 0` | -0.018 % | -0.043 % | 0.08 % |

**Read `V_out`, not `V_tot`.** The raw total is not a mass-conservation
diagnostic in this fork — up to 14 % of it is fluid buried inside the particle,
where it is not physical. Subtracting it recovers a quantity conserved to a few
tenths of a percent.

The last row is the control: with the relaxation switched off nothing is
injected, so `extend.f90` is unambiguously the source. What it injects lands
almost entirely inside the solid (at theta = 30: 4.836 injected, `V_in(T)` =
4.809).

### Mechanism

Not overshoot — `psi` stays exactly in `[0,1]` (see behaviour 4 above). The
source is that `advect_vof_upwind` uses the **advective** form
`psi -= dtau*(u.grad psi)` on a **masked** band: the upwind stencil reads
neighbours outside `alphac > alpha_min` that are never debited in return, so
nothing telescopes and the update has no discrete conservation property. The
strong theta-dependence follows the wetted fraction of the band — `cot(pi-theta)`
swings `u_ext` tangential, and a wetting drop covers far more of the shell.

### Bounded vs unbounded

Out to `t = 40` (theta = 30 and 90):

- **`V_in` saturates** — 4.81 -> 5.14 (theta = 30), 1.64 -> 1.92 (theta = 90).
  Once the diffuse shell has filled it is a one-off offset, not a runaway.
- **`V_out` does not.** It drifts linearly at `+0.0099 %/t` (theta = 30) and
  `+0.0098 %/t` (theta = 90) — essentially theta-independent asymptotically,
  against a `-0.0011 %/t` baseline with the relaxation off. This is the only
  genuinely unbounded error, and at `t = 100` it is ~1 % of the drop.

### `alpha_min` is not the knob *for the drift*

At theta = 150, sweeping `alpha_min` over 0.1 -> 0.9 moves the residual drift
only from `-0.011` to `-0.003 %/t` and `V_in` from 0.096 to 0.090. It sets the
band width; **theta** sets the injection. (It *is* the knob for the near-wall
voids — see the next section. The two statements are consistent: `alpha_min`
barely moves the integrated drift while strongly moving where the error is
deposited.)

All of the above is at 16 points per diameter and sigma = 1000. For the
resolution scaling, and for why `alpha_min` *is* the knob for the near-wall
voids even though it barely touches the drift, see the next section.

---

## Near-wall nucleation: the same defect, seen locally

Measured 2026-09-09; runs and scripts in `studies/contact-line-2026-09-09/`
(gitignored). Spurious depressions in `psi` appear in the fluid on the
drop–sphere surface — reported as "nucleation", clearest at high resolution.

They are **not a separate problem**. They are the *donor cells*: the fluid cells
immediately outside `alphac > alpha_min`, which `advect_vof_upwind`'s upwind
stencil reads but never updates. The mass budget sees the net transfer as a
volume source; the pitting is the same transfer seen locally.

The pairing is explicit in the raw field — at a void, `psi = 0.88-0.90` with
`alphac = 0.00-0.31` (outside the band), while the band cells directly beneath
at `alphac = 0.60-1.00` sit full at `psi ~ 1.00`.

Sessile_Drop, theta = 30, sigma = 100, particle fixed, `t ~ 12`. "voids" counts
majority-fluid cells (`alphac < 0.5`) with `psi < 0.95` that still have bulk drop
(`psi > 0.98`) further out — the qualifier matters, or the drop's own free
surface and the solid interior both register as voids:

| D/delta | `alpha_min` | voids | worst dip | dV_out | dV_tot | V_in |
|---|---|---|---|---|---|---|
| 16 | 0.5 | 188 | 0.177 | +0.256 % | +6.64 % | 2.34 |
| 32 | 0.5 | 82 | 0.256 | +0.291 % | +2.70 % | 0.91 |
| 16 | 0.1 | **0** | — | +2.77 % | +13.65 % | 4.18 |
| 32 | 0.1 | **338** | 0.888 | +1.26 % | +5.54 % | 1.61 |

### `alpha_min` moves both, in opposite directions

Narrow the band and the injection falls but the donor cells pit; widen it and
the pitting goes but the solid floods. **`max_pseudo_iter` (1 vs 5) and
`dtau_cfl` (0.1 vs 0.3) change the void integral by under 5 %** — it is the
band-edge geometry, not the strength of the relaxation. This is the sharpest
evidence that the two artefacts are one defect: one knob, two symptoms,
opposite signs.

### Refinement concentrates the error, it does not remove it

At `alpha_min = 0.5`, going 16 -> 32 points makes the worst dip *deeper*
(0.177 -> 0.256) and `V_out` slightly *worse* (+0.256 -> +0.291 %). Only `V_in`
improves (2.34 -> 0.91), and only because the band is a fixed number of *cells*,
so its physical volume halves. Hence the phenomenon reads more clearly at
D/delta = 32: the same error packed into a thinner ring.

### Do not "fix" it with a low `alpha_min`

`alpha_min = 0.1` gives exactly zero voids at 16 points out to `t = 38.5`, and
holds at 32 points only to `t ~ 10` — then the drop bridges and traps a gas
pocket (`psi` down to 0.11, void integral 2.23) worse than the pitting it
replaced. Low `alpha_min` relocates the defect from pitting to bridging; which
one shows up depends on resolution and run length. There is no setting of
`alpha_min` that removes it.

### Answered, 2026-09-10: the depressions are capillary-mediated

The question was whether the relaxation writes them into `psi` directly, or
whether they are driven through spurious curvature (`cmpt_norm_curv` runs on
the relaxed field, so a non-physical profile feeds `sigma*kappa*grad psi`).

The `psi_cl` build (`f45d41a`) settles it. There the relaxation **cannot write
`psi` at all** — and the user reports the nucleation is unchanged, same as
before. So the depressions are not written; they are **advected in by the flow
that spurious `kappa` drives**.

That relocates the defect. It is not the advective form, not the masked update,
not the missing donor/receiver pairing — none of which can act on `psi` in that
build. What survives is that **`kappa` is garbage at the band edge**: the
relaxation writes inside `alphac > alpha_min` and not outside, so the relaxed
field has a kink at the band boundary, and `cmpt_norm_curv` differentiates
across it twice. The resulting capillary force pulls fluid into the band and
pits the cells just outside it — which is exactly where the voids are observed,
and exactly why `alpha_min` (which moves the band edge) is the knob that moves
them while `max_pseudo_iter` and `dtau_cfl` (which change the strength, not the
edge) do not.

**Consequence for the fix:** it is a smoothness problem at the band edge, not a
conservation problem. Ramping the relaxation strength continuously to zero as
`alphac -> alpha_min`, instead of the current hard on/off mask in
`compute_uextend`/`advect_vof_upwind` (`extend.f90:48,114`), would remove the
kink and with it the spurious curvature. Not attempted.

### The fix, for both

**Not** the flux form. `psi -= dtau*div(u_ext*psi)` with a paired band edge
conserves `sum_band psi`, which is the wrong integral: `u_ext` points into the
solid, so the relaxation moves `psi` toward higher `alphac` where the weight
`(1-alphac)` is smaller, and closing the band edge means the debit comes from
the fluid-side band cells. `V_out` would then *decrease* at comparable
magnitude — the sign flips, the drift does not go away. That earlier proposal
is superseded; do not implement it.

**The `psi_cl` boundary-condition reformulation — tried 2026-09-10, reverted,
not rejected in principle.** That version ran the relaxation on a persistent
copy so the transported `psi` was never written, making the extension exactly
volume-preserving (verified bit-for-bit). Commit `f45d41a`, reverted in
`HEAD`. Two results from it, both from the user:

1. **The nucleation was unchanged.** That is the diagnostic result above — it
   proves the pitting is capillary-mediated, not written. Worth the experiment
   on its own.
2. **It imposed the contact angle less well.** Expected, and the reason it is
   not in `HEAD`: with `psi` untouched, theta is enforced only *dynamically*
   through `kappa` in the momentum equation, never *kinematically*. The
   kinematic drag turns out to matter for how tightly the angle is held.

So `psi_cl` fixes volume but weakens the angle, and it does not touch the
nucleation. It is a live option to revisit — most plausibly combined with a
smoothed band edge (see above), which would fix the curvature that is the
actual cause of the pitting — but on its own it trades one problem for another.

**Nothing is in the code now.** The band-edge smoothing (`alpha_ramp`) and the
volume projection (`crrct_vout`) were both implemented and both reverted; each
is documented below with what it measured. The `extend.f90` / `rotnorm.f90`
hygiene fixes from `f45d41a` were reverted too, so the source is identical to
`db3dace`: the volume drift and the nucleation are both still present. Next steps are in `planned-changes.md` items 2 and 3.

---

## `alpha_ramp`: tried and reverted — and the rule it taught us

Implemented `dfc55b1`, reverted. It weighted the relaxation update by a quintic
smootherstep in `alphac`, ramping up from the outer band edge, on the theory
that the hard on/off mask left a kink for `cmpt_norm_curv` to differentiate.

**It did not change the nucleation at `alpha_ramp = 1`.** User-verified.

It could not have. The update is `psi -= w*dtau*(u_ext.grad psi)`, and any
positive scalar `w` multiplying the whole right-hand side has fixed point
`u_ext.grad psi = 0` — **the same fixed point as `w = 1`**. `psi` persists
across timesteps, so the relaxation is near-converged and `w` changes only the
convergence rate, never the converged field.

### The rule, which covers the whole 2026-09-09 knob sweep

| knob | what it is | effect on the converged field |
|---|---|---|
| `max_pseudo_iter` | rate | none |
| `dtau_cfl` | rate | none |
| `alpha_ramp` | rate | none — measured |
| `alpha_min` | geometry | the only one that moves the voids |

This is why 1 vs 5 iterations and 0.1 vs 0.3 CFL moved the void integral by
under 5 %, and why `alpha_min` was the only knob that ever did anything. The
sweep was measuring convergence rate, not the defect.

**Do not propose anything that only scales the relaxation.** To change the
converged field you have to change either the geometry of the region it acts on
or the *direction* of `u_ext` — at steady state `u_ext.grad psi = 0` means the
`psi` isosurfaces contain `u_ext`, so direction is what sets the interface
orientation. See `planned-changes.md` item 3.

### What is now the leading suspect

The relaxation loop discards the SDF-based curvature domain-wide and replaces
it with a Youngs difference of the raw VOF field, so `_SDF_NORMALS` is
effectively off whenever the contact line is active. Since `psi_cl` proved the
voids are driven *through* `kappa`, that is where to look next. Full write-up
and the proposed change in `planned-changes.md` item 2.

---

## `crrct_vout`: restoring V_out by projection — implemented, then reverted

Implemented in `3ab76ac`, reverted at the user's request when the work moved to
the band edge. **The code is not in `HEAD`**; this section records what it did
and what it measured, because the approach is still valid for the volume
problem and may come back.

It left the relaxation exactly as it was and removed the volume it injects
afterwards. `cmpt_massbal` already brackets the loop, so
`dvol = vol_ext(2) - vol_adv(2)` is the V_out just created, and it is taken out
by solving

```
sum dpsi*(1-alphac)*dV = -dvol      with   dpsi = -c*g,  g = psi*(1-psi)*(1-alphac)
```

— one scalar `c = dvol / sum g*(1-alphac)*dV`, one `MPI_ALLREDUCE`. The weight
`g` does three jobs: `psi*(1-psi)` confines the correction to interface cells,
`(1-alphac)` keeps it out of the solid interior and debits each cell in
proportion to what it contributes to `V_out`, and together they make the update
bound-preserving for any `|c| < 1` with no clipping.

Measured: per-step residual exactly zero once the seeding debt clears; at
theta = 150 the drift over the first 25 steps (`-0.034 -> -0.117 %`) became a
bounded wobble ~8x smaller. At theta = 30 the two were comparable — no
demonstrated win. It does **not** touch the pitting.

### Why the local numbers stop there — read before trusting any of this

The shipped `examples/Solid_Particles/Three_Phase/Sessile_Drop` at `sigma = 1000` on 64x64x48 is
**marginal at t ~ 0.09**: `dt_cfl` collapses to ~1e-8 and the run aborts on the
divergence check, in some configurations and not others, independently of any
switch. Past step ~25 every number above is contaminated, `dt` has collapsed so
`%/t` is meaningless, and nothing here can judge stability or long-time drift.

**Evaluation has to happen on the cluster with the study configuration** (the
2026-09-09 runs reached t = 40; those inputs are in the gitignored
`studies/contact-line-2026-09-09/`, not in the examples).

---

## Tuning knobs, in order of usefulness

| knob | where | effect |
|---|---|---|
| `theta` | `input.nml` `&two_fluid` | the prescribed contact angle, **degrees** |
| `eps_sol` | `input.nml` `&particle_euler` | width of the diffuse solid shell in cells → width of the contact-line band |
| `max_pseudo_iter` | `input.nml` `&contact_line` | default `5`. More = stronger enforcement, more round-off |
| `dtau_cfl` | `input.nml` `&contact_line` | default `0.3`; `dtau = dtau_cfl/maxval(dli)`, a CFL number on the smallest cell |
| `alpha_min` | `input.nml` `&contact_line` | default `0.5`, the relaxation band threshold. **Do not tune to chase the near-wall voids** — it trades them for flooding/bridging, see above |

All five are runtime inputs. The last three used to be hard-coded — in
`main.f90` (`max_pseudo_iter`, `dtau`) and `extend.f90` (`alpha_min`) — and were
promoted to the new `&contact_line` namelist following the same pattern as
commits `183391a` and `7fd5ccb`. Their defaults reproduce the old hard-coded
values exactly, so an `input.nml` without a `&contact_line` group is unchanged
bit for bit; see [`input-namelists.md`](input-namelists.md) for the full table.

Promoting `alpha_min` also folded in the *second* hard-coded `0.5` in
`extend.f90`, the band test at the top of `compute_uextend`. The two must move
together: building `u_ext` on one band and advecting `psi` on another is
incoherent. This does **not** touch `rotnorm.f90`'s `alphac > 0`, so the band
mismatch in the section above still stands — lowering `alpha_min` narrows it.
