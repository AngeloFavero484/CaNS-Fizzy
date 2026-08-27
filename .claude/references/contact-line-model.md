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

## Known numerical behaviours

1. **`psi ≈ 1e-16` noise around the particle.** Machine epsilon injected by the
   repeated normalisations (`/ (norm + epsilon(1._rp))`) across the
   `max_pseudo_iter` (default 5) iterations per step. Appears as faint concentric "levels" following the `alphac` shells.
   Physically meaningless — 15 orders below any real volume fraction. Already
   diagnosed and dismissed with the user. Clip with
   `where(abs(psi) < 1e-12) psi = 0._rp` only if it pollutes a diagnostic.
   *Would stop being negligible under `SINGLE_PRECISION=1`* (`~1e-7`).

2. **Band mismatch.** `extend.f90` uses `alphac > alpha_min` (default `0.5`),
   `rotnorm.f90` uses `alphac > 0`. The force is integrated over a wider shell
   than the one the interface is relaxed on. Lowering `alpha_min` narrows the
   gap; `rotnorm.f90`'s threshold is still hard-coded.

3. **Normals inconsistency.** The main phase-field step computes normals from
   `phi` (the SDF) under `_SDF_NORMALS`, but the pseudo-loop recomputes them
   from `psi` directly (`main.f90:649`). The contact-line normals are therefore
   noisier than the bulk ones.

4. **No `psi` clipping inside the loop.** `advect_vof_upwind` can push `psi`
   slightly outside `[0,1]`; `clip_field` (in `two_fluid.f90`) is only applied
   inside `rk_2fl`, before the loop runs.

---

## Tuning knobs, in order of usefulness

| knob | where | effect |
|---|---|---|
| `theta` | `input.nml` `&two_fluid` | the prescribed contact angle, **degrees** |
| `eps_sol` | `input.nml` `&particle_euler` | width of the diffuse solid shell in cells → width of the contact-line band |
| `max_pseudo_iter` | `input.nml` `&contact_line` | default `5`. More = stronger enforcement, more round-off |
| `dtau_cfl` | `input.nml` `&contact_line` | default `0.3`; `dtau = dtau_cfl/maxval(dli)`, a CFL number on the smallest cell |
| `alpha_min` | `input.nml` `&contact_line` | default `0.5`, the relaxation band threshold |

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
