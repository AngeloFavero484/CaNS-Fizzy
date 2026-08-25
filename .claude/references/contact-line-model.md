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

Runs on cells with `alphac > 0.5 .and. alphac < 1` (note: **0.5**, a narrower
band than `rot_norm` uses).

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

First-order upwind advection of `psi` along `u_ext`, on cells with
`alphac > 0.5 .and. alphac < 1`:

```fortran
psi(i,j,k) = psi(i,j,k) - dtau * (u*dpsidx + v*dpsidy + w*dpsidz)
```

> **Latent bug — the gradient is divided, not multiplied, by the spacing.**
> ```fortran
> dpsidx = (psi(i,j,k) - psi(i-1,j,k)) / dli(1)
> ```
> `dli` is the **inverse** spacing (`dli = 1/dl`, set in `param.f90`). A correct
> derivative is `(psi_i - psi_{i-1}) * dli(1)`. As written this computes
> `dpsi * dx` instead of `dpsi / dx` — off by `dx²`. Combined with
> `dtau = 0.3 * minval(dli)` (also inverted relative to a normal CFL scaling),
> the two errors partially compensate at the grid spacings currently in use,
> which is why the model behaves acceptably. **Do not "fix" one without the
> other**, and re-tune `dtau`/`max_pseudo_iter` if either is touched.
> Everywhere else in the codebase the idiom is `*dli(1)` — compare
> `rotnorm.f90:99-107`, which does it correctly.

Driven from `main.f90` as 5 iterations with
`dtau = 0.3_rp * minval(dli(1:3))`, both hard-coded
(`main.f90:637–638`, and again at `502–503` for the IC).

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
- **NOT applied to the particle** — the momentum terms
  `½dt(Fstot+Fstot_old)/(vol·rho_s)` are commented out at
  `prt_intgr_nwtn_eulr.f90:675, 685, 694, 725, 735, 744`.

So in the current state the contact-line model **shapes the interface** (via
`extend.f90`) and **measures** the capillary force, but does not feed it back
into the particle dynamics. Uncommenting those six lines is the switch.

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
   repeated normalisations (`/ (norm + epsilon(1._rp))`) across 5 iterations per
   step. Appears as faint concentric "levels" following the `alphac` shells.
   Physically meaningless — 15 orders below any real volume fraction. Already
   diagnosed and dismissed with the user. Clip with
   `where(abs(psi) < 1e-12) psi = 0._rp` only if it pollutes a diagnostic.
   *Would stop being negligible under `SINGLE_PRECISION=1`* (`~1e-7`).

2. **Band mismatch.** `extend.f90` uses `alphac > 0.5`, `rotnorm.f90` uses
   `alphac > 0`. The force is integrated over a wider shell than the one the
   interface is relaxed on.

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
| `max_pseudo_iter` | `main.f90:638` (and `503`) | hard-coded `5`. More = stronger enforcement, more round-off |
| `dtau` | `main.f90:637` (and `502`) | hard-coded `0.3*minval(dli)`. See the scaling caveat above |
| `alpha_min` | `extend.f90:101` | hard-coded `0.5`, the relaxation band threshold |

The three hard-coded values are prime candidates for promotion to `input.nml` if
the user wants to sweep them — that would follow the same pattern as commits
`183391a` and `7fd5ccb`, which moved particle and collision constants to the
namelist.
