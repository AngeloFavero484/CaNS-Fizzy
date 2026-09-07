# Reference: planned changes — decided, not yet implemented

Modifications the user has settled on but has **not** asked to be written yet.
**Do not implement anything on this page unless explicitly asked.** It exists so a
future session does not have to re-derive the analysis.

---

## 1. Replace `ratiorho` with the live solid / inner-fluid density ratio

*Raised 2026-09-07. Code untouched as of that date.*

### The complaint

`ratiorho` is a constant read from `&particle`. With two fluids the density of
the fluid occupying the particle volume changes in time as the interface sweeps
across the sphere, so a constant `ρ_s/ρ_f` is meaningless. The intent is to
replace it with the **instantaneous** ratio

```
ratiorho(t) = rho_s / rho_f(t),     rho_f(t) = (mass of fluid inside the particle) / V_p
```

### What `ratiorho` actually does today

It is **not** a general density ratio — it appears in exactly two places, and
both are about the *collision* path:

| site | use |
|---|---|
| `prt_param.f90:246,258` | `meffn_ss = ratiorho*volp/2`, `meffn_sw = ratiorho*volp` → the DEM effective mass, from which `kn_*`, `kt_*`, `etan_*`, `etat_*` are derived **once**, at input-read time |
| `prt_intgr_nwtn_eulr.f90:689,699,708` and `:734,744,753` | denominator of the collision-force term in the velocity update: `(ep(p)%colfx+op(p)%colfx)/(ep(p)%vol*ep(p)%ratiorho)` |

`ep(p)%ratiorho` is filled from the namelist value in two places
(`prt_initparticles.f90:545`, `prt_coordsfp.f90:115`), both **before** the time
loop — `coordsfp` is called once, at `main.f90:573`. Nothing rewrites it later.

Everywhere *else* in the momentum update the mass is the honest `ep(p)%vol*rho_s`.

### Why the current code is self-consistent, and what that implies

Two contributions are summed into the same `ep(p)%colf*` accumulator, and they
scale differently:

**(a) The DEM spring/dashpot** (`prt_collisions.f90:108`):
`fxn = -kn·δ·nx - etan·v·nx` with `kn = (π² + ln²e)·meffn/(N·dt_estim)²` and
`meffn ∝ ratiorho·volp`. Dividing by `volp·ratiorho` gives

```
a_spring = -(π² + ln²e)·δ / (Nstretch·dt_estim)²
```

— **`ratiorho` cancels exactly.** The contact time and the restitution `en` are
therefore completely insensitive to its value. This is why nobody has noticed
that `ratiorho = 5` sits next to `rho_s = 1320` in the same namelist.

**(b) The lubrication correction** (`prt_collisions.f90:466`) with
`coeff_f = 6π·(mu12(1)/rho12(1))·radius` (`prt_common.f90:198`). Here `colf` is a
force *divided by a fluid density*, so dividing by `volp·(ρ_s/ρ_f)` restores

```
a_lub = 6π·μ₁·(…)·R / (V_p·ρ_s)
```

— correct **only because** `ratiorho` means `ρ_s/ρ₁` and `coeff_f` was divided by
that same `ρ₁`. The source even admits it: `! mu12(1) is a temporary solution to
be used instead of visc`.

So `ratiorho` is a real physical quantity for the lubrication term and a pure
dummy for the spring term. **Any change must respect that split**, and there are
two traps that follow directly from it:

> **Trap 1 — the DEM constants go stale.** `meffn/kn/kt/etan/etat` are computed
> once in `read_particle_input`. Make `ratiorho` time-dependent without touching
> them and the cancellation in (a) breaks: the spring is tuned for the namelist
> mass but the acceleration is divided by the live mass, so the realised
> restitution is no longer `en`. For a sphere straddling a water/air interface
> `ρ_f(t)` swings by ~10³, so this is not a small error.
>
> **Trap 2 — the lubrication force silently changes.** Substituting `ρ_s/ρ_f(t)`
> while `coeff_f` still carries `ρ₁` leaves a spurious factor `ρ_f(t)/ρ₁`, i.e. an
> effective viscosity `μ₁·ρ_f(t)/ρ₁` that is not the viscosity of the fluid in
> the gap. Fixing `ratiorho` alone makes the lubrication model *worse*, not
> better.

### The quantity to use — it already exists

`intgr_over_sphere(4,…)` integrates `ρ = ρ₂ + (ρ₁-ρ₂)ψ` over the particle and
stores the result in `ep(p)%intrhox/intrhoy/intrhoz`
(`prt_intgr_over_sphere.f90:853-855`; the commented print at :891 calls it
"Mass of the fluid occupied by the particle"). It is already called every
timestep at `prt_intgr_nwtn_eulr.f90:180`, and the buoyancy term already consumes
it as `gacc(i)·(1 - intrho_i/(vol·rho_s))`. So

```
rho_f(t) = intrho / V_p
```

needs no new integral and no new field in the halo exchange — but **not in the
form it is stored in today**. See the next section.

### Why `intrhox/y/z` are not three copies of one number — add `intrhoc`

*(User correction, 2026-09-07. The earlier version of this note got this wrong.)*

The three components are **not** the same integral differing by round-off. They
are integrals over **three different staggered control volumes**:

- each sweep samples the sphere on its own lattice — u-faces
  (`prt_intgr_over_sphere.f90:374–450`), v-faces (`:496–570`), w-faces
  (`:617–690`) — so `dist2`, and with it the digitised `alpha_eul`, is evaluated
  at a **different set of points relative to the particle centre** in each case;
- the box swept in each direction runs out to `R + eps_sol·dl`
  (`:318–323`, with `eps_sol = 1.5`), i.e. over the whole **diffuse** shell, and
  `alpha` tapers through it (`prt_digitiser.f90:26`). The enclosed volume is
  therefore a digitised, shell-smeared volume, not `4/3πR³`, and it comes out
  slightly different on each of the three lattices;
- on top of that, `ρ` is interpolated **onto the face** in each sweep —
  `rhox = ρ₂ + Δρ·0.5(psi(i,j,k)+psi(i+1,j,k))` (`:374`), and likewise `rhoy`
  (`:496`), `rhoz` (`:617`) — a second, direction-dependent averaging of ψ.

That is exactly what makes them right for their present job: each momentum
component is divided by the mass on **its own** control volume. It also makes
them the wrong basis for a single scalar "density of the fluid inside the
particle", and averaging the three would just blend three different volumes.

**Planned fix: add a cell-centred `intrhoc`** — a new `cas` in
`intgr_over_sphere` sampling at cell centres, where ψ natively lives:

```
rho_c(i,j,k) = rho12(2) + (rho12(1)-rho12(2))*psi(i,j,k)   ! no face interpolation
intrhoc      = sum over the box of  rho_c * alpha_c * dV
intvolc      = sum over the box of           alpha_c * dV   ! same sweep
rho_f(t)     = intrhoc / intvolc
```

Two things this buys, both of which retire caveats from the earlier draft:

- **ψ is not interpolated.** The cell-centred sweep reads `psi(i,j,k)` directly,
  so `rho_f` carries none of the extra face-averaging smoothing that
  `intrhox/y/z` do.
- **`rho_f` is bounded by construction.** Dividing by `intvolc` from the *same*
  sweep makes `rho_f` an `alpha`-weighted mean of `rho_c`, with `alpha ≥ 0`.
  It therefore lies in `[min(rho12), max(rho12)]` exactly, for any digitiser
  error and any shell width. Dividing `intrho` by the analytic `volp` does not
  — that mixes a digitised numerator with an analytic denominator and can drift
  outside the physical range. Compute `intvolc` in the same sweep; do not reach
  for `volp`, and do not revive the staggered `cas == 1`
  (`prt_intgr_nwtn_eulr.f90:177`) for this.

**Do not shortcut it through `alphac`.** The cell-centred solid indicator already
exists as a field and is rebuilt every step (`prt_initeul.f90:340–348`, called
from `main.f90:663,685`), but it is unusable here for two reasons: it is
digitised on a **retracted** radius `radius-retrac` (`:339`; `retrac = 0`
today at `prt_common.f90:179`, but it is a live knob), and it is **accumulated
over all particles and capped at 1** (`:344–348`), so it cannot be attributed to
one particle when `np > 1`. A new `cas` inside `intgr_over_sphere` inherits the
master/slave halo bookkeeping and the neighbour reduction for free; a raw field
sum would need its own MPI reduction over the particle's box.

**MPI/ownership caveat, unchanged:** `intrho*` is set **only on the master rank**
(`:857–859` zeroes it everywhere else), and `ratiorho` sits *outside* the
`send_real = 24+10*nqmax` contiguous block that gets exchanged
(`prt_common.f90:79` — it is after the communicated region). `intrhoc` must be
declared in the same non-communicated part. Since the master both computes and
consumes the value inside one timestep this is fine — but only if it is
recomputed every step, never stored and carried across a mastership change.

### Recommended implementation

Do **not** just swap the number. Remove the mass-scale role from `ratiorho`
entirely and let the density ratio appear only where physics asks for it:

1. **Split the accumulator.** Keep `ep%colf*` for the DEM contact force (a true
   force, in N) and add `ep%lubf*` for the lubrication force. They are two
   different scalings sharing one variable today.
2. **Make lubrication a true force.** `coeff_f = 6π·μ_f·R`, `coeff_t = 8π·μ_f·R²`
   with `μ_f` the viscosity of the fluid in the gap — computable exactly like
   `intrho` by adding a `cas` that integrates `μ = μ₂ + (μ₁-μ₂)ψ`, or by
   evaluating `μ` from the local `ψ` at the contact point.
3. **One denominator for everything:** `ep(p)%vol*rho_s`, the true particle mass,
   which is already what every other term in the update uses.
4. **Rebuild the DEM constants on `rho_s*volp`** in `read_particle_input`. The
   cancellation of (a) is then exact and permanent, and `en`/`et` are honoured
   for free — no per-step recomputation of `kn`/`etan`, no per-particle copies.
5. `ratiorho` now appears nowhere. Delete it from `&particle`, from
   `prt_mod_param`, from the `particle` type (`prt_common.f90:79`), and from its
   two assignment sites. `ρ_f(t)` still enters the physics — through buoyancy
   (already via `intrho`) and through `μ_f` in the lubrication model.

This is strictly more than "substitute the variable", but it is the version where
the answer is right, and it removes a parameter instead of adding a subtler one.

### Minimal variant, if the literal substitution is preferred

Keep the structure and make `ep(p)%ratiorho` live:

- compute it from the new cell-centred `intrhoc/intvolc` right after the
  `intgr_over_sphere` calls (`prt_intgr_nwtn_eulr.f90:178-180`) and **before**
  the `do r = 1,r_dtcol` loop at `:205`, so it stays frozen over the macro step
  like the other `int*` quantities (see the "new --> old" comment near `:860`);
- guard it: skip non-masters, and clamp `intvolc` away from zero — the `int*`
  are identically 0 on a non-master and on the first step before the first
  integration (a zero *volume*, not a zero density, is the division hazard once
  `rho_f = intrhoc/intvolc` is bounded by construction);
- **still** fix Trap 1 (move `meffn/kn/etan` per-particle and recompute each
  step) and Trap 2 (`coeff_f`), otherwise the restitution and the lubrication
  force are both silently wrong.

### Touch points checklist

```
src/prt_param.f90            :44,175,187        declaration / namelist / default
src/prt_param.f90            :246-263           meffn,kn,kt,etan,etat derivation
src/prt_common.f90           :79                ratiorho in the particle type
src/prt_common.f90           :198-199           coeff_f, coeff_t  (mu12(1)/rho12(1))
src/prt_collisions.f90       :108-110,235       DEM spring -> colf
src/prt_collisions.f90       :454-466           lubrication -> colf
src/prt_initparticles.f90    :545               ep%ratiorho = ratiorho
src/prt_coordsfp.f90         :115               ep%ratiorho = ratiorho  (once, pre-loop)
src/prt_intgr_nwtn_eulr.f90  :178-180           intgr_over_sphere calls; add the new cas
src/prt_intgr_over_sphere.f90:18,374,496,617    cas dispatch + the three staggered sweeps
src/prt_intgr_over_sphere.f90                   NEW: cell-centred sweep -> intrhoc, intvolc
src/prt_common.f90           :79                NEW: intrhoc,intvolc in the particle type,
                                                outside the send_real block
src/prt_intgr_nwtn_eulr.f90  :689,699,708       collision term, colliding branch
src/prt_intgr_nwtn_eulr.f90  :734,744,753       collision term, non-colliding branch
src/prt_intgr_over_sphere.f90:853-859           intrho, master-only
src/input.nml + examples/*/input.nml            drop/replace ratiorho
.claude/references/input-namelists.md           the &particle table
```

Restart files carry no `ratiorho` (it is outside the communicated block and is
re-derived), so `prt_loadpart.f90` needs nothing — but a derived `ρ_f` must be
recomputed before its first use on a restarted step, not read back.

Worth adding a `rho_f` column to `forces_data.csv` while doing this: it is the
one number that says whether the sphere currently thinks it is in the liquid,
in the gas, or straddling. `intrhoc/intvolc` is also the honest thing to plot
against `F_buoy` — note that the buoyancy term keeps using the staggered
`intrhox/y/z`, and should: each momentum component wants the mass on its own
control volume. `intrhoc` is for the scalar density ratio only, not a
replacement for them.
