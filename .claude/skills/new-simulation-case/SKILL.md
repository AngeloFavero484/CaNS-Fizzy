---
name: new-simulation-case
description: Create or modify a simulation case (input.nml) for this three-phase particle solver. Use when the user wants a new example case, a new configuration, to change contact angle / density ratio / particle position / grid, or asks "how do I set up X".
---

# Setting up a simulation case

A case is **one directory under `examples/`** containing `input.nml`, plus
`spheres.in` when the phase-field IC needs it. Nothing else — no build.conf
override unless the case genuinely needs different compile flags.

```
examples/Three_Phase/<CaseName>/
├── input.nml
└── spheres.in        # only for inipsi = 'bub3'/'drp3'/'bub2'/'drp1'/...
```

To run: copy `input.nml` (and `spheres.in`) into `run/`, then
`cd run && mpirun -n N ./cans`.

## Start from the closest existing case

| case | what it is |
|---|---|
| `Bouncing_Sphere` | dense sphere falling onto a flat liquid film, bounces. `inipsi='flm'`, `theta=154`, `w_ini=-28.78` |
| `Sinking_Sphere` | same setup, lower surface tension → sphere penetrates |
| `Wall_Collision` | sphere impacting a wall, no interface (`inipsi='zer'`), fine grid `192×192×480` |
| `Head_On` | sphere fired at a droplet, zero gravity, `theta=60` |
| `Sessile_Drop` | droplet resting on a sphere, zero gravity, `theta=150` |
| `Particle_Capture` | **`np=500`** particles + bubble, the only multi-particle case |

Copy the nearest one and edit — do not write an `input.nml` from scratch.

## The namelist groups

Full parameter tables: `.claude/references/input-namelists.md`. The groups are
`&dns`, `&scalar`, `&two_fluid`, `&cudecomp`, `&particle`,
`&collision_parameters`, `&particle_euler` — all in the **same file**, read by
two different routines (`param.f90` and `prt_param.f90`).

## Checklist for a new case

1. **Grid and domain.** `ng(1:3)` and `l(1:3)`. Resolution across the particle is
   `2*radius / (l(1)/ng(1))` — aim for **at least ~16 cells per diameter**, and
   remember the diffuse solid shell eats `eps_sol` (~1.5) cells on each side.

2. **Particle.** `radius`, `rho_s`, `ratiorho`, `u/v/w_ini`, and for `np=1`
   the explicit `x_ini,y_ini,z_ini`. Leaving positions unset falls back to the
   legacy default (`l/2, l/2, 0.755*l(3)`) — the recent commits `5044311`/`344a400`
   exist specifically to stop relying on that, so **always set them explicitly**.
   Keep `rho_s` and `ratiorho` physically consistent; nothing checks them.

3. **Fluids.** `rho12(1:2)`, `mu12(1:2)`, `sigma`. Index 1 is the phase where
   `psi = 1`.

4. **Contact angle.** `theta`, **in degrees**. `< 90` wetting, `> 90` non-wetting.

5. **Phase-field IC.** `inipsi`:
   - `'zer'` / `'uni'` — uniform, no interface (single-fluid + particle)
   - `'flm'` — flat horizontal film
   - `'drp3'` / `'bub3'` — spheres, **requires `spheres.in`**
   - full list in `.claude/references/input-namelists.md`

6. **Boundary conditions.** `cbcvel`/`cbcpre` must be consistent: `'D'`+`'N'` at
   a wall, `'P'`+`'P'` periodic. `test_sanity_input` catches most errors at startup.

7. **Gravity.** `gacc(1:3)`, signed — `0.,0.,-9.81`, not positive.

8. **Collision parameters.** `dt_estim` must be **close to the actual running
   `dt`** (printed each `icheck`). If it is off by an order of magnitude, the
   collision is either numerically rigid (blows up) or mushy (particle sinks
   into the wall). Run a few steps, read the reported `dt`, then set `dt_estim`.

9. **Decomposition.** `dims(1:2) = 0,0` is fine. But
   `radius + offset` must be `< l(1)/dims(1)` and `< l(2)/dims(2)`, or
   `prt_InitMemo` aborts. With a big particle, **use fewer ranks**.

10. **Output cadence.** `iout0d` also controls `forces_data.csv`.
    `iout3d`/`isave` dominate disk use — `iout3d = 0` disables 3-D output.

## `spheres.in` format

One line per sphere, read by `two_fluid.f90::read_sphere_file`. Copy the format
from an existing case (`Head_On/spheres.in`) rather than guessing.

## Sanity-check the physics before running

`main.f90:280-308` prints a dimensionless-parameter banner at startup —
`Ar`, `Eo`, `Oh`, `Mo`, capillary length, and `dt_cap` (the capillary timestep
limit `sqrt((rho1+rho2)*dV / (4*pi*sigma))`). **Read those numbers on the first
run** and confirm they match the regime the user intends. A capillary `dt_cap`
far below the CFL `dt` means surface tension will dominate the timestep cost.

## Registering the case

Add the directory, then commit and push to `myfork` (never `origin`):

```bash
git add examples/Three_Phase/<CaseName>/
git commit -m "Add <CaseName> example case"
git push myfork main
```

Match the existing commit-message style — short imperative, no body needed for a
case addition.
