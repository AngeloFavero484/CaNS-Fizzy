# Reference: `input.nml` — every runtime parameter

Read by two routines:
- `src/param.f90::read_input` → `&dns`, `&scalar`, `&two_fluid`, `&contact_line`, `&cudecomp`
- `src/prt_param.f90::read_particle_input` → `&particle`, `&collision_parameters`, `&particle_euler`

Both open the **same file** `input.nml`, in the working directory of the
executable (normally `run/`). A missing namelist group is not fatal — the
defaults listed below apply. Upstream docs: [`docs/INFO_INPUT.md`](../../docs/INFO_INPUT.md).

---

## `&dns` — grid, time, boundary conditions

| parameter | default | meaning |
|---|---|---|
| `ng(1:3)` | `128,128,128` | global grid points |
| `l(1:3)` | `1.,1.,1.` | domain size (physical units — this fork uses **dimensional** SI-like values) |
| `gtype` | `1` | z-grid clustering type; `1` = uniform |
| `gr` | `0.` | clustering strength; `0.` = uniform |
| `cfl` | `0.95` | safety factor on the computed stability limit |
| `dtmax` | `1.e9` | hard cap on `dt` |
| `dt_f` | `-1.` | **fixed** timestep; `<0` = adaptive. If `dt_f` exceeds the stability limit the run **aborts** at the next `icheck`. |
| `is_solve_ns` | `T` | solve Navier–Stokes; `F` re-imposes `initflow` each step (prescribed-velocity mode) |
| `is_track_interface` | `T` | advance `psi`. **`F` also disables the whole contact-line block.** |
| `inivel` | `'zer'` | velocity IC (`'zer'`, `'poi'`, `'log'`, …) |
| `is_wallturb` | `F` | superimpose wall-turbulence perturbations |
| `is_forced_hit` | `F` | large-scale forcing for homogeneous isotropic turbulence |
| `nstep` | `1000` | max timesteps |
| `time_max` | `1.` | max physical time |
| `tw_max` | `0.5` | max wall-clock **hours** |
| `stop_type(1:3)` | `T,F,F` | which of `nstep` / `time_max` / `tw_max` are active |
| `restart` | `F` | resume from `data/fld_*.bin` + particle restart |
| `is_overwrite_save` | `T` | overwrite the single checkpoint instead of numbering them |
| `nsaves_max` | `0` | cap on retained numbered checkpoints; `0` = unlimited |
| `icheck` | `10` | steps between `chkdt` + `chkdiv` |
| `iout0d` | `10` | steps between scalar output (`time.out`, **`forces_data.csv`**) |
| `iout1d` | `100` | steps between 1-D profiles |
| `iout2d` | `1000` | steps between 2-D slices |
| `iout3d` | `500` | steps between 3-D fields |
| `isave` | `1000` | steps between checkpoints |
| `cbcvel(0:1,1:3,1:3)` | `'P'` | velocity BC **type**: `[lower,upper] × [x,y,z] × [u,v,w]`. `'P'` periodic, `'D'` Dirichlet, `'N'` Neumann |
| `cbcpre(0:1,1:3)` | `'P'` | pressure BC type |
| `bcvel`, `bcpre` | `0.` | BC **values**, same index layout |
| `bforce(1:3)` | `0.` | uniform body force |
| `gacc(1:3)` | `0.` | gravity vector — e.g. `0.,0.,-9.81` |
| `dims(1:2)` | `0,0` | MPI pencil decomposition in x,y. `0,0` = let MPI decide. Product must equal the rank count. |

**BC gotcha:** `cbcpre` must be `'N'` where `cbcvel` is `'D'` (solid wall), and
`'P'`/`'P'` together for periodic. `test_sanity_input` catches most mistakes.

---

## `&scalar` — passive/active scalar (needs `-D_SCALAR`, currently OFF)

| parameter | default | meaning |
|---|---|---|
| `inisca` | `'zer'` | scalar IC |
| `cbcsca(0:1,1:3)`, `bcsca` | `'P'`, `0.` | scalar BCs |
| `ssource` | `0.` | uniform source term |

---

## `&two_fluid` — the two fluid phases and the interface

| parameter | default | meaning |
|---|---|---|
| `inipsi` | `'uni'` | `psi` IC. See table below. |
| `cbcpsi(0:1,1:3)`, `bcpsi` | `'P'`, `0.` | phase-field BCs |
| `cbcnor(0:1,1:3,1:3)`, `bcnor` | `'P'`, `0.` | interface-normal BCs (one set per component) |
| **`theta`** | `90°` (`2*atan(1)` rad… stored in **degrees** at input) | **contact angle in DEGREES.** Fork-specific. Consumed by `extend.f90` and `rotnorm.f90`, both of which do `theta*pi/180`. |
| `sigma` | `0.` | surface tension coefficient |
| `rho12(1:2)` | `1.,1.` | densities of phase 1 (`psi=1`) and phase 2 (`psi=0`) |
| `mu12(1:2)` | `0.01` | dynamic viscosities, same ordering |
| `ka12`, `cp12`, `beta12` | `0.01`,`1.`,`1.` | conductivity / heat capacity / expansion (scalar runs only) |
| `psi_thickness_factor` | `0.51`, but **`0.50` under `_INTERFACE_CAPTURING_VOF`** | interface thickness. **Active in both schemes.** Under ACDI it sets `seps = max(dl)*factor`. Under VOF it sets the THINC sharpness: `vof_thinc_beta = 1/(2*psi_thickness_factor)` (`param.f90:183`) — so *smaller* factor means *larger* beta means a *sharper* interface. |

### `inipsi` options (`two_fluid.f90::init2fl`)

| value | shape |
|---|---|
| `'uni'` | uniform `psi = 1` everywhere |
| `'zer'` | uniform `psi = 0` everywhere |
| `'bub3'` / `'drp3'` / `'dis3'` | spheres from `spheres.in` (bubble = phase 2 inside, drop = phase 1 inside) |
| `'bub2'` / `'drp2'` / `'dis2'` | cylinders (2-D) from `spheres.in` |
| `'bub1'` / `'drp1'` / `'dis1'` | planar films from `spheres.in` |
| `'flm'` | a flat film / horizontal interface — used by the bouncing/sinking-sphere cases |
| `'cap-wav-1d'` | capillary-wave validation |
| `'zalesak-disk'` | Zalesak rotating-slotted-disk advection test |

Cases using a `*3`/`*2`/`*1` option **require a `spheres.in`** next to `input.nml`.

---

## `&contact_line` — numerics of the extended contact-line model

Fork-only. These three drive the pseudo-time relaxation in `main.f90` that
imposes `theta` at the particle surface (see
[`contact-line-model.md`](contact-line-model.md)). They were hard-coded until
they were promoted to the namelist; **the defaults are exactly the old
hard-coded values**, so an `input.nml` without this group behaves as before.

| parameter | default | meaning |
|---|---|---|
| `max_pseudo_iter` | `5` | relaxation iterations per timestep. More = stronger enforcement of `theta`, but more of the machine-epsilon `psi` round-off noted in `contact-line-model.md`. |
| `dtau_cfl` | `0.3` | pseudo-timestep as a CFL number on the smallest cell: `dtau = dtau_cfl/maxval(dli)`. `u_ext` is a unit vector, so this is a true CFL. Raising it past ~0.5 risks the upwind advection going unstable. |
| `alpha_min` | `0.5` | lower edge of the `alphac` band the relaxation acts on (band is `alpha_min < alphac < 1`). Applies to **both** `compute_uextend` and `advect_vof_upwind` in `extend.f90`, which must agree. |

Note `alpha_min` does **not** move the band `rotnorm.f90` integrates the
capillary force over — that one is still `alphac > 0`. The band mismatch
documented in `contact-line-model.md` is therefore unchanged, and lowering
`alpha_min` towards 0 narrows the gap.

The whole block is skipped when `is_track_interface = F`, so these knobs have
no effect in that mode.

---

## `&particle` — the rigid particle(s) (`-D_PARTICLE`)

| parameter | default | meaning |
|---|---|---|
| `np` | `1` | number of particles |
| `radius` | `1.` | sphere radius (same units as `l`) |
| `rho_s` | `1320.` | **solid density**, used in the Newton–Euler force balance |
| `ratiorho` | `5.` | density ratio used by the **collision** effective masses (`meffn_ss = ratiorho*volp/2`). Note it is *independent* of `rho_s` — keep them physically consistent by hand. |
| `u_ini,v_ini,w_ini` | `0,0,-28.78` | initial particle velocity (applied to **all** particles) |
| `x_ini,y_ini,z_ini` | unset sentinel | **only honoured when `np == 1`.** Left unset → legacy default `x=l(1)/2`, `y=l(2)/2`, `z=0.755*l(3)`. Added by commit `7e46c21`. |

`np > 1` ignores `*_ini` positions and places particles pseudo-randomly in
`prt_initparticles.f90` with overlap rejection.

Derived automatically: `volp = 4/3 πr³`, `mominert = 2/5 volp r²`.

---

## `&collision_parameters` — soft-sphere DEM

| parameter | default | meaning |
|---|---|---|
| `Nstretch` | `8.0` | collision time stretched over this many `dt_estim` — sets spring stiffness |
| `dt_estim` | `0.003` | the timestep the collision model is *tuned* for. **Must be kept near the actual `dt`** or the collision becomes either rigid (unstable) or mushy. |
| `r_dtcol` | `50` | ratio `dt / dt_collision` — collision substepping |
| `en` | `0.97` | normal restitution coefficient |
| `et` | `0.10` | tangential restitution coefficient |
| `muc` | `0.0` | Coulomb friction coefficient |

Everything else (`kn_ss`, `kt_ss`, `etan_ss`, `etat_ss`, `kn_sw`, …) is
**derived** in `read_particle_input`:
```
kn   = (π² + |ln e|²) · m_eff / (Nstretch·dt_estim)²
etan = -2·ln(e) · m_eff / (Nstretch·dt_estim)
```
with `m_eff = ratiorho·volp/2` (sphere–sphere) or `ratiorho·volp` (sphere–wall).

The lubrication-model coefficients (`a11_ini_pp`, `a22_sat_pw`, …) are
compile-time `parameter`s in `prt_param.f90`, keyed off
`eps_ini_pp=0.025`, `eps_sat_pp=0.001`, `eps_ini_pw=0.075`, `eps_sat_pw=0.001`.

---

## `&particle_euler` — Eulerian coupling (`-D_EULER`)

| parameter | default | meaning |
|---|---|---|
| `eps_sol` | `1.5` | **width of the diffuse solid shell in cells.** `digitiser` returns `alpha = 0` beyond `eps_sol*delta` from the surface. Larger = smoother `alphac`, thicker contact-line band, more `psi` noise; smaller = sharper but more grid-locking. |

---

## `&cudecomp` — only read under `-D_OPENACC` (GPU=0 here, so ignored)

`cudecomp_t_comm_backend`, `cudecomp_is_t_enable_nccl`,
`cudecomp_is_t_enable_nvshmem`, and the matching `_h_` halo variants.
Harmless to leave in the file on a CPU build.

---

## Cross-parameter consistency rules

1. `product(dims) == number of MPI ranks` (unless `dims = 0,0`).
2. `radius + offset < l(1)/dims(1)` **and** `< l(2)/dims(2)` — otherwise
   `prt_InitMemo` aborts with *"Radius spheres larger than x-dimension of processes"*.
   Fewer ranks, or a bigger domain, fixes it.
3. `dt_estim` ≈ the actual running `dt` (printed each `icheck` as `dt = ...`).
4. `theta` is in **degrees**, `gacc` is signed (gravity is `-9.81`, not `9.81`).
5. `rho12(1)` is the density where `psi = 1`. With `inipsi = 'flm'` and
   `rho12 = 935., 1.`, phase 1 is the liquid pool and phase 2 the gas above it.
