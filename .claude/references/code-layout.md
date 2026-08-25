# Reference: file-by-file map of `src/`

Every source file, what it owns, and whether it is live in the current build
(`INTERFACE_CAPTURING_VOF=1`, `SDF_NORMALS=1`, `CONSTANT_COEFFS_POISSON=1`,
`PARTICLE=1`, `EULER=1`, `GPU=0`).

Legend: **LIVE** = compiled and executed · *dead* = compiled out or unreachable
in this configuration.

---

## Driver

| file | module | notes |
|---|---|---|
| `main.f90` | `program cans` | **LIVE**. Everything: allocation, init, timeloop, output, finalisation. ~950 lines. No subroutines — it is one long program body with `#include` of the `out*.h90` files. |

---

## Core solver (inherited from CaNS / CaNS-Fizzy)

| file | module | status | role |
|---|---|---|---|
| `types.f90` | `mod_types` | **LIVE** | defines `rp` (real kind), `MPI_REAL_RP`. Fork added a few aliases. |
| `param.f90` | `mod_param` | **LIVE** | all `&dns`/`&scalar`/`&two_fluid` namelist variables + `read_input`. Fork added `theta`, `nh_wide`, `is_ibm`. |
| `common_mpi.f90` | `mod_common_mpi` | **LIVE** | rank, ierr, neighbour ranks, `prt_comm_cart`, `boundleftmyid/boundfrontmyid`. |
| `common_cudecomp.f90` | `mod_common_cudecomp` | *dead* (GPU=0) | cuDecomp handles/streams. |
| `initmpi.f90` | `mod_initmpi` | **LIVE** | pencil decomposition, `lo/hi/n`, `is_bound`, neighbour topology. |
| `initgrid.f90` | `mod_initgrid` | **LIVE** | z-grid (uniform or clustered per `gtype`,`gr`). |
| `initflow.f90` | `mod_initflow` | **LIVE** | `initflow` (velocity IC per `inivel`), `initscal`. |
| `initsolver.f90` | `mod_initsolver` | **LIVE** | eigenvalues + tridiagonal coefficients for the FFT Poisson solver. |
| `sanity.f90` | `mod_sanity` | **LIVE** | `test_sanity_input` — BC/decomposition consistency checks. |
| `bound.f90` | `mod_bound` | **LIVE** | `boundp` (scalars), `bounduvw` (velocity), `updt_rhs_b`. All halo exchange. Fork extended it. |
| `rk.f90` | `mod_rk` | **LIVE** | `rk` (momentum), `rk_scal`, `rk_2fl` (phase field). Called as `tm`/`tm_scal`/`tm_2fl`. |
| `mom.f90` | `mod_mom` | **LIVE** | advective + diffusive momentum operators, surface tension, gravity. |
| `two_fluid.f90` | `mod_two_fluid` | **LIVE** | `init2fl` (psi IC per `inipsi`), `cmpt_norm_curv_youngs`, `clip_field`, `read_sphere_file`. |
| `vof_thinc_qq.f90` | `mod_vof_thinc_qq` | **LIVE** | THINC/QQ VOF: `vof_thinc_transport_psi`, `vof_thinc_cmpt_phi`. **The active interface-capturing scheme.** |
| `acdi.f90` | `mod_acdi` | *dead* | ACDI diffuse interface. Only compiled when `INTERFACE_CAPTURING_VOF=0`. **Do not diagnose interface behaviour from this file.** |
| `fillps.f90` | `mod_fillps` | **LIVE** | builds the Poisson RHS from the velocity divergence. |
| `solver.f90` | `mod_solver` | **LIVE** | FFT-based constant-coefficient Poisson solver (CPU). |
| `solver_gpu.f90` | `mod_solver_gpu` | *dead* (GPU=0) | cuFFT version. |
| `solver_vc.f90` | `mod_solver_vc` | *dead* | variable-coefficient (HYPRE) solver; only when `CONSTANT_COEFFS_POISSON=0`. |
| `correc.f90` | `mod_correc` | **LIVE** | velocity correction with the density-weighted pressure gradient. |
| `updatep.f90` | `mod_updatep` | **LIVE** | pressure update for the split/incremental scheme. |
| `fft.f90` / `fftw.f90` | `mod_fft` / `mod_fftw_param` | **LIVE** | FFTW guru plans, transform-type selection from BCs. |
| `workspaces.f90` | `mod_workspaces` | *dead* (GPU=0) | GPU workspace sizing. |
| `chkdt.f90` | `mod_chkdt` | **LIVE** | stability limit: CFL, viscous, **capillary**, and interface-scheme constraints. |
| `chkdiv.f90` | `mod_chkdiv` | **LIVE** | divergence monitor; aborts the run when `divmax > small`. |
| `forcing.f90` | `mod_forcing` | conditional | `lscale_forcing` for forced HIT (`is_forced_hit`). |
| `scal.f90` | `mod_scal` | *dead* | scalar transport; needs `-D_SCALAR`, not set. |
| `output.f90` | `mod_output` | **LIVE** | `out0d/1d/2d/3d`, `write_visu_2d/3d`, `gen_alias`. |
| `load.f90` | `mod_load` | **LIVE** | `load_one` — MPI-IO checkpoint read/write. |
| `post.f90` | `mod_post` | **LIVE** | statistics helpers used by `out1d.h90`. |
| `utils.f90` | `mod_utils` | **LIVE** | `bulk_mean_12_stag`, `device_memory_footprint`. |
| `debug.f90` | `mod_debug` | conditional | checksum/dump helpers under `-D_DEBUG`. |
| `timer.f90` | `mod_timer` | **LIVE** | `-D_TIMING` wall-clock instrumentation. |
| `nvtx.f90` | `mod_nvtx` | *dead* | NVIDIA profiler ranges. |
| `rotnorm.f90` | `mod_rotnorm` | **LIVE** | **fork addition** — capillary force at the contact line. See `contact-line-model.md`. |
| `extend.f90` | `mod_extend` | **LIVE** | **fork addition** — contact-angle extension velocity + upwind relaxation. |

---

## Particle module (`prt_*.f90`, all behind `-D_PARTICLE`)

| file | module | status | role |
|---|---|---|---|
| `prt_param.f90` | `prt_mod_param` | **LIVE** | `&particle`, `&collision_parameters`, `&particle_euler` namelists. Derives `volp`, `mominert`, spring/damper constants `kn_ss/kt_ss/etan_ss/...`, and the full table of lubrication coefficients. |
| `prt_common.f90` | `prt_mod_common` | **LIVE** | `type particle` / `particle_old` / `particle_sumrk`; arrays `ep`, `op`, `tp`, `rkp`; fields `alphac`, `norm_partx/y/z`, `uphase/vphase/wphase`; `prt_InitMemo` allocates everything and builds `prt_comm_cart`. |
| `prt_initparticles.f90` | `prt_mod_initparticles` | **LIVE** | initial placement. `np==1` uses `x_ini,y_ini,z_ini` from input; `np>1` places pseudo-randomly with overlap rejection. |
| `prt_digitiser.f90` | `prt_mod_digitiser` | **LIVE** | `digitiser(ds, n, alpha)` — the **diffuse solid indicator**: `alpha = 0.5*(1-tanh(ds/(sigma*lambda*delta)))` inside `eps_sol*delta`, else 0. `lambda = |nx|+|ny|+|nz|`, `sigma = 0.05*(1-lambda^2)+0.3`. |
| `prt_initeul.f90` | `prt_mod_initeul` | **LIVE** | rebuilds `alphac` **and** `norm_partx/y/z` at cell centres for the current particle position. Called once per step before the contact-line loop. |
| `prt_initvof.f90` | `prt_mod_initvof` | **LIVE** | builds `uphase/vphase/wphase = (1-alpha)*u_fluid + alpha*u_rigidbody` on the **staggered** faces, for phase-field advection. |
| `prt_eulint.f90` | `prt_mod_eulint` | **LIVE** | the Eulerian IBM: forces `u,v,w` toward rigid-body motion inside the particle and accumulates the reaction force. Also carries the surface-tension term at lines ~428–441. |
| `prt_intgr_over_sphere.f90` | `prt_mod_intgr_over_sphere` | **LIVE** | volume integrals over the sphere: `intu/intv/intw` (momentum) and `intrhox/intrhoy/intrhoz` (density-weighted). `cas` = 1/2/3 selects the component. |
| `prt_intgr_nwtn_eulr.f90` | `prt_mod_intgr_nwtn_eulr` | **LIVE** | **Newton–Euler integration** of every particle: assembles IBM force, gravity/buoyancy, collision and lubrication contributions, updates `x,y,z,u,v,w,om*`, handles master/slave migration, writes `forces_data.csv`. ~1660 lines, the largest file. |
| `prt_collisions.f90` | `prt_mod_collisions` | **LIVE** | `collisions` (soft-sphere DEM, normal+tangential springs/dashpots, Coulomb cap) and `lubrication` (Stokes amplification, saturated below `eps_sat`). |
| `prt_output.f90` | `prt_mod_output` | **LIVE** | `outpart(istep)` — parallel MPI-IO of particle state. |
| `prt_loadpart.f90` | `prt_mod_loadpart` | **LIVE** | `loadpart('r'/'w')` — particle restart. |
| `prt_phase_indicator.f90` | `prt_mod_phase_indicator` | **LIVE** (unused) | sharp indicator `gamu/gamv/gamw/gamp` + rigid-body velocity fields. Only called from the **commented-out** block in `prt_out.h90`. |
| `prt_coordsfp.f90` | `prt_mod_coordsfp` | **LIVE** (vestigial) | Lagrangian forcing-point coordinates. Still `call`ed from `main.f90:551`, but the arrays it fills are only consumed by the `!_EULER` path. |
| `prt_interp_spread.f90` | `prt_mod_interp_spread` | *dead* | `eulr2lagr`/`lagr2eulr` — Lagrangian interpolation/spreading. `#if !defined(_EULER)`. |
| `prt_forcing.f90` | `prt_mod_forcing` | *dead* | `complagrforces`/`updtlagrforces`. `#if !defined(_EULER)`. |
| `prt_kernel.f90` | `prt_mod_kernel` | *dead* | regularised delta kernel + `kerneltest`. `#if !defined(_EULER)`. |

---

## Include files (`.h90`) — compiled into `main.f90`

| file | purpose |
|---|---|
| `out1d.h90` | which 1-D profiles/statistics to write at `iout1d` |
| `out2d.h90` | which 2-D slices to write at `iout2d` — fork added the `Alpha_C` slice |
| `out3d.h90` | which 3-D fields to write at `iout3d` — fork swapped `Kappa` for `Alpha_C` |
| `prt_out.h90` | particle output hook; body is mostly commented out, ends with `call outpart(istep)` |

**Editing any `.h90` requires a full recompile of `main.f90`.**

---

## Other input files shipped in `src/`

| file | purpose |
|---|---|
| `input.nml` | the default/template namelist, copied next to the executable |
| `spheres.in` | sphere list for `inipsi = 'bub3'/'drp3'/...` — read by `read_sphere_file` |
| `cudecomp.in` | cuDecomp autotuning options (GPU only) |
| `dns.in` | legacy CaNS input, copied by the Makefile |

---

## Outside `src/`

| path | purpose |
|---|---|
| `src/poslfp/` | standalone program: post-processing of Lagrangian forcing points. Own Makefile. |
| `PostPrt/PrtPos/` | standalone: reads particle trajectory binaries. Own Makefile + `param.dat`. |
| `PostPrt/SphereData/` | standalone: force decomposition and interpolation on the sphere surface. Own Makefile + `param.f90` + `run.sh`. |
| `utils/visualize_fields/` | Python/XDMF generation from the `log_visu_*.out` logs. |
| `utils/read_binary_data/` | readers for the `.bin` field format. |
