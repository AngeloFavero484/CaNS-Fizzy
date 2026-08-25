# Reference: the timeloop, in call order

All line numbers are in `src/main.f90` as of commit `26b6dfb`. They drift —
re-grep if they do not match.

---

## Time integration scheme

```fortran
! main.f90:595-602
! Runge-Kutta
!    do irk=1,3
!      tm_coeff(:) = rkcoeff(:,irk)
!      dtrk = sum(rkcoeff(:,irk))*dt
!Adams-Bashfort
    do irk=1,1
      tm_coeff(:) = [2._rp+dt/dto,-dt/dto]/2._rp
      dtrk = dt
```

**The scheme is second-order Adams–Bashforth with a variable-timestep
correction**, not RK3. The `do irk=1,1` loop is a vestigial wrapper kept so the
RK3 body below it still compiles unchanged.

- `tm_coeff(1) = (2 + dt/dto)/2`, `tm_coeff(2) = -(dt/dto)/2`
- With `dt == dto` this reduces to the classic `[3/2, -1/2]`.
- `dto` is the *previous* step's `dt`, updated at `main.f90:789`.
- Every routine taking `rkpar`/`tm_coeff` uses `rkcoeffab = tm_coeff(1)+tm_coeff(2)`
  as the effective full-step weight — which equals `1` exactly. Multiplications
  by `rkcoeffab` throughout `prt_intgr_nwtn_eulr.f90` are therefore no-ops in
  the AB configuration, left over from RK3.

Consequence: **anything in this codebase written assuming three substeps is
stale.** `rkp` (the `particle_sumrk` accumulator in `prt_common.f90`) exists for
that reason and is essentially unused now.

---

## One timestep, in order

### 0. Bookkeeping — `main.f90:592–594`
```
istep = istep + 1 ;  time = time + dt
```
Note `time` is incremented *before* the step is taken.

### 1. Phase-field + contact line — `main.f90:609–667`
Only when `is_track_interface = T`. Otherwise the `psiflx_*` arrays are zeroed
and the whole block is skipped.

```fortran
psio(:,:,:) = psi(:,:,:)                                   ! 610  save old psi (needed by rk)

uphase = u ; vphase = v ; wphase = w                       ! 613-615
call initvof(n,u,v,w,uphase,vphase,wphase)                 ! 616  blend in rigid-body velocity
call bounduvw(...uphase,vphase,wphase)                     ! 617

call tm_2fl(tm_coeff,...,psi,psiflx_x,psiflx_y,psiflx_z)   ! 618  THINC/QQ transport of psi
call bounduvw(... psiflx_* ...)                            ! 620
call boundp(cbcpsi,...,psi)                                ! 621

call vof_thinc_cmpt_phi(n,vof_thinc_beta,psi,phi)          ! 625  rebuild the SDF
call cmpt_norm_curv(n,dli,dzci,dzfi,phi,...)               ! 628  normals+curvature FROM phi
call boundp(cbccur,...,kappa) ; boundp(cbcnor(:,:,1..3))   ! 632-635

dtau = 0.3_rp * minval(dli(1:3))                           ! 637  <-- see note below
max_pseudo_iter = 5                                        ! 638

call initeul(n)                                            ! 640  rebuild alphac + norm_part*
call boundp(cbcpsi,...,alphac)                             ! 641

do iter = 1, 5                                             ! 642  CONTACT-LINE RELAXATION
  u_ext = 0 ; v_ext = 0 ; w_ext = 0
  call compute_uextend(n,theta,normx,normy,normz,u_ext,v_ext,w_ext)   ! extend.f90
  call advect_vof_upwind(n,dli,dtau,u_ext,v_ext,w_ext,psi)            ! extend.f90
  call boundp(cbcpsi,...,psi)
  call cmpt_norm_curv(n,dli,dzci,dzfi,psi,...)             ! 649  NOTE: on psi, not phi
  call boundp(cbcnor(:,:,1..3)) ; call boundp(...,kappa)
end do

call rot_norm(n,dli,dzci,psi,theta,is_bound,...,Fs)        ! 655  rotnorm.f90 -> capillary force
Fstot_old = Fstot
call MPI_ALLREDUCE(Fs,Fstot,3,MPI_REAL_RP,MPI_SUM,...)     ! 657
```

**Two things to flag about `dtau`:**
`dtau = 0.3 * minval(dli(1:3))` — `dli` is *inverse* spacing, so this is
`0.3 / max(dl)`, i.e. `dtau` grows as the grid is refined. A pseudo-time step
that scales like `1/dx` combined with an upwind advection of `psi` is the
opposite of the usual CFL scaling. It is stable in practice only because
`u_ext` is a unit vector and the loop runs 5 iterations, but it is worth
knowing if the contact line behaves oddly under grid refinement.

The same block is duplicated verbatim at `main.f90:502–522` for the initial
condition, before the timeloop starts.

### 2. Scalar transport — `main.f90:668–671`
`#if defined(_SCALAR)`. **Not compiled** in the current build.

### 3. Momentum predictor — `main.f90:672–697`

```fortran
if (.not. is_solve_ns) then
  call initflow(...)                     ! re-impose the prescribed velocity field
else
  is_cmpt_rho_av = (abs(gacc) > 0) .and. cbcpre(0,:)//cbcpre(1,:) == 'PP'
  call bulk_mean_12_stag(...)            ! 681-685 mean density, only for periodic+gravity dirs
  call tm(tm_coeff,...,psi,kappa,p,pn,po,s,psio,psiflx_*,u,v,w)    ! 687  mom.f90 via rk.f90
  if (is_forced_hit) call lscale_forcing(...)
  call bounduvw(...)
```

`tm` (= `mod_rk::rk`) computes advection, diffusion, surface tension from
`sigma*kappa*grad(psi)`, gravity relative to `rho_av`, and applies the AB
weights to `dudtrko` etc.

### 4. Particle coupling (IBM) — `main.f90:698–749`

```fortran
!$acc update self(u,v,w)
is_ibm = .true.
call eulint(tm_coeff,dt,irk,n,psi,psio,kappa,fx_old,fy_old,fz_old,u,v,w)   ! 734
call boundp(cbcpsi,...,alphac)
call bounduvw(...,u,v,w)
is_ibm = .false.
```

`eulint` (`prt_eulint.f90`) is the **direct-forcing Eulerian IBM**: wherever
`alphac > 0` it drives `u,v,w` toward the local rigid-body velocity and stores
the reaction force in `ep(p)%fxltot/fyltot/fzltot` (and torques). It also carries
a surface-tension contribution at lines ~428–441 using the local density
`rho + drho*psi`.

`fx_old/fy_old/fz_old` persist the forcing between steps (they are checkpointed
as `fexts(6..8)`).

### 5. Pressure projection — `main.f90:750–764`

```fortran
pp = p
call fillps(n,dli,dzfi,dtrki,rho0,u,v,w,p)         ! divergence -> RHS
call updt_rhs_b(['c','c','c'],cbcpre,...)          ! BC contributions
call solver(n,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],p)   ! FFT Poisson
call boundp(cbcpre,...,p)
call correc(n,dli,dzci,rho0,rho12,dtrk,p,psi,u,v,w)   ! velocity correction
call bounduvw(...,.true.,...)
call updatep(pp,p)
call boundp(cbcpre,...,p)
```

`rho0` is the constant reference density for the split; `correc` applies the
density-weighted gradient using the local `psi`.

### 6. Particle motion — `main.f90:765–767`

```fortran
call intgr_nwtn_eulr(n,l,dl,dli,dt,tm_coeff,istep,psi,u,v,w, &
                     Fstot,Fstot_old,F_ibm,F_inertia,F_w,F_buoy,F_cap)
```

Inside (`prt_intgr_nwtn_eulr.f90`), per particle:

```
u^{n+1} = u^n
        + (1-colflg)·[ -½dt(fltot^{n+1}+fltot^n)/(vol·rho_s)      <- IBM reaction
                       + (int^{n+1} - int^n)/(vol·rho_s) ]         <- fluid momentum inside
        + dt·gacc·(1 - intrho/(vol·rho_s))                         <- gravity + buoyancy
        + ½dt(colf^{n+1}+colf^n)/(vol·ratiorho)                    <- collisions
!       + ½dt(Fstot+Fstot_old)/(vol·rho_s)                         <- CAPILLARY: COMMENTED OUT
x^{n+1} = x^n + ½dt(u^{n+1}+u^n)
```

**The capillary feedback `Fstot` is commented out** at lines 675, 685, 694, 725,
735, 744 for both the collision and no-collision branches. `Fstot` is still
*computed*, *reduced*, and *logged* as `F_cap`, but does not act on the particle.
Check this before claiming the contact-line force is coupled.

Also inside: collision detection against other particles and walls, the
`lubrication` correction, master/slave re-assignment as particles cross
subdomain boundaries, and the `forces_data.csv` write every `iout0d` steps.

### 7. Pressure history — `main.f90:770–775`
```
po = pn ; pn = p
```
Needed by the pressure-splitting scheme.

### 8. Stop criteria — `main.f90:779–789`
`nstep` / `time_max` / `tw_max` per `stop_type(1:3)`. Then `dto = dt`.

### 9. Stability + divergence check — `main.f90:809–841`
Every `icheck` steps:
- `chkdt` → new `dt_cfl`; `dt = min(cfl*dt_cfl, dtmax)`, or the fixed `dt_f`
- **aborts** if `dt_f` exceeds the stability limit
- **aborts** if `dt_cfl < small`
- `chkdiv` → **aborts** if `divmax > small` or `divtot` is NaN
  (unless `-D_MASK_DIVERGENCE_CHECK`)

### 10. Output — `main.f90:845–900`
`iout0d` → `time.out` (and `forces_data.csv`, written from inside
`intgr_nwtn_eulr`) · `iout1d` → `out1d.h90` · `iout2d` → `out2d.h90` ·
`iout3d` → `out3d.h90` + `prt_out.h90` · `isave` → checkpoints.

Each output block does `!$acc update self(u,v,w,p,psi,kappa,s)` first.

---

## Initialisation order (before the loop)

```
MPI_INIT → read_input → read_particle_input → initmpi → prt_InitMemo
→ [allocations] → initgrid → grid.bin/grid.out/geometry.out
→ test_sanity_input → initsolver
→ restart? { loadpart + 8× load_one }  :  { initflow, initparticles, initeul, init2fl }
→ bounduvw, boundp, boundp(alphac)
→ vof_thinc_cmpt_phi → cmpt_norm_curv → boundp(kappa,norm*)
→ [contact-line pseudo-loop ×5]  → rot_norm → MPI_ALLREDUCE → Fstot
→ initial output → coordsfp → intgr_over_sphere(1,2,3)
→ chkdt → dt, dto, dti
→ MAIN LOOP
```

`forces_data.csv` is opened with `status='replace'` at `main.f90:181`, **before
`MPI_INIT`** — so every rank truncates and holds the same filename. Only rank 0
writes to it in practice (guarded inside `intgr_nwtn_eulr`), but the open itself
is unguarded.

---

## Restart contents

`restart = T` reads, in this order:
`fld_u`, `fld_v`, `fld_w`, `fld_p`, `fld_psi`, `fld_fx_old`, `fld_fy_old`,
`fld_fz_old`, `fld_alphac`, plus the particle state via `loadpart('r')`.

`kappa`, `normx/y/z` and `phi` are **not** checkpointed — they are recomputed
from `psi` immediately after loading. `Fstot`/`Fstot_old` are **not** restored
(they restart from zero), which introduces a one-step transient in `F_cap`
after every restart.
