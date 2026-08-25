# Reference: data structures, array conventions, MPI layout

---

## Real kind

Everything uses `rp` from `mod_types`:

```fortran
use mod_types            ! brings in rp, sp, dp, MPI_REAL_RP
real(rp) :: x
x = 0.5_rp               ! CORRECT
x = 0.5d0                ! WRONG - breaks SINGLE_PRECISION builds
```

`rp = dp` normally; `rp = sp` under `-D_SINGLE_PRECISION`. `MPI_REAL_RP` is the
matching MPI datatype — always use it, never `MPI_DOUBLE_PRECISION` directly.

---

## Field arrays

### Core solver fields — one halo cell

```fortran
allocate(u(0:n(1)+1, 0:n(2)+1, 0:n(3)+1))
```
`u, v, w, p, pp, pn, po, psi, psio, phi, kappa, normx, normy, normz,
psiflx_x/y/z, fx_old/fy_old/fz_old, u_ext/v_ext/w_ext, s`

Interior is `1:n(1), 1:n(2), 1:n(3)`. Index `0` and `n+1` are halo, filled by
`boundp` (scalars) or `bounduvw` (velocity).

### Particle fields — `nh_wide` halo

```fortran
allocate(alphac(1-nh_wide:n(1)+nh_wide, ...))
```
`alphac, norm_partx/y/z, uphase/vphase/wphase` (and `uf/vf/wf` in the
non-EULER build).

`nh_wide = 1` (`param.f90:24`), so numerically identical bounds to the core
fields — but they are declared as `dimension(1-nh_wide:,1-nh_wide:,1-nh_wide:)`
in every `prt_*` routine signature. **Do not mix the two declaration styles when
passing arrays across the boundary**; assumed-shape dummies with different lower
bounds will silently reindex.

### Staggering

Standard MAC arrangement:
- `p`, `psi`, `alphac`, `kappa` at cell **centres**
- `u` at `x`-faces, `v` at `y`-faces, `w` at `z`-faces
- index `i` for `u` means the face between cells `i` and `i+1`

`prt_initvof.f90` is the clearest example of handling this — it computes
`coorx_cent = boundleftmyid + (i-0.5)*dx` and
`coorx_stag = boundleftmyid + i*dx`, then digitises `alphac` **separately** for
each of the three staggered positions (`alpha_eulx`, `alpha_euly`, `alpha_eulz`).

### Non-uniform z

`dzc`/`dzf` (spacing at centres/faces) and `dzci`/`dzfi` (inverses) are 1-D
arrays over `0:n(3)+1`. The `_g` suffixed versions are the **global** z-grid.
`dl(1:2)`/`dli(1:2)` are scalars because x and y are always uniform.

**`dli` and `dzci`/`dzfi` are INVERSE spacings.** A derivative is
`(f(i) - f(i-1)) * dli(1)`, never `/ dli(1)`.
(See the caveat in `contact-line-model.md` — `extend.f90` gets this wrong.)

---

## Particle derived types (`prt_common.f90`)

Three parallel arrays, all `dimension(1:npmax)`:

| array | type | meaning |
|---|---|---|
| `ep(:)` | `particle` | **current** state |
| `op(:)` | `particle_old` | **previous** step's state (for the AB/trapezoidal update) |
| `tp(:)` | `particle_old` | temporary, used during master/slave migration |
| `rkp(:)` | `particle_sumrk` | RK3 substep accumulator — **vestigial** under Adams–Bashforth |

### `type particle` field ordering is load-bearing

```fortran
type particle
   real(rp) :: x,y,z,theta,phi, &
               u,v,w, &
               omx,omy,omz,omtheta, &
               intu,intv,intw, &
               intomx,intomy,intomz, &
               intrhox,intrhoy,intrhoz, &
               colfx,colfy,colfz, &
               coltx,colty,coltz            ! 24 contiguous reals
   real(rp), dimension(nqmax) :: dx,dy,dz, ...
   ...
end type
```

The comment block above the type in `prt_common.f90` states the rule:

1. reals that must be communicated on a master change — **contiguous, first**
2. then integers that must be communicated — contiguous
3. then reals that need not be communicated
4. then integers that need not be communicated

`prt_param.f90` encodes the counts:
```fortran
integer, parameter :: send_real = 24 + 10*nqmax
integer, parameter :: send_int  = 1  + nqmax
```

**Adding a field to `type particle` in the wrong position silently corrupts MPI
transfers.** Add non-communicated fields at the end; if a field must be
communicated, insert it inside the contiguous real block *and* bump `send_real`.

`nqmax = 15` is a compile-time cap on simultaneous collision partners + walls per
particle. It sizes fixed-shape components, so it cannot come from the runtime `np`.

### Key `particle` members

| member | meaning |
|---|---|
| `x,y,z` | centre position |
| `u,v,w` | translational velocity |
| `omx,omy,omz` | angular velocity; `theta,phi,omtheta` track orientation |
| `intu,intv,intw` | ∫ fluid momentum over the sphere volume (`intgr_over_sphere`, `cas=1/2/3`) |
| `intomx/y/z` | ∫ angular momentum over the sphere |
| `intrhox/y/z` | ∫ density over the sphere → buoyancy |
| `fxltot,fyltot,fzltot` | total IBM reaction force |
| `torqxltot,...` | total IBM torque |
| `colfx,colfy,colfz` / `coltx,...` | collision force / torque |
| `fcapx,fcapy,fcapz` | capillary force slot (**currently unused in the update**) |
| `vol`, `mominert`, `ratiorho` | per-particle derived constants |
| `mslv` | **master/slave flag** — see below |
| `nb(1:8)` | which of the 8 neighbours this particle overlaps |

---

## MPI decomposition and the master/slave protocol

### Decomposition

2-D pencil in x–y; **z is never split**. `dims(1:2)` from `input.nml`
(`0,0` = MPI chooses). Set up in `initmpi.f90`, giving each rank
`lo(1:3)`, `hi(1:3)`, `n(1:3)`, `is_bound(0:1,1:3)`, and neighbour ranks.

The particle module builds its **own** communicator in `prt_InitMemo`:

```fortran
call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims(1:2), (/.true.,.true./), .false., prt_comm_cart, ierr)
```

Note it is created **periodic in both directions regardless of `cbcvel`** — the
periodic-image corrections in the `prt_*` routines rely on that, and handle
non-periodic physical BCs separately via `is_bound`.

`neighbor(0:8)`: index 0 is self, 1–8 are
right, right-front, front, left-front, left, left-back, back, right-back.
`boundleftmyid`/`boundfrontmyid` are this rank's physical x/y origin.

### Master / slave

A particle larger than one cell straddles subdomains. Exactly one rank owns it:

| `ep(p)%mslv` | meaning |
|---|---|
| `> 0` | this rank is **master** of particle `mslv` |
| `< 0` | this rank is **slave** of particle `-mslv` |
| `= 0` | slot unused |

`ep(p)%nb(1:8)` records which neighbours the sphere reaches into.

**Every particle routine opens with the same ~150-line boilerplate:**

1. pack `ep(p)%{x,y,z,...}` into a local `type pneighbor` array `anb(0:8,1:npmax)`
2. masters `MPI_ISEND` to the slave ranks in `nb`; slaves `MPI_IRECV`
   (with a self-send special case when `neighbor(nbsend) == myid`, which happens
   under periodicity with few ranks)
3. `MPI_WAITALL`
4. **periodic-image correction**: shift `anb(nbrecv,p)%x/y` by `±l(1)`/`±l(2)`
   so the received centre is expressed in the *receiver's* coordinates
5. compute the local index box `ilow..ihigh, jlow..jhigh, klow..khigh`,
   clamp it to `1..n`, and loop

This is duplicated verbatim in `prt_initeul`, `prt_initvof`,
`prt_intgr_over_sphere`, `prt_eulint`, `prt_phase_indicator`, and
`prt_intgr_nwtn_eulr`. **When changing the protocol, change all of them.**

`npmax` is the per-rank slot count:
```fortran
npmax = nint(min(1.*np, max(1., 10.*np/(1.*dims(1)*dims(2)))))
```
i.e. ~10× the average share, capped at `np`.

### Hard size constraint

`prt_InitMemo` aborts if
```
radius + offset > l(1)/dims(1)   or   radius + offset > l(2)/dims(2)
```
with `offset = 1.01/dli(1)` under `_EULER`. A particle may span at most one
neighbour in each direction. **Fix by using fewer ranks**, not by editing the check.

---

## The two indicator fields, side by side

| | `psi` | `alphac` |
|---|---|---|
| meaning | fluid-1 volume fraction | solid indicator |
| range | `[0,1]` | `[0,1]` |
| built by | `vof_thinc_qq.f90` (transport) | `prt_digitiser.f90` (geometric, rebuilt each step) |
| shape | evolves with the flow | always a sphere at `ep(p)%{x,y,z}` |
| thickness | THINC `beta` | `eps_sol` cells (`tanh` profile) |
| checkpointed | yes (`fexts(5)`) | yes (`fexts(9)`) |
| halo | `0:n+1` | `1-nh_wide:n+nh_wide` |

`alphac` is **recomputed from scratch** every step by `initeul` — it is never
advected. That is why the particle can move without any interface-capturing
error accumulating on the solid.

---

## Output field naming

`write_visu_3d(datadir, 'xxx_fld_NNNNNNN.bin', 'log_visu_3d.out', 'VarName', ...)`

Current 3-D outputs include `alx_fld_*` (`Alpha_C`) — the fork replaced the
upstream `cur_fld_*` (`Kappa`) with it in `out3d.h90`. 2-D adds
`alx_slice_fld*`. The `log_visu_*.out` files are consumed by
`utils/visualize_fields/` to generate XDMF.
