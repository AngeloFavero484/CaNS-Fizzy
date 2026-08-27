---
name: build-solver
description: Compile CaNS-Fizzy (this three-phase particle fork) — locally or on the CINECA cluster. Use when the user asks to build, compile, recompile, or reports a make/compiler error, missing module files, "No rule to make target", or link failures.
---

# Building the solver

## The short version

```bash
make libs && make        # first time, or after changing build.conf
make                     # incremental
make clean && make       # after ANY build.conf change (see below)
```

The executable lands in `run/cans`, together with copies of `src/dns.in` and
`src/cudecomp.in`.

## Before anything else: check the submodules

`dependencies/2decomp-fft` and `dependencies/cuDecomp` are **git submodules**.
A plain `git clone` leaves them empty, and the failure looks unrelated:

```
make[2]: *** Nessuna regola per generare l'obiettivo "clean".  Arresto.
make[2]: *** No rule to make target 'clean'.  Stop.
make[1]: *** [dependencies/external.mk:15: libsclean] Error 2
```

Diagnose and fix:

```bash
ls dependencies/2decomp-fft          # empty or missing Makefile => submodules not init'd
git submodule update --init --recursive
```

For a fresh clone, prefer `git clone --recursive`.

## `make clean` is required after editing `build.conf`

The Makefile does **not** treat `-D` defines as dependencies. Changing
`INTERFACE_CAPTURING_VOF`, `PARTICLE`, `EULER`, `SDF_NORMALS`,
`SINGLE_PRECISION`, `CONSTANT_COEFFS_POISSON` etc. and running a plain `make`
produces an inconsistently-compiled binary that may link and then behave
nonsensically. Always `make clean && make`.

Same applies after editing any `.h90` include file — those are pulled into
`main.f90` by the preprocessor.

## Reading `build.conf`

Current settings and what they mean are documented in
`.claude/references/preprocessor-flags.md`. The ones that change physics:

| setting | current | effect |
|---|---|---|
| `INTERFACE_CAPTURING_VOF` | `1` | THINC/QQ VOF active; **`acdi.f90` is dead code** |
| `SDF_NORMALS` | `1` | normals from the SDF `phi`, not `psi` |
| `PARTICLE` | `1` | particle module on — **cannot be 0** in this fork (`extend.f90`/`rotnorm.f90` need `alphac`) |
| `EULER` | `1` | Eulerian IBM; the Lagrangian path is dead *and* commented out |
| `CONSTANT_COEFFS_POISSON` | `1` | FFT Poisson solver (fast); `0` needs HYPRE |
| `GPU` | `0` | CPU build; `1` requires `FCOMP=NVIDIA` |

## Prerequisites

- MPI (compiler wrapper — the Makefile uses `mpifort`/`mpif90` via `FCOMP`)
- FFTW3 (CPU) or cuFFT (GPU)
- `nvfortran` only for `GPU=1`
- HYPRE only for `CONSTANT_COEFFS_POISSON=0`

## On the CINECA cluster (Galileo100)

Build in **`$WORK/CaNS-Fizzy`** (`/g100_work/IscrC_TP-PBR/CaNS-Fizzy`), never in
a run directory — see [`run-on-galileo`](../run-on-galileo/SKILL.md). Runs are
staged copies under `$SCRATCH/runs/`, so rebuilding never disturbs a live job.

```bash
cd $WORK/CaNS-Fizzy && module purge && module load openmpi fftw && make libs && make
```

`make libs` is needed once per fresh clone; without it the build stops at
`Cannot open module file 'decomp_2d.mod'`, which does not mention the missing
library. Modules must be loaded before `make`; `openmpi` + `fftw` is verified
working with `FCOMP=GNU`.
If a build fails on the cluster but works locally, **suspect the module
environment first**, not the code.

To ask the user to run something interactive there, tell them to type
`! <command>` in the prompt so the output lands in the conversation.

## Debug builds

Set exactly one profile in `build.conf`:

```
FFLAGS_OPT=0
FFLAGS_DEBUG=1
```

Then `make clean && make`. Two things to expect, both normal:

1. **FPE traps that do not fire at `-O3`.** `-ffpe-trap=invalid` will catch the
   `0/0` cases that the `+ epsilon(1._rp)` guards in `extend.f90` and
   `rotnorm.f90` mitigate but do not eliminate. This is the debug build being
   correct.
2. **`-std=f2018` warnings/errors** in the `prt_*` files, which were ported from
   an older code and use constructs `-O3` accepts silently.

## Verifying a build

```bash
cd run && mpirun -n 1 ./cans
```

A good start prints, in order: the compiler/MPI banner (`DEBUG=1`), the
dimensionless-parameter block (`Ar`, `Eo`, `Oh`, `Mo`, capillary length and
`dt_cap`), grid setup, `*** Particle initialization ***`, the solid volume
fraction, `Fstot`, then `dt_cfl = ... dt = ...` and the first timestep.

If it aborts immediately with *"Radius spheres larger than x-dimension of
processes"*, that is `prt_InitMemo`'s size check — **use fewer MPI ranks** or
enlarge the domain; do not edit the check.
