## Synopsis

This repository is a **research fork of [CaNS-Fizzy](https://github.com/CaNS-World/CaNS-Fizzy)**
extending it from two-phase flow to **three-phase flow: two immiscible incompressible
fluids plus finite-size rigid spherical particles**, coupled through an
**extended contact-line model** that imposes a prescribed contact angle where the
fluid–fluid interface meets the particle surface.

The base solver is unchanged in spirit: a massively-parallel, second-order
finite-difference Navier–Stokes solver on a staggered 3D Cartesian grid, using a
pressure-splitting technique so that the fast FFT-based constant-coefficient
Poisson solver of [CaNS](https://github.com/CaNS-World/CaNS) can be used even at
high density contrasts. On top of it this fork adds:

* a **finite-size particle module** — rigid spheres resolved on the Eulerian grid
  via a diffuse solid indicator and a direct-forcing immersed-boundary method,
  with Newton–Euler rigid-body dynamics, soft-sphere DEM collisions and a
  lubrication correction;
* an **extended contact-line model** — the interface is relaxed each timestep
  along an extension velocity built from the prescribed contact angle, and the
  resulting capillary force on the particle is integrated along the numerical
  contact line.

The interface between the two fluid phases may be captured with either the
THINC/QQ volume-of-fluid method (the default here) or the Accurate Conservative
Diffuse Interface (ACDI) method.

**References**

P. Costa. *A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows.* *Computers & Mathematics with Applications* 76: 1853--1862 (2018). [doi:10.1016/j.camwa.2018.07.034](https://doi.org/10.1016/j.camwa.2018.07.034) [[arXiv preprint]](https://arxiv.org/abs/1802.10323)

G. Frantzis & D. Grigoriadis. *An efficient method for two-fluid incompressible flows appropriate for the immersed boundary method.* *Journal of Computational Physics* 376 (2019): 28-53. [doi.org/10.1016/j.jcp.2018.09.035](https://doi.org/10.1016/j.jcp.2018.09.035).

P. Cifani. *Analysis of a constant-coefficient pressure equation method for fast computations of two-phase flows at high density ratios.* *Journal of Computational Physics* 398 (2019): 108904. [doi:10.1016/j.jcp.2019.108904](https://doi.org/10.1016/j.jcp.2019.108904).

S. Jain. *Accurate conservative phase-field method for simulation of two-phase flows.* *Journal of Computational Physics* 469 (2022): 111529. [doi.org/10.1016/j.jcp.2022.111529](https://doi.org/10.1016/j.jcp.2022.111529)

Bin & Xiao. *Toward efficient and accurate interface capturing on arbitrary hybrid unstructured grids: The THINC method with quadratic surface representation and Gaussian quadrature* *Journal of Computational Physics* 349 (2017): 415-440. [doi.org/10.1016/j.jcp.2017.08.028](https://doi.org/10.1016/j.jcp.2017.08.028)

## What the code currently does

* Solves the incompressible Navier–Stokes equations for two immiscible Newtonian
  fluid phases in a one-fluid formulation on a staggered Cartesian grid.
* Advances time with a **second-order Adams–Bashforth scheme** with a
  variable-timestep correction (the three-step Runge–Kutta driver is retained in
  `src/main.f90` but commented out).
* Captures the fluid–fluid interface with the **THINC/QQ volume-of-fluid method**
  (`INTERFACE_CAPTURING_VOF=1`, the current default) or with ACDI.
* Computes interface normals and curvature from a reconstructed signed-distance
  field (`SDF_NORMALS=1`).
* Represents rigid spheres with a **diffuse solid indicator** `alphac`, rebuilt
  geometrically every timestep from a hyperbolic-tangent profile of the signed
  distance to the sphere surface, over a shell of `eps_sol` cells.
* Couples particles to the fluid with an **Eulerian direct-forcing immersed
  boundary method**, driving the velocity toward rigid-body motion inside the
  particle and accumulating the reaction force and torque.
* Integrates **Newton–Euler rigid-body dynamics** for each particle, including
  the fluid momentum and buoyancy integrals over the sphere volume, soft-sphere
  DEM collisions (particle–particle and particle–wall, with normal/tangential
  springs and dashpots and a Coulomb friction cap) and a lubrication correction
  with Stokes amplification factors.
* Enforces a **prescribed contact angle** `theta` at the fluid–fluid–solid
  contact line by relaxing the volume fraction along an extension velocity in a
  short pseudo-time loop each timestep (`src/extend.f90`).
* Integrates the **capillary force** acting on the particle along the numerical
  contact line and reduces it across ranks (`src/rotnorm.f90`). This force is
  currently reported and logged; feeding it back into the particle momentum
  equation is available but commented out in `src/prt_intgr_nwtn_eulr.f90`.
* Reports a dimensionless-parameter banner at startup (Archimedes, Eötvös,
  Ohnesorge, Morton numbers, capillary length and the capillary timestep limit),
  and writes a per-step force breakdown to `forces_data.csv`.

## Features

 * MPI parallelization with 2D pencil decomposition
 * FFTW guru interface / cuFFT used for computing multi-dimensional vectors of 1D transforms
 * The right type of transformation (Fourier, cosine, sine, etc) is automatically determined form the input file
 * [cuDecomp](https://github.com/NVIDIA/cuDecomp) pencil decomposition library for _hardware-adaptive_ distributed memory calculations on _many GPUs_
 * [2DECOMP&FFT](https://github.com/xcompact3d/2decomp-fft) library used for performing global data transpositions on CPUs and some of the data I/O
 * GPU acceleration using OpenACC directives (core solver; the particle module has not been validated on GPU)
 * A different canonical flow can be simulated just by changing the input files
 * Mass transport-consistent discretization of advection terms that enable simulations at high density contrasts
 * Two options for interface-capturing algorithms: ACDI (diffuse interface method) and THINC/QQ (volume-of-fluid method)
 * Finite-size resolved rigid particles with collisions and lubrication
 * Prescribed contact angle at the particle surface, with contact-line force diagnostics

## Repository layout

A full description is in [`ARCHITECTURE.md`](ARCHITECTURE.md).

```
src/               solver sources (flat; the Makefile globs it)
  main.f90         driver and timeloop
  *.f90            core two-phase Navier-Stokes solver
  prt_*.f90        particle module
  extend.f90       contact-line relaxation
  rotnorm.f90      capillary force at the contact line
examples/
  Three_Phase/     particle + two-fluid cases (this fork)
  Two_Phase/       upstream Fizzy validation cases
  _CaNS-example-files/  single-phase CaNS heritage cases
PostPrt/           standalone post-processing programs
utils/             visualisation and binary-reading helpers
docs/              INFO_COMPILING / INFO_INPUT / INFO_VISU
```

## Method

The two-phase flow is described by a one-fluid formulation and solved with a
second-order finite-difference incremental pressure-correction scheme on a
staggered grid. The interface between the two fluid phases is advected with the
THINC/QQ volume-of-fluid method (or represented as a diffuse interface of
specified thickness under ACDI). A pressure-splitting technique converts the
variable-coefficients Poisson equation into a constant-coefficients one, so the
fast Poisson solver of CaNS applies.

Rigid particles are superimposed on this formulation through a second indicator
field. The solid indicator is built geometrically each timestep, so no
interface-capturing error accumulates on the particle surface as it moves. The
fluid inside the particle is driven toward rigid-body motion by direct forcing,
and the reaction force closes the Newton–Euler equations for the particle.

At the three-phase contact line, the fluid–fluid interface is relaxed toward the
prescribed contact angle by advecting the volume fraction along an extension
velocity constructed from the solid normal, the interface normal and the contact
angle. The capillary force is then obtained by integrating
`sigma * |grad(psi) x grad(alphac)|` along the tangent to the interface over the
contact-line cells.

## Usage

### Downloading

The external pencil decomposition libraries are Git Submodules, so clone
recursively:
```bash
git clone --recursive <repository-url>
```
If the repository has already been cloned without them (i.e. `dependencies/cuDecomp`
and `dependencies/2decomp-fft` are empty), run:
```bash
git submodule update --init --recursive
```
A missing submodule shows up at build time as `No rule to make target 'clean'`,
which does not otherwise hint at the cause.

### Compilation

#### Prerequisites

 * MPI
 * FFTW3/cuFFT library for CPU/GPU runs
 * The `nvfortran` compiler (for GPU runs)
 * NCCL and NVSHMEM (optional, may be exploited by the cuDecomp library)
 * HYPRE library in case the variable-coefficients Poisson equation is solved without the pressure splitting technique (only available for CPU runs)

#### In short

```bash
make libs && make
```
compiles the 2DECOMP&FFT/cuDecomp libraries and then the solver. The executable
is placed in `run/`.

#### Detailed instructions

The `Makefile` in the root directory is expected to work out-of-the-box on most
systems. `build.conf` selects the Fortran compiler (MPI wrapper), a build profile
(production vs debugging) and the pre-processing options. General options:

 * `DEBUG`                    : performs some basic checks for debugging purposes
 * `TIMING`                   : wall-clock time per time step is computed
 * `PENCIL_AXIS`              : sets the default pencil direction, one of [1,2,3] for [X,Y,Z]-aligned pencils
 * `SINGLE_PRECISION`         : calculation will be carried out in single precision (the default precision is double)
 * `GPU`                      : enable GPU-accelerated runs

Options specific to this fork:

 * `INTERFACE_CAPTURING_VOF`  : `1` uses THINC/QQ volume-of-fluid, `0` uses ACDI
 * `SDF_NORMALS`              : `1` computes normals/curvature from the reconstructed signed-distance field
 * `PARTICLE`                 : `1` activates the particle module
 * `EULER`                    : `1` selects the Eulerian immersed-boundary coupling (`0` selects the Lagrangian forcing-point path, which is not maintained)

**`build.conf` changes are not tracked as build dependencies — run
`make clean && make` after editing it.**

See [`INFO_COMPILING.md`](docs/INFO_COMPILING.md) for the full list of options.

### Input file

`input.nml` sets the physical and computational parameters. It contains the
`&dns`, `&scalar`, `&two_fluid` and `&cudecomp` groups read by the core solver,
plus the `&particle`, `&collision_parameters` and `&particle_euler` groups read
by the particle module. Examples for several configurations are in `examples/`;
see [`INFO_INPUT.md`](docs/INFO_INPUT.md) for the base parameters.

Files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/` set which data are
written in 1-, 2- and 3-dimensional output files, and `prt_out.h90` the particle
output. *The code must be recompiled after editing these files.*

### Example cases

Under `examples/Three_Phase/`:

| case | description |
|---|---|
| `Bouncing_Sphere` | dense sphere impacting a flat liquid film and rebounding |
| `Sinking_Sphere`  | same configuration at lower surface tension, the sphere penetrates |
| `Wall_Collision`  | sphere impacting a wall in a single fluid |
| `Head_On`         | sphere fired head-on at a droplet in zero gravity |
| `Sessile_Drop`    | droplet resting on a sphere at a prescribed contact angle |
| `Particle_Capture`| 500 particles interacting with a rising bubble |

Cases whose `inipsi` refers to spheres, cylinders or films (`bub3`, `drp3`, …)
additionally require a `spheres.in` file.

### Running the code

Run the executable with `mpirun` using a number of tasks consistent with `dims`
in `input.nml`. Data are written by default into `data/`, which must exist where
the executable is run (by default the `run/` folder). Note that the particle
module requires `radius + offset` to be smaller than the subdomain size in `x`
and `y`, so a large particle constrains the maximum number of ranks.

### Output

| file | content |
|---|---|
| `data/fld_*.bin` | checkpoint fields (velocity, pressure, `psi`, `alphac`, forcing history) |
| `data/*_fld_*.bin` + `log_visu_3d.out` | 3D field series for visualisation |
| `data/*_slice_fld*.bin` + `log_visu_2d_*.out` | 2D slices |
| `data/time.out` | time step number, `dt`, physical time |
| `forces_data.csv` | per-step force breakdown on the particle: `F_drag,F_ibm,F_inertia,F_w,F_bouy,F_cap,ep_z,ep_w` |

### Visualizing field data

See [`INFO_VISU.md`](docs/INFO_VISU.md). Post-processing programs for particle
trajectories and surface force decomposition are in `PostPrt/`.

## Contributing

This is a research fork. For contributions to the upstream solver, please refer
to the [CaNS-Fizzy](https://github.com/CaNS-World/CaNS-Fizzy) repository.

## Final notes

Please read the `ACKNOWLEDGEMENTS`, `LICENSE` files.
