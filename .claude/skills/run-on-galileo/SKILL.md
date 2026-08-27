---
name: run-on-galileo
description: Stage and launch a production run on CINECA Galileo100 — the $WORK/$SCRATCH split, sizing the grid to a target particle resolution, choosing dims/ntasks, and validating a configuration before burning wall time. Use when the user wants to run a case on the cluster, asks "set up a run", picks a resolution in cells per diameter, or hits a decomposition/resource error at job start.
---

# Running a case on Galileo100

Configuring `input.nml` itself is [`new-simulation-case`](../new-simulation-case/SKILL.md);
compiling is [`build-solver`](../build-solver/SKILL.md). This skill is about
where things live on the cluster and how to get a job onto the queue safely.

## Storage layout — code in `$WORK`, output in `$SCRATCH`

| area | path | quota | lifetime |
|---|---|---|---|
| `$HOME` | `/g100/home/userexternal/afavero0` | 50 GB | permanent, backed up |
| `$WORK` | `/g100_work/IscrC_TP-PBR` | **1 TB** | **dies with the project** |
| `$SCRATCH` | `/g100_scratch/userexternal/afavero0` | none | **40-day cleanup** |

```
$WORK/CaNS-Fizzy/          the git clone; pull and build HERE
$WORK/bin/stage-run        deployed copy of utils/cluster/stage-run
$SCRATCH/runs/<case>/      cans + input.nml + jobfile.slurm + data/ + PROVENANCE
```

The scratch cleanup is **not** a threat to the code — the code is in `myfork`,
and `git clone --recursive` restores it in minutes. It threatens *output*, and
output is exactly what cannot live in `$WORK` because of the 1 TB quota. Hence
the split.

**`$WORK` is not an archive.** It is bound to the project grant
(`IscrC_TP-PBR`, ends **2026-12-02** — check with `saldo -b`). Anything that
must outlive the grant has to leave the cluster entirely; the user's habit is
`rsync … data.g100.cineca.it → ~/Dottorato/CINECA/`.

## Staging a run

```bash
$WORK/bin/stage-run <case-name> [input.nml]
```

Defaults to `src/input.nml`; pass an `examples/Three_Phase/*/input.nml` to start
from a real case. It creates the run directory, copies the **binary** (not a
symlink, so a later rebuild never rewrites the provenance of an old run), writes
a `PROVENANCE` file (commit, branch, dirty flag, `build.conf`), and generates a
`jobfile.slurm`. It refuses to overwrite an existing run directory.

If `$WORK/bin/stage-run` is missing or stale, redeploy it from the repo:
`cp $WORK/CaNS-Fizzy/utils/cluster/stage-run $WORK/bin/ && chmod +x $WORK/bin/stage-run`.

## Sizing the grid to a particle resolution

The user asks for resolution as **cells per sphere diameter**, not as `ng`.
The grid is uniform and isotropic (`gtype = 1, gr = 0.`), so:

```
dx = l(1)/ng(1)          D = 2*radius          D/dx = cells per diameter
ng(i) = l(i) * (D_cells / D) = l(i) * D_cells / (2*radius)
```

Check the example you started from — it is often coarser than assumed.
`examples/Three_Phase/Bouncing_Sphere` ships at `ng = 64,64,48`, `l = 16,16,12`,
`radius = 1`, i.e. `dx = 0.25` and **D = 8 cells**. Going to D = 32 is a 4×
refinement in every direction: `ng = 256,256,192`.

**Refining costs more than cells.** `dt` falls roughly linearly with `dx`, so the
same `nstep` covers proportionally less physical time. Say so explicitly and let
the user decide `nstep` — do not silently change it. Measure the real `dt` from
a short validation run (below) rather than predicting it.

## Choosing `ntasks` and `dims`

`dims(1)*dims(2)` must equal `ntasks`, and **`ng(1)` must be divisible by
`dims(1)`, `ng(2)` by `dims(2)`**. With `PENCIL_AXIS=3` only x and y are
decomposed; `ng(3)` is not.

This bites: for `ng = 256,256,192`, `ntasks = 48` has no valid factorisation
(256 is not divisible by 6 or 12). Workable: 16 (`4,4`), 32 (`8,4`), 64 (`8,8`).

Do not leave `dims(1:2) = 0,0` for particle runs — set it explicitly so the
decomposition is reproducible and the particle pencil ownership is predictable.

Memory is rarely the constraint: ~25 arrays × 8 B × cells. 12.6 M cells ≈ 2.5 GB
total. Nodes have 48 cores / 375 GB.

## Validate before submitting

A 24 h job that dies in the first second on a `dims` mismatch is a wasted day.
Copy the staged directory, shrink it, and run it on the debug partition:

```bash
cd $SCRATCH/runs && cp -r <case> _check && cd _check
sed -i "s|^nstep = .*|nstep = 3, time_max = 10., tw_max = 0.05|" input.nml
sed -i "s|iout2d = [0-9]*|iout2d = 1|" input.nml
sed -i "s|g100_usr_prod|g100_usr_dbg|; s|--time=24:00:00|--time=00:15:00|" jobfile.slurm
sbatch jobfile.slurm
```

Then confirm, and **delete `_check`** so the staged directory stays pristine:

- job `COMPLETED`, exit `0:0`
- 2D slice size == `ng(1)*ng(3)*8` bytes for the default `y = l(2)/2` plane
- `data/log_visu_2d_slice_1.out` exists (XDMF needs it)
- `forces_data.csv` header is `F_cap_ibm,F_ibm,F_inertia,F_w,F_bouy,F_cap,ep_z,ep_w`
- read `dt` from `Time = …` at step 3 → extrapolate wall time for the real `nstep`

Never run the solver on a login node.

## Gotchas

- **`tw_max` is in HOURS.** The examples ship `tw_max = 0.1` — six minutes. It is
  inert only because `stop_type(1:3) = T,F,F` selects `nstep`. Flipping the third
  flag silently caps the run.
- **Never rebuild while a job is running out of that directory.** Relinking
  `run/cans` under a live process fails with `ETXTBSY` or corrupts it. This is the
  main reason the build lives in `$WORK` and runs are copies on `$SCRATCH`.
- **`login.g100.cineca.it` round-robins** over `130.186.24.{4,5,6}`. All three
  share one ECDSA host key but each has its own ED25519 key, so ssh can report
  "REMOTE HOST IDENTIFICATION HAS CHANGED" purely from landing on a different
  node. Verify with `ssh-keyscan` before touching `known_hosts` — pinning one IP
  with `-o HostKeyAlias=login.g100.cineca.it` keeps strict checking intact.
- Job budget: `saldo -b`. Wall time on `g100_usr_prod` maxes at 24 h;
  `g100_usr_dbg` is for the validation runs.
