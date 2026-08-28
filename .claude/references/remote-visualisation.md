# Remote visualisation — ParaView on CINECA Galileo100

How to look at simulation output **without pulling it off the cluster**. The
production runs are 384x384x960 (~142 M cells); a single 3D field is ~1.1 GB in
double precision and `run/data` is already 14 GB. `rsync`-ing that to the laptop
for every look is the thing this setup replaces.

## What G100 actually has (checked 2026-08-27)

| module | what it gives | notes |
|---|---|---|
| `paraview/5.9.1-osmesa` | `pvserver`, `pvbatch`, `pvpython` | **no GUI binary** — OSMesa build, software off-screen rendering |
| `rcm/02` | TurboVNC remote desktop | *not useful here*: there is no ParaView GUI installed on G100 to run inside it |

So the only workable route is **client–server**: the GUI runs on the laptop, the
data stays on `/g100_scratch` and is read and rendered by `pvserver` on a
compute node, which ships rendered images back over an SSH tunnel.

`pvpython` from this module bundles **numpy 1.19.2**, which matters because the
system `python/*` modules on G100 have no numpy.

## Version lock

ParaView refuses a client/server pair whose versions differ. The cluster is
pinned at **5.9.1**, so the laptop needs a 5.9.1 client:

    ~/opt/ParaView-5.9.1/bin/paraview

The distro ParaView (`/usr/bin/paraview`, 5.11.0) is untouched and still fine
for local work on downloaded data — it just cannot talk to G100.

## Everyday use

    g100-paraview

That one command submits the server job, waits for the allocation, opens the
tunnel, and starts the local GUI already connected. On exit it cancels the job
and closes the tunnel, so nothing keeps burning budget.

Options:

    g100-paraview -n 16 -m 200GB -t 08:00:00     # bigger session
    g100-paraview -p g100_usr_bmem -m 350GB      # very large 3D fields
    g100-paraview --no-client                    # job + tunnel only

Defaults: `g100_usr_interactive`, 1 node, 8 ranks, 120 GB, 4 h, account
`IscrC_TP-PBR`. That partition allows at most 2 nodes and 8 h, which is the
right shape for interactive work; use `g100_usr_prod` for longer sessions.

### Ports are negotiated, not fixed

SLURM packs interactive jobs onto the same node, so a hard-coded 11111 collides
with another `pvserver` — and **pvserver still exits 0 when the bind fails**, so
the job just ends a couple of seconds after starting with only this in the
`.err` file:

    ERR| vtkServerSocket: Socket error in call to bind. Address already in use.

The job script therefore probes upward from 11111 for a free port, confirms
`pvserver` is really listening on it, and publishes the result in
`~/pvserver/job-<jobid>.info`. The launcher reads that back and tunnels to the
port that was actually taken. It also steps the *local* port up if 11111 is
already in use on the laptop. Nothing needs to be chosen by hand.

### Cleanup

Quitting the client (or Ctrl-C on `--no-client`) fires a trap that closes the
tunnel and `scancel`s the job — verified, so a forgotten session does not sit
on the allocation. The walltime is the backstop. `pvserver` also exits by
itself when the client disconnects.

Once the GUI is up, **File > Open browses the cluster filesystem**, not the
laptop. Data lives under `/g100_scratch/userexternal/afavero0/`.

## The XDMF step (required)

ParaView cannot open the raw `.bin` fields. It needs the `.xmf` metadata that
describes the grid and the time series, generated from the `log_visu_*.out`
file the solver writes. Do this **on the cluster**, once per run:

    ~/pvserver/make-xdmf.sh <data-dir> [logfile] [grid-suffix] [out.xmf]

For 3D fields (`iout3d` > 0) the defaults are right:

    ~/pvserver/make-xdmf.sh /g100_scratch/userexternal/afavero0/<case>/run/data

For the 2D slices (`iout2d`, the current Bouncing_Sphere setup):

    ~/pvserver/make-xdmf.sh /g100_scratch/.../run/data \
        log_visu_2d_slice_1.out 2d viewfld_DNS_2d.xmf

It wraps `utils/visualize_fields/gen_xdmf_easy/write_xdmf.py`, answering its
three prompts, and runs it under `pvpython` for numpy. Then open the resulting
`.xmf` from the connected client.

Note `iprecision = 8` in `write_xdmf.py` — it must match the solver's output
precision or the fields come out as garbage.

## Files involved

| where | path | role |
|---|---|---|
| laptop | `~/bin/g100-paraview` | launcher: sbatch + tunnel + client |
| laptop | `~/opt/ParaView-5.9.1/` | version-matched client |
| G100 | `~/pvserver/pvserver.slurm` | the server job |
| G100 | `~/pvserver/make-xdmf.sh` | XDMF generation |
| G100 | `~/pvserver/job-<jobid>.info` | node/port the server actually got |

## Authentication

CINECA certificates last 12 h. When `g100-paraview` reports it cannot ssh:

    step ssh login 'angeloraimondo.favero@studenti.unipd.it' --provisioner cineca-hpc

## Batch rendering, no GUI

For producing figures/animations unattended, `pvbatch` on the cluster runs a
ParaView Python script with the same OSMesa rendering and writes PNGs:

    module load paraview/5.9.1-osmesa
    mpiexec -n 8 pvbatch myscript.py

Trace a session in the GUI (Tools > Start Trace) to get a script to adapt.
