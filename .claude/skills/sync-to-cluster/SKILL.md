---
name: sync-to-cluster
description: Commit and push work to the user's fork, or move changes between the local machine and the CINECA supercomputer. Use when the user says "commit and push", "push my changes", "sync to the cluster", or hits git/SSH/submodule problems on the HPC side.
---

# Syncing work between local and the cluster

## Remotes — get this right

| remote | URL | role |
|---|---|---|
| `origin` | `https://github.com/CaNS-World/CaNS-Fizzy` | **upstream**, read-only |
| `myfork` | `git@github.com:AngeloFavero484/CaNS-Fizzy.git` | **the user's fork** |

**Always `git push myfork main`.** A bare `git push` targets `origin`, which is
the upstream project — wrong, and it will be rejected anyway.

The local `main` sits many commits ahead of `origin/main`. That is the permanent,
intended state of this fork, not drift to be corrected.

## Standard commit-and-push

```bash
git status                       # always look first
git add <specific paths>         # prefer explicit paths over -A
git status                       # review what got staged
git commit -m "<short imperative summary>"
git push myfork main
```

Commit-message style in this repo is a short imperative line, no body:
`Add Wall_Collision example case`, `Make initial particle velocity a runtime input`,
`Fix lubrication coefficient and duplicate torque assignment typos`.

Before pushing, scan the staged diff for anything that should not be published —
simulation output, `data/`, `*.o`, `*.mod`, `forces_data.csv`, large binaries.
`.gitignore` covers most of it, but check.

## First-time setup on the supercomputer

Working clone: **`$WORK/CaNS-Fizzy`** (`/g100_work/IscrC_TP-PBR/CaNS-Fizzy`) on
CINECA Galileo100 — code in `$WORK`, run output in `$SCRATCH/runs/`. See
[`run-on-galileo`](../run-on-galileo/SKILL.md) for why, and for staging runs.

(An older clone under `/g100_scratch/.../New_Version/CaNS-Fizzy` is superseded
but may still hold live output — check `squeue` before touching it.)

### 1. SSH key on the cluster

```bash
ssh-keygen -t ed25519 -C "angelorai.favero@gmail.com" -f ~/.ssh/id_ed25519_github
```

`~/.ssh/config`:
```
Host github.com
  HostName github.com
  User git
  IdentityFile ~/.ssh/id_ed25519_github
```

Add the `.pub` key at GitHub → Settings → SSH and GPG keys, then `ssh -T git@github.com`.

**Many HPC login nodes block outbound port 22.** If that test hangs, use GitHub
over 443 instead:
```
Host github.com
  Hostname ssh.github.com
  Port 443
  User git
  IdentityFile ~/.ssh/id_ed25519_github
```

### 2. Clone with submodules

```bash
git clone --recursive git@github.com:AngeloFavero484/CaNS-Fizzy.git
cd CaNS-Fizzy
git remote add origin https://github.com/CaNS-World/CaNS-Fizzy
```

The clone's default remote is named `origin` and points at the **fork**. Rename
it to `myfork` so both machines use the same vocabulary:
```bash
git remote rename origin myfork
git remote add origin https://github.com/CaNS-World/CaNS-Fizzy
```

### 3. Identity

```bash
git config --global user.name  "Angelo Raimondo Favero"
git config --global user.email "angelorai.favero@gmail.com"
```

## Submodules — the recurring failure

`dependencies/2decomp-fft` and `dependencies/cuDecomp` are submodules. A
non-recursive clone leaves them empty and `make` fails with a message that does
not mention submodules at all:

```
make[2]: *** Nessuna regola per generare l'obiettivo "clean".  Arresto.
make[1]: *** [dependencies/external.mk:15: libsclean] Error 2
```

Fix:
```bash
git submodule update --init --recursive
```

## Day-to-day workflow

Editing happens locally; production runs happen on the cluster.

```
local:    edit → commit → git push myfork main
cluster:  git pull myfork main → make → sbatch
```

Pull before starting new work on either machine. Both clones track the same
branch (`main`), so **divergent commits on the two machines mean a manual merge**
— avoid by always pulling first.

If the user edited files directly on the cluster, commit and push from there,
then pull locally. Do not copy files over `scp` between the two clones; it
desynchronises git state.

## Connecting to the cluster yourself

**You can usually SSH to Galileo100 directly — check before asking the user to.**

```bash
ssh afavero0@login.g100.cineca.it            # compute / login nodes
rsync -PravzHS afavero0@data.g100.cineca.it:$SCRATCH/runs/<case>/forces_data.csv .   # data-mover host, for run output
```

There is **no `~/.ssh/config` entry** for CINECA and no long-lived key. Auth is a
short-lived SSH certificate issued by SmallStep, which the user obtains with

```bash
step ssh login angeloraimondo.favero@studenti.unipd.it --provisioner cineca-hpc
```

That command is genuinely interactive (browser SSO + 2FA) — **the user must run
it**, via `! <command>` in the prompt. But it deposits the certificate in the
user's ssh-agent (`SSH_AUTH_SOCK=/run/user/1000/keyring/ssh`), which this session
inherits, so while the certificate is valid every `ssh` and `rsync` works
non-interactively from the Bash tool.

Check for a live certificate before concluding anything:

```bash
ssh-add -l | grep CERT          # an ECDSA-CERT line = a step certificate is loaded
```

and confirm with a cheap probe rather than trusting the agent listing:

```bash
ssh -o BatchMode=yes -o ConnectTimeout=20 afavero0@login.g100.cineca.it 'hostname; echo $WORK'
```

`BatchMode=yes` matters: it makes an expired certificate fail fast instead of
hanging on a password prompt. If the probe fails, the certificate has expired —
ask the user to re-run `step ssh login`, then retry.

### What to run unattended, and what to hand back

Once connected, read-only and idempotent commands are yours to run: `git
pull`/`status`/`log`, `squeue`, `saldo -b`, `ls`, tailing `time.out` or
`forces_data.csv`.

Hand back — or confirm first — anything that spends or destroys:

- `sbatch` and `scancel` — burn (or waste) the project budget
- `make` in `$WORK/CaNS-Fizzy` — never while a job is running out of it
- any `rm` under `$SCRATCH/runs/`

Note that these are not blocked by being "interactive" — they run fine over
`ssh`. They need a decision from the user, which is a different thing.

## Merging upstream changes

```bash
git fetch origin
git merge origin/main
```

Expect conflicts **almost exclusively in `src/main.f90`**. The fork-only files
(`prt_*.f90`, `extend.f90`, `rotnorm.f90`) never conflict. The ordering
invariants that must survive a merge are listed in
`.claude/references/upstream-divergence.md`.
