# PyAutoLens Base Project — Claude Instructions

This is the **base template** for PyAutoLens science projects. When asked to create a
new project from this template, follow the conventions below.

---

## Creating a New Project

New projects live outside this repo (e.g. `/mnt/c/Users/Jammy/Results/<project_name>/`).
Use `rsync` to copy the template, excluding what isn't needed.

The HPC folder contains one submit script per script type (`submit_imaging`,
`submit_interferometer`, `submit_group`) in both `batch_gpu/` and `batch_cpu/`.
Use the rsync exclusions below to copy **only** the submit scripts that match
the chosen SLaM pipeline(s) — exclude everything else.

To run a single dataset as a test, just put one entry in the `datasets=()` array
in the submit script; no separate template file is needed.

### Imaging-only project (most common)

```bash
rsync -av \
  --exclude='scripts/interferometer.py' \
  --exclude='scripts/group.py' \
  --exclude='hpc/batch_gpu/submit_interferometer' \
  --exclude='hpc/batch_gpu/submit_group' \
  --exclude='hpc/batch_gpu/submit' \
  --exclude='hpc/batch_cpu/submit_interferometer' \
  --exclude='hpc/batch_cpu/submit_group' \
  --exclude='hpc/batch_cpu/submit' \
  --exclude='hpc/batch_cpu/template' \
  --exclude='dataset/' \
  --exclude='output/' \
  --exclude='simulators/' \
  --exclude='__pycache__/' \
  --exclude='*.pyc' \
  /mnt/c/Users/Jammy/Code/PyAutoLabs/autolens_base_project/ \
  /path/to/new/project/
```

### Interferometer-only project

```bash
rsync -av \
  --exclude='scripts/imaging.py' \
  --exclude='scripts/group.py' \
  --exclude='hpc/batch_gpu/submit_imaging' \
  --exclude='hpc/batch_gpu/submit_group' \
  --exclude='hpc/batch_gpu/submit' \
  --exclude='hpc/batch_cpu/submit_imaging' \
  --exclude='hpc/batch_cpu/submit_group' \
  --exclude='hpc/batch_cpu/submit' \
  --exclude='hpc/batch_cpu/template' \
  --exclude='dataset/' \
  --exclude='output/' \
  --exclude='simulators/' \
  --exclude='__pycache__/' \
  --exclude='*.pyc' \
  /mnt/c/Users/Jammy/Code/PyAutoLabs/autolens_base_project/ \
  /path/to/new/project/
```

### Group-only project

```bash
rsync -av \
  --exclude='scripts/imaging.py' \
  --exclude='scripts/interferometer.py' \
  --exclude='hpc/batch_gpu/submit_imaging' \
  --exclude='hpc/batch_gpu/submit_interferometer' \
  --exclude='hpc/batch_gpu/submit' \
  --exclude='hpc/batch_cpu/submit_imaging' \
  --exclude='hpc/batch_cpu/submit_interferometer' \
  --exclude='hpc/batch_cpu/submit' \
  --exclude='hpc/batch_cpu/template' \
  --exclude='dataset/' \
  --exclude='output/' \
  --exclude='simulators/' \
  --exclude='__pycache__/' \
  --exclude='*.pyc' \
  /mnt/c/Users/Jammy/Code/PyAutoLabs/autolens_base_project/ \
  /path/to/new/project/
```

### Multiple data types

Omit the exclusions for any script types you need; keep all others.
Always exclude `submit` and `template` (no suffix) — those are legacy files.

### What to always exclude

- `dataset/` — add real datasets separately (see below)
- `output/` — never pre-populate; written by PyAutoFit at runtime
- `simulators/` — only needed when generating simulated data

## Codex / sandboxed runs

When running Python from Codex or any restricted environment, set writable cache directories so `numba` and `matplotlib` do not fail on unwritable home or source-tree paths:

```bash
NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/matplotlib python scripts/imaging.py
```

This workspace is often imported from `/mnt/c/...` and Codex may not be able to write to module `__pycache__` directories or `/home/jammy/.cache`, which can cause import-time `numba` caching failures without this override.

---

## Dataset Handling

### Directory layout

```
dataset/
└── <sample_name>/
    └── <dataset_name>/
        ├── data.fits
        ├── noise_map.fits
        ├── psf.fits
        ├── positions.json
        └── info.json
```

`sample_name` is the survey/batch name (e.g. `slacs`, `bells`).
`dataset_name` is the individual lens name (e.g. `slacs0737+3216`).

### Copying vs symlinking

**Copy** the dataset when the source may be deleted or reorganised (e.g. copying
from a `Results/old_project/dataset/` folder before deleting it).

**Symlink** only when the source is stable and permanent (e.g. a shared NFS mount
on the HPC, or a dedicated raw-data archive that will never move).

```bash
# Copy (preferred for real-data projects that will be archived)
cp -r /path/to/source/slacs  dataset/slacs

# Symlink (only when source is stable)
ln -s /path/to/source/slacs  dataset/slacs
```

---

## info.json Fields

Every dataset directory needs an `info.json`. The imaging script reads all values
via `info.get(key, default)` so fields can be omitted when the default is correct.

| Field | Default | Notes |
|---|---|---|
| `pixel_scale` | `0.05` | Arcsec/pixel. HST ≈ 0.05, Euclid ≈ 0.1 |
| `n_batch` | `50` | Pixelization batch size. Lower for high-res data |
| `mask_radius` | `3.5` | Circular mask radius in arcsec |
| `subhalo_grid_dimensions_arcsec` | `3.0` | Grid search half-width for subhalo pipeline |
| `redshift_lens` | `0.5` | Used by all SLAM stages |
| `redshift_source` | `1.0` | Used by all SLAM stages |

Interferometer datasets additionally support `real_space_shape` ([256,256]) and the
same `mask_radius`.

---

## HPC (`hpc/`)

The `hpc/` directory contains everything needed to submit, sync, and monitor SLURM
jobs on the HPC cluster.

### Directory Structure

```
hpc/
├── batch_gpu/                  # GPU submit scripts + SLURM log dirs
│   ├── submit_imaging          # SLURM batch script for imaging pipeline
│   ├── submit_interferometer   # SLURM batch script for interferometer pipeline
│   ├── submit_group            # SLURM batch script for group pipeline
│   ├── submit                  # LEGACY — do not use, do not modify
│   ├── output/                 # SLURM stdout logs (*.out)
│   └── error/                  # SLURM stderr logs (*.err)
├── batch_cpu/                  # CPU submit scripts + SLURM log dirs
│   ├── submit_imaging
│   ├── submit_interferometer
│   ├── submit_group
│   ├── submit                  # LEGACY — do not use, do not modify
│   ├── template                # LEGACY — do not use, do not modify
│   ├── output/
│   └── error/
├── sync                        # Bidirectional sync script (local ↔ HPC)
├── sync.conf.example           # Template config for sync
├── sync_jump                   # Two-hop relay sync (local → jump → build → cosma → local)
├── sync_jump.conf.example      # Template config for sync_jump
├── .gitignore                  # Ignores sync.conf, sync_jump.conf, subhalo/
└── __init__.py
```

### Submit Scripts — GPU vs CPU

Each script type (`imaging`, `interferometer`, `group`) has a submit script in both
`batch_gpu/` and `batch_cpu/`. The key differences:

| | GPU (`batch_gpu/`) | CPU (`batch_cpu/`) |
|---|---|---|
| Partition | `--partition=gpu` | `--partition=cpu` |
| GPU | `--gres=gpu:1` | none |
| CPUs | `--cpus-per-task=1` | `--cpus-per-task=8` |
| Memory | `--mem=32gb` | `--mem=64gb` |
| Wall time | `--time=08:00:00` | `-t 18:00:00` |
| JAX platform | (uses GPU by default) | Forces `JAX_PLATFORM_NAME=cpu` |
| Thread pinning | none | Sets `OPENBLAS/MKL/OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK` |
| Echo block | Includes `nvidia-smi` | No `nvidia-smi` |

**CPU scripts set these environment variables** to pin threads and force CPU-only JAX:

```bash
export JAX_PLATFORM_NAME=cpu
export JAX_PLATFORMS=cpu
THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS
export VECLIB_MAXIMUM_THREADS=$THREADS
export NUMEXPR_NUM_THREADS=$THREADS
export NPROC=$THREADS
```

### Submit Script Structure

All submit scripts follow the same pattern:

1. **SBATCH headers** — job name, partition, resources, array range, log paths, email
2. **Environment** — `source $PROJECT_PATH/activate.sh` (must set `PROJECT_PATH` before submitting)
3. **Sample** — `sample=<sample_name>` matches subdirectory under `dataset/`
4. **Dataset list** — `datasets=(...)` array, one dataset name per line
5. **Array task selection** — `dataset="${datasets[$SLURM_ARRAY_TASK_ID]}"`
6. **Run** — `python3 $PROJECT_PATH/scripts/<type>.py --sample=$sample --dataset=$dataset`

### HPC Script Checklist (after copying)

For each script type present in the project (`imaging`, `interferometer`, `group`),
update these fields in both `hpc/batch_gpu/submit_<type>` and
`hpc/batch_cpu/submit_<type>`:

1. `#SBATCH -J <job_name>` — descriptive name for the SLURM queue
2. `#SBATCH --array=0-N` — set N = number of datasets minus 1
3. `sample=<sample_name>` — matches the subdirectory under `dataset/`
4. `datasets=(...)` — one dataset name per line, in the same order as the array indices

The GPU submit scripts have `nvidia-smi` in the echo block — leave it in place.

To test a single lens, temporarily set `--array=0-0` and put just that lens in
`datasets=(...)` — no separate template file is needed.

### Legacy Files

`batch_gpu/submit`, `batch_cpu/submit`, and `batch_cpu/template` (no suffix) are
**legacy files**. Do not use, modify, or copy them to new projects. Always exclude
them in rsync commands. They use an older format without `--sample` / `--dataset`
separation and reference `scripts/base.py` instead of the typed scripts.

### `hpc/sync` — Bidirectional Project Sync

A single script that handles all transfer and job management between your local
machine and the HPC. Run from the project root or `hpc/` directory.

**Setup:**
```bash
cp hpc/sync.conf.example hpc/sync.conf
# Edit hpc/sync.conf with your HPC host, base path, and project name.
# sync.conf is gitignored — it stays on your local machine only.
```

**Config fields** (`sync.conf`):
- `HPC_HOST` — SSH host alias or `user@hostname`
- `HPC_BASE` — base directory on the HPC (e.g. `/mnt/ral/jnightin`)
- `PROJECT_NAME` — defaults to local folder name if unset

The remote path is `$HPC_HOST:$HPC_BASE/$PROJECT_NAME`.

**Transfer commands:**

| Command | Description |
|---|---|
| `hpc/sync push` | Upload code, config, and data to the HPC |
| `hpc/sync push --no-data` | Upload code only (skip `dataset/`) |
| `hpc/sync pull` | Download SLURM logs then results from the HPC |
| `hpc/sync logs` | Download SLURM output/error logs only (fast, use mid-run) |
| `hpc/sync sync` | Push then pull (default if no command given) |
| `hpc/sync sync --no-data` | Push code only, then pull |
| `hpc/sync push-data-init` | First-time dataset upload via tar pipe (faster for initial large uploads) |
| `hpc/sync pull-full` | Full output download via tar pipe (avoids per-file rsync overhead) |
| `hpc/sync status` | Dry run — show what would transfer without transferring |

**Job commands** (no manual SSH required):

| Command | Description |
|---|---|
| `hpc/sync submit [gpu\|cpu] <script>` | Submit a SLURM job (e.g. `submit gpu submit_imaging`) |
| `hpc/sync push-submit [gpu\|cpu] <script>` | Push code then submit in one step |
| `hpc/sync jobs` | Show queued/running jobs (`squeue`) |
| `hpc/sync sacct` | Show job history and exit codes |
| `hpc/sync cancel <job_id>` | Cancel a job by ID |
| `hpc/sync wait-and-pull [secs]` | Poll until all jobs finish, then pull (default: 60s interval) |

**Inspect commands:**

| Command | Description |
|---|---|
| `hpc/sync tail [gpu\|cpu]` | Stream live SLURM log output (Ctrl+C to stop; default: gpu) |
| `hpc/sync du` | Show remote disk usage |
| `hpc/sync check` | Verify SSH connection and remote paths |
| `hpc/sync clear-logs [gpu\|cpu]` | Delete SLURM output/error log files (local + remote) |

**What gets synced:**

- **Push (code):** `config/`, `hpc/`, `scripts/`, `slam_pipeline/`, `simulators/` + root files (`activate.sh`, `util.py`, `__init__.py`, `README.rst`, `LICENSE`). Changed files are updated normally.
- **Push (data):** `dataset/` — uses `--ignore-existing` so FITS files already on the HPC are never re-transferred.
- **Pull (logs):** `hpc/batch_gpu/output/`, `hpc/batch_gpu/error/`, `hpc/batch_cpu/output/`, `hpc/batch_cpu/error/`
- **Pull (results):** `output/` — excludes `search_internal/` (large sampler state not needed locally).
- **Always excluded:** `__pycache__/`, `*.pyc`, `.git/`, `*.egg-info/`, `sync.conf`

**rsync options:** archive mode, compression (skipping FITS/gz/bz2/xz/zst), partial resume, SSH ControlMaster connection reuse.

### `hpc/sync_jump` — Two-Hop Relay Sync

For topologies where results live on a build server only reachable via a jump host,
and must be staged through an intermediate server before reaching your local machine:

```
local  ──ssh──►  JUMP_HOST  ──ssh──►  BUILD_HOST
                                           │
                                      rsync/tar
                                           │
                                           ▼
local  ◄──rsync──  COSMA_HOST  ◄──rsync──  BUILD_HOST
```

**Setup:**
```bash
cp hpc/sync_jump.conf.example hpc/sync_jump.conf
# Edit sync_jump.conf with your host names and paths.
```

**Config fields** (`sync_jump.conf`):
- `JUMP_HOST` — first hop from local machine (e.g. `euclid_jump`)
- `BUILD_HOST` — build server, only reachable via JUMP_HOST (e.g. `euclid-ral-build`)
- `BUILD_OUTPUT_PATH` — path to `output/` on the build server
- `COSMA_HOST` — intermediate staging server (e.g. `cosma8`)
- `COSMA_STAGING_DIR` — staging directory on COSMA_HOST
- `PROJECT_NAME` — used to namespace archives; defaults to local folder name

**Commands:**

| Command | Description |
|---|---|
| `hpc/sync_jump push` | Relay output from build server → cosma staging (rsync) |
| `hpc/sync_jump push --zip` | Same, but via a single tar.gz archive |
| `hpc/sync_jump pull` | Download cosma staging → local `output/` |
| `hpc/sync_jump pull --zip` | Download archive then extract locally |
| `hpc/sync_jump sync [--zip]` | Push then pull (default) |
| `hpc/sync_jump status` | Dry run — show what would transfer |

**Options:**
- `--zip` — transfer a single tar.gz instead of many small files (faster when `output/` has thousands of files)
- `--include-search-internal` — include `search_internal/` dirs (excluded by default)

**Requires:** `ssh-agent` running with your key loaded (`ssh-add -L` should list a key).

### `.gitignore`

The `hpc/.gitignore` ignores:
- `subhalo/` — subhalo grid search output (generated at runtime)
- `sync.conf` — local HPC connection config (contains host-specific paths)
- `sync_jump.conf` — local jump-host connection config

---

## Scripts and info.json

`scripts/imaging.py` reads all dataset-specific values from `info.json` using
`info.get(key, default)`. Hard-coded values for `mask_radius`,
`subhalo_grid_dimensions_arcsec`, `pixel_scale`, and `n_batch` should never appear
in the script body — always source them from info.json.

`scripts/interferometer.py` similarly reads `pixel_scale`, `n_batch`,
`real_space_shape`, and `mask_radius` from info.json.

---

## slam_pipeline/ — Do Not Modify

`slam_pipeline/` is dataset-type agnostic. Never modify these files when setting up
a new project. Project-specific changes belong in `scripts/imaging.py` or
`scripts/interferometer.py`.

---

## Line Endings — Always Unix (LF)

All files in this project **must use Unix line endings (LF, `\n`)**. Windows/DOS
line endings (CRLF, `\r\n`) will break shell scripts and Python files on the HPC.

**When writing or editing any file**, always produce Unix line endings. Never write
`\r\n` line endings.

After creating or copying files, verify and convert if needed:

```bash
# Check for DOS line endings
file hpc/batch_gpu/submit_imaging   # should say "ASCII text", not "CRLF"

# Convert a single file
dos2unix hpc/batch_gpu/submit_imaging

# Convert all shell scripts and Python files in the project
find . -type f \( -name "*.py" -o -name "*.sh" -o -name "submit*" -o -name "template*" \) \
  | xargs dos2unix
```

---

## Test Runs

A "test run" means running a script with `PYAUTOFIT_TEST_MODE=1`, which makes all
non-linear searches complete almost instantly with a trivial number of samples. Use
this to verify the full pipeline executes without errors before submitting to the HPC.

```bash
# Imaging
PYAUTOFIT_TEST_MODE=1 python3 scripts/imaging.py --sample=<sample> --dataset=<dataset>

# Interferometer
PYAUTOFIT_TEST_MODE=1 python3 scripts/interferometer.py --sample=<sample> --dataset=<dataset>

# Group
PYAUTOFIT_TEST_MODE=1 python3 scripts/group.py --sample=<sample> --dataset=<dataset>
```

Example datasets for each script type live at:
- Imaging: `dataset/sample_imaging/example_imaging/`
- Interferometer: `dataset/sample_interferometer/example_interferometer/`
- Group: `dataset/sample_group/102021990_NEG650312660474055399/`

---

## Bash Project Alias

Every new project gets a shell function in `~/.bashrc` that activates the venv and
`cd`s into the project directory. Add it immediately after creating the project,
grouped with the other `Project*` functions:

```bash
Project<ProjectName>() {
  source ~/venv/PyAuto/bin/activate
  cd /mnt/c/Users/Jammy/Results/<project_name>
}
```

Use the `PyAuto` venv unless the project requires a different one.

---

## Context (`context/`)

The `context/` folder provides AI agents with the scientific and technical background
needed to work on the project. It contains files copied from the software's workspace
repository (e.g. `autolens_workspace`) — tutorials, feature examples, guide scripts,
and reference material that explain how the APIs, modeling conventions, and science
cases work.

**Purpose:** When an agent needs to understand how a feature works, interpret modeling
results, or make decisions about pipeline configuration, `context/` is where it looks
first. The files bridge the gap between the raw API and the science — they explain
*why* certain choices are made, not just *how* to call the code.

**When to read:** Always read relevant files in `context/` before modifying scripts,
interpreting results, or advising on modeling choices. For example:
- Before changing pixelization settings → read the pixelization tutorial/example
- Before interpreting subhalo results → read the subhalo modeling guide
- Before adjusting mass model configuration → read the relevant feature example

**Typical contents** (varies per project):
- Feature examples (e.g., pixelization, subhalo modeling, multi-Gaussian expansion)
- Guide scripts explaining API usage or modeling conventions
- Reference outputs or worked examples relevant to the science case
- Tutorials from `autolens_workspace/notebooks/` converted or copied as `.py` scripts

**Population:** The `context/` folder is empty in a fresh rsync of the base template.
It is populated manually per-project by copying relevant files from the workspace
repository (e.g. `autolens_workspace/scripts/`, `autolens_workspace/notebooks/`).
There is no automated population step — the user selects which files are relevant
to the specific science case and copies them in.

**Source repositories** (common):
- `autolens_workspace` — PyAutoLens tutorials, feature examples, guides
- `autofit_workspace` — PyAutoFit non-linear search tutorials, analysis examples
- `autogalaxy_workspace` — galaxy modeling tutorials (light profiles, mass profiles)

---

## Modeling Scripts (`Scripts/`)

The `Scripts/` folder is empty in the base template. After rsync-ing the template into a new
project, use the `/init-slam` skill to select and copy the appropriate SLaM pipeline script(s)
from `autolens_workspace`. The skill presents categorized options, copies the chosen script(s),
and creates `Scripts/slam_claude.md` with full SLaM context for future AI sessions.

The skill is defined at `autolens_base_project/skills/init-slam/SKILL.md`. Install it once from there.

See `Scripts/CLAUDE.md` for the full list of available pipeline options.

---

## Typical New-Project Workflow

1. `rsync` the template (with appropriate exclusions)
2. **Run `/init-slam`** to select and copy SLaM pipeline script(s) into `Scripts/`
3. Copy or symlink the dataset into `dataset/<sample_name>/`
4. Verify every lens has an `info.json` with at least `pixel_scale` and `n_batch`
   (or confirm the defaults in `Scripts/imaging.py` are correct for the instrument)
5. Update `hpc/batch_gpu/submit_<type>` and `hpc/batch_cpu/submit_<type>` for each
   chosen script type: job name, `--array`, `sample=`, `datasets=(...)`
7. **Run `dos2unix` on all shell scripts and Python files** to ensure Unix line endings
8. **Add a `Project<Name>()` function to `~/.bashrc`** (see Bash Project Alias above)
9. Test locally on one lens before submitting the full array:
   ```bash
   python3 Scripts/imaging.py --sample=<sample> --dataset=<one_lens>
   ```
