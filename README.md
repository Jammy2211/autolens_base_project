# PyAutoLens Science Project

This folder is the base template for a PyAutoLens gravitational lens analysis project.
Copy and adapt it as the starting point for a new science project.

---

## Project Structure

```
project/
├── config/           # PyAutoLens configuration (priors, non-linear samplers, visualisation)
├── dataset/          # Imaging and interferometer data (see Dataset Layout below)
├── hpc/              # HPC batch submission scripts
│   ├── batch_cpu/    # CPU job scripts + SLURM output/error logs
│   └── batch_gpu/    # GPU job scripts + SLURM output/error logs
├── output/           # Analysis results (written automatically by PyAutoFit)
├── scripts/          # Analysis scripts — run locally or on the HPC unchanged
│   ├── imaging.py    # SLAM pipeline for imaging data
│   └── interferometer.py  # SLAM pipeline for interferometer data
├── simulators/       # Scripts for generating simulated datasets
│   └── base.py
└── slam_pipeline/    # SLAM pipeline stage definitions (dataset-type agnostic)
```

---

## Dataset Layout

Datasets live inside `dataset/` and are organised into **samples**. A sample is a named
subdirectory that groups related datasets (e.g. all lenses from a survey). Each dataset
is then a subdirectory inside its sample.

```
dataset/
└── <sample>/
    ├── <dataset_1>/
    │   ├── data.fits
    │   ├── noise_map.fits
    │   ├── psf.fits          # imaging only
    │   ├── uv_wavelengths.fits  # interferometer only
    │   ├── positions.json
    │   └── info.json
    └── <dataset_2>/
        └── ...
```

The included example datasets use separate sample folders per data type:

```
dataset/
├── sample_imaging/
│   └── example_imaging/
│       ├── data.fits
│       ├── noise_map.fits
│       ├── psf.fits
│       ├── positions.json
│       └── info.json
└── sample_interferometer/
    └── example_interferometer/
        ├── data.fits
        ├── noise_map.fits
        ├── uv_wavelengths.fits
        ├── positions.json
        └── info.json
```

### info.json

Every dataset directory must contain an `info.json` file. This is the single source of
truth for dataset-specific properties used by analysis scripts. It removes any need to
hard-code or pass these values as arguments.

Required fields:

```json
{
    "pixel_scale": 0.05,
    "n_batch": 40
}
```

| Field         | Type  | Description                                                                 |
|---------------|-------|-----------------------------------------------------------------------------|
| `pixel_scale` | float | Arcseconds per pixel. Varies by instrument (e.g. HST ≈ 0.05, Euclid ≈ 0.1) |
| `n_batch`     | int   | Pixelization batch size. Use lower values for higher-resolution data to reduce VRAM usage (e.g. 40 for HST, 8 for AO) |

Optional fields (used by the SLAM pipeline if present):

```json
{
    "pixel_scale": 0.05,
    "n_batch": 40,
    "redshift_lens": 0.5,
    "redshift_source": 1.0
}
```

For interferometer datasets, two additional optional fields are supported:

```json
{
    "pixel_scale": 0.1,
    "n_batch": 25,
    "real_space_shape": [256, 256],
    "mask_radius": 3.5
}
```

| Field              | Type       | Description                                                                    |
|--------------------|------------|--------------------------------------------------------------------------------|
| `real_space_shape` | [int, int] | Height × width of the real-space reconstruction grid (default `[256, 256]`)    |
| `mask_radius`      | float      | Circular mask radius in arcseconds (default `3.5`)                             |

When generating a simulated dataset with `simulators/base.py`, `info.json` is written
automatically alongside the data.

For real observational data, create `info.json` manually or with a preprocessing script.

---

## Running Scripts

### Locally

Run from anywhere — paths are resolved relative to the script's location:

```bash
python3 scripts/imaging.py --sample=<sample> --dataset=<dataset>
python3 scripts/interferometer.py --sample=<sample> --dataset=<dataset>
```

The included examples:

```bash
python3 scripts/imaging.py --sample=sample_imaging --dataset=example_imaging
python3 scripts/interferometer.py --sample=sample_interferometer --dataset=example_interferometer
```

Both `--sample` and `--dataset` are optional. Output paths are organised as
`output/<sample>/<dataset>/<pipeline_stage>/`.

### On the HPC

The batch scripts in `hpc/batch_gpu/` and `hpc/batch_cpu/` handle all
HPC-specific concerns (SLURM directives, environment activation, paths).
The Python script itself requires no modification between local and HPC runs.

All batch scripts use a `$PROJECT_PATH` environment variable so no paths are
hard-coded in the scripts. Set it once before submitting:

```bash
export PROJECT_PATH=/path/to/your/project
```

Each dataset type has its own set of batch scripts. Imaging scripts call
`scripts/imaging.py`; interferometer scripts call `scripts/interferometer.py`.

**GPU (recommended):**

| Script | Purpose |
|--------|---------|
| `hpc/batch_gpu/submit_imaging` | Array job for imaging datasets |
| `hpc/batch_gpu/submit_interferometer` | Array job for interferometer datasets |

1. Edit the appropriate submit script:
   - Update `--mail-user` for your email
   - Set `sample=` to the sample subdirectory name
   - Populate the `datasets` array with the dataset names to run
   - Update `--array=0-N` to match the number of datasets
   - Adjust `--mem` and `--time` as needed

2. Submit from the `hpc/batch_gpu/` directory:

```bash
cd hpc/batch_gpu
export PROJECT_PATH=/path/to/your/project
sbatch submit_imaging          # imaging
sbatch submit_interferometer   # interferometer
```

**CPU:**

| Script | Purpose |
|--------|---------|
| `hpc/batch_cpu/submit_imaging` | CPU array job for imaging datasets |
| `hpc/batch_cpu/submit_interferometer` | CPU array job for interferometer datasets |
| `hpc/batch_cpu/template_imaging` | Single-dataset imaging job template |
| `hpc/batch_cpu/template_interferometer` | Single-dataset interferometer job template |

```bash
cd hpc/batch_cpu
export PROJECT_PATH=/path/to/your/project
sbatch submit_imaging          # imaging
sbatch submit_interferometer   # interferometer
```

SLURM logs are written to the `output/` and `error/` subdirectories inside each batch folder.

---

## Syncing with the HPC

`hpc/sync` is a single script that handles all data movement between your local
machine and the HPC. It wraps `rsync` with sensible defaults and transfers only
what has actually changed.

### First-time setup

```bash
cp hpc/sync.conf.example hpc/sync.conf
# Edit hpc/sync.conf — set HPC_HOST, HPC_BASE, and PROJECT_NAME
```

`sync.conf` is gitignored and stays on your local machine only.

### Commands

```bash
hpc/sync push     # Upload code, config, and data to the HPC
hpc/sync pull     # Download results from the HPC
hpc/sync sync     # Push then pull (default)
hpc/sync status   # Dry run — see what would transfer without moving anything
```

### What gets transferred

| Direction | Folders | Strategy |
|-----------|---------|----------|
| push | `config/` `hpc/` `scripts/` `slam_pipeline/` `simulators/` | Normal sync — only changed files |
| push | `dataset/` | `--ignore-existing` — skips files already on HPC, avoiding re-checksumming large FITS archives |
| pull | `output/` | `--update --exclude=search_internal` — only downloads files newer than local copies, omits large sampler internals |

The `--ignore-existing` flag on dataset is the key optimisation for large projects:
once a FITS file is on the HPC, it is never re-examined on subsequent syncs.

### Connection to HPC batch scripts

`$HPC_BASE/$PROJECT_NAME` in `sync.conf` is the same path as `$PROJECT_PATH`
used inside the SLURM batch scripts, so activation paths and script calls stay
consistent across local, push, and job submission steps.

---

## Configuration

`config/` contains all PyAutoLens configuration files. The HPC jobs use the same
`config/` as local runs — there is no separate HPC config.

Key config files:

| File | Purpose |
|------|---------|
| `config/general.yaml` | Global settings |
| `config/non_linear/nest.yaml` | Nested sampling settings (Nautilus / MultiNest) |
| `config/priors/` | Prior distributions for all model components |
| `config/visualize/` | Matplotlib output settings |

---

## SLAM Pipeline

`slam_pipeline/` contains the modular pipeline stages:

| Module | Stage |
|--------|-------|
| `source_lp.py` | Parametric source (light profile) |
| `source_pix.py` | Pixelised source (mesh + regularization) |
| `light_lp.py` | Lens light |
| `mass_total.py` | Total mass |
| `subhalo/detection.py` | Dark matter subhalo detection |

---

## Simulating Data

`simulators/base.py` generates a synthetic imaging dataset. Edit the dataset properties
at the top of the file (`pixel_scale`, `shape_native`, `n_batch`) then run:

```bash
# Single simulated dataset
python3 simulators/base.py

# Named subdirectory
python3 simulators/base.py my_dataset
```

The simulator writes `info.json` automatically, so analysis scripts will pick up the
correct `pixel_scale` and `n_batch` without any further configuration.
