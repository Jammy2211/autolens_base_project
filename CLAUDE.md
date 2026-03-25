# PyAutoLens Base Project — Claude Instructions

This is the **base template** for PyAutoLens science projects. When asked to create a
new project from this template, follow the conventions below.

---

## Creating a New Project

New projects live outside this repo (e.g. `/mnt/c/Users/Jammy/Results/<project_name>/`).
Use `rsync` to copy the template, excluding what isn't needed.

### Imaging-only project (most common)

```bash
rsync -av \
  --exclude='scripts/interferometer.py' \
  --exclude='hpc/batch_gpu/submit_interferometer' \
  --exclude='hpc/batch_gpu/submit' \
  --exclude='hpc/batch_cpu/submit_interferometer' \
  --exclude='hpc/batch_cpu/template_interferometer' \
  --exclude='hpc/batch_cpu/submit' \
  --exclude='hpc/batch_cpu/template' \
  --exclude='dataset/' \
  --exclude='output/' \
  --exclude='simulators/' \
  --exclude='__pycache__/' \
  --exclude='*.pyc' \
  /mnt/c/Users/Jammy/Code/PyAutoJAX/base_project/ \
  /path/to/new/project/
```

### Interferometer-only project

Same rsync but swap the exclusions: exclude `scripts/imaging.py` and the
`submit_imaging` / `template_imaging` HPC scripts instead.

### Both data types

Omit the script/HPC exclusions entirely.

### What to always exclude

- `dataset/` — add real datasets separately (see below)
- `output/` — never pre-populate; written by PyAutoFit at runtime
- `simulators/` — only needed when generating simulated data

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

## HPC Script Checklist (after copying)

For each new project, update these fields in `hpc/batch_gpu/submit_imaging` and
`hpc/batch_cpu/submit_imaging`:

1. `#SBATCH -J <job_name>` — descriptive name for the SLURM queue
2. `#SBATCH --array=0-N` — set N = number of datasets minus 1
3. `sample=<sample_name>` — matches the subdirectory under `dataset/`
4. `datasets=(...)` — one dataset name per line, in the same order as the array indices
5. `hpc/batch_cpu/template_imaging`: update `sample=` and `dataset=` to a real example lens

The GPU submit script also has `nvidia-smi` in the echo block — leave it in place.

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

## Typical New-Project Workflow

1. `rsync` the template (with appropriate exclusions)
2. Copy or symlink the dataset into `dataset/<sample_name>/`
3. Verify every lens has an `info.json` with at least `pixel_scale` and `n_batch`
   (or confirm the defaults in `scripts/imaging.py` are correct for the instrument)
4. Update `hpc/batch_gpu/submit_imaging` and `hpc/batch_cpu/submit_imaging`:
   job name, `--array`, `sample=`, `datasets=(...)`
5. Update `hpc/batch_cpu/template_imaging`: `sample=` and `dataset=`
6. **Run `dos2unix` on all shell scripts and Python files** to ensure Unix line endings
7. **Add a `Project<Name>()` function to `~/.bashrc`** (see Bash Project Alias above)
8. Test locally on one lens before submitting the full array:
   ```bash
   python3 scripts/imaging.py --sample=<sample> --dataset=<one_lens>
   ```
