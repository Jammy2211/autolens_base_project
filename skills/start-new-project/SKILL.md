---
name: start-new-project
description: Create a new science project from the autolens_base_project template. Use this skill when the user says "start a new project", "create a project", "set up a new project", "new project from template", or wants to initialize a fresh project directory with datasets, scripts, simulators, and context. Walks through project name, description, datasets, modeling scripts, simulators, and context folder step by step.
user-invocable: true
---

# Start New Project

Interactive setup for a new science project from the `autolens_base_project` template. Guides the user
through naming, describing, and populating a project step by step.

## Intro — Explain the base project

Before asking any questions, explain the following to the user:

> This repository (`autolens_base_project`) is a **template** for scientific modeling projects. It contains
> the standard directory structure, HPC submission scripts, sync tools, and pipeline scaffolding
> that every project needs.
>
> When you create a new project, this template is copied to a fresh directory (under
> `/mnt/c/Users/Jammy/Results/<project_name>/`). The new project starts mostly empty — you then
> populate it with your datasets, modeling scripts, simulators, and reference context.
>
> **Modeling scripts, simulators, and context files** are normally adapted from the **workspace
> repository** of the software you're using. For example, if you're using PyAutoLens, the
> `autolens_workspace` repository contains template scripts for different modeling pipelines,
> simulator examples, and tutorials that you can copy and customize for your science case.
> Other workspaces include `autofit_workspace`, `autogalaxy_workspace`, etc.
>
> I'll walk you through six steps to set everything up. Most steps are optional — you can
> always add things later.

Then proceed to Step 1.

## Step 1 — Project name

Ask the user:

> **What should this project be called?**
>
> This will be the folder name, created at `/mnt/c/Users/Jammy/Results/<project_name>/`.
> Use a short, descriptive name (e.g. `slacs_subhalo`, `euclid_pilot`, `bells_mge`).

Store the answer as `PROJECT_NAME`. Confirm the full path will be
`/mnt/c/Users/Jammy/Results/<PROJECT_NAME>/`.

## Step 2 — Project description

Ask the user:

> **Give a short description of this project's scientific goals.**
>
> This will be written into the project's `CLAUDE.md` so that future AI sessions understand
> the purpose and scope of the project. A sentence or two is fine — e.g. "Subhalo detection
> in SLACS lenses using pixelized source reconstruction" or "MGE modeling of Euclid strong
> lens candidates for the initial catalog."

Store the answer as `PROJECT_DESCRIPTION`.

## Step 3 — Datasets

Ask the user:

> **Do you have datasets to include in this project?**
>
> Datasets live under `dataset/<sample_name>/<dataset_name>/`, where `sample_name` is the
> survey or batch (e.g. `slacs`, `bells`) and `dataset_name` is the individual lens (e.g.
> `slacs0737+3216`). Each dataset directory should contain at minimum `data.fits`,
> `noise_map.fits`, and an `info.json` with fields like `pixel_scale`, `redshift_lens`, etc.
>
> You can either:
> - **Point me to dataset paths** and I'll copy them into the project
> - **Skip for now** and add datasets later
>
> If pointing to paths, give me the source directory and tell me the sample name to file them
> under.

If the user provides paths, validate they exist and note the sample/dataset structure. If they
skip, move on.

## Step 4 — Modeling scripts

Ask the user:

> **Do you have modeling scripts for this project?**
>
> Modeling scripts define the analysis pipeline (e.g. SLaM stages for source, light, mass,
> subhalo detection). They live in `scripts/` in the project.
>
> These are normally **adapted from the workspace** of the software you're using. For example,
> `autolens_workspace/scripts/` contains template pipelines for imaging, interferometer, and
> group-scale lensing, organized by feature (pixelization, MGE, subhalo, etc.).
>
> You can either:
> - **Point me to specific script paths** and I'll copy them in
> - **Describe your science case** and I'll scan the workspace to suggest relevant templates
> - **Skip for now** — you can add scripts later or use the `/init-slam` skill
>
> If you want me to scan a workspace, tell me which workspace to look in (e.g.
> `autolens_workspace`) and describe what kind of modeling you need.

If the user asks for a workspace scan, use the Explore agent to search the workspace's
`scripts/` directory for relevant pipelines based on the user's science description. Present
the findings and let the user pick which scripts to copy.

If the user points to specific files, copy them to `scripts/` in the project.

If they skip, move on.

## Step 5 — Simulators

Ask the user:

> **Do you want to set up any simulators?**
>
> The `simulators/` folder contains scripts that generate mock/simulated data for your project.
> This is useful for testing your pipeline on known inputs before running on real data, or for
> sensitivity analyses (e.g. simulating subhalos at different masses and positions).
>
> Simulator templates are normally found in the workspace repository — for example,
> `autolens_workspace/scripts/simulators/` has examples for imaging, interferometer, and
> group-scale data.
>
> You can either:
> - **Point me to specific simulator scripts** and I'll copy them in
> - **Describe what you want to simulate** and I'll scan the workspace for relevant examples
> - **Skip for now** — simulators can be added any time

If the user asks for a workspace scan, search `simulators/` in the workspace and present
options. Otherwise copy specified files or skip.

## Step 6 — Context

Ask the user:

> **Do you want to add any context files?**
>
> The `context/` folder holds tutorials, examples, and reference material from workspace
> repositories. These files give AI agents the scientific and technical background needed to
> work on the project — for example, explaining how pixelization works, what MGE profiles are,
> or how subhalo sensitivity mapping is done.
>
> Context files are copied from workspace repos like `autolens_workspace/scripts/` or
> `autolens_workspace/notebooks/`.
>
> You can either:
> - **Point me to specific files** to include
> - **Describe topics you'll need context for** and I'll scan the workspace for relevant
>   tutorials and examples
> - **Skip for now** — context can be added at any time as the project evolves
>
> If you're not sure yet, that's fine — it's easy to add context later as you discover what
> background the AI needs.

If the user asks for a scan, search the workspace for relevant tutorials, guides, and feature
examples. Present findings and let the user pick. Otherwise copy specified files or skip.

## Step 7 — Execute

Once all questions are answered, perform these actions:

### 7a — rsync the template

Copy the base template to the new project directory. Always exclude `dataset/`, `output/`,
`simulators/`, `__pycache__/`, `*.pyc`, and legacy files (`submit` and `template` without
suffix):

```bash
rsync -av \
  --exclude='dataset/' \
  --exclude='output/' \
  --exclude='simulators/' \
  --exclude='__pycache__/' \
  --exclude='*.pyc' \
  --exclude='hpc/batch_gpu/submit' \
  --exclude='hpc/batch_cpu/submit' \
  --exclude='hpc/batch_cpu/template' \
  --exclude='skills/' \
  /mnt/c/Users/Jammy/Code/PyAutoLabs/autolens_base_project/ \
  /mnt/c/Users/Jammy/Results/<PROJECT_NAME>/
```

### 7b — Copy datasets

If the user provided dataset paths, copy them into
`/mnt/c/Users/Jammy/Results/<PROJECT_NAME>/dataset/<sample_name>/`. Verify each dataset
directory has at least `data.fits` and `info.json`.

### 7c — Copy modeling scripts

If the user provided or selected scripts, copy them into
`/mnt/c/Users/Jammy/Results/<PROJECT_NAME>/scripts/`.

### 7d — Copy simulator files

If the user provided or selected simulators, copy them into
`/mnt/c/Users/Jammy/Results/<PROJECT_NAME>/simulators/`.

### 7e — Copy context files

If the user provided or selected context files, copy them into
`/mnt/c/Users/Jammy/Results/<PROJECT_NAME>/context/`.

### 7f — Update CLAUDE.md

Open `/mnt/c/Users/Jammy/Results/<PROJECT_NAME>/CLAUDE.md`. Add a project-specific section
at the very top, before the existing template instructions:

```markdown
# <PROJECT_NAME>

<PROJECT_DESCRIPTION>

---

```

### 7g — Run dos2unix

Convert all shell scripts and Python files to Unix line endings:

```bash
find /mnt/c/Users/Jammy/Results/<PROJECT_NAME>/ \
  -type f \( -name "*.py" -o -name "*.sh" -o -name "submit*" -o -name "template*" \) \
  | xargs dos2unix
```

### 7h — Add bash alias

Append a `Project<Name>()` function to `~/.bashrc`:

```bash
Project<ProjectName>() {
  source ~/venv/PyAuto/bin/activate
  cd /mnt/c/Users/Jammy/Results/<PROJECT_NAME>
}
```

Use PascalCase for the function name (e.g. `SlacsSubhalo` for project `slacs_subhalo`).

## Step 8 — GitHub repository

Ask the user:

> **Do you want to create a GitHub repository for this project?**
>
> A GitHub repo lets you manage the project remotely — edit scripts, push changes, and
> trigger HPC jobs from your phone or any device. Only **code** goes in the repo (scripts,
> config, HPC submit files, context). Large data files (`dataset/`, `output/`) stay excluded
> via `.gitignore` — data is already on the HPC from `hpc/sync push` and doesn't need to
> live in GitHub.
>
> Options:
> - **Private repo** — visible only to you (recommended for unpublished science)
> - **Public repo** — open to everyone
> - **Skip** — no GitHub repo for now

If the user wants a repo:

### 8a — Create .gitignore

Write a `.gitignore` at the project root with these entries:

```
# Data — too large for GitHub; lives on HPC via hpc/sync
dataset/
output/

# Simulator output (generated data)
simulators/output/

# HPC local config (host-specific paths, never committed)
hpc/sync.conf
hpc/sync_jump.conf

# SLURM logs (pulled from HPC, not version-controlled)
hpc/batch_gpu/output/
hpc/batch_gpu/error/
hpc/batch_cpu/output/
hpc/batch_cpu/error/

# Python
__pycache__/
*.pyc
*.egg-info/

# OS
.DS_Store
Thumbs.db
```

### 8b — Initialize git and create the repo

```bash
cd /mnt/c/Users/Jammy/Results/<PROJECT_NAME>
git init
git add -A
git commit -m "Initial project setup from autolens_base_project template"
```

Then determine the user's GitHub username:

```bash
gh api user --jq .login
```

Store the result as `GH_USER`. Create the repo:

```bash
# Private (default):
gh repo create <GH_USER>/<PROJECT_NAME> --private --source=. --push

# Public:
gh repo create <GH_USER>/<PROJECT_NAME> --public --source=. --push
```

Use `--private` or `--public` based on the user's choice.

### 8c — Confirm

Report:
- Repo URL: `https://github.com/<GH_USER>/<PROJECT_NAME>`
- What's tracked: scripts, config, HPC submit files, slam_pipeline, context
- What's excluded: dataset/, output/, SLURM logs, sync configs
- Remind them: data stays on the HPC. To work from your phone, edit code via GitHub,
  then SSH to the HPC to pull and submit jobs.

## Step 9 — Final summary

Report a summary:
- Project path: `/mnt/c/Users/Jammy/Results/<PROJECT_NAME>/`
- What was copied: datasets, scripts, simulators, context (or "none" for each)
- Bash alias added
- GitHub repo (if created): URL and visibility
- Remind the user of next steps:
  - Verify `info.json` in each dataset if datasets were copied
  - Update HPC submit scripts (`hpc/batch_gpu/submit_*`, `hpc/batch_cpu/submit_*`) with job
    name, `--array` range, `sample=`, and `datasets=(...)`
  - Run a test: `PYAUTOFIT_TEST_MODE=1 python3 scripts/<type>.py --sample=<sample> --dataset=<dataset>`
  - Use `/init-slam` if they skipped modeling scripts and want to set up SLaM pipelines later
  - If GitHub repo was created: push changes with `git add -A && git commit -m "..." && git push`
    after making edits
