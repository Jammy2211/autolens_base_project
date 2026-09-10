# skills/

Procedural how-to-do-X skills for the PyAuto\* lensing stack. Each skill is a single
Markdown file with YAML frontmatter; the body teaches an agent (and through them, the
user) how to write Python that accomplishes one lensing task.

Skills are also exposed at `.claude/skills/` (Claude Code) and `~/.codex/skills/` (when
configured) via symlinks; the canonical files live here.

## Conventions

- File names use the `al_<task>` convention for lensing-API skills, e.g. `al_run_search.md`.
- Project-workflow skills (repo-level operations, template manipulation) use a plain
  kebab-case name, e.g. `init-slam.md`, `start-new-project.md`.
- Domain-mode skills that pair an external pipeline repo to the assistant use that
  domain as prefix — `euclid_<task>.md` for the Euclid pipeline (euclid mode).
- Meta-skills (writing guide, bootstrap protocol) start with `_`.
- Every lensing-API skill is **python-first**: the deliverable is a runnable `.py` script
  + the understanding to evolve it. Project-workflow skills may instead drive `rsync`,
  `cp`, or other repo-level operations.
- Source citations use the project-name + repo-relative-path form,
  e.g. `PyAutoFit:autofit/non_linear/search/nest/nautilus/`, resolved via
  [`../sources.yaml`](../sources.yaml).
- Wiki references use workspace-relative paths,
  e.g. `wiki/core/concepts/non_linear_search.md`.

## Index

Every skill below is a **complete recipe** unless marked `(stub)` — the stubs are gathered
under "Pending — stubbed" at the end, with a queue of catalogued-but-unstubbed topics after
them.

### Meta

- [`_style.md`](./_style.md) — writing guide every skill is authored against. Read first
  before adding or editing any skill.
- [`_bootstrap_skill.md`](./_bootstrap_skill.md) — protocol for authoring a new skill on
  demand when a user requests a capability not yet covered.

### Setup & maintenance

- [`al_setup_environment.md`](./al_setup_environment.md) — detect absent or broken PyAuto\*
  environments, install via pip or editable clones when needed, configure caches, verify imports.
- [`al_update_wiki.md`](./al_update_wiki.md) — refresh `wiki/core/` pages whose pinned
  source commits have moved; surface new public APIs for review.
- [`al_audit_skill_apis.md`](./al_audit_skill_apis.md) — verify every PyAuto\* symbol
  cited in `skills/` and `wiki/core/api+stack/` resolves in the installed stack;
  report stale references with suggested replacements.
- [`al_refresh_api_docs.md`](./al_refresh_api_docs.md) — orchestrate a full maintenance
  sweep across skill recipes, wiki API pages, and pinned-source drift after a PyAuto\*
  upgrade or source refresh.
- [`al_ingest_paper.md`](./al_ingest_paper.md) — add a strong-lensing paper (local PDF
  or arxiv URL): project-local `wiki/project/bibliography.md` by default inside a science
  project; shared `wiki/literature/` in the assistant clone or on explicit promotion.

### Project workflow

- [`start-new-project.md`](./start-new-project.md) — the single bridge to a standalone
  **science project** and its full lifecycle (Create → Work → Collaborate → Publish):
  scaffold a lean repo that copies the reproducible science and refers back to the assistant
  for skills/wiki, run modelling with reproducibility manifests + the `wiki/project/` journal,
  build collaborator summaries, and harden for an open-science release (CITATION/license/Zenodo).
  Optional HPC folder.
- [`contribute-upstream.md`](./contribute-upstream.md) — prepare a scoped change,
  push it either to your collaborator branch on `PyAutoLabs/autolens_assistant`
  or to your fork, and open a draft PR into `PyAutoLabs/autolens_assistant`.
- [`init-slam.md`](./init-slam.md) — populate an empty `scripts/` folder with SLaM
  pipeline script(s) copied from `autolens_workspace` and tailored to the chosen data
  type.

### Euclid mode (pipeline-paired)

Skills pairing the collaboration's
[`euclid_strong_lens_modeling_pipeline`](https://github.com/PyAutoLabs/euclid_strong_lens_modeling_pipeline)
to the assistant; the paired literature context is the dedicated
[`wiki/euclid/`](../wiki/euclid/index.md) sub-wiki.

- [`euclid_setup_pipeline.md`](./euclid_setup_pipeline.md) — clone/install the pipeline,
  the `dataset/<sample>/<dataset>/` layout, and the black-box `start_here.py` run.
- [`euclid_prepare_data.md`](./euclid_prepare_data.md) — segmentation validation, binary
  mask tuning, extra-galaxy masks/centres via the GUI tools, hand-assembling a dataset
  folder.
- [`euclid_model_lens.md`](./euclid_model_lens.md) — choose and run the staged
  pipelines: initial MGE+SIE, Sersic photometry, lens-only subtraction, multi-waveband,
  full SLaM.
- [`euclid_workflow_products.md`](./euclid_workflow_products.md) — aggregate many fits
  into .csv catalogues, .fits stacks, and one-line .png summaries via `workflow/`.
- [`euclid_hpc_runs.md`](./euclid_hpc_runs.md) — `hpc/sync` configuration and the
  SLURM batch-array templates for sample-scale runs.

### Data preparation

- [`al_prepare_imaging_data.md`](./al_prepare_imaging_data.md) — load and preprocess
  FITS imaging, decide masking for real data, measure noise, prepare PSF.
- [`al_simulate_dataset.md`](./al_simulate_dataset.md) — synthesise a lens dataset
  (imaging or interferometer) from a ground-truth model.

### Model building

- [`al_build_imaging_model.md`](./al_build_imaging_model.md) — compose a `Tracer` from
  light + mass profiles and wrap it in an `AnalysisImaging`.
- [`al_build_interferometer_model.md`](./al_build_interferometer_model.md) — same, but
  for visibility-plane data.
- [`al_custom_profile.md`](./al_custom_profile.md) — write a new light or mass profile
  subclass and register it for use in models.

### Fitting

- [`al_configure_search.md`](./al_configure_search.md) — pick and tune a non-linear
  search (Nautilus, Dynesty, Emcee, Zeus, …) for your problem.
- [`al_run_search.md`](./al_run_search.md) — execute `search.fit(model=..., analysis=...)`
  and monitor convergence.
- [`al_chain_searches.md`](./al_chain_searches.md) — sequence searches so a later phase
  inherits priors from an earlier one.
- [`al_run_slam_pipeline.md`](./al_run_slam_pipeline.md) — run a Source-Light-Mass
  pipeline (the canonical automated lensing workflow).
- [`al_debug_fit_failure.md`](./al_debug_fit_failure.md) — diagnose a fit that didn't
  converge or produced unphysical results.

### Results & visualisation

- [`al_load_results.md`](./al_load_results.md) — load a completed fit's `Tracer`,
  `Samples`, dataset and FITS products from its output folder.
- [`al_inspect_results_mcp.md`](./al_inspect_results_mcp.md) — the read-only
  results-inspector MCP server: browse fits, summaries, result images and bulk
  subplot/FITS extraction from MCP clients without running Python (Claude
  Desktop first).
- [`al_plot_tracer.md`](./al_plot_tracer.md) — plot ray tracing, critical curves,
  caustics, magnification maps.
- [`al_plot_fit_residuals.md`](./al_plot_fit_residuals.md) — plot model image,
  residuals, normalised residuals, chi-squared map.
- [`al_inspect_source_reconstruction.md`](./al_inspect_source_reconstruction.md) —
  inspect a pixelised inversion: regularisation, source-plane image, reconstruction
  diagnostics.
- [`al_to_notebook.md`](./al_to_notebook.md) — convert a generated narrative-docstring
  script to a Jupyter notebook (docstrings → markdown cells, code → code cells) via the
  stdlib-only `autoassistant/to_notebook.py`.

### Pending — stubbed (need full recipes)

Drafted as scaffolds during the 2026-05-22 coverage audit against
`autolens_workspace/scripts/`. Each has frontmatter + Orient/Ask/Branch/Combine
structure + `Further reading`; the `Branch` recipes are TODO markers. Fill in
one at a time, paired with their companion wiki/core stub.

**Data types and regimes**

- [`al_point_source.md`](./al_point_source.md) (stub) — quasar / multi-image
  position fits, flux ratios, point-source deblending.
- [`al_time_delay_cosmography.md`](./al_time_delay_cosmography.md) (stub) — H0
  from time-delay strong lenses.
- [`al_group_lensing.md`](./al_group_lensing.md) (stub) — extra galaxies,
  scaling-relation members.
- [`al_cluster_csv_api.md`](./al_cluster_csv_api.md) (stub) — cluster-scale
  CSV-driven model composition.
- [`al_multi_dataset.md`](./al_multi_dataset.md) (stub) — joint imaging +
  interferometer, multi-band, wavelength-dependent sources.
- [`al_weak_lensing.md`](./al_weak_lensing.md) (stub) — shear catalogue fits
  (`WeakDataset` / `AnalysisWeak`).
- [`al_datacube_modeling.md`](./al_datacube_modeling.md) (stub) — interferometer
  spectral cubes.

**Dark-matter substructure**

- [`al_subhalo_detect.md`](./al_subhalo_detect.md) (stub) — Bayesian-evidence
  grid search for perturbing subhaloes.
- [`al_sensitivity_mapping.md`](./al_sensitivity_mapping.md) (stub) — quantitative
  detectability calibration.
- [`al_potential_correction.md`](./al_potential_correction.md) — gravitational
  imaging: pixelized corrections to the lensing potential (`al.pc`), reconstructed
  jointly with the source, for substructure whose form you do not want to assume.

**Advanced techniques**

- [`al_hierarchical_inference.md`](./al_hierarchical_inference.md) (stub) —
  population-level / graphical models, expectation propagation.
- [`al_aggregator_bulk_analysis.md`](./al_aggregator_bulk_analysis.md) (stub) —
  bulk operations across many completed fits, optional result database.
- [`al_adaptive_pixelization.md`](./al_adaptive_pixelization.md) (stub) — adaptive
  mesh + adaptive regularisation source reconstructions.
- [`al_mge_decomposition.md`](./al_mge_decomposition.md) (stub) — Multi-Gaussian
  Expansion workflows for lens / source.
- [`al_custom_analysis.md`](./al_custom_analysis.md) (stub) — subclassing
  `Analysis` to add custom likelihood terms.

**Queue (catalogued, not yet stubbed):** `al_multi_plane`, `al_los_halos`,
`al_over_sampling`, `al_workflow_outputs`, `al_data_prep_interactive`,
`al_bayesian_model_comparison`.
