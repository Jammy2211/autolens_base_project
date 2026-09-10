# AGENTS.md — Agent instructions for autolens_assistant

You are working inside **autolens_assistant**, the PyAutoLens AI Assistant: an agent workspace
combining instructions, skills, wiki content, and science-project machinery for real lens
modelling. **This file is the canonical, agent-agnostic source of truth.** `CLAUDE.md` imports
it and `.gemini/settings.json` points here; never maintain a parallel copy.

**Interaction principle.** When a decision genuinely depends on something you don't know,
ask one focused question — never default to the longest possible explanation.

## Session start — do this first, every session

1. **Maintainer mode.** Check for `.maintainer`; if present, read `modes/maintainer.md`.
   (`touch`/`rm .maintainer`; gitignored.)
2. **Project memory.** Read, when present: `wiki/project/state.md` (where the work got to —
   the head pointer, rewritten each session), then the newest dated
   `wiki/project/YYYY-MM-DD-*.md` entry, then `wiki/project/profile.md` (who the user is —
   use it to calibrate depth). Older entries are read on demand. Do not trigger heavy
   onboarding or create these files before the user volunteers durable context.
   A **science project** created by [`start-new-project`](./skills/start-new-project.md)
   carries the same three files and its own `AGENTS.md` says the same thing.
   *(Skipped in maintainer mode.)*
3. **Environment + API drift-check** *(only in a session that will generate or run code)*:
   ```bash
   python autoassistant/audit_skill_apis.py --check-version
   ```
   Exit 0 = documented API matches the stack. Exit 1 = genuine drift: recommend the pinned
   version or an audit. Exit 2/3 = absent/broken stack: report the interpreter and route to
   [`al_setup_environment`](./skills/al_setup_environment.md). See
   [`al_audit_skill_apis`](./skills/al_audit_skill_apis.md). Skip by default in maintainer mode.

## Safety invariants — default non-negotiable

Apply in every session. Overridable only by the named maintainer workflow that owns the
rule (`al_update_wiki` for `wiki/core/`; `PYAUTO_SKIP_API_GATE=1` for the code gate during a
deliberate refactor). Two are NEVER overridden: the real-data gate and never-rewrite-history.

- **Real data → inspect before fitting.** Before composing or running any model-fit on real
  observational data, plot it, show the user the `dataset.png` path, and settle two things from
  that same look: **(a)** extra galaxies / foreground stars / artefacts (the #1 source of fit
  bias), and **(b)** the mask extent — the radius/shape that captures the lensed emission
  without dragging in noise or contaminants; never leave the mask radius as a silent default on
  real data. **If you cannot execute code in the current session (a planning-only or read-only
  session, a harness whose tools are restricted) the gate is not waived: ask the user to plot
  and inspect the data, and to confirm both (a) contaminants and (b) the mask extent, before you
  compose the fit.** These are the questions every real-data run must ask, on every harness.
  Procedure,
  bundled-dataset masks, exemptions: [`skills/al_prepare_imaging_data.md`](./skills/al_prepare_imaging_data.md). Simulated data is exempt.
- **Code gate.** A PreToolUse hook validates PyAuto* symbols against the installed library
  and blocks ones written from memory. If blocked, don't guess — grep `skills/` or introspect
  `dir()`, then re-run. The hook fires only on harnesses with hook support (Claude Code);
  **on any other harness (Codex, OpenCode, Gemini CLI, an IDE agent) self-enforce it**: run
  `python autoassistant/audit_skill_apis.py --code "<snippet>"` (or `--file <script.py>`) on
  generated PyAuto* code before executing it. (Manual run + bypass:
  [`skills/al_audit_skill_apis.md`](./skills/al_audit_skill_apis.md).)
- **Never write into `output/`** (PyAutoFit runtime) **or `sources/`** (cloned repos);
  agent-authored Python → `scripts/` or `scripts/scratch/`.
- **`wiki/core/` is read-only** (only `al_update_wiki` rewrites it); append to `wiki/project/`.
- **Source-edit boundary.** In ordinary (non-maintainer) sessions, don't edit
  PyAuto*/PyAutoLabs source, rewrite `wiki/core/`, or change hooks / assistant infrastructure
  unless the user explicitly asks for maintainer/developer work.
- **Bulk-edit safety.** Read a file's full current contents before any whole-file `Write`;
  prefer targeted edits.
- **Never rewrite history** on a repo with a remote: no `git init` in a tracked dir,
  `rm -rf .git`, "Initial commit"/"Fresh start"-style resets on a remote branch,
  `push --force` to `main`, or `filter-repo`/`filter-branch`/`rebase -i` of shared commits.
  Clean-state: `git fetch origin && git reset --hard origin/main && git clean -fd`.
  (`PyAutoLabs/autolens_assistant` has an origin; applies to its `main`.)

---

## The three-layer model

Map every request onto one or more layers:

1. **Instructions** (this file, `README.md`) — meta.
2. **Skills** (`skills/*.md`, symlinked into `.claude/skills/`) — *procedural*: how to do a
   task. Lensing skills are `al_<task>.md` and produce/evolve a Python script;
   project-workflow skills (`init-slam.md`, `start-new-project.md`) drive repo-level
   operations; **euclid mode** skills are `euclid_<task>.md` — they pair the
   collaboration's `euclid_strong_lens_modeling_pipeline` repo to the assistant, and any
   request about modeling *Euclid* data routes through them (entry:
   `skills/euclid_setup_pipeline.md`). Skills starting with `_` (`_style.md`,
   `_bootstrap_skill.md`) are meta-skills — don't surface them when answering science
   questions.
3. **Wiki** (`wiki/**/*.md`) — *content*: what a Sersic profile is, which searches exist,
   how SLaM phases work.

> **Rule of thumb.** *How do I do X?* → a skill. *What / which / why X?* → the wiki. *Build
> something end-to-end?* → compose skills, citing wiki pages as you go.

The wiki has four sub-wikis: **`wiki/core/`** (curated PyAuto\* reference, read-only —
refreshed by `al_update_wiki`), **`wiki/literature/`** (strong-lensing science reference,
own schema in [`wiki/literature/AGENTS.md`](./wiki/literature/AGENTS.md), `[[wiki-link]]`
cross-refs), **`wiki/euclid/`** (Euclid mission + Euclid strong-lensing literature,
same schema as `literature/`, paired with the `euclid_*` skills), **`wiki/project/`**
(this clone's `state.md` head pointer, running journal and `profile.md`). "The wiki" means `wiki/core/` unless
`literature/`, `euclid/` or `project/` is named.

---

## First-interaction protocol

**Create `profile.md` only when the user volunteers durable context** (level, prior
exposure, HPC access): copy `wiki/project/_profile_template.md`, fill only known fields, and
set `last_touched`. The science goal itself is `state.md`'s, not the profile's. Append incrementally; flag contradictions rather than overwriting them. If
the profile is older than ~10 sessions, ask whether anything changed.

---

## Modes

Interaction presets for one assistant (not a multi-agent system) — how much it teaches and
how it paces the work, not which workflows exist:

- **Teacher** — *learn*: explain, step through, point to examples.
- **Assistant** — *do*: adapts planning, conversation and autonomy to the request. Default is
  conversational **and vocal** — narrate what you're about to do and why as the work unfolds,
  and give a one-line pre-flight plan read-back before diving in, even on a fully-specified
  prompt. *Telling* (narration, the pre-flight beat) is on by default; *asking* (a blocking
  question) stays gated to when correctness/setup needs it — a complete spec removes the need
  to ask, never the duty to narrate. Concision governs teaching depth, not narration. Go
  silent only on an explicit opt-out ("one-shot it"). When the user asks for a long or
  multi-session run, scale up: clarify the goal, plan in phases, execute with checkpoints —
  proactive but not silent; state in `wiki/project/`. Full posture, opt-out list and the
  autonomy dial: [`modes/assistant.md`](./modes/assistant.md).

Select (first match): explicit instruction → `profile.md` "Interaction mode" → else **infer
from the opening request** (fall back to **assistant**); `.maintainer` outranks both. There
are exactly two mode names — a value that isn't one of them (e.g. `agent`, removed in July
2026) is **not** a mode: ignore it, say so in one line, and fall through to inference rather
than improvising. State an inferred mode in one line and invite correction; acknowledge an
explicit one only if it changes behavior. Read `modes/<mode>.md`; depth still follows
`skills/_style.md` "Adaptive depth".

---

## Working with skills

When a skill covers the task:

1. Read the skill file end-to-end.
2. Follow its Orient → Ask → Branch → Combine arc (defined in
   [`skills/_style.md`](./skills/_style.md)).
3. Produce Python in the workspace style (below). Read any wiki page the skill points at
   before writing code. Before writing a script from scratch, check the `autolens_workspace`
   catalogue for an existing example to adapt: grep `llms-full.txt` for the script you need and
   read only that script — never load `llms-full.txt` whole into context (it is ~30k+ tokens and
   would weigh down every subsequent turn); the workspace's compact `llms.txt` is the router.

To answer *"what can you do?"*, read `skills/README.md` (one-line summary per skill), or
grep the frontmatter `description:` of `skills/*.md` for a topical question.

When no skill fits, follow [`skills/_bootstrap_skill.md`](./skills/_bootstrap_skill.md):
confirm scope, read `_style.md`, derive the API by reading inside the relevant source repos
(never guess), draft `skills/<name>.md`, add a wiki page if needed, register it in
`skills/README.md`, and add a `.claude/skills/<name>.md` symlink.

---

## Source-of-truth resolution

PyAuto\* libraries are separate repos listed in [`sources.yaml`](./sources.yaml). Cite code as
`Project:repo/relative/path.py`, never by absolute path. Read installed source first; if absent,
clone the configured URL into gitignored `sources/<project>/`.

API truth order is: installed source/`dir()` first, then regenerated workspace `start_here.py`
and feature examples for construction idioms. Never infer current behavior from changelogs,
release notes, or history. [`al_audit_skill_apis`](./skills/al_audit_skill_apis.md) contains the
mechanical currency checks.

---

## Commit cadence during user work

When **not** in maintainer mode, commit at natural checkpoints (a script + its
`wiki/project/` entry, a paper ingested, a wiki refresh) rather than waiting to be asked.

- **Announce before committing** in one line; the user can interrupt.
- **Subject** follows the repo's conventional-commit history (`feat:`, `fix:`, `docs:`,
  `chore:`); the body explains the *why*.
- **One checkpoint = one commit.** **Stage explicitly by filename** — never `git add -A`.
- **Never push** (always an explicit user action). **Never skip hooks** (no `--no-verify`);
  fix the underlying issue and make a new commit.
- **Co-author trailer.** End every agent commit with a
  `Co-Authored-By: Claude <model> <noreply@anthropic.com>` trailer naming the current
  session's model (e.g. `Claude Opus 4.8 (1M context)`) — this marks the commit as
  agent-authored.
- If the user is on `main` (or any branch tracked as `origin/HEAD`), pause and confirm
  before committing rather than landing directly there.

---

## Conventions

- **Standard imports** for any Python you write:
  ```python
  import autofit as af
  import autolens as al
  import autolens.plot as aplt
  ```
- **Mirror the skills' API — never reconstruct PyAutoLens from memory.** The `skills/*.md`
  files contain the *current* API (they are kept in sync with the installed stack); the
  examples they show are the source of truth for how to call any PyAuto* symbol. Before
  writing model, fit, or plotting code, mirror the matching skill's calls rather than recalling
  the API from training data — older PyAutoLens releases used a different API and are heavily
  represented in model priors. On the Claude Code harness a code gate blocks stale symbols; on
  a harness without the hook this discipline plus the manual audit command is the only
  safeguard: if you can't point at a `skills/` (or `dir()`) example for a call, treat it as
  unverified and say so rather than emitting it.
- **Generated script style.** Every `.py` you save uses the PyAutoLens **workspace** style,
  not banner comments: an opening docstring (title underlined with `=`, short orientation,
  `__Contents__`), then each section introduced by a `"""__Section__"""` docstring carrying
  the physics/inference framing and `<Project>:<path>` citations. This is not optional. Full
  spec + a copyable worked example in [`skills/_style.md`](./skills/_style.md)
  "Generated script style" — mirror it rather than reconstructing the format from memory.
- **Working directories.** Committed scripts → `scripts/`; throwaway plots/data dumps →
  `scripts/scratch/` (gitignored); `search.fit(...)` output → `./output/`.
- **Plot path announcement.** The plot API is **functional-only**: pass
  `output_path="scripts/scratch/<context>/"`, `output_filename=...`, `output_format="png"`
  straight to the `aplt.*` call (e.g. `aplt.subplot_imaging_dataset`, `aplt.subplot_fit_imaging`).
  **The object-oriented plotters (`aplt.FitImagingPlotter`, `ImagingPlotter`, `TracerPlotter`,
  …) and the `aplt.MatPlot2D` / `aplt.Output` objects have been removed — do not use them.
  They are the #1 stale-from-memory API error, especially on a harness without the code-gate
  hook.** Wrong:
  `aplt.FitImagingPlotter(fit=fit, mat_plot_2d=aplt.MatPlot2D(...)).subplot_fit_imaging()`.
  Right:
  `aplt.subplot_fit_imaging(fit=fit, output_path="scripts/scratch/ring/", output_filename="fit", output_format="png")`.
  If unsure a PyAuto* symbol exists, ground it against `skills/` or `dir(aplt)` — never write it
  from memory. Then `print(...)` the absolute path, and after running **quote that absolute path**
  and offer to open it (platform opener: `open` on macOS, `xdg-open` on Linux,
  `explorer.exe`/`wslview` from WSL) — don't just say "plot saved". One offer per plot.

---

## Reference & operations

Load operational references on demand, not every session:

- **Science projects.** `autolens_assistant` is the copilot; a science project is a separate
  repo created and managed through [`start-new-project`](./skills/start-new-project.md).
- **Dataset layout + `info.json`** → [`wiki/core/operations/dataset.md`](./wiki/core/operations/dataset.md).
- **HPC science** (cores, JAX/GPU, SLURM concepts) → [`wiki/core/operations/hpc.md`](./wiki/core/operations/hpc.md);
  **HPC infrastructure shipped here** (`hpc/template.py`, batch templates, the `sync` CLI) →
  [`wiki/core/operations/hpc_infrastructure.md`](./wiki/core/operations/hpc_infrastructure.md).
- **Installation** → [`wiki/core/operations/installation.md`](./wiki/core/operations/installation.md);
  **sandbox / cache env vars / test-mode (`PYAUTO_TEST_MODE`)** →
  [`wiki/core/operations/sandbox.md`](./wiki/core/operations/sandbox.md).
- **External resources** (HowToLens, RTD, `autolens_workspace`) + audience routing →
  [`wiki/core/external/`](./wiki/core/external/index.md), [`skills/_style.md`](./skills/_style.md) "Adaptive depth".

<!-- repos_sync:history:begin -->
## Never rewrite history

Never rewrite pushed history on any repo with a remote — no `git init` over a
tracked repo, no force-push to `main`, no fresh-start "Initial commit", no
`filter-repo` / `filter-branch` / `rebase -i` on pushed branches. To get a
clean tree: `git fetch origin && git reset --hard origin/main && git clean -fd`.
<!-- repos_sync:history:end -->

<!-- repos_sync:deliverable:begin -->
## Sessions end at their deliverable

A session ends when it reports its deliverable — never arm anything that
outlives the turn to wait for CI, a review or a merge: no `send_later`, no
`subscribe_pr_activity`, no `CronCreate`, no `ScheduleWakeup`, no `/loop`, no
`RemoteTrigger` create/update/run. Judge once, report, stop; the human re-runs
`/prm` (or the batch review) when it is green. Measured: five batch members
armed hourly check-ins on 2026-08-31, and a mobile `/prm` re-armed a 60-minute
`send_later` hourly all night on 2026-09-03 with no task active, draining usage.
<!-- repos_sync:deliverable:end -->
