---
id: harness-smoke
version: 1
mode: assistant
difficulty: easy
datasets:
  - dataset/imaging/cosmos_web_ring
workspace_packages:
  - imaging
added: 2026-09-10
---

# Benchmark: harness smoke (assistant · easy)

A short, operator-driven card whose purpose is **to qualify an agent or model**, not to test
lensing skill: can this harness ground its answers in the repository, run a small fit and look
at the result, and recover from a stale-API error? It is the evidence behind every support
statement in the README; the protocol and the promotion rule are in
[`docs/evaluation/agent_evaluation.md`](../../docs/evaluation/agent_evaluation.md). Unlike the
science cards it is a three-turn conversation, so the operator's three messages are all part of
the frozen prompt.

## Prompt

Three messages, pasted verbatim in order as a fresh session progresses (see
[`../AGENTS.md`](../AGENTS.md) for the run protocol). Wait for the agent to finish each before
sending the next.

```
Message 1 — grounded answering:
Which non-linear searches does the assistant recommend for a first lens model of the
bundled COSMOS-Web Ring, and where in this repository is that documented? Point me at
the exact file(s) you read before answering.

Message 2 — a small fit with inspected outputs:
Now fit the F277W imaging in dataset/imaging/cosmos_web_ring with the simplest sensible
lens light, mass and source model, running with PYAUTO_TEST_MODE=1 so it finishes in
minutes. When it is done, open the fit subplot and tell me what you see in the image.

Message 3 — recovery from an API error:
I edited your script to plot with aplt.FitImagingPlotter(fit=fit).subplot_fit_imaging().
Run it again and sort out whatever happens.
```

Before sending message 3 the operator makes exactly that edit to the agent's script.

## What this measures

- Grounding: repository files read and cited rather than answered from memory.
- The real-data gate and figure inspection: the data is plotted and discussed **before** the
  fit, and the fit image is actually looked at afterwards (a model without image input cannot
  pass this, which is the point).
- Error recovery without hallucination: the stale plotter is diagnosed against the installed
  library or `skills/`, fixed with the functional API, and the script re-run.

## Success rubric (100 points)

### Machine-checkable (50)

| # | Check | Pts |
|---|-------|-----|
| M1 | Message 1 answer names at least one real repository file that documents search choice (e.g. `skills/al_configure_search.md`), and that file exists | 10 |
| M2 | A script under `scripts/` performing the message-2 fit exists and uses only functional `aplt.*` plotting | 10 |
| M3 | A completed search result exists under `output/` for the message-2 fit (test mode is allowed for this card) | 10 |
| M4 | A dataset plot was written **before** the fit was run (file timestamp precedes the search output) | 10 |
| M5 | After message 3 the script no longer contains `FitImagingPlotter` / `MatPlot2D` and was re-run successfully | 10 |

### Judged (50)

| # | Criterion | Pts |
|---|-----------|-----|
| J1 | Message 1 answer is correct for the current stack and consistent with the cited file; no stale symbols | 10 |
| J2 | Real-data gate honoured: contaminants and mask extent raised with the operator before the fit | 10 |
| J3 | Figure inspection is genuine: the description of the fit subplot refers to what is visible (residual structure, source, arc) rather than to the code | 15 |
| J4 | Error recovery: diagnosis names the real cause (removed object-oriented plotters), fix mirrors the plotting skill, no new invented API | 10 |
| J5 | Conduct: concise, no fabricated numbers, self-enforced audit run on harnesses without the code-gate hook | 5 |

## Operator notes

- Expected wall-clock: 10–30 minutes. `PYAUTO_TEST_MODE=1` keeps the search short; do not use
  `=2` or `=3`, which skip the fit output the card checks.
- Record the model, harness, provider (for OpenCode) and date in `meta.yaml` exactly — this card
  exists to be compared across agents.
- A fail on J3 (cannot see the figure) is the single most useful negative result to record: it
  disqualifies a configuration from the full workflow regardless of its other scores.
