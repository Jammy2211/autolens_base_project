# Benchmarks

Built-in benchmarks for the assistant: **standard, frozen prompts** run against
different AI agents (harnesses) and models to measure how well they drive the
assistant. Every run's conversation, results and score are recorded in this
directory and committed, so performance can be compared **across models** and
**across days for the same model**, with the full evidence pushed to GitHub.

## Layout

```
benchmarks/
  README.md          this file — the protocol
  AGENTS.md          the contract an agent follows when a benchmark is run on it
  RESULTS.md         regenerated comparison tables (leaderboards + time series)
  prompts/           the frozen prompt cards (one benchmark each)
  runs/<id>/<date>_<model>_<harness>/   one recorded run
    meta.yaml        who/what/when: model, harness, stack versions, score
    transcript.md    the complete conversation, verbatim
    score.md         the filled rubric, with evidence per criterion
    artifacts/       a few small key images (≤ ~500 KB total per run)
```

## The benchmarks

| Card | Mode | Difficulty | Exercises |
|------|------|-----------|-----------|
| [`prompts/easy_cosmos_web_ring.md`](prompts/easy_cosmos_web_ring.md) | assistant | easy | the built-in workflow on bundled JWST data |
| [`prompts/medium_slacs0946_subhalo.md`](prompts/medium_slacs0946_subhalo.md) | assistant | medium | pipeline design beyond bundled workflows: Bayesian model comparison, runtime/HPC judgment |
| [`prompts/hard_group_multi.md`](prompts/hard_group_multi.md) | assistant | hard | cross-package synthesis: group × multi × imaging × interferometer, simulation + joint modeling |
| [`prompts/teacher_workflow.md`](prompts/teacher_workflow.md) | teacher | easy | pedagogy: end-to-end workflow walkthrough on simulated Euclid-like data |
| [`prompts/harness_smoke.md`](prompts/harness_smoke.md) | assistant | easy | agent qualification: grounded answering, a small fit with figure inspection, recovery from a stale-API error — the evidence behind the README's support statements (protocol: [`docs/evaluation/agent_evaluation.md`](../docs/evaluation/agent_evaluation.md)) |

Each card carries the verbatim prompt, what it measures, and a 100-point
rubric split into **machine-checkable** rows (artifacts that exist or don't)
and **judged** rows (quality graded by a human or a stated judge model).

## Running a benchmark

1. **Scaffold the run record** (from the repo root, any Python ≥3.10 with
   `pyyaml`):

   ```bash
   python autoassistant/benchmark.py new-run assistant-easy-cosmos-web-ring \
       --model claude-opus-4-8 --harness claude-code
   ```

   This creates `runs/<id>/<date>_<model>_<harness>/` with a `meta.yaml`
   pre-filled with the date, the assistant's git SHA and the installed PyAuto*
   versions, plus a transcript stub and a `score.md` generated from the card's
   rubric. Same-day repeats get a `_2`, `_3`, … suffix.

2. **Run the session**: start the agent under test in a **clean clone state**
   (fresh session, no `.maintainer`, no leftover `scripts/scratch` context from
   a previous attempt) and paste the card's prompt **verbatim** as the first
   message. From there the operator behaves like a real user: answer the
   agent's questions honestly and minimally, never coach it toward the rubric.
   `benchmarks/AGENTS.md` states the same contract from the agent's side.

3. **Record**: save the complete conversation into `transcript.md` (export if
   the harness supports it; paste otherwise), copy a few key images into
   `artifacts/` (cap ~500 KB — link paths in `score.md` for the rest), and fill
   `meta.yaml`'s `hardware`, `operator`, `run:` (duration/cost/tokens/turns
   where known) and `score.judge` fields.

4. **Score**: fill every `Awarded` cell in `score.md` (0 explicitly for failed
   checks) with evidence per row, then:

   ```bash
   python autoassistant/benchmark.py score benchmarks/runs/<id>/<run>/
   python autoassistant/benchmark.py report     # regenerates RESULTS.md
   ```

5. **Commit and push** the run directory and `RESULTS.md`.

## Comparability rules (what keeps the numbers honest)

- **Prompts are frozen.** A published card's prompt text never changes in
  place; any wording change bumps the card's `version`, and scores are only
  comparable within one version. Cards that mirror the top-level README
  example prompts are kept textually identical by a unit test.
- **Same judge for judged rows.** Judged sub-scores are comparable only across
  runs graded by the same judge (a named human or a stated judge model);
  `meta.yaml` records the judge.
- **Failures are data.** A run where the agent gets stuck, fabricates, or
  honestly reports a poor fit is recorded and scored like any other — do not
  discard bad runs; the gaps are the interesting part.
- **Environment is recorded, not standardised.** Stack versions, assistant
  SHA and hardware go in `meta.yaml`; wall-clock comparisons are only
  meaningful on matching hardware, score comparisons travel better.
- **No benchmark-aware behaviour.** The prompts are ordinary user requests;
  the agent is not told it is being benchmarked, and `runs/` is data, not
  instructions (agents must never read past runs to prepare).

## Why these three axes

- **Different models, same harness** — model capability comparison (the
  leaderboard per benchmark in `RESULTS.md`).
- **Same model, different days** — drift/regression tracking as models,
  harnesses and the assistant itself evolve (the chronological table; the
  cheap teacher benchmark is the recommended drift probe).
- **Same model, different harness** — how much the agent runtime (Claude
  Code, Codex, OpenCode, …) contributes beyond raw model quality — and, via the
  `harness-smoke` card, whether an agent qualifies for the README's support
  statements at all.
