# OpenCode (coding agent — experimental alternative)

[OpenCode](https://github.com/sst/opencode) is an **open-source** coding agent. The client is
free; you connect it to a model provider of your choice. It reads `AGENTS.md` from the
repository root, so the assistant's instructions, skills and wiki all apply, and it can install
PyAutoLens, run fits and inspect results **if the model you connect can drive the workflow**.

This page is here because OpenCode is the route people ask about when they have no paid Claude
Code or Codex access. Read the caveats before spending time on it.

## What "free" does and does not mean

- **The client is free; the model is not OpenCode's.** Cost, capability, rate limits and
  availability all belong to the provider you configure. OpenCode's own docs recommend new users
  start with its hosted provider (OpenCode Zen), which is pay-as-you-go with a handful of models
  offered at no cost — and those free models are described by OpenCode as available *for a
  limited time* while feedback is collected. Do not plan a project around a free model staying
  free.
- **Not every model can run the assistant.** The workflow needs a model that (a) follows
  multi-step tool use reliably over a long session, (b) does not hallucinate PyAutoLens API from
  training data faster than the self-enforced audit can catch it, and (c) can **look at figures**
  — the real-data gate and every fit check involve inspecting a PNG. Many free models have no
  vision input, which rules them out for the full workflow even if they can write code.
- **Nothing here is validated yet.** As of 2026-09-10 **no provider/model configuration has been
  tested against this assistant's benchmarks**, so there is no default free model to recommend.
  When one is validated it will be recorded here with its test date, the model and provider, and
  its limitations; until then treat OpenCode as *compatible* (it loads the instructions and runs
  code) rather than *supported*. The evaluation protocol maintainers use is in
  [`../evaluation/agent_evaluation.md`](../evaluation/agent_evaluation.md).

## Setup

1. Install OpenCode and configure a model provider — follow the official instructions at
   [github.com/sst/opencode](https://github.com/sst/opencode). Pick a model that supports
   **tool use and image input**; check the provider's model card, not the model's name.
2. Clone this repository and start the agent inside it (running from the repository root is
   what lets it discover the assistant's instructions in `AGENTS.md`):

```bash
git clone https://github.com/PyAutoLabs/autolens_assistant.git
cd autolens_assistant
opencode
```

3. OpenCode has no PreToolUse hook, so the code gate is self-enforced: the assistant should run
   `python autoassistant/audit_skill_apis.py --file <script.py>` on generated PyAutoLens code
   before executing it. If the model skips this, ask it to.

## Your first prompt

```text
Find the data on the Cosmos-Web ring, give me a short script to plot it in PyAutoLens
and then given that I'm a new user give me an overview of the different ways we can
perform strong lens modeling of this system.
```

Then a quick capability check before trusting it with real work: ask it to plot the dataset
and *describe what it sees in the image*. A model that cannot view the figure will describe the
code instead — that is the signal it cannot run the full workflow.

## Reporting back

If you run the assistant on OpenCode, please
[open an issue](https://github.com/PyAutoLabs/autolens_assistant/issues) with the provider,
model, date, and what worked or failed (the three checks in the evaluation protocol are a good
template). User reports are how a configuration graduates from "compatible" to "tested".
