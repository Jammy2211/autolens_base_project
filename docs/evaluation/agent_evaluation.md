# Evaluating coding agents against the assistant (maintainers)

This page is the maintainer-facing protocol for deciding what the README may claim about an
agent. It exists because the support statements are now deliberately narrow — **Claude Code
and Codex recommended; OpenCode experimental; nothing else supported** — and every widening of
that list should be backed by a recorded run, not an impression.

Vocabulary used throughout the assistants:

| Word | Meaning |
|---|---|
| **Recommended / tested** | The agent the assistant is developed against day to day, and on which the three checks below pass repeatedly. Claude Code and Codex. |
| **Experimental / compatible** | The agent loads `AGENTS.md`, runs code and follows the skills in principle, but no configuration has passed the checks below on record. OpenCode today. |
| **Candidate for evaluation** | Worth a look; not offered to users at all until it has a record. Gemini CLI today. |
| **Unsupported** | Cannot meet the requirements (no execution, no data access) — the retired chat routes. |

"Compatible" is not a performance claim. Do not write "works brilliantly", "reliable" or
"benchmarked" anywhere unless `benchmarks/RESULTS.md` carries the run that justifies it; as of
2026-09-10 it carries none.

## The three checks

Small enough to run in an afternoon on a laptop, large enough to expose the failure modes that
matter. They are shipped as the frozen benchmark card
[`benchmarks/prompts/harness_smoke.md`](../../benchmarks/prompts/harness_smoke.md) so a run is
recorded exactly like any other benchmark (`python autoassistant/benchmark.py new-run
harness-smoke --model <model> --harness <harness>`, then score and `report`).

1. **Grounded answering.** Ask a scientific/API question whose answer lives in a specific
   repository file (the card names one). Pass = the agent opens the right `skills/` or
   `wiki/core/` file, cites it, and the answer matches; fail = an answer from memory, or a
   stale API symbol.
2. **A small fit with inspected outputs.** Fit the bundled COSMOS-Web Ring with the assistant's
   own quick workflow under `PYAUTO_TEST_MODE=1` (reduced iterations, minutes not hours). Pass =
   the real-data gate is honoured *before* the fit (data plotted, contaminants and mask settled
   with the operator), the search completes, and the agent **opens the fit subplot and describes
   what it shows** — this is the figure-inspection test; an agent that describes the code
   instead of the image cannot run the full workflow.
3. **Recovery from an API or execution error.** The operator introduces one stale symbol
   (e.g. `aplt.FitImagingPlotter`) into a script and asks the agent to run it. Pass = the agent
   diagnoses it against the installed library or `skills/`, fixes it without inventing another
   symbol, and re-runs; on harnesses without the PreToolUse hook it should run
   `audit_skill_apis.py --file` before re-executing.

Record failures as scored runs too — a negative result is the evidence that keeps an agent off
the supported list.

## Status of the outstanding evaluations (2026-09-10)

None of these runs could be performed in the session that adopted the agentic-only policy (a
remote container with no access to the agents themselves), so the table below is the gap, not
the result. Nothing in it should be copied into user-facing docs as a support claim.

| Agent / configuration | Access verified from official docs? | Runs recorded | Notes |
|---|---|---|---|
| Claude Code | Yes — subscription (Pro/Max, Team/Enterprise seats) or Claude Console usage billing, or Bedrock / Google Cloud / Foundry ([costs page](https://code.claude.com/docs/en/costs)). No free tier described. | 0 scored (developed against daily) | Primary harness; the only one with the PreToolUse code gate. |
| Codex | Partially — the Codex README states sign-in with Plus, Pro, Business, Edu or Enterprise, or an API key. The pricing page (`developers.openai.com/codex/pricing`) was unreachable from the adopting session, so **free-plan Codex access is unverified**. | 0 scored (developed against daily) | Re-check the pricing page before stating anything about free access. |
| OpenCode + a named free model | OpenCode's own docs: client open source; hosted provider (Zen) is pay-as-you-go with a few no-cost models "available for a limited time"; only one Zen model is documented as vision-capable, and it is not one of the free ones. | 0 | **No default free configuration is validated.** A candidate must pass check 2's figure inspection, which needs image input. |
| Gemini CLI | Conflicting: the assistant's setup page recorded on 2026-06-18 that free/individual access was retired, while the Gemini CLI quota document currently lists a free individual tier (1,000 requests/day) alongside Google AI Pro/Ultra, Code Assist Standard/Enterprise, API-key and Vertex routes. Resolve at evaluation time. | 0 | Candidate only. `.gemini/settings.json` is kept so it loads `AGENTS.md`; do not add a user-facing setup page until it has a record. |

Promotion rule: an agent (or an OpenCode provider/model pair) moves to *experimental with a
tested configuration* after one recorded pass of all three checks, and to *recommended* only
after repeated passes across dates plus at least one of the full science benchmark cards. Every
promotion is a dated edit to this table and to the README.
