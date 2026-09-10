# Troubleshooting (coding agents)

Failure modes we know about when running the assistant inside a coding agent. If you hit
something not listed here, [open a GitHub issue](https://github.com/PyAutoLabs/autolens_assistant/issues)
— agents and models change fast and user reports are how this page stays current.

**It wrote `aplt.FitImagingPlotter(...)` or `aplt.MatPlot2D(...)`.**
That is the stale-API failure — those classes were removed and older PyAutoLens is heavily
represented in model training data. On Claude Code the PreToolUse code gate blocks it; on any
other agent reply: *"Plotting is functional now — re-check the plotting skill and rewrite using
`aplt.subplot_fit_imaging(...)`, then run `python autoassistant/audit_skill_apis.py --file
<script>` before executing."* If it keeps happening, the agent has stopped reading `skills/`
— start a fresh session in the repository root so `AGENTS.md` loads again.

**It answered a PyAutoLens API question without opening any file.**
Ask it to name the `skills/` or `wiki/core/` page it read. The assistant's value is that it
grounds answers in the repository and the installed library; an answer from memory is the
failure mode the whole design exists to prevent.

**It started composing a fit on real data without plotting it.**
The real-data gate is not optional: the assistant must plot the dataset, show you the path, and
settle contaminants and the mask extent with you first. Tell it to stop and do that; if it
happens repeatedly on a non-Claude-Code agent, please report the agent and model.

**`audit_skill_apis.py --check-version` exits 1 or 2/3 at session start.**
Exit 1 is genuine API drift between the documented and installed stack — follow the
`al_audit_skill_apis` skill (usually: install the pinned version, or run the audit). Exit 2/3
means PyAutoLens is absent or broken in the interpreter the agent is using; route to
`al_setup_environment` (the agent will offer to install it).

**Cache or permission errors from numba / matplotlib in a sandboxed agent.**
Prefix runs with `NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/matplotlib`. Details in
`wiki/core/operations/sandbox.md`.

**It ran out of context / got slow and vague.**
Something large was pulled in — `wiki/literature/` and `llms-full.txt` are big enough to crowd
out a session on their own. Start a fresh session and tell it to read single pages only (grep
`llms-full.txt`, never read it whole).

**It stopped mid-task, or the agent switched to a weaker model.**
You hit the plan's usage window. On a subscription the agent tells you when it resets; on
usage-based billing check your spend controls. Long modelling sessions are where paid access
earns its keep — see the access notes on [Claude Code](claude_code.md) and [Codex](codex_cli.md).

**On OpenCode the model cannot see figures, or ignores the audit step.**
That is a model limitation, not a bug in the assistant: see
[OpenCode — what "free" does and does not mean](opencode_cli.md). Switch to a model with image
input and reliable tool use, and report the configuration you tried.

**You were using a browser chat (ChatGPT / Claude with a GitHub connector, the custom GPT, a
pasted bundle).**
Those routes are no longer supported (retired 2026-09-10). Their old pages are archived under
[`docs/archive/`](../archive/README.md); the supported path is a coding agent, above.
