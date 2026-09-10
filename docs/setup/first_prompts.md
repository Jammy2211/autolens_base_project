# First prompts to try

Once [Claude Code](claude_code.md), [Codex](codex_cli.md) or (experimentally)
[OpenCode](opencode_cli.md) is open inside the repository, these all work immediately — the
COSMOS-Web Ring data ships with the repository:

<sub><b>Example Natural Language Prompt for Claude Code, Codex or other AI coding agent</b></sub>

```text
Find the data on the COSMOS-Web ring, give me a short script to plot it in PyAutoLens,
and then, given that I'm a new user, give me an overview of the different ways we can
perform strong lens modeling of this system.
```

<sub><b>Example Natural Language Prompt for Claude Code, Codex or other AI coding agent</b></sub>

```text
Teacher mode.

I'm new to PyAutoLens and want to learn the basic workflow end-to-end. Walk me through
simulating Euclid-like imaging of a simple strong lens, plotting it, and fitting it.
```

<sub><b>Example Natural Language Prompt for Claude Code, Codex or other AI coding agent</b></sub>

```text
I have HST imaging of a galaxy-scale lens. Help me plan the model: lens light, mass, and
source. Ask me what you need to know about the data first — don't run anything yet.
```

The last one exercises two things that matter. First, you can hold a planning discussion inside
the agent without it executing anything — say so, as the prompt does. Second, on **real data**
the assistant is required to make you look at the image before it composes a fit, and to settle
two things with you: whether there are extra galaxies or artefacts in the frame, and how big the
mask should be. It will plot the data itself and show you the file; that is the rule working,
not the assistant being slow.

More ambitious examples — dark-matter subhalo detection, joint imaging + interferometer +
weak-lensing fits — are in the [README](../../README.md).
