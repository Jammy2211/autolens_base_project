# Claude Code (coding agent — recommended)

The **primary, most thoroughly exercised** harness for the assistant. Claude Code is Anthropic's
coding agent (terminal, desktop app, and IDE extensions): it reads the repository's instructions
through `CLAUDE.md`, runs the code-gate hook that blocks stale PyAutoLens API written from memory,
and can install PyAutoLens, run fits end-to-end, inspect results and drive HPC workflows.

**Access.** Claude Code is a paid product. Anthropic's
[cost page](https://code.claude.com/docs/en/costs) describes the routes: a Claude subscription
(Pro / Max, or a Team / Enterprise seat through your institution), usage-based billing through
the Claude Console (API), or a cloud provider your institution already uses (Amazon Bedrock,
Google Cloud, Microsoft Foundry). Which is cheapest depends on how much you use it; sustained
lens modelling is heavy use, so budget for it rather than assuming a personal subscription is
the only or best option.

## Setup

1. Install Claude Code — follow the official instructions at
   [code.claude.com/docs](https://code.claude.com/docs).
2. Clone this repository and start the agent inside it:

```bash
git clone https://github.com/PyAutoLabs/autolens_assistant.git
cd autolens_assistant
claude
```

The assistant configures itself on your first prompt, and will install PyAutoLens for you if
it isn't already installed. The desktop app and the IDE extensions work the same way once they
are opened on this folder.

## Your first prompt

You're set up — copy and paste this to start (the COSMOS-Web Ring data ships with the
repository, so it works immediately):

<sub><b>Example Natural Language Prompt for Claude Code, Codex or other AI coding agent</b></sub>

```text
Find the data on the Cosmos-Web ring, give me a short script to plot it in PyAutoLens
and then given that I'm a new user give me an overview of the different ways we can
perform strong lens modeling of this system.
```

You do not have to let it run anything: start with `Teacher mode.` to learn, or ask it to plan
an analysis and discuss the choices before any fit is submitted.

More examples: [First prompts to try](first_prompts.md). If anything misbehaves:
[Troubleshooting](troubleshooting.md).
