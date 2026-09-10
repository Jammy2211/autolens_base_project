<!-- ARCHIVED 2026-09-10 — UNSUPPORTED. Conversational chat routes (ChatGPT / Claude chat with a GitHub
     connector, the custom GPT, pasted bundles, Claude Project knowledge packs) were retired when the assistants
     moved to an agentic-only support policy. This page is kept for reference only: it is not maintained, the
     bundles it refers to (`chat_pack/`, `llms-chat.txt`, `AGENTS_CHAT.md`) no longer exist in the repository, and
     the setup it describes is not a supported way to use the assistant. Use a coding agent instead — see the
     top-level README "Getting Started". -->

> **Archived and unsupported (2026-09-10).** This describes a retired conversational-chat route. See the
> [README](../../../README.md) for the supported agentic setup.

# Claude chat — Free plan (Project + knowledge pack)

Run the assistant inside claude.ai on the **Free plan**: a curated pack of this repository is
uploaded straight into a Project's knowledge. On Free, project knowledge is put in **full
context** (~200K tokens) rather than retrieved in fragments, so the assistant sees all of it
from the first message.

**How well it works:** well — the assistant plans models, writes current-API scripts and
teaches concepts faithfully — **but it goes through the free token allowance quickly**. Free
usage limits are unpublished and vary with load (in the region of a few dozen messages per 5
hours), and every message carries the pack. Fine for planning sessions and short questions;
long modelling conversations will hit the ceiling — that is where a
[paid plan](claude_chat_paid.md) or a free coding agent such as [OpenCode](../../setup/opencode_cli.md) pays for itself.

## Setup (~5 minutes, once)

1. Download the files in `chat_pack/` (removed — see the archive notice) (11 files, ~61k tokens total).
   Easiest way without git: download the repository ZIP from the green **Code** button on
   GitHub and take that folder.
2. In Claude, create a **Project** (Free allows up to 5), then upload every file from
   `chat_pack/` into its **knowledge**.
3. Paste this into the project's **custom instructions**:

```text
You are the PyAutoLens Assistant. Follow 00_instructions.md in your project knowledge
exactly — especially: never write PyAutoLens from memory, check symbols against
01_api_surface.md, and use only the functional plotting API (aplt.subplot_*), never the
removed object-oriented plotters. Lead by engaging with my science goal; raise the
handoff to a local coding agent only when actually running code is the blocker.
```

Every chat in that project is now configured. You can add your own papers or data notes to
the same project knowledge.

> **Notes.** Claude Projects can't be shared on Free or Pro (sharing is a Team/Enterprise
> feature), so each person does this once for themselves. The Free tier's GitHub connector is
> missing features which hurt performance, so use the Project + knowledge pack here (see
> [Troubleshooting](../../setup/troubleshooting.md)); on a paid plan the connector works excellently —
> see [Claude chat — paid plans](claude_chat_paid.md).

## What it can and can't do

In chat the assistant does the thinking work — planning models, writing scripts for you to
run, explaining concepts, reviewing errors and figures you paste in. It **cannot** run fits,
read your `.fits` files or inspect your results folder; when that is the blocker it will say
so and point you to a coding agent.

## Your first prompt

You're set up — copy and paste this to start (the COSMOS-Web Ring data ships with the
repository, so it works immediately):

```text
Find the data on the Cosmos-Web ring, give me a short script to plot it in PyAutoLens
and then given that I'm a new user give me an overview of the different ways we can
perform strong lens modeling of this system.
```

More examples: [First prompts to try](../../setup/first_prompts.md). If anything misbehaves:
[Troubleshooting](../../setup/troubleshooting.md).
