<!-- ARCHIVED 2026-09-10 — UNSUPPORTED. Conversational chat routes (ChatGPT / Claude chat with a GitHub
     connector, the custom GPT, pasted bundles, Claude Project knowledge packs) were retired when the assistants
     moved to an agentic-only support policy. This page is kept for reference only: it is not maintained, the
     bundles it refers to (`chat_pack/`, `llms-chat.txt`, `AGENTS_CHAT.md`) no longer exist in the repository, and
     the setup it describes is not a supported way to use the assistant. Use a coding agent instead — see the
     top-level README "Getting Started". -->

> **Archived and unsupported (2026-09-10).** This describes a retired conversational-chat route. See the
> [README](../../../README.md) for the supported agentic setup.

# Claude chat — paid plans (GitHub connector)

Run the assistant inside claude.ai on a **paid plan** (Pro/Max/Team) by connecting this
repository through Claude's GitHub connector. **It works excellently:** the assistant reads
the repository directly, always sees current content, and grounds every PyAutoLens symbol in
the live API surface. The higher usage limits also remove the token squeeze that constrains
the [Free-plan route](claude_chat_free.md), so long modelling conversations, deep concept
dives and multi-script planning sessions all run without hitting the ceiling.

On Claude's Free plan the connector is missing features which hurt performance, which is why
that route uses a Project + knowledge pack instead.

## Setup (~2 minutes)

1. In Claude, open **Settings → Connectors** and connect **GitHub**. Authorise access to
   public repositories (this repo is public — you don't need to grant anything of your own),
   then attach `https://github.com/PyAutoLabs/autolens_assistant`.
2. Start a chat with this bootstrap prompt:

```text
Use the autolens_assistant repository: https://github.com/PyAutoLabs/autolens_assistant

Start by reading its front door:
https://raw.githubusercontent.com/PyAutoLabs/autolens_assistant/main/llms.txt

Follow its read order (AGENTS_CHAT.md → the relevant skill → wiki) and its API rules.
First tell me whether you can actually read llms.txt — if you can't, say so plainly
and don't answer from memory.
```

**Naming `llms.txt` explicitly matters** — connectors do not reliably find it on their own,
and answers are markedly better when it is pointed there first. The final sentence is a
deliberate honesty check: if the assistant can't read the file, you want to know immediately,
not after it writes you a script from a 2023 API.

**Prefer uploading to connecting?** Project knowledge remains a fine alternative: follow
[Claude chat — Free plan](claude_chat_free.md) steps 1–3. On a paid plan the higher limits
leave room for your own papers, data notes and analysis logs in the same project knowledge
too.

## What it can and can't do

In chat the assistant does the thinking work — planning models, writing current-API scripts
for you to run, explaining concepts, reviewing errors and figures. It **cannot** run fits or
read your `.fits` files; for that, pair it with a coding agent such as
[Claude Code](../../setup/claude_code.md), which shares your Claude subscription.

## Your first prompt

You're set up — copy and paste this to start (the COSMOS-Web Ring data ships with the
repository, so it works immediately):

```text
Use the autolens_assistant repository: https://github.com/PyAutoLabs/autolens_assistant

Start by reading its front door:
https://raw.githubusercontent.com/PyAutoLabs/autolens_assistant/main/llms.txt

Follow its read order (AGENTS_CHAT.md → the relevant skill → wiki) and its API rules.
First tell me whether you can actually read llms.txt — if you can't, say so plainly
and don't answer from memory.

Find the data on the Cosmos-Web ring, give me a short script to plot it in PyAutoLens
and then given that I'm a new user give me an overview of the different ways we can
perform strong lens modeling of this system.
```

More examples: [First prompts to try](../../setup/first_prompts.md). If anything misbehaves:
[Troubleshooting](../../setup/troubleshooting.md).
