<!-- ARCHIVED 2026-09-10 — UNSUPPORTED. Conversational chat routes (ChatGPT / Claude chat with a GitHub
     connector, the custom GPT, pasted bundles, Claude Project knowledge packs) were retired when the assistants
     moved to an agentic-only support policy. This page is kept for reference only: it is not maintained, the
     bundles it refers to (`chat_pack/`, `llms-chat.txt`, `AGENTS_CHAT.md`) no longer exist in the repository, and
     the setup it describes is not a supported way to use the assistant. Use a coding agent instead — see the
     top-level README "Getting Started". -->

> **Archived and unsupported (2026-09-10).** This describes a retired conversational-chat route. See the
> [README](../../../README.md) for the supported agentic setup.

# Paste the bundle (works in any chat, free or paid)

The universal fallback: a single self-contained paste that configures **any** conversational
AI — ChatGPT Free, Gemini, Copilot, anything with a text box. No connectors, no uploads, no
accounts beyond the chat itself. It is deliberately compact (~6k tokens) so it fits free-tier
context windows and still leaves room for an actual conversation.

## Setup (~30 seconds, per chat)

1. Open **`llms-chat.txt` (removed — see the archive notice)** and copy the whole file.
2. Paste it as your first message in a new chat, followed by your question.

That file is self-contained: chat-mode instructions, the complete generated list of public
PyAuto\* symbols, and a routing table of raw URLs. If your chat has browsing, it can fetch
any skill or wiki page from those URLs; if not, it still has the rules and the API surface,
and it is instructed to tell you when an answer would need a page it can't reach.

Two things to know:

- **Re-paste it in each new chat.** Nothing persists between conversations on a free plan
  unless you put it in custom instructions or a project.
- **ChatGPT tip:** paste the bundle and add your question in the *same* message. A first
  message that is only context sometimes gets a "what would you like to do?" reply that
  wastes one of your limited flagship-model turns.

## Your first prompt

You're set up — copy and paste this to start (the COSMOS-Web Ring data ships with the
repository, so it works immediately):

```text
[paste the full contents of llms-chat.txt here]

Find the data on the Cosmos-Web ring, give me a short script to plot it in PyAutoLens
and then given that I'm a new user give me an overview of the different ways we can
perform strong lens modeling of this system.
```

More examples: [First prompts to try](../../setup/first_prompts.md). If anything misbehaves:
[Troubleshooting](../../setup/troubleshooting.md).
