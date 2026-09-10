<!-- ARCHIVED 2026-09-10 — UNSUPPORTED. Conversational chat routes (ChatGPT / Claude chat with a GitHub
     connector, the custom GPT, pasted bundles, Claude Project knowledge packs) were retired when the assistants
     moved to an agentic-only support policy. This page is kept for reference only: it is not maintained, the
     bundles it refers to (`chat_pack/`, `llms-chat.txt`, `AGENTS_CHAT.md`) no longer exist in the repository, and
     the setup it describes is not a supported way to use the assistant. Use a coding agent instead — see the
     top-level README "Getting Started". -->

> **Archived and unsupported (2026-09-10).** This describes a retired conversational-chat route. See the
> [README](../../../README.md) for the supported agentic setup.

# ChatGPT — the custom GPT (any plan, experimental)

**The GPT is live:**
**[PyAutoLens AI Assistant](https://chatgpt.com/g/g-6a74c33c58c48191b8cd353e7b46f18b-pyautolens-ai-assistant)**
(built 2026-08-06 against the 2026.8.4.1 `chat_pack/`). Open the link and start chatting —
no setup, on any ChatGPT plan including Free.

**How well it works: it works, but is not yet able to do all tasks.** In testing (paid and
free accounts, 2026-08-06) the GPT's uploaded knowledge was not reliably retrievable in
published chats (a known custom-GPT platform issue), and on free accounts browsing could not
fetch raw GitHub pages either. The GPT is honest about this — it will *decline* to write code
it can't verify rather than fabricate a stale API, which is the designed behaviour. It
remains good for planning a lens model, explaining concepts, and reviewing your scripts and
errors. **If it says it can't retrieve its documentation, or you need a verified full
script, [paste the bundle](paste_bundle.md) instead** — that puts everything in-context and
is immune to both problems.

Free ChatGPT users **can use** custom GPTs; they just can't create them. So one person with a
Plus/Pro account builds the GPT once and shares the link — the only first-class way to hand a
configured assistant to free ChatGPT users.

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

## Build recipe (maintainers)

For rebuilding the GPT after a `chat_pack/` regeneration — regenerate the inputs first (see
"Maintaining the chat bundles" in [`modes/maintainer.md`](../../../modes/maintainer.md)):

1. **Create a GPT** → *Configure*.
2. **Instructions**: paste the contents of `chat_pack/00_instructions.md`.
3. **Knowledge**: upload all files from `chat_pack/`. They are pre-split into ≤20 topic files
   with distinctive headings, because GPT knowledge is retrieved by **RAG chunking** rather
   than read whole — one topic per file is what makes retrieval land on the right chunk.
4. **Capabilities**: enable *Web Browsing* (so it can fetch skills and wiki pages beyond the
   pack). Code Interpreter is optional and does **not** give it PyAutoLens — the library is
   not installed in that sandbox, so it cannot run a fit there.
5. *(Optional)* Add an **Action** against `raw.githubusercontent.com` for always-current
   content. Publishing a GPT with an Action requires a privacy-policy URL.
6. Publish **Anyone with the link**, and record the link here.
