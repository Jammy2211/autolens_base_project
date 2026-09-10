# docs/archive — retired material (unsupported)

Everything in this folder is kept **for reference only**. It is not maintained, it is not linked from active
onboarding, and none of it describes a supported way to use the assistant.

## Conversational chat routes (retired 2026-09-10)

The assistant is supported **only inside an AI coding agent** (Claude Code or Codex; OpenCode as an experimental
alternative). Between July and September 2026 it also offered browser-chat routes; they were retired because a chat
that cannot execute code or read data cannot enforce the scientific safeguards the assistant depends on
(current-API verification, inspecting the data before fitting, checking a fit's outputs), and because keeping a
parallel chat-mode instruction set and generated knowledge bundles current was ongoing work that agentic use never
needed.

| Archived page | What it was |
|---|---|
| [`CHOOSING_YOUR_AI_TOOL.md`](CHOOSING_YOUR_AI_TOOL.md) | The old chooser between chat and coding agents, free and paid (last updated 2026-08-06) |
| [`README_choosing_your_ai_tool_section.md`](README_choosing_your_ai_tool_section.md), [`README_ai_chat_assistant_section.md`](README_ai_chat_assistant_section.md) | README sections removed on 2026-09-03 |
| [`chat/chatgpt_paid_connector.md`](chat/chatgpt_paid_connector.md) | ChatGPT paid plans via the GitHub connector |
| [`chat/claude_chat_paid.md`](chat/claude_chat_paid.md) | Claude chat paid plans via the GitHub connector |
| [`chat/claude_chat_free.md`](chat/claude_chat_free.md) | Claude Free via a Project + uploaded knowledge pack |
| [`chat/chatgpt_custom_gpt.md`](chat/chatgpt_custom_gpt.md) | The experimental "PyAutoLens AI Assistant" custom GPT (built 2026-08-06) |
| [`chat/paste_bundle.md`](chat/paste_bundle.md) | Pasting a generated bundle into any chat |

**Removed rather than archived** (recoverable from git history before this retirement): `AGENTS_CHAT.md` (the
chat-mode subset of `AGENTS.md`), the generated `chat_pack/` knowledge pack and `llms-chat.txt` paste bundle, their
generator `autoassistant/chat_bundle.py` with its tests and `make chat-bundle{,-check}` targets, and the
`FREE_TIER_SETUP.md` redirect stub. The custom GPT link in the archived page may still resolve, but the GPT is not
maintained.

The dated observations these pages record (connector behaviour on 2026-08-06, free-tier limits) are left as written;
they were true when measured and are not being re-verified.
