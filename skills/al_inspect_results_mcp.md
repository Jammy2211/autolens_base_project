---
name: al_inspect_results_mcp
description: Run and configure the read-only results-inspector MCP server, which lets an MCP client (Claude Desktop, Claude Code) inspect PyAutoLens fit results without running Python — list fits ranked by evidence, read model and posterior summaries, view result images inline in chat, and extract/combine subplot panels and FITS HDUs across many fits. Use when the user wants to "browse/inspect my results from chat", asks about the MCP server, or wants Claude Desktop wired to their output folder. Not for loading results in Python (that is `al_load_results` / `al_aggregator_bulk_analysis`) and not for running fits.
user-invocable: true
---

# The results-inspector MCP server

`autoassistant/mcp/` is a read-only MCP (Model Context Protocol) stdio server over
PyAutoFit/PyAutoLens output directories. It exists for harnesses that cannot execute
code — a Claude Desktop chat gets tools to list fits, read summaries and display result
images inline, against the same `output/` folder a script-based session works with.

It is deliberately **not** a fitting interface: composing models and running searches
stay python-first through the other skills. Exposing `search.fit` through a JSON tool
schema flattens the compositional API and is out of scope by design.

## Orient

- Server: `autoassistant/mcp/server.py`, run as `python -m autoassistant.mcp` from the
  repo root (stdio; nothing listens on a port).
- The **core tools are mirrored verbatim** from
  `autofit_assistant:autoassistant/mcp/` (`tools.py` and `__main__.py` byte-identical
  — `diff` them when syncing; graduation to an `autofit[mcp]` extra is the recorded
  de-dup path). `lens_tools.py` and the lens registrations in `server.py` are this
  repo's own layer.
- The tool core ships as an optional PyAutoFit extra: `pip install autofit[mcp]`
  (the lens image/FITS tools then layer on top, needing PyAutoLens installed)
  (assistant-environment dependency only — never a library requirement).
- All tools are read-only against `output/`; `extract_lens_fits` writes only to an
  explicit destination outside it.

## Tools

Core (mirrored): `list_searches`, `get_model`, `get_result_summary`,
`get_samples_summary`, `get_search_info`, `list_images`, `fetch_image` — see
`autofit_assistant:skills/af_inspect_results_mcp.md` for the full table.

Lens layer:

| Tool | Returns |
|------|---------|
| `list_extractable_images()` | The `al.agg` enum groups and names the two tools below accept, as "group.name" specs (e.g. `subplot_fit.data`, `fits_tracer.convergence`). |
| `combine_lens_images(directory, subplots, subplot_width=0)` | Named panels extracted from every fit under `directory`, combined into one image rendered inline in chat (`af.AggregateImages`). |
| `extract_lens_fits(directory, hdus, destination_path, overwrite=False)` | A single `.fits` of the named HDUs from every fit, written to `destination_path` (refused inside the output directory); returns the path (`af.AggregateFITS`). |

Wrappers over `PyAutoFit:autofit/aggregator/` (`Aggregator.from_directory`,
`SearchOutput`, `AggregateImages`, `AggregateFITS`) and the
`PyAutoLens:autolens/aggregator/subplot.py` enum groups exposed as `al.agg.*`.

## Configure a client

**Claude Desktop** (`claude_desktop_config.json`, Settings → Developer → Edit Config):

```json
{
  "mcpServers": {
    "pyauto-results-inspector": {
      "command": "/absolute/path/to/venv/bin/python",
      "args": ["-m", "autoassistant.mcp"],
      "env": { "PYTHONPATH": "/absolute/path/to/autolens_assistant" }
    }
  }
}
```

Use the interpreter that has PyAutoLens + `mcp` installed. Ask about fits by absolute
path ("list the fits under /home/me/project/output ranked by evidence, then show me
the data and residuals of the best one side by side").

MCP clients spawn the server with a **minimal environment** — nothing from your shell
propagates. Anything the stack needs (`PYTHONPATH` for editable/source checkouts,
`NUMBA_CACHE_DIR`/`MPLCONFIGDIR` in restricted setups) must be declared in the
config's `env` block.

**Windows (Claude Desktop → WSL):** when the interpreter lives in WSL, launch it
through `wsl.exe`. Nothing beyond `PYTHONPATH` is required — the server pins its own
config directory and forces JAX onto CPU before importing autofit, so it does not
depend on the launch directory or extra environment:

```json
{
  "mcpServers": {
    "pyauto-results-inspector": {
      "command": "wsl.exe",
      "args": ["-e", "bash", "-c",
        "PYTHONPATH=/home/you/autolens_assistant /home/you/venv/bin/python -m autoassistant.mcp"]
    }
  }
}
```

For the Microsoft Store build of Claude Desktop, the config and the failure log
(`logs/mcp-server-pyauto-results-inspector.log`) live under
`%LOCALAPPDATA%\Packages\Claude_*\LocalCache\Roaming\Claude\`, not `%APPDATA%\Claude\`.

**Claude Code**: the repo-root `.mcp.json` registers the same server automatically for
sessions opened in this repo.

## Deployment tiers

1. **Local stdio (built, above)** — Claude Desktop / Claude Code on the machine that
   holds `output/`.
2. **Remote (documented only)** — claude.ai web/mobile custom connectors and ChatGPT
   developer mode speak MCP but only to servers reachable over the public internet
   (no stdio): expose the server via `mcp.run(transport="streamable-http")` behind an
   ngrok/cloudflared tunnel. Not built or hardened here; do not tunnel a machine you
   care about without thinking about auth.
3. **Hosted (future)** — a collaboration-scale deployment next to shared outputs
   (e.g. Euclid sample-wide triage). Same tools; hosting, auth and scale are its own
   task.

## Design rules (maintainers)

- **Glue, not code.** Every tool is argument parsing + one existing public
  PyAutoFit/PyAutoLens call + serialization. If a tool needs more, add the method to
  the library first.
- **Read-only.** No fit-running, no compute, no writes into `output/`.
- **stdout is the protocol.** Autofit calls run under `tools._stdout_to_stderr()` and
  `server._route_logging_to_stderr()` rebinds stdout log handlers — keep both when
  adding tools; a single stray print corrupts the JSON-RPC channel.
- **Anti-drift.** `autoassistant/mcp/*.py` is scanned by the symbol audit via
  `autoassistant/audit_skill_apis.py` (`--scope scripts`); tests
  (`autoassistant/tests/test_mcp_tools.py`) build their fixture by running a real
  tiny fit so format drift fails loudly.
