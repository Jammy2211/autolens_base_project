"""Tests for the benchmark harness (autoassistant/benchmark.py)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from autoassistant import benchmark

REPO_ROOT = Path(__file__).resolve().parents[2]

CARD = """\
---
id: test-bench
version: 2
mode: assistant
difficulty: easy
datasets: []
workspace_packages:
  - imaging
added: 2026-07-10
---

# Benchmark: test

## Prompt

```
Assistant mode. Do the thing.
```

## Success rubric (30 points)

### Machine-checkable (20)

| # | Check | Pts |
|---|-------|-----|
| M1 | A thing exists | 10 |
| M2 | Another thing exists | 10 |

### Judged (10)

| # | Criterion | Pts |
|---|-----------|-----|
| J1 | The thing was done well | 10 |
"""


@pytest.fixture
def root(tmp_path):
    prompts = tmp_path / "benchmarks" / "prompts"
    prompts.mkdir(parents=True)
    (prompts / "test_bench.md").write_text(CARD)
    return tmp_path


def fill_score(run_dir, awards):
    path = run_dir / "score.md"
    lines = []
    for line in path.read_text().splitlines():
        m = benchmark.SCORE_ROW.match(line)
        if m and m.group(1) in awards:
            code, criterion, max_points = m.group(1), m.group(2), m.group(3)
            lines.append(
                f"| {code} | {criterion} | {max_points} | {awards[code]} | evidence |"
            )
        else:
            lines.append(line)
    path.write_text("\n".join(lines) + "\n")


def test_load_card_parses_frontmatter_and_rubric(root):
    cards = benchmark.load_cards(root)
    card = cards["test-bench"]
    assert card.version == 2
    assert [r.code for r in card.rubric] == ["M1", "M2", "J1"]
    assert card.rubric[0].max_points == 10
    assert card.rubric[0].machine and not card.rubric[2].machine


def test_new_run_scaffolds_run_directory(root):
    run_dir = benchmark.new_run(root, "test-bench", "Claude Opus 4.8", "claude-code")
    assert run_dir.name.endswith("_claude-opus-4.8_claude-code")
    assert (run_dir / "transcript.md").exists()
    assert (run_dir / "artifacts").is_dir()
    meta = yaml.safe_load((run_dir / "meta.yaml").read_text())
    assert meta["benchmark"] == "test-bench"
    assert meta["prompt_version"] == 2
    assert meta["status"] == "pending"
    assert set(meta["stack"]) == set(benchmark.STACK_PACKAGES)
    score_text = (run_dir / "score.md").read_text()
    assert "| M1 |" in score_text and "| J1 |" in score_text


def test_new_run_same_day_repeat_gets_suffix(root):
    first = benchmark.new_run(root, "test-bench", "m", "h")
    second = benchmark.new_run(root, "test-bench", "m", "h")
    assert second.name == f"{first.name}_2"


def test_new_run_unknown_benchmark_fails(root):
    with pytest.raises(SystemExit):
        benchmark.new_run(root, "nope", "m", "h")


def test_score_run_totals_and_updates_meta(root):
    run_dir = benchmark.new_run(root, "test-bench", "m", "h")
    fill_score(run_dir, {"M1": 10, "M2": 0, "J1": 7.5})
    score = benchmark.score_run(run_dir)
    assert score.machine == 10 and score.judged == 7.5 and score.total == 17.5
    meta = yaml.safe_load((run_dir / "meta.yaml").read_text())
    assert meta["score"]["total"] == 17.5
    assert meta["status"] == "complete"


def test_score_run_rejects_unfilled_rows(root):
    run_dir = benchmark.new_run(root, "test-bench", "m", "h")
    fill_score(run_dir, {"M1": 10})
    with pytest.raises(SystemExit, match="unfilled"):
        benchmark.score_run(run_dir)


def test_score_run_rejects_over_max_award(root):
    run_dir = benchmark.new_run(root, "test-bench", "m", "h")
    fill_score(run_dir, {"M1": 11, "M2": 0, "J1": 0})
    with pytest.raises(SystemExit, match="outside"):
        benchmark.score_run(run_dir)


def test_report_builds_leaderboard_and_pending(root):
    scored = benchmark.new_run(root, "test-bench", "model-a", "h")
    fill_score(scored, {"M1": 10, "M2": 10, "J1": 5})
    benchmark.score_run(scored)
    pending = benchmark.new_run(root, "test-bench", "model-b", "h")

    text = benchmark.report(root)
    assert "## test-bench" in text
    assert "| model-a | h | 1 | 25 | 25 |" in text
    assert f"`{pending.relative_to(root / 'benchmarks')}`" in text

    path = benchmark.write_report(root)
    assert path == root / "benchmarks" / "RESULTS.md"
    assert path.read_text() == text


def test_repo_prompt_cards_parse():
    """Every committed prompt card must load: unique ids, frontmatter, rubric."""
    cards = benchmark.load_cards(REPO_ROOT)
    assert set(cards) == {
        "assistant-easy-cosmos-web-ring",
        "assistant-medium-slacs0946-subhalo",
        "assistant-hard-group-multi",
        "teacher-basic-workflow",
        "harness-smoke",
    }
    for card in cards.values():
        machine = sum(r.max_points for r in card.rubric if r.machine)
        judged = sum(r.max_points for r in card.rubric if not r.machine)
        assert machine + judged == 100, f"{card.id}: rubric totals {machine + judged}"


def test_repo_card_datasets_exist():
    """Bundled datasets a card declares must exist — a missing one is a stale card."""
    for card in benchmark.load_cards(REPO_ROOT).values():
        for dataset in card.meta.get("datasets", []):
            assert (REPO_ROOT / dataset).is_dir(), f"{card.id}: missing {dataset}"


def test_repo_readme_prompts_match_cards():
    """Cards that mirror README example prompts must stay textually identical."""
    readme = (REPO_ROOT / "README.md").read_text()
    for name in ("easy_cosmos_web_ring", "medium_slacs0946_subhalo", "teacher_workflow"):
        card_text = (REPO_ROOT / "benchmarks" / "prompts" / f"{name}.md").read_text()
        prompt = card_text.split("```\n", 1)[1].split("```", 1)[0]
        assert prompt.strip() in readme, f"{name}: prompt text diverges from README.md"
