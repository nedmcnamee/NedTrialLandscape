"""Post-build checks. Runs in CI after build.py and fails the job loudly.

Silent drift is the failure mode that matters here: upstream adds a new cancer
type or drug class, it quietly lands in "Unmapped", and the figure is wrong
without anyone noticing. These checks make that a build failure instead.
"""
from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "docs" / "data.json"

# Tunable guardrails
MIN_TRIALS = 300           # a collapse below this means a source broke
MAX_UNMAPPED_PCT = 10.0    # tissue grouping coverage floor
MIN_REGISTRIES = 1


def main() -> int:
    if not DATA.exists():
        print("FAIL: docs/data.json missing -- build.py did not run")
        return 1

    d = json.loads(DATA.read_text())
    trials = d["trials"]
    problems, warnings = [], []

    if len(trials) < MIN_TRIALS:
        problems.append(f"only {len(trials)} trials (expected >= {MIN_TRIALS}) "
                        f"-- a source probably failed silently")

    registries = {t["registry"] for t in trials}
    if len(registries) < MIN_REGISTRIES:
        problems.append(f"no registries returned data")
    for r in sorted(registries):
        n = sum(1 for t in trials if t["registry"] == r)
        print(f"  {r}: {n} trials")

    unmapped = sorted({c for t in trials if "Unmapped" in t["tissue_groups"]
                       for c in t["cancer_types"]})
    n_unmapped = sum(1 for t in trials if t["tissue_groups"] == ["Unmapped"])
    pct = 100 * n_unmapped / max(len(trials), 1)
    print(f"  tissue grouping: {pct:.1f}% of trials fully unmapped")
    if pct > MAX_UNMAPPED_PCT:
        problems.append(f"{pct:.1f}% of trials have no tissue group "
                        f"(limit {MAX_UNMAPPED_PCT}%)")
    if unmapped:
        warnings.append("cancer types with no tissue group -- add them to "
                        "config.yaml tissue_groups:\n    " +
                        "\n    ".join(unmapped[:40]))

    no_target = sum(1 for t in trials if not t["classes"])
    if no_target:
        warnings.append(f"{no_target} trials have no drug class assigned")

    for w in warnings:
        print(f"\nWARNING: {w}")
    for p in problems:
        print(f"\nFAIL: {p}")

    if problems:
        return 1
    print("\nvalidate: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
