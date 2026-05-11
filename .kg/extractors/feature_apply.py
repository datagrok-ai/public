"""Re-apply LLM-derived Feature outputs from `data/enrichment/*/outputs/`.

Bridges the enricher pipeline into the deterministic extractor pipeline so
that `build.py --clean` doesn't lose Features. The clean step preserves
`data/enrichment/` (the precious LLM outputs); this extractor reads them
and emits `Feature` / `HasFeature` / `PartOfFeature` JSONL slices.

Two enrichers are applied:
  - `feature_clustering` -> Feature nodes + HasFeature/PartOfFeature edges
  - `feature_implementation` -> additional IsImplementedIn/IsTestedIn edges

Both `apply()` functions are no-ops when their `outputs/` directory has
nothing new. Safe to run on every build.
"""

from __future__ import annotations

from pathlib import Path

from . import _common as c

EXTRACTOR_NAME = "feature_apply"


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    # Side-step the Bundle: the enrichers write directly to JSONL because they
    # produce thousands of rows and want their own `_merge_into` path. The
    # Bundle's role is to dedupe entities/edges within one extractor's
    # output; the enrichers handle dedupe themselves.
    try:
        from enrichers.feature_clustering import apply as fc_apply
    except ImportError as e:
        print(f"[{EXTRACTOR_NAME}] feature_clustering enricher not importable ({e})")
        return
    try:
        from enrichers.feature_implementation import apply as fi_apply
    except ImportError:
        fi_apply = None

    pkgs = list(package_filter) if package_filter else None
    fc_apply(packages=pkgs)
    if fi_apply is not None:
        try:
            fi_apply(packages=pkgs)
        except TypeError:
            # Older signature — no packages kwarg
            fi_apply()
