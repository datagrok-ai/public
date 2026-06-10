"""Enrichment CLI.

Three commands:
  py enrich.py prepare --enricher feature_clustering --packages Chem,Chembl,PowerPack
  py enrich.py apply   --enricher feature_clustering [--packages ...]
  py enrich.py rebuild-db                # rebuild Kuzu after apply

Step 2 (the LLM call itself) is intentionally external — it can be Claude
Code agents, the Anthropic SDK, a batch job, or hand-curation.
"""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

ENRICHERS = {
    "feature_clustering":     "enrichers.feature_clustering",
    "doc_linking":            "enrichers.doc_linking",
    "feature_implementation": "enrichers.feature_implementation",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("prepare")
    sp.add_argument("--enricher", required=True, choices=list(ENRICHERS))
    sp.add_argument("--packages", default="")

    sa = sub.add_parser("apply")
    sa.add_argument("--enricher", required=True, choices=list(ENRICHERS))
    sa.add_argument("--packages", default="")

    sr = sub.add_parser("rebuild-db")

    args = ap.parse_args()
    if args.cmd in ("prepare", "apply"):
        mod = importlib.import_module(ENRICHERS[args.enricher])
        pkgs = [p for p in (args.packages or "").split(",") if p] or None
        if args.cmd == "prepare":
            mod.prepare(pkgs)
        else:
            mod.apply(pkgs)
        return 0

    if args.cmd == "rebuild-db":
        # Just delegate to build.py's DB-load step by re-running it
        import build
        sys.argv = ["build.py"]   # default: full build, with DB
        return build.main()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
