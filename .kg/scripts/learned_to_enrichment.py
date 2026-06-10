"""Drain `.kg/.learned/*.md` into the canonical enrichment store.

Each `.learned/*.md` file holds tab-separated facts the agents found
during a /dg-task session that the graph didn't already know:

    feature:Bio:atomic-level-structure  IS_IMPLEMENTED_IN  packages/Bio/src/utils/foo.ts  reason: ...
    func:Bio:To Atomic Level            IS_TESTED_IN       packages/Bio/src/tests/x.ts    reason: ...

Currently supported predicates: `IS_IMPLEMENTED_IN`, `IS_TESTED_IN`.
We append directly to the matching JSONL edge slice with
`derived_by: manual` and `extracted_by: learned` so the next Kuzu
rebuild absorbs them. After successful drain the source `.md` files
are moved into `.learned/applied/` (kept for audit, never deleted).

Idempotent: if the same edge already exists in the slice, it's
skipped (Bundle.merge_jsonl dedupes by (from_id, to_id, predicate)).

Run by the Stop hook in `.claude/hooks/refresh-kg.sh`. Safe
to invoke manually: `py .kg/scripts/learned_to_enrichment.py`.
"""

from __future__ import annotations

import json
import re
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEARNED_DIR = ROOT / ".learned"
APPLIED_DIR = LEARNED_DIR / "applied"
EDGES_DIR = ROOT / "data" / "edges"

SUPPORTED = {"IS_IMPLEMENTED_IN", "IS_TESTED_IN"}

LINE_RE = re.compile(
    r"""^
        (?P<subject>[^\s]+)\s+
        (?P<predicate>[A-Z_]+)\s+
        (?P<object>[^\s]+)
        (?:\s+reason:\s*(?P<reason>.+))?
        \s*$
    """, re.VERBOSE,
)


def _ensure_file_id(obj: str) -> str:
    """Allow bare paths in the .md (gets `file:` prefix added)."""
    if obj.startswith("file:") or obj.startswith("func:") or obj.startswith("feature:"):
        return obj
    if obj.endswith((".ts", ".tsx", ".js", ".py", ".sql", ".md", ".json", ".yaml", ".css")):
        return f"file:{obj}"
    return obj


def _parse(path: Path) -> list[dict]:
    out: list[dict] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        # Allow tabs OR runs of spaces as separator
        line = re.sub(r"\s+", " ", line)
        m = LINE_RE.match(line)
        if not m:
            print(f"  ! {path.name}: skip unparsable line: {raw[:120]}", file=sys.stderr)
            continue
        pred = m.group("predicate")
        if pred not in SUPPORTED:
            print(f"  ! {path.name}: unsupported predicate {pred} (ignored)", file=sys.stderr)
            continue
        out.append({
            "from_id": m.group("subject"),
            "to_id": _ensure_file_id(m.group("object")),
            "predicate": pred,
            "reason": (m.group("reason") or "").strip() or None,
        })
    return out


def _existing_keys(jsonl: Path) -> set[tuple[str, str, str]]:
    if not jsonl.is_file():
        return set()
    keys: set[tuple[str, str, str]] = set()
    for line in jsonl.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            continue
        keys.add((obj.get("from_id", ""), obj.get("to_id", ""), obj.get("predicate", "")))
    return keys


def _append(jsonl: Path, rows: list[dict]) -> int:
    if not rows:
        return 0
    EDGES_DIR.mkdir(parents=True, exist_ok=True)
    existing = _existing_keys(jsonl)
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    written = 0
    with jsonl.open("a", encoding="utf-8", newline="\n") as fh:
        for row in rows:
            key = (row["from_id"], row["to_id"], row["predicate"])
            if key in existing:
                continue
            edge = {
                "from_id": row["from_id"],
                "to_id": row["to_id"],
                "predicate": row["predicate"],
                "evidence": [],
                "derived_by": "manual",
                "confidence": 0.9,
                "extras": {"reason": row["reason"], "source": "learned"},
                "extracted_by": "learned",
                "extracted_at": now,
                "member_count": 1,
            }
            fh.write(json.dumps(edge, separators=(",", ":"), ensure_ascii=False) + "\n")
            existing.add(key)
            written += 1
    return written


def main() -> int:
    if not LEARNED_DIR.is_dir():
        return 0
    files = sorted(p for p in LEARNED_DIR.glob("*.md") if p.is_file())
    if not files:
        return 0
    APPLIED_DIR.mkdir(parents=True, exist_ok=True)

    grouped: dict[str, list[dict]] = {}
    for f in files:
        for row in _parse(f):
            grouped.setdefault(row["predicate"], []).append(row)

    total = 0
    for predicate, rows in grouped.items():
        jsonl = EDGES_DIR / f"{predicate}.jsonl"
        n = _append(jsonl, rows)
        total += n
        print(f"  {predicate}: appended {n} new edge(s)")

    # Move applied files to applied/ (preserve, never delete)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    for f in files:
        dst = APPLIED_DIR / f"{stamp}_{f.name}"
        shutil.move(str(f), str(dst))

    print(f"[learned] drained {len(files)} file(s), {total} new edge(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
