"""One-shot migration: strip the `public/` prefix from every stored path
inside the canonical JSONL store and the LLM enrichment outputs.

Run **once** after moving `.kg/` from `reddata/` to `reddata/public/`. After
this runs and `build.py` rebuilds Kuzu, the graph will use repo-relative
paths anchored at `public/` (e.g. `packages/Bio/src/foo.ts`,
`file:packages/Bio/src/foo.ts`).

Idempotent: rerunning produces no change once the prefix is gone.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
ENRICH = DATA / "enrichment"

PREFIX = "public/"
FILE_ID_PREFIX = "file:public/"


def _strip(s: str) -> str:
    if s.startswith(PREFIX):
        return s[len(PREFIX):]
    if s.startswith(FILE_ID_PREFIX):
        return "file:" + s[len(FILE_ID_PREFIX):]
    return s


def _walk_obj(obj):
    if isinstance(obj, dict):
        return {k: _walk_obj(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_walk_obj(v) for v in obj]
    if isinstance(obj, str):
        return _strip(obj)
    return obj


def _migrate_jsonl(path: Path) -> int:
    if not path.is_file():
        return 0
    changed = 0
    out_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            out_lines.append(line)
            continue
        obj = json.loads(line)
        new_obj = _walk_obj(obj)
        new_line = json.dumps(new_obj, separators=(",", ":"), ensure_ascii=False)
        if new_line != line:
            changed += 1
        out_lines.append(new_line)
    path.write_text("\n".join(out_lines) + "\n", encoding="utf-8", newline="\n")
    return changed


def _migrate_json(path: Path) -> int:
    if not path.is_file():
        return 0
    obj = json.loads(path.read_text(encoding="utf-8"))
    new_obj = _walk_obj(obj)
    if obj == new_obj:
        return 0
    path.write_text(json.dumps(new_obj, indent=2, ensure_ascii=False),
                    encoding="utf-8")
    return 1


def main() -> int:
    total = 0
    for sub in ("nodes", "edges"):
        for jl in sorted((DATA / sub).glob("*.jsonl")):
            n = _migrate_jsonl(jl)
            if n:
                print(f"  {jl.relative_to(ROOT)}: {n} rows updated")
                total += n
    if ENRICH.is_dir():
        for jf in sorted(ENRICH.rglob("*.json")):
            n = _migrate_json(jf)
            if n:
                print(f"  {jf.relative_to(ROOT)}: rewritten")
                total += n
    print(f"done. {total} units rewritten.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
