"""Restore the materialized Kuzu DB from the committed distributable artifact.

Decompresses `public/.kg-dist/kg.kuzu.xz` into `.kg/kg.kuzu` so a fresh
checkout has a queryable graph in ~1s, without running `build.py`.

Called by bootstrap when `.kg/kg.kuzu` is missing. Exits non-zero (without
touching anything) when the artifact is absent, so the caller can fall back
to a full build.

Usage:
    python .kg/scripts/unpack.py           # restore if kg.kuzu missing
    python .kg/scripts/unpack.py --force   # overwrite an existing kg.kuzu
"""

from __future__ import annotations

import argparse
import lzma
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
KG_ROOT = SCRIPTS_DIR.parent
PUBLIC_ROOT = KG_ROOT.parent

DB = KG_ROOT / "kg.kuzu"
ARTIFACT = PUBLIC_ROOT / ".kg-dist" / "kg.kuzu.xz"

CHUNK = 1 << 20  # 1 MiB


def main() -> int:
    ap = argparse.ArgumentParser(description="Restore kg.kuzu from .kg-dist/kg.kuzu.xz")
    ap.add_argument("--force", action="store_true", help="overwrite an existing kg.kuzu")
    args = ap.parse_args()

    if not ARTIFACT.exists():
        print(f"[unpack] no prebuilt artifact at {ARTIFACT}", file=sys.stderr)
        return 2
    if DB.exists() and not args.force:
        print(f"[unpack] {DB} already exists — skipping (use --force to overwrite)")
        return 0

    tmp = DB.with_suffix(".kuzu.tmp")
    print(f"[unpack] restoring {DB} from {ARTIFACT} ({ARTIFACT.stat().st_size / 1e6:.2f} MB) ...")
    with lzma.open(ARTIFACT, "rb") as src, tmp.open("wb") as out:
        while True:
            chunk = src.read(CHUNK)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(DB)  # atomic — never leave a half-written kg.kuzu

    print(f"[unpack] done: {DB} ({DB.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
