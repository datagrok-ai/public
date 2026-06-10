"""Pack the materialized Kuzu DB into the committed, distributable artifact.

Compresses `.kg/kg.kuzu` (≈560 MB, mostly zero-padded pages) into
`public/.kg-dist/kg.kuzu.xz` (≈3 MB) using stdlib lzma — no external
binary, so it runs the same on Windows, macOS, and Linux.

The `.xz` lives in the *public* repo (outside the `.kg` submodule), so a
plain clone ships a ready-to-use graph that `unpack.py` restores in ~1s,
skipping the multi-minute `build.py`.

Usage:
    python .kg/scripts/pack.py            # pack kg.kuzu -> .kg-dist/kg.kuzu.xz
    python .kg/scripts/pack.py --preset 9 # stronger compression (slower)

Pack first CHECKPOINTs the DB (folding any write-ahead log into the main
file and removing the .wal), so the resulting .xz is a self-contained
snapshot regardless of whether the query server was running.
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
WAL = KG_ROOT / "kg.kuzu.wal"
DIST_DIR = PUBLIC_ROOT / ".kg-dist"
ARTIFACT = DIST_DIR / "kg.kuzu.xz"

CHUNK = 1 << 20  # 1 MiB


def _checkpoint() -> None:
    """Fold the WAL into kg.kuzu and drop it, so we pack a clean snapshot."""
    import kuzu
    db = kuzu.Database(str(DB))
    conn = kuzu.Connection(db)
    conn.execute("CHECKPOINT")
    conn.close()
    db.close()


def main() -> int:
    ap = argparse.ArgumentParser(description="Pack kg.kuzu into .kg-dist/kg.kuzu.xz")
    ap.add_argument("--preset", type=int, default=6, help="lzma preset 0-9 (default 6)")
    args = ap.parse_args()

    if not DB.exists():
        print(f"[pack] {DB} not found — build the graph first (build.py)", file=sys.stderr)
        return 1

    print("[pack] checkpointing (folding WAL into kg.kuzu) ...")
    _checkpoint()
    if WAL.exists():
        print(f"[pack] {WAL} still present after checkpoint — stop the query server "
              f"(qq.py --stop) and retry", file=sys.stderr)
        return 1

    DIST_DIR.mkdir(parents=True, exist_ok=True)
    src_size = DB.stat().st_size
    print(f"[pack] compressing {DB} ({src_size / 1e6:.1f} MB) -> {ARTIFACT} ...")

    with DB.open("rb") as src, lzma.open(ARTIFACT, "wb", preset=args.preset) as out:
        while True:
            chunk = src.read(CHUNK)
            if not chunk:
                break
            out.write(chunk)

    out_size = ARTIFACT.stat().st_size
    print(f"[pack] done: {out_size / 1e6:.2f} MB ({src_size / out_size:.0f}x smaller). "
          f"Commit {ARTIFACT.relative_to(PUBLIC_ROOT)} in the public repo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
