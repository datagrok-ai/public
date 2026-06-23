"""Ad-hoc Cypher CLI for the local Kuzu DB.

Examples:
    py query.py "MATCH (p:Package) RETURN p.name, p.version LIMIT 10"
    py query.py --table "MATCH (p:Package)-[:EXPORTS]->(f:RegisteredFunction) WHERE p.name = 'Chem' RETURN f.name, f.role LIMIT 20"
    py query.py --demos       # run a few canned diagnostic queries
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent   # .kg/ (this file lives in .kg/scripts/)
DB_PATH = ROOT / "kg.kuzu"


DEMO_QUERIES: list[tuple[str, str]] = [
    ("Counts by entity kind",
     "MATCH (n) RETURN label(n) AS kind, count(*) AS n ORDER BY n DESC;"),
    ("Counts by edge predicate",
     "MATCH ()-[r]->() RETURN label(r) AS predicate, count(*) AS n ORDER BY n DESC;"),
    ("Top 10 functions by tag count",
     "MATCH (f:RegisteredFunction) WHERE f.tags IS NOT NULL "
     "RETURN f.name, size(f.tags) AS tag_count ORDER BY tag_count DESC LIMIT 10;"),
    ("All apps in Chem",
     "MATCH (p:Package {name: 'Chem'})-[:EXPORTS]->(f:RegisteredFunction {role: 'app'}) "
     "RETURN f.name, f.friendly_name;"),
    ("Functions per role",
     "MATCH (f:RegisteredFunction) WHERE f.role IS NOT NULL "
     "RETURN f.role AS role, count(*) AS n ORDER BY n DESC;"),
]


def render(rows: list[tuple], headers: list[str], *, table: bool = True) -> None:
    if not rows:
        print("(no rows)")
        return
    if not table:
        for r in rows:
            print(" | ".join(_fmt(v) for v in r))
        return
    cols = [list(map(_fmt, [headers[i]] + [r[i] for r in rows])) for i in range(len(headers))]
    widths = [max(len(s) for s in c) for c in cols]
    line = "  ".join(c[0].ljust(w) for c, w in zip(cols, widths))
    print(line)
    print("  ".join("-" * w for w in widths))
    for i in range(len(rows)):
        print("  ".join(c[i + 1].ljust(w) for c, w in zip(cols, widths)))


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, list):
        return ",".join(_fmt(x) for x in v)
    return str(v)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cypher", nargs="?", help="Cypher query")
    ap.add_argument("--demos", action="store_true", help="Run canned demo queries")
    ap.add_argument("--no-table", action="store_true", help="Plain pipe-separated output")
    args = ap.parse_args()

    if not DB_PATH.exists():
        print(f"DB not found at {DB_PATH}. Run qq.py first — it restores kg.kuzu from kg.kuzu.xz.", file=sys.stderr)
        return 2

    import kuzu
    db = kuzu.Database(str(DB_PATH))
    conn = kuzu.Connection(db)

    queries: list[tuple[str, str]] = []
    if args.demos:
        queries = DEMO_QUERIES
    elif args.cypher:
        queries = [("query", args.cypher)]
    else:
        ap.print_help()
        return 0

    for label, q in queries:
        print(f"\n# {label}")
        print(f"  {q}")
        print()
        try:
            res = conn.execute(q)
            headers = res.get_column_names()
            rows: list[tuple] = []
            while res.has_next():
                rows.append(tuple(res.get_next()))
            render(rows, headers, table=not args.no_table)
        except Exception as e:
            print(f"  ! {e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
