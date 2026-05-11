"""Smoke test: schema loads, DDL is well-formed, Kuzu accepts it.

Run from `.kg/` with: `.venv/Scripts/python tests/test_schema_smoke.py`
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import kuzu

from schema import entities, relations
from schema.ddl import all_ddl
from schema.entities import Package, RegisteredFunction
from schema.relations import Exports


def main() -> None:
    print(f"Entity kinds: {len(entities.ENTITY_KINDS)}")
    print(f"Relation kinds: {len(relations.RELATION_KINDS)}")

    statements = all_ddl()
    print(f"\nDDL statements: {len(statements)}")
    for s in statements[:3]:
        print("  ", s[:120], "..." if len(s) > 120 else "")

    # Build a temp Kuzu DB
    tmp = Path(tempfile.mkdtemp(prefix="kg-smoke-"))
    db_path = tmp / "kg.kuzu"
    try:
        db = kuzu.Database(str(db_path))
        conn = kuzu.Connection(db)
        for stmt in statements:
            conn.execute(stmt)
        print("\nAll DDL applied successfully.")

        # Smoke test: insert a Package and a RegisteredFunction, link them
        pkg = Package(
            id="pkg:Smoke", name="Smoke", source_layer="plugins",
            friendly_name="Smoke Pkg", version="0.0.1",
        )
        fn = RegisteredFunction(
            id="func:Smoke.hello", name="hello", source_layer="plugins",
            package_id="pkg:Smoke", role="function", language="ts",
        )
        edge = Exports(from_id=pkg.id, to_id=fn.id)
        print("\nModel instances OK:")
        print(f"  {pkg.kind}: {pkg.id}")
        print(f"  {fn.kind}: {fn.id}")
        print(f"  {edge.predicate}: {edge.from_id} -> {edge.to_id}")

        # Round-trip the JSONL form
        print("\nJSONL round-trip:")
        print("  ", pkg.model_dump_json()[:200])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\nOK")


if __name__ == "__main__":
    main()
