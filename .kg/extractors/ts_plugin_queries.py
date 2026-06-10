"""Extract DataQuery entities from `packages/*/queries/*.sql`.

A `.sql` file may contain multiple queries, each delimited by `--end`.
Each query block has a `--name:` header. Annotations follow the same
`--key: value` mini-format as scripts.

Emits:
  - DataQuery
  - HasQuery (Package -> DataQuery)
  - UsesConnection (DataQuery -> DataConnection)  (resolved by name)
  - MentionsTicket if the description mentions a GROK ticket
"""

from __future__ import annotations

import re
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import DataQuery
from schema.relations import HasQuery, UsesConnection

from . import _common as c

EXTRACTOR_NAME = "ts_plugin_queries"

# Split a `.sql` file on `--end` lines (one per query).
END_LINE = re.compile(r"^\s*--\s*end\s*$", re.IGNORECASE | re.MULTILINE)


def _split_blocks(text: str) -> list[str]:
    parts = END_LINE.split(text)
    return [p for p in (s.strip() for s in parts) if p]


def _extract_query(block: str, file_path: Path, pkg_name: str, pkg_node_id: str, bundle) -> bool:
    # The annotation header: contiguous leading lines starting with `--`
    lines = block.splitlines()
    head: list[str] = []
    for ln in lines:
        if ln.strip().startswith("--"):
            head.append(ln)
        else:
            break
    if not head:
        return False
    ann = c.parse_annotation_block(head)
    name = ann.get("name")
    if not name:
        return False
    qid = c.query_id(pkg_name, name)
    conn_name = ann.get("raw", {}).get("connection")
    q = DataQuery(
        id=qid, name=name, source_layer=SourceLayer.PLUGINS,
        package_id=pkg_node_id,
        connection_name=conn_name,
        inputs=ann.get("inputs", []),
        tags=ann.get("tags", []),
        cache=str(ann.get("meta", {}).get("cache", "")).lower() in {"true", "client", "server"},
        description=ann.get("description"),
        description_provenance=Provenance.ANNOTATION if ann.get("description") else None,
        paths=[c.file_ref(file_path, role="definition")],
        extras={"friendly_name": ann.get("friendly_name"),
                "top_menu": ann.get("top_menu"),
                "meta": ann.get("meta") or None},
    )
    bundle.add(q)
    bundle.add(HasQuery(from_id=pkg_node_id, to_id=qid,
                        derived_by=Provenance.ANNOTATION))
    if conn_name:
        # Two patterns: bare 'Foo' (same package) or 'Pkg:Foo' (cross-package)
        if ":" in conn_name:
            tgt_pkg, tgt = conn_name.split(":", 1)
            bundle.add(UsesConnection(
                from_id=qid, to_id=c.connection_id(tgt_pkg, tgt),
                derived_by=Provenance.ANNOTATION))
        else:
            bundle.add(UsesConnection(
                from_id=qid, to_id=c.connection_id(pkg_name, conn_name),
                derived_by=Provenance.ANNOTATION))
    return True


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return
    total = 0
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        qdir = pkg_dir / "queries"
        if not qdir.is_dir():
            continue
        for sql_file in sorted(qdir.glob("*.sql")):
            try:
                text = sql_file.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            for blk in _split_blocks(text):
                if _extract_query(blk, sql_file, pkg_dir.name, c.pkg_id(pkg_dir.name), bundle):
                    total += 1
    print(f"[{EXTRACTOR_NAME}] {total} queries")
