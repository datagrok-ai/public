"""Promote function tags to first-class nodes.

Reads `RegisteredFunction.jsonl` (already produced by `ts_plugin_package`),
collects every distinct string in `tags[]` AND in `meta.role` (when role
matches a known platform-dispatched function role), emits one `Tag` node
per tag, and one `HasTag(RegisteredFunction → Tag)` edge per (function, tag).

Why both `tags[]` and `meta.role`? Plugin authors are inconsistent. Helm's
`editMoleculeCell` declares `tags: ['cellEditor']`, while Chem's declares
`meta.role: 'cellEditor'` only. A query that filters by `:Tag {name:'cellEditor'}`
must hit BOTH. We canonicalize: every dispatched-role string (cellEditor,
cellRenderer, panel, ...) is reflected as a Tag, regardless of which field
the source declares it in.

Lets queries like
  MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'cellEditor'})
  RETURN rf.id, rf.package_id, rf.meta
work without scanning every RF's tag array OR meta dict.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import Tag
from schema.relations import HasTag

from . import _common as c

EXTRACTOR_NAME = "tags"

# Platform-dispatched roles. Whenever `meta.role` is one of these, also
# emit a Tag of the same name — so role-based and tag-based queries
# converge on the same answer.
_DISPATCHED_ROLES = {
    "app", "init", "autostart", "panel", "viewer", "widget", "filter",
    "cellRenderer", "cellEditor", "fileViewer", "fileExporter", "fileHandler",
    "semTypeDetector", "valueEditor", "dashboard", "transform",
    "function", "converter", "notationProvider", "notationRefiner",
    "packageSettingsEditor", "appTreeBrowser", "model", "demo",
}


def _tag_id(name: str) -> str:
    return f"tag:{name}"


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    rf_path = c.REPO_ROOT / ".kg" / "data" / "nodes" / "RegisteredFunction.jsonl"
    if not rf_path.is_file():
        print(f"[{EXTRACTOR_NAME}] no RegisteredFunction.jsonl yet")
        return

    tag_to_funcs: dict[str, set[str]] = defaultdict(set)
    n_role_synth = 0
    for line in rf_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if package_filter:
            pkg_id = row.get("package_id") or ""
            if pkg_id.split(":", 1)[-1] not in package_filter:
                continue
        fid = row.get("id")
        if not fid:
            continue
        # 1. Tags as declared in source.
        for tag in (row.get("tags") or []):
            t = (tag or "").strip()
            if t:
                tag_to_funcs[t].add(fid)
        # 2. Synthesize from role (top-level + meta.role).
        meta = row.get("meta") or {}
        for src in (row.get("role"), meta.get("role")):
            if not isinstance(src, str):
                continue
            r = src.strip()
            if r and r in _DISPATCHED_ROLES and fid not in tag_to_funcs[r]:
                tag_to_funcs[r].add(fid)
                n_role_synth += 1

    n_tags = n_edges = 0
    for tag_name in sorted(tag_to_funcs):
        tid = _tag_id(tag_name)
        bundle.add(Tag(
            id=tid, name=tag_name,
            source_layer=SourceLayer.SYNTHETIC,
            function_count=len(tag_to_funcs[tag_name]),
        ))
        n_tags += 1
        for fid in sorted(tag_to_funcs[tag_name]):
            bundle.add(HasTag(from_id=fid, to_id=tid,
                              derived_by=Provenance.ANNOTATION))
            n_edges += 1

    print(f"[{EXTRACTOR_NAME}] {n_tags} tags, {n_edges} HasTag edges "
          f"({n_role_synth} synthesized from role)")
