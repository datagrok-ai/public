"""Promote Datagrok semantic types to first-class nodes and link them to
the functions that detect, consume, or produce them.

Sources scanned (per RegisteredFunction):
  - `meta.semType: <X>`         → Detects (when role==semTypeDetector) OR Produces (otherwise)
  - `inputs[i].options.semType` → ConsumesSemtype
  - `outputs[i].options.semType` → ProducesSemtype
  - `meta.columnTags`           → RequiresColumnTag (one per `key=value` pair).
                                  When key is `quality`, the pair also implies
                                  ConsumesSemtype (since `quality=X` means the
                                  function only fires for columns of semType X).

Cardinality: ~30 distinct semantic types, ~700+ edges (post-2026-05).
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import SemanticType
from schema.relations import (
    ConsumesSemtype, DetectsSemtype, ProducesSemtype, RequiresColumnTag,
)

from . import _common as c

EXTRACTOR_NAME = "semtypes"


def _semtype_id(value: str) -> str:
    return f"semtype:{value}"


def _normalize_semtype(raw: str | None) -> str | None:
    """Strip `; description: ...`, `; units: ...`, `; caption: ...` pollution
    that crept in when `_parse_inline_options` accepted the rest of the
    `//input:` line as part of the value. The semType is always the first
    semicolon-separated token, trimmed."""
    if not raw or not isinstance(raw, str):
        return None
    s = raw.strip()
    # Cut on `;` (description/units/caption) and on `,` (rare - aliasing)
    s = re.split(r"[;,]", s, maxsplit=1)[0].strip()
    if not s:
        return None
    return s


_TAG_KV_RE = re.compile(r"\s*([A-Za-z][\w.\-]*)\s*=\s*([^,]+?)\s*(?:,|$)")


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    nodes_dir = c.REPO_ROOT / ".kg" / "data" / "nodes"
    rf_path = nodes_dir / "RegisteredFunction.jsonl"
    if not rf_path.is_file():
        print(f"[{EXTRACTOR_NAME}] no RegisteredFunction.jsonl yet")
        return

    detected: dict[str, set[str]] = defaultdict(set)   # semtype → func_ids
    consumed: dict[str, set[str]] = defaultdict(set)
    produced: dict[str, set[str]] = defaultdict(set)
    # (fid, key, value) tuples for RequiresColumnTag — the value side becomes
    # a SemanticType node so equality queries on quality= work uniformly.
    requires: list[tuple[str, str, str]] = []

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
        fid = row["id"]
        meta = row.get("meta") or {}
        role = (row.get("role") or "").lower()

        # meta.semType
        m_st = _normalize_semtype(meta.get("semType") or meta.get("semtype"))
        if m_st:
            if role == "semtypedetector":
                detected[m_st].add(fid)
            else:
                produced[m_st].add(fid)

        # inputs / outputs
        for inp in row.get("inputs") or []:
            opts = inp.get("options") or {}
            st = _normalize_semtype(opts.get("semType") or opts.get("semtype"))
            if st:
                consumed[st].add(fid)
        for out in row.get("outputs") or []:
            opts = out.get("options") or {}
            st = _normalize_semtype(opts.get("semType") or opts.get("semtype"))
            if st:
                produced[st].add(fid)

        # meta.columnTags — `quality=Macromolecule, units=helm` → one
        # RequiresColumnTag per pair; quality= also implies ConsumesSemtype
        col_tags_raw = meta.get("columnTags") or meta.get("columntags")
        if col_tags_raw and isinstance(col_tags_raw, str):
            for km in _TAG_KV_RE.finditer(col_tags_raw):
                k = km.group(1).strip()
                v = km.group(2).strip().rstrip(",").strip()
                if not (k and v):
                    continue
                requires.append((fid, k, v))
                if k.lower() == "quality":
                    consumed[v].add(fid)

    all_st = set(detected) | set(consumed) | set(produced)
    # Also seed nodes for tag values so RequiresColumnTag has typed targets.
    for _fid, _k, v in requires:
        all_st.add(v)

    n_edges = 0
    for st in sorted(all_st):
        sid = _semtype_id(st)
        bundle.add(SemanticType(
            id=sid, name=st,
            source_layer=SourceLayer.SYNTHETIC,
            semtype_value=st,
            detected_by_count=len(detected.get(st, ())),
            consumed_by_count=len(consumed.get(st, ())),
            description=None,
        ))
        for fid in detected.get(st, ()):
            bundle.add(DetectsSemtype(from_id=fid, to_id=sid,
                                      derived_by=Provenance.ANNOTATION))
            n_edges += 1
        for fid in consumed.get(st, ()):
            bundle.add(ConsumesSemtype(from_id=fid, to_id=sid,
                                       derived_by=Provenance.ANNOTATION))
            n_edges += 1
        for fid in produced.get(st, ()):
            bundle.add(ProducesSemtype(from_id=fid, to_id=sid,
                                       derived_by=Provenance.ANNOTATION))
            n_edges += 1
    n_req = 0
    for fid, k, v in requires:
        bundle.add(RequiresColumnTag(
            from_id=fid, to_id=_semtype_id(v),
            tag_key=k, tag_value=v,
            derived_by=Provenance.ANNOTATION,
        ))
        n_req += 1

    print(f"[{EXTRACTOR_NAME}] {len(all_st)} semantic types, {n_edges} sem edges, "
          f"{n_req} RequiresColumnTag edges")
