"""Classify each Feature deterministically:

  - `is_user_facing`   — True iff any PART_OF_FEATURE member is a
                         RegisteredFunction with a user-facing role
                         (app, viewer, widget, panel, dashboard, fileViewer,
                         fileHandler, fileExporter, transform, editor,
                         semTypeDetector, cellRenderer, valueEditor,
                         searchProvider) OR the Feature has any DocPage /
                         Tutorial member.

  - `interaction_kind` — single category derived from the dominant member role:
                         'app' | 'viewer' | 'widget-panel' | 'transform' |
                         'data-import' | 'api' | 'infra' | 'data-query' |
                         'script' | 'mixed'

This is a synthesis pass — reads existing JSONL, writes a Feature.jsonl
overlay that the build script merges. Run AFTER `feature_clustering` has
emitted Feature + PART_OF_FEATURE edges.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from . import _common as c

EXTRACTOR_NAME = "feature_classify"

NODES_DIR = c.REPO_ROOT / ".kg" / "data" / "nodes"
EDGES_DIR = c.REPO_ROOT / ".kg" / "data" / "edges"

USER_FACING_ROLES = {
    "app", "viewer", "widget", "panel", "dashboard", "fileViewer",
    "fileHandler", "fileExporter", "transform", "editor", "semTypeDetector",
    "cellRenderer", "valueEditor", "searchProvider", "tooltip",
    "moleculeSketcher", "converter", "filter", "hitTriageDataSource",
    "hitTriageFunction",
}
INFRA_ROLES = {"init", "autostart"}


def _read_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _classify(roles: Counter, has_doc: bool, has_tutorial: bool,
              has_query: bool, has_script: bool) -> tuple[bool, str]:
    """Return (is_user_facing, interaction_kind)."""
    user_facing = (
        any(r in USER_FACING_ROLES for r in roles)
        or has_doc or has_tutorial
    )

    # Pick interaction_kind by precedence
    if "app" in roles:                 kind = "app"
    elif "viewer" in roles:            kind = "viewer"
    elif "widget" in roles or "panel" in roles or "dashboard" in roles or "tooltip" in roles:
        kind = "widget-panel"
    elif "transform" in roles or "editor" in roles:  kind = "transform"
    elif "fileHandler" in roles or "fileViewer" in roles or "fileExporter" in roles:
        kind = "data-import"
    elif "cellRenderer" in roles or "semTypeDetector" in roles or "converter" in roles or "valueEditor" in roles:
        kind = "api"
    elif INFRA_ROLES & set(roles) and not user_facing:
        kind = "infra"
    elif has_query and not roles:      kind = "data-query"
    elif has_script and not roles:     kind = "script"
    else:                              kind = "mixed"

    return user_facing, kind


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    feats = _read_jsonl(NODES_DIR / "Feature.jsonl")
    if not feats:
        print(f"[{EXTRACTOR_NAME}] no features yet")
        return
    rfs = {n["id"]: n for n in _read_jsonl(NODES_DIR / "RegisteredFunction.jsonl")}
    docs_index = {n["id"] for n in _read_jsonl(NODES_DIR / "DocPage.jsonl")}
    tut_index = {n["id"] for n in _read_jsonl(NODES_DIR / "Tutorial.jsonl")}
    query_index = {n["id"] for n in _read_jsonl(NODES_DIR / "DataQuery.jsonl")}
    script_index = {n["id"] for n in _read_jsonl(NODES_DIR / "Script.jsonl")}

    # Walk PART_OF_FEATURE: feat_id -> [(member_id, role-or-None)]
    pof_by_feat: dict[str, list[tuple[str, str]]] = {}
    for row in _read_jsonl(EDGES_DIR / "PART_OF_FEATURE.jsonl"):
        f = row["to_id"]
        m = row["from_id"]
        pof_by_feat.setdefault(f, []).append((m, row.get("role")))

    overlay: list[dict] = []
    counter = Counter()
    for feat in feats:
        members = pof_by_feat.get(feat["id"], [])
        roles_seen: Counter = Counter()
        has_doc = has_tutorial = has_query = has_script = False
        for mid, edge_role in members:
            if mid in docs_index:    has_doc = True
            elif mid in tut_index:   has_tutorial = True
            elif mid in query_index: has_query = True
            elif mid in script_index: has_script = True
            elif mid in rfs:
                role = rfs[mid].get("role")
                # tags like 'widgets,panel' — split on comma
                if role:
                    for r in role.split(","):
                        roles_seen[r.strip()] += 1
        is_uf, kind = _classify(roles_seen, has_doc, has_tutorial, has_query, has_script)
        counter[(is_uf, kind)] += 1

        # Overlay: re-emit the Feature node with the new fields populated.
        # The Bundle merges by id, so this overwrites the prior row.
        out = dict(feat)
        out["is_user_facing"] = is_uf
        out["interaction_kind"] = kind
        # Keep the rest as-is. Bundle.add expects an Entity instance though,
        # so we hand-roll a JSONL append below.
        overlay.append(out)

    # Write overlay directly (bypassing Bundle since we're updating, not adding)
    feat_path = NODES_DIR / "Feature.jsonl"
    by_id = {row["id"]: row for row in feats}
    for row in overlay:
        by_id[row["id"]] = row
    with feat_path.open("w", encoding="utf-8", newline="\n") as f:
        for row in by_id.values():
            f.write(json.dumps(row, separators=(",", ":"), ensure_ascii=False) + "\n")

    print(f"[{EXTRACTOR_NAME}] classified {len(overlay)} features:")
    for (uf, k), n in sorted(counter.items(), key=lambda x: -x[1])[:12]:
        print(f"  {('user' if uf else 'infra'):<5}  {k:<14}  {n:>4}")
