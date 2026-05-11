"""Feature clustering enricher.

For each Package, gather every entity that belongs to it (functions, scripts,
queries, tutorials, package-level docs) plus its CHANGELOG summaries, write
a single context JSON, and ask an LLM to group them into Features.

Three-step contract:
  1. prepare()  -> .kg/data/enrichment/feature_clustering/inputs/<Package>.json
  2. (Claude Code agent or any LLM) -> outputs/<Package>.json with Feature clusters
  3. apply()    -> emit Feature, HasFeature, PartOfFeature JSONL slices
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from schema.base import Provenance, SourceLayer
from schema.entities import Feature
from schema.relations import HasFeature, PartOfFeature, Documents

from extractors import _common as c

NAME = "feature_clustering"

DATA_DIR = c.REPO_ROOT / ".kg" / "data"
ENRICH_DIR = DATA_DIR / "enrichment" / NAME
INPUTS_DIR = ENRICH_DIR / "inputs"
OUTPUTS_DIR = ENRICH_DIR / "outputs"
INSTRUCTIONS_PATH = ENRICH_DIR / "PROMPT.md"


# ---------------------------------------------------------------------------
# JSONL loaders (read canonical store)
# ---------------------------------------------------------------------------

def _load_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    out: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            out.append(json.loads(line))
    return out


def _load_all() -> dict[str, list[dict]]:
    """Read every node JSONL into {kind: [rows]}."""
    nodes_dir = DATA_DIR / "nodes"
    out: dict[str, list[dict]] = {}
    if not nodes_dir.is_dir():
        return out
    for f in nodes_dir.glob("*.jsonl"):
        out[f.stem] = _load_jsonl(f)
    return out


def _load_edges(predicate: str) -> list[dict]:
    return _load_jsonl(DATA_DIR / "edges" / f"{predicate}.jsonl")


# ---------------------------------------------------------------------------
# prepare(): build per-package context
# ---------------------------------------------------------------------------

def _truncate(s: str | None, n: int) -> str | None:
    if not s:
        return s
    s = s.strip().replace("\r", "")
    if len(s) <= n:
        return s
    return s[:n].rsplit(" ", 1)[0] + "…"


def _scrub_function(f: dict) -> dict:
    """Compact a RegisteredFunction row for the LLM context."""
    return {
        "id": f["id"],
        "name": f.get("name"),
        "friendly_name": f.get("friendly_name"),
        "role": f.get("role"),
        "tags": f.get("tags") or [],
        "description": _truncate(f.get("description"), 220),
        "input_types": [i.get("type") for i in (f.get("inputs") or [])],
        "output_types": [o.get("type") for o in (f.get("outputs") or [])],
        "meta": {k: v for k, v in (f.get("meta") or {}).items()
                 if k in {"role", "domain", "browsePath", "fileViewer", "ext", "semType"}},
        "help_url": f.get("help_url"),
        "language": f.get("language"),
    }


def _scrub_script(s: dict) -> dict:
    return {"id": s["id"], "name": s.get("name"),
            "language": s.get("language"),
            "description": _truncate(s.get("description"), 200)}


def _scrub_query(q: dict) -> dict:
    return {"id": q["id"], "name": q.get("name"),
            "connection": q.get("connection_name"),
            "description": _truncate(q.get("description"), 200),
            "tags": q.get("tags") or []}


def _scrub_tutorial(t: dict) -> dict:
    return {"id": t["id"], "name": t.get("name"),
            "class_name": t.get("class_name"),
            "description": _truncate(t.get("description"), 250),
            "step_count": t.get("step_count"),
            "help_url": t.get("help_url")}


def _scrub_doc(d: dict) -> dict:
    return {"id": d["id"], "title": d.get("title"),
            "url": (d.get("extras") or {}).get("url"),
            "kind": d.get("doc_kind"), "audience": d.get("audience")}


def _scrub_changelog_entry(e: dict, *, max_per_pkg: int) -> dict | None:
    return {"id": e["id"], "version": e.get("version"),
            "released_at": e.get("released_at"),
            "ticket": e.get("ticket_id"),
            "text": _truncate(e.get("description"), 200)}


def prepare(packages: list[str] | None = None,
            *, max_changelog_per_pkg: int = 30,
            max_docs_per_pkg: int = 25) -> list[Path]:
    """Write per-package input JSONs. Returns the list of input file paths."""
    INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    nodes = _load_all()
    pkgs = nodes.get("Package", [])

    by_pkg_funcs: dict[str, list[dict]] = defaultdict(list)
    for f in nodes.get("RegisteredFunction", []):
        by_pkg_funcs[f.get("package_id")].append(f)
    by_pkg_scripts: dict[str, list[dict]] = defaultdict(list)
    for s in nodes.get("Script", []):
        by_pkg_scripts[s.get("package_id")].append(s)
    by_pkg_queries: dict[str, list[dict]] = defaultdict(list)
    for q in nodes.get("DataQuery", []):
        by_pkg_queries[q.get("package_id")].append(q)
    by_pkg_tutorials: dict[str, list[dict]] = defaultdict(list)
    for t in nodes.get("Tutorial", []):
        by_pkg_tutorials[t.get("package_id")].append(t)
    by_pkg_changelog: dict[str, list[dict]] = defaultdict(list)
    for e in nodes.get("ChangelogEntry", []):
        by_pkg_changelog[e.get("package_id")].append(e)
    by_pkg_props: dict[str, list[dict]] = defaultdict(list)
    for pp in nodes.get("PackageProperty", []):
        by_pkg_props[pp.get("package_id")].append(pp)

    # ---- DocPages relevant to each package ---------------------------------
    # Three signals, summed to a relevance score:
    #   (a) DocPage reachable from any package-owned entity via DOCUMENTS edges  (+10)
    #   (b) DocPage URL contains the package folder name OR its category slug   (+5)
    #   (c) DocPage title shares non-trivial words with any package function    (+1 each, capped at +5)
    # Plus: README/META/CLAUDE.md docs always attached to their package.
    by_pkg_docs = _select_docs_per_package(
        nodes.get("DocPage", []),
        documents_edges=_load_edges("DOCUMENTS"),
        packages=pkgs,
        by_pkg_funcs=by_pkg_funcs,
        max_per_pkg=max_docs_per_pkg,
    )

    written: list[Path] = []
    for pkg in pkgs:
        folder = pkg["name"]
        if packages and folder not in packages:
            continue
        funcs = by_pkg_funcs.get(pkg["id"], [])
        scripts = by_pkg_scripts.get(pkg["id"], [])
        queries = by_pkg_queries.get(pkg["id"], [])
        tutorials = by_pkg_tutorials.get(pkg["id"], [])
        changelog = sorted(by_pkg_changelog.get(pkg["id"], []),
                           key=lambda e: (e.get("version") != "v.next",
                                          e.get("released_at") or ""),
                           reverse=True)[:max_changelog_per_pkg]
        # Distinct DocPages, capped
        seen_doc_ids: set[str] = set()
        docs: list[dict] = []
        for d in by_pkg_docs.get(pkg["id"], []):
            if d["id"] in seen_doc_ids:
                continue
            seen_doc_ids.add(d["id"])
            docs.append(d)
            if len(docs) >= max_docs_per_pkg:
                break

        ctx = {
            "package": {
                "id": pkg["id"],
                "name": pkg["name"],
                "friendly_name": pkg.get("friendly_name"),
                "version": pkg.get("version"),
                "category": pkg.get("category"),
                "description": _truncate(pkg.get("description"), 600),
                "extras": {k: v for k, v in (pkg.get("extras") or {}).items()
                           if k in {"npm_name", "fullName", "browserFeatures"}},
            },
            "properties": [{"name": pp.get("name"), "type": pp.get("type"),
                            "choices": pp.get("choices"),
                            "default": pp.get("default_value")}
                           for pp in by_pkg_props.get(pkg["id"], [])],
            "functions":  [_scrub_function(f) for f in funcs],
            "scripts":    [_scrub_script(s) for s in scripts],
            "queries":    [_scrub_query(q) for q in queries],
            "tutorials":  [_scrub_tutorial(t) for t in tutorials],
            "docs":       [_scrub_doc(d) for d in docs],
            "changelog":  [_scrub_changelog_entry(e, max_per_pkg=max_changelog_per_pkg)
                           for e in changelog],
            "stats": {
                "functions":  len(funcs),
                "scripts":    len(scripts),
                "queries":    len(queries),
                "tutorials":  len(tutorials),
                "docs":       len(docs),
                "changelog":  len(by_pkg_changelog.get(pkg["id"], [])),
            },
        }
        out_path = INPUTS_DIR / f"{folder}.json"
        out_path.write_text(json.dumps(ctx, indent=2, ensure_ascii=False),
                            encoding="utf-8")
        written.append(out_path)
    print(f"[{NAME}] wrote {len(written)} input contexts to {INPUTS_DIR}")
    write_prompt()
    return written


# ---------------------------------------------------------------------------
# Prompt template (written once next to the inputs)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# DocPage relevance scoring
# ---------------------------------------------------------------------------

# Map common Datagrok package categories to URL substrings that hint at
# matching help content. Add as we discover more.
CATEGORY_URL_HINTS: dict[str, list[str]] = {
    "Cheminformatics":   ["/chem/", "/chemistry"],
    "Bioinformatics":    ["/bio/", "/peptides", "/macromolecule"],
    "Visualization":     ["/visualize/", "/viewers/"],
    "Database":          ["/access/databases", "/db-explorer"],
    "Compute":           ["/compute/", "/diff-studio"],
    "Connectors":        ["/access/databases"],
    "Tutorials":         ["/tutorial"],
    "Statistics":        ["/explore/"],
    "Demo":              ["/datagrok/"],
}

_STOPWORDS = {
    "the", "a", "an", "and", "or", "of", "to", "for", "with", "in", "on", "by",
    "as", "is", "be", "at", "this", "that", "from", "into", "via", "data",
    "table", "column", "value", "viewer", "widget", "panel", "function",
    "type", "show", "get", "set", "new", "all", "one", "two", "three",
    "demo", "test", "init", "load", "save", "list", "find", "result", "obj",
    "object", "name", "items", "options", "default",
}


def _word_tokens(s: str | None) -> set[str]:
    if not s:
        return set()
    raw = re.findall(r"[A-Za-z][A-Za-z0-9]+", s)
    return {w.lower() for w in raw if len(w) > 3 and w.lower() not in _STOPWORDS}


def _select_docs_per_package(
    all_docs: list[dict],
    *,
    documents_edges: list[dict],
    packages: list[dict],
    by_pkg_funcs: dict[str, list[dict]],
    max_per_pkg: int,
) -> dict[str, list[dict]]:
    """Score each (Package, DocPage) pair, return top-N docs per package.

    A DocPage's `extras.url` (canonical /help/... or /packages/<X>/...) is
    the matching surface; `extras.package` (set on README docs) gives a hard
    +50 bonus to ensure README always lands.
    """
    docs_by_id = {d["id"]: d for d in all_docs}

    # Per-package signals
    pkg_by_id: dict[str, dict] = {p["id"]: p for p in packages}
    pkg_keywords: dict[str, set[str]] = {}
    pkg_url_hints: dict[str, list[str]] = {}
    for p in packages:
        kws = {p["name"].lower()}
        if p.get("friendly_name"):
            kws |= _word_tokens(p["friendly_name"])
        if p.get("category"):
            kws.add(p["category"].lower())
        for f in by_pkg_funcs.get(p["id"], []):
            kws |= _word_tokens(f.get("name"))
            kws |= _word_tokens(f.get("friendly_name"))
        pkg_keywords[p["id"]] = {k for k in kws if k}
        hints = CATEGORY_URL_HINTS.get(p.get("category") or "", [])
        # Add the package folder name as a hint
        hints = list(hints) + [f"/{p['name'].lower()}", f"-{p['name'].lower()}"]
        pkg_url_hints[p["id"]] = hints

    # Signal 1: explicit Documents edges
    explicit: dict[str, set[str]] = defaultdict(set)
    for e in documents_edges:
        tgt = e.get("to_id", "")
        m = re.match(r"^(?:pkg|func|script|query|tutorial):([^:]+)", tgt)
        if not m:
            continue
        explicit[c.pkg_id(m.group(1))].add(e.get("from_id", ""))

    # Score every (pkg, doc) pair, then top-N
    out: dict[str, list[dict]] = defaultdict(list)
    for pkg_id, pkg in pkg_by_id.items():
        scores: list[tuple[int, dict]] = []
        kws = pkg_keywords[pkg_id]
        url_hints = pkg_url_hints[pkg_id]
        for d in all_docs:
            url = (d.get("extras") or {}).get("url", "") or ""
            url_l = url.lower()
            score = 0
            # README/META/CLAUDE: hard match by extras.package
            if (d.get("extras") or {}).get("package") == pkg["name"]:
                score += 50
            # Explicit Documents edge
            if d["id"] in explicit.get(pkg_id, ()):
                score += 10
            # URL-substring hint
            for h in url_hints:
                if h.lower() in url_l:
                    score += 5
                    break
            # Title-word overlap with package keywords
            title_tokens = _word_tokens(d.get("title")) | _word_tokens(d.get("name"))
            overlap = len(title_tokens & kws)
            if overlap:
                score += min(overlap, 5)
            if score > 0:
                scores.append((score, d))
        scores.sort(key=lambda x: (-x[0], x[1].get("id", "")))
        seen: set[str] = set()
        for sc, d in scores:
            if d["id"] in seen:
                continue
            seen.add(d["id"])
            out[pkg_id].append(d)
            if len(out[pkg_id]) >= max_per_pkg:
                break
    return out


PROMPT_BODY = """\
# Feature clustering for one Datagrok plugin

You are analyzing a single Datagrok plugin (package). Your job is to identify
**Features** — coherent, user-facing capabilities that this plugin provides.

A Feature is NOT one function. It is a *theme* that may be implemented by:
- one or more `RegisteredFunction`s (the user-callable functions)
- supporting `Script`s (server-side compute)
- backing `DataQuery`s (SQL / data sources)
- a `Tutorial` that walks a user through it
- one or more `DocPage`s that explain it
- changelog entries describing its evolution

## Input

You will receive a JSON object describing the package: its metadata, its
functions (with role/tags/inputs/outputs/description/help_url/meta), scripts,
queries, tutorials, doc pages, properties, and a recent changelog. Every entity
has an `id`. Use those IDs verbatim in your output.

## Your output (STRICT JSON, no other text)

```json
{
  "package_id": "pkg:<Folder>",
  "features": [
    {
      "id": "feature:<Folder>:<slug>",
      "name": "Human-readable feature name",
      "description": "1-2 sentences describing what this feature gives the user.",
      "category": "<short tag, e.g. chem-search | chem-rendering | chem-mpo | viz | ml | data-import | infra>",
      "confidence": 0.0,
      "members": [
        {"id": "func:Pkg:foo",       "weight": 1.0,  "role": "core"},
        {"id": "script:Pkg:bar",     "weight": 0.7,  "role": "core"},
        {"id": "doc:/help/.../X",    "weight": 0.9,  "role": "doc"},
        {"id": "tutorial:Pkg:Foo",   "weight": 0.95, "role": "doc"}
      ],
      "evidence_summary": "1 sentence on why these belong together (shared name root, shared menu path, shared tutorial, etc.)"
    }
  ],
  "feature_relations": [
    { "from": "feature:Pkg:foo", "to": "feature:Pkg:bar", "kind": "subfeature|sibling|related" }
  ],
  "unassigned_member_ids": ["..."]
}
```

## Rules

- `id` slugs are **kebab-case**, derived from the feature name. Example: `"Activity Cliffs"` → `feature:Chem:activity-cliffs`.
- A feature may contain entities of any kind in `members`. Use only IDs from the input.
- `weight` ∈ [0, 1]. Use 1.0 for clearly core members, 0.6–0.9 for supporting members, lower for loose associations.
- `role` ∈ `{"core", "doc", "test", "sample", "related"}`.
- `confidence` ∈ [0, 1]. Use higher when multiple naming/tag/menu signals agree, lower when you're inferring.
- Be **conservative** about cross-feature relations. Emit `feature_relations` only when it's clearly subfeature/sibling/related.
- It is OK to have many features — typical plugin has 5–25 features. Don't lump unrelated functions together.
- Anything genuinely uncategorized goes in `unassigned_member_ids`. Be honest about uncertainty.
- Aim for features that a developer (or a user) would actually search for: "Substructure search", "MPO profiles", "Reaction enumeration", "Activity cliffs", "Scaffold tree", not "Chem widget #3".

Return ONLY the JSON object. No markdown fences, no commentary.
"""


def write_prompt() -> None:
    INSTRUCTIONS_PATH.parent.mkdir(parents=True, exist_ok=True)
    INSTRUCTIONS_PATH.write_text(PROMPT_BODY, encoding="utf-8")


# ---------------------------------------------------------------------------
# apply(): validate agent outputs, emit JSONL
# ---------------------------------------------------------------------------

def _validate_output(obj: dict, valid_ids: set[str]) -> tuple[list[dict], list[str]]:
    """Return (valid_features, problems)."""
    problems: list[str] = []
    if not isinstance(obj, dict):
        problems.append("output is not a dict")
        return [], problems
    feats = obj.get("features") or []
    if not isinstance(feats, list):
        problems.append("'features' is not a list")
        return [], problems
    valid: list[dict] = []
    for i, f in enumerate(feats):
        fid = f.get("id", "")
        if not isinstance(f, dict) or not fid.startswith("feature:"):
            problems.append(f"feature[{i}].id missing or wrong prefix: {fid!r}")
            continue
        members = f.get("members") or []
        if not isinstance(members, list) or not members:
            problems.append(f"feature {fid} has no members")
            continue
        f_clean = dict(f)
        f_clean["members"] = []
        for m in members:
            mid = m.get("id") if isinstance(m, dict) else None
            if not mid or mid not in valid_ids:
                problems.append(f"feature {fid} unknown member id: {mid}")
                continue
            f_clean["members"].append(m)
        if not f_clean["members"]:
            problems.append(f"feature {fid} has zero valid members after filtering")
            continue
        valid.append(f_clean)
    return valid, problems


def apply(packages: list[str] | None = None) -> dict[str, int]:
    """Read agent outputs from OUTPUTS_DIR, write Feature/HasFeature/PartOfFeature
    JSONL into the canonical store, return counts."""
    nodes_dir = DATA_DIR / "nodes"
    edges_dir = DATA_DIR / "edges"
    nodes_dir.mkdir(parents=True, exist_ok=True)
    edges_dir.mkdir(parents=True, exist_ok=True)

    # Build a known-id set from existing JSONL so we validate references.
    valid_ids: set[str] = set()
    for jl in nodes_dir.glob("*.jsonl"):
        for line in jl.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                valid_ids.add(json.loads(line)["id"])
            except Exception:
                continue

    feature_rows: list[dict] = []
    has_feature_rows: list[dict] = []
    part_of_rows: list[dict] = []
    feat_ids: set[str] = set()

    counts = {"packages_applied": 0, "features": 0, "members": 0, "problems": 0}
    for out_path in sorted(OUTPUTS_DIR.glob("*.json")):
        pkg_folder = out_path.stem
        if packages and pkg_folder not in packages:
            continue
        try:
            obj = json.loads(out_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            print(f"  ! {out_path.name}: bad JSON ({e})")
            counts["problems"] += 1
            continue

        valid_feats, problems = _validate_output(obj, valid_ids)
        if problems:
            for p in problems[:8]:
                print(f"  ~ {pkg_folder}: {p}")
            counts["problems"] += len(problems)
        if not valid_feats:
            continue

        pkg_id = obj.get("package_id") or c.pkg_id(pkg_folder)
        for f in valid_feats:
            fid = f["id"]
            if fid in feat_ids:
                # rename per-package collision
                fid = f"{fid}__{pkg_folder}"
                f["id"] = fid
            feat_ids.add(fid)
            feature = Feature(
                id=fid,
                name=f.get("name") or fid,
                source_layer=SourceLayer.SYNTHETIC,
                cluster_score=float(f.get("confidence") or 0.0),
                member_ids=[m["id"] for m in f["members"]],
                description=f.get("description"),
                description_provenance=Provenance.LLM,
                extras={"category": f.get("category"),
                        "evidence_summary": f.get("evidence_summary"),
                        "package_id": pkg_id},
                extracted_by=NAME,
            )
            feature_rows.append(feature.model_dump(mode="json"))

            has_feature_rows.append(HasFeature(
                from_id=pkg_id, to_id=fid,
                derived_by=Provenance.LLM,
                confidence=float(f.get("confidence") or 0.0),
                extracted_by=NAME,
            ).model_dump(mode="json"))

            for m in f["members"]:
                part_of_rows.append(PartOfFeature(
                    from_id=m["id"], to_id=fid,
                    derived_by=Provenance.LLM,
                    confidence=float(m.get("weight") or 1.0),
                    weight=float(m.get("weight") or 1.0),
                    role=m.get("role"),
                    extracted_by=NAME,
                ).model_dump(mode="json"))
            counts["members"] += len(f["members"])
            counts["features"] += 1
        counts["packages_applied"] += 1

    # Merge into canonical JSONL slices (dedupe by id / tuple)
    _merge_into(nodes_dir / "Feature.jsonl", feature_rows, "id")
    _merge_into(edges_dir / "HAS_FEATURE.jsonl", has_feature_rows, ("from_id", "to_id"))
    _merge_into(edges_dir / "PART_OF_FEATURE.jsonl", part_of_rows, ("from_id", "to_id"))

    print(f"[{NAME}] applied {counts['packages_applied']} packages -> "
          f"{counts['features']} features, {counts['members']} memberships, "
          f"{counts['problems']} problems")
    return counts


def _merge_into(path: Path, new_rows: list[dict], key: str | tuple[str, ...]) -> None:
    existing: dict[Any, str] = {}
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            obj = json.loads(line)
            k = obj[key] if isinstance(key, str) else tuple(obj[x] for x in key)
            existing[k] = line
    for r in new_rows:
        k = r[key] if isinstance(key, str) else tuple(r[x] for x in key)
        existing[k] = json.dumps(r, separators=(",", ":"), ensure_ascii=False)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for line in existing.values():
            f.write(line + "\n")
