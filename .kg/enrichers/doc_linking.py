"""LLM enricher: link help/ DocPages to Packages, Features, and Functions.

The deterministic `docs` extractor only finds DocPages that code explicitly
points at via `//help-url:` annotations or `this.helpUrl = ...` assignments
— that's ~200 edges total. But `help/` has 695 pages, most of which
describe specific package features yet have no machine-readable back-pointer.
This enricher fills the gap: an agent reads each subtree of help/ and emits
`Documents` edges to the relevant code entities.

Three-step contract (same as `feature_clustering`):
  1. prepare()  — slice DocPages by help/ subtree, write input JSON per slice
  2. (Claude Code agent) — read input + PROMPT.md, emit output JSON
  3. apply()    — validate, dedupe vs deterministic edges, emit Documents JSONL

Slicing strategy: one slice per second-level path under `/help/`
(e.g. `/help/visualize`, `/help/datagrok`, `/help/explore`). The agent for
each slice sees ALL packages + ALL features as candidate targets, so
cross-domain docs (e.g. a Compute doc that mentions Chem) can still link.

Output edges always carry `derivation: "llm_doc_linking"` and `derived_by:
LLM` so they're distinguishable from deterministic ones.
"""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from schema.base import Provenance, SourceLayer
from schema.relations import Documents, PartOfFeature

from extractors import _common as c

NAME = "doc_linking"

DATA_DIR = ROOT / "data"
ENRICH_DIR = DATA_DIR / "enrichment" / NAME
INPUTS_DIR = ENRICH_DIR / "inputs"
OUTPUTS_DIR = ENRICH_DIR / "outputs"
PROMPT_PATH = ENRICH_DIR / "PROMPT.md"


# ---------------------------------------------------------------------------
# Loaders (same shape as feature_clustering)
# ---------------------------------------------------------------------------

def _load_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _load_all_nodes() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for f in sorted((DATA_DIR / "nodes").glob("*.jsonl")):
        out[f.stem] = _load_jsonl(f)
    return out


def _load_edges(predicate: str) -> list[dict]:
    return _load_jsonl(DATA_DIR / "edges" / f"{predicate}.jsonl")


def _truncate(s: str | None, n: int) -> str | None:
    if not s:
        return s
    s = s.strip().replace("\n\n", " · ").replace("\n", " ")
    return s if len(s) <= n else s[:n].rsplit(" ", 1)[0] + "…"


def _read_doc_snippet(doc: dict, char_budget: int = 350) -> str:
    """Read the first paragraph or two of the source markdown for context."""
    paths = doc.get("paths") or []
    if not paths:
        return ""
    src = paths[0].get("path")
    if not src:
        return ""
    f = c.REPO_ROOT / src
    try:
        text = f.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    # Strip frontmatter and HTML title comments
    text = re.sub(r"^---\n.*?\n---\n", "", text, count=1, flags=re.DOTALL)
    text = re.sub(r"<!--[^>]*-->", "", text)
    # Drop heading lines so the snippet is prose
    lines = [l for l in text.splitlines() if not l.lstrip().startswith("#")]
    body = " ".join(l.strip() for l in lines if l.strip())
    return _truncate(body, char_budget) or ""


# ---------------------------------------------------------------------------
# Slicing
# ---------------------------------------------------------------------------

def _slice_for(doc: dict) -> str | None:
    """Return the slice id for a DocPage, or None to skip."""
    url = (doc.get("extras") or {}).get("url", "") or ""
    if not url.startswith("/help/"):
        return None
    parts = [p for p in url.split("/") if p]   # ['help', 'visualize', 'viewers', 'scatter-plot']
    if len(parts) < 2:
        return None
    sub = parts[1]
    if sub.startswith("_") or sub == "uploads":
        return None
    return f"help-{sub}"


def _scrub_function(f: dict) -> dict | None:
    """Compact RegisteredFunction for the candidates list. Skip non-public roles."""
    role = f.get("role")
    if not role:
        return None
    if role in ("init", "autostart"):
        return None
    return {
        "id": f["id"],
        "name": f.get("name"),
        "friendly_name": f.get("friendly_name"),
        "role": role,
        "package_id": f.get("package_id"),
        "tags": f.get("tags") or [],
        "description": _truncate(f.get("description"), 140),
    }


def _scrub_package(p: dict) -> dict:
    e = p.get("extras") or {}
    return {
        "id": p["id"], "name": p["name"],
        "category": p.get("category"),
        "friendly_name": p.get("friendly_name"),
        "description": _truncate(p.get("description"), 200),
        "browserFeatures": e.get("browserFeatures"),
    }


def _scrub_feature(f: dict, by_pkg: dict[str, dict]) -> dict:
    e = f.get("extras") or {}
    pkg_id = e.get("package_id")
    return {
        "id": f["id"],
        "name": f.get("name"),
        "package": pkg_id,
        "package_name": (by_pkg.get(pkg_id) or {}).get("name"),
        "category": e.get("category"),
        "description": _truncate(f.get("description"), 200),
    }


def prepare(slices_filter: list[str] | None = None) -> list[Path]:
    """Build per-slice input JSON. Returns list of input paths written."""
    INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

    nodes = _load_all_nodes()
    docs = nodes.get("DocPage", [])
    pkgs = nodes.get("Package", [])
    feats = nodes.get("Feature", [])
    funcs = nodes.get("RegisteredFunction", [])

    by_pkg = {p["id"]: p for p in pkgs}

    # Group docs by slice
    docs_by_slice: dict[str, list[dict]] = defaultdict(list)
    for d in docs:
        sl = _slice_for(d)
        if sl is None:
            continue
        docs_by_slice[sl].append(d)

    # Compact candidate sets — same in every slice (cross-domain links matter)
    pkg_cands = [_scrub_package(p) for p in pkgs]
    feat_cands = [_scrub_feature(f, by_pkg) for f in feats]
    func_cands_all = [s for s in (_scrub_function(f) for f in funcs) if s]

    written: list[Path] = []
    for sl, ds in sorted(docs_by_slice.items()):
        if slices_filter and sl not in slices_filter:
            continue
        # Gather doc snippets
        slim_docs = []
        for d in sorted(ds, key=lambda x: ((x.get("extras") or {}).get("url") or "")):
            extras = d.get("extras") or {}
            slim_docs.append({
                "id": d["id"],
                "url": extras.get("url"),
                "title": d.get("title") or d.get("name"),
                "audience": d.get("audience"),
                "doc_kind": d.get("doc_kind"),
                "snippet": _read_doc_snippet(d),
            })

        # Per-slice function candidate trim: keep all role-tagged user-facing
        # ones; cap at 600 to stay under context budget.
        func_cands = func_cands_all[:600]

        ctx = {
            "slice_id": sl,
            "doc_count": len(slim_docs),
            "docs": slim_docs,
            "candidate_packages": pkg_cands,
            "candidate_features": feat_cands,
            "candidate_functions": func_cands,
            "stats": {
                "docs": len(slim_docs),
                "packages": len(pkg_cands),
                "features": len(feat_cands),
                "functions": len(func_cands),
            },
        }
        out = INPUTS_DIR / f"{sl}.json"
        out.write_text(json.dumps(ctx, indent=2, ensure_ascii=False), encoding="utf-8")
        written.append(out)

    PROMPT_PATH.write_text(PROMPT_BODY, encoding="utf-8")
    print(f"[{NAME}] wrote {len(written)} slice inputs to {INPUTS_DIR}")
    for p in written:
        size_kb = p.stat().st_size / 1024
        print(f"  {p.name:<32} {size_kb:>7.0f} KB")
    return written


# ---------------------------------------------------------------------------
# Prompt
# ---------------------------------------------------------------------------

PROMPT_BODY = """\
# Doc → code linking for one help/ subtree

You are linking documentation pages in `help/` to the code entities
they describe. Your output drives `Documents` edges in the Datagrok
knowledge graph.

## Input

A JSON object with:
- `slice_id` — short id for this batch (e.g. `help-visualize`).
- `docs` — DocPage entries (`id`, `url`, `title`, `snippet` of body).
- `candidate_packages` — every Package in the corpus.
- `candidate_features` — every Feature (LLM-clustered groups of functions/scripts/docs that already hang together).
- `candidate_functions` — user-facing role-tagged Functions (role in {app, viewer, widget, panel, fileViewer, dashboard, semTypeDetector, transform, editor, …}).

Use IDs verbatim. Do NOT invent ids.

## Your output (STRICT JSON, no markdown fences, no commentary)

```json
{
  "slice_id": "help-visualize",
  "links": [
    {
      "doc_id": "doc:/help/visualize/viewers/scatter-plot",
      "target_id": "func:Charts:Sankey",         // or feature:..., or pkg:...
      "confidence": 0.92,
      "evidence": "title and content describe Sankey diagrams"
    }
  ]
}
```

## Rules

1. **Prefer Feature → DocPage links over Function → DocPage** when a feature
   exists. Features already group their member functions, so the indirect
   chain `Function ←PART_OF_FEATURE— Feature ←DOCUMENTS— DocPage` is more
   meaningful than spamming direct edges.

2. **Emit Function-level links only when the doc is unambiguously about a
   single function** (e.g. `/help/visualize/viewers/scatter-plot` → the
   scatter-plot viewer function specifically).

3. **Package-level links are useful for** package overviews, READMEs, or
   getting-started pages that cover many features at once
   (e.g. `/help/datagrok/solutions/domains/chem/overview` → `pkg:Chem`).

4. **A doc may link to multiple targets.** Common case: a "Cheminformatics
   solutions" doc links to `pkg:Chem` AND several `feature:Chem:*` features.
   Emit one row per target.

5. **Confidence**:
   - `0.95+` — title or URL slug exactly matches a Feature/Function name.
   - `0.80-0.94` — strong textual match (snippet describes the entity unambiguously).
   - `0.60-0.79` — domain match without specific naming alignment.
   - Below `0.60` — don't emit. Skip.

6. **Cross-domain links are fine.** A "Add new column" doc in `/help/transform/`
   may legitimately reference a Power Pack feature. Cross-package linking is
   the whole point of this enricher.

7. **One doc → many targets is normal.** Don't try to make it 1:1.

8. **Be exhaustive but not noisy.** Aim to link every doc to ≥ 1 target if
   any plausible match exists, but don't reach. Better to skip an unclear
   doc than to invent a wrong link.

9. **Reserve viewer-class slug matching** (e.g. `scatter-plot.md` → scatter-plot
   viewer). Datagrok core viewers (Bar Chart, Scatter Plot, Box Plot, …) are
   built into the platform itself — they may NOT have a corresponding
   RegisteredFunction in the input. In that case, link to the most relevant
   Charts/PowerGrid feature instead, with a slightly lower confidence (0.7-0.8).

Return ONLY the JSON. No prose around it.
"""


# ---------------------------------------------------------------------------
# Apply
# ---------------------------------------------------------------------------

def apply(slices_filter: list[str] | None = None) -> dict[str, int]:
    """Read agent outputs, emit Documents edges to the canonical store."""
    nodes_dir = DATA_DIR / "nodes"
    edges_dir = DATA_DIR / "edges"
    edges_dir.mkdir(parents=True, exist_ok=True)

    # Build set of valid IDs
    valid_ids: set[str] = set()
    for jl in nodes_dir.glob("*.jsonl"):
        for line in jl.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                valid_ids.add(json.loads(line)["id"])
            except Exception:
                continue

    # Load existing Documents edges so we can preserve deterministic provenance
    existing: dict[tuple[str, str], dict] = {}
    edges_path = edges_dir / "DOCUMENTS.jsonl"
    if edges_path.is_file():
        for line in edges_path.read_text(encoding="utf-8").splitlines():
            if not line.strip(): continue
            row = json.loads(line)
            existing[(row["from_id"], row["to_id"])] = row

    counts = {"slices_applied": 0, "links_added": 0, "skipped_dup": 0,
              "skipped_invalid": 0, "skipped_low_conf": 0,
              "part_of_feature_added": 0}

    # Existing PartOfFeature edges so we don't dupe membership
    pof_path = edges_dir / "PART_OF_FEATURE.jsonl"
    pof_existing: dict[tuple[str, str], dict] = {}
    if pof_path.is_file():
        for line in pof_path.read_text(encoding="utf-8").splitlines():
            if not line.strip(): continue
            row = json.loads(line)
            pof_existing[(row["from_id"], row["to_id"])] = row
    new_pof_rows: list[dict] = []

    new_rows: list[dict] = []
    for out_path in sorted(OUTPUTS_DIR.glob("*.json")):
        sl = out_path.stem
        if slices_filter and sl not in slices_filter:
            continue
        try:
            obj = json.loads(out_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            print(f"  ! {out_path.name}: {e}")
            continue
        for ln in obj.get("links") or []:
            d = ln.get("doc_id"); t = ln.get("target_id")
            if not d or not t:
                counts["skipped_invalid"] += 1; continue
            if d not in valid_ids or t not in valid_ids:
                counts["skipped_invalid"] += 1; continue
            try:
                conf = float(ln.get("confidence", 0))
            except (TypeError, ValueError):
                conf = 0.0
            if conf < 0.6:
                counts["skipped_low_conf"] += 1; continue
            key = (d, t)
            prior = existing.get(key)
            should_emit_doc = True
            if prior:
                prior_der = prior.get("derived_by") or "llm"
                prior_conf = float(prior.get("confidence") or 0)
                if prior_der != Provenance.LLM.value:
                    should_emit_doc = False
                elif prior_conf >= conf:
                    should_emit_doc = False
            if should_emit_doc:
                edge = Documents(
                    from_id=d, to_id=t,
                    derivation="llm_doc_linking",
                    derived_by=Provenance.LLM,
                    confidence=conf,
                    extras={"slice": sl, "evidence": ln.get("evidence")},
                    extracted_by=NAME,
                )
                new_rows.append(edge.model_dump(mode="json"))
                counts["links_added"] += 1
            else:
                counts["skipped_dup"] += 1

            # Seed a PART_OF_FEATURE membership edge if the target is a Feature
            # — independent of whether the DOCUMENTS edge is new or a dup,
            # because the Feature membership may not yet exist even when the
            # Documents edge does (e.g. on a re-apply pass).
            if t.startswith("feature:") and (d, t) not in pof_existing:
                pof = PartOfFeature(
                    from_id=d, to_id=t,
                    weight=conf,
                    role="doc",
                    derived_by=Provenance.LLM,
                    confidence=conf,
                    extras={"source_enricher": NAME},
                    extracted_by=NAME,
                ).model_dump(mode="json")
                new_pof_rows.append(pof)
                pof_existing[(d, t)] = pof   # avoid double-add within one run
                counts["part_of_feature_added"] += 1
        counts["slices_applied"] += 1

    # Merge into the existing DOCUMENTS JSONL by (from_id, to_id)
    for r in new_rows:
        existing[(r["from_id"], r["to_id"])] = r
    with edges_path.open("w", encoding="utf-8", newline="\n") as f:
        for line in existing.values():
            if isinstance(line, dict):
                f.write(json.dumps(line, separators=(",", ":"), ensure_ascii=False) + "\n")
            else:
                f.write(line if line.endswith("\n") else line + "\n")

    # Merge new PART_OF_FEATURE seeds (only when key is new)
    for r in new_pof_rows:
        pof_existing[(r["from_id"], r["to_id"])] = r
    with pof_path.open("w", encoding="utf-8", newline="\n") as f:
        for row in pof_existing.values():
            f.write(json.dumps(row, separators=(",", ":"), ensure_ascii=False) + "\n")

    print(f"[{NAME}] applied {counts['slices_applied']} slices -> "
          f"{counts['links_added']} DOCUMENTS edges, "
          f"{counts['part_of_feature_added']} PART_OF_FEATURE edges "
          f"(dups: {counts['skipped_dup']}, "
          f"invalid: {counts['skipped_invalid']}, "
          f"low-conf: {counts['skipped_low_conf']})")
    return counts
