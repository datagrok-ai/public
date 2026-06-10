"""LLM enricher for IS_IMPLEMENTED_IN / IS_TESTED_IN refinement.

The deterministic extractor (`extractors.feature_implementation`) handles
the obvious case: parse `package.ts` decorator wrappers, walk imports,
attribute files. But many features have implementations that the import
chain doesn't reach — e.g. cell renderers in `viewers/`, JSON config in
`files/`, Python scripts that get registered separately, helper modules
called only via `grok.functions.eval`. Tests are similarly under-covered
because not every test file mentions the function name verbatim.

This enricher dispatches one Claude Code agent per package. The agent
reads the package source tree and assigns implementation/test files to
each Feature.

Three-step contract (same as the other enrichers):
  1. prepare()  -> data/enrichment/feature_implementation/inputs/<Pkg>.json
  2. (Claude Code agents)  -> outputs/<Pkg>.json
  3. apply()    -> emit IS_IMPLEMENTED_IN / IS_TESTED_IN edges (derived_by=LLM)
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from schema.base import Provenance, SourceLayer
from schema.relations import IsImplementedIn, IsTestedIn

from extractors import _common as c

NAME = "feature_implementation"

DATA_DIR = c.REPO_ROOT / ".kg" / "data"
ENRICH_DIR = DATA_DIR / "enrichment" / NAME
INPUTS_DIR = ENRICH_DIR / "inputs"
OUTPUTS_DIR = ENRICH_DIR / "outputs"
INSTRUCTIONS_PATH = ENRICH_DIR / "PROMPT.md"

# Cap context so input JSONs stay reasonable for the agent prompt.
_MAX_FILES_LISTED = 600
_MAX_FUNCTIONS_LISTED = 400


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def _load_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()]


# ---------------------------------------------------------------------------
# prepare()
# ---------------------------------------------------------------------------

def _truncate(s: str | None, n: int) -> str | None:
    if not s:
        return s
    s = s.strip().replace("\r", "")
    return s if len(s) <= n else s[:n].rsplit(" ", 1)[0] + "…"


def _pkg_src_files(pkg_dir: Path) -> dict[str, list[str]]:
    """Walk the package tree and bucket files by purpose. Returns dict with
    keys: src, tests, scripts, queries, docs, files. Paths are repo-relative
    POSIX, sorted, capped per bucket."""
    buckets: dict[str, list[str]] = {
        "src":     [],
        "tests":   [],
        "scripts": [],
        "queries": [],
        "docs":    [],
        "config":  [],
    }
    skip_exts = {".png", ".jpg", ".gif", ".ico", ".svg", ".webp",
                 ".wasm", ".woff", ".woff2"}
    for f in pkg_dir.rglob("*"):
        if not f.is_file():
            continue
        if any(p in {"node_modules", "dist", "build", ".cache", ".git"}
               for p in f.parts):
            continue
        if f.suffix.lower() in skip_exts:
            continue
        rel = c.repo_rel(f)
        if "/src/tests/" in rel or rel.endswith("package-test.ts"):
            buckets["tests"].append(rel)
        elif rel.endswith(".g.ts") or rel.endswith("package-api.ts"):
            continue                                # generated — excluded
        elif "/src/" in rel:
            buckets["src"].append(rel)
        elif "/scripts/" in rel:
            buckets["scripts"].append(rel)
        elif "/queries/" in rel:
            buckets["queries"].append(rel)
        elif rel.endswith(".md"):
            buckets["docs"].append(rel)
        elif "/connections/" in rel or "/environments/" in rel:
            buckets["config"].append(rel)
    for k in buckets:
        buckets[k] = sorted(buckets[k])[:_MAX_FILES_LISTED]
    return buckets


def _features_for_package(pkg_id: str,
                          features: list[dict],
                          part_of: list[dict]) -> list[dict]:
    """Return features whose package_id (in extras or via HAS_FEATURE) match."""
    feat_index = {f["id"]: f for f in features}
    member_to_feat: dict[str, set[str]] = defaultdict(set)
    for row in part_of:
        member_to_feat[row["from_id"]].add(row["to_id"])
    out: list[dict] = []
    for f in features:
        extras = f.get("extras") or {}
        if extras.get("package_id") == pkg_id:
            out.append(f)
    return out


def _ast_edges_for_features(feat_ids: set[str], edges: list[dict]) -> dict[str, list[str]]:
    """Return {feature_id: [file_paths_already_attributed_by_AST]}.
    Helps the agent build on top of the deterministic baseline."""
    by_feat: dict[str, set[str]] = defaultdict(set)
    for e in edges:
        if e.get("from_id") in feat_ids and e.get("derived_by") == "ast":
            tid = e.get("to_id", "")
            if tid.startswith("file:"):
                by_feat[e["from_id"]].add(tid.removeprefix("file:"))
    return {k: sorted(v) for k, v in by_feat.items()}


def prepare(packages: list[str] | None = None) -> list[Path]:
    """Write per-package input JSONs."""
    INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    write_prompt()

    nodes_dir = DATA_DIR / "nodes"
    edges_dir = DATA_DIR / "edges"

    pkgs = _load_jsonl(nodes_dir / "Package.jsonl")
    features = _load_jsonl(nodes_dir / "Feature.jsonl")
    rfs = _load_jsonl(nodes_dir / "RegisteredFunction.jsonl")
    scripts = _load_jsonl(nodes_dir / "Script.jsonl")
    queries = _load_jsonl(nodes_dir / "DataQuery.jsonl")
    docs = _load_jsonl(nodes_dir / "DocPage.jsonl")
    pof = _load_jsonl(edges_dir / "PART_OF_FEATURE.jsonl")
    impl_edges = _load_jsonl(edges_dir / "IS_IMPLEMENTED_IN.jsonl")
    test_edges = _load_jsonl(edges_dir / "IS_TESTED_IN.jsonl")

    rf_by_id = {rf["id"]: rf for rf in rfs}
    script_by_id = {s["id"]: s for s in scripts}
    query_by_id = {q["id"]: q for q in queries}
    doc_by_id = {d["id"]: d for d in docs}

    # Map feature -> [member entity descriptions]
    members_by_feat: dict[str, list[dict]] = defaultdict(list)
    for row in pof:
        feat = row["to_id"]
        mid = row["from_id"]
        if mid in rf_by_id:
            rf = rf_by_id[mid]
            members_by_feat[feat].append({
                "id": mid, "kind": "RegisteredFunction",
                "name": rf.get("name"),
                "friendly_name": rf.get("friendly_name"),
                "role": rf.get("role"),
                "description": _truncate(rf.get("description"), 160),
                "tags": rf.get("tags") or [],
            })
        elif mid in script_by_id:
            s = script_by_id[mid]
            members_by_feat[feat].append({
                "id": mid, "kind": "Script",
                "name": s.get("name"), "language": s.get("language"),
                "description": _truncate(s.get("description"), 160),
            })
        elif mid in query_by_id:
            q = query_by_id[mid]
            members_by_feat[feat].append({
                "id": mid, "kind": "DataQuery",
                "name": q.get("name"),
                "description": _truncate(q.get("description"), 160),
            })
        elif mid in doc_by_id:
            d = doc_by_id[mid]
            members_by_feat[feat].append({
                "id": mid, "kind": "DocPage",
                "title": d.get("title"),
                "url": (d.get("extras") or {}).get("url"),
            })
        else:
            members_by_feat[feat].append({"id": mid, "kind": "other"})

    pkgs_root = c.REPO_ROOT / "packages"
    written: list[Path] = []
    for pkg in pkgs:
        folder = pkg["name"]
        if packages and folder not in packages:
            continue
        pkg_dir = pkgs_root / folder
        if not pkg_dir.is_dir():
            continue

        pkg_features = _features_for_package(pkg["id"], features, pof)
        if not pkg_features:
            continue

        feat_ids = {f["id"] for f in pkg_features}
        ast_impl = _ast_edges_for_features(feat_ids, impl_edges)
        ast_test = _ast_edges_for_features(feat_ids, test_edges)

        ctx = {
            "package": {
                "id": pkg["id"],
                "name": folder,
                "friendly_name": pkg.get("friendly_name"),
                "category": pkg.get("category"),
                "description": _truncate(pkg.get("description"), 400),
            },
            "files": _pkg_src_files(pkg_dir),
            "features": [
                {
                    "id": f["id"],
                    "name": f.get("name"),
                    "description": _truncate(f.get("description"), 220),
                    "category": (f.get("extras") or {}).get("category"),
                    "members": members_by_feat.get(f["id"], [])[:_MAX_FUNCTIONS_LISTED],
                    "ast_baseline": {
                        "implements": ast_impl.get(f["id"], []),
                        "tests":      ast_test.get(f["id"], []),
                    },
                }
                for f in pkg_features
            ],
        }
        out_path = INPUTS_DIR / f"{folder}.json"
        out_path.write_text(json.dumps(ctx, indent=2, ensure_ascii=False),
                            encoding="utf-8")
        written.append(out_path)

    print(f"[{NAME}] wrote {len(written)} input contexts to {INPUTS_DIR}")
    return written


# ---------------------------------------------------------------------------
# Prompt
# ---------------------------------------------------------------------------

PROMPT_BODY = """\
# Feature implementation/test mapping

You are mapping each Feature in a single Datagrok plugin to the source
files that **implement** it and the source files that **test** it.

## Why this exists

Datagrok plugins follow a wrapper pattern:
- `src/package.g.ts` is auto-generated and just contains `//name:` annotation
  blocks + thin `export async function foo() { return PackageFunctions.foo() }`
  wrappers.
- `src/package.ts` declares `class PackageFunctions` whose decorated `static`
  methods call into the real implementation in `src/widgets/`, `src/utils/`,
  `src/viewers/`, `src/analysis/` etc.

A naive walker only sees `package.g.ts`. We have a deterministic AST walker
that follows the `package.ts` import chain — that result is in each feature's
`ast_baseline.implements`. Your job is to **extend and correct** that baseline.

## Your input

A JSON object:

```
{
  "package": { "id", "name", "friendly_name", "category", "description" },
  "files": {
    "src":     [paths under <pkg>/src/ excluding tests/, .g.ts, package-api.ts],
    "tests":   [paths under <pkg>/src/tests/ + package-test.ts],
    "scripts": [paths under <pkg>/scripts/],
    "queries": [paths under <pkg>/queries/],
    "docs":    [.md files in pkg],
    "config":  [paths under connections/ + environments/]
  },
  "features": [
    {
      "id":          "feature:<Pkg>:<slug>",
      "name":        "Human-readable feature name",
      "description": "What it does for the user",
      "members":     [ {kind, id, name, friendly_name, role, ...}, ... ],
      "ast_baseline": {
        "implements": [...repo-relative paths...],
        "tests":      [...repo-relative paths...]
      }
    }
  ]
}
```

## What to do

For each feature, decide which **files** in the package implement it and
which **test files** test it. Use the `members` list (it tells you the
RegisteredFunctions, Scripts, etc. that comprise the feature) plus the
file listings to find the real implementation. You may use `Read`, `Glob`,
and `Grep` tools to inspect any file under the package directory.

## Rules

- Only emit paths that appear in the `files` listing — exact match.
- Skip `package.g.ts` and any `*-api.ts` / `*.api.g.ts` (auto-generated).
- `package.ts` IS valid — it's the entry point that wires everything up.
- A feature typically has 2–10 implementation files and 0–3 test files.
- Be conservative: if a file is referenced but is just a one-liner type
  re-export, skip it. If you're not sure, skip.
- Tests must be tests of *this* feature (file mentions the feature's name,
  one of its functions, a related class, etc.) — not just any test in the
  package.
- For demo files (`src/demo*/`, `src/apps/`-like), if they specifically
  demo this feature, include them in `implements` with `role: "demo"`.

## Your output (STRICT JSON, no other text)

```json
{
  "package_id": "pkg:<Pkg>",
  "features": [
    {
      "id": "feature:<Pkg>:<slug>",
      "implements": [
        {"path": "packages/<Pkg>/src/widgets/foo.ts", "role": "core",  "confidence": 0.95},
        {"path": "packages/<Pkg>/src/utils/bar.ts",   "role": "support", "confidence": 0.85}
      ],
      "tests": [
        {"path": "packages/<Pkg>/src/tests/foo-tests.ts", "confidence": 0.9}
      ],
      "notes": "1-sentence justification, optional"
    }
  ]
}
```

`role` ∈ `{"core", "support", "demo", "config"}`. `confidence` ∈ [0, 1].

Return ONLY the JSON object. No markdown fences, no commentary.
"""


def write_prompt() -> None:
    INSTRUCTIONS_PATH.parent.mkdir(parents=True, exist_ok=True)
    INSTRUCTIONS_PATH.write_text(PROMPT_BODY, encoding="utf-8")


# ---------------------------------------------------------------------------
# apply()
# ---------------------------------------------------------------------------

def _validate(obj: dict, valid_files: set[str], valid_feats: set[str]
              ) -> tuple[list[dict], list[str]]:
    problems: list[str] = []
    if not isinstance(obj, dict):
        return [], ["output is not a dict"]
    feats = obj.get("features") or []
    if not isinstance(feats, list):
        return [], ["'features' is not a list"]
    out: list[dict] = []
    for f in feats:
        fid = f.get("id", "")
        if fid not in valid_feats:
            problems.append(f"unknown feature id: {fid}")
            continue
        impl = [x for x in (f.get("implements") or [])
                if isinstance(x, dict) and x.get("path") in valid_files]
        tests = [x for x in (f.get("tests") or [])
                 if isinstance(x, dict) and x.get("path") in valid_files]
        out.append({"id": fid, "implements": impl, "tests": tests})
    return out, problems


def apply(packages: list[str] | None = None) -> dict[str, int]:
    nodes_dir = DATA_DIR / "nodes"
    edges_dir = DATA_DIR / "edges"
    edges_dir.mkdir(parents=True, exist_ok=True)

    valid_files: set[str] = set()
    for f in _load_jsonl(nodes_dir / "File.jsonl"):
        rp = f.get("relative_path")
        if rp:
            valid_files.add(rp)
    valid_feats: set[str] = {f["id"] for f in _load_jsonl(nodes_dir / "Feature.jsonl")}

    impl_rows: list[dict] = []
    test_rows: list[dict] = []
    counts = {"packages_applied": 0, "impl_edges": 0, "test_edges": 0,
              "problems": 0, "skipped_files": 0}

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

        feats, problems = _validate(obj, valid_files, valid_feats)
        if problems:
            for p in problems[:6]:
                print(f"  ~ {pkg_folder}: {p}")
            counts["problems"] += len(problems)

        for f in feats:
            for x in f["implements"]:
                impl_rows.append(IsImplementedIn(
                    from_id=f["id"], to_id=f"file:{x['path']}",
                    member_count=1,
                    derived_by=Provenance.LLM,
                    confidence=float(x.get("confidence") or 0.7),
                    extracted_by=NAME,
                    extras={"role": x.get("role")},
                ).model_dump(mode="json"))
            for x in f["tests"]:
                test_rows.append(IsTestedIn(
                    from_id=f["id"], to_id=f"file:{x['path']}",
                    member_count=1,
                    derived_by=Provenance.LLM,
                    confidence=float(x.get("confidence") or 0.7),
                    extracted_by=NAME,
                ).model_dump(mode="json"))
        counts["packages_applied"] += 1
        counts["impl_edges"] += sum(len(f["implements"]) for f in feats)
        counts["test_edges"] += sum(len(f["tests"]) for f in feats)

    _merge_into(edges_dir / "IS_IMPLEMENTED_IN.jsonl", impl_rows,
                ("from_id", "to_id"))
    _merge_into(edges_dir / "IS_TESTED_IN.jsonl", test_rows,
                ("from_id", "to_id"))

    print(f"[{NAME}] applied {counts['packages_applied']} packages -> "
          f"{counts['impl_edges']} impl + {counts['test_edges']} test edges, "
          f"{counts['problems']} problems")
    return counts


def _merge_into(path: Path, new_rows: list[dict], key: tuple[str, ...]) -> None:
    existing: dict[tuple, str] = {}
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            obj = json.loads(line)
            k = tuple(obj[x] for x in key)
            existing[k] = line
    for r in new_rows:
        k = tuple(r[x] for x in key)
        existing[k] = json.dumps(r, separators=(",", ":"), ensure_ascii=False)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for line in existing.values():
            f.write(line + "\n")
