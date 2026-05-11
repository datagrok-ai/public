"""Build orchestrator: extractors -> JSONL -> Kuzu.

Usage:
    py build.py                          # run everything (extract + Kuzu)
    py build.py --extractors ts_plugin_package
    py build.py --packages Chem,Chembl   # restrict packages-layer extractors
    py build.py --no-db                  # JSONL only, skip Kuzu rebuild
    py build.py --clean                  # wipe JSONL slices first

The extractor registry (`EXTRACTORS`) maps short names to module paths.
Adding a new extractor: import it, register it, done.
"""

from __future__ import annotations

import argparse
import importlib
import json
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from extractors._common import Bundle, REPO_ROOT
from schema.ddl import all_ddl, _entity_columns, _edge_columns
from schema.entities import ENTITY_KINDS
from schema.relations import RELATION_KINDS


EXTRACTORS: dict[str, str] = {
    # name -> module path. Order matters — earlier ones populate IDs that
    # later ones may reference. Library runs before plugin so that
    # Package -> Library DependsOn edges have a target to bind to.
    "library":                 "extractors.library",
    "ts_plugin_package":       "extractors.ts_plugin_package",
    "ts_plugin_scripts":       "extractors.ts_plugin_scripts",
    "ts_plugin_queries":       "extractors.ts_plugin_queries",
    "ts_plugin_connections":   "extractors.ts_plugin_connections",
    "ts_plugin_environments":  "extractors.ts_plugin_environments",
    "ts_plugin_dockerfiles":   "extractors.ts_plugin_dockerfiles",
    "ts_plugin_changelog":     "extractors.ts_plugin_changelog",
    "semtype_detectors":       "extractors.semtype_detectors",  # detectors.js → semType detectors as RegisteredFunction(role=semTypeDetector)
    "docs":                    "extractors.docs",
    "cross_package":           "extractors.cross_package",
    "js_api":                  "extractors.js_api",         # canonical JS API surface (TsClass etc. with is_jsapi=True). MUST run before ts_code_structure (resolves DG.X parents).
    "ts_code_structure":       "extractors.ts_code_structure",  # TS classes/interfaces/enums/funcs in packages + libraries
    "files":                   "extractors.files",          # promote files to nodes; emits DefinedIn for ALL definitional kinds + DefinedByFunction bridge
    "package_tests":           "extractors.package_tests",  # category()/test() blocks → PackageTest nodes
    "ts_imports":              "extractors.ts_imports",     # File→File / File→LibraryModule / File→JsApiNamespace import edges (depends on files extractor)
    "tags":                    "extractors.tags",           # Tag entity + HasTag edges (depends on RegisteredFunction.jsonl)
    "semtypes":                "extractors.semtypes",       # extract SemanticType + Consumes/Produces/Detects + RequiresColumnTag (depends on RegisteredFunction.jsonl)
    "feature_apply":           "extractors.feature_apply",  # re-emit Feature/HasFeature/PartOfFeature from preserved LLM outputs in data/enrichment/
    "feature_files":           "extractors.feature_files",  # post-pass: derive shallow IsImpl/IsTest (kept as fallback signal)
    "feature_implementation":  "extractors.feature_implementation",  # deeper AST-based IsImpl/IsTest from package.ts wrapper chains
    "js_api_usage":            "extractors.js_api_usage",   # File / TsMethod / TsFunction / Package -> JsApi edges (multi-granularity since 2026-05)
    "feature_classify":        "extractors.feature_classify",  # derive is_user_facing + interaction_kind
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extractors", default="", help="Comma list; default = all")
    ap.add_argument("--packages", default="", help="Comma list of package folder names to restrict to")
    ap.add_argument("--no-db", action="store_true", help="Skip Kuzu rebuild")
    ap.add_argument("--clean", action="store_true", help="Wipe data/ and kg.kuzu/ first")
    args = ap.parse_args()

    data_dir = ROOT / "data"
    db_path = ROOT / "kg.kuzu"

    if args.clean:
        # PRESERVE data/enrichment/ — the LLM agent outputs are precious and
        # not easily reproducible. Only wipe nodes/, edges/, and the DB.
        for sub in ("nodes", "edges"):
            d = data_dir / sub
            if d.exists() and d.is_dir():
                shutil.rmtree(d)
        if db_path.exists():
            if db_path.is_dir():
                shutil.rmtree(db_path)
            else:
                db_path.unlink()
            for sib in db_path.parent.glob(db_path.name + ".*"):
                try:
                    sib.unlink()
                except OSError:
                    pass
        print(f"cleaned {data_dir}/nodes, {data_dir}/edges, {db_path} "
              "(preserved data/enrichment/)")

    # ---- Run extractors -----------------------------------------------------
    selected = [s for s in (args.extractors or "").split(",") if s] or list(EXTRACTORS)
    pkg_filter = [s for s in (args.packages or "").split(",") if s] or None

    total_counts: dict[str, int] = {}
    t0 = time.time()
    for name in selected:
        mod_path = EXTRACTORS.get(name)
        if not mod_path:
            print(f"unknown extractor: {name}", file=sys.stderr)
            return 2
        mod = importlib.import_module(mod_path)
        bundle = Bundle(data_dir, extractor_name=name)
        kwargs = {}
        if pkg_filter and "package_filter" in mod.run.__code__.co_varnames:
            kwargs["package_filter"] = pkg_filter
        mod.run(REPO_ROOT, bundle, **kwargs)
        counts = bundle.flush()
        for k, v in counts.items():
            total_counts[k] = total_counts.get(k, 0) + v
    print(f"\nExtraction took {time.time() - t0:.1f}s")
    print(f"JSONL slices: {sum(1 for _ in (data_dir / 'nodes').glob('*.jsonl'))} node files, "
          f"{sum(1 for _ in (data_dir / 'edges').glob('*.jsonl'))} edge files")
    print("Counts:")
    for k, v in sorted(total_counts.items()):
        print(f"  {k:<40} {v:>6}")

    # ---- Build Kuzu ---------------------------------------------------------
    if args.no_db:
        return 0

    import kuzu
    if db_path.exists():
        if db_path.is_dir():
            shutil.rmtree(db_path)
        else:
            db_path.unlink()
        # Also remove sibling WAL/lock files Kuzu may have written next to it
        for sib in db_path.parent.glob(db_path.name + ".*"):
            try:
                sib.unlink()
            except OSError:
                pass

    print(f"\nBuilding Kuzu DB at {db_path} ...")
    t0 = time.time()
    db = kuzu.Database(str(db_path))
    conn = kuzu.Connection(db)

    # 1. DDL
    for stmt in all_ddl():
        conn.execute(stmt)

    # 2. Load nodes first
    nodes_dir = data_dir / "nodes"
    edges_dir = data_dir / "edges"
    n_nodes = 0
    if nodes_dir.is_dir():
        for jl in sorted(nodes_dir.glob("*.jsonl")):
            kind = jl.stem
            lines = [json.loads(line) for line in jl.read_text(encoding="utf-8").splitlines() if line.strip()]
            n_nodes += _insert_nodes(conn, kind, lines)

    # 3. Build id->kind index, then load edges
    id_to_kind = _index_node_ids(conn, list(ENTITY_KINDS.keys()))

    # 3a. Resolve `pkg-npm:<scoped-name>` -> `pkg:<Folder>` based on each
    # Package node's extras.npm_name. Lets cross-package devDeps land.
    npm_to_folder = _build_npm_to_folder(conn)

    n_edges = 0
    if edges_dir.is_dir():
        for jl in sorted(edges_dir.glob("*.jsonl")):
            predicate = jl.stem
            lines = [json.loads(line) for line in jl.read_text(encoding="utf-8").splitlines() if line.strip()]
            # Rewrite pkg-npm: targets to pkg: targets where we can resolve them
            for row in lines:
                tgt = row.get("to_id", "")
                if isinstance(tgt, str) and tgt.startswith("pkg-npm:"):
                    folder = npm_to_folder.get(tgt.removeprefix("pkg-npm:"))
                    if folder:
                        row["to_id"] = f"pkg:{folder}"
            n_edges += _insert_edges(conn, predicate, lines, id_to_kind)

    print(f"Loaded {n_nodes} nodes + {n_edges} edges in {time.time() - t0:.1f}s")
    return 0


def _build_npm_to_folder(conn) -> dict[str, str]:
    """Walk every Package node and read its extras.npm_name, returning
    `{<scoped-tail-after-/>: folder_name}` so we can resolve `pkg-npm:`
    placeholders emitted before all packages were known."""
    out: dict[str, str] = {}
    try:
        res = conn.execute("MATCH (p:Package) RETURN p.name, p.extras;")
    except Exception:
        return out
    while res.has_next():
        folder, extras_json = res.get_next()
        if not extras_json:
            continue
        try:
            extras = json.loads(extras_json)
        except (TypeError, json.JSONDecodeError):
            continue
        npm = extras.get("npm_name")
        if isinstance(npm, str) and npm.startswith("@datagrok/"):
            out[npm.split("/", 1)[1]] = folder
    return out


def _index_node_ids(conn, kinds: list[str]) -> dict[str, str]:
    """Walk every loaded node table once and build {id -> kind} for edge lookup."""
    id_to_kind: dict[str, str] = {}
    for kind in kinds:
        try:
            res = conn.execute(f"MATCH (n:`{kind}`) RETURN n.id;")
        except Exception:
            continue
        while res.has_next():
            (node_id,) = res.get_next()
            if node_id:
                id_to_kind[node_id] = kind
    return id_to_kind


def _insert_nodes(conn, kind: str, rows: list[dict]) -> int:
    """Insert all rows for one node kind via parameterized CREATE.
    Returns the count actually inserted."""
    if not rows:
        return 0
    model = ENTITY_KINDS.get(kind)
    if model is None:
        print(f"  ! unknown node kind: {kind}")
        return 0
    cols = [c for c, _t, _pk in _entity_columns(model)]
    cypher = f"CREATE (n:`{kind}` {{{', '.join(f'`{c}`: ${c}' for c in cols)}}});"
    n = 0
    for row in rows:
        params = {c: _coerce(row.get(c)) for c in cols}
        try:
            conn.execute(cypher, parameters=params)
            n += 1
        except Exception as e:
            # Common cause: duplicate primary key. Skip silently.
            if "duplicate" not in str(e).lower():
                print(f"  ! insert {kind} {row.get('id')}: {e}")
    return n


def _insert_edges(conn, predicate: str, rows: list[dict],
                  id_to_kind: dict[str, str]) -> int:
    """Insert edges. Look up the node-kind for both endpoints via the
    pre-built id_to_kind index, so we can pick the correct labeled MATCH."""
    if not rows:
        return 0
    model = RELATION_KINDS.get(predicate)
    if model is None:
        print(f"  ! unknown predicate: {predicate}")
        return 0
    extra_cols = [c for c, _t in _edge_columns(model)]

    n = 0
    skipped = 0
    for row in rows:
        from_id = row.get("from_id")
        to_id = row.get("to_id")
        if not from_id or not to_id:
            skipped += 1
            continue
        src_kind = id_to_kind.get(from_id)
        dst_kind = id_to_kind.get(to_id)
        if not src_kind or not dst_kind:
            skipped += 1   # orphan; either dep on Library that hasn't been extracted yet, or cross-package id
            continue
        if src_kind not in model.SUBJECT_KINDS or dst_kind not in model.OBJECT_KINDS:
            skipped += 1
            continue
        col_set = ", ".join(f"`{c}`: ${c}" for c in extra_cols)
        cypher = (
            f"MATCH (a:`{src_kind}` {{id: $from_id}}), "
            f"      (b:`{dst_kind}` {{id: $to_id}}) "
            f"CREATE (a)-[r:`{predicate}` {{{col_set}}}]->(b);"
        )
        params = {"from_id": from_id, "to_id": to_id}
        for c in extra_cols:
            params[c] = _coerce(row.get(c))
        try:
            conn.execute(cypher, parameters=params)
            n += 1
        except Exception as e:
            print(f"  ! insert {predicate} {from_id}->{to_id}: {e}")
    if skipped:
        print(f"    ({predicate}: {skipped} edges skipped — orphan endpoints or kind mismatch)")
    return n


def _coerce(v):
    """Coerce a JSONL value to something Kuzu can store. Lists of primitives
    pass through; complex things go to JSON strings. Empty lists/dicts are
    NULL because Kuzu cannot infer the element type of an empty list during
    parameter binding."""
    if v is None:
        return None
    if isinstance(v, (str, bool, int, float)):
        return v
    if isinstance(v, list):
        if not v:
            return None
        if all(isinstance(x, (str, int, float, bool)) for x in v):
            return v
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, dict):
        if not v:
            return None
        return json.dumps(v, ensure_ascii=False)
    return str(v)


if __name__ == "__main__":
    raise SystemExit(main())
