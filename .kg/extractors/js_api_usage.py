"""Scan plugin (and library) TS source for JS API references and emit
edges at THREE granularities: File, containing Method/Function, and Package.

The same regex pass that previously aggregated to package level now also
attributes each match to the source File and (where the line falls inside
a method body) the containing TsMethod or top-level TsFunction. This lets
queries like

  MATCH (m:TsMethod {name:'editMoleculeCell'})-[u:USES_UI_COMPONENT]->(c:UiComponent)
  RETURN c.name, u.use_count

answer "what UI components does *this method* build", not just "which package".

Containment lookup is line-based: build `(file_id, sorted [(line, method_id)])`
from already-extracted TsMethod / TsFunction `paths[0].line_start`. For each
match at line N, pick the entity with the largest line_start ≤ N. TS doesn't
have nested function declarations of any depth, so a flat last-wins lookup
is correct ~99% of the time. False attributions land on the immediately
preceding sibling, which is recoverable; method count is correct because we
de-duplicate on (subject, target).

Patterns:
  - `import * as DG from 'datagrok-api/dg'`            → ImportsNamespace(DG)
  - `import * as grok from 'datagrok-api/grok'`        → ImportsNamespace(grok)
  - `import * as ui from 'datagrok-api/ui'`            → ImportsNamespace(ui)
  - `DG.<JsApiClass>`     (resolved to known classes)  → UsesApiClass
  - `DG.<JsApiEnum>.X`    (resolved to known enums)    → UsesApiEnum
  - `DG.<JsApiConstant>.X`(resolved to known consts)   → UsesApiEnum
  - `grok.dapi.<endpoint>`(any segment in the chain)   → CallsDapiEndpoint
  - `grok.events.<event>` (resolved to known streams)  → SubscribesToEvent
  - `ui.<factory>(`       (resolved to known ui funcs) → UsesUiComponent

Resolution requires the JS API extractor to have run first.
"""

from __future__ import annotations

import bisect
import json
import re
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.relations import (
    CallsDapiEndpoint, ImportsNamespace, SubscribesToEvent,
    UsesApiClass, UsesApiEnum, UsesUiComponent,
)

from . import _common as c

EXTRACTOR_NAME = "js_api_usage"

NODES_DIR = c.REPO_ROOT / ".kg" / "data" / "nodes"

IMPORT_RE = re.compile(
    r"""import\s+\*\s+as\s+(?P<alias>[A-Za-z_]\w*)\s+from\s+
        ['"]datagrok-api/(?P<which>dg|grok|ui)['"]
    """, re.VERBOSE,
)

# `DG.<Identifier>` or `DG.<ns>.<Identifier>`
DG_REF_RE = re.compile(r"\bDG\.(?P<chain>[A-Za-z_][\w.]*)")
# Walk the full dotted chain after `grok.dapi.` so multi-segment endpoints
# like `grok.dapi.docker.dockerContainers` register both segments.
GROK_DAPI_RE = re.compile(r"\bgrok\.dapi\.(?P<chain>[A-Za-z_][\w.]*)")
GROK_EVENTS_RE = re.compile(r"\bgrok\.events\.(?P<name>on[A-Z]\w*)")
UI_CALL_RE = re.compile(r"\bui\.(?P<name>[a-z_$][\w$]*)\s*\(")


def _read_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _build_indexes() -> tuple[dict, dict, dict, dict, dict, dict]:
    """Return name→id maps for the JS API resolution targets."""
    def jsapi_only(rows: list[dict]) -> dict:
        return {n["name"]: n["id"] for n in rows if n.get("is_jsapi")}

    classes = jsapi_only(_read_jsonl(NODES_DIR / "TsClass.jsonl"))
    enums = jsapi_only(_read_jsonl(NODES_DIR / "TsEnum.jsonl"))
    consts = jsapi_only(_read_jsonl(NODES_DIR / "TsConstant.jsonl"))
    dapi = {n["name"]: n["id"] for n in _read_jsonl(NODES_DIR / "DapiEndpoint.jsonl")}
    events = {n["name"]: n["id"] for n in _read_jsonl(NODES_DIR / "JsEventStream.jsonl")}
    ui = {n["name"]: n["id"] for n in _read_jsonl(NODES_DIR / "UiComponent.jsonl")}
    return classes, enums, consts, dapi, events, ui


def _build_containing_index() -> dict[str, list[tuple[int, str]]]:
    """{file_relative_path: sorted list of (line_start, entity_id)} for every
    TsMethod and TsFunction. Used to look up the containing entity at a
    given match line via bisect. Methods come before standalone functions
    when both share a line_start (rare; resolved arbitrarily by sort)."""
    out: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for kind in ("TsMethod.jsonl", "TsFunction.jsonl"):
        for row in _read_jsonl(NODES_DIR / kind):
            paths = row.get("paths") or []
            if not paths:
                continue
            p = paths[0]
            rel = p.get("path"); line = p.get("line_start")
            mid = row.get("id")
            if rel and isinstance(line, int) and mid:
                out[rel].append((line, mid))
    for rel, lst in out.items():
        lst.sort(key=lambda t: t[0])
    return out


def _containing(file_rel: str, line: int,
                idx: dict[str, list[tuple[int, str]]]) -> str | None:
    arr = idx.get(file_rel)
    if not arr:
        return None
    # Find the largest line_start <= line
    pos = bisect.bisect_right(arr, (line, "￿")) - 1
    if pos < 0:
        return None
    return arr[pos][1]


def _scan_pkg(pkg_dir: Path, owner_id: str,
              classes: dict, enums: dict, consts: dict,
              dapi: dict, events: dict, ui: dict,
              containing_idx: dict[str, list[tuple[int, str]]],
              # aggregations keyed by (subject_id, target_id)
              imports_ns: dict, uses_class: dict, uses_enum: dict,
              calls_dapi: dict, subs_event: dict, uses_ui: dict) -> int:
    src_dirs = [pkg_dir / "src"]
    n_files = 0
    for sd in src_dirs:
        if not sd.is_dir():
            continue
        for ts in sd.rglob("*.ts"):
            if "node_modules" in ts.parts:
                continue
            try:
                text = ts.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            n_files += 1
            rel = c.repo_rel(ts)
            file_id = f"file:{rel}"

            def _bump(agg: dict, target: str, line: int):
                # Three granularities: package, file, containing method
                agg[(owner_id, target)] += 1
                agg[(file_id, target)] += 1
                m = _containing(rel, line, containing_idx)
                if m:
                    agg[(m, target)] += 1

            # Imports — emit at namespace alias level
            for m in IMPORT_RE.finditer(text):
                alias = m.group("which")
                ns_key = "DG" if alias == "dg" else alias
                line = text[:m.start()].count("\n") + 1
                # Imports use a 2-tuple-arg helper since they have a different
                # value space (namespace IDs are strings, not entity IDs)
                imports_ns[(owner_id, ns_key)] += 1
                imports_ns[(file_id, ns_key)] += 1
                # Containing method/function for an import is the file itself,
                # not a method, so we skip the third granularity here.

            # DG.<X>...
            for m in DG_REF_RE.finditer(text):
                chain = m.group("chain").split(".", 1)[0]
                line = text[:m.start()].count("\n") + 1
                if chain in classes:
                    _bump(uses_class, classes[chain], line)
                elif chain in enums:
                    _bump(uses_enum, enums[chain], line)
                elif chain in consts:
                    _bump(uses_enum, consts[chain], line)

            # grok.dapi.<seg1>.<seg2>...
            for m in GROK_DAPI_RE.finditer(text):
                chain = m.group("chain")
                line = text[:m.start()].count("\n") + 1
                for seg in chain.split("."):
                    if seg in dapi:
                        _bump(calls_dapi, dapi[seg], line)

            # grok.events.<stream>
            for m in GROK_EVENTS_RE.finditer(text):
                name = m.group("name")
                line = text[:m.start()].count("\n") + 1
                if name in events:
                    _bump(subs_event, events[name], line)

            # ui.<factory>(
            for m in UI_CALL_RE.finditer(text):
                name = m.group("name")
                line = text[:m.start()].count("\n") + 1
                if name in ui:
                    _bump(uses_ui, ui[name], line)
    return n_files


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    classes, enums, consts, dapi, events, ui = _build_indexes()
    if not classes:
        print(f"[{EXTRACTOR_NAME}] no TsClass(is_jsapi=true) nodes - run `js_api` extractor first")
        return
    containing_idx = _build_containing_index()

    ns_id_for = {"DG": "jsns:DG", "grok": "jsns:grok", "ui": "jsns:ui"}

    imports_ns: dict = defaultdict(int)
    uses_class: dict = defaultdict(int)
    uses_enum: dict = defaultdict(int)
    calls_dapi: dict = defaultdict(int)
    subs_event: dict = defaultdict(int)
    uses_ui: dict = defaultdict(int)

    n_files = 0
    pkgs = repo_root / "packages"
    if pkgs.is_dir():
        for pkg_dir in sorted(pkgs.iterdir()):
            if not pkg_dir.is_dir(): continue
            if package_filter and pkg_dir.name not in package_filter: continue
            n_files += _scan_pkg(
                pkg_dir, c.pkg_id(pkg_dir.name),
                classes, enums, consts, dapi, events, ui,
                containing_idx,
                imports_ns, uses_class, uses_enum, calls_dapi, subs_event, uses_ui,
            )

    libs = repo_root / "libraries"
    if libs.is_dir():
        for lib_dir in sorted(libs.iterdir()):
            if not lib_dir.is_dir(): continue
            n_files += _scan_pkg(
                lib_dir, c.lib_id(lib_dir.name),
                classes, enums, consts, dapi, events, ui,
                containing_idx,
                imports_ns, uses_class, uses_enum, calls_dapi, subs_event, uses_ui,
            )

    n_emit = 0
    by_kind = defaultdict(int)

    def _kind_of(subject_id: str) -> str:
        if subject_id.startswith("pkg:"):    return "pkg"
        if subject_id.startswith("lib:"):    return "lib"
        if subject_id.startswith("file:"):   return "file"
        if subject_id.startswith("ts:method:"): return "method"
        if subject_id.startswith("ts:fn:"):  return "fn"
        return "other"

    for (caller, alias), count in imports_ns.items():
        ns = ns_id_for.get(alias)
        if not ns: continue
        bundle.add(ImportsNamespace(
            from_id=caller, to_id=ns, import_count=count,
            derived_by=Provenance.AST,
        ))
        n_emit += 1; by_kind[("imports", _kind_of(caller))] += 1
    for (caller, target), count in uses_class.items():
        bundle.add(UsesApiClass(from_id=caller, to_id=target, use_count=count,
                                derived_by=Provenance.AST))
        n_emit += 1; by_kind[("class", _kind_of(caller))] += 1
    for (caller, target), count in uses_enum.items():
        bundle.add(UsesApiEnum(from_id=caller, to_id=target, use_count=count,
                               derived_by=Provenance.AST))
        n_emit += 1; by_kind[("enum", _kind_of(caller))] += 1
    for (caller, target), count in calls_dapi.items():
        bundle.add(CallsDapiEndpoint(from_id=caller, to_id=target, call_count=count,
                                     derived_by=Provenance.AST))
        n_emit += 1; by_kind[("dapi", _kind_of(caller))] += 1
    for (caller, target), count in subs_event.items():
        bundle.add(SubscribesToEvent(from_id=caller, to_id=target, subscribe_count=count,
                                     derived_by=Provenance.AST))
        n_emit += 1; by_kind[("event", _kind_of(caller))] += 1
    for (caller, target), count in uses_ui.items():
        bundle.add(UsesUiComponent(from_id=caller, to_id=target, use_count=count,
                                   derived_by=Provenance.AST))
        n_emit += 1; by_kind[("ui", _kind_of(caller))] += 1

    breakdown = ", ".join(f"{kind}/{level}: {n}"
                          for (kind, level), n in sorted(by_kind.items()))
    print(f"[{EXTRACTOR_NAME}] scanned {n_files} files, emitted {n_emit} usage edges; "
          f"{breakdown}")
