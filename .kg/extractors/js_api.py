"""Extract the canonical JS API surface from `js-api/src/`.

Three passes:

1. **Hand-written TS** (everything except `*.g.ts`): regex over each file to
   pull out `export class`, `export enum`, `export const`, `export interface`,
   `export type`, plus methods/getters per class. Light-weight — we don't run
   the full TS compiler.

2. **`Dapi` class** in `dapi.ts`: parse getter declarations of the form
   `get <name>(): <DataSourceClass>` to enumerate `grok.dapi.<endpoint>`
   accessors and the entity class each one is typed for.

3. **`api/grok_api.g.ts` IDartApi interface**: every line is a `grok_*`
   stub — emit `GeneratedBinding` nodes with parsed Dart class + member.

Namespaces (`grok`, `ui`, `DG`, plus nested) are seeded from the source-file
location and from any explicit `namespace`-style usage observed.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import (
    DapiEndpoint, GeneratedBinding, JsApiNamespace, JsEventStream,
    TsClass, TsConstant, TsEnum, TsInterface, TsMethod, UiComponent,
)
from schema.relations import (
    DefinesSemtype, ExtendsClass, HasMethod, InNamespace, TypedFor,
)

from . import _common as c

EXTRACTOR_NAME = "js_api"

JS_API_ROOT = c.REPO_ROOT / "js-api"

# ---------------------------------------------------------------------------
# Regexes (light-weight — works for ~90% of declarations per the survey)
# ---------------------------------------------------------------------------

CLASS_RE = re.compile(
    r"^export\s+(?P<abstract>abstract\s+)?class\s+(?P<name>[A-Z][\w]*)"
    r"(?:\s+extends\s+(?P<extends>[A-Z][\w.<>, ]*?))?"
    r"\s*[{<]",
    re.MULTILINE,
)
ENUM_RE = re.compile(r"^export\s+enum\s+(?P<name>[A-Z][A-Z0-9_]*)\s*\{",
                     re.MULTILINE)
CONST_OBJ_RE = re.compile(
    r"^export\s+const\s+(?P<name>[A-Z][A-Z0-9_]*)\s*[:=]?\s*(?:[^=]*=\s*)?\{",
    re.MULTILINE,
)
TYPE_RE = re.compile(r"^export\s+type\s+(?P<name>[A-Z][\w]*)\s*=", re.MULTILINE)
INTERFACE_RE = re.compile(r"^export\s+interface\s+(?P<name>[A-Z][\w]*)\s*[<{]",
                          re.MULTILINE)

# Method / getter inside a class body. Permissive about modifiers and
# allows parameter lists to span multiple lines with nested parens (e.g.
# decorator-annotated params).  Two levels of paren-balancing is enough.
METHOD_RE = re.compile(
    r"""^[\t ]+
        (?:public\s+|private\s+|protected\s+|override\s+|readonly\s+)*
        (?P<static>static\s+)?
        (?:public\s+|private\s+|protected\s+|override\s+|readonly\s+)*
        (?P<async>async\s+)?
        (?P<getter>get\s+|set\s+)?
        (?!if\b|while\b|for\b|switch\b|return\b|catch\b|else\b|new\b|class\b|function\b|interface\b|throw\b|do\b|try\b|delete\b|typeof\b)
        (?P<name>[a-z_$][\w$]*)
        \s*\(
        (?P<params>(?:[^()]|\((?:[^()]|\([^()]*\))*\))*)
        \)
        (?:\s*:\s*(?P<ret>[^{=;]+?))?
        \s*\{
    """, re.VERBOSE | re.MULTILINE)

# `api.grok_*(...)` call inside a method body
DELEGATES_RE = re.compile(r"\bapi\.(grok_[A-Za-z0-9_]+)\s*\(")

# Dapi getter: `get users(): UsersDataSource {`
DAPI_GET_RE = re.compile(
    r"^\s*get\s+(?P<name>\w+)\s*\(\)\s*:\s*(?P<ds>[A-Z][\w]*DataSource)",
    re.MULTILINE,
)
# `class UsersDataSource extends HttpDataSource<User>`
DS_DECL_RE = re.compile(
    r"^export\s+class\s+(?P<ds>[A-Z][\w]*DataSource)\s+extends\s+HttpDataSource<(?P<entity>[A-Z][\w]*)>",
    re.MULTILINE,
)

# `IDartApi` interface body — each line `grok_<class>_<member>(...)`
DART_API_RE = re.compile(r"^\s+(?P<name>grok_[A-Za-z0-9_]+)\s*\(",
                         re.MULTILINE)

# Events.ts getter: `get onTableAdded(): Observable<…> {` (with greedy
# generic-parameter capture so `EventData<DataFrameArgs>` and similar
# nested generics survive). Also accept getters whose return type is just
# `rxjs.Observable<...>` or `: Observable<…>;` (no body).
EVENT_RE = re.compile(
    r"""^\s*
        get\s+
        (?P<name>on[A-Z][\w]*)\s*
        \(\)\s*:\s*
        (?:rxjs\.)?Observable\s*<\s*(?P<args>[^>]+(?:<[^>]*>[^>]*)*?)\s*>
    """,
    re.MULTILINE | re.VERBOSE,
)

# ui.ts: `export function dialog(...)` and arrow `export const button = ...`
UI_FUNC_RE = re.compile(
    r"^export\s+function\s+(?P<name>[a-z_$][\w$]*)\s*\((?P<params>[^)]*)\)"
    r"(?:\s*:\s*(?P<ret>[^{=;]+?))?\s*\{",
    re.MULTILINE,
)
UI_CONST_RE = re.compile(
    r"^export\s+const\s+(?P<name>[a-z_$][\w$]*)\s*=\s*\(",
    re.MULTILINE,
)

# `EVENT_TYPE.X = '...'` not directly used here, but SEMTYPE values matter:
SEMTYPE_VAL_RE = re.compile(
    r"^\s*([A-Z][A-Z0-9_]*)\s*[:=]\s*['\"]([^'\"]+)['\"]",
    re.MULTILINE,
)


# ---------------------------------------------------------------------------
# Namespace inference
# ---------------------------------------------------------------------------

# File path → top-level namespace heuristic
def _ns_for_path(rel: str) -> str:
    """Crude heuristic: 'src/dataframe/...' → DG (DG re-exports), 'src/dapi.ts'
    → grok.dapi-related, 'src/ui.ts' → ui."""
    if "/dapi.ts" in rel:    return "grok.dapi"
    if "/events.ts" in rel:  return "grok.events"
    if "/functions.ts" in rel: return "grok.functions"
    if "/shell.ts" in rel:   return "grok.shell"
    if "/ui.ts" in rel:      return "ui"
    return "DG"


def _ns_id(path: str) -> str:
    return f"jsns:{path}"


# ---------------------------------------------------------------------------
# Pass A — walk hand-written TS
# ---------------------------------------------------------------------------

def _emit_namespaces(bundle, seen: set[str]) -> None:
    """Seed the standard namespaces."""
    for path in ("grok", "ui", "DG",
                 "grok.dapi", "grok.events", "grok.functions",
                 "grok.shell", "DG.chem", "DG.bio", "ui.dialog",
                 "ui.input"):
        if path in seen: continue
        seen.add(path)
        parent = ".".join(path.split(".")[:-1]) if "." in path else None
        bundle.add(JsApiNamespace(
            id=_ns_id(path),
            name=path,
            source_layer=SourceLayer.SYNTHETIC,
            dotted_path=path,
            parent=_ns_id(parent) if parent else None,
        ))


def _process_class(text: str, file_path: Path, ns: str, bundle,
                    classes_index: dict[str, str]) -> int:
    """Find each `export class X extends Y` and emit JsApiClass + ExtendsClass."""
    n = 0
    # Build per-class slice so we can attribute methods to the right class
    matches = list(CLASS_RE.finditer(text))
    for i, m in enumerate(matches):
        name = m.group("name")
        cid = f"ts:class:DG.{name}"
        classes_index[name] = cid
        extends_raw = m.group("extends") or ""
        extends_name = extends_raw.split("<", 1)[0].strip().split(",")[0].strip() or None

        # Slice out the class body up to the next class declaration
        body_start = m.end()
        body_end = matches[i+1].start() if i + 1 < len(matches) else len(text)
        body = text[body_start:body_end]

        # Count methods within this body
        method_count = sum(1 for _ in METHOD_RE.finditer(body))

        bundle.add(TsClass(
            id=cid, name=name, source_layer=SourceLayer.SYNTHETIC,
            namespace=ns, dotted_name=f"DG.{name}",
            is_abstract=bool(m.group("abstract")),
            extends_class=extends_name,
            method_count=method_count,
            is_jsapi=True,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1,
                              role="definition")],
        ))
        bundle.add(InNamespace(from_id=cid, to_id=_ns_id(ns),
                               derived_by=Provenance.AST))

        # Methods + delegates
        for mm in METHOD_RE.finditer(body):
            mname = mm.group("name")
            if mname in {"constructor", "if", "while", "for", "switch", "return"}:
                continue
            mid = f"ts:method:DG.{name}.{mname}"
            ret = (mm.group("ret") or "").strip()
            params = mm.group("params") or ""
            param_count = 0 if not params.strip() else (params.count(",") + 1)

            # Look at the method body for `api.grok_*` delegate
            method_text = body[mm.end(): mm.end() + 800]   # 800 chars look-ahead
            delegate = None
            dm = DELEGATES_RE.search(method_text)
            if dm:
                delegate = dm.group(1)

            bundle.add(TsMethod(
                id=mid, name=mname, source_layer=SourceLayer.SYNTHETIC,
                class_id=cid,
                is_static=bool(mm.group("static")),
                is_getter=bool(mm.group("getter")),
                return_type=ret or None,
                parameter_count=param_count,
                delegates_to=delegate,
                paths=[c.file_ref(file_path,
                                  line_start=body[:mm.start()].count("\n") + text[:body_start].count("\n") + 1)],
            ))
            bundle.add(HasMethod(from_id=cid, to_id=mid,
                                 derived_by=Provenance.AST))
        n += 1

    # Resolve ExtendsClass edges (forward-reference safe via classes_index)
    return n


def _emit_extends_pass(text: str, classes_index: dict[str, str], bundle) -> None:
    """Second pass to wire ExtendsClass once all classes are known."""
    for m in CLASS_RE.finditer(text):
        name = m.group("name")
        ext = m.group("extends")
        if not ext: continue
        ext_clean = ext.split("<", 1)[0].strip().split(",")[0].strip()
        cid = classes_index.get(name)
        eid = classes_index.get(ext_clean)
        if cid and eid:
            bundle.add(ExtendsClass(from_id=cid, to_id=eid,
                                    derived_by=Provenance.AST))


def _process_enums_consts(text: str, file_path: Path, ns: str, bundle) -> int:
    n = 0
    for m in ENUM_RE.finditer(text):
        name = m.group("name")
        eid = f"ts:enum:DG.{name}"
        bundle.add(TsEnum(
            id=eid, name=name, source_layer=SourceLayer.SYNTHETIC,
            namespace=ns, is_jsapi=True,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(InNamespace(from_id=eid, to_id=_ns_id(ns),
                               derived_by=Provenance.AST))
        n += 1
    for m in CONST_OBJ_RE.finditer(text):
        name = m.group("name")
        cid = f"ts:const:DG.{name}"
        bundle.add(TsConstant(
            id=cid, name=name, source_layer=SourceLayer.SYNTHETIC,
            namespace=ns, is_jsapi=True,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(InNamespace(from_id=cid, to_id=_ns_id(ns),
                               derived_by=Provenance.AST))
        n += 1
        # If this is SEMTYPE / UNITS / TAGS, parse member values to bridge
        if name in {"SEMTYPE", "UNITS", "TAGS"}:
            # Slice out body up to closing brace
            body_start = m.end()
            depth = 1
            i = body_start
            while i < len(text) and depth > 0:
                if text[i] == "{": depth += 1
                elif text[i] == "}": depth -= 1
                i += 1
            body = text[body_start:i]
            for vm in SEMTYPE_VAL_RE.finditer(body):
                value = vm.group(2)
                if name == "SEMTYPE":
                    bundle.add(DefinesSemtype(
                        from_id=cid, to_id=f"semtype:{value}",
                        derived_by=Provenance.AST,
                    ))
    return n


def _process_interfaces_types(text: str, file_path: Path, ns: str, bundle) -> int:
    n = 0
    for m in INTERFACE_RE.finditer(text):
        name = m.group("name")
        iid = f"ts:iface:DG.{name}"
        bundle.add(TsInterface(
            id=iid, name=name, source_layer=SourceLayer.SYNTHETIC,
            namespace=ns, is_type_alias=False, is_jsapi=True,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(InNamespace(from_id=iid, to_id=_ns_id(ns),
                               derived_by=Provenance.AST))
        n += 1
    for m in TYPE_RE.finditer(text):
        name = m.group("name")
        iid = f"ts:iface:DG.{name}"
        bundle.add(TsInterface(
            id=iid, name=name, source_layer=SourceLayer.SYNTHETIC,
            namespace=ns, is_type_alias=True, is_jsapi=True,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(InNamespace(from_id=iid, to_id=_ns_id(ns),
                               derived_by=Provenance.AST))
        n += 1
    return n


# ---------------------------------------------------------------------------
# Pass B — Dapi endpoints
# ---------------------------------------------------------------------------

def _process_dapi(repo_root: Path, bundle, classes_index: dict[str, str]) -> int:
    dapi_file = repo_root / "js-api" / "src" / "dapi.ts"
    if not dapi_file.is_file():
        return 0
    text = dapi_file.read_text(encoding="utf-8", errors="replace")

    # Map DataSource class -> served entity (T in HttpDataSource<T>)
    ds_to_entity: dict[str, str] = {}
    for m in DS_DECL_RE.finditer(text):
        ds_to_entity[m.group("ds")] = m.group("entity")

    n = 0
    for m in DAPI_GET_RE.finditer(text):
        name = m.group("name")
        ds = m.group("ds")
        eid = f"jsapi:dapi:{name}"
        served = ds_to_entity.get(ds)
        bundle.add(DapiEndpoint(
            id=eid, name=name, source_layer=SourceLayer.SYNTHETIC,
            accessor_path=f"grok.dapi.{name}",
            data_source_class=ds, served_entity=served,
            paths=[c.file_ref(dapi_file, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(InNamespace(from_id=eid, to_id=_ns_id("grok.dapi"),
                               derived_by=Provenance.AST))
        if served and (cid := classes_index.get(served)):
            bundle.add(TypedFor(from_id=eid, to_id=cid,
                                derived_by=Provenance.AST))
        n += 1
    return n


# ---------------------------------------------------------------------------
# Pass C — Generated bindings (IDartApi)
# ---------------------------------------------------------------------------

def _process_generated_bindings(repo_root: Path, bundle) -> int:
    g_file = repo_root / "js-api" / "src" / "api" / "grok_api.g.ts"
    if not g_file.is_file():
        return 0
    text = g_file.read_text(encoding="utf-8", errors="replace")
    n = 0
    for m in DART_API_RE.finditer(text):
        name = m.group("name")     # grok_DataFrame_Columns
        # Split into Dart class + member by underscore (best-effort)
        parts = name.split("_", 2)
        dart_class = parts[1] if len(parts) > 1 else None
        dart_member = parts[2] if len(parts) > 2 else None
        gid = f"jsapi:bind:{name}"
        bundle.add(GeneratedBinding(
            id=gid, name=name, source_layer=SourceLayer.SYNTHETIC,
            dart_class=dart_class, dart_member=dart_member,
            paths=[c.file_ref(g_file, line_start=text[:m.start()].count("\n") + 1)],
        ))
        n += 1
    return n


# ---------------------------------------------------------------------------
# Pass D — events
# ---------------------------------------------------------------------------

def _process_events(repo_root: Path, bundle) -> int:
    f = repo_root / "js-api" / "src" / "events.ts"
    if not f.is_file(): return 0
    text = f.read_text(encoding="utf-8", errors="replace")
    n = 0
    for m in EVENT_RE.finditer(text):
        name = m.group("name")
        eid = f"jsapi:event:{name}"
        bundle.add(JsEventStream(
            id=eid, name=name, source_layer=SourceLayer.SYNTHETIC,
            accessor_path=f"grok.events.{name}",
            event_args_type=m.group("args").strip(),
            paths=[c.file_ref(f, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(InNamespace(from_id=eid, to_id=_ns_id("grok.events"),
                               derived_by=Provenance.AST))
        n += 1
    return n


# ---------------------------------------------------------------------------
# Pass E — UI components
# ---------------------------------------------------------------------------

def _process_ui(repo_root: Path, bundle) -> int:
    # ui.ts lives at js-api/ui.ts (top-level), not under src/. Other public
    # modules (dg.ts, grok.ts, ui.ts, datagrok.ts) follow the same pattern.
    for candidate in (repo_root / "js-api" / "ui.ts",
                      repo_root / "js-api" / "src" / "ui.ts"):
        if candidate.is_file():
            f = candidate
            break
    else:
        return 0
    text = f.read_text(encoding="utf-8", errors="replace")
    n = 0
    for m in UI_FUNC_RE.finditer(text):
        name = m.group("name")
        if name.startswith("_"):
            continue
        cid = f"jsapi:ui:{name}"
        params = m.group("params") or ""
        ret = (m.group("ret") or "").strip()
        bundle.add(UiComponent(
            id=cid, name=name, source_layer=SourceLayer.SYNTHETIC,
            factory_name=name, return_type=ret or None,
            parameter_count=0 if not params.strip() else params.count(",") + 1,
            paths=[c.file_ref(f, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(InNamespace(from_id=cid, to_id=_ns_id("ui"),
                               derived_by=Provenance.AST))
        n += 1
    return n


# ---------------------------------------------------------------------------
# Entry
# ---------------------------------------------------------------------------

def _iter_jsapi_ts(repo_root: Path):
    """Yield all hand-written js-api TS files: js-api/src/**.ts plus the
    top-level entries (dg.ts, grok.ts, ui.ts, datagrok.ts). Skips
    src/api/*.g.ts (Dart bindings) and any .d.ts."""
    src = repo_root / "js-api" / "src"
    if src.is_dir():
        for ts in sorted(src.rglob("*.ts")):
            if "/api/" in ts.as_posix() or ts.name.endswith(".g.ts"):
                continue
            if ts.name.endswith(".d.ts"):
                continue
            yield ts
    root = repo_root / "js-api"
    if root.is_dir():
        for ts in sorted(root.glob("*.ts")):
            if ts.name.endswith(".d.ts") or ts.name.endswith(".g.ts"):
                continue
            yield ts


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    if not (repo_root / "js-api").is_dir():
        print(f"[{EXTRACTOR_NAME}] no js-api/ at {repo_root}")
        return

    seen_ns: set[str] = set()
    _emit_namespaces(bundle, seen_ns)

    classes_index: dict[str, str] = {}
    n_class = n_enum = n_iface = 0

    # Pass A: hand-written TS — src/**.ts + top-level entries
    for ts_file in _iter_jsapi_ts(repo_root):
        try:
            text = ts_file.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        rel = c.repo_rel(ts_file)
        ns = _ns_for_path(rel)
        n_class += _process_class(text, ts_file, ns, bundle, classes_index)
        n_enum += _process_enums_consts(text, ts_file, ns, bundle)
        n_iface += _process_interfaces_types(text, ts_file, ns, bundle)

    # Resolve ExtendsClass after all classes known
    for ts_file in _iter_jsapi_ts(repo_root):
        try:
            text = ts_file.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        _emit_extends_pass(text, classes_index, bundle)

    n_dapi = _process_dapi(repo_root, bundle, classes_index)
    n_bind = _process_generated_bindings(repo_root, bundle)
    n_event = _process_events(repo_root, bundle)
    n_ui = _process_ui(repo_root, bundle)

    print(f"[{EXTRACTOR_NAME}] {n_class} classes, {n_enum} enums/consts, "
          f"{n_iface} interfaces/types, {n_dapi} dapi endpoints, "
          f"{n_bind} generated bindings, {n_event} event streams, "
          f"{n_ui} ui components")
