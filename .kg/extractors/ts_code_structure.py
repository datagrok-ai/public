"""Extract TS code structure (classes, methods, interfaces, enums, exported
functions) from `packages/*/src/` and `libraries/*/src/`.

The same shape as `js_api.py` but uses `is_jsapi=False` and `source_layer
= PLUGINS | LIBRARIES`. Lets you ask 'what classes does Chem declare?'
'which interfaces does Bio define?' 'where are top-level helper functions
that aren't platform-registered?' alongside the existing js-api queries.

Edges emitted:
  - HasClass     : Package/Library -> TsClass | TsInterface | TsEnum | TsConstant
  - HasMethod    : TsClass -> TsMethod
  - HasFunction  : Package/Library -> TsFunction
  - ExtendsClass : TsClass -> TsClass (when both ends are extracted)
"""

from __future__ import annotations

import re
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import (
    TsClass, TsConstant, TsEnum, TsFunction, TsInterface, TsMethod,
)
from schema.relations import (
    ExtendsClass, HasClass, HasFunction, HasMethod, ImplementsInterface,
)

from . import _common as c

EXTRACTOR_NAME = "ts_code_structure"

# Reuse the same regexes from js_api.py — duplicated here so the modules are independent.
CLASS_RE = re.compile(
    r"^export\s+(?P<abstract>abstract\s+)?class\s+(?P<name>[A-Za-z_][\w]*)"
    r"(?:\s+extends\s+(?P<extends>[A-Za-z_][\w.<>, ]*?))?"
    r"(?:\s+implements\s+(?P<implements>[A-Za-z_][\w.<>, ]*?))?"
    r"\s*[{<]",
    re.MULTILINE,
)
ENUM_RE = re.compile(r"^export\s+(?:const\s+)?enum\s+(?P<name>[A-Za-z_][\w]*)\s*\{",
                     re.MULTILINE)
CONST_OBJ_RE = re.compile(
    r"^export\s+const\s+(?P<name>[A-Za-z_][\w]*)\s*[:=]?\s*(?:[^=]*=\s*)?\{",
    re.MULTILINE,
)
TYPE_RE = re.compile(r"^export\s+type\s+(?P<name>[A-Za-z_][\w]*)\s*=", re.MULTILINE)
INTERFACE_RE = re.compile(r"^export\s+interface\s+(?P<name>[A-Za-z_][\w]*)\s*[<{]",
                          re.MULTILINE)
FUNCTION_RE = re.compile(
    r"^export\s+(?P<async>async\s+)?function\s+(?P<name>[A-Za-z_][\w]*)\s*"
    r"\((?P<params>[^)]*)\)"
    r"(?:\s*:\s*(?P<ret>[^{=;]+?))?\s*\{",
    re.MULTILINE,
)
# Method or getter inside a class body. Permissive about modifiers (static,
# async, public/private/protected, override, get/set) and allows the
# parameter list to span multiple lines and contain nested parens — so
# decorator-annotated params like `@grok.decorators.param({...}) x: T`
# parse correctly. Two levels of paren-balancing in the params capture is
# enough for everything in this codebase (decorators with object literals).
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


def _ts_id(kind: str, namespace: str, name: str) -> str:
    """Stable ID; namespace is `pkg:Chem` or `lib:utils`."""
    return f"ts:{kind}:{namespace}.{name}"


def _process_file(text: str, file_path: Path,
                   ns_id: str, ns_dotted: str, owner_kind: str,
                   bundle, classes_index: dict[str, str]) -> dict[str, int]:
    counts = {"class": 0, "method": 0, "interface": 0, "enum": 0,
              "constant": 0, "function": 0}
    layer = SourceLayer.LIBRARIES if owner_kind == "lib" else SourceLayer.PLUGINS

    # Classes (and methods inside)
    matches = list(CLASS_RE.finditer(text))
    for i, m in enumerate(matches):
        name = m.group("name")
        cid = _ts_id("class", ns_dotted, name)
        classes_index[f"{ns_dotted}.{name}"] = cid
        classes_index[name] = cid              # also unprefixed for in-package extends

        body_start = m.end()
        body_end = matches[i+1].start() if i + 1 < len(matches) else len(text)
        body = text[body_start:body_end]
        method_count = sum(1 for _ in METHOD_RE.finditer(body))

        bundle.add(TsClass(
            id=cid, name=name, source_layer=layer,
            namespace=ns_dotted, dotted_name=f"{ns_dotted}.{name}",
            is_abstract=bool(m.group("abstract")),
            extends_class=(m.group("extends") or "").split("<", 1)[0].strip() or None,
            method_count=method_count,
            is_jsapi=False,
            paths=[c.file_ref(file_path,
                              line_start=text[:m.start()].count("\n") + 1,
                              role="definition")],
        ))
        bundle.add(HasClass(from_id=ns_id, to_id=cid,
                            derived_by=Provenance.AST))
        counts["class"] += 1

        for mm in METHOD_RE.finditer(body):
            mname = mm.group("name")
            if mname in {"constructor", "if", "while", "for", "switch", "return",
                          "catch", "else"}:
                continue
            mid = _ts_id("method", ns_dotted, f"{name}.{mname}")
            params = mm.group("params") or ""
            param_count = 0 if not params.strip() else (params.count(",") + 1)
            line_start = body[:mm.start()].count("\n") + text[:body_start].count("\n") + 1
            bundle.add(TsMethod(
                id=mid, name=mname, source_layer=layer,
                class_id=cid,
                is_static=bool(mm.group("static")),
                is_getter=bool(mm.group("getter")),
                return_type=(mm.group("ret") or "").strip() or None,
                parameter_count=param_count,
                paths=[c.file_ref(file_path, line_start=line_start)],
            ))
            bundle.add(HasMethod(from_id=cid, to_id=mid,
                                 derived_by=Provenance.AST))
            counts["method"] += 1

    # Interfaces / type aliases
    for m in INTERFACE_RE.finditer(text):
        name = m.group("name")
        iid = _ts_id("iface", ns_dotted, name)
        bundle.add(TsInterface(
            id=iid, name=name, source_layer=layer,
            namespace=ns_dotted, is_type_alias=False, is_jsapi=False,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(HasClass(from_id=ns_id, to_id=iid,
                            derived_by=Provenance.AST))
        counts["interface"] += 1
    for m in TYPE_RE.finditer(text):
        name = m.group("name")
        iid = _ts_id("iface", ns_dotted, name)
        bundle.add(TsInterface(
            id=iid, name=name, source_layer=layer,
            namespace=ns_dotted, is_type_alias=True, is_jsapi=False,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(HasClass(from_id=ns_id, to_id=iid,
                            derived_by=Provenance.AST))
        counts["interface"] += 1

    # Enums
    for m in ENUM_RE.finditer(text):
        name = m.group("name")
        eid = _ts_id("enum", ns_dotted, name)
        bundle.add(TsEnum(
            id=eid, name=name, source_layer=layer,
            namespace=ns_dotted, is_jsapi=False,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(HasClass(from_id=ns_id, to_id=eid,
                            derived_by=Provenance.AST))
        counts["enum"] += 1

    # Constants (object-shaped)
    for m in CONST_OBJ_RE.finditer(text):
        name = m.group("name")
        cnid = _ts_id("const", ns_dotted, name)
        bundle.add(TsConstant(
            id=cnid, name=name, source_layer=layer,
            namespace=ns_dotted, is_jsapi=False,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(HasClass(from_id=ns_id, to_id=cnid,
                            derived_by=Provenance.AST))
        counts["constant"] += 1

    # Top-level exported functions (not class methods)
    for m in FUNCTION_RE.finditer(text):
        name = m.group("name")
        if name in {"if", "while", "for", "switch"}:
            continue
        fid = _ts_id("fn", ns_dotted, name)
        params = m.group("params") or ""
        param_count = 0 if not params.strip() else (params.count(",") + 1)
        bundle.add(TsFunction(
            id=fid, name=name, source_layer=layer,
            namespace=ns_dotted,
            return_type=(m.group("ret") or "").strip() or None,
            parameter_count=param_count,
            is_async=bool(m.group("async")),
            is_jsapi=False,
            paths=[c.file_ref(file_path, line_start=text[:m.start()].count("\n") + 1)],
        ))
        bundle.add(HasFunction(from_id=ns_id, to_id=fid,
                               derived_by=Provenance.AST))
        counts["function"] += 1

    return counts


def _wire_extends(text: str, classes_index: dict[str, str],
                   ns_dotted: str, bundle,
                   jsapi_classes: dict[str, str] | None = None,
                   jsapi_interfaces: dict[str, str] | None = None) -> None:
    """Second pass once everything is indexed.

    Resolves the parent class via, in order:
      1. `<ns>.<Name>` in the same package/library (intra-namespace)
      2. unprefixed `<Name>` in any indexed namespace (in-package)
      3. `DG.<Name>` aliases against the canonical js-api index — strip
         `DG.` prefix, look up in `jsapi_classes`. This is what makes
         "all subclasses of DG.GridCellRenderer" answerable.
    """
    for m in CLASS_RE.finditer(text):
        name = m.group("name")
        cid = classes_index.get(f"{ns_dotted}.{name}")
        if not cid: continue
        ext = (m.group("extends") or "").split("<", 1)[0].strip()
        if ext:
            ext_clean = ext.split(",")[0].strip()
            base = ext_clean.removeprefix("DG.").removeprefix("dg.")
            tgt = (
                classes_index.get(f"{ns_dotted}.{ext_clean}")
                or classes_index.get(ext_clean)
                or (jsapi_classes or {}).get(base)
            )
            if tgt:
                bundle.add(ExtendsClass(from_id=cid, to_id=tgt,
                                        derived_by=Provenance.AST))
        # implements I, J, K  — emit ImplementsInterface for each
        impls_raw = (m.group("implements") or "").strip()
        if impls_raw:
            for raw_iface in impls_raw.split(","):
                iname = raw_iface.split("<", 1)[0].strip()
                if not iname:
                    continue
                base = iname.removeprefix("DG.").removeprefix("dg.")
                tgt = (
                    classes_index.get(f"{ns_dotted}.{iname}")
                    or classes_index.get(iname)
                    or (jsapi_interfaces or {}).get(base)
                )
                if tgt:
                    bundle.add(ImplementsInterface(
                        from_id=cid, to_id=tgt,
                        derived_by=Provenance.AST,
                    ))


def _load_jsapi_index() -> tuple[dict[str, str], dict[str, str]]:
    """{name: id} maps for js-api TsClass / TsInterface, used to resolve
    `DG.<X>` parents from plugin TS to the canonical js-api node."""
    import json
    nodes_dir = c.REPO_ROOT / ".kg" / "data" / "nodes"
    classes: dict[str, str] = {}
    ifaces: dict[str, str] = {}
    for kind, dest in [("TsClass.jsonl", classes), ("TsInterface.jsonl", ifaces)]:
        p = nodes_dir / kind
        if not p.is_file():
            continue
        for line in p.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not row.get("is_jsapi"):
                continue
            n = row.get("name"); rid = row.get("id")
            if n and rid:
                dest.setdefault(n, rid)
    return classes, ifaces


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    classes_index: dict[str, str] = {}
    totals = {"class": 0, "method": 0, "interface": 0, "enum": 0,
              "constant": 0, "function": 0}
    jsapi_classes, jsapi_ifaces = _load_jsapi_index()

    pkgs_root = repo_root / "packages"
    if pkgs_root.is_dir():
        for pkg_dir in sorted(pkgs_root.iterdir()):
            if not pkg_dir.is_dir():
                continue
            if package_filter and pkg_dir.name not in package_filter:
                continue
            ns_dotted = pkg_dir.name           # 'Chem'
            ns_id = c.pkg_id(pkg_dir.name)
            src = pkg_dir / "src"
            if not src.is_dir():
                continue
            for ts in src.rglob("*.ts"):
                if ts.name.endswith(".g.ts") or "node_modules" in ts.parts:
                    continue
                try:
                    text = ts.read_text(encoding="utf-8", errors="replace")
                except OSError:
                    continue
                for k, n in _process_file(text, ts, ns_id, ns_dotted,
                                           "pkg", bundle, classes_index).items():
                    totals[k] += n
            # Second pass for extends
            for ts in src.rglob("*.ts"):
                if ts.name.endswith(".g.ts") or "node_modules" in ts.parts:
                    continue
                try:
                    text = ts.read_text(encoding="utf-8", errors="replace")
                except OSError:
                    continue
                _wire_extends(text, classes_index, ns_dotted, bundle,
                              jsapi_classes, jsapi_ifaces)

    libs_root = repo_root / "libraries"
    if libs_root.is_dir():
        for lib_dir in sorted(libs_root.iterdir()):
            if not lib_dir.is_dir():
                continue
            ns_dotted = lib_dir.name
            ns_id = c.lib_id(lib_dir.name)
            src = lib_dir / "src"
            if not src.is_dir():
                continue
            for ts in src.rglob("*.ts"):
                if ts.name.endswith(".g.ts") or "node_modules" in ts.parts:
                    continue
                try:
                    text = ts.read_text(encoding="utf-8", errors="replace")
                except OSError:
                    continue
                for k, n in _process_file(text, ts, ns_id, ns_dotted,
                                           "lib", bundle, classes_index).items():
                    totals[k] += n
            for ts in src.rglob("*.ts"):
                if ts.name.endswith(".g.ts") or "node_modules" in ts.parts:
                    continue
                try:
                    text = ts.read_text(encoding="utf-8", errors="replace")
                except OSError:
                    continue
                _wire_extends(text, classes_index, ns_dotted, bundle,
                              jsapi_classes, jsapi_ifaces)

    print(f"[{EXTRACTOR_NAME}] {totals['class']} classes, {totals['method']} methods, "
          f"{totals['interface']} interfaces, {totals['enum']} enums, "
          f"{totals['constant']} constants, {totals['function']} functions")
