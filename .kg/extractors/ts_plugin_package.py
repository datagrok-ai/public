"""Extract Package + RegisteredFunction + PackageProperty entities from
`packages/*/`.

Inputs read:
  - `package.json` (canonical metadata)
  - `src/package.g.ts` (canonical generated annotation form — preferred)
  - `src/package.ts` (fallback if .g.ts missing)
  - `detectors.js` / `detectors.ts` at package root or `src/`

Emits:
  - Package
  - RegisteredFunction
  - PackageProperty
  - Exports (Package -> RegisteredFunction)
  - HasProperty (Package -> PackageProperty)
  - DependsOn / DependsOn(is_dev=True) (Package -> Library | Package)

Note on dependency edges: target IDs may not yet exist when this extractor
runs (the library extractor populates them later). That's fine — Kuzu's
COPY skips edges with missing endpoints, and the JSONL persists for the
next build pass. The build orchestrator runs library extraction first
when both are scheduled.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from schema.base import SourceLayer, Provenance
from schema.entities import Package, RegisteredFunction, PackageProperty
from schema.relations import Exports, HasProperty, DependsOn

from . import _common as c

EXTRACTOR_NAME = "ts_plugin_package"


def _parse_dep_target_id(npm_name: str) -> str | None:
    """Map an npm dep name to a KG entity ID, or None to skip.

    @datagrok-libraries/utils  -> lib:utils
    @datagrok/chem             -> pkg:Chem  (case fixed by lookup later)
    datagrok-api               -> None (it's the JS API; we model it elsewhere)
    """
    if npm_name.startswith("@datagrok-libraries/"):
        lib = npm_name.split("/", 1)[1]
        return c.lib_id(lib)
    if npm_name.startswith("@datagrok/"):
        # The folder name is PascalCase; npm name is kebab-case scoped.
        # We can't recover the folder here without scanning. Use the npm tail
        # as the ID and let a post-pass rewrite to the actual pkg:<Folder>.
        tail = npm_name.split("/", 1)[1]
        return f"pkg-npm:{tail}"   # marker; resolved by the build pass
    return None


def _emit_dep_edges(pkg: Package, pkg_json: dict[str, Any], deps_field: str,
                    is_dev: bool) -> list[DependsOn]:
    out: list[DependsOn] = []
    deps = pkg_json.get(deps_field) or {}
    for dep_name, semver in deps.items():
        target = _parse_dep_target_id(dep_name)
        if not target:
            continue
        is_test_harness = dep_name == "@datagrok-libraries/test"
        out.append(DependsOn(
            from_id=pkg.id, to_id=target,
            semver_range=str(semver), is_dev=is_dev,
            is_test_harness=is_test_harness,
            derived_by=Provenance.PACKAGE_JSON,
        ))
    return out


def _emit_property(pkg: Package, prop: dict[str, Any]) -> tuple[PackageProperty, HasProperty]:
    name = prop.get("name", "?")
    pid = c.prop_id(pkg.name, name)
    pe = PackageProperty(
        id=pid, name=name, source_layer=SourceLayer.PLUGINS,
        package_id=pkg.id,
        type=prop.get("propertyType"),
        choices=list(prop.get("choices") or []),
        default_value=prop.get("defaultValue"),
        extras={k: v for k, v in prop.items()
                if k not in {"name", "propertyType", "choices", "defaultValue"}},
    )
    edge = HasProperty(from_id=pkg.id, to_id=pid,
                       derived_by=Provenance.PACKAGE_JSON)
    return pe, edge


def _emit_function(pkg: Package, ann: dict[str, Any], paths: list,
                   src_path: Path, lang: str = "ts") -> tuple[RegisteredFunction, Exports]:
    name = ann.get("name") or _infer_name_from_anchor(ann)
    if not name:
        name = "<anonymous>"
    role = ann.get("meta", {}).get("role")
    if not role:
        # Infer from output type / tags
        outs = ann.get("outputs") or []
        if outs and outs[0].get("type") == "viewer":
            role = "viewer"
        elif outs and outs[0].get("type") == "widget":
            role = "widget"
        elif outs and outs[0].get("type") == "filter":
            role = "filter"
        elif "app" in ann.get("tags", []):
            role = "app"
    fid = c.func_id(pkg.name, name)
    rf = RegisteredFunction(
        id=fid, name=name, source_layer=SourceLayer.PLUGINS,
        package_id=pkg.id,
        friendly_name=ann.get("friendly_name"),
        role=role,
        tags=ann.get("tags", []),
        inputs=ann.get("inputs", []),
        outputs=ann.get("outputs", []),
        meta=ann.get("meta", {}),
        help_url=ann.get("help_url"),
        language=lang,
        description=ann.get("description"),
        description_provenance=Provenance.ANNOTATION if ann.get("description") else None,
        paths=paths,
    )
    edge = Exports(from_id=pkg.id, to_id=fid, derived_by=Provenance.ANNOTATION)
    return rf, edge


_EXPORT_RE = re.compile(
    r"^\s*export\s+(?:async\s+)?(?:function|const|class|var|let)\s+([a-zA-Z_$][\w$]*)"
)


def _infer_name_from_anchor(ann: dict[str, Any]) -> str | None:
    """If the annotation block didn't carry an explicit //name, infer from
    the export declaration on the next line (TS file)."""
    return None  # placeholder; the caller can pass the anchor line


def _scan_ts_file(pkg: Package, src_path: Path) -> tuple[list[RegisteredFunction], list[Exports]]:
    text = src_path.read_text(encoding="utf-8", errors="replace")
    funcs: list[RegisteredFunction] = []
    edges: list[Exports] = []
    for start, end, block_lines, anchor in c.find_annotation_blocks(text):
        ann = c.parse_annotation_block(block_lines)
        if not ann.get("name"):
            mexp = _EXPORT_RE.match(anchor)
            if mexp:
                ann["name"] = mexp.group(1)
        if not ann.get("name") and not ann.get("meta"):
            continue
        paths = [c.file_ref(src_path, line_start=start, line_end=end + 1,
                            role="definition")]
        rf, edge = _emit_function(pkg, ann, paths, src_path,
                                  lang="js" if src_path.suffix == ".js" else "ts")
        funcs.append(rf)
        edges.append(edge)
    return funcs, edges


def _build_package(pkg_dir: Path, pkg_json: dict[str, Any]) -> Package:
    folder = pkg_dir.name
    return Package(
        id=c.pkg_id(folder),
        name=folder,                       # folder name is the canonical identifier
        source_layer=SourceLayer.PLUGINS,
        friendly_name=pkg_json.get("friendlyName"),
        version=pkg_json.get("version"),
        category=pkg_json.get("category"),
        author=(pkg_json.get("author") or {}).get("name") if isinstance(pkg_json.get("author"), dict) else pkg_json.get("author"),
        description=pkg_json.get("description"),
        description_provenance=Provenance.PACKAGE_JSON if pkg_json.get("description") else None,
        properties_count=len(pkg_json.get("properties") or []),
        dependencies=list((pkg_json.get("dependencies") or {}).keys()),
        dev_dependencies=list((pkg_json.get("devDependencies") or {}).keys()),
        paths=[c.file_ref(pkg_dir, role="root"),
               c.file_ref(pkg_dir / "package.json", role="manifest")],
        extras={
            "npm_name": pkg_json.get("name"),
            "fullName": pkg_json.get("fullName"),
            "canEdit": pkg_json.get("canEdit"),
            "canView": pkg_json.get("canView"),
            "browserFeatures": (pkg_json.get("meta") or {}).get("browserFeatures"),
        },
    )


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    """Run the extractor over `packages/*/`.

    `package_filter` (optional) — a list of folder names to restrict to; useful
    for incremental rebuilds and testing.
    """
    packages_dir = repo_root / "packages"
    if not packages_dir.is_dir():
        print(f"[{EXTRACTOR_NAME}] no packages dir at {packages_dir}")
        return

    n_pkg = n_func = n_prop = n_dep = 0
    for pkg_dir in sorted(packages_dir.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        pkg_json_path = pkg_dir / "package.json"
        if not pkg_json_path.is_file():
            continue

        try:
            pkg_json = json.loads(pkg_json_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            print(f"[{EXTRACTOR_NAME}] {pkg_dir.name}: bad package.json ({e})")
            continue

        pkg = _build_package(pkg_dir, pkg_json)
        bundle.add(pkg)
        n_pkg += 1

        # Properties
        for prop in pkg_json.get("properties") or []:
            pe, edge = _emit_property(pkg, prop)
            bundle.add(pe)
            bundle.add(edge)
            n_prop += 1

        # Dependency edges
        for e in _emit_dep_edges(pkg, pkg_json, "dependencies", is_dev=False):
            bundle.add(e); n_dep += 1
        for e in _emit_dep_edges(pkg, pkg_json, "devDependencies", is_dev=True):
            bundle.add(e); n_dep += 1

        # Functions: prefer package.g.ts (canonical), fall back to package.ts
        candidates = [
            pkg_dir / "src" / "package.g.ts",
            pkg_dir / "src" / "package.ts",
            pkg_dir / "src" / "package.js",     # older JS-only packages (Notebooks, TensorFlow.js)
            pkg_dir / "src" / "detectors.ts",
            pkg_dir / "src" / "detectors.js",
            pkg_dir / "detectors.ts",
            pkg_dir / "detectors.js",
        ]
        seen_funcs: set[str] = set()
        for src in candidates:
            if not src.is_file():
                continue
            funcs, edges = _scan_ts_file(pkg, src)
            for rf, edge in zip(funcs, edges):
                # `package.g.ts` is canonical; once we have a function ID,
                # don't re-emit from `package.ts`.
                if rf.id in seen_funcs:
                    continue
                seen_funcs.add(rf.id)
                bundle.add(rf)
                bundle.add(edge)
                n_func += 1

    print(f"[{EXTRACTOR_NAME}] {n_pkg} packages, {n_func} functions, "
          f"{n_prop} properties, {n_dep} dep edges")
