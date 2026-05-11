"""Extract Library + LibraryModule entities from `libraries/*`,
plus `Package → LibraryModule` import edges by scanning package TS for
`from '@datagrok-libraries/<lib>/src/<path>'`.

Emits:
  - Library
  - LibraryModule (one per `.ts` under each library's `src/`)
  - DependsOn (Library → Library)
  - ImportsFromModule (Package → LibraryModule)
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import Library, LibraryModule
from schema.relations import DependsOn, ImportsFromModule

from . import _common as c

EXTRACTOR_NAME = "library"

DEEP_IMPORT_RE = re.compile(
    r"""from\s+['"]@datagrok-libraries/([a-zA-Z0-9._-]+)/(src/[^'"]+?)['"]"""
)
NAMED_IMPORT_RE = re.compile(
    r"""import\s+(?:type\s+)?(?:\*\s+as\s+\w+|\{([^}]+)\}|\w+)\s+from""", re.MULTILINE
)
TYPE_ONLY_RE = re.compile(r"""import\s+type\s""")


def _classify_role(pkg_json: dict, src_count: int) -> str:
    """Heuristic role classification (LLM enricher can refine later)."""
    name = (pkg_json.get("name") or "").split("/")[-1]
    deps = pkg_json.get("dependencies") or {}
    has_dg_api = "datagrok-api" in deps
    if name == "test":
        return "test-harness"
    if name in {"utils", "math", "statistics"}:
        return "utility"
    if name in {"compute-utils", "compute-api", "tutorials", "cruddy", "db-explorer"}:
        return "compute-framework" if "compute" in name else "feature"
    if name in {"webcomponents", "webcomponents-vue", "gridext"}:
        return "ui-component"
    if name in {"chem-meta"}:
        return "domain-wrapper"
    if name in {"bio", "ml"}:
        return "domain-framework"
    return "feature" if has_dg_api else "domain-wrapper"


def _build_library(lib_dir: Path) -> Library | None:
    pkj_path = lib_dir / "package.json"
    if not pkj_path.is_file():
        return None
    try:
        pkj = json.loads(pkj_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    npm_name = pkj.get("name", "")
    if not npm_name.startswith("@datagrok-libraries/"):
        return None
    short = npm_name.split("/", 1)[1]
    src_files = list((lib_dir / "src").rglob("*.ts")) if (lib_dir / "src").is_dir() else []
    deps = pkj.get("dependencies") or {}
    return Library(
        id=c.lib_id(short),
        name=short,
        source_layer=SourceLayer.LIBRARIES,
        npm_name=npm_name,
        version=pkj.get("version"),
        role=_classify_role(pkj, len(src_files)),
        main=pkj.get("main"),
        is_platform_agnostic="datagrok-api" not in deps,
        ts_file_count=len(src_files),
        description=pkj.get("description"),
        description_provenance=Provenance.PACKAGE_JSON if pkj.get("description") else None,
        paths=[c.file_ref(lib_dir, role="root"),
               c.file_ref(pkj_path, role="manifest")],
        extras={"dependencies": list(deps.keys())},
    )


def _emit_lib_modules(lib: Library, lib_dir: Path, bundle) -> int:
    src = lib_dir / "src"
    if not src.is_dir():
        return 0
    n = 0
    for ts_file in sorted(src.rglob("*.ts")):
        rel = ts_file.relative_to(lib_dir).as_posix()    # 'src/typed-metrics/consts.ts'
        mid = c.lib_mod_id(lib.name, rel)
        bundle.add(LibraryModule(
            id=mid, name=rel, source_layer=SourceLayer.LIBRARIES,
            library_id=lib.id, relative_path=rel,
            paths=[c.file_ref(ts_file, role="definition")],
        ))
        n += 1
    return n


def _emit_lib_deps(lib: Library, lib_dir: Path, bundle) -> int:
    pkj = json.loads((lib_dir / "package.json").read_text(encoding="utf-8"))
    n = 0
    for field, is_dev in (("dependencies", False), ("devDependencies", True)):
        for npm, semver in (pkj.get(field) or {}).items():
            if not npm.startswith("@datagrok-libraries/"):
                continue
            target = c.lib_id(npm.split("/", 1)[1])
            bundle.add(DependsOn(
                from_id=lib.id, to_id=target,
                semver_range=str(semver), is_dev=is_dev,
                is_test_harness=npm == "@datagrok-libraries/test",
                derived_by=Provenance.PACKAGE_JSON,
            ))
            n += 1
    return n


def _scan_package_imports(pkg_dir: Path, pkg_node_id: str,
                          known_modules: dict[tuple[str, str], str],
                          bundle) -> int:
    """For each `import ... from '@datagrok-libraries/<lib>/src/<path>'`
    emit an ImportsFromModule edge. Symbol extraction is best-effort regex
    (a real TS AST pass would be more accurate)."""
    src_dirs = [pkg_dir / "src"]
    n = 0
    grouped: dict[tuple[str, str], dict] = defaultdict(lambda: {"symbols": set(), "count": 0, "type_only": True})
    for sd in src_dirs:
        if not sd.is_dir():
            continue
        for ts in sd.rglob("*.ts"):
            try:
                text = ts.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            for m in re.finditer(
                r"""import\s+(?P<head>(?:type\s+)?(?:\{[^}]+\}|\*\s+as\s+\w+|\w+))?\s*(?:,\s*\{[^}]+\})?\s*from\s+['"](?P<spec>@datagrok-libraries/(?P<lib>[a-zA-Z0-9._-]+)/src/[^'"]+?)['"]""",
                text,
            ):
                head = (m.group("head") or "").strip()
                spec = m.group("spec")
                lib_name = m.group("lib")
                rel = "src/" + spec.split("/src/", 1)[1]
                # Some imports omit the .ts extension; try both
                key = (lib_name, rel)
                key_ts = (lib_name, rel + ".ts") if not rel.endswith(".ts") else key
                tgt = known_modules.get(key) or known_modules.get(key_ts)
                if not tgt:
                    # try with /index.ts
                    tgt = known_modules.get((lib_name, rel + "/index.ts"))
                if not tgt:
                    continue
                bucket = grouped[(pkg_node_id, tgt)]
                bucket["count"] += 1
                if head.startswith("type "):
                    pass    # keep type_only True
                else:
                    bucket["type_only"] = False
                # Best-effort symbol extraction from `{ a, b as c }`
                inner = re.search(r"\{([^}]+)\}", head)
                if inner:
                    for sym in (s.strip() for s in inner.group(1).split(",")):
                        if sym:
                            bucket["symbols"].add(sym.split(" as ")[0].strip())
    for (src, dst), info in grouped.items():
        bundle.add(ImportsFromModule(
            from_id=src, to_id=dst,
            imported_symbols=sorted(info["symbols"])[:20],
            import_count=info["count"],
            is_type_only=info["type_only"],
            derived_by=Provenance.AST,
        ))
        n += 1
    return n


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    libs_dir = repo_root / "libraries"
    if not libs_dir.is_dir():
        print(f"[{EXTRACTOR_NAME}] no libraries dir")
        return
    libs: list[Library] = []
    known_modules: dict[tuple[str, str], str] = {}     # (lib_name, rel_path) -> mod_id

    for lib_dir in sorted(libs_dir.iterdir()):
        if not lib_dir.is_dir():
            continue
        lib = _build_library(lib_dir)
        if not lib:
            continue
        bundle.add(lib)
        libs.append(lib)
        n_mods = _emit_lib_modules(lib, lib_dir, bundle)
        # Index
        if (lib_dir / "src").is_dir():
            for ts in (lib_dir / "src").rglob("*.ts"):
                rel = ts.relative_to(lib_dir).as_posix()
                known_modules[(lib.name, rel)] = c.lib_mod_id(lib.name, rel)
        _emit_lib_deps(lib, lib_dir, bundle)

    # Scan package imports
    pkgs_dir = repo_root / "packages"
    n_imports = 0
    if pkgs_dir.is_dir():
        for pkg_dir in sorted(pkgs_dir.iterdir()):
            if not pkg_dir.is_dir():
                continue
            if package_filter and pkg_dir.name not in package_filter:
                continue
            n_imports += _scan_package_imports(
                pkg_dir, c.pkg_id(pkg_dir.name), known_modules, bundle,
            )
    print(f"[{EXTRACTOR_NAME}] {len(libs)} libraries, "
          f"{sum(l.ts_file_count for l in libs)} library modules, "
          f"{n_imports} import-from-module edges")
