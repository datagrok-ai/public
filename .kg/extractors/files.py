"""Promote source files to first-class graph nodes.

Two passes:

1. **Walk the filesystem.** For every file under the modeled roots
   (`packages/*/{src,scripts,queries,connections,environments,
   dockerfiles,css,detectors.*,package*.{ts,js},README.md,CHANGELOG.md,
   META.md,CLAUDE.md,*.json}` plus `libraries/*/src/`), emit a
   `File` node + a `ContainsFile` edge from the owning Package/Library.
   Compute language, size, line count, and `is_test` / `is_generated`
   heuristics so gap queries can split impl from tests cheaply.

2. **Index existing nodes.** For every entity's `paths[]`, emit a
   `DefinedIn` edge to the corresponding `File` node. This creates the
   `Function/Script/Query → File` traversal that `feature_files.py`
   needs to derive `IsImplementedIn` / `IsTestedIn`.

We deliberately don't include `node_modules`, `dist`, `build`, `.git`,
or any binary blob (>2 MB).
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from schema.base import Provenance, SourceLayer
from schema.entities import File
from schema.relations import ContainsFile, DefinedByFunction, DefinedIn

from . import _common as c

EXTRACTOR_NAME = "files"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Per-package directories we scan recursively
PACKAGE_SUBDIRS = ("src", "scripts", "queries", "connections", "environments",
                   "dockerfiles", "css", "tests", "test")
# Per-package single files we always include if present
PACKAGE_FILES = ("package.json", "README.md", "CHANGELOG.md", "META.md",
                 "CLAUDE.md", "CONTRIB.md", "CREDITS.md", "detectors.js",
                 "detectors.ts", "package.png", "webpack.config.js")

# Per-library: same idea, narrower
LIBRARY_SUBDIRS = ("src", "tests", "test")
LIBRARY_FILES = ("package.json", "README.md", "CHANGELOG.md", "CLAUDE.md")

EXCLUDE_DIR_NAMES = {"node_modules", "dist", "build", ".cache", ".git",
                     "__pycache__", ".idea", ".vscode"}

MAX_BYTES = 2 * 1024 * 1024     # skip blobs > 2 MB (rdkit wasm etc.)

LANG_BY_EXT: dict[str, str] = {
    ".ts": "ts", ".tsx": "ts", ".js": "js", ".mjs": "js", ".jsx": "js",
    ".py": "python", ".r": "r", ".R": "r", ".jl": "julia", ".m": "octave",
    ".sql": "sql", ".md": "markdown", ".mdx": "markdown",
    ".json": "json", ".yaml": "yaml", ".yml": "yaml",
    ".css": "css", ".scss": "css", ".less": "css",
    ".html": "html", ".dart": "dart", ".java": "java",
    ".sh": "shell", ".cmd": "batch", ".ps1": "powershell",
    ".dockerfile": "dockerfile",
    ".png": "image", ".jpg": "image", ".jpeg": "image", ".gif": "image",
    ".svg": "image", ".ico": "image", ".webp": "image",
    ".wasm": "wasm",
    ".csv": "data", ".tsv": "data",
}


def _id_for(rel_path: str) -> str:
    return f"file:{rel_path}"


# ---------------------------------------------------------------------------
# Heuristics
# ---------------------------------------------------------------------------

_TEST_PATH_RE = re.compile(r"(^|/)(tests?|playwright-public)/", re.IGNORECASE)
_TEST_NAME_RE = re.compile(r"(^|[-._])(test|spec)([-._]|s?\.(?:t|j|p)y?$|$)", re.IGNORECASE)
_GEN_NAME_RE = re.compile(r"\.(?:g|api\.g|xcmd\.g|dsc\.g)\.(?:ts|js|dart)$|^package-api\.ts$")


def _is_test(rel: str) -> bool:
    if _TEST_PATH_RE.search(rel):
        return True
    name = rel.rsplit("/", 1)[-1]
    if _TEST_NAME_RE.search(name):
        return True
    if name in {"package-test.ts", "package-test.js"}:
        return True
    return False


def _is_generated(rel: str) -> bool:
    name = rel.rsplit("/", 1)[-1]
    return bool(_GEN_NAME_RE.search(name))


def _language_for(path: Path) -> str:
    if path.name.lower() == "dockerfile":
        return "dockerfile"
    return LANG_BY_EXT.get(path.suffix.lower(), "other")


def _line_count(path: Path) -> int:
    if path.suffix.lower() in {".png", ".jpg", ".jpeg", ".gif", ".ico",
                                ".webp", ".wasm", ".woff", ".woff2", ".ttf"}:
        return 0
    try:
        with path.open("rb") as f:
            return sum(1 for _ in f)
    except OSError:
        return 0


# ---------------------------------------------------------------------------
# Walk one root
# ---------------------------------------------------------------------------

def _iter_files(root: Path, subdirs: Iterable[str], single_files: Iterable[str]):
    """Yield `Path` objects under `root` according to the directory and
    file allow-lists. Skips excluded dir names."""
    for sub in subdirs:
        d = root / sub
        if not d.is_dir():
            continue
        for f in d.rglob("*"):
            if not f.is_file():
                continue
            if any(part in EXCLUDE_DIR_NAMES for part in f.parts):
                continue
            yield f
    for fname in single_files:
        f = root / fname
        if f.is_file():
            yield f


def _emit_file(path: Path, owner_pkg: str | None, owner_lib: str | None,
               bundle, seen: set[str]) -> str | None:
    """Add a File node + ContainsFile edge. Returns the file id."""
    rel = c.repo_rel(path)
    fid = _id_for(rel)
    if fid in seen:
        return fid
    seen.add(fid)

    try:
        size = path.stat().st_size
    except OSError:
        return None
    if size > MAX_BYTES:
        return None
    lang = _language_for(path)
    lines = _line_count(path) if lang not in {"image", "wasm"} else 0
    is_test = _is_test(rel)
    is_gen = _is_generated(rel)
    is_lib = owner_lib is not None
    layer = SourceLayer.LIBRARIES if is_lib else SourceLayer.PLUGINS

    node = File(
        id=fid, name=path.name,
        source_layer=layer,
        package_id=c.pkg_id(owner_pkg) if owner_pkg else None,
        library_id=c.lib_id(owner_lib) if owner_lib else None,
        language=lang, size_bytes=size, line_count=lines,
        is_test=is_test, is_generated=is_gen,
        relative_path=rel,
        paths=[c.file_ref(path, role="self")],
    )
    bundle.add(node)
    if owner_pkg:
        bundle.add(ContainsFile(from_id=c.pkg_id(owner_pkg), to_id=fid,
                                derived_by=Provenance.FILESYSTEM))
    elif owner_lib:
        bundle.add(ContainsFile(from_id=c.lib_id(owner_lib), to_id=fid,
                                derived_by=Provenance.FILESYSTEM))
    return fid


# ---------------------------------------------------------------------------
# Pass 2: emit DefinedIn from existing entities
# ---------------------------------------------------------------------------

# Which entity kinds carry meaningful `paths[0]` to a definition file.
# TS code kinds (TsClass/TsMethod/TsFunction/TsInterface/TsEnum/TsConstant)
# were added 2026-05 — every TS extractor stamps `paths[0]` with the source
# location, so the same DefinedIn pass attributes them all.
_DEFINEDIN_KINDS = {
    "RegisteredFunction", "Script", "DataQuery", "DataConnection",
    "Tutorial", "TutorialTrack", "DocPage", "ChangelogEntry",
    "TsClass", "TsMethod", "TsFunction", "TsInterface", "TsEnum", "TsConstant",
    "PackageTest",
}


def _emit_defined_in(bundle, seen_files: set[str]) -> int:
    """Read every node JSONL and emit DefinedIn for the first `definition`-
    role path on each entity. We deliberately restrict to definitional kinds
    so we don't spam edges from every leaf node."""
    nodes_dir = c.REPO_ROOT / ".kg" / "data" / "nodes"
    if not nodes_dir.is_dir():
        return 0
    n = 0
    for jl in sorted(nodes_dir.glob("*.jsonl")):
        kind = jl.stem
        if kind not in _DEFINEDIN_KINDS:
            continue
        for line in jl.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            paths = row.get("paths") or []
            # Pick the definition path (or first path otherwise)
            chosen = next((p for p in paths if (p.get("role") or "") == "definition"), None)
            if chosen is None and paths:
                chosen = paths[0]
            if not chosen or not chosen.get("path"):
                continue
            fid = _id_for(chosen["path"])
            if fid not in seen_files:
                continue       # path points outside our tracked tree
            bundle.add(DefinedIn(
                from_id=row["id"], to_id=fid,
                line_start=chosen.get("line_start"),
                line_end=chosen.get("line_end"),
                role=chosen.get("role"),
                derived_by=Provenance.FILESYSTEM,
            ))
            n += 1
    return n


def _emit_defined_by_function(bundle) -> int:
    """Bridge `RegisteredFunction` → `TsMethod`/`TsFunction` of its impl.

    The wrapper emitted by `grok api` lives in `package.g.ts` and has
    `name` matching the source method name. The actual method lives in
    `PackageFunctions.<methodName>` of `package.ts`, modeled as
    `TsMethod` with id `ts:method:<Pkg>.PackageFunctions.<name>`.
    Falls back to a top-level `TsFunction` of the same name when the
    package uses the older `//name:`-comment style.

    Pure JOIN over already-extracted JSONL — no parsing.
    """
    nodes_dir = c.REPO_ROOT / ".kg" / "data" / "nodes"

    # Index TsMethod / TsFunction by (package_lc, name_lc)
    method_idx: dict[tuple[str, str], list[str]] = {}
    function_idx: dict[tuple[str, str], list[str]] = {}

    def _ns_to_pkg(ns: str | None) -> str | None:
        # `Chem` → `Chem` (TsMethod.namespace is the package folder for plugin TS)
        if not ns: return None
        return ns

    m_path = nodes_dir / "TsMethod.jsonl"
    if m_path.is_file():
        for line in m_path.read_text(encoding="utf-8").splitlines():
            if not line.strip(): continue
            try: row = json.loads(line)
            except json.JSONDecodeError: continue
            mid = row.get("id"); mname = row.get("name")
            if not (mid and mname): continue
            # mid like ts:method:<ns>.<Class>.<method>
            parts = (mid or "").split(":", 2)
            if len(parts) < 3: continue
            tail = parts[2]   # `<ns>.<Class>.<method>`
            chunks = tail.rsplit(".", 2)
            if len(chunks) != 3: continue
            ns, cls, mn = chunks
            # Prefer methods on `PackageFunctions` (the canonical class for
            # decorator-style plugin packages) but accept any class.
            pkg = _ns_to_pkg(ns)
            if not pkg: continue
            key = (pkg.lower(), mn.lower())
            method_idx.setdefault(key, []).append((mid, cls))

    f_path = nodes_dir / "TsFunction.jsonl"
    if f_path.is_file():
        for line in f_path.read_text(encoding="utf-8").splitlines():
            if not line.strip(): continue
            try: row = json.loads(line)
            except json.JSONDecodeError: continue
            fid = row.get("id"); fname = row.get("name")
            ns = row.get("namespace")
            if not (fid and fname and ns): continue
            function_idx.setdefault((ns.lower(), fname.lower()), []).append(fid)

    n = 0
    rf_path = nodes_dir / "RegisteredFunction.jsonl"
    if not rf_path.is_file():
        return 0
    for line in rf_path.read_text(encoding="utf-8").splitlines():
        if not line.strip(): continue
        try: row = json.loads(line)
        except json.JSONDecodeError: continue
        rid = row.get("id"); rname = row.get("name") or ""
        pkg_id = row.get("package_id") or ""
        if not (rid and rname and pkg_id.startswith("pkg:")): continue
        pkg = pkg_id.removeprefix("pkg:")
        key = (pkg.lower(), rname.lower())
        # Prefer `PackageFunctions.<rname>`
        candidates = method_idx.get(key) or []
        chosen = None
        for mid, cls in candidates:
            if cls == "PackageFunctions":
                chosen = mid; break
        if not chosen and candidates:
            chosen = candidates[0][0]
        if not chosen:
            funcs = function_idx.get(key) or []
            if funcs:
                chosen = funcs[0]
        if chosen:
            bundle.add(DefinedByFunction(
                from_id=rid, to_id=chosen,
                derived_by=Provenance.AST, confidence=0.95,
            ))
            n += 1
    return n


# ---------------------------------------------------------------------------
# Entry
# ---------------------------------------------------------------------------

def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    seen: set[str] = set()
    n_files_pkg = n_files_lib = 0

    pkgs = repo_root / "packages"
    if pkgs.is_dir():
        for pkg_dir in sorted(pkgs.iterdir()):
            if not pkg_dir.is_dir():
                continue
            if package_filter and pkg_dir.name not in package_filter:
                continue
            for f in _iter_files(pkg_dir, PACKAGE_SUBDIRS, PACKAGE_FILES):
                if _emit_file(f, pkg_dir.name, None, bundle, seen):
                    n_files_pkg += 1

    libs = repo_root / "libraries"
    if libs.is_dir():
        for lib_dir in sorted(libs.iterdir()):
            if not lib_dir.is_dir():
                continue
            for f in _iter_files(lib_dir, LIBRARY_SUBDIRS, LIBRARY_FILES):
                if _emit_file(f, None, lib_dir.name, bundle, seen):
                    n_files_lib += 1

    # Pass 2: entity → File DefinedIn
    n_defined = _emit_defined_in(bundle, seen)
    # Pass 3: RegisteredFunction → TsMethod/TsFunction bridge
    n_bridge = _emit_defined_by_function(bundle)

    print(f"[{EXTRACTOR_NAME}] {n_files_pkg + n_files_lib} files "
          f"({n_files_pkg} pkg + {n_files_lib} lib), {n_defined} DefinedIn edges, "
          f"{n_bridge} DefinedByFunction edges")
