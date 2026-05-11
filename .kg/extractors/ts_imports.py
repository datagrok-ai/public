"""File-level TypeScript import edges.

For every TS/JS file under `packages/*/src/` and `libraries/*/src/`, parse
every ES module `import` statement (including multi-line and type-only forms)
plus dynamic `import('spec')` calls, and emit one `Imports(File → target)`
edge per (source_file, resolved_target) pair.

Resolution rules:
  - relative `./foo` / `../foo`     → File node in the same package/library
                                      (resolved with TS extension fallbacks)
  - `@datagrok-libraries/<lib>/...` → LibraryModule node
  - `@datagrok-libraries/<lib>`     → LibraryModule for the lib's main entry
  - `@datagrok/<pkg>/<rel>`         → File node in the target package's src/
  - `datagrok-api/{dg,grok,ui}`     → JsApiNamespace node
  - bare specifiers (`react`, `wu`) → skipped (external npm)

Symbols across multiple statements (rare) are unioned. Named, default,
namespace, and type-only flags are recorded for query-time filtering.

This is the only extractor that uses a stateful tokenizer rather than a
single regex — multi-line `import { a, b, c, ... } from '...'` and
brace-balanced `import { Foo, Bar } from '...';` need it.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.relations import Imports

from . import _common as c

EXTRACTOR_NAME = "ts_imports"

# Static `import ... from '...'` — multiline-friendly.  Captures everything
# between `import` and `from '...'` so we can post-process it for symbols.
STATIC_IMPORT_RE = re.compile(
    r"""\bimport\s+
        (?P<head>[\s\S]*?)
        \bfrom\s+
        ['"](?P<spec>[^'"]+)['"]
    """,
    re.VERBOSE,
)
# `import 'side-effect'`
SIDE_EFFECT_RE = re.compile(
    r"""\bimport\s+['"](?P<spec>[^'"]+)['"](?:\s*;)?\s*$""",
    re.MULTILINE,
)
# `await import('...')`, `import('...')`, `dynamic.then(import('...'))`
DYNAMIC_RE = re.compile(
    r"""\bimport\s*\(\s*['"](?P<spec>[^'"]+)['"]"""
)

# Strip out string literals and comments before scanning, so `'import x...'`
# and `// import` are ignored. Returns sanitized text + a function to recover
# original byte offsets if we cared (we don't — line numbers from the original).
_STRING_RE = re.compile(r"""(['"`])(?:[^'"`\\]|\\.|\$\{[^{}]*\})*?\1|/\*[\s\S]*?\*/|//[^\n]*""")

TS_EXTS = (".ts", ".tsx", ".d.ts", ".js", ".jsx", ".mjs", ".cjs")


def _strip_strings_and_comments(text: str) -> str:
    """Replace strings, template literals, and comments with same-length
    blanks so regex offsets still line up with the original text."""
    def repl(m):
        s = m.group(0)
        # Preserve newlines so line counts stay accurate
        return "".join(ch if ch == "\n" else " " for ch in s)
    # Note: this is fine because import specifiers are themselves *strings*,
    # and we run our import regex AGAINST THE ORIGINAL TEXT, not the stripped
    # one. The stripped version is only used as a sanity check for "is this
    # 'import' actually code or a comment?" — handled separately below.
    return _STRING_RE.sub(repl, text)


def _is_real_import(text: str, sanitized: str, match_start: int) -> bool:
    """Reject matches inside comments or strings by checking the sanitized
    copy at the same offset."""
    return sanitized[match_start:match_start + 6] == "import"


def _parse_head(head: str) -> dict:
    """Return {symbols: set[str], is_default: bool, is_namespace: bool,
    is_type_only: bool} from the bit between `import` and `from`."""
    info = {"symbols": set(), "is_default": False,
            "is_namespace": False, "is_type_only": False}
    h = head.strip()
    if h.startswith("type "):
        info["is_type_only"] = True
        h = h[5:].lstrip()
    # `* as X`
    if h.startswith("*"):
        info["is_namespace"] = True
        return info
    # Default + named: `Foo, { a, b }` or `Foo`
    # Pull out the named part if present
    brace_open = h.find("{")
    brace_close = h.find("}", brace_open + 1) if brace_open >= 0 else -1
    default_part = h[:brace_open] if brace_open >= 0 else h
    default_part = default_part.strip().rstrip(",").strip()
    if default_part and not default_part.startswith("*"):
        # `Default` or `Default, ...` — may be empty if only `{ a }`
        # First identifier in default_part is the default
        m = re.match(r"^([A-Za-z_$][\w$]*)", default_part)
        if m:
            info["is_default"] = True
            info["symbols"].add(m.group(1))
    if brace_open >= 0 and brace_close > brace_open:
        named = h[brace_open + 1: brace_close]
        for raw in named.split(","):
            sym = raw.strip().split(" as ")[0].strip()
            sym = sym.lstrip().rstrip().lstrip("*").strip()
            if sym.startswith("type "):
                sym = sym[5:].strip()
            if sym:
                info["symbols"].add(sym)
    return info


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------

def _resolve_relative(spec: str, importer: Path, file_index: set[str]) -> str | None:
    """Resolve `./foo` against the importer's directory, with TS extension
    fallbacks (.ts, .tsx, .d.ts, .js, /index.ts, /index.tsx, /index.js).
    Returns the repo-rel POSIX path or None."""
    base = (importer.parent / spec).resolve()
    candidates = [base.with_suffix(ext) for ext in TS_EXTS]
    candidates += [base / f"index{ext}" for ext in TS_EXTS]
    candidates.append(base)   # exact path (already had extension)
    for cand in candidates:
        try:
            rel = c.repo_rel(cand)
        except Exception:
            continue
        if rel in file_index:
            return rel
    return None


def _resolve_lib_module(spec: str, lib_index: dict[str, dict[str, str]]) -> str | None:
    """`@datagrok-libraries/lib/src/path` → lib_mod:lib:src/path.ts (best-match).
    `@datagrok-libraries/lib`            → lib's main entry."""
    if not spec.startswith("@datagrok-libraries/"):
        return None
    tail = spec[len("@datagrok-libraries/"):]
    parts = tail.split("/", 1)
    lib = parts[0]
    if lib not in lib_index:
        return None
    rest = parts[1] if len(parts) > 1 else None
    modules = lib_index[lib]
    if rest is None:
        # main entry — pick `src/index.ts` if present, else first
        for cand in ("src/index.ts", "src/main.ts", "index.ts"):
            if cand in modules:
                return modules[cand]
        # fall back to any module
        return next(iter(modules.values()), None)
    # Try rest with TS extensions
    for ext in ("", ".ts", ".tsx", ".d.ts", ".js"):
        key = rest + ext
        if key in modules:
            return modules[key]
        key2 = f"{rest}/index.ts" if not rest.endswith(".ts") else rest
        if key2 in modules:
            return modules[key2]
    # Permissive: last component match
    last = rest.rsplit("/", 1)[-1]
    for path, mid in modules.items():
        if path.endswith(f"/{last}.ts") or path.endswith(f"/{last}/index.ts"):
            return mid
    return None


def _resolve_pkg_file(spec: str, pkg_to_folder: dict[str, str],
                      file_index: set[str]) -> str | None:
    """`@datagrok/pkg/path` → file:packages/<Folder>/path resolution."""
    if not spec.startswith("@datagrok/"):
        return None
    tail = spec[len("@datagrok/"):]
    parts = tail.split("/", 1)
    npm = parts[0]
    folder = pkg_to_folder.get(npm)
    if not folder:
        return None
    rest = parts[1] if len(parts) > 1 else "src/package.ts"
    base = f"packages/{folder}/{rest}"
    for cand in (base, f"{base}.ts", f"{base}/index.ts"):
        if cand in file_index:
            return f"file:{cand}"
    return None


def _resolve_jsapi(spec: str) -> str | None:
    if spec == "datagrok-api/dg":   return "jsns:DG"
    if spec == "datagrok-api/grok": return "jsns:grok"
    if spec == "datagrok-api/ui":   return "jsns:ui"
    return None


# ---------------------------------------------------------------------------
# Indexes (built from existing JSONL)
# ---------------------------------------------------------------------------

def _load_file_index() -> set[str]:
    p = c.REPO_ROOT / ".kg" / "data" / "nodes" / "File.jsonl"
    if not p.is_file():
        return set()
    out: set[str] = set()
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        rel = row.get("relative_path")
        if rel:
            out.add(rel)
    return out


def _load_lib_module_index() -> dict[str, dict[str, str]]:
    """{lib_folder: {relative_path: lib_mod_id}}"""
    p = c.REPO_ROOT / ".kg" / "data" / "nodes" / "LibraryModule.jsonl"
    out: dict[str, dict[str, str]] = defaultdict(dict)
    if not p.is_file():
        return out
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        mid = row.get("id"); lib_id = row.get("library_id") or ""
        rel = row.get("relative_path")
        if not (mid and rel):
            continue
        if not lib_id.startswith("lib:"):
            continue
        lib = lib_id.removeprefix("lib:")
        out[lib][rel] = mid
    return out


def _load_pkg_npm_lookup() -> dict[str, str]:
    """`{npm_name_after_slash: folder}` from Package extras.npm_name."""
    p = c.REPO_ROOT / ".kg" / "data" / "nodes" / "Package.jsonl"
    out: dict[str, str] = {}
    if not p.is_file():
        return out
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        folder = row.get("name")
        extras = row.get("extras") or {}
        npm = extras.get("npm_name") or ""
        if folder and isinstance(npm, str) and npm.startswith("@datagrok/"):
            out[npm.split("/", 1)[1]] = folder
    return out


# ---------------------------------------------------------------------------
# Per-file scan
# ---------------------------------------------------------------------------

def _scan_file(ts_path: Path, source_fid: str, file_index: set[str],
               lib_index: dict[str, dict[str, str]],
               pkg_lookup: dict[str, str],
               agg: dict) -> int:
    try:
        text = ts_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return 0
    if "import" not in text:
        return 0
    sanitized = _strip_strings_and_comments(text)
    n = 0

    # Static `import ... from '...'`
    for m in STATIC_IMPORT_RE.finditer(text):
        if not _is_real_import(text, sanitized, m.start()):
            continue
        spec = m.group("spec")
        info = _parse_head(m.group("head") or "")
        target = _resolve(spec, ts_path, file_index, lib_index, pkg_lookup)
        if not target:
            continue
        _add_edge(agg, source_fid, target, info, dynamic=False)
        n += 1

    # `import 'side-effect'` (no `from`)
    for m in SIDE_EFFECT_RE.finditer(text):
        if not _is_real_import(text, sanitized, m.start()):
            continue
        spec = m.group("spec")
        target = _resolve(spec, ts_path, file_index, lib_index, pkg_lookup)
        if not target:
            continue
        _add_edge(agg, source_fid, target, {"symbols": set(),
            "is_default": False, "is_namespace": False, "is_type_only": False},
            dynamic=False)
        n += 1

    # `await import('spec')` / `import('spec')`
    for m in DYNAMIC_RE.finditer(text):
        if not _is_real_import(text, sanitized, m.start()):
            continue
        spec = m.group("spec")
        target = _resolve(spec, ts_path, file_index, lib_index, pkg_lookup)
        if not target:
            continue
        _add_edge(agg, source_fid, target, {"symbols": set(),
            "is_default": False, "is_namespace": False, "is_type_only": False},
            dynamic=True)
        n += 1

    return n


def _resolve(spec: str, importer: Path, file_index: set[str],
             lib_index: dict[str, dict[str, str]],
             pkg_lookup: dict[str, str]) -> str | None:
    if spec.startswith("./") or spec.startswith("../"):
        rel = _resolve_relative(spec, importer, file_index)
        return f"file:{rel}" if rel else None
    if spec.startswith("@datagrok-libraries/"):
        return _resolve_lib_module(spec, lib_index)
    if spec.startswith("@datagrok/"):
        return _resolve_pkg_file(spec, pkg_lookup, file_index)
    if spec.startswith("datagrok-api/"):
        return _resolve_jsapi(spec)
    return None    # external npm or unknown


def _add_edge(agg: dict, src: str, tgt: str, info: dict, *, dynamic: bool) -> None:
    key = (src, tgt)
    cur = agg.get(key)
    if cur is None:
        cur = {
            "symbols": set(),
            "is_default": False, "is_namespace": False,
            "is_type_only": True, "is_dynamic": False, "count": 0,
        }
        agg[key] = cur
    cur["symbols"].update(info.get("symbols") or ())
    cur["is_default"] = cur["is_default"] or info.get("is_default", False)
    cur["is_namespace"] = cur["is_namespace"] or info.get("is_namespace", False)
    # type_only only stays true if EVERY statement is type-only
    if not info.get("is_type_only", False):
        cur["is_type_only"] = False
    if dynamic:
        cur["is_dynamic"] = True
    cur["count"] += 1


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    file_index = _load_file_index()
    lib_index = _load_lib_module_index()
    pkg_lookup = _load_pkg_npm_lookup()
    if not file_index:
        print(f"[{EXTRACTOR_NAME}] no File.jsonl yet — run files extractor first")
        return
    agg: dict = {}

    n_files = n_skipped = 0
    pkgs = repo_root / "packages"
    if pkgs.is_dir():
        for pkg_dir in sorted(pkgs.iterdir()):
            if not pkg_dir.is_dir():
                continue
            if package_filter and pkg_dir.name not in package_filter:
                continue
            src = pkg_dir / "src"
            if not src.is_dir():
                continue
            for ts in src.rglob("*"):
                if not ts.is_file() or ts.suffix.lower() not in TS_EXTS:
                    continue
                if "node_modules" in ts.parts:
                    continue
                rel = c.repo_rel(ts)
                if rel not in file_index:
                    n_skipped += 1
                    continue
                fid = f"file:{rel}"
                if _scan_file(ts, fid, file_index, lib_index, pkg_lookup, agg):
                    n_files += 1

    libs = repo_root / "libraries"
    if libs.is_dir():
        for lib_dir in sorted(libs.iterdir()):
            if not lib_dir.is_dir():
                continue
            src = lib_dir / "src"
            if not src.is_dir():
                continue
            for ts in src.rglob("*"):
                if not ts.is_file() or ts.suffix.lower() not in TS_EXTS:
                    continue
                if "node_modules" in ts.parts:
                    continue
                rel = c.repo_rel(ts)
                if rel not in file_index:
                    n_skipped += 1
                    continue
                fid = f"file:{rel}"
                if _scan_file(ts, fid, file_index, lib_index, pkg_lookup, agg):
                    n_files += 1

    n_edges = 0
    by_kind = defaultdict(int)
    for (src_fid, tgt), info in agg.items():
        # Classify target for stats
        if tgt.startswith("file:"): by_kind["file"] += 1
        elif tgt.startswith("lib_mod:"): by_kind["lib_mod"] += 1
        elif tgt.startswith("jsns:"): by_kind["jsns"] += 1
        bundle.add(Imports(
            from_id=src_fid, to_id=tgt,
            imported_symbols=sorted(info["symbols"])[:30],
            is_type_only=info["is_type_only"],
            is_default=info["is_default"],
            is_namespace=info["is_namespace"],
            is_dynamic=info["is_dynamic"],
            import_count=info["count"],
            derived_by=Provenance.AST,
        ))
        n_edges += 1

    print(f"[{EXTRACTOR_NAME}] scanned {n_files} files ({n_skipped} skipped - not in File index), "
          f"{n_edges} Imports edges (file->file: {by_kind['file']}, file->lib_mod: {by_kind['lib_mod']}, "
          f"file->jsns: {by_kind['jsns']})")
