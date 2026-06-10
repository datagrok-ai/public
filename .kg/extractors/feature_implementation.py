"""Deterministic IS_IMPLEMENTED_IN / IS_TESTED_IN extractor.

Replaces the shallow `feature_files.py` pass. For every plugin under
`packages/`, it:

1. Reads `src/package.ts` and parses its `import` statements, resolving
   each relative path (`./utils/foo`, `./widgets/bar`) to a real File node.
2. Locates each `@grok.decorators.<role>({ name: 'X', ... })`-annotated
   `static` method on `PackageFunctions`, extracts its body, and finds
   which imported symbols it references.
3. For each `RegisteredFunction` named 'X', the implementation files are
   `package.ts` + every directly-referenced file + every file transitively
   imported by those (one extra hop, in-package only).
4. Aggregates per Feature: union of impl files across all PART_OF_FEATURE
   members. Emits `IsImplementedIn` edges with `derived_by=AST`.

For tests:

5. Walks every test file under `src/tests/`, `tests/`, `package-test.ts`,
   and matches each `RegisteredFunction` by name / friendly_name / TS
   method name appearing in the file. Emits `IsTestedIn` edges.

Provenance: deterministic (AST-style), so `derived_by=AST`. The shallow
`package.g.ts` edges from `feature_files.py` are kept (still true) but
this extractor's edges are the load-bearing ones for "where is X
implemented". The downstream LLM enricher refines the residual cases.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.relations import IsImplementedIn, IsTestedIn

from . import _common as c

EXTRACTOR_NAME = "feature_implementation"

# Recurse this many extra hops from package.ts when collecting impl files.
# 1 = direct imports of package.ts; 2 = imports-of-imports too.
_IMPORT_HOP_DEPTH = 2

# Cap files per RegisteredFunction so a single run-away decorator method
# doesn't spam thousands of edges via the second hop.
_MAX_FILES_PER_RF = 60


# ---------------------------------------------------------------------------
# Regex toolkit
# ---------------------------------------------------------------------------

# Captures: 1=symbols-block (with braces) or default-name, 2=source path
_IMPORT_RE = re.compile(
    r"""^\s*import\s+
        (?:type\s+)?
        (?:
            (?P<default>[A-Za-z_$][\w$]*)            # default
            (?:\s*,\s*\{(?P<named1>[^}]*)\})?        # default + named
          | \{(?P<named2>[^}]*)\}                    # named only
          | \*\s+as\s+(?P<star>[A-Za-z_$][\w$]*)     # * as ns
        )
        \s+from\s+
        ['"](?P<src>[^'"]+)['"]\s*;?\s*$
    """,
    re.VERBOSE,
)

# Match `@grok.decorators.<role>({...})` followed (possibly across lines)
# by the static method declaration `static (async )? <methodName>(`.
_DECORATOR_RE = re.compile(
    r"""@grok\.decorators\.([A-Za-z_]\w*)\s*\(""",
    re.VERBOSE,
)

# Identifier scan inside method bodies — JS identifiers including $.
_IDENT_RE = re.compile(r"\b([A-Za-z_$][\w$]*)\b")

# Resolves common TS extensions for import paths.
_RESOLVE_EXTS = (".ts", ".tsx", ".js", ".jsx")
_RESOLVE_INDEX = ("/index.ts", "/index.tsx", "/index.js")


# ---------------------------------------------------------------------------
# JSONL load helpers
# ---------------------------------------------------------------------------

def _read_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


# ---------------------------------------------------------------------------
# TS source parsing
# ---------------------------------------------------------------------------

def _parse_imports(src: str, *, src_file: Path, pkg_src_root: Path) -> dict[str, Path]:
    """Return {local_symbol: resolved_file_path} for every relative import.

    External imports (`datagrok-api/...`, `@datagrok-libraries/...`, plain
    npm names) are skipped — we can't add IS_IMPLEMENTED_IN edges to them
    because they don't live in the same package's File nodes. The library
    layer is handled separately.
    """
    out: dict[str, Path] = {}
    for line in src.splitlines():
        m = _IMPORT_RE.match(line)
        if not m:
            continue
        path_str = m.group("src")
        if not path_str.startswith("."):
            continue                                # external / npm
        target = _resolve_relative(path_str, src_file, pkg_src_root)
        if target is None:
            continue
        # Default
        if m.group("default"):
            out[m.group("default")] = target
        # Named (with possible `X as Y`)
        named = m.group("named1") or m.group("named2") or ""
        if named.strip():
            for sym in named.split(","):
                sym = sym.strip()
                if not sym:
                    continue
                if " as " in sym:
                    _, _, alias = sym.partition(" as ")
                    out[alias.strip()] = target
                else:
                    out[sym] = target
        # Star
        if m.group("star"):
            out[m.group("star")] = target
    return out


def _resolve_relative(import_path: str, importer: Path, pkg_src_root: Path) -> Path | None:
    """Resolve `'./foo'` / `'../bar/baz'` relative to importer; constrain to
    paths under pkg_src_root. Returns None if it lands outside the package."""
    base = (importer.parent / import_path).resolve()
    # Try exact, then with extensions, then as folder/index.
    candidates: list[Path] = [base]
    for ext in _RESOLVE_EXTS:
        candidates.append(base.with_suffix(ext))
    for idx in _RESOLVE_INDEX:
        candidates.append(Path(str(base) + idx))
    for cand in candidates:
        if not cand.is_file():
            continue
        try:
            cand.resolve().relative_to(pkg_src_root)
        except ValueError:
            return None                             # escaped the package
        return cand
    return None


def _find_legacy_top_level_functions(src: str) -> list[tuple[str, str, int, int]]:
    """Older `//name: X\\n//description: ...\\nexport (async) function fn(...)`
    style — used by ApiSamples-derivatives and a few legacy plugins. Return
    the same shape as `_find_decorator_methods`.
    """
    results: list[tuple[str, str, int, int]] = []
    lines = src.splitlines(keepends=True)
    line_offsets: list[int] = [0]
    for ln in lines:
        line_offsets.append(line_offsets[-1] + len(ln))
    i = 0
    n_lines = len(lines)
    while i < n_lines:
        ln = lines[i].lstrip()
        if not ln.startswith("//"):
            i += 1
            continue
        # Collect contiguous //... block.
        block_start = i
        block: list[str] = []
        while i < n_lines and lines[i].lstrip().startswith("//"):
            block.append(lines[i].lstrip()[2:].strip())
            i += 1
        # Look for `name:` in block.
        reg_name: str | None = None
        for b in block:
            m = re.match(r"name\s*:\s*(.+)", b)
            if m:
                reg_name = m.group(1).strip().strip('"').strip("'")
                break
        if not reg_name:
            continue
        # Anchor line should be `export (async )? function <fn>(`.
        if i >= n_lines:
            continue
        anchor = lines[i].lstrip()
        m = re.match(r"export\s+(?:async\s+)?function\s+([A-Za-z_$][\w$]*)\s*\(", anchor)
        if not m:
            continue
        method_name = m.group(1)
        # Find the opening `{` after the parameter list.
        anchor_offset = line_offsets[i]
        paren_open = src.find("(", anchor_offset)
        if paren_open < 0:
            continue
        params_end = _find_matching(src, paren_open, "(", ")")
        if params_end < 0:
            continue
        j = params_end + 1
        while j < len(src) and src[j] != "{":
            if src[j] == ";":
                j = -1
                break
            j += 1
        if j < 0 or j >= len(src):
            continue
        body_start = j
        body_end = _find_matching(src, body_start, "{", "}")
        if body_end < 0:
            continue
        results.append((reg_name, method_name, body_start, body_end))
    return results


def _find_decorator_methods(src: str) -> list[tuple[str, str, int, int]]:
    """Find every `@grok.decorators.<role>({...name: 'X'...}) static <fn>(...)`
    and return [(registered_name, method_name, body_start, body_end)] using
    character offsets (slice src[body_start:body_end] for the brace body).

    `name:` is optional — falls back to method name when missing.
    """
    results: list[tuple[str, str, int, int]] = []
    for m in _DECORATOR_RE.finditer(src):
        # Walk forward from the opening `(` after `@grok.decorators.<role>`
        # to find the matching `)`. Then look ahead for the static method.
        decorator_start = m.start()
        paren_open = m.end() - 1                    # the `(` itself
        paren_end = _find_matching(src, paren_open, "(", ")")
        if paren_end < 0:
            continue
        decorator_body = src[paren_open + 1:paren_end]

        # Optional explicit `name: '...'` — falls back to method name below.
        name_m = re.search(r"""name\s*:\s*['"]([^'"]+)['"]""", decorator_body)

        # After the decorator's closing `)` we may have whitespace/comments,
        # possibly more decorators on the same target (rare), then
        # `static [async] <methodName>(`. Walk past nested @decorator()
        # blocks via brace-aware matching, since their args may span lines.
        i = paren_end + 1
        while True:
            j = i
            while j < len(src) and src[j] in " \t\r\n":
                j += 1
            if j < len(src) - 1 and src[j] == "/" and src[j + 1] == "/":
                k = src.find("\n", j)
                i = k + 1 if k >= 0 else len(src)
                continue
            if j < len(src) - 1 and src[j] == "/" and src[j + 1] == "*":
                k = src.find("*/", j + 2)
                i = k + 2 if k >= 0 else len(src)
                continue
            if j < len(src) and src[j] == "@":
                # Skip another decorator: jump to its `(`, find its `)`.
                pop = src.find("(", j)
                if pop < 0:
                    break
                pcl = _find_matching(src, pop, "(", ")")
                if pcl < 0:
                    break
                i = pcl + 1
                continue
            i = j
            break
        sig = re.match(
            r"""static\s+(?:async\s+)?([A-Za-z_$][\w$]*)\s*\(""",
            src[i:], re.VERBOSE,
        )
        if not sig:
            continue
        method_name = sig.group(1)
        registered_name = name_m.group(1) if name_m else method_name

        # Walk to the opening `{` of the method body.
        method_sig_end = i + sig.end()              # at the `(` after method name
        # Find the `)` that ends the parameter list.
        params_end = _find_matching(src, method_sig_end - 1, "(", ")")
        if params_end < 0:
            continue
        # Skip optional return type `: Promise<void>` etc., then `{`.
        i = params_end + 1
        while i < len(src) and src[i] != "{":
            if src[i] == ";":                       # abstract / overload
                i = -1
                break
            i += 1
        if i < 0 or i >= len(src):
            continue
        body_start = i
        body_end = _find_matching(src, body_start, "{", "}")
        if body_end < 0:
            continue
        results.append((registered_name, method_name, body_start, body_end))
    return results


def _find_matching(src: str, open_idx: int, open_ch: str, close_ch: str) -> int:
    """Return the index of the close char matching the open at `open_idx`,
    accounting for strings/regex/comments. -1 if not found."""
    if src[open_idx] != open_ch:
        return -1
    depth = 0
    i = open_idx
    n = len(src)
    while i < n:
        ch = src[i]
        # Skip line comments
        if ch == "/" and i + 1 < n and src[i + 1] == "/":
            j = src.find("\n", i)
            i = j + 1 if j >= 0 else n
            continue
        # Skip block comments
        if ch == "/" and i + 1 < n and src[i + 1] == "*":
            j = src.find("*/", i + 2)
            i = j + 2 if j >= 0 else n
            continue
        # Skip string / template literal
        if ch in ("'", '"', "`"):
            quote = ch
            i += 1
            while i < n:
                if src[i] == "\\":
                    i += 2
                    continue
                if src[i] == quote:
                    i += 1
                    break
                # Template-literal interpolation: `...${ ... }...`
                if quote == "`" and src[i] == "$" and i + 1 < n and src[i + 1] == "{":
                    inner_end = _find_matching(src, i + 1, "{", "}")
                    if inner_end < 0:
                        return -1
                    i = inner_end + 1
                    continue
                i += 1
            continue
        if ch == open_ch:
            depth += 1
        elif ch == close_ch:
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1


def _identifiers_in(body: str) -> set[str]:
    """Tokenize identifiers in body, ignoring those inside strings/comments."""
    # Strip strings & comments first so we don't pick up identifiers there.
    cleaned = _strip_strings_and_comments(body)
    return set(_IDENT_RE.findall(cleaned))


def _strip_strings_and_comments(src: str) -> str:
    out: list[str] = []
    i = 0
    n = len(src)
    while i < n:
        ch = src[i]
        if ch == "/" and i + 1 < n and src[i + 1] == "/":
            j = src.find("\n", i)
            i = j if j >= 0 else n
            continue
        if ch == "/" and i + 1 < n and src[i + 1] == "*":
            j = src.find("*/", i + 2)
            i = j + 2 if j >= 0 else n
            continue
        if ch in ("'", '"', "`"):
            quote = ch
            i += 1
            while i < n:
                if src[i] == "\\":
                    i += 2
                    continue
                if src[i] == quote:
                    i += 1
                    break
                if quote == "`" and src[i] == "$" and i + 1 < n and src[i + 1] == "{":
                    inner = _find_matching(src, i + 1, "{", "}")
                    if inner < 0:
                        i = n
                        break
                    out.append(_strip_strings_and_comments(src[i + 2:inner]))
                    i = inner + 1
                    continue
                i += 1
            continue
        out.append(ch)
        i += 1
    return "".join(out)


# ---------------------------------------------------------------------------
# Per-package extraction
# ---------------------------------------------------------------------------

def _file_id(rel_path: str) -> str:
    return f"file:{rel_path}"


def _safe_read(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


def _gather_pkg_impl(pkg_dir: Path, repo_root: Path) -> dict[str, set[Path]]:
    """For one package, return {registered_name: set[impl_file_path]}.

    The set always includes `package.ts` itself for any decorated method
    found there. For each decorated method we add the files of imported
    symbols referenced in its body, then recurse one extra hop.
    """
    pkg_src = (pkg_dir / "src").resolve()
    pkg_ts = pkg_src / "package.ts"
    if not pkg_ts.is_file():
        return {}

    src_text = _safe_read(pkg_ts)
    if not src_text:
        return {}

    imports = _parse_imports(src_text, src_file=pkg_ts, pkg_src_root=pkg_src)
    methods = _find_decorator_methods(src_text) + _find_legacy_top_level_functions(src_text)
    if not methods:
        return {}

    # Cache of nested-import maps so we don't re-parse the same file.
    nested_imports_cache: dict[Path, dict[str, Path]] = {}

    def nested_files(file_path: Path) -> set[Path]:
        if file_path in nested_imports_cache:
            return set(nested_imports_cache[file_path].values())
        text = _safe_read(file_path)
        if not text:
            nested_imports_cache[file_path] = {}
            return set()
        sub = _parse_imports(text, src_file=file_path, pkg_src_root=pkg_src)
        nested_imports_cache[file_path] = sub
        return set(sub.values())

    out: dict[str, set[Path]] = {}
    for registered_name, _method_name, body_start, body_end in methods:
        body = src_text[body_start:body_end + 1]
        idents = _identifiers_in(body)
        impl: set[Path] = {pkg_ts}
        first_hop: set[Path] = set()
        for sym in idents:
            target = imports.get(sym)
            if target is not None:
                first_hop.add(target)
        impl.update(first_hop)

        if _IMPORT_HOP_DEPTH >= 2:
            for fp in list(first_hop):
                for nf in nested_files(fp):
                    impl.add(nf)
                    if len(impl) >= _MAX_FILES_PER_RF:
                        break
                if len(impl) >= _MAX_FILES_PER_RF:
                    break

        out.setdefault(registered_name, set()).update(impl)
    return out


def _gather_pkg_tests(pkg_dir: Path, registered_names: dict[str, dict],
                      pkg_src: Path) -> dict[str, set[Path]]:
    """Scan test files; return {registered_name: set[test_file]} when
    the file mentions the function's name / friendly_name / method name.

    Matching is conservative: we look for the literal string forms most
    commonly used in Datagrok tests:
      - 'PkgPrefix:RegisteredName' (used in grok.functions.call)
      - The registered name as a quoted string
      - The friendly_name as a quoted string (when present)
      - The TS method name as a bare identifier
    """
    tests_dir = (pkg_dir / "src" / "tests").resolve()
    extra_test_files = [
        (pkg_dir / "src" / "package-test.ts").resolve(),
        (pkg_dir / "tests").resolve(),
    ]
    test_files: list[Path] = []
    if tests_dir.is_dir():
        for f in tests_dir.rglob("*"):
            if f.is_file() and f.suffix.lower() in (".ts", ".tsx", ".js"):
                test_files.append(f)
    for ef in extra_test_files:
        if ef.is_file():
            test_files.append(ef)
        elif ef.is_dir():
            for f in ef.rglob("*"):
                if f.is_file() and f.suffix.lower() in (".ts", ".tsx", ".js"):
                    test_files.append(f)
    if not test_files:
        return {}

    # Pre-build patterns per registered function.
    # Keys we already know we need: name, friendly_name, method name.
    out: dict[str, set[Path]] = defaultdict(set)
    for tf in test_files:
        text = _safe_read(tf)
        if not text:
            continue
        for reg_name, info in registered_names.items():
            method = info.get("method_name")
            friendly = info.get("friendly_name")
            # Match a quoted form first — most reliable.
            patterns = [f"'{reg_name}'", f'"{reg_name}"']
            if friendly and friendly != reg_name:
                patterns.append(f"'{friendly}'")
                patterns.append(f'"{friendly}"')
            hit = any(p in text for p in patterns)
            # Method name as a bare identifier — only if it's distinctive
            # (length >= 5 and not a common word) to avoid false positives.
            if not hit and method and len(method) >= 5:
                pat = re.compile(rf"\b{re.escape(method)}\b")
                if pat.search(text):
                    hit = True
            if hit:
                out[reg_name].add(tf)
    return out


# ---------------------------------------------------------------------------
# Entry
# ---------------------------------------------------------------------------

def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    nodes_dir = c.REPO_ROOT / ".kg" / "data" / "nodes"
    edges_dir = c.REPO_ROOT / ".kg" / "data" / "edges"
    if not (nodes_dir.is_dir() and edges_dir.is_dir()):
        print(f"[{EXTRACTOR_NAME}] no data/ yet — run extractors first")
        return

    rfs = _read_jsonl(nodes_dir / "RegisteredFunction.jsonl")
    files = _read_jsonl(nodes_dir / "File.jsonl")
    pof = _read_jsonl(edges_dir / "PART_OF_FEATURE.jsonl")
    if not rfs:
        print(f"[{EXTRACTOR_NAME}] no RegisteredFunction nodes — bail")
        return

    # File set for fast membership check.
    valid_file_ids: set[str] = {f["id"] for f in files}
    file_is_test: dict[str, bool] = {f["id"]: bool(f.get("is_test")) for f in files}

    # Index RFs by package and by name.
    rfs_by_pkg: dict[str, list[dict]] = defaultdict(list)
    for rf in rfs:
        pkg_id = rf.get("package_id")
        if pkg_id:
            rfs_by_pkg[pkg_id].append(rf)

    # PART_OF_FEATURE: rf_id -> [feature_ids]
    rf_to_features: dict[str, set[str]] = defaultdict(set)
    for row in pof:
        rf_to_features[row["from_id"]].add(row["to_id"])

    # Walk every package.
    pkgs_root = repo_root / "packages"
    if not pkgs_root.is_dir():
        print(f"[{EXTRACTOR_NAME}] no packages/")
        return

    n_impl = n_test = 0
    pkgs_processed = 0
    pkgs_skipped = 0

    for pkg_dir in sorted(pkgs_root.iterdir()):
        if not pkg_dir.is_dir():
            continue
        pkg_name = pkg_dir.name
        if package_filter and pkg_name not in package_filter:
            continue
        pkg_id = c.pkg_id(pkg_name)
        rfs_for_pkg = rfs_by_pkg.get(pkg_id, [])
        if not rfs_for_pkg:
            continue

        impl_map = _gather_pkg_impl(pkg_dir, repo_root)
        if not impl_map:
            pkgs_skipped += 1
            continue
        pkgs_processed += 1

        # Per-RF info for test scanning (method name, friendly).
        # Reverse-lookup: registered_name -> (method_name, friendly_name)
        # from impl_map's keys + RF rows.
        rf_by_name = {rf.get("name"): rf for rf in rfs_for_pkg if rf.get("name")}
        # Also discover method names from package.ts.
        pkg_ts = (pkg_dir / "src" / "package.ts").resolve()
        ts_text = _safe_read(pkg_ts) or ""
        methods = []
        if ts_text:
            methods = _find_decorator_methods(ts_text) + _find_legacy_top_level_functions(ts_text)
        method_by_name = {reg: meth for reg, meth, _bs, _be in methods}

        registered_info: dict[str, dict] = {}
        for reg_name, rf in rf_by_name.items():
            registered_info[reg_name] = {
                "rf_id": rf["id"],
                "method_name": method_by_name.get(reg_name),
                "friendly_name": rf.get("friendly_name") or None,
            }

        test_map = _gather_pkg_tests(pkg_dir, registered_info, pkg_ts.parent)

        # Aggregate per Feature: union of all member impl files.
        feat_impl: dict[str, set[str]] = defaultdict(set)
        feat_test: dict[str, set[str]] = defaultdict(set)
        feat_member_count: dict[tuple[str, str], int] = defaultdict(int)

        for reg_name, info in registered_info.items():
            rf_id = info["rf_id"]
            feat_ids = rf_to_features.get(rf_id, set())
            if not feat_ids:
                continue
            # Implementation files for this RF
            for fp in impl_map.get(reg_name, ()):
                rel = c.repo_rel(fp)
                fid = _file_id(rel)
                if fid not in valid_file_ids:
                    continue
                if file_is_test.get(fid, False):
                    continue
                for feat in feat_ids:
                    feat_impl[feat].add(fid)
                    feat_member_count[(feat, fid)] += 1
            # Test files for this RF
            for tf in test_map.get(reg_name, ()):
                rel = c.repo_rel(tf)
                fid = _file_id(rel)
                if fid not in valid_file_ids:
                    continue
                if not file_is_test.get(fid, False):
                    continue
                for feat in feat_ids:
                    feat_test[feat].add(fid)
                    feat_member_count[(feat, fid)] += 1

        for feat_id, file_ids in feat_impl.items():
            for fid in file_ids:
                bundle.add(IsImplementedIn(
                    from_id=feat_id, to_id=fid,
                    member_count=feat_member_count.get((feat_id, fid), 1),
                    derived_by=Provenance.AST,
                ))
                n_impl += 1
        for feat_id, file_ids in feat_test.items():
            for fid in file_ids:
                bundle.add(IsTestedIn(
                    from_id=feat_id, to_id=fid,
                    member_count=feat_member_count.get((feat_id, fid), 1),
                    derived_by=Provenance.AST,
                ))
                n_test += 1

    print(f"[{EXTRACTOR_NAME}] {pkgs_processed} packages, "
          f"{pkgs_skipped} skipped (no decorators)")
    print(f"  {n_impl} IsImplementedIn (AST), {n_test} IsTestedIn (AST)")
