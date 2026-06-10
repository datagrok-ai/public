"""Extract `category('Cat', () => { test('name', async () => { ... }) })` blocks
from `packages/<Pkg>/src/tests/**/*.ts` (and `package-test.ts`) as
`PackageTest` nodes.

Each `test(...)` call becomes one node; `category` and `package_id` give the
grouping. Emits `DefinedIn(PackageTest → File)` so click-through to the source
file works. Also emits `Covers(PackageTest → RegisteredFunction)` when the
test name (or category) matches a known RF's name or friendly_name — fragile
heuristic, but strictly additive (low confidence).

Doesn't try to parse the test body — that's a separate concern handled by
test-coverage tooling, not the KG.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import PackageTest
from schema.relations import Covers, DefinedIn

from . import _common as c

EXTRACTOR_NAME = "package_tests"

# `category('Cat', () => {` or `category("Cat", () => {`
CATEGORY_RE = re.compile(
    r"""\bcategory\s*\(\s*
        (?P<q>['"`])(?P<name>(?:(?!(?P=q)).)*?)(?P=q)
        \s*,
    """, re.VERBOSE,
)
# `test('name', ...)`, `test.skip('name', ...)`, `test.only('name', ...)`,
# also accept multi-line params (parameters: { ... })
TEST_RE = re.compile(
    r"""\btest(?P<modifier>\.skip|\.only)?\s*\(\s*
        (?P<q>['"`])(?P<name>(?:(?!(?P=q)).)*?)(?P=q)
    """, re.VERBOSE,
)


def _walk_blocks(text: str) -> list[tuple[int, int, str | None]]:
    """For each `category(...)` opening, find its matching `})` close so we
    can attribute nested `test(...)` calls to the right category. Brace-aware,
    handles strings and template literals."""
    out: list[tuple[int, int, str | None]] = []
    for m in CATEGORY_RE.finditer(text):
        name = m.group("name")
        # Find the start of the arrow function body `=> {`
        i = m.end()
        # Scan forward to the first `{` that opens the block (skipping the arrow)
        depth = 0
        in_q: str | None = None
        body_open = -1
        # Allow up to ~200 chars between `,` and `{` (parameters, comments)
        scan_limit = min(len(text), i + 600)
        while i < scan_limit:
            ch = text[i]
            if in_q:
                if ch == "\\" and i + 1 < len(text):
                    i += 2; continue
                if ch == in_q:
                    in_q = None
                i += 1; continue
            if ch in "'\"`":
                in_q = ch; i += 1; continue
            if ch == "{":
                body_open = i; break
            i += 1
        if body_open < 0:
            continue
        # Now scan from body_open until matching close `}`
        depth = 1
        j = body_open + 1
        while j < len(text) and depth > 0:
            ch = text[j]
            if in_q:
                if ch == "\\" and j + 1 < len(text):
                    j += 2; continue
                if ch == in_q:
                    in_q = None
                j += 1; continue
            if ch in "'\"`":
                in_q = ch; j += 1; continue
            if ch == "{": depth += 1
            elif ch == "}": depth -= 1
            j += 1
        out.append((body_open, j, name))
    return out


def _scan_file(file_path: Path, pkg_name: str, bundle,
               known_funcs: dict[tuple[str, str], str]) -> int:
    try:
        text = file_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return 0
    if "test(" not in text and "test.skip(" not in text and "test.only(" not in text:
        return 0

    rel = c.repo_rel(file_path)
    fid = f"file:{rel}"
    blocks = _walk_blocks(text)

    n = 0
    for tm in TEST_RE.finditer(text):
        # Find which category block this test sits in (innermost wins)
        pos = tm.start()
        category = None
        for start, end, name in blocks:
            if start < pos < end:
                category = name
                # don't break — pick last matching (deepest)
        test_name = tm.group("name")
        modifier = tm.group("modifier") or ""
        # ID: include file path so collisions across files don't merge
        tid = f"pkgtest:{pkg_name}:{rel}:{test_name}"
        line_no = text[:pos].count("\n") + 1
        bundle.add(PackageTest(
            id=tid, name=test_name,
            source_layer=SourceLayer.PLUGINS,
            package_id=c.pkg_id(pkg_name),
            category=category,
            test_name=test_name,
            is_skipped=(modifier == ".skip"),
            is_only=(modifier == ".only"),
            file_path=rel,
            paths=[c.file_ref(file_path, line_start=line_no, role="test")],
        ))
        bundle.add(DefinedIn(
            from_id=tid, to_id=fid,
            line_start=line_no, role="test",
            derived_by=Provenance.AST,
        ))
        # Best-effort Covers: name match (case-insensitive substring)
        tn = test_name.lower()
        cn = (category or "").lower()
        for (pkg_lc, key_lc), rf_id in known_funcs.items():
            if pkg_lc != pkg_name.lower():
                continue
            if len(key_lc) < 4:
                continue   # avoid `init`/`run`/`get` matching everything
            if key_lc in tn or key_lc in cn:
                bundle.add(Covers(
                    from_id=tid, to_id=rf_id,
                    derived_by=Provenance.AST, confidence=0.5,
                ))
                # one match per test is enough
                break
        n += 1
    return n


def _build_function_lookup() -> dict[tuple[str, str], str]:
    """Same shape as cross_package._build_function_lookup but kept local
    to avoid the dependency."""
    rf_path = c.REPO_ROOT / ".kg" / "data" / "nodes" / "RegisteredFunction.jsonl"
    out: dict[tuple[str, str], str] = {}
    if not rf_path.is_file():
        return out
    for line in rf_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        rid = row.get("id"); pkg_id = row.get("package_id") or ""
        if not rid or not pkg_id.startswith("pkg:"):
            continue
        pkg = pkg_id.removeprefix("pkg:").lower()
        for key in (row.get("name"), row.get("friendly_name")):
            if key:
                out.setdefault((pkg, key.lower()), rid)
    return out


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return
    known_funcs = _build_function_lookup()
    n_files = n_tests = 0
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        # Scan src/tests/**.ts and src/package-test.ts
        for sub in (pkg_dir / "src" / "tests", pkg_dir / "src" / "test"):
            if not sub.is_dir():
                continue
            for ts in sub.rglob("*.ts"):
                if ts.name.endswith(".g.ts") or "node_modules" in ts.parts:
                    continue
                added = _scan_file(ts, pkg_dir.name, bundle, known_funcs)
                if added:
                    n_files += 1
                    n_tests += added
        for fname in ("package-test.ts", "package-test.js"):
            f = pkg_dir / "src" / fname
            if f.is_file():
                added = _scan_file(f, pkg_dir.name, bundle, known_funcs)
                if added:
                    n_files += 1
                    n_tests += added

    print(f"[{EXTRACTOR_NAME}] {n_tests} test() blocks across {n_files} files")
