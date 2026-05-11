"""Extract cross-package code-level edges:

  - `ImportsFromPackage`: Package -> Package, from `import ... from '@datagrok/<other>'`
  - `Calls`              : Package -> RegisteredFunction, from
                            `grok.functions.eval('<OtherPkg>:<fn>')`
                            `grok.functions.call('<OtherPkg>:<fn>')`
                            `grok.functions.find('<OtherPkg>:<fn>')`

Both signals are stronger than the package.json `DEPENDS_ON` edges because
they prove the dependency is *actually used* in source. Use them to find
real coupling, dead deps, and runtime call graphs.

Function-level call resolution: when the target is `<Pkg>:<fn>` and the
matching `func:<Pkg>:<fn>` node exists, the edge resolves to that function.
Otherwise it falls back to the Package node (so we still record the call
even if the function spec is dynamic / typo / older platform name).
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.relations import Calls, ImportsFromPackage

from . import _common as c

EXTRACTOR_NAME = "cross_package"

# `import ... from '@datagrok/<other>/...'` — covers all 4 forms:
#   import { A, B } from '@datagrok/charts'
#   import type { X } from '@datagrok/charts/src/foo'
#   import * as Charts from '@datagrok/charts'
#   import Charts from '@datagrok/charts'
IMPORT_RE = re.compile(
    r"""import\s+
        (?:type\s+)?
        (?:
            (?P<braced>\{[^}]+\})            # named imports
          | (?P<star>\*\s+as\s+\w+)          # namespace
          | (?P<default>\w+)                 # default
        )?
        \s*(?:,\s*\{[^}]+\})?\s*
        from\s+['"](?P<spec>@datagrok/(?P<pkg>[a-z0-9._-]+)(?:/[^'"]*)?)['"]
    """, re.VERBOSE)

TYPE_ONLY_RE = re.compile(r"import\s+type\s")

# Many ways to invoke another package's function — all yield the same Calls edge.
#
#   grok.functions.eval('Pkg:fn', ...)
#   grok.functions.call('Pkg:fn', ...)
#   grok.functions.find('Pkg:fn')                 (some older paths)
#   _package.functions.eval('Pkg:fn', ...)
#   DG.Func.byName('Pkg:fn')
#   DG.Func.find({name:'fn', package:'Pkg'})       — separate regex
#   DG.Func.find({name:'My Function', package:'Pkg'}) — friendly_name with spaces
#   await import('@datagrok/charts')               — handled by IMPORT_RE
#
# `fn` accepts: identifier, friendly-name with spaces / dashes / dots.
# `pkg` is a single token (no whitespace).
CALL_STRING_SPEC = re.compile(
    r"""(?:
            (?:grok|_package)\.functions\.(?:eval|call|find)
          | DG\.Func\.byName
        )\s*\(\s*
        ['"]
        (?P<spec>
            (?P<pkg>[A-Za-z][\w.\-]*)
            \s*:\s*
            (?P<fn>[A-Za-z_$][\w$ \.\-]*?)
        )
        ['"]
    """, re.VERBOSE)

# DG.Func.find({name:'fn', package:'Pkg'}) — order of keys can vary; both
# required. `fn` allows friendly-name spelling (spaces, dashes, dots).
CALL_FIND_OBJ = re.compile(
    r"""DG\.Func\.find\s*\(\s*\{
        (?=[^}]*?\bname\s*:\s*['"](?P<fn>[^'"]+)['"])
        (?=[^}]*?\bpackage\s*:\s*['"](?P<pkg>[A-Za-z][\w.\-]*)['"])
        [^}]*?\}\s*\)
    """, re.VERBOSE)


def _kebab_to_folder(npm_kebab: str, folder_lookup: dict[str, str]) -> str | None:
    """`charts` -> `Charts` etc. Falls back to a placeholder ID that
    build.py rewrites once all Package nodes are loaded."""
    return folder_lookup.get(npm_kebab)


def _build_folder_lookup(repo_root: Path) -> dict[str, str]:
    """Read every `packages/*/package.json` and map npm-kebab → folder."""
    out: dict[str, str] = {}
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return out
    import json as _json
    for d in sorted(pkgs.iterdir()):
        pj = d / "package.json"
        if not pj.is_file():
            continue
        try:
            obj = _json.loads(pj.read_text(encoding="utf-8"))
        except Exception:
            continue
        name = obj.get("name", "")
        if name.startswith("@datagrok/"):
            out[name.split("/", 1)[1]] = d.name
    return out


def _build_function_lookup() -> dict[tuple[str, str], str]:
    """Read RegisteredFunction.jsonl and map (pkg_folder_lower, key_lower) → id,
    where `key` is both `name` and `friendly_name`. Lets Calls resolve when the
    code spells the function as either."""
    import json as _json
    path = c.REPO_ROOT / ".kg" / "data" / "nodes" / "RegisteredFunction.jsonl"
    out: dict[tuple[str, str], str] = {}
    if not path.is_file():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            row = _json.loads(line)
        except _json.JSONDecodeError:
            continue
        rid = row.get("id")
        pkg_id = row.get("package_id") or ""
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

    folder_lookup = _build_folder_lookup(repo_root)
    fn_lookup = _build_function_lookup()

    # (caller_pkg, target_pkg) -> {symbols, count, type_only_all}
    imports: dict[tuple[str, str], dict] = defaultdict(
        lambda: {"symbols": set(), "count": 0, "type_only_all": True, "evidence": set()})
    # (caller_pkg, target_pkg, target_fn) -> {count, evidence_files}
    calls: dict[tuple[str, str, str], dict] = defaultdict(
        lambda: {"count": 0, "evidence": set()})

    n_files = 0
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        caller = pkg_dir.name

        for src_dir in (pkg_dir / "src", pkg_dir):     # search src/ first, then root .js
            if not src_dir.is_dir():
                continue
            for f in src_dir.rglob("*.ts"):
                if "node_modules" in f.parts: continue
                _scan_file(f, caller, folder_lookup, imports, calls)
                n_files += 1
            for f in src_dir.rglob("*.js"):
                if "node_modules" in f.parts: continue
                _scan_file(f, caller, folder_lookup, imports, calls)
                n_files += 1

    # Emit ImportsFromPackage edges
    n_imp = 0
    for (caller, target_folder), info in imports.items():
        if caller == target_folder:        # self-import via npm name (rare, ignore)
            continue
        bundle.add(ImportsFromPackage(
            from_id=c.pkg_id(caller),
            to_id=c.pkg_id(target_folder),
            imported_symbols=sorted(info["symbols"])[:20],
            import_count=info["count"],
            is_type_only=info["type_only_all"],
            evidence=[c.file_ref(p) for p in list(info["evidence"])[:5]],
            derived_by=Provenance.AST,
        ))
        n_imp += 1

    # Emit Calls edges (function-resolved when possible, else Package-targeted)
    n_call = n_resolved_friendly = 0
    for (caller, target_pkg, target_fn), info in calls.items():
        if caller == target_pkg and not info.get("count", 0):
            continue
        target_folder = folder_lookup.get(target_pkg.lower(), target_pkg)
        # Resolve via name first, then friendly_name (covers `Pkg:My Function`
        # spellings whose RegisteredFunction.name is `myFunction`).
        fid = c.func_id(target_folder, target_fn)
        lookup_key = (target_folder.lower(), target_fn.lower())
        if lookup_key in fn_lookup and fn_lookup[lookup_key] != fid:
            fid = fn_lookup[lookup_key]
            n_resolved_friendly += 1
        bundle.add(Calls(
            from_id=c.pkg_id(caller),
            to_id=fid,
            call_count=info["count"],
            callee_spec=f"{target_pkg}:{target_fn}",
            evidence=[c.file_ref(p) for p in list(info["evidence"])[:5]],
            derived_by=Provenance.AST,
        ))
        n_call += 1

    print(f"[{EXTRACTOR_NAME}] scanned {n_files} files, "
          f"{n_imp} ImportsFromPackage, {n_call} Calls "
          f"({n_resolved_friendly} resolved via friendly_name fallback)")


def _scan_file(path: Path, caller: str, folder_lookup: dict[str, str],
               imports: dict, calls: dict) -> None:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return
    if ("@datagrok/" not in text
            and "grok.functions." not in text
            and "_package.functions." not in text
            and "DG.Func." not in text):
        return

    # ImportsFromPackage
    for m in IMPORT_RE.finditer(text):
        npm_kebab = m.group("pkg")
        target = folder_lookup.get(npm_kebab)
        if not target:
            continue
        key = (caller, target)
        info = imports[key]
        info["count"] += 1
        info["evidence"].add(c.repo_rel(path))
        if not TYPE_ONLY_RE.match(m.group(0)):
            info["type_only_all"] = False
        # Pull symbols from {a, b as c}
        braced = m.group("braced")
        if braced:
            for sym in braced.strip("{}").split(","):
                sym = sym.strip().split(" as ")[0].strip()
                if sym:
                    info["symbols"].add(sym)

    # Calls — string-form: grok.functions.eval/call/find('Pkg:fn'), DG.Func.byName('Pkg:fn')
    for m in CALL_STRING_SPEC.finditer(text):
        _record_call(calls, caller, m.group("pkg"), m.group("fn"), path)

    # Calls — DG.Func.find({name:'fn', package:'Pkg'})
    for m in CALL_FIND_OBJ.finditer(text):
        _record_call(calls, caller, m.group("pkg"), m.group("fn"), path)


def _record_call(calls: dict, caller: str, target_pkg: str, target_fn: str, path: Path) -> None:
    if target_pkg.lower() == caller.lower():
        return
    key = (caller, target_pkg, target_fn)
    info = calls[key]
    info["count"] += 1
    info["evidence"].add(c.repo_rel(path))
