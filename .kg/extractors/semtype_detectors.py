"""Extract Datagrok semantic-type detectors from `detectors.js`/`detectors.ts`.

The Datagrok platform discovers semType detectors by walking the
`detectors.js` file in each package. Inside, any class extending
`DG.Package` whose method names match `detect<SemType>(col)` becomes a
synchronous detector for that semType. The platform doesn't require an
annotation comment.

This extractor mirrors that discovery: each `detect<X>` method becomes a
synthetic `RegisteredFunction` (role=`semTypeDetector`) and emits a
`DetectsSemtype` edge to the SemanticType node.

Why synthetic? Because `ts_plugin_package` already creates
RegisteredFunction rows from annotated blocks; detectors are unannotated.
We give them a stable `func:<Pkg>:detect<Sem>` ID and let `semtypes.py`
pick up the connection on its next pass — but to avoid an extractor-order
dance we emit DetectsSemtype directly here.
"""

from __future__ import annotations

import re
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import RegisteredFunction, SemanticType
from schema.relations import DetectsSemtype, Exports

from . import _common as c

EXTRACTOR_NAME = "semtype_detectors"

# Match `class XxxPackageDetectors extends DG.Package { ... }` — body bracket
CLASS_RE = re.compile(
    r"class\s+(?P<class>[A-Za-z_]\w*)\s+extends\s+DG\.Package\s*\{",
    re.MULTILINE,
)
# Match `detectMolecule(col) { ... }` or `detectFoo (col, row) {` — JS ES2015
# class-body method form. We don't try to parse arrow-style class fields.
DETECT_METHOD_RE = re.compile(
    r"^\s+(?:async\s+)?detect(?P<sem>[A-Z][\w]*)\s*\([^)]*\)\s*\{",
    re.MULTILINE,
)


def _find_class_body(text: str, open_brace_pos: int) -> tuple[int, int]:
    """Return (start, end) of the brace-balanced body that begins at
    open_brace_pos (which must point at `{`). Quote-aware; handles `/.../`
    regex by skipping over `//`-comment-prone characters with a tiny
    inline tokenizer."""
    depth = 1
    i = open_brace_pos + 1
    in_q: str | None = None
    in_c: str | None = None    # 'line' | 'block'
    while i < len(text) and depth > 0:
        ch = text[i]
        if in_c == "line":
            if ch == "\n": in_c = None
            i += 1; continue
        if in_c == "block":
            if ch == "*" and i + 1 < len(text) and text[i + 1] == "/":
                in_c = None; i += 2; continue
            i += 1; continue
        if in_q:
            if ch == "\\" and i + 1 < len(text):
                i += 2; continue
            if ch == in_q:
                in_q = None
            i += 1; continue
        if ch in "'\"`":
            in_q = ch; i += 1; continue
        if ch == "/" and i + 1 < len(text):
            nxt = text[i + 1]
            if nxt == "/":
                in_c = "line"; i += 2; continue
            if nxt == "*":
                in_c = "block"; i += 2; continue
        if ch == "{": depth += 1
        elif ch == "}": depth -= 1
        i += 1
    return (open_brace_pos, i)


def _scan_detectors_file(text: str, file_path: Path, pkg_name: str,
                          bundle) -> int:
    n = 0
    for cm in CLASS_RE.finditer(text):
        # Find body
        brace = text.find("{", cm.end() - 1)
        if brace < 0: continue
        body_start, body_end = _find_class_body(text, brace)
        body = text[body_start: body_end]
        for dm in DETECT_METHOD_RE.finditer(body):
            sem_raw = dm.group("sem")
            # Convert PascalCase/camelCase to a plausible semType label.
            # detectMacromolecule → Macromolecule; detectChemicalReaction →
            # ChemicalReaction. These match Datagrok's canonical names.
            sem = sem_raw
            sid = f"semtype:{sem}"
            # Stable per-detector ID
            fname = f"detect{sem_raw}"
            fid = f"func:{pkg_name}:{fname}"
            line_no = (text[:body_start + dm.start()]).count("\n") + 1
            bundle.add(RegisteredFunction(
                id=fid, name=fname,
                source_layer=SourceLayer.PLUGINS,
                package_id=c.pkg_id(pkg_name),
                friendly_name=f"detect{sem_raw}",
                role="semTypeDetector",
                tags=["semTypeDetector"],
                meta={"semType": sem},
                language="js" if file_path.suffix == ".js" else "ts",
                description=f"Auto-discovered detector for semType '{sem}' (detectors.js method).",
                description_provenance=Provenance.AST,
                paths=[c.file_ref(file_path, line_start=line_no, role="definition")],
            ))
            bundle.add(Exports(
                from_id=c.pkg_id(pkg_name), to_id=fid,
                derived_by=Provenance.AST,
            ))
            bundle.add(SemanticType(
                id=sid, name=sem, source_layer=SourceLayer.SYNTHETIC,
                semtype_value=sem,
                detected_by_count=0, consumed_by_count=0,
                description=None,
            ))
            bundle.add(DetectsSemtype(
                from_id=fid, to_id=sid,
                derived_by=Provenance.AST,
            ))
            n += 1
    return n


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return
    n_pkgs = n_detectors = 0
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        for cand in (pkg_dir / "detectors.js", pkg_dir / "detectors.ts",
                     pkg_dir / "src" / "detectors.js",
                     pkg_dir / "src" / "detectors.ts"):
            if not cand.is_file():
                continue
            try:
                text = cand.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            added = _scan_detectors_file(text, cand, pkg_dir.name, bundle)
            if added:
                n_pkgs += 1
                n_detectors += added
    print(f"[{EXTRACTOR_NAME}] {n_detectors} semType detectors across {n_pkgs} packages")
