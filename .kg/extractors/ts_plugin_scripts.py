"""Extract Script entities from `packages/*/scripts/`.

Supports Python (`.py`), R (`.r`/`.R`), JavaScript (`.js`/`.mjs`), Julia
(`.jl`), Octave (`.m`). Each script begins with a comment header in the
language-appropriate marker (`#` for Py/R/Julia/Octave, `//` for JS).

Emits:
  - Script
  - HasScript (Package -> Script)
  - RequiresEnvironment (Script -> ScriptEnvironment)  [populated when env name matches]
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from schema.base import Provenance, SourceLayer
from schema.entities import Script
from schema.relations import HasScript, RequiresEnvironment

from . import _common as c

EXTRACTOR_NAME = "ts_plugin_scripts"

LANG_BY_EXT = {
    ".py": "python", ".r": "r", ".R": "r", ".jl": "julia",
    ".m": "octave", ".js": "js", ".mjs": "js",
}


def _read_header(path: Path) -> list[str]:
    """Read until the first non-comment, non-blank line."""
    head: list[str] = []
    try:
        with path.open(encoding="utf-8", errors="replace") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith("#") or s.startswith("//"):
                    head.append(line.rstrip("\n"))
                    continue
                break
    except OSError:
        return []
    return head


def _scan_dir(scripts_dir: Path, pkg_name: str, pkg_node_id: str, bundle) -> int:
    n = 0
    for path in sorted(scripts_dir.rglob("*")):
        if not path.is_file():
            continue
        lang = LANG_BY_EXT.get(path.suffix)
        if lang is None:
            continue
        head = _read_header(path)
        if not head:
            continue
        ann = c.parse_annotation_block(head)
        name = ann.get("name") or path.stem
        sid = c.script_id(pkg_name, name)
        s = Script(
            id=sid, name=name, source_layer=SourceLayer.PLUGINS,
            package_id=pkg_node_id,
            language=ann.get("language") or lang,
            inputs=ann.get("inputs", []),
            outputs=ann.get("outputs", []),
            environment=ann.get("raw", {}).get("environment"),
            description=ann.get("description"),
            description_provenance=Provenance.ANNOTATION if ann.get("description") else None,
            paths=[c.file_ref(path, role="definition")],
            extras={"top_menu": ann.get("top_menu"), "tags": ann.get("tags"),
                    "meta": ann.get("meta") or None},
        )
        bundle.add(s)
        bundle.add(HasScript(from_id=pkg_node_id, to_id=sid,
                             derived_by=Provenance.ANNOTATION))
        env = s.environment
        if env:
            bundle.add(RequiresEnvironment(
                from_id=sid, to_id=c.env_id(pkg_name, env),
                derived_by=Provenance.ANNOTATION))
        n += 1
    return n


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    pkgs = repo_root / "packages"
    if not pkgs.is_dir():
        return
    total = 0
    for pkg_dir in sorted(pkgs.iterdir()):
        if not pkg_dir.is_dir():
            continue
        if package_filter and pkg_dir.name not in package_filter:
            continue
        scripts_dir = pkg_dir / "scripts"
        if not scripts_dir.is_dir():
            continue
        total += _scan_dir(scripts_dir, pkg_dir.name, c.pkg_id(pkg_dir.name), bundle)
    print(f"[{EXTRACTOR_NAME}] {total} scripts")
