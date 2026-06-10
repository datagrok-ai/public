"""Extract ScriptEnvironment entities from `packages/*/environments/*.yaml`.

Conda-style YAML with `name:`, `channels:`, `dependencies:`.

Emits:
  - ScriptEnvironment
  - HasEnvironment (Package -> ScriptEnvironment)
"""

from __future__ import annotations

from pathlib import Path

import yaml

from schema.base import Provenance, SourceLayer
from schema.entities import ScriptEnvironment
from schema.relations import HasEnvironment

from . import _common as c

EXTRACTOR_NAME = "ts_plugin_environments"


def _flatten_deps(deps) -> list[str]:
    out: list[str] = []
    for d in deps or []:
        if isinstance(d, str):
            out.append(d)
        elif isinstance(d, dict):
            # e.g. {"pip": ["a", "b"]}
            for k, v in d.items():
                for it in v or []:
                    out.append(f"{k}:{it}")
    return out


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
        edir = pkg_dir / "environments"
        if not edir.is_dir():
            continue
        for yf in sorted(list(edir.glob("*.yaml")) + list(edir.glob("*.yml"))):
            try:
                obj = yaml.safe_load(yf.read_text(encoding="utf-8"))
            except (OSError, yaml.YAMLError):
                continue
            if not isinstance(obj, dict):
                continue
            name = obj.get("name") or yf.stem
            eid = c.env_id(pkg_dir.name, name)
            se = ScriptEnvironment(
                id=eid, name=name, source_layer=SourceLayer.PLUGINS,
                package_id=c.pkg_id(pkg_dir.name),
                channels=list(obj.get("channels") or []),
                deps=_flatten_deps(obj.get("dependencies")),
                paths=[c.file_ref(yf, role="definition")],
            )
            bundle.add(se)
            bundle.add(HasEnvironment(
                from_id=c.pkg_id(pkg_dir.name), to_id=eid,
                derived_by=Provenance.PACKAGE_JSON))
            total += 1
    print(f"[{EXTRACTOR_NAME}] {total} environments")
