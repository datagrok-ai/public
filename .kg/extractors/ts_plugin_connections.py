"""Extract DataConnection entities from `packages/*/connections/*.json`.

Each JSON file is a single connection with `"#type": "DataConnection"`.

Emits:
  - DataConnection
  - HasConnection (Package -> DataConnection)
"""

from __future__ import annotations

import json
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import DataConnection
from schema.relations import HasConnection

from . import _common as c

EXTRACTOR_NAME = "ts_plugin_connections"


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
        cdir = pkg_dir / "connections"
        if not cdir.is_dir():
            continue
        for jf in sorted(cdir.glob("*.json")):
            try:
                obj = json.loads(jf.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if obj.get("#type") != "DataConnection":
                continue
            name = obj.get("name") or jf.stem
            cid = c.connection_id(pkg_dir.name, name)
            dc = DataConnection(
                id=cid, name=name, source_layer=SourceLayer.PLUGINS,
                package_id=c.pkg_id(pkg_dir.name),
                data_source=obj.get("dataSource"),
                description=obj.get("description"),
                description_provenance=Provenance.PACKAGE_JSON,
                paths=[c.file_ref(jf, role="definition")],
                extras={"friendlyName": obj.get("friendlyName"),
                        "tags": obj.get("tags"),
                        "parameters": obj.get("parameters")},
            )
            bundle.add(dc)
            bundle.add(HasConnection(
                from_id=c.pkg_id(pkg_dir.name), to_id=cid,
                derived_by=Provenance.PACKAGE_JSON))
            total += 1
    print(f"[{EXTRACTOR_NAME}] {total} connections")
