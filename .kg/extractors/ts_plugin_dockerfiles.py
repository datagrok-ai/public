"""Extract DockerContainer entities from `packages/*/dockerfiles/<svc>/`.

Each subdir is one container, with a `Dockerfile` and optional
`container.json` (cpu, memory, gpu, on_demand).
"""

from __future__ import annotations

import json
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.entities import DockerContainer
from schema.relations import HasContainer

from . import _common as c

EXTRACTOR_NAME = "ts_plugin_dockerfiles"


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
        ddir = pkg_dir / "dockerfiles"
        if not ddir.is_dir():
            continue
        for sub in sorted(ddir.iterdir()):
            if not sub.is_dir():
                continue
            df = sub / "Dockerfile"
            if not df.is_file():
                continue
            cfg: dict = {}
            cj = sub / "container.json"
            if cj.is_file():
                try:
                    cfg = json.loads(cj.read_text(encoding="utf-8"))
                except json.JSONDecodeError:
                    cfg = {}
            cid = c.docker_id(pkg_dir.name, sub.name)
            dc = DockerContainer(
                id=cid, name=sub.name, source_layer=SourceLayer.PLUGINS,
                package_id=c.pkg_id(pkg_dir.name),
                cpu=str(cfg.get("cpu")) if cfg.get("cpu") is not None else None,
                memory=str(cfg.get("memory")) if cfg.get("memory") is not None else None,
                gpu=bool(cfg.get("gpu") or 0),
                on_demand=bool(cfg.get("on_demand") or False),
                paths=[c.file_ref(df, role="dockerfile")] +
                      ([c.file_ref(cj, role="config")] if cj.is_file() else []),
                extras=cfg,
            )
            bundle.add(dc)
            bundle.add(HasContainer(
                from_id=c.pkg_id(pkg_dir.name), to_id=cid,
                derived_by=Provenance.PACKAGE_JSON))
            total += 1
    print(f"[{EXTRACTOR_NAME}] {total} containers")
