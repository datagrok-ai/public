"""Derive `IsImplementedIn` and `IsTestedIn` from Feature → PART_OF_FEATURE
member entity → the file each member is defined in.

This is a pure synthesis step — it reads existing JSONL, doesn't touch the
filesystem. Run AFTER `files.py` (so File nodes + DefinedIn edges exist) AND
after the feature_clustering enricher has been applied.

Algorithm per Feature:
  1. Walk PART_OF_FEATURE → member ids
  2. For each member, look up its DefinedIn → File
  3. Group by file, count members
  4. Classify each File as test (`is_test == True`) or implementation
  5. Emit one IsImplementedIn or IsTestedIn edge per (Feature, File) pair

A Feature with zero IsTestedIn edges is a test gap.
A Feature with zero IsImplementedIn but ≥ 1 doc member is a docs-only feature.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from schema.base import Provenance, SourceLayer
from schema.relations import IsImplementedIn, IsTestedIn

from . import _common as c

EXTRACTOR_NAME = "feature_files"


def _read_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def run(repo_root: Path, bundle, package_filter: list[str] | None = None) -> None:
    nodes_dir = c.REPO_ROOT / ".kg" / "data" / "nodes"
    edges_dir = c.REPO_ROOT / ".kg" / "data" / "edges"
    if not (nodes_dir.is_dir() and edges_dir.is_dir()):
        print(f"[{EXTRACTOR_NAME}] no data/ yet — run extractors first")
        return

    # Build member → file map from DefinedIn
    member_to_file: dict[str, str] = {}
    for row in _read_jsonl(edges_dir / "DEFINED_IN.jsonl"):
        member_to_file[row["from_id"]] = row["to_id"]

    if not member_to_file:
        print(f"[{EXTRACTOR_NAME}] no DefinedIn edges yet — run files.py first")
        return

    # Build file → is_test map
    file_is_test: dict[str, bool] = {}
    for row in _read_jsonl(nodes_dir / "File.jsonl"):
        file_is_test[row["id"]] = bool(row.get("is_test", False))

    # Walk Feature → PART_OF_FEATURE → member
    pof = _read_jsonl(edges_dir / "PART_OF_FEATURE.jsonl")
    feat_files: dict[tuple[str, str], int] = defaultdict(int)
    for row in pof:
        member = row["from_id"]
        feat = row["to_id"]
        fid = member_to_file.get(member)
        if not fid:
            continue
        feat_files[(feat, fid)] += 1

    n_impl = n_test = 0
    for (feat_id, file_id), count in feat_files.items():
        is_test = file_is_test.get(file_id, False)
        if is_test:
            bundle.add(IsTestedIn(
                from_id=feat_id, to_id=file_id,
                member_count=count,
                derived_by=Provenance.FILESYSTEM,
            ))
            n_test += 1
        else:
            bundle.add(IsImplementedIn(
                from_id=feat_id, to_id=file_id,
                member_count=count,
                derived_by=Provenance.FILESYSTEM,
            ))
            n_impl += 1

    # Coverage stats
    feats = {row["id"] for row in _read_jsonl(nodes_dir / "Feature.jsonl")}
    feats_with_impl = {f for (f, _file), _n in feat_files.items()
                       if not file_is_test.get(_file, False)}
    feats_with_tests = {f for (f, _file), _n in feat_files.items()
                        if file_is_test.get(_file, False)}

    print(f"[{EXTRACTOR_NAME}] {n_impl} IsImplementedIn, {n_test} IsTestedIn")
    print(f"  {len(feats_with_impl)}/{len(feats)} features have impl files; "
          f"{len(feats_with_tests)}/{len(feats)} have test files; "
          f"{len(feats - feats_with_tests)} are test gaps.")
