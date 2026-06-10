"""Shared utilities for all extractors.

- `Bundle` — accumulates emitted entities/edges and writes JSONL slices.
- `parse_annotation_block` — parses Datagrok's `//key: value` comment header
  used in TS, Python (`#`), R (`#`), Julia (`#`), JS (`//`), SQL (`--`).
- ID helpers — stable global IDs for every entity kind.
- Path helpers — POSIX-form repo-relative strings on any OS.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from schema.base import Edge, Entity, FileRef


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]


def repo_rel(p: Path | str) -> str:
    """Return a POSIX-form path relative to the repo root."""
    pp = Path(p).resolve()
    try:
        rel = pp.relative_to(REPO_ROOT)
    except ValueError:
        return str(pp).replace("\\", "/")
    return rel.as_posix()


def file_ref(p: Path | str, *, line_start: int | None = None,
             line_end: int | None = None, role: str | None = None) -> FileRef:
    return FileRef(path=repo_rel(p), line_start=line_start,
                   line_end=line_end, role=role)


# ---------------------------------------------------------------------------
# IDs
# ---------------------------------------------------------------------------

def pkg_id(folder: str) -> str:               return f"pkg:{folder}"
def lib_id(folder: str) -> str:               return f"lib:{folder}"
def lib_mod_id(lib: str, rel_path: str) -> str:
    return f"lib_mod:{lib}:{rel_path}"
def func_id(pkg: str, name: str) -> str:      return f"func:{pkg}:{name}"
def script_id(pkg: str, name: str) -> str:    return f"script:{pkg}:{name}"
def query_id(pkg: str, name: str) -> str:     return f"query:{pkg}:{name}"
def connection_id(pkg: str, name: str) -> str:return f"conn:{pkg}:{name}"
def env_id(pkg: str, name: str) -> str:       return f"env:{pkg}:{name}"
def docker_id(pkg: str, svc: str) -> str:     return f"docker:{pkg}:{svc}"
def wasm_id(pkg: str, fname: str) -> str:     return f"wasm:{pkg}:{fname}"
def clog_id(pkg: str, version: str, idx: int) -> str:
    return f"clog:{pkg}:{version}:{idx}"
def prop_id(pkg: str, name: str) -> str:      return f"prop:{pkg}:{name}"
def doc_id(rel_path: str) -> str:             return f"doc:{rel_path}"
def jira_id(key: str) -> str:                 return f"jira:{key}"
def commit_id(repo: str, sha: str) -> str:    return f"commit:{repo}:{sha[:12]}"
def feature_id(slug: str) -> str:             return f"feature:{slug}"


# ---------------------------------------------------------------------------
# Annotation block parser
# ---------------------------------------------------------------------------

# Match e.g.   //meta.role: panel    or    --input: string target = "X"
_ANN_LINE = re.compile(
    r"""^\s*
        (?:\/\/|\#|--)            # comment marker
        \s*
        ([a-zA-Z_][\w.\-]*)       # key (may include . for meta.foo)
        \s*:\s*
        (.*?)\s*$
    """, re.VERBOSE)

# `//input: <type> <name> {<json-ish>}` or `//input: <type> <name> = <default>`
_INPUT_LINE = re.compile(
    r"""^
        ([a-zA-Z_][\w]*)             # type
        \s+
        ([a-zA-Z_$][\w$]*)           # name
        (?:\s*=\s*(.+?))?            # = default
        (?:\s*\{(.*)\})?             # { options }
        \s*$
    """, re.VERBOSE)

_OUTPUT_LINE = re.compile(
    r"""^
        ([a-zA-Z_][\w]*)             # type
        (?:\s+([a-zA-Z_$][\w$]*))?   # optional name
        (?:\s*\{(.*)\})?
        \s*$
    """, re.VERBOSE)


def _parse_inline_options(s: str) -> dict[str, Any]:
    """Parse the `{key: value, key2: "v2"}` mini-DSL used in //input options."""
    s = s.strip()
    if not s:
        return {}
    out: dict[str, Any] = {}
    # Split on commas but not inside quotes/brackets. Cheap implementation.
    depth = 0
    parts: list[str] = []
    cur: list[str] = []
    in_q: str | None = None
    for ch in s:
        if in_q:
            cur.append(ch)
            if ch == in_q:
                in_q = None
            continue
        if ch in "\"'":
            in_q = ch
            cur.append(ch)
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
            continue
        cur.append(ch)
    if cur:
        parts.append("".join(cur))
    for part in parts:
        if ":" not in part:
            continue
        k, _, v = part.partition(":")
        k = k.strip()
        v = v.strip()
        if v.startswith('"') and v.endswith('"'):
            v = v[1:-1]
        elif v == "true":
            v = True  # type: ignore[assignment]
        elif v == "false":
            v = False  # type: ignore[assignment]
        elif v.startswith("[") and v.endswith("]"):
            try:
                v = json.loads(v)  # type: ignore[assignment]
            except json.JSONDecodeError:
                pass
        out[k] = v
    return out


def parse_annotation_block(lines: Iterable[str]) -> dict[str, Any]:
    """Parse a contiguous comment block into a structured dict.

    Returns:
        {
          "name": str | None,
          "description": str | None,
          "language": str | None,
          "tags": list[str],
          "inputs": [{name, type, default?, options}],
          "outputs": [{name?, type, options}],
          "meta": {dotted-key: value},
          "top_menu": str | None,
          "raw": {key: last-value}     # everything else
        }
    """
    out: dict[str, Any] = {
        "name": None, "description": None, "language": None,
        "friendly_name": None, "help_url": None,
        "tags": [], "inputs": [], "outputs": [], "meta": {},
        "top_menu": None, "raw": {},
    }
    for raw_line in lines:
        m = _ANN_LINE.match(raw_line)
        if not m:
            continue
        key, value = m.group(1).strip(), m.group(2).strip()

        if key.startswith("meta."):
            out["meta"][key[5:]] = value
            continue
        kl = key.lower()
        if kl == "name":
            out["name"] = value
        elif kl == "description":
            out["description"] = value
        elif kl in ("friendlyname", "friendly-name"):
            out["friendly_name"] = value
        elif kl == "language":
            out["language"] = value.lower()
        elif kl == "help-url":
            out["help_url"] = value
        elif kl == "top-menu":
            out["top_menu"] = value
        elif kl == "tags":
            out["tags"] = [t.strip() for t in value.split(",") if t.strip()]
        elif kl == "input":
            mi = _INPUT_LINE.match(value)
            if mi:
                t, n, default, opts_s = mi.groups()
                inp: dict[str, Any] = {"name": n, "type": t}
                if default is not None:
                    inp["default"] = default.strip().strip('"')
                if opts_s:
                    inp["options"] = _parse_inline_options(opts_s)
                out["inputs"].append(inp)
            else:
                out["inputs"].append({"raw": value})
        elif kl == "output":
            mo = _OUTPUT_LINE.match(value)
            if mo:
                t, n, opts_s = mo.groups()
                outp: dict[str, Any] = {"type": t}
                if n:
                    outp["name"] = n
                if opts_s:
                    outp["options"] = _parse_inline_options(opts_s)
                out["outputs"].append(outp)
            else:
                out["outputs"].append({"raw": value})
        else:
            out["raw"][kl] = value

    return out


def find_annotation_blocks(text: str) -> list[tuple[int, int, list[str], str]]:
    """Find contiguous //, #, or -- comment blocks plus the line that follows.

    Returns list of (start_line_1based, end_line_1based, comment_lines, anchor_line).
    Comment markers must be uniform within a block (no mixing).
    """
    lines = text.splitlines()
    blocks: list[tuple[int, int, list[str], str]] = []
    i = 0
    while i < len(lines):
        ln = lines[i]
        m = re.match(r"\s*(\/\/|\#|--)", ln)
        if not m:
            i += 1
            continue
        marker = m.group(1)
        start = i
        block: list[str] = []
        while i < len(lines):
            ln2 = lines[i]
            m2 = re.match(r"\s*(\/\/|\#|--)", ln2)
            if not m2 or m2.group(1) != marker:
                break
            block.append(ln2)
            i += 1
        anchor = lines[i] if i < len(lines) else ""
        # Only keep blocks that contain at least one structured `key:` line
        if any(_ANN_LINE.match(b) for b in block):
            blocks.append((start + 1, i, block, anchor))
        # If the anchor wasn't a comment we don't advance; loop will move on
        if i == start:
            i += 1
    return blocks


# ---------------------------------------------------------------------------
# Bundle: collect entities/edges, write JSONL slices.
# ---------------------------------------------------------------------------

class Bundle:
    """Accumulate emitted entities/edges and flush them to JSONL slices.

    JSONL files are sharded by entity-kind (one file per kind) and by
    relation-predicate (one file per predicate). This keeps the canonical
    layout simple and predictable for both Kuzu COPY and code review.
    """

    def __init__(self, data_dir: Path, *, extractor_name: str):
        self.data_dir = data_dir
        self.extractor_name = extractor_name
        self.entities: dict[str, list[Entity]] = defaultdict(list)
        self.edges: dict[str, list[Edge]] = defaultdict(list)
        self._now = datetime.now(timezone.utc).isoformat(timespec="seconds")

    def add(self, item: Entity | Edge) -> None:
        if isinstance(item, Entity):
            item.extracted_by = self.extractor_name
            item.extracted_at = self._now
            self.entities[item.kind].append(item)
        elif isinstance(item, Edge):
            item.extracted_by = self.extractor_name
            item.extracted_at = self._now
            self.edges[item.predicate].append(item)
        else:
            raise TypeError(f"Bundle.add: not an Entity/Edge: {type(item)}")

    def add_many(self, items: Iterable[Entity | Edge]) -> None:
        for it in items:
            self.add(it)

    def flush(self) -> dict[str, int]:
        """Write JSONL slices. Returns counts per kind/predicate."""
        nodes_dir = self.data_dir / "nodes"
        edges_dir = self.data_dir / "edges"
        nodes_dir.mkdir(parents=True, exist_ok=True)
        edges_dir.mkdir(parents=True, exist_ok=True)

        counts: dict[str, int] = {}
        for kind, items in self.entities.items():
            self._merge_jsonl(nodes_dir / f"{kind}.jsonl", items)
            counts[f"node:{kind}"] = len(items)
        for predicate, items in self.edges.items():
            self._merge_jsonl(edges_dir / f"{predicate}.jsonl", items)
            counts[f"edge:{predicate}"] = len(items)
        return counts

    def _merge_jsonl(self, path: Path, new_items: list[Entity | Edge]) -> None:
        """Merge new items into an existing JSONL slice, deduped by ID
        (entities) or (from_id, to_id, predicate) (edges).

        This makes extractors idempotent and lets multiple extractors emit
        into the same file without overwriting each other's contributions.
        """
        existing: dict[str, str] = {}
        if path.exists():
            for line in path.read_text(encoding="utf-8").splitlines():
                if not line.strip():
                    continue
                obj = json.loads(line)
                key = obj.get("id") or f"{obj.get('from_id')}|{obj.get('to_id')}|{obj.get('predicate')}"
                existing[key] = line
        for it in new_items:
            d = it.model_dump(mode="json")
            key = d.get("id") or f"{d.get('from_id')}|{d.get('to_id')}|{d.get('predicate')}"
            existing[key] = json.dumps(d, separators=(",", ":"), ensure_ascii=False)
        with path.open("w", encoding="utf-8", newline="\n") as f:
            for line in existing.values():
                f.write(line + "\n")


__all__ = [
    "Bundle", "REPO_ROOT", "repo_rel", "file_ref",
    "pkg_id", "lib_id", "lib_mod_id", "func_id", "script_id", "query_id",
    "connection_id", "env_id", "docker_id", "wasm_id", "clog_id",
    "prop_id", "doc_id", "jira_id", "commit_id", "feature_id",
    "parse_annotation_block", "find_annotation_blocks",
]
