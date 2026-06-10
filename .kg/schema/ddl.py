"""Derive Kuzu DDL from Pydantic models.

Maps Python types to Kuzu types. Lists become Kuzu LIST, dicts become
JSON-serialized STRING (queryable as text; Kuzu's MAP support is still
maturing). Optional fields are simply nullable.

We re-derive DDL on every build. If the schema changed, we drop and
re-create the affected tables — safe because the JSONL files are
authoritative, the DB is just a cache.
"""

from __future__ import annotations

import types
from datetime import date, datetime
from enum import Enum
from typing import Any, Union, get_args, get_origin

from pydantic import BaseModel
from pydantic.fields import FieldInfo

from .base import Edge, Entity
from .entities import ENTITY_KINDS
from .relations import RELATION_KINDS


# Fields stored in `extras` JSON blob and not unpacked into columns
_BASE_ENTITY_FIELDS = {"id", "name", "kind", "source_layer", "paths", "docs",
                       "description", "description_provenance", "extras",
                       "extracted_by", "extracted_at"}
_BASE_EDGE_FIELDS = {"from_id", "to_id", "predicate", "evidence", "derived_by",
                     "confidence", "extras", "extracted_by", "extracted_at"}


def _is_optional(annotation: Any) -> tuple[bool, Any]:
    """Return (is_optional, inner_type) for `T | None` / `Optional[T]`."""
    origin = get_origin(annotation)
    if origin is Union or origin is types.UnionType:
        args = [a for a in get_args(annotation) if a is not type(None)]
        if len(args) == 1:
            return True, args[0]
    return False, annotation


def _kuzu_type(annotation: Any) -> str:
    """Map a Python annotation to a Kuzu column type."""
    optional, inner = _is_optional(annotation)
    origin = get_origin(inner)

    if inner is str:
        return "STRING"
    if inner is bool:
        return "BOOLEAN"
    if inner is int:
        return "INT64"
    if inner is float:
        return "DOUBLE"
    if inner in (datetime,):
        return "TIMESTAMP"
    if inner is date:
        return "DATE"
    if isinstance(inner, type) and issubclass(inner, Enum):
        return "STRING"

    if origin is list:
        item_t = get_args(inner)[0] if get_args(inner) else str
        # Lists of primitives → typed LIST; lists of dicts/objects → STRING (JSON)
        if item_t in (str, int, float, bool):
            return f"{_kuzu_type(item_t)}[]"
        return "STRING"   # JSON-serialized

    if origin is dict or inner is dict:
        return "STRING"   # JSON-serialized

    if isinstance(inner, type) and issubclass(inner, BaseModel):
        return "STRING"   # JSON-serialized

    # Fallback — unknown types stored as JSON string
    return "STRING"


def _entity_columns(model: type[Entity]) -> list[tuple[str, str, bool]]:
    """Return [(column_name, kuzu_type, is_primary_key), ...]"""
    cols: list[tuple[str, str, bool]] = [
        ("id", "STRING", True),
        ("name", "STRING", False),
        ("kind", "STRING", False),
        ("source_layer", "STRING", False),
        ("paths", "STRING", False),                # JSON
        ("docs", "STRING", False),                 # JSON
        ("description", "STRING", False),
        ("description_provenance", "STRING", False),
        ("extras", "STRING", False),               # JSON
        ("extracted_by", "STRING", False),
        ("extracted_at", "STRING", False),
    ]
    for fname, finfo in model.model_fields.items():
        if fname in _BASE_ENTITY_FIELDS:
            continue
        cols.append((fname, _kuzu_type(finfo.annotation), False))
    return cols


def _edge_columns(model: type[Edge]) -> list[tuple[str, str]]:
    cols: list[tuple[str, str]] = [
        ("evidence", "STRING"),                    # JSON
        ("derived_by", "STRING"),
        ("confidence", "DOUBLE"),
        ("extras", "STRING"),                      # JSON
        ("extracted_by", "STRING"),
        ("extracted_at", "STRING"),
    ]
    for fname, finfo in model.model_fields.items():
        if fname in _BASE_EDGE_FIELDS:
            continue
        cols.append((fname, _kuzu_type(finfo.annotation)))
    return cols


def entity_ddl(model: type[Entity]) -> str:
    cols = _entity_columns(model)
    parts = [
        f"`{n}` {t}{' PRIMARY KEY' if pk else ''}" for n, t, pk in cols
    ]
    # Move PRIMARY KEY out of column list per Kuzu syntax: PRIMARY KEY(col)
    # Build properly:
    col_decls = [f"`{n}` {t}" for n, t, _pk in cols]
    pk_col = next((n for n, _t, pk in cols if pk), "id")
    body = ", ".join(col_decls) + f", PRIMARY KEY (`{pk_col}`)"
    return f"CREATE NODE TABLE IF NOT EXISTS `{model.__name__}` ({body});"


def edge_ddl(model: type[Edge]) -> str:
    """Kuzu requires (FROM X TO Y) pairs; we emit one rel table that lists
    every allowed (subject, object) combination."""
    cols = _edge_columns(model)
    col_decls = ", ".join(f"`{n}` {t}" for n, t in cols)
    pairs = []
    for s in model.SUBJECT_KINDS:
        for o in model.OBJECT_KINDS:
            pairs.append(f"FROM `{s}` TO `{o}`")
    pairs_str = ", ".join(pairs)
    body = pairs_str + (f", {col_decls}" if col_decls else "")
    return f"CREATE REL TABLE IF NOT EXISTS `{model.predicate_name()}` ({body});"


def all_ddl() -> list[str]:
    """All CREATE statements, node tables first."""
    statements: list[str] = []
    for model in ENTITY_KINDS.values():
        statements.append(entity_ddl(model))
    for model in RELATION_KINDS.values():
        statements.append(edge_ddl(model))
    return statements


# Add helper on Edge subclasses to fetch the predicate name (the auto-set default)
def _predicate_name(cls: type[Edge]) -> str:
    return cls.model_fields["predicate"].default


Edge.predicate_name = classmethod(_predicate_name)  # type: ignore[attr-defined]


__all__ = ["entity_ddl", "edge_ddl", "all_ddl"]
