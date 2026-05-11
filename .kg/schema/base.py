"""Base models for all entities and edges.

Every node in the graph is an `Entity`; every relationship is an `Edge`.
Both carry provenance so we always know where a fact came from and how
confident we are.

These models are layer-agnostic on purpose — the same `Entity` shape
works for a TS plugin function, a Dart viewer, a Jenkins job, or a
help-doc page. Layer-specific shape lives in subclasses in `entities.py`.
"""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class SourceLayer(str, Enum):
    PLUGINS = "plugins"
    LIBRARIES = "libraries"
    CORE_CLIENT = "core_client"
    CORE_SERVER = "core_server"
    CORE_SHARED = "core_shared"
    DOCS = "docs"
    PROCESS = "process"           # git, jira, changelog, jenkins
    INFRA = "infra"
    SYNTHETIC = "synthetic"       # LLM-derived (Feature, Concern, etc.)


class Provenance(str, Enum):
    """How a fact was derived."""
    ANNOTATION = "annotation"     # //name:, @Prop, @CmdFunc, frontmatter
    PACKAGE_JSON = "package_json"
    AST = "ast"
    REGISTRY = "registry"         # ViewerDescriptor.all, initEntities()
    FILESYSTEM = "filesystem"     # path-based inference
    GIT = "git"
    EXTERNAL_API = "external_api" # JIRA, GitHub, Jenkins
    LLM = "llm"
    MANUAL = "manual"


class FileRef(BaseModel):
    """Pointer to a file (or a span within it). All paths are repo-relative
    POSIX form (forward slashes), regardless of host OS."""
    model_config = ConfigDict(frozen=True)

    path: str
    line_start: int | None = None
    line_end: int | None = None
    role: str | None = None        # "definition", "usage", "test", "css", "doc"


class DocRef(BaseModel):
    """Pointer to documentation. May be a help/ markdown page, a README,
    a JSDoc block, a tutorial, or an external URL."""
    model_config = ConfigDict(frozen=True)

    kind: str                       # "help", "readme", "jsdoc", "tutorial", "external"
    target: str                     # path or URL
    anchor: str | None = None       # heading anchor within the doc
    title: str | None = None


class Entity(BaseModel):
    """Base node. Subclass per kind in `entities.py`. The `kind` field is
    auto-set from the subclass name and becomes the Kuzu node-table name."""
    model_config = ConfigDict(extra="forbid")

    id: str = Field(..., description="Globally unique, e.g. 'pkg:Chem' or 'func:Chem.detectMolecules'")
    name: str = Field(..., description="Human-readable display name")
    kind: str = Field(..., description="Auto-set from subclass name; matches the Kuzu node table")
    source_layer: SourceLayer

    paths: list[FileRef] = Field(default_factory=list, description="Code locations where this entity lives or is registered")
    docs: list[DocRef] = Field(default_factory=list, description="Documentation references")
    description: str | None = None
    description_provenance: Provenance | None = None

    # Free-form per-kind extras. Anything that doesn't deserve a typed field
    # but is worth keeping. Stored as JSON in Kuzu, queryable but not indexed.
    extras: dict[str, Any] = Field(default_factory=dict)

    # Tracking
    extracted_by: str | None = None    # extractor module name
    extracted_at: str | None = None    # ISO timestamp


class Edge(BaseModel):
    """Base edge. Subclass per predicate in `relations.py`. The `predicate`
    field becomes the Kuzu rel-table name."""
    model_config = ConfigDict(extra="forbid")

    from_id: str
    to_id: str
    predicate: str                  # auto-set from subclass name

    evidence: list[FileRef] = Field(default_factory=list, description="Where this relationship is detectable")
    derived_by: Provenance = Provenance.AST
    confidence: float = 1.0          # 0..1; 1.0 for deterministic facts

    extras: dict[str, Any] = Field(default_factory=dict)

    extracted_by: str | None = None
    extracted_at: str | None = None
