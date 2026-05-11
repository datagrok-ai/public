"""Concrete edge types.

Each subclass of `Edge` becomes a Kuzu rel table. The `predicate` field
is auto-populated from the subclass name (UPPER_SNAKE_CASE).

Each relation declares the allowed source/target entity kinds via
`SUBJECT_KINDS` and `OBJECT_KINDS`. This makes Kuzu's typed REL tables
declarative and lets us catch wiring mistakes at extraction time.
"""

from __future__ import annotations

import re
from typing import Any, ClassVar

from pydantic import Field

from .base import Edge


RELATION_KINDS: dict[str, type[Edge]] = {}


def _camel_to_upper_snake(name: str) -> str:
    s1 = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s1).upper()


class _RelBase(Edge):
    SUBJECT_KINDS: ClassVar[list[str]] = []
    OBJECT_KINDS: ClassVar[list[str]] = []

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        predicate = _camel_to_upper_snake(cls.__name__)
        cls.model_fields["predicate"].default = predicate
        RELATION_KINDS[predicate] = cls


# ============================================================================
# Structural — plugins
# ============================================================================

class Exports(_RelBase):
    """Package exports a registered function."""
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["RegisteredFunction"]


class HasScript(_RelBase):
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["Script"]


class HasQuery(_RelBase):
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["DataQuery"]


class HasConnection(_RelBase):
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["DataConnection"]


class HasEnvironment(_RelBase):
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["ScriptEnvironment"]


class HasContainer(_RelBase):
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["DockerContainer"]


class HasWasmModule(_RelBase):
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["WasmModule"]


class HasProperty(_RelBase):
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["PackageProperty"]


class HasChangelogEntry(_RelBase):
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["ChangelogEntry"]


# ============================================================================
# Cross-package and library coupling
# ============================================================================

class DependsOn(_RelBase):
    """Package depends on another Package or a Library (npm-style dep).
    Carries `semver_range` and `is_dev` flag for distinguishing runtime vs dev."""
    SUBJECT_KINDS = ["Package", "Library"]
    OBJECT_KINDS = ["Package", "Library"]
    semver_range: str | None = None
    is_dev: bool = False
    is_test_harness: bool = False    # @datagrok-libraries/test special-case


class ImportsFromModule(_RelBase):
    """Package imports symbols from a specific LibraryModule via deep path
    `@datagrok-libraries/<lib>/src/<path>`."""
    SUBJECT_KINDS = ["Package", "Library"]
    OBJECT_KINDS = ["LibraryModule"]
    imported_symbols: list[str] = Field(default_factory=list)
    import_count: int = 0
    is_type_only: bool = False       # `import type { ... }`


class ImportsFromPackage(_RelBase):
    """Cross-package TS import: one Package imports symbols from another via
    `import ... from '@datagrok/<other>'`. Stronger signal than package.json
    DEPENDS_ON because it proves the dep is actually used at the source level."""
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["Package"]
    imported_symbols: list[str] = Field(default_factory=list)
    import_count: int = 0
    is_type_only: bool = False


class Calls(_RelBase):
    """Cross-package runtime function call detected via
    `grok.functions.eval('<Pkg>:<fn>')` or `grok.functions.call('<Pkg>:<fn>')`.
    Subject is the calling Package; object is the called RegisteredFunction
    (or Package if the function isn't in the graph — dynamic name, etc.)."""
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["RegisteredFunction", "Package"]
    call_count: int = 0
    callee_spec: str | None = None    # raw "Pkg:fn" string


# ============================================================================
# Multi-language: scripts/queries → environments/connections
# ============================================================================

class RequiresEnvironment(_RelBase):
    SUBJECT_KINDS = ["Script", "RegisteredFunction"]
    OBJECT_KINDS = ["ScriptEnvironment"]


class UsesConnection(_RelBase):
    SUBJECT_KINDS = ["DataQuery"]
    OBJECT_KINDS = ["DataConnection"]


class UsesContainer(_RelBase):
    SUBJECT_KINDS = ["RegisteredFunction", "Script"]
    OBJECT_KINDS = ["DockerContainer"]


# ============================================================================
# Tags / typing
# ============================================================================

class HasTag(_RelBase):
    """RegisteredFunction has a tag (`cellEditor`, `app`, `panel`, `init`,
    `autostart`, …). One edge per (function, tag). Promoted from
    `RegisteredFunction.tags[]` so queries like
    `MATCH (rf)-[:HAS_TAG]->(:Tag {name:'cellEditor'}) RETURN rf` work
    without scanning every RF's tag array."""
    SUBJECT_KINDS = ["RegisteredFunction"]
    OBJECT_KINDS = ["Tag"]


# ============================================================================
# Documentation
# ============================================================================

class Documents(_RelBase):
    """Doc page documents a function/script/query/tutorial/feature/package.

    Direction is doc → code so a query like 'show me what scatter-plot.md
    documents' reads naturally as `MATCH (d:DocPage {url:'…'})-[:DOCUMENTS]->(x)`.
    To go the other way: `MATCH (f:RegisteredFunction)<-[:DOCUMENTS]-(d:DocPage)`.
    """
    SUBJECT_KINDS = ["DocPage", "HelpAnchor"]
    OBJECT_KINDS = ["RegisteredFunction", "Script", "DataQuery", "Tutorial",
                    "Package", "Library", "Feature"]
    derivation: str | None = None    # 'help_url' | 'this_help_url' | 'fuzzy_name' | 'llm_prose'


class Demonstrates(_RelBase):
    """Sample/tutorial demonstrates a function/feature."""
    SUBJECT_KINDS = ["ApiSample", "PlaywrightScenario", "DocPage", "Tutorial"]
    OBJECT_KINDS = ["RegisteredFunction", "Feature"]


class HasAnchor(_RelBase):
    SUBJECT_KINDS = ["DocPage"]
    OBJECT_KINDS = ["HelpAnchor"]


class LinksTo(_RelBase):
    """Doc page links to another doc page (relative .md link)."""
    SUBJECT_KINDS = ["DocPage"]
    OBJECT_KINDS = ["DocPage", "HelpAnchor"]


class HasTutorial(_RelBase):
    """Track contains tutorials; package owns the track."""
    SUBJECT_KINDS = ["Package", "TutorialTrack"]
    OBJECT_KINDS = ["Tutorial", "TutorialTrack"]


# ============================================================================
# Tests
# ============================================================================

class Covers(_RelBase):
    """Test covers a function/feature/package."""
    SUBJECT_KINDS = ["ApiTest", "PlaywrightScenario", "PackageTest"]
    OBJECT_KINDS = ["RegisteredFunction", "Feature", "Package"]


# ============================================================================
# Process
# ============================================================================

class FixedIn(_RelBase):
    """Jira ticket fixed in a commit."""
    SUBJECT_KINDS = ["JiraTicket"]
    OBJECT_KINDS = ["Commit"]


class MentionsTicket(_RelBase):
    """Changelog/commit/doc mentions a Jira ticket."""
    SUBJECT_KINDS = ["ChangelogEntry", "Commit", "DocPage"]
    OBJECT_KINDS = ["JiraTicket"]


# ============================================================================
# Synthetic — LLM-derived
# ============================================================================

class PartOfFeature(_RelBase):
    """Entity is part of a Feature cluster (LLM-derived)."""
    SUBJECT_KINDS = ["RegisteredFunction", "Script", "DataQuery", "DocPage",
                     "Tutorial", "ApiTest", "ApiSample", "ChangelogEntry"]
    OBJECT_KINDS = ["Feature"]
    weight: float = 1.0
    role: str | None = None             # 'core' | 'doc' | 'test' | 'sample' | 'related'


class HasFeature(_RelBase):
    """Package owns this Feature (primary). LLM-derived; one Feature can
    belong to one Package primarily, even if it spans multiple."""
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["Feature"]


class RelatesToFeature(_RelBase):
    """Cross-feature relationships within a package or domain."""
    SUBJECT_KINDS = ["Feature"]
    OBJECT_KINDS = ["Feature"]
    relation_kind: str | None = None    # 'subfeature' | 'sibling' | 'depends_on' | 'related'


# ============================================================================
# Files (cross-cutting)
# ============================================================================

class ContainsFile(_RelBase):
    """A Package or Library owns this File (file lives under that root)."""
    SUBJECT_KINDS = ["Package", "Library"]
    OBJECT_KINDS = ["File"]


class DefinedIn(_RelBase):
    """A code entity is defined in this file. Carries the line range when known."""
    SUBJECT_KINDS = ["RegisteredFunction", "Script", "DataQuery", "DataConnection",
                     "Tutorial", "TutorialTrack", "DocPage", "ChangelogEntry",
                     "TsClass", "TsMethod", "TsFunction", "TsInterface",
                     "TsEnum", "TsConstant", "PackageTest"]
    OBJECT_KINDS = ["File"]
    line_start: int | None = None
    line_end: int | None = None
    role: str | None = None              # 'definition' | 'usage' | 'test' | …


class IsImplementedIn(_RelBase):
    """Feature is implemented (in part) in this File. Derived from PART_OF_FEATURE
    members → their DefinedIn files, restricted to non-test files. A Feature
    typically has many implementation files; the `member_count` field counts
    how many of the Feature's members live in this file."""
    SUBJECT_KINDS = ["Feature"]
    OBJECT_KINDS = ["File"]
    member_count: int = 0


class IsTestedIn(_RelBase):
    """Feature is tested in this File. Same derivation as IsImplementedIn but
    restricted to files classified as test (path contains /tests/, ends with
    .test.ts, etc.). If a Feature has zero IsTestedIn edges, it's a test gap."""
    SUBJECT_KINDS = ["Feature"]
    OBJECT_KINDS = ["File"]
    member_count: int = 0


class Imports(_RelBase):
    """A source file imports from another file or library module via an
    ES module `import` statement (or dynamic `import()`). Aggregated to one
    edge per (from_file, to_target); `imported_symbols` lists every named
    symbol seen across all matching statements.

    Resolution:
      - relative paths (./foo, ../bar) → File node in the same package/library
      - `@datagrok-libraries/<lib>/...` → LibraryModule node
      - `datagrok-api/{dg,grok,ui}` → JsApiNamespace node
      - `@datagrok/<pkg>/...` → File in target package (rare in practice)

    Reverse query: `MATCH (other)-[:IMPORTS]->(target) WHERE target.id = '...'`
    answers "who imports this file" — this is the file-level back-reference
    that was missing before.
    """
    SUBJECT_KINDS = ["File"]
    OBJECT_KINDS = ["File", "LibraryModule", "JsApiNamespace"]
    imported_symbols: list[str] = Field(default_factory=list)
    is_type_only: bool = False
    is_default: bool = False
    is_namespace: bool = False
    is_dynamic: bool = False             # `await import(...)`
    import_count: int = 1                # number of distinct import statements


class RequiresColumnTag(_RelBase):
    """A RegisteredFunction's `meta.columnTags` constrains which columns the
    platform will dispatch it for. Each `key=value` pair becomes one edge to
    a SemanticType-style node; `tag_key` records the key (e.g. `quality`,
    `units`, `cell.renderer`) and the SemanticType node holds the value.

    Today the platform's most-used keys are `quality` (semType-equivalent)
    and `units` (notation). Promoting these to edges lets queries like
    "find all cellEditors valid for an OligoNucleotide column" work without
    string-grepping `meta`.
    """
    SUBJECT_KINDS = ["RegisteredFunction"]
    OBJECT_KINDS = ["SemanticType"]
    tag_key: str | None = None           # 'quality' | 'units' | 'cell.renderer' | ...
    tag_value: str | None = None         # raw value string


# ============================================================================
# Semantic types (cross-cutting)
# ============================================================================

class DetectsSemtype(_RelBase):
    """A `semTypeDetector` function declares this semantic type."""
    SUBJECT_KINDS = ["RegisteredFunction"]
    OBJECT_KINDS = ["SemanticType"]


class ConsumesSemtype(_RelBase):
    """A function's input column has `{semType: 'X'}` — function expects this type."""
    SUBJECT_KINDS = ["RegisteredFunction"]
    OBJECT_KINDS = ["SemanticType"]


class ProducesSemtype(_RelBase):
    """A function's output column has `{semType: 'X'}` — function emits this type."""
    SUBJECT_KINDS = ["RegisteredFunction"]
    OBJECT_KINDS = ["SemanticType"]


# ============================================================================
# JS API structure
# ============================================================================

class HasMethod(_RelBase):
    """TsClass owns this TsMethod."""
    SUBJECT_KINDS = ["TsClass"]
    OBJECT_KINDS = ["TsMethod"]


class HasFunction(_RelBase):
    """Package or Library exports a top-level TsFunction (helper fn,
    not the platform-registered //name: kind which is RegisteredFunction)."""
    SUBJECT_KINDS = ["Package", "Library"]
    OBJECT_KINDS = ["TsFunction"]


class HasClass(_RelBase):
    """Package or Library declares this TsClass / TsInterface / TsEnum / TsConstant
    in its source. One edge-kind for all four — target's `kind` field discriminates."""
    SUBJECT_KINDS = ["Package", "Library"]
    OBJECT_KINDS = ["TsClass", "TsInterface", "TsEnum", "TsConstant"]


class ExtendsClass(_RelBase):
    """One TsClass extends another (TS `extends`)."""
    SUBJECT_KINDS = ["TsClass"]
    OBJECT_KINDS = ["TsClass"]


class ImplementsInterface(_RelBase):
    """TsClass implements a TsInterface (TS `implements`)."""
    SUBJECT_KINDS = ["TsClass"]
    OBJECT_KINDS = ["TsInterface"]


class InNamespace(_RelBase):
    """A JS API entity belongs to a JS API namespace (DG, grok, ui, nested).
    Only used for js-api entities."""
    SUBJECT_KINDS = ["TsClass", "TsEnum", "TsConstant", "TsInterface", "TsFunction",
                     "DapiEndpoint", "JsEventStream", "UiComponent", "JsApiNamespace"]
    OBJECT_KINDS = ["JsApiNamespace"]


class DelegatesTo(_RelBase):
    """JsApiMethod calls a GeneratedBinding (Dart-side handler) at runtime."""
    SUBJECT_KINDS = ["TsMethod"]
    OBJECT_KINDS = ["GeneratedBinding"]


class TypedFor(_RelBase):
    """DapiEndpoint serves entities of a typed kind (the T in HttpDataSource<T>)."""
    SUBJECT_KINDS = ["DapiEndpoint"]
    OBJECT_KINDS = ["TsClass"]      # the entity class — User, Project, etc.


class DefinedByFunction(_RelBase):
    """RegisteredFunction is implemented by this top-level TsFunction.
    Bridges the //name:-annotated wrapper to its actual TS function."""
    SUBJECT_KINDS = ["RegisteredFunction"]
    OBJECT_KINDS = ["TsFunction", "TsMethod"]


class DefinesSemtype(_RelBase):
    """A JsApiConstant value defines a Datagrok semantic type — bridges
    `DG.SEMTYPE.MOLECULE` to the SemanticType node 'Molecule'."""
    SUBJECT_KINDS = ["TsConstant"]
    OBJECT_KINDS = ["SemanticType"]


# ============================================================================
# JS API usage (Package -> JS API)
# ============================================================================

class ImportsNamespace(_RelBase):
    """A code unit imports `* as DG from 'datagrok-api/dg'` (or grok / ui).
    Emitted at three granularities so queries can attribute use precisely:
    - `File`        → the importing TS file (one edge per import statement)
    - `Package`     → aggregate count per package (existing)
    - `Library`     → aggregate count per library (existing)
    """
    SUBJECT_KINDS = ["Package", "Library", "File"]
    OBJECT_KINDS = ["JsApiNamespace"]
    import_count: int = 0


class UsesApiClass(_RelBase):
    """A code unit references `DG.<JsApiClass>` (type, instantiation, static).
    Three granularities — File / containing TsMethod or TsFunction / Package
    — so you can answer 'which method uses DG.Grid?' not just 'which package'.
    """
    SUBJECT_KINDS = ["Package", "Library", "File", "TsMethod", "TsFunction"]
    OBJECT_KINDS = ["TsClass"]
    use_count: int = 0


class UsesApiEnum(_RelBase):
    """A code unit references `DG.<JsApiEnum>.VALUE` or `DG.<JsApiConstant>.KEY`.
    See UsesApiClass for granularity notes."""
    SUBJECT_KINDS = ["Package", "Library", "File", "TsMethod", "TsFunction"]
    OBJECT_KINDS = ["TsEnum", "TsConstant"]
    use_count: int = 0


class CallsDapiEndpoint(_RelBase):
    """A code unit calls `grok.dapi.<endpoint>` for server access.
    Three granularities — File / TsMethod or TsFunction / Package."""
    SUBJECT_KINDS = ["Package", "Library", "File", "TsMethod", "TsFunction"]
    OBJECT_KINDS = ["DapiEndpoint"]
    call_count: int = 0


class SubscribesToEvent(_RelBase):
    """A code unit subscribes to `grok.events.<stream>`.
    Three granularities — File / TsMethod or TsFunction / Package."""
    SUBJECT_KINDS = ["Package", "Library", "File", "TsMethod", "TsFunction"]
    OBJECT_KINDS = ["JsEventStream"]
    subscribe_count: int = 0


class UsesUiComponent(_RelBase):
    """A code unit calls `ui.<factory>(...)` to build UI.
    Three granularities — File / TsMethod or TsFunction / Package. Lets you
    answer 'which functions/methods build a ui.dialog' end-to-end."""
    SUBJECT_KINDS = ["Package", "Library", "File", "TsMethod", "TsFunction"]
    OBJECT_KINDS = ["UiComponent"]
    use_count: int = 0


__all__ = [
    "RELATION_KINDS",
    "Exports", "HasScript", "HasQuery", "HasConnection", "HasEnvironment",
    "HasContainer", "HasWasmModule", "HasProperty", "HasChangelogEntry",
    "DependsOn", "ImportsFromModule", "ImportsFromPackage", "Calls",
    "RequiresEnvironment", "UsesConnection", "UsesContainer",
    "HasTag",
    "Documents", "Demonstrates", "HasAnchor", "LinksTo", "HasTutorial",
    "Covers",
    "FixedIn", "MentionsTicket",
    "PartOfFeature", "HasFeature", "RelatesToFeature",
    # files
    "ContainsFile", "DefinedIn", "IsImplementedIn", "IsTestedIn", "Imports",
    # semtypes
    "DetectsSemtype", "ConsumesSemtype", "ProducesSemtype", "RequiresColumnTag",
    # JS API + TS structure
    "HasMethod", "HasFunction", "HasClass", "ExtendsClass", "ImplementsInterface",
    "InNamespace", "DelegatesTo", "TypedFor", "DefinesSemtype", "DefinedByFunction",
    # JS API usage
    "ImportsNamespace", "UsesApiClass", "UsesApiEnum",
    "CallsDapiEndpoint", "SubscribesToEvent", "UsesUiComponent",
]
