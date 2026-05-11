"""Concrete entity types.

Each subclass of `Entity` becomes a Kuzu node table. The `kind` field is
auto-populated from the subclass name (set in `__init_subclass__`).

Add new entity kinds here. Update extractors to emit them. Run `build.py`.

Layering note: many of these kinds (`Function`, `Test`, `DocPage`, `Commit`)
apply across plugins / core / libraries — the `source_layer` field
distinguishes them. Plugin-only kinds are clearly named (`Package`,
`RegisteredFunction`, `DataConnection`).
"""

from __future__ import annotations

from typing import Any, ClassVar

from pydantic import Field

from .base import Entity


def _set_kind(cls: type[Entity]) -> None:
    """Stamp the kind on a subclass and add it to the registry."""
    cls.model_fields["kind"].default = cls.__name__
    ENTITY_KINDS[cls.__name__] = cls


ENTITY_KINDS: dict[str, type[Entity]] = {}


class _EntityBase(Entity):
    """Internal helper that auto-stamps `kind` on subclassing."""

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        _set_kind(cls)


# ============================================================================
# PLUGINS LAYER (`packages/`)
# ============================================================================

class Package(_EntityBase):
    """A Datagrok plugin published as a package under `packages/`."""
    friendly_name: str | None = None
    version: str | None = None
    category: str | None = None
    author: str | None = None
    server_compatibility: str | None = None  # serverDependencies
    properties_count: int = 0                 # length of package.json["properties"]
    dependencies: list[str] = Field(default_factory=list)
    dev_dependencies: list[str] = Field(default_factory=list)


class Library(_EntityBase):
    """Shared TypeScript library under `libraries/`, published as
    `@datagrok-libraries/<name>` and consumed by packages."""
    npm_name: str | None = None              # @datagrok-libraries/<x>
    version: str | None = None
    role: str | None = None                  # utility | domain-framework | compute-framework | ui-component | domain-wrapper | feature | test-harness
    main: str | None = None                  # package.json main entrypoint
    is_platform_agnostic: bool = False       # True if no datagrok-api dep (e.g. chem-meta)
    ts_file_count: int = 0


class LibraryModule(_EntityBase):
    """A `.ts` file inside a library's `src/`. Referenced by package imports
    via deep paths (`@datagrok-libraries/ml/src/typed-metrics/consts`)."""
    library_id: str
    relative_path: str                       # e.g. 'src/typed-metrics/consts.ts'


class RegisteredFunction(_EntityBase):
    """A function exposed to the platform via the `//name: //tags: //meta.role:`
    annotation block (or `@grok.decorators.func` decorator). Lives in
    `package.g.ts` (canonical, generated form) or `package.ts`."""
    package_id: str | None = None            # owning Package id
    friendly_name: str | None = None
    role: str | None = None                  # app | viewer | widget | panel | semTypeDetector | cellRenderer | fileViewer | fileHandler | searchProvider | valueEditor | dashboard | transform | function
    tags: list[str] = Field(default_factory=list)
    inputs: list[dict[str, Any]] = Field(default_factory=list)   # [{name, type, options}]
    outputs: list[dict[str, Any]] = Field(default_factory=list)
    meta: dict[str, Any] = Field(default_factory=dict)            # //meta.* k/v
    help_url: str | None = None
    language: str = "ts"                                          # ts | js | py | r | sql | julia


class Script(_EntityBase):
    """A standalone script (Python/R/JS/Julia/SQL) under a package's
    `scripts/` folder, exposed to the platform via comment header."""
    package_id: str | None = None
    language: str                                                 # python | r | js | julia
    inputs: list[dict[str, Any]] = Field(default_factory=list)
    outputs: list[dict[str, Any]] = Field(default_factory=list)
    environment: str | None = None                                # //environment: name


class DataQuery(_EntityBase):
    """A named SQL query under a package's `queries/` folder."""
    package_id: str | None = None
    connection_name: str | None = None
    inputs: list[dict[str, Any]] = Field(default_factory=list)
    tags: list[str] = Field(default_factory=list)
    cache: bool = False


class DataConnection(_EntityBase):
    """A connection definition under `connections/`."""
    package_id: str | None = None
    data_source: str | None = None                                # Postgres, REST, etc.


class ScriptEnvironment(_EntityBase):
    """A conda/pip environment YAML under `environments/`."""
    package_id: str | None = None
    channels: list[str] = Field(default_factory=list)
    deps: list[str] = Field(default_factory=list)


class DockerContainer(_EntityBase):
    """A Dockerfile under a package's `dockerfiles/<service>/`."""
    package_id: str | None = None
    cpu: str | None = None
    memory: str | None = None
    gpu: bool = False
    on_demand: bool = False


class WasmModule(_EntityBase):
    """A WebAssembly artifact bundled with a package."""
    package_id: str | None = None
    bytes: int | None = None


class ChangelogEntry(_EntityBase):
    """One bullet under a `## <version>` section in `CHANGELOG.md`."""
    package_id: str | None = None
    version: str | None = None                                    # "v.next" or semver
    released_at: str | None = None                                # YYYY-MM-DD if versioned
    ticket_id: str | None = None                                  # GROK-NNNNN if present


class PackageProperty(_EntityBase):
    """A user-configurable setting declared in package.json `properties[]`."""
    package_id: str | None = None
    type: str | None = None
    choices: list[str] = Field(default_factory=list)
    default_value: Any = None


# ============================================================================
# DOCS LAYER (`help/`, `public/docusaurus/`, plugin READMEs)
# ============================================================================

class DocPage(_EntityBase):
    """A markdown documentation page."""
    title: str | None = None
    audience: str | None = None                                   # user | developer | admin | internal
    doc_kind: str | None = None                                   # guide | reference | tutorial | release-note | concept
    public: bool = True                                           # vs `_internal/`
    frontmatter: dict[str, Any] = Field(default_factory=dict)


class HelpAnchor(_EntityBase):
    """A named heading inside a DocPage (cross-reference target)."""
    doc_id: str
    anchor: str                                                   # slugified
    level: int = 2


class Tutorial(_EntityBase):
    """An interactive tutorial — a class extending the Tutorial base in the
    Tutorials package (or in a plugin's `tutorials/`). Drives users through
    a sequence of UI steps demonstrating a feature."""
    package_id: str | None = None
    class_name: str | None = None
    track: str | None = None                  # owning TutorialTrack name
    help_url: str | None = None               # normalized — points at DocPage/HelpAnchor
    prerequisites: list[str] = Field(default_factory=list)
    step_count: int | None = None


class TutorialTrack(_EntityBase):
    """A grouping of Tutorials registered via `new Track(name, [tutorials], helpUrl)`."""
    package_id: str | None = None
    help_url: str | None = None


# ============================================================================
# PROCESS LAYER (git, jira, jenkins, builds, tests)
# ============================================================================

class JiraTicket(_EntityBase):
    """A GROK-N ticket in JIRA."""
    ticket_key: str                                               # GROK-12345
    status: str | None = None
    assignee: str | None = None
    reporter: str | None = None
    client: str | None = None
    summary: str | None = None


class Commit(_EntityBase):
    """A git commit (in either repo)."""
    sha: str
    repo: str                                                     # "reddata" | "public"
    author: str | None = None
    committed_at: str | None = None
    subject: str | None = None
    ticket_id: str | None = None                                  # parsed GROK-N from message


# ============================================================================
# TESTS
# ============================================================================

class ApiTest(_EntityBase):
    """An ApiTests test() block under `packages/ApiTests/`."""
    category: str | None = None
    test_name: str | None = None
    owner: str | None = None
    stress_test: bool = False


class ApiSample(_EntityBase):
    """A canonical usage sample under `packages/ApiSamples/scripts/`."""
    help_url: str | None = None


class PlaywrightScenario(_EntityBase):
    """A browser e2e test under `public/playwright-public/`."""
    manual_scenario_ref: str | None = None                        # filename from leading comment


class PackageTest(_EntityBase):
    """A `test()` block inside `packages/<Pkg>/src/tests/**/*.ts`, run via
    `grok test`. Each `test('name', async () => { ... })` becomes one node;
    its enclosing `category('Cat', ...)` is the `category` field."""
    package_id: str | None = None
    category: str | None = None
    test_name: str | None = None
    is_skipped: bool = False
    is_only: bool = False
    file_path: str | None = None     # POSIX repo-rel; for click-through


# ============================================================================
# TAGS — first-class strings (not LLM)
# ============================================================================

class Tag(_EntityBase):
    """A function tag — `cellEditor`, `app`, `panel`, `init`, `autostart`,
    `cellRenderer`, etc. Promoted from `RegisteredFunction.tags[]`. Lets
    queries like `(:Tag {name:'cellEditor'})<-[:HAS_TAG]-(:RegisteredFunction)`
    enumerate without scanning every RF's tag array."""
    function_count: int = 0


# ============================================================================
# SYNTHETIC LAYER (LLM-derived; populated by enrichers)
# ============================================================================

class Feature(_EntityBase):
    """A user-facing feature, derived by clustering names and references
    across help docs, tests, code, and changelogs. Always `LLM` provenance.

    A Feature can live anywhere — not just `<viewer>/features/`. Each
    Feature carries the file refs and doc refs of its constituents.
    """
    cluster_score: float = 0.0                                    # how cohesive the cluster is
    member_ids: list[str] = Field(default_factory=list)           # entity ids that comprise it

    # Classification (populated deterministically by feature_classify.py)
    is_user_facing: bool = False        # True if any member is app/viewer/widget/panel/transform/etc.
    interaction_kind: str | None = None # 'app'|'viewer'|'widget-panel'|'transform'|'data-import'|
                                        # 'api'|'infra'|'data-query'|'script'|'mixed'


# ============================================================================
# CODE / FILESYSTEM (cross-cutting)
# ============================================================================

class File(_EntityBase):
    """A source file in the repo. Modeled at file granularity (not per line).
    Lets you ask 'where is this Feature implemented?', 'which files have no
    code entity?', 'how many test files cover this Function?'.
    """
    package_id: str | None = None
    library_id: str | None = None
    language: str | None = None       # ts | js | py | r | sql | md | json | yaml | dockerfile | css | other
    size_bytes: int = 0
    line_count: int = 0
    is_test: bool = False             # heuristic: path/filename signals test code
    is_generated: bool = False        # *.g.ts, *.api.g.ts, *.xcmd.g.dart
    relative_path: str | None = None  # repo-relative POSIX (also stored as `paths[0]`)


class SemanticType(_EntityBase):
    """A Datagrok semantic type identifier (e.g. 'Molecule', 'ChEMBL_ID',
    'Macromolecule'). Detected by `semTypeDetector` functions; consumed by
    other functions whose inputs declare `{semType: 'X'}`. Cross-cutting —
    a semtype produced by Chem may be consumed by Bio, etc.
    """
    semtype_value: str | None = None       # canonical string as used in code
    detected_by_count: int = 0             # how many detectors yield this
    consumed_by_count: int = 0             # how many functions filter on this


# ============================================================================
# JS API LAYER (`js-api/`) — used by every plugin
# ============================================================================

class JsApiNamespace(_EntityBase):
    """A logical namespace exposed by the JS API: 'grok', 'ui', 'DG' at
    the top, plus nested ones like 'grok.dapi', 'grok.events', 'DG.chem',
    'ui.dialog'. Plugins import these and reference members through them."""
    dotted_path: str | None = None         # e.g. 'grok.dapi'
    parent: str | None = None              # parent namespace id, or None for top-level


class TsClass(_EntityBase):
    """A TypeScript class declared in `public/{js-api,libraries,packages}/`.
    `source_layer` distinguishes where it lives. Use the `is_jsapi` field
    when you specifically want the public API surface (the `DG.*` classes)."""
    namespace: str | None = None           # e.g. 'DG.chem', 'pkg:Chem', 'lib:utils'
    dotted_name: str | None = None         # 'DG.chem.Sketcher' or 'Chem.ScaffoldTreeViewer'
    is_abstract: bool = False
    extends_class: str | None = None       # parent class id
    method_count: int = 0
    is_jsapi: bool = False                 # in js-api/?
    is_exported: bool = True


class TsMethod(_EntityBase):
    """A public instance/static method or getter on a TsClass."""
    class_id: str | None = None
    is_static: bool = False
    is_getter: bool = False
    return_type: str | None = None
    parameter_count: int = 0
    delegates_to: str | None = None        # name of the IDartApi `grok_*` stub if delegated


class TsFunction(_EntityBase):
    """A top-level `export function` declaration in any TS file (not a class
    method). For platform-registered functions (the `//name:`-annotated
    ones), `RegisteredFunction` is the richer node — but plenty of helper
    functions exist that aren't registered."""
    namespace: str | None = None
    return_type: str | None = None
    parameter_count: int = 0
    is_async: bool = False
    is_jsapi: bool = False


class TsEnum(_EntityBase):
    """An `export enum` declaration."""
    namespace: str | None = None
    value_count: int = 0
    values: list[str] = Field(default_factory=list)
    is_jsapi: bool = False


class TsConstant(_EntityBase):
    """An `export const` holding a string-keyed plain object / registry
    (e.g. `SEMTYPE`, `UNITS`). Plain primitive constants are skipped."""
    namespace: str | None = None
    value_count: int = 0
    is_jsapi: bool = False


class TsInterface(_EntityBase):
    """An `export interface` or `export type` declaration."""
    namespace: str | None = None
    member_count: int = 0
    is_type_alias: bool = False
    is_jsapi: bool = False


class DapiEndpoint(_EntityBase):
    """A typed `HttpDataSource<T>` accessor on `Dapi` — the typed REST clients.
    e.g. `grok.dapi.users` returns `UsersDataSource` typed for `User`."""
    accessor_path: str | None = None       # 'grok.dapi.users'
    data_source_class: str | None = None   # 'UsersDataSource'
    served_entity: str | None = None       # 'User' (the T in HttpDataSource<T>)


class JsEventStream(_EntityBase):
    """An RxJS Observable property on the `Events` class, accessed as
    `grok.events.onTableAdded` etc. ~115 streams."""
    accessor_path: str | None = None       # 'grok.events.onTableAdded'
    event_args_type: str | None = None     # 'EventData<DataFrameArgs>'
    event_id: str | None = None            # if backed by an EVENT_TYPE constant


class GeneratedBinding(_EntityBase):
    """A `grok_<DartClass>_<DartMember>` stub on the `IDartApi` interface —
    auto-generated bridge to Dart. ~900 of these. They are how every JS API
    method actually crosses into Dart land at runtime."""
    dart_class: str | None = None          # parsed from name segments
    dart_member: str | None = None
    parameter_count: int = 0


class UiComponent(_EntityBase):
    """A factory function exported from `js-api/src/ui.ts` that
    produces a DOM element or Widget. e.g. `ui.button`, `ui.dialog`."""
    factory_name: str | None = None
    return_type: str | None = None
    parameter_count: int = 0


# ============================================================================
# Public re-export
# ============================================================================

__all__ = [
    "ENTITY_KINDS",
    # plugin
    "Package", "Library", "LibraryModule", "RegisteredFunction", "Script",
    "DataQuery", "DataConnection", "ScriptEnvironment", "DockerContainer",
    "WasmModule", "ChangelogEntry", "PackageProperty",
    # docs
    "DocPage", "HelpAnchor", "Tutorial", "TutorialTrack",
    # process
    "JiraTicket", "Commit",
    # tests
    "ApiTest", "ApiSample", "PlaywrightScenario", "PackageTest",
    # tags
    "Tag",
    # synthetic
    "Feature",
    # cross-cutting
    "File", "SemanticType",
    # TS code structure (any layer)
    "TsClass", "TsMethod", "TsFunction", "TsEnum", "TsConstant", "TsInterface",
    # JS API specific
    "JsApiNamespace", "DapiEndpoint", "JsEventStream", "GeneratedBinding", "UiComponent",
]
