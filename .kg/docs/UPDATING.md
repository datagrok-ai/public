# Updating the knowledge graph

This is the maintainer's guide. It explains when to rebuild, how to add new entity kinds, how to add or modify extractors and enrichers, how schema migrations work, and the test discipline that keeps the graph honest.

## Mental model (recap)

```
schema/      Pydantic models — single source of truth for shape.
extractors/  Deterministic parsers (annotations, package.json, regex on TS, …).
enrichers/   LLM-driven enrichment (Feature clustering, doc↔code bridges).
data/        CANONICAL JSONL — committed to git, the truth.
kg.kuzu/     Materialized graph DB — gitignored, rebuilt from JSONL.
build.py     extractors + DDL + load → Kuzu.
enrich.py    prepare/apply lifecycle for LLM enrichers.
viz.py       generates kg.html.
```

The DB is a cache. Edit JSONL (or extractors that produce it), rebuild. Never write to `kg.kuzu/` directly.

## Rebuild commands

```powershell
.kg/.venv/Scripts/python.exe .kg/build.py                 # full extract → JSONL → Kuzu
.kg/.venv/Scripts/python.exe .kg/build.py --packages Chem # restrict plugins-layer to one
.kg/.venv/Scripts/python.exe .kg/build.py --extractors library  # run just one extractor
.kg/.venv/Scripts/python.exe .kg/build.py --no-db         # JSONL only
.kg/.venv/Scripts/python.exe .kg/build.py --clean         # wipe data/ and kg.kuzu/ first
```

After build:

```powershell
.kg/.venv/Scripts/python.exe .kg/tests/test_kg_queries.py
.kg/.venv/Scripts/python.exe .kg/viz.py                   # regenerate kg.html
```

## When to rebuild

| Trigger | Action |
|---|---|
| Schema change (new entity kind, renamed field) | Full `--clean` rebuild + tests |
| New extractor added | Run `--extractors <name>` then full `build.py` to refresh dependents |
| Plugin source change in `packages/` | Re-run plugin extractors (most are idempotent and merge by id) |
| New help docs | Re-run `docs` extractor |
| Updated CHANGELOGs after release | Re-run `ts_plugin_changelog` |
| LLM clustering produced new outputs | `enrich.py apply --enricher feature_clustering` then `build.py` |
| Test failure after change | Inspect failure, fix the root cause, never weaken the assertion blindly |

## Adding a new entity kind

Concrete example: imagine you want to add `Connector` (Java JDBC connectors under `public/connectors/`).

### 1. Define the Pydantic model

In [schema/entities.py](../schema/entities.py):

```python
class Connector(_EntityBase):
    """A Java JDBC connector under `public/connectors/`."""
    db_family: str | None = None      # postgres | mysql | oracle | mongo | …
    java_class: str | None = None     # e.g. "grok.connect.providers.PostgresProvider"
    supports_streaming: bool = False
```

Then add `"Connector"` to the `__all__` list at the bottom.

The base class auto-stamps `kind = "Connector"` and registers it in `ENTITY_KINDS`. The DDL deriver picks it up automatically — no DDL to write by hand.

### 2. Define the relation(s)

In [schema/relations.py](../schema/relations.py):

```python
class HasConnector(_RelBase):
    """Package owns / ships a Connector (rare; usually connectors are global)."""
    SUBJECT_KINDS = ["Package"]
    OBJECT_KINDS = ["Connector"]


class ProvidesDataSource(_RelBase):
    """A Connector provides a `dataSource` value referenced by DataConnection."""
    SUBJECT_KINDS = ["Connector"]
    OBJECT_KINDS = ["DataConnection"]
```

Add to `__all__`.

### 3. Write the extractor

Create [extractors/connectors.py](../extractors/connectors.py):

```python
from pathlib import Path
from schema.base import Provenance, SourceLayer
from schema.entities import Connector
from schema.relations import ProvidesDataSource
from . import _common as c

EXTRACTOR_NAME = "connectors"

def run(repo_root: Path, bundle, package_filter=None):
    root = repo_root / "public" / "connectors"
    if not root.is_dir():
        return
    n = 0
    for provider_dir in sorted(root.iterdir()):
        if not provider_dir.is_dir():
            continue
        # …pick out the Java class and dataSource string…
        cid = f"connector:{provider_dir.name}"
        bundle.add(Connector(
            id=cid, name=provider_dir.name,
            source_layer=SourceLayer.CORE_SHARED,    # or whichever fits
            db_family=…, java_class=…,
            paths=[c.file_ref(provider_dir, role="root")],
        ))
        n += 1
    print(f"[{EXTRACTOR_NAME}] {n} connectors")
```

### 4. Register the extractor

In [build.py](../build.py), add to `EXTRACTORS`:

```python
EXTRACTORS = {
    …,
    "connectors": "extractors.connectors",
}
```

Order matters when one extractor's nodes are referenced as edge endpoints by another. Library runs before plugin so cross-package `DependsOn` edges land. Add `connectors` near the top if other extractors will reference its IDs.

### 5. Build, test, visualize

```powershell
.kg/.venv/Scripts/python.exe .kg/build.py --extractors connectors
.kg/.venv/Scripts/python.exe .kg/build.py            # full rebuild to load into Kuzu
.kg/.venv/Scripts/python.exe .kg/tests/test_kg_queries.py
.kg/.venv/Scripts/python.exe .kg/viz.py              # regen viewer (now shows Connectors)
```

If you add new visual style for a kind, edit `KIND_STYLE` in [viz.py](../viz.py).

### 6. Add coverage tests

In [tests/test_kg_queries.py](../tests/test_kg_queries.py):

```python
@test("coverage", "At least 25 Java connectors")
def test_connectors_count(conn):
    n = q_one(conn, "MATCH (c:Connector) RETURN count(c);")
    assert n >= 25, f"Only {n} connectors"

@test("crosslayer", "Every DataConnection's dataSource resolves to a known Connector")
def test_connection_resolves_to_connector(conn):
    bad = q(conn, """
        MATCH (dc:DataConnection)
        WHERE dc.data_source IS NOT NULL
          AND NOT EXISTS { MATCH (:Connector {db_family: dc.data_source}) }
        RETURN count(dc);
    """)
    assert bad[0][0] == 0, "DataConnections reference unknown connectors"
```

Tests double as the schema's executable spec. Adding tests at the same time as the extractor catches typos and missed cases on day one.

## Adding a new extractor for an existing kind

Most common case. Suppose you want to extract `Documents` edges from `ApiSamples` script headers (currently the `docs` extractor does this for `script_help_url` derivation, but only for `.js`):

1. Add (or edit) the relevant logic in [extractors/docs.py](../extractors/docs.py).
2. Run `--extractors docs` (or full build).
3. Verify edge counts go up (`MATCH ()-[r:DOCUMENTS {derivation:'script_help_url'}]->() RETURN count(r);`).
4. Update or add a test asserting the new count is non-zero.

Extractors are **idempotent** — they merge into the canonical JSONL by primary key. You can rerun any number of times without dupes. The orchestrator clears slices only on `--clean`.

## Adding a new LLM enricher

Concrete example: a `Concern` enricher that flags features touching security, performance, or accessibility.

### 1. Add the entity (if not present)

```python
# schema/entities.py
class Concern(_EntityBase):
    concern_type: str | None = None   # security | perf | a11y | privacy
    severity: str | None = None        # high | medium | low
```

### 2. Add the relation

```python
# schema/relations.py
class FlaggedAs(_RelBase):
    SUBJECT_KINDS = ["Feature", "RegisteredFunction"]
    OBJECT_KINDS = ["Concern"]
```

### 3. Write the enricher

Follow the `prepare → external LLM → apply` contract used by [enrichers/feature_clustering.py](../enrichers/feature_clustering.py):

```python
# enrichers/concern_flagging.py
NAME = "concern_flagging"
ENRICH_DIR = ROOT / "data" / "enrichment" / NAME
INPUTS_DIR = ENRICH_DIR / "inputs"
OUTPUTS_DIR = ENRICH_DIR / "outputs"

def prepare(packages=None):
    # For each Feature, write an input JSON with its description, member names,
    # related changelog entries. Write data/enrichment/concern_flagging/inputs/<feature>.json
    ...

PROMPT_BODY = "..."   # written to PROMPT.md alongside inputs/

def apply(packages=None):
    # Read outputs/<feature>.json, validate, emit Concern + FlaggedAs JSONL
    ...
```

### 4. Register in `enrich.py`

```python
ENRICHERS = {
    "feature_clustering": "enrichers.feature_clustering",
    "concern_flagging":   "enrichers.concern_flagging",
}
```

### 5. Run

```powershell
.kg/.venv/Scripts/python.exe .kg/enrich.py prepare --enricher concern_flagging --packages Chem
# Then in this Claude Code session: spawn Agent tool calls per input file
# (one per feature/package), each agent reads input + PROMPT.md, writes output JSON.
.kg/.venv/Scripts/python.exe .kg/enrich.py apply --enricher concern_flagging
.kg/.venv/Scripts/python.exe .kg/build.py
```

**Important:** if you don't have an Anthropic API key (e.g. Claude Max subscription only), do not write enricher code that calls the Anthropic SDK directly — it will fail. Use the Claude Code `Agent` tool to do the LLM work, one agent per input file. The lifecycle is `prepare()` writes per-unit input JSON → orchestrator launches one Agent per input → `apply()` validates outputs and merges into JSONL.

## Currently registered enrichers

| Name | Purpose | Output edges |
|---|---|---|
| `feature_clustering` | LLM clusters every Package's entities into user-facing Features | `Feature` nodes, `HAS_FEATURE`, `PART_OF_FEATURE`, `RELATES_TO_FEATURE` |
| `doc_linking` | LLM links DocPages under `help/` to the Packages, Features, and Functions they describe | `DOCUMENTS` edges with `derivation: "llm_doc_linking"` |

Run them like this:

```powershell
.kg/.venv/Scripts/python.exe .kg/enrich.py prepare --enricher doc_linking
# (Then in this Claude Code session, spawn one Agent per slice, pointed at
# data/enrichment/doc_linking/inputs/<slice>.json with the contract in PROMPT.md)
.kg/.venv/Scripts/python.exe .kg/enrich.py apply --enricher doc_linking
.kg/.venv/Scripts/python.exe .kg/qq.py --stop          # release the DB lock
.kg/.venv/Scripts/python.exe .kg/build.py              # rebuild Kuzu
```

The `apply()` step is idempotent and merges by `(from_id, to_id)`. It will NOT overwrite a deterministic Documents edge (`derivation: "help_url_annotation"`, `"viewer_class_help_url"`, `"readme"`, etc.) — those win over LLM-derived ones. LLM edges with higher confidence overwrite older LLM edges of the same key.

## Schema migrations

The schema lives in code (Pydantic), so most changes are git diffs:

| Change | What to do |
|---|---|
| **Add a field** to an existing entity | Edit the model, rebuild. Old JSONL rows missing the field default to None — this is fine. |
| **Rename a field** | Edit the model, write a migration script in `schema/migrations/0001_rename_X_to_Y.py` that rewrites every JSONL row. Run it once. Bump `schema/version.txt`. |
| **Drop a field** | Edit the model, optionally delete the column from JSONL via a migration. |
| **Rename an entity kind** | Migration must rewrite every JSONL row's `kind`, every edge's `from_id` / `to_id` if those embed the kind in IDs, plus the JSONL filename. |
| **Add an entity kind** | Just add the model — no migration needed. |
| **Drop an entity kind** | Delete the model + the JSONL slice + any dangling edges. |
| **Add a relation kind** | Just add the model — no migration needed. |
| **Drop a relation kind** | Delete the model + the JSONL slice. |

Always run `tests/test_kg_queries.py` after a migration. The test suite is the regression net.

## ID conventions

Stable IDs are essential because they're how every edge points to its endpoints. Conventions live in [extractors/_common.py](../extractors/_common.py):

| Kind | ID pattern | Example |
|---|---|---|
| Package | `pkg:<Folder>` | `pkg:Chem` |
| Library | `lib:<short-name>` | `lib:utils` |
| LibraryModule | `lib_mod:<lib>:<rel-path>` | `lib_mod:ml:src/typed-metrics/consts.ts` |
| RegisteredFunction | `func:<Pkg>:<name>` | `func:Chem:scaffoldTreeViewer` |
| Script | `script:<Pkg>:<name>` | `script:Chem:descriptors` |
| DataQuery | `query:<Pkg>:<name>` | `query:Chembl:patternSimilaritySearch` |
| DataConnection | `conn:<Pkg>:<name>` | `conn:Chembl:Chembl` |
| ScriptEnvironment | `env:<Pkg>:<name>` |  |
| DockerContainer | `docker:<Pkg>:<svc>` |  |
| WasmModule | `wasm:<Pkg>:<file>` |  |
| ChangelogEntry | `clog:<Pkg>:<version>:<idx>` |  |
| PackageProperty | `prop:<Pkg>:<name>` |  |
| DocPage | `doc:<canonical-url>` | `doc:/help/visualize/viewers/scatter-plot` |
| HelpAnchor | `anchor:<doc-url>#<slug>` |  |
| Tutorial | `tutorial:<Pkg>:<class>` |  |
| TutorialTrack | `track:<Pkg>:<slug>` |  |
| JiraTicket | `jira:GROK-NNNNN` |  |
| Commit | `commit:<repo>:<sha12>` |  |
| Feature | `feature:<Pkg>:<slug>` |  |

If you change an ID convention, you must (1) update `_common.py`, (2) write a migration to rewrite every JSONL row including edge endpoints, (3) bump `schema/version.txt`.

## CI / pre-commit

Recommended (not yet wired up):

- Run `tests/test_kg_queries.py` on every push to `master`.
- Pre-commit hook: if any `packages/*/package.json` changed, re-run `--extractors ts_plugin_package` before committing — keeps JSONL in lockstep with code.

## Visualizer maintenance

The visualizer is a single self-contained HTML, regenerated from the current data:

```powershell
.kg/.venv/Scripts/python.exe .kg/viz.py
.kg/.venv/Scripts/python.exe .kg/viz.py --exclude-kinds HelpAnchor,ChangelogEntry,LibraryModule,JiraTicket
```

By default the noisy high-cardinality kinds are excluded so the canvas stays usable. Pass `--exclude-kinds ""` to include everything (slow).

If you add a new entity kind, also add a style entry to `KIND_STYLE` in [viz.py](../viz.py).

## Persistent query server

`kg_server.py` keeps the Kuzu DB loaded so `qq.py` queries return in tens of milliseconds (vs ~1.5s cold-start for `query.py`). The server is auto-started by `qq.py` if it isn't running. **Important:** Kuzu's local file lock means **only one process can have the DB open at a time** — that includes `tests/test_kg_queries.py` and any `python -c "import kuzu; kuzu.Database(...)"` you run. If you're debugging:

```powershell
.kg/.venv/Scripts/python.exe .kg/qq.py --stop          # release the lock
.kg/.venv/Scripts/python.exe .kg/tests/test_kg_queries.py
.kg/.venv/Scripts/python.exe .kg/qq.py "..."           # auto-restart
```

After a `build.py` rebuild, the server's in-memory copy is stale until you restart it: `qq.py --restart`.

## Common gotchas

- **Edge endpoints not loaded**: Kuzu skips edges whose endpoints don't exist. The build prints `(<predicate>: N edges skipped)` per predicate. Most are legitimate (e.g. `DEPENDS_ON` to a third-party npm package we don't model). Investigate if the count grows after a change.
- **Empty lists become NULL**: A function with `tags: []` is stored as `tags = NULL` (Kuzu can't infer empty-list element type during parameter binding). Filter with `WHERE tags IS NOT NULL` before `size(tags)`.
- **Reserved Cypher words as identifiers**: `default`, `match`, `from`, `to`, `value`. Use backticks: `` `default` ``.
- **`kg.kuzu` is a *file* not a directory** in current Kuzu (≥ 0.10). The build script handles both — but if you write tooling that walks it, check `Path.is_file()` first.
- **WSL/native path mismatches** on Windows: `_common.py::repo_rel` always emits POSIX paths. Use those when constructing IDs that embed paths.
