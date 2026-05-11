# Datagrok knowledge graph — agent entry point

**This file is for AI agents (Claude Code, Copilot, etc.) working in the `public` repo. Read it before you grep.**

There is a queryable knowledge graph of every plugin, function, script, query, library, doc page, tutorial, changelog entry, and Jira ticket in this repository. It lives at `.kg/` (inside the `public` git submodule, so it ships with the open-source repo). Stored file paths are anchored at `./` — e.g. a graph result `packages/Bio/src/foo.ts` is on disk at `packages/Bio/src/foo.ts`. Use it as your **first lookup** when you need to:

- Find which plugin owns a feature: "what implements activity cliffs?"
- Find documentation for a function: "where is `searchSubstructure` documented?"
- Find related code: "what other features does Bio depend on?"
- Find coverage gaps: "which Charts viewers have help docs but no `helpUrl` in code?"
- Trace a Jira ticket to commits, packages, and docs.
- Cluster work that touches a Feature.

Use it **instead of** filesystem grep when those questions can be answered structurally — the graph already knows.

## How to query — use `qq.py` (auto-starts a persistent server)

Each `query.py` invocation pays ~1.5s of Python + Kuzu cold-start. **Use [`qq.py`](qq.py) instead** — it talks to a long-running HTTP server (`kg_server.py`) on `127.0.0.1:7475`. First call starts the server in the background; every subsequent call returns in ~50ms (Kuzu queries themselves are 1-5ms).

```powershell
# First call: spawns the server in the background, then queries
.kg/.venv/Scripts/python.exe .kg/qq.py "MATCH (p:Package {name:'Chem'})-[:HAS_FEATURE]->(f:Feature) RETURN f.name LIMIT 10"

# All subsequent calls hit the running server directly
.kg/.venv/Scripts/python.exe .kg/qq.py "MATCH (n) RETURN label(n), count(n) ORDER BY 2 DESC LIMIT 5"

# Server status
.kg/.venv/Scripts/python.exe .kg/qq.py                    # show node/edge counts + uptime

# Raw JSON (for piping into another script)
.kg/.venv/Scripts/python.exe .kg/qq.py --json "MATCH (f:Feature) WHERE f.name CONTAINS 'Activity' RETURN f.id, f.name"

# Restart (e.g. after a rebuild)
.kg/.venv/Scripts/python.exe .kg/qq.py --restart
.kg/.venv/Scripts/python.exe .kg/qq.py --stop             # release the DB lock
```

The server binds to **127.0.0.1 only** — it's never network-exposed. If you'd rather not auto-spawn, pass `--no-start` and qq will fail loudly.

`query.py` still works for one-offs, demos, and pre-server scripts — but for any iterative work, prefer `qq.py`.

Common queries are catalogued in [docs/QUERYING.md](docs/QUERYING.md). When unsure, start with `--demos`:

```powershell
.kg/.venv/Scripts/python.exe .kg/query.py --demos
```

## What the graph contains (current state)

| Layer | Entity kinds | Counts |
|---|---|---:|
| **Plugins** (`packages/`) | Package, RegisteredFunction, Script, DataQuery, DataConnection, ScriptEnvironment, DockerContainer, PackageProperty, WasmModule | 79 / 1073 / 545 / 508 / 29 / 11 / 15 / 29 / 0 |
| **Libraries** (`libraries/`) | Library, LibraryModule | 18 / 515 |
| **Docs** (`help/`, plugin READMEs, Tutorials package) | DocPage, HelpAnchor, Tutorial, TutorialTrack | 695 / 3089 / 22 / 7 |
| **Process** (CHANGELOGs, Jira) | ChangelogEntry, JiraTicket | 2920 / 306 |
| **Tests** (each `test()` block in `src/tests/`) | PackageTest | 1701 |
| **TS code** (per-file extraction) | TsClass, TsMethod, TsFunction, TsInterface, TsEnum, TsConstant | 1181 / 9325 / 2559 / 1256 / 282 / 232 |
| **Files** (any modeled root) | File | 4597 |
| **JS API** | JsApiNamespace, DapiEndpoint, JsEventStream, UiComponent, GeneratedBinding | 11 / 31 / 47 / 68 / 1985 |
| **Tags & semtypes** | Tag, SemanticType | 31 / 82 |
| **Synthetic** (LLM-derived) | Feature | 458 |

**Edges** (~95,000 total). The high-leverage ones for an agent:

| Edge | Subject → Object | What it answers |
|---|---|---|
| `IMPORTS` | File → File / LibraryModule / JsApiNamespace | per-file ES-module import graph (11K+ edges) |
| `DEFINED_IN` | TsClass/TsMethod/TsFunction/RegisteredFunction/PackageTest → File | "where is X defined" with line number |
| `DEFINED_BY_FUNCTION` | RegisteredFunction → TsMethod/TsFunction | bridges `package.g.ts` wrapper to its real impl |
| `HAS_TAG` | RegisteredFunction → Tag | `cellEditor`, `app`, `panel`, etc. — synthesized from both `tags[]` AND `meta.role` |
| `REQUIRES_COLUMN_TAG` | RegisteredFunction → SemanticType | parsed from `meta.columnTags` (`quality=X, units=Y`) — what column shape a function dispatches for |
| `USES_API_CLASS` / `USES_API_ENUM` / `USES_UI_COMPONENT` / `CALLS_DAPI_ENDPOINT` / `SUBSCRIBES_TO_EVENT` | **Package + File + TsMethod + TsFunction** → js-api target | per-method usage of every JS API surface (multi-granularity since 2026-05) |
| `EXTENDS_CLASS` / `IMPLEMENTS_INTERFACE` | plugin TsClass → js-api TsClass / TsInterface | DG.X aliases resolved — "all subclasses of DG.GridCellRenderer" |
| `DETECTS_SEMTYPE` | RegisteredFunction → SemanticType | from `meta.role: semTypeDetector` AND `detectors.js` `detect<X>()` methods |
| `COVERS` | PackageTest → RegisteredFunction | best-effort name-match coverage |
| `PART_OF_FEATURE` / `HAS_FEATURE` / `IS_IMPLEMENTED_IN` / `IS_TESTED_IN` | LLM-derived Feature edges | feature-level rollups |

What's **NOT** yet in the graph (fall back to file reads):

- Dart client/server code (`core/`)
- Git commits (other than via changelog `MENTIONS_TICKET`)
- Live Jira state (status/assignee/etc — only the ticket-key node exists)
- Datagrok server runtime state (use `grok s` for that)
- Function-body call graph (`RegisteredFunction → RegisteredFunction` direct calls — only cross-package `Calls` via `grok.functions.call`)
- Runtime preconditions (`if (col.semType !== X) throw`, `seqHelper.getSeqHandler(col)`) — partially covered by `REQUIRES_COLUMN_TAG` for declared `meta.columnTags`, but body-level guards aren't extracted
- ApiTest / PlaywrightScenario nodes (extractor pending — but `PackageTest` covers `src/tests/`)

If you need any of the above, read files directly. The graph will surface what's there but never lies about what isn't.

## What you can ask now — query patterns by intent

For the full cookbook see [`docs/QUERYING.md`](docs/QUERYING.md). The patterns below cover the highest-leverage agent questions.

### "Show me a donor pattern for X"

> Find all `cellEditor`s with their column-scope and resolved impl file:

```cypher
MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'cellEditor'})
OPTIONAL MATCH (rf)-[:REQUIRES_COLUMN_TAG]->(s:SemanticType)
OPTIONAL MATCH (rf)-[:DEFINED_BY_FUNCTION]->(:TsMethod)-[:DEFINED_IN]->(f:File)
RETURN rf.id, rf.package_id, collect(DISTINCT s.name) AS column_scope, f.relative_path;
```

`HAS_TAG` is synthesized from both `tags[]` AND `meta.role`, so this captures both styles. Replace `cellEditor` with `panel` / `cellRenderer` / `fileViewer` / `valueEditor` etc. for the analogous question.

### "Which method uses X JS API surface?"

> Methods that open a `ui.dialog`:

```cypher
MATCH (m:TsMethod)-[u:USES_UI_COMPONENT]->(c:UiComponent {name:'dialog'})
RETURN m.id, u.use_count ORDER BY u.use_count DESC LIMIT 20;
```

> Methods that subscribe to `grok.events.onTableAdded`:

```cypher
MATCH (m:TsMethod)-[s:SUBSCRIBES_TO_EVENT]->(:JsEventStream {name:'onTableAdded'})
RETURN m.id, s.subscribe_count;
```

> Top files that call `grok.dapi.docker`:

```cypher
MATCH (f:File)-[c:CALLS_DAPI_ENDPOINT]->(:DapiEndpoint {name:'docker'})
RETURN f.relative_path, c.call_count ORDER BY c.call_count DESC LIMIT 10;
```

The same `USES_API_CLASS` / `USES_API_ENUM` / `USES_UI_COMPONENT` / `CALLS_DAPI_ENDPOINT` / `SUBSCRIBES_TO_EVENT` edges are emitted at four levels — `Package`, `File`, `TsMethod`, `TsFunction` — pick whichever granularity your question needs.

### "Where is X defined?" (TS code, not just RegisteredFunction)

```cypher
MATCH (c:TsClass {name:'MacromoleculeSequenceCellRenderer'})-[d:DEFINED_IN]->(f:File)
RETURN c.id, f.relative_path, d.line_start;
```

Same for `TsMethod`, `TsFunction`, `TsInterface`, `TsEnum`, `TsConstant`.

### "Which subclasses extend a JS API base class?"

```cypher
MATCH (sub:TsClass)-[:EXTENDS_CLASS]->(:TsClass {name:'GridCellRenderer', is_jsapi:true})
RETURN sub.id, sub.namespace LIMIT 20;
```

Plugin classes that say `extends DG.GridCellRenderer` resolve to the canonical js-api class via the alias map.

### "Which files import a specific symbol?"

> Per-file granularity for `getHelmHelper`:

```cypher
MATCH (f:File)-[i:IMPORTS]->(target)
WHERE 'getHelmHelper' IN i.imported_symbols
RETURN f.relative_path, target.id, i.imported_symbols LIMIT 10;
```

> What does this file actually import?

```cypher
MATCH (:File {id:'file:packages/Bio/src/utils/cell-renderer.ts'})-[i:IMPORTS]->(target)
RETURN target.id, i.imported_symbols, i.is_type_only, i.import_count;
```

> Reverse — who imports this file?

```cypher
MATCH (other:File)-[i:IMPORTS]->(:File {id:'file:packages/Bio/src/utils/seq-helper/seq-helper.ts'})
RETURN other.relative_path, i.imported_symbols ORDER BY other.relative_path;
```

### "What semType does a function dispatch for?" (donor verification)

```cypher
MATCH (rf:RegisteredFunction {id:'func:Helm:editMoleculeCell'})-[:REQUIRES_COLUMN_TAG]->(s:SemanticType)
RETURN s.name;
```

Returns `Macromolecule` and `helm` — the `quality=` and `units=` clauses from `meta.columnTags`. **Use this before delegating to a function from another package** to confirm preconditions match.

### "Test coverage at the test-block level"

> Every `test()` in Bio matching "renderer":

```cypher
MATCH (t:PackageTest)
WHERE t.package_id = 'pkg:Bio'
  AND (toLower(t.name) CONTAINS 'render' OR toLower(t.category) CONTAINS 'render')
RETURN t.category, t.name, t.file_path LIMIT 20;
```

> Tests that cover a specific function (best-effort name match):

```cypher
MATCH (t:PackageTest)-[:COVERS]->(:RegisteredFunction {name:'toAtomicLevel'})
RETURN t.category, t.name, t.file_path;
```

### "Feature → all the way down to JS API usage"

This is the chain that was previously broken:

```cypher
MATCH (feat:Feature)<-[:PART_OF_FEATURE]-(rf:RegisteredFunction)
      -[:DEFINED_BY_FUNCTION]->(m:TsMethod)
      -[u:USES_UI_COMPONENT]->(c:UiComponent {name:'dialog'})
RETURN feat.name, sum(u.use_count) AS dialog_calls
ORDER BY dialog_calls DESC LIMIT 10;
```

Returns features like *Substructure Search*, *Molecule Sketcher*, *Hierarchical Clustering* with their `ui.dialog` usage counts.

## Reliability — when to trust the graph

| Source of fact | Trust |
|---|---|
| Annotation comments (`//name: //tags: //help-url:`) | Very high — direct parse |
| `package.json`, JSON connections, YAML envs, Dockerfile dirs | Very high |
| CHANGELOG entries, Jira ticket IDs in commits | High |
| Library module imports (regex on TS source) | High |
| `Documents` edge from `//help-url:` annotation | High (when present — most are missing in code, see below) |
| `Documents` edge from `this.helpUrl=` regex match | Medium (class-name fuzzy match) |
| `Feature` clusters and `PART_OF_FEATURE` edges | Medium-high (LLM-derived, with `confidence` field) |
| `Documents` edges from `Tutorial.helpUrl` | Medium (URL normalization may miss anchors) |

Every node carries `description_provenance`; every edge carries `derived_by` and `confidence`. Filter on those if you want only deterministic facts.

## Verify before recommending

The graph is a snapshot from the last `py build.py` run. If a recent commit added/removed/renamed an entity, the graph won't know yet. Before recommending an answer to the user based on the graph:

- If the question is about **recent** state (e.g. "what was just added?"), use `git log` instead.
- If the answer names a specific file/function, **also check the file exists** before passing it to the user as actionable advice.

## How to extend

If you find a query the graph can't answer well — that's a signal the schema or an extractor needs work. Don't paper over it; capture the gap:

1. If a new entity kind is needed → add it in [schema/entities.py](schema/entities.py).
2. If a new relation is needed → add it in [schema/relations.py](schema/relations.py).
3. If extraction is missing → add an extractor in [extractors/](extractors/).
4. If it's an LLM-clusterable concept → add an enricher in [enrichers/](enrichers/).
5. Run `py build.py` and `py tests/test_kg_queries.py`.

Full guide: [docs/UPDATING.md](docs/UPDATING.md).

## Visual exploration

```powershell
.kg/.venv/Scripts/python.exe .kg/viz.py
```

Generates `.kg/kg.html` — open it in a browser. Default view shows packages + libraries. Pick a package to drill in; click any node for full detail; double-click to expand neighborhood.

## Don't break the graph

- Don't write to `kg.kuzu/` directly. It's a derived cache — edit JSONL, run `build.py`.
- Don't rename entity kinds without a migration. Many extractors and queries will hard-break.
- Don't commit `kg.kuzu/` (it's gitignored). Do commit `data/**.jsonl` so the graph travels with the repo.

## Tests

```powershell
.kg/.venv/Scripts/python.exe .kg/tests/test_kg_queries.py
```

33 tests covering schema integrity, coverage minima, connectivity, cross-layer traversals, feature invariants, and known coverage gaps. **Always run these after any schema or extractor change.**
