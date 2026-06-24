# Querying the Datagrok knowledge graph

A cookbook of useful Cypher queries against the local Kuzu graph at `.kg/kg.kuzu/`. Run them via:

```powershell
.kg/.venv/Scripts/python.exe .kg/scripts/query.py "<your cypher>"
.kg/.venv/Scripts/python.exe .kg/scripts/query.py --demos        # canned diagnostic queries
```

Or from Python directly:

```python
import kuzu
db = kuzu.Database(".kg/kg.kuzu")
conn = kuzu.Connection(db)
res = conn.execute("MATCH (p:Package) RETURN p.name LIMIT 5;")
while res.has_next():
    print(res.get_next())
```

## Schema cheat-sheet

### Entity kinds (node tables)

```
Package                     Library                  LibraryModule
RegisteredFunction          Script                   DataQuery
DataConnection              ScriptEnvironment        DockerContainer
PackageProperty             WasmModule               ChangelogEntry
DocPage                     HelpAnchor
Tutorial                    TutorialTrack
JiraTicket                  Commit
ApiTest                     ApiSample                PlaywrightScenario
Feature                                  ← LLM-derived
File                        SemanticType             ← cross-cutting
```

Every node has these base fields: `id` (PK), `name`, `kind`, `source_layer`, `description`, `description_provenance`, `paths` (JSON of `[{path, line_start, line_end, role}]`), `extras` (JSON), `extracted_by`, `extracted_at`. Kind-specific fields live alongside (e.g. `RegisteredFunction.role`, `Library.is_platform_agnostic`).

### Relation kinds (rel tables)

| Predicate | Subject → Object | Notes |
|---|---|---|
| `EXPORTS` | Package → RegisteredFunction | Always present |
| `HAS_SCRIPT` / `HAS_QUERY` / `HAS_CONNECTION` / `HAS_ENVIRONMENT` / `HAS_CONTAINER` / `HAS_PROPERTY` / `HAS_CHANGELOG_ENTRY` / `HAS_TUTORIAL` | Package → … | Containment |
| `DEPENDS_ON` | Package/Library → Package/Library | npm-style; carries `semver_range`, `is_dev`, `is_test_harness` |
| `IMPORTS_FROM_MODULE` | Package/Library → LibraryModule | Carries `imported_symbols`, `import_count`, `is_type_only` |
| `USES_CONNECTION` | DataQuery → DataConnection | from `--connection:` SQL annotation |
| `REQUIRES_ENVIRONMENT` | Script → ScriptEnvironment | from `#environment:` script header |
| `USES_CONTAINER` | RegisteredFunction/Script → DockerContainer | (sparse) |
| `DOCUMENTS` | DocPage/HelpAnchor → RegisteredFunction/Script/DataQuery/Tutorial/Package/Feature | Carries `derivation` (`help_url_annotation` / `script_help_url` / `viewer_class_help_url` / `tutorial_help_url` / `readme`) |
| `LINKS_TO` | DocPage → DocPage/HelpAnchor | Markdown cross-link |
| `HAS_ANCHOR` | DocPage → HelpAnchor | h2/h3 inside a page |
| `MENTIONS_TICKET` | ChangelogEntry/Commit/DocPage → JiraTicket | `GROK-NNNNN` mentions |
| `FIXED_IN` | JiraTicket → Commit | (extractor pending) |
| `HAS_FEATURE` | Package → Feature | Primary owner; LLM-derived |
| `PART_OF_FEATURE` | RegisteredFunction/Script/DataQuery/DocPage/Tutorial/ApiTest/ApiSample/ChangelogEntry → Feature | Carries `weight` (0-1) and `role` (`core`/`doc`/`test`/`sample`/`related`) |
| `RELATES_TO_FEATURE` | Feature → Feature | Carries `relation_kind` (`subfeature`/`sibling`/`related`) |
| `CONTAINS_FILE` | Package/Library → File | Containment |
| `DEFINED_IN` | RegisteredFunction/Script/DataQuery/DocPage/Tutorial/… → File | Where the entity lives; carries `line_start`, `line_end`, `role` |
| `IS_IMPLEMENTED_IN` | Feature → File | Synthesized: file holds ≥1 implementation member of the Feature; `member_count` |
| `IS_TESTED_IN` | Feature → File | Same but file is `is_test=true` |
| `DETECTS_SEMTYPE` | RegisteredFunction → SemanticType | Function with `meta.role: semTypeDetector` + `meta.semType: X` |
| `CONSUMES_SEMTYPE` | RegisteredFunction → SemanticType | Function input has `{semType: X}` |
| `PRODUCES_SEMTYPE` | RegisteredFunction → SemanticType | Function output has `{semType: X}` |
| `CALLS` | Package → RegisteredFunction | Cross-package runtime call (`grok.functions.eval/call/find`, `DG.Func.byName/find`); carries `call_count`, `callee_spec` |
| `IMPORTS_FROM_PACKAGE` | Package → Package | Direct TS import — empty in this codebase by convention |
| `Demonstrates` / `Covers` | … → … | Used by Tutorial / ApiTest extractors |

## Recipe book

### 1. Find what implements a feature

```cypher
MATCH (p:Package {name:'Chem'})-[:HAS_FEATURE]->(f:Feature)
WHERE f.name CONTAINS 'Activity'
MATCH (f)<-[r:PART_OF_FEATURE]-(m)
RETURN f.name, label(m) AS kind, m.name, r.weight, r.role
ORDER BY r.weight DESC;
```

### 2. Find docs for a function

```cypher
MATCH (f:RegisteredFunction {name:'Recalculate Coordinates'})<-[d:DOCUMENTS]-(doc:DocPage)
RETURN doc.title, doc.extras, d.derivation;
```

### 3. Find functions with no docs (coverage gap)

```cypher
MATCH (p:Package {name:'Charts'})-[:EXPORTS]->(f:RegisteredFunction {role:'viewer'})
WHERE NOT EXISTS { MATCH (f)<-[:DOCUMENTS]-() }
RETURN f.name;
```

### 4. Find features with no doc members

```cypher
MATCH (p:Package {name:'Chem'})-[:HAS_FEATURE]->(f:Feature)
WHERE NOT EXISTS { MATCH (f)<-[:PART_OF_FEATURE]-(:DocPage) }
RETURN f.name;
```

### 5. Most-used libraries

```cypher
MATCH (p:Package)-[:DEPENDS_ON]->(l:Library)
RETURN l.name, l.role, count(p) AS users
ORDER BY users DESC LIMIT 15;
```

### 6. Which Bio modules does Chem actually import?

```cypher
MATCH (p:Package {name:'Chem'})-[i:IMPORTS_FROM_MODULE]->(m:LibraryModule)
WHERE m.library_id = 'lib:bio'
RETURN m.relative_path, i.imported_symbols, i.import_count;
```

### 7. All Jira tickets fixed by Chem changelog

```cypher
MATCH (p:Package {name:'Chem'})-[:HAS_CHANGELOG_ENTRY]->(e:ChangelogEntry)-[:MENTIONS_TICKET]->(t:JiraTicket)
RETURN DISTINCT t.ticket_key, e.version, e.released_at, e.description
ORDER BY e.released_at DESC LIMIT 30;
```

### 8. Find the package that owns a Jira ticket

```cypher
MATCH (t:JiraTicket {ticket_key:'GROK-18139'})<-[:MENTIONS_TICKET]-(e:ChangelogEntry)<-[:HAS_CHANGELOG_ENTRY]-(p:Package)
RETURN p.name, e.version, e.description;
```

### 9. Cross-package: who depends on this package as a devDependency?

```cypher
MATCH (p:Package)-[d:DEPENDS_ON {is_dev:true}]->(:Package {name:'PowerGrid'})
RETURN p.name;
```

### 10. SQL queries by connection

```cypher
MATCH (q:DataQuery)-[:USES_CONNECTION]->(c:DataConnection)
RETURN c.name, c.data_source, count(q) AS query_count
ORDER BY query_count DESC;
```

### 11. Scripts by language

```cypher
MATCH (s:Script)
RETURN s.language, count(s) AS n
ORDER BY n DESC;
```

### 12. Find a feature by fuzzy name (across all packages)

```cypher
MATCH (f:Feature)
WHERE toLower(f.name) CONTAINS 'substructure'
   OR toLower(f.id) CONTAINS 'substructure'
RETURN f.id, f.name, f.cluster_score;
```

### 13. Tutorials that target a specific help anchor

```cypher
MATCH (t:Tutorial)-[:DEMONSTRATES]->(target)
RETURN t.name, t.help_url, target.id;
```

(Note: `DEMONSTRATES` from Tutorials is currently sparse — populated from `Tutorial.helpUrl`. LLM enrichment to extract `grok.functions.call('Pkg:Foo')` from `_run()` bodies is planned.)

### 14. Orphan documentation (docs nothing links to)

```cypher
MATCH (d:DocPage)
WHERE NOT EXISTS { MATCH (d)<-[:LINKS_TO]-() }
  AND NOT EXISTS { MATCH (d)-[:DOCUMENTS]->() }
RETURN d.audience, count(d) AS orphan_docs ORDER BY orphan_docs DESC;
```

### 15. Stale `//help-url:` (URL doesn't resolve to any DocPage)

```cypher
MATCH (f:RegisteredFunction) WHERE f.help_url IS NOT NULL
  AND NOT EXISTS { MATCH (f)<-[:DOCUMENTS]-(:DocPage) }
RETURN f.name, f.help_url, f.package_id LIMIT 30;
```

### 16. Heavy library imports (potential coupling smell)

```cypher
MATCH (p:Package)-[i:IMPORTS_FROM_MODULE]->(m:LibraryModule)
WITH p, m.library_id AS lib, count(distinct m) AS modules
WHERE modules >= 10
RETURN p.name, lib, modules ORDER BY modules DESC;
```

### 17. WASM-using packages

```cypher
MATCH (p:Package) WHERE p.extras CONTAINS 'wasm'
RETURN p.name;
```

(`extras` is a JSON-string column; `CONTAINS 'wasm'` is a cheap text search on it.)

### 18. Functions per role

```cypher
MATCH (f:RegisteredFunction) WHERE f.role IS NOT NULL
RETURN f.role, count(f) AS n ORDER BY n DESC;
```

### 19. Sibling features within a package

```cypher
MATCH (p:Package {name:'Bio'})-[:HAS_FEATURE]->(f1:Feature)
      -[r:RELATES_TO_FEATURE {relation_kind:'sibling'}]->(f2:Feature)
RETURN f1.name, f2.name;
```

### 20. Walk: from a doc page → its feature → all sibling features → all their docs

```cypher
MATCH (d:DocPage)-[:DOCUMENTS]->()<-[:PART_OF_FEATURE]-(:Feature)
      -[:RELATES_TO_FEATURE]->(f2:Feature)<-[:PART_OF_FEATURE]-(d2:DocPage)
WHERE d.title = 'Scaffold tree'
RETURN DISTINCT d2.title;
```

## Files & implementation/test gaps

### 21. Where is a Feature implemented?

```cypher
MATCH (f:Feature {id:'feature:Chem:scaffold-tree'})-[r:IS_IMPLEMENTED_IN]->(file:File)
RETURN file.relative_path, file.language, file.line_count, r.member_count
ORDER BY r.member_count DESC, file.line_count DESC;
```

### 22. Where is a Feature tested?

```cypher
MATCH (f:Feature {id:'feature:Chem:scaffold-tree'})-[:IS_TESTED_IN]->(file:File)
RETURN file.relative_path, file.line_count;
```

### 23. **Test gaps** — Features with implementation but no tests

```cypher
MATCH (p:Package)-[:HAS_FEATURE]->(f:Feature)
WHERE EXISTS  { MATCH (f)-[:IS_IMPLEMENTED_IN]->(:File) }
  AND NOT EXISTS { MATCH (f)-[:IS_TESTED_IN]->(:File) }
RETURN p.name, f.name ORDER BY p.name LIMIT 30;
```

### 24. Most-tested Features

```cypher
MATCH (f:Feature)-[:IS_TESTED_IN]->(file:File)
RETURN f.id, count(file) AS test_files
ORDER BY test_files DESC LIMIT 15;
```

### 25. **Doc gaps** — Features with implementation but no docs

```cypher
MATCH (p:Package)-[:HAS_FEATURE]->(f:Feature)
WHERE EXISTS  { MATCH (f)-[:IS_IMPLEMENTED_IN]->(:File) }
  AND NOT EXISTS { MATCH (f)<-[:DOCUMENTS]-(:DocPage) }
RETURN p.name, f.name ORDER BY p.name LIMIT 30;
```

### 26. **Orphan files** — files with no DefinedIn entity (helper/utility surface)

```cypher
MATCH (file:File)
WHERE NOT EXISTS { MATCH ()-[:DEFINED_IN]->(file) }
  AND file.language IN ['ts', 'js']
  AND file.is_test = false
  AND file.is_generated = false
RETURN file.package_id, file.relative_path, file.line_count
ORDER BY file.line_count DESC LIMIT 20;
```

### 27. Files with the most line count per package (heaviest implementation)

```cypher
MATCH (p:Package)-[:CONTAINS_FILE]->(file:File)
WHERE file.is_test = false AND file.is_generated = false
WITH p, file ORDER BY file.line_count DESC
WITH p, collect(file)[0] AS biggest
RETURN p.name, biggest.relative_path, biggest.line_count
ORDER BY biggest.line_count DESC LIMIT 15;
```

### 28. Test-to-impl ratio per Feature

```cypher
MATCH (f:Feature)-[:IS_IMPLEMENTED_IN]->(impl:File)
OPTIONAL MATCH (f)-[:IS_TESTED_IN]->(test:File)
WITH f, count(distinct impl) AS impl_count, count(distinct test) AS test_count
WHERE impl_count > 0
RETURN f.id, impl_count, test_count
ORDER BY impl_count DESC LIMIT 15;
```

## Semantic types

### 29. Most-consumed semtype across the platform

```cypher
MATCH (s:SemanticType)<-[:CONSUMES_SEMTYPE]-(f:RegisteredFunction)
RETURN s.name, count(distinct f) AS consumers
ORDER BY consumers DESC LIMIT 10;
```

### 30. Functions that detect a given semtype

```cypher
MATCH (f:RegisteredFunction)-[:DETECTS_SEMTYPE]->(s:SemanticType {name:'Molecule'})
RETURN f.package_id, f.name;
```

### 31. **Cross-package semtype flow** — who produces what others consume

```cypher
MATCH (producer:Package)-[:EXPORTS]->(pf:RegisteredFunction)-[:PRODUCES_SEMTYPE]->(s:SemanticType)
      <-[:CONSUMES_SEMTYPE]-(cf:RegisteredFunction)<-[:EXPORTS]-(consumer:Package)
WHERE producer.name <> consumer.name
RETURN producer.name, s.name, consumer.name, count(*) AS edges
ORDER BY edges DESC LIMIT 15;
```

## Cross-package code edges

### 32. Top runtime callers (Package → another Package's function)

```cypher
MATCH (a:Package)-[:CALLS]->(f:RegisteredFunction)<-[:EXPORTS]-(b:Package)
RETURN a.name AS caller, b.name AS target, count(distinct f) AS funcs
ORDER BY funcs DESC LIMIT 15;
```

### 33. Functions called by ≥ 2 other packages (popular API surface)

```cypher
MATCH (caller:Package)-[:CALLS]->(f:RegisteredFunction)<-[:EXPORTS]-(owner:Package)
WHERE caller.name <> owner.name
WITH f, owner, count(distinct caller) AS callers
WHERE callers >= 2
RETURN owner.name, f.name, callers ORDER BY callers DESC LIMIT 15;
```

### 34. Two packages that share the most libraries (similarity proxy)

```cypher
MATCH (a:Package)-[:DEPENDS_ON]->(l:Library)<-[:DEPENDS_ON]-(b:Package)
WHERE a.name < b.name
WITH a, b, count(distinct l) AS shared
RETURN a.name, b.name, shared ORDER BY shared DESC LIMIT 10;
```

## TS code structure (TsClass / TsMethod / TsFunction across layers)

These entity kinds exist in **all three layers** — js-api (`is_jsapi=true`), packages, libraries — distinguished by `source_layer` and `is_jsapi`.

### 35. Class count per package (top declarations)

```cypher
MATCH (p:Package)-[:HAS_CLASS]->(c:TsClass)
RETURN p.name, count(c) AS classes
ORDER BY classes DESC LIMIT 15;
```

### 36. Largest classes by method count (anywhere)

```cypher
MATCH (c:TsClass)
WHERE c.is_jsapi = false AND c.method_count > 0
RETURN c.namespace, c.name, c.method_count
ORDER BY c.method_count DESC LIMIT 15;
```

### 37. Top-level helper functions a Package exports (TsFunction, NOT registered)

```cypher
MATCH (p:Package {name:'Chem'})-[:HAS_FUNCTION]->(f:TsFunction)
RETURN f.name, f.return_type, f.is_async, f.parameter_count
ORDER BY f.name LIMIT 20;
```

### 38. Plugins that subclass a given JS API class

```cypher
MATCH (sub:TsClass)-[:EXTENDS_CLASS]->(parent:TsClass {name:'Widget', is_jsapi:true})
RETURN sub.namespace, sub.name;
```

## JS API & DapiEndpoint

### 39. Most-used JS API class across plugins

```cypher
MATCH (p:Package)-[u:USES_API_CLASS]->(c:TsClass {is_jsapi:true})
WITH c, sum(u.use_count) AS total_uses, count(distinct p) AS users
RETURN c.name, total_uses, users ORDER BY total_uses DESC LIMIT 15;
```

### 40. Which DapiEndpoint each Package calls

```cypher
MATCH (p:Package)-[c:CALLS_DAPI_ENDPOINT]->(e:DapiEndpoint)
RETURN p.name, e.accessor_path, c.call_count
ORDER BY c.call_count DESC LIMIT 20;
```

### 41. The Dart-side method behind a JS API method (DelegatesTo bridge)

```cypher
MATCH (m:TsMethod)-[:DELEGATES_TO]->(b:GeneratedBinding)
WHERE m.class_id = 'ts:class:DG.DataFrame'
RETURN m.name, b.name, b.dart_class, b.dart_member LIMIT 15;
```

### 42. Function returns DapiEndpoint typed for which entity

```cypher
MATCH (e:DapiEndpoint)-[:TYPED_FOR]->(c:TsClass)
RETURN e.accessor_path, c.name AS served_entity LIMIT 20;
```

### 43. JS API enum/constant → SemanticType bridge

```cypher
MATCH (k:TsConstant)-[:DEFINES_SEMTYPE]->(s:SemanticType)
RETURN k.name, s.name LIMIT 15;
```

### 44. **Dead deps** — packages that DEPENDS_ON another but never CALLS into it

```cypher
MATCH (a:Package)-[:DEPENDS_ON]->(b:Package)
WHERE NOT EXISTS { MATCH (a)-[:CALLS]->(:RegisteredFunction)<-[:EXPORTS]-(b) }
RETURN a.name AS depender, b.name AS unused_dep
ORDER BY a.name LIMIT 20;
```

(Caveat: also unused if the dep is consumed via another mechanism we don't yet model — Dart-side, npm transitive, etc. Treat as a hint.)

## Useful patterns

### Filter on provenance / confidence

```cypher
MATCH (f:Feature)<-[r:PART_OF_FEATURE]-(m)
WHERE r.confidence >= 0.9 AND r.role = 'core'
RETURN f.name, m.name, r.weight;
```

### Anti-join (find missing relationship)

```cypher
MATCH (p:Package)
WHERE NOT EXISTS { MATCH (p)-[:HAS_FEATURE]->() }
RETURN p.name;     -- packages that didn't get clustered
```

### Aggregate after group

```cypher
MATCH (p:Package)-[:HAS_FEATURE]->(f:Feature)
WITH p, count(f) AS feature_count
WHERE feature_count >= 5
RETURN p.name, feature_count
ORDER BY feature_count DESC;
```

### Multi-hop existence (Kuzu syntax)

```cypher
-- "Does there exist a path from A to B through Documents and PartOfFeature?"
MATCH (f:RegisteredFunction {name:'Activity Cliffs'})
MATCH (d:DocPage {title:'Activity Cliffs'})
WHERE EXISTS { MATCH (d)-[:DOCUMENTS]->()-[:PART_OF_FEATURE]->(:Feature)<-[:PART_OF_FEATURE]-(f) }
RETURN d.id, f.id;
```

## Kuzu Cypher gotchas

- **No pattern comprehensions**: `[(a)-[:R]->(b) | b]` is not supported. Rewrite with `OPTIONAL MATCH … WITH … collect(b)`.
- **No label wildcards**: `MATCH (a)-[r]->(b)` works, but `MATCH (a:*)` does not. Use `label(a)` to inspect.
- **`extras` and `paths` are STRING (JSON-encoded)**: use `CONTAINS` for substring; for typed access, parse JSON in your script.
- **Empty arrays are NULL**: when extracting, an empty `tags: []` becomes NULL in the column. Use `WHERE x IS NOT NULL`.
- **Reserved words**: `default`, `match`, `from`, `to`, `value` need backticks if used as identifiers.

## Performance

The Kuzu DB is ~150 MB and indexes id (primary key) per node table. Most single-hop queries return in &lt; 50 ms. Multi-hop traversals over `LINKS_TO` (5K edges) or `PART_OF_FEATURE` (2.8K) take 100-500 ms. Whole-graph aggregates (e.g. count by edge predicate) are sub-second.

## Multi-granularity JS API usage (since 2026-05)

`USES_API_CLASS`, `USES_API_ENUM`, `USES_UI_COMPONENT`, `CALLS_DAPI_ENDPOINT`, `SUBSCRIBES_TO_EVENT`, and `IMPORTS_NAMESPACE` are now emitted at **four** granularities for the same usage event: `Package`, `File`, `TsMethod`, `TsFunction`. Pick the one that fits the question — the cardinalities differ by an order of magnitude (Package edges aggregate, File/Method edges are per-site).

### 45. Method-level: which method opens a HELM editor dialog?

```cypher
MATCH (m:TsMethod)-[u:USES_UI_COMPONENT]->(c:UiComponent {name:'dialog'})
WHERE m.id CONTAINS 'Helm' OR m.id CONTAINS 'Sequence'
RETURN m.id, c.name, u.use_count;
```

### 46. File-level: top files using `ui.dialog`

```cypher
MATCH (f:File)-[u:USES_UI_COMPONENT]->(:UiComponent {name:'dialog'})
RETURN f.relative_path, u.use_count ORDER BY u.use_count DESC LIMIT 10;
```

### 47. Method-level: which methods subscribe to `grok.events.onTableAdded`?

```cypher
MATCH (m:TsMethod)-[s:SUBSCRIBES_TO_EVENT]->(:JsEventStream {name:'onTableAdded'})
RETURN m.id, s.subscribe_count ORDER BY s.subscribe_count DESC;
```

### 48. Feature → method → JS API surface (the chain)

> "Which Features end up calling `ui.dialog` somewhere in their implementation?"

```cypher
MATCH (feat:Feature)<-[:PART_OF_FEATURE]-(rf:RegisteredFunction)
      -[:DEFINED_BY_FUNCTION]->(m:TsMethod)
      -[u:USES_UI_COMPONENT]->(:UiComponent {name:'dialog'})
RETURN feat.name, sum(u.use_count) AS dialog_calls
ORDER BY dialog_calls DESC LIMIT 10;
```

### 49. Per-file dapi usage breakdown

```cypher
MATCH (f:File)-[c:CALLS_DAPI_ENDPOINT]->(e:DapiEndpoint)
WHERE f.relative_path STARTS WITH 'packages/Chem/'
RETURN f.relative_path, e.name, c.call_count
ORDER BY c.call_count DESC LIMIT 20;
```

## File-level imports (since 2026-05)

The `Imports` edge resolves every `import ... from '...'` in TS source to its actual target node — `File` (relative paths), `LibraryModule` (`@datagrok-libraries/X/...`), or `JsApiNamespace` (`datagrok-api/{dg,grok,ui}`). Aggregated to one edge per `(source_file, target)` with the union of imported symbols.

### 50. Per-file imports out of one file

```cypher
MATCH (:File {id:'file:packages/Bio/src/utils/cell-renderer.ts'})-[i:IMPORTS]->(t)
RETURN t.id, i.imported_symbols, i.is_type_only, i.is_namespace, i.import_count;
```

### 51. **Reverse import** — who imports this file?

```cypher
MATCH (other:File)-[i:IMPORTS]->(:File {id:'file:packages/Bio/src/utils/seq-helper/seq-helper.ts'})
RETURN other.relative_path, i.imported_symbols ORDER BY other.relative_path;
```

This is the question grep can't answer cleanly — you'd grep for the file's basename plus all its imports, then manually filter false positives.

### 52. Top files reaching into a library

```cypher
MATCH (f:File)-[i:IMPORTS]->(m:LibraryModule)
WHERE m.id STARTS WITH 'lib_mod:bio:'
RETURN f.relative_path, count(DISTINCT m) AS modules_imported,
       sum(i.import_count) AS total_imports
ORDER BY modules_imported DESC LIMIT 10;
```

### 53. Files that import a specific symbol

```cypher
MATCH (f:File)-[i:IMPORTS]->(t)
WHERE 'getHelmHelper' IN i.imported_symbols
RETURN f.relative_path, t.id LIMIT 20;
```

### 54. Files that namespace-import the JS API

```cypher
MATCH (f:File)-[i:IMPORTS]->(ns:JsApiNamespace)
RETURN f.relative_path, ns.name, i.import_count
ORDER BY i.import_count DESC LIMIT 10;
```

### 55. Type-only imports (compile-time deps that webpack drops)

```cypher
MATCH (f:File)-[i:IMPORTS {is_type_only:true}]->(t)
RETURN f.relative_path, t.id, i.imported_symbols LIMIT 20;
```

## Function tags (since 2026-05)

Tags are first-class nodes. `HasTag` is synthesized from BOTH the source `tags[]` array AND `meta.role` (when `role` matches a known platform-dispatched function role like `cellEditor`, `panel`, `viewer`, `init`, `autostart`, `semTypeDetector`, `fileViewer`, etc.). This bridges the inconsistency where Helm declares `tags: ['cellEditor']` while Chem declares only `meta.role: 'cellEditor'`.

### 56. All functions with a given tag (replaces scanning `tags[]` arrays)

```cypher
MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'cellEditor'})
RETURN rf.id, rf.package_id, rf.role;
```

### 57. Most-used tags

```cypher
MATCH (t:Tag)
RETURN t.name, t.function_count ORDER BY t.function_count DESC LIMIT 15;
```

## Donor verification with `RequiresColumnTag` (since 2026-05)

`meta.columnTags: 'quality=Macromolecule, units=helm'` is parsed into one `RequiresColumnTag` edge per `key=value` pair, with the value side as a `SemanticType` node. When key is `quality`, it also emits `ConsumesSemtype` (since `quality=X` is the platform's way of saying "this function only fires for columns whose semType is X").

### 58. **Donor verification** — does X function require a specific column shape?

```cypher
MATCH (rf:RegisteredFunction {id:'func:Helm:editMoleculeCell'})-[r:REQUIRES_COLUMN_TAG]->(s:SemanticType)
RETURN s.name AS required_value, r.tag_key AS required_key;
```

Returns `Macromolecule` (key=`quality`) and `helm` (key=`units`). **If your caller's data shape doesn't satisfy these, do NOT delegate to this function** — it will throw.

### 59. All cellEditors with their column scope (the donor-lookup query)

```cypher
MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'cellEditor'})
OPTIONAL MATCH (rf)-[r:REQUIRES_COLUMN_TAG]->(s:SemanticType)
RETURN rf.id, rf.package_id,
       collect({key: r.tag_key, value: s.name}) AS column_scope;
```

## PackageTest — per-test-block coverage (since 2026-05)

Every `test('name', async () => {...})` block under `packages/<Pkg>/src/tests/**/*.ts` is a `PackageTest` node with category, file_path, line, and `is_skipped`/`is_only` flags. Best-effort `Covers(PackageTest → RegisteredFunction)` edges link tests to functions when names overlap.

### 60. Tests in a category

```cypher
MATCH (t:PackageTest)
WHERE t.package_id = 'pkg:Chem' AND t.category CONTAINS 'similarity'
RETURN t.name, t.file_path LIMIT 30;
```

### 61. Skipped tests across all packages

```cypher
MATCH (t:PackageTest {is_skipped:true})
RETURN t.package_id, t.category, t.name LIMIT 30;
```

### 62. Tests that cover a function (name-match)

```cypher
MATCH (t:PackageTest)-[:COVERS]->(:RegisteredFunction {name:'toAtomicLevel'})
RETURN t.category, t.name, t.file_path;
```

### 63. RegisteredFunctions with no PackageTest covering them

```cypher
MATCH (p:Package {name:'Bio'})-[:EXPORTS]->(rf:RegisteredFunction)
WHERE NOT EXISTS { MATCH (rf)<-[:COVERS]-(:PackageTest) }
  AND rf.role IN ['app','viewer','panel','widget','transform','cellEditor','cellRenderer']
RETURN rf.id, rf.role LIMIT 20;
```

## Where the graph isn't enough yet (next round)

These remain grep territory until the next round of extractor work:

- **Function-body call graph at the symbol level** — `RegisteredFunction → RegisteredFunction` direct calls. Today only cross-package `Calls` (via `grok.functions.call`) is modeled.
- **Runtime preconditions in code** — `if (col.semType !== X) throw`, `seqHelper.getSeqHandler(col)`. Partially covered by `REQUIRES_COLUMN_TAG` for declared `meta.columnTags`, but body-level guards aren't extracted.
- **CSS class definitions** — `.d4-scatter-plot { ... }` selectors aren't parsed; use grep.
- **Recent code** — anything added since the last `build.py` run.

## When the graph isn't enough

The graph deliberately omits transient and runtime state. Fall back to:

| Question | Use instead |
|---|---|
| What did someone just commit? | `git log` |
| What's the live status of GROK-12345? | JIRA UI / `mcp__atlassian` |
| What's running on dev right now? | `grok s` CLI |
| Which CSS class styles `.d4-scatter-plot`? | grep / IDE find-in-files |
| How do I implement X in Dart? | `core/client/CLAUDE.md`, then read source |
