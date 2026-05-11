## Knowledge graph — query before you grep

A queryable graph of every plugin, function, script, doc page, tutorial,
JS API class, changelog entry, Jira ticket, **TS file, class, method, and
JS API usage at file/method granularity** lives at `.kg/`. ~33K nodes,
~95K edges. **Use it as your first lookup** for structural questions:

| Question | One-shot Cypher |
|---|---|
| "What implements activity cliffs?" | `MATCH (f:Feature)-[:IS_IMPLEMENTED_IN]->(file:File) WHERE toLower(f.name) CONTAINS 'activity cliff' RETURN file.relative_path` |
| "Where is `MacromoleculeSequenceCellRenderer` defined?" | `MATCH (c:TsClass {name:'MacromoleculeSequenceCellRenderer'})-[d:DEFINED_IN]->(f:File) RETURN f.relative_path, d.line_start` |
| "All cellEditors in the codebase, with their column scope" | `MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'cellEditor'}) OPTIONAL MATCH (rf)-[:REQUIRES_COLUMN_TAG]->(s:SemanticType) RETURN rf.id, collect(s.name)` |
| "What does Helm:editMoleculeCell require of its caller's column?" | `MATCH (rf:RegisteredFunction {id:'func:Helm:editMoleculeCell'})-[r:REQUIRES_COLUMN_TAG]->(s:SemanticType) RETURN r.tag_key, s.name` |
| "Which method opens a `ui.dialog`?" | `MATCH (m:TsMethod)-[u:USES_UI_COMPONENT]->(:UiComponent {name:'dialog'}) RETURN m.id, u.use_count` |
| "Who imports this file?" (reverse import lookup) | `MATCH (other:File)-[i:IMPORTS]->(:File {id:'file:packages/Bio/src/utils/seq-helper/seq-helper.ts'}) RETURN other.relative_path, i.imported_symbols` |
| "Which files import `getHelmHelper`?" | `MATCH (f:File)-[i:IMPORTS]->(t) WHERE 'getHelmHelper' IN i.imported_symbols RETURN f.relative_path` |
| "All subclasses of `DG.GridCellRenderer`" | `MATCH (sub:TsClass)-[:EXTENDS_CLASS]->(:TsClass {name:'GridCellRenderer', is_jsapi:true}) RETURN sub.id` |
| "Which packages subscribe to `grok.events.onTableAdded`?" | `MATCH (p:Package)-[s:SUBSCRIBES_TO_EVENT]->(:JsEventStream {name:'onTableAdded'}) RETURN p.name, s.subscribe_count` |
| "Test gaps in Bio" | `MATCH (p:Package {name:'Bio'})-[:HAS_FEATURE]->(f:Feature) WHERE NOT EXISTS { MATCH (f)-[:IS_TESTED_IN]->() } RETURN f.name` |
| "Per-test-block listing for Chem similarity tests" | `MATCH (t:PackageTest) WHERE t.package_id='pkg:Chem' AND t.category CONTAINS 'similarity' RETURN t.name, t.file_path` |

The graph almost always returns *the answer set* faster than ripgrep
returns *candidate matches*. Stored paths are repo-root-relative — a
graph result `packages/Bio/src/foo.ts` is on disk at the same path
relative to the public repo root.

**Multi-granularity** since 2026-05: every JS-API-usage edge
(`USES_API_CLASS`, `USES_UI_COMPONENT`, `CALLS_DAPI_ENDPOINT`,
`SUBSCRIBES_TO_EVENT`, `USES_API_ENUM`, `IMPORTS_NAMESPACE`) is emitted
at four levels — `Package`, `File`, `TsMethod`, `TsFunction`. Pick the
level that fits your question. Methods → traverse via
`(rf)-[:DEFINED_BY_FUNCTION]->(m:TsMethod)` to bridge to RegisteredFunction
or to `Feature` via `(rf)-[:PART_OF_FEATURE]->(f:Feature)`.

### How to query

```powershell
.kg/.venv/Scripts/python.exe .kg/qq.py "MATCH (p:Package {name:'Chem'})-[:HAS_FEATURE]->(f:Feature) RETURN f.name LIMIT 10"
```

First call auto-spawns a persistent server (~0.8s warmup); every
subsequent call is 5–50 ms. Cookbook: [`docs/QUERYING.md`](../../.kg/docs/QUERYING.md).

If the venv is missing run [`.kg/scripts/bootstrap.ps1`](../../.kg/scripts/bootstrap.ps1)
or `bash .kg/scripts/bootstrap.sh` once.

### When the graph is wrong or insufficient — fall back to grep

The graph is a **snapshot**. Trust it for structural facts; verify any
specific path/symbol you're about to act on with a quick `Read` or
`Grep`. Fall back to filesystem search when:

- The KG returns 0 rows for a query you have strong reason to believe
  should match (the snapshot may be stale, or this slice may not be
  modeled).
- You're asking about **recent code** — code added/renamed since the
  last `build.py` won't be in the graph. Use `git log` / `Grep`.
- You're asking about **Dart `core/`**, runtime server state, or live
  Jira fields — not in the graph. Use the source / `grok s` / Jira UI.
- The query errors out (e.g. `qq.py` can't reach the server). Restart
  the server (`qq.py --restart`) or fall back; don't block.

### When you discover something the graph doesn't know — write it down

If during your work you find a real fact the graph is missing — a
feature implementation in a file the KG didn't link, a cross-package
runtime call, an undocumented test that exercises a function — append
a one-line note to `.kg/.learned/<YYYY-MM-DD>-<short-slug>.md`:

```
feature:Bio:atomic-level-structure  IS_IMPLEMENTED_IN  packages/Bio/src/utils/helm-to-molfile/converter/connection-list.ts  reason: handles loop closure, found while debugging GROK-12345
func:Bio:To Atomic Level            IS_TESTED_IN       packages/Bio/src/tests/cyclic-helm-tests.ts                            reason: covers branched HELM
```

Format: one tab-separated line per fact: `subject_id  PREDICATE  object_id  reason: ...`. The Stop hook drains these into the next
enrichment cycle so subsequent agents see the new facts.

### Never

- Don't `grep -r` across `packages/` to "find a feature" before
  trying `MATCH (f:Feature {...})-[:IS_IMPLEMENTED_IN]->(file)`. The
  graph already has the answer.
- Don't write to `.kg/kg.kuzu/` directly. It's a derived cache —
  edit JSONL, run `build.py`.
