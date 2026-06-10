---
name: datagrok-knowledge-graph
description: Use for structural questions about the Datagrok codebase — "what implements X", "where is class/function Y defined", "which tests/docs cover Z", "who imports this file", "what does package A depend on", coverage/gap analysis. Query the knowledge graph FIRST, before grepping the source. Do NOT use for runtime/live state (use the Datagrok MCP), recent uncommitted code (use git/grep), or Dart core.
---

# datagrok-knowledge-graph

A queryable Kuzu graph of the whole public codebase lives at `workspace/.kg`. It
answers structural questions far faster than grep returns candidate matches.

## How to query

```
workspace/.kg/.venv/bin/python workspace/.kg/qq.py "<cypher>"
```

First call warms a persistent local server (~1-3s); later calls return in tens of ms.
Stored paths are repo-root-relative — a result `packages/Bio/src/foo.ts` is on disk at
`workspace/packages/Bio/src/foo.ts`. Always use `workspace/.kg/...` relative paths; never
the absolute `/workspace` form (blocked at the hook level).

Full schema + 60-query cookbook: `workspace/.kg/docs/QUERYING.md`. Agent guide:
`workspace/.kg/CLAUDE.md`.

## Node kinds (node tables)

`Package`, `RegisteredFunction`, `Script`, `DataQuery`, `DataConnection`, `Library`,
`LibraryModule`, `DockerContainer`, `Feature` (LLM-derived), `File`, `SemanticType`,
`TsClass`, `TsMethod`, `TsFunction`, `DocPage`, `HelpAnchor`, `Tutorial`, `ChangelogEntry`,
`JiraTicket`, `Tag`, `PackageTest`, `DapiEndpoint`, `UiComponent`, `JsEventStream`.

Every node has `id` (PK), `name`, `kind`, `description`, plus kind-specific fields.

## Key edges

`EXPORTS` (Package→RegisteredFunction), `DEPENDS_ON` (Package/Library→…),
`IMPORTS` (File→File/LibraryModule/JsApiNamespace), `DEFINED_IN` (entity→File),
`HAS_FEATURE` / `PART_OF_FEATURE`, `IS_IMPLEMENTED_IN` / `IS_TESTED_IN` (Feature→File),
`DOCUMENTS` (DocPage→entity), `HAS_TAG` (→Tag), `COVERS` (PackageTest→RegisteredFunction),
`USES_API_CLASS` / `USES_UI_COMPONENT` / `CALLS_DAPI_ENDPOINT` / `SUBSCRIBES_TO_EVENT`
(at Package / File / TsMethod / TsFunction granularity), `EXTENDS_CLASS`,
`CONSUMES_SEMTYPE` / `PRODUCES_SEMTYPE` / `DETECTS_SEMTYPE`.

## Example queries

What implements a feature:
```cypher
MATCH (f:Feature)-[:IS_IMPLEMENTED_IN]->(file:File)
WHERE toLower(f.name) CONTAINS 'activity cliff'
RETURN f.name, file.relative_path;
```

Where a class is defined:
```cypher
MATCH (c:TsClass {name:'MacromoleculeSequenceCellRenderer'})-[d:DEFINED_IN]->(f:File)
RETURN f.relative_path, d.line_start;
```

Reverse import — who imports this file:
```cypher
MATCH (other:File)-[i:IMPORTS]->(:File {id:'file:packages/Bio/src/utils/seq-helper/seq-helper.ts'})
RETURN other.relative_path, i.imported_symbols;
```

Files importing a specific symbol:
```cypher
MATCH (f:File)-[i:IMPORTS]->(t)
WHERE 'getHelmHelper' IN i.imported_symbols
RETURN f.relative_path;
```

All functions with a given role/tag:
```cypher
MATCH (rf:RegisteredFunction)-[:HAS_TAG]->(:Tag {name:'cellEditor'})
RETURN rf.id, rf.package_id, rf.role;
```

Test gaps — features implemented but not tested:
```cypher
MATCH (p:Package {name:'Bio'})-[:HAS_FEATURE]->(f:Feature)
WHERE EXISTS { MATCH (f)-[:IS_IMPLEMENTED_IN]->(:File) }
  AND NOT EXISTS { MATCH (f)-[:IS_TESTED_IN]->(:File) }
RETURN f.name;
```

Which packages subscribe to an event:
```cypher
MATCH (p:Package)-[s:SUBSCRIBES_TO_EVENT]->(:JsEventStream {name:'onTableAdded'})
RETURN p.name, s.subscribe_count;
```

Docs for a function:
```cypher
MATCH (rf:RegisteredFunction {name:'Recalculate Coordinates'})<-[d:DOCUMENTS]-(doc:DocPage)
RETURN doc.title, d.derivation;
```

## Kuzu gotchas

- No pattern comprehensions; rewrite with `OPTIONAL MATCH … collect()`.
- `extras` / `paths` are JSON strings — use `CONTAINS` for substring matching.
- Reserved words (`default`, `match`, `from`, `to`, `value`) need backticks.

## When the graph isn't enough

Fall back to grep/git for **recent uncommitted code**, to the Datagrok MCP for **live
server state**, and to `workspace/core/` source for **Dart**. The graph is a snapshot —
verify a specific path before acting on it.
