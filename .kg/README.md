# Datagrok Knowledge Graph (`.kg/`)

A queryable knowledge graph of the Datagrok codebase: code, docs, tests, infra, and process artifacts — connected so that "what does this bug touch?", "which features lack tests?", "show me everything related to scatter plot" become single Cypher queries.

Lives inside the `public/` git submodule and ships with the open-source repo. **Stored paths are anchored at `public/`** — a graph result like `packages/Bio/src/foo.ts` maps to `packages/Bio/src/foo.ts` on disk.

## Three-link entry points

- **For agents** (Claude Code, Copilot, etc.) → [CLAUDE.md](CLAUDE.md). Read this before grepping.
- **For users + developers** → [docs/QUERYING.md](docs/QUERYING.md). Cookbook of 20+ useful Cypher recipes.
- **For maintainers** → [docs/UPDATING.md](docs/UPDATING.md). Adding entities, extractors, enrichers; schema migrations.

## Quick start

```powershell
# One-time setup
py -m venv .venv
.\.venv\Scripts\Activate.ps1
py -m pip install -r requirements.txt

# Build the whole graph (~40s for the first run)
py build.py

# Query
py query.py --demos
py query.py "MATCH (p:Package {name:'Chem'})-[:HAS_FEATURE]->(f:Feature) RETURN f.name LIMIT 10"

# Visualize (open kg.html in a browser, no server needed)
py viz.py

# Tests (33 tests across schema/coverage/connectivity/crosslayer/features/gaps)
py tests/test_kg_queries.py
```

## Architecture (one-paragraph version)

**JSONL-canonical, Kuzu-derived.** The truth lives in `data/{nodes,edges}/*.jsonl` (committed to git, diffable, code-reviewable). The `kg.kuzu/` graph DB is an index, rebuilt from JSONL by `build.py` in seconds. Edit JSONL or extractors that produce it — never write to the DB directly. The schema is a set of Pydantic models in `schema/`; Kuzu DDL is derived from those. Adding a new entity kind = edit one file, rerun build.

```
.kg/
├── schema/        Pydantic models — single source of truth for shape
├── extractors/    Deterministic parsers (annotations, package.json, regex)
├── enrichers/     LLM-driven enrichment (Feature clustering, doc bridges)
├── data/          CANONICAL JSONL — committed
│   ├── nodes/     One file per entity kind: Package.jsonl, Feature.jsonl, …
│   └── edges/     One file per <Subject>__<PREDICATE>__<Object> triple
├── kg.kuzu/       MATERIALIZED graph DB — gitignored, rebuildable
├── build.py       extractors → JSONL → Kuzu COPY
├── enrich.py      prepare/apply contract for LLM enrichers
├── query.py       ad-hoc Cypher CLI
├── viz.py         standalone HTML viewer generator
├── kg.html        ← open in browser (regenerate with `py viz.py`)
├── CLAUDE.md      ← agent-facing entry point
├── docs/
│   ├── QUERYING.md   ← user/developer cookbook
│   └── UPDATING.md   ← maintainer guide
└── tests/
    ├── test_schema_smoke.py
    └── test_kg_queries.py
```

## What's in the graph (current state, 2026-05-10)

| Layer | Counts |
|---|---|
| Packages × Libraries × LibraryModules | 79 × 18 × 515 |
| RegisteredFunction / Script / DataQuery | 1043 / 545 / 508 |
| TsClass / TsMethod / TsFunction / TsInterface (across js-api, packages, libraries) | 1177 / 5153 / 2556 / 1399 |
| DocPage / HelpAnchor / Tutorial / TutorialTrack | 695 / 3089 / 22 / 7 |
| ChangelogEntry / JiraTicket | 2918 / 306 |
| **Feature** (LLM-derived) / **PartOfFeature** memberships | **458 / 2810** |
| **IS_IMPLEMENTED_IN / IS_TESTED_IN edges** | **3358 / 518** |
| **Total nodes / edges** | **27,433 / 42,947** |
| Build time (clean rebuild) | ~5 min (LLM enrichment cached in `data/enrichment/`) |

Top packages by feature richness: Chem (27) · Bio (21) · ApiSamples (22) · DBTests (16) · Samples (15) · UsageAnalysis (11) · PowerPack (12) · Chembl (12) · Charts (12) · EDA (12).

## Layer-agnostic by design

Every `Entity` has `source_layer` (`plugins | libraries | core_client | core_server | core_shared | docs | process | infra | synthetic`). The extractors built so far cover **plugins**, **libraries**, **docs**, **process** (changelog + Jira), and **synthetic** (Feature). Extending into `core_client` / `core_server` requires no schema changes — just write the analogous extractors.

## Provenance everywhere

Every node carries `description_provenance`; every edge carries `derived_by` and `confidence`. So you can always answer "how do we know this?" — whether the fact came from an annotation comment, a JSON file, an AST scan, an LLM cluster, or a manual edit.

## Tests

`tests/test_kg_queries.py` is 33 assertions in 6 categories — schema integrity, coverage minima, connectivity, cross-layer traversals, feature invariants, known coverage gaps. **Run them after every schema or extractor change.** Two real bugs were caught and fixed this way (CHANGELOG paragraph parsing, `package.js` fallback).

## What's not in the graph (yet)

- Dart client (`core/client/`) and Dart server (`core/server/`) — schema is layer-agnostic, just need extractors. Note: these live one level above `public/`, so any future Dart extractor needs to walk `../core/...` from `REPO_ROOT` (which now points at `public/`).
- Git commits and live Jira state.
- ApiTest / PlaywrightScenario nodes (extractors pending).
- Live Datagrok server state (use `grok s` for that).

## License

Internal Datagrok artifact, not for redistribution.
