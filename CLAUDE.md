# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is the **Datagrok public repository** - the open-source monorepo for the Datagrok data analytics and visualization platform. It contains the JavaScript/TypeScript API, shared libraries, CLI tools, and 76+ extension packages.

## Knowledge graph — query it before you grep

**Read this before you reach for grep/find/rg.** This repo ships a **queryable knowledge graph** of every package, registered function, script, query, library, TS class/method, doc page, tutorial, changelog entry, and Jira ticket — plus the edges between them. It lives at [`.kg/`](.kg/) as a ~2.4 MB committed binary.

For any **structural** question, your FIRST action must be a `.kg/scripts/qq.py` query — **not** a filesystem search. Use it whenever you need to:

- Find which plugin owns a feature: *"what implements substructure search?"*
- Find documentation for a function: *"where is `searchSubstructure` documented?"*
- Find who imports / tests / extends / calls / depends on something.
- Find coverage gaps, or trace a Jira ticket to commits, packages, and docs.

```bash
python3 .kg/scripts/qq.py "MATCH (p:Package {name:'Chem'})-[:HAS_FEATURE]->(f:Feature) RETURN f.name LIMIT 10"
```

Run it with plain system Python — **no setup required**. The first call self-installs its engine into `.kg/.venv` (~30s) and unpacks the DB; every call after returns in ~5ms. Cookbook and node/edge schema are in [`.kg/CLAUDE.md`](.kg/CLAUDE.md).

Reach for grep **only** when the question isn't structural (e.g. a raw string literal) or the graph comes up empty. Don't grep `packages/` or `libraries/` to answer a "who implements / owns / uses X" question — that's what the graph is for.

## Repository Structure

```
public/
├── js-api/          # Core JavaScript/TypeScript API (datagrok-api)
├── libraries/       # Shared TypeScript libraries (@datagrok-libraries/*)
├── packages/        # Extension packages (viewers, connectors, scientific tools)
├── tools/           # CLI tool "grok" (datagrok-tools)
├── python-api/      # Python API bindings
├── help/            # Documentation (Docusaurus)
├── connectors/      # GrokConnect Java/Maven JDBC connectors
├── docker/          # Docker deployment configs
├── environments/    # Environment configurations
└── misc/            # ESLint config, utilities
```

## js-api (datagrok-api)

The core TypeScript API providing bindings to the Datagrok Dart backend. Packages import three namespaces:

- **`grok`** - High-level APIs: shell, functions, events, data operations, server API, AI
- **`ui`** - UI components: elements, dialogs, inputs, menus, viewers
- **`DG`** (from `datagrok-api/dg`) - Re-exports all types, constants, and classes

Key modules in `js-api/src/`: `dataframe.ts` (DataFrame, Column, BitSet), `entities.ts` (User, Project, Package), `widgets.ts` (Dialog, Menu, InputBase), `viewer.ts`, `grid.ts`, `dapi.ts` (server HTTP API), `shell.ts`, `functions.ts`, `events.ts`, `const.ts` (enums).

Files ending in `.g.ts` or `.api.g.ts` are auto-generated from the Dart codebase — do not edit manually.

## Libraries (@datagrok-libraries/*)

Reusable TypeScript libraries published under `@datagrok-libraries` scope in `libraries/`:

| Library | Purpose |
|---------|---------|
| **bio** | Macromolecule & Molecule3D data support |
| **chem-meta** | RDKit JS API, molfile parsing |
| **compute-api** | Compute global API |
| **compute-utils** | Compute-related utilities |
| **cruddy** | CRUD app framework |
| **db-explorer** | Database exploration |
| **math** | High-performance computation |
| **ml** | Machine learning utilities |
| **statistics** | Statistics & aggregation |
| **utils** | Common utilities (hashing, encoding) |

## Tools (datagrok-tools / `grok` CLI)

The `grok` CLI manages the full package lifecycle: `grok create`, `grok publish`, `grok test`, `grok link`, `grok check`, `grok api`. Install via `npm install -g datagrok-tools` or link locally from `tools/`.

Server config is stored in `~/.grok/config.yaml`.

### Managing a running Datagrok server — use `grok s`

For any task that creates, updates, shares, or inspects entities on a running server
(users, groups, connections, queries, scripts, packages, files, function calls, raw
API), use `grok s` — **not** Playwright / Selenium / hand-written HTTP. It hits the
public REST API directly, is idempotent, and scripts cleanly from the shell.

Quick orientation:

```bash
grok s users list                                    # list entities
grok s groups get Admins --output json               # introspect the JSON shape
grok s users save --json user.json                   # create / update from JSON
grok s groups add-members Chemists alice bob         # membership (idempotent)
grok s shares add "Me:MyConn" Chemists --access Edit # sharing
grok s functions run 'Chem:smilesToMw("CCO")'        # call a function
grok s raw GET /api/users/current                    # any endpoint
```

Full reference with JSON shapes, batch operations, and scripting patterns:
[`tools/GROK_S.md`](tools/GROK_S.md).

## Code Style

- 2-space indentation
- Single quotes for strings
- Semicolons required
- Windows line endings (CRLF)
- TypeScript strict mode
- File naming: kebab-case (`my-component.ts`)
- Package naming: letters, numbers, underscores, hyphens
- Semantic versioning: `X.Y.Z` or `X.Y.Z-rc` or `X.Y.Z-rc.N`
- Never edit `.g.ts` or `.api.g.ts` files (auto-generated by `grok api`)

## Before Implementing

Before creating or adding any package artifact, check `.claude/skills/` for a matching skill and follow it.

When building UI components (viewers, file viewers, dialogs, layouts), use the `/ui` skill for authoritative guidelines.
