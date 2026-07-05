# Grit changelog

## v.next

* GROK-20308: Introduced Grit (GRok Issue Tracker) — reference app for entity-mapped domain schemas: projects/issues/comments/labels manifest, minimal UI over `grok.dapi.domains`, issue CRUD tests
* GROK-20319: Moved the app and tests to the `grok api`-generated typed clients (`src/generated/db.ts`); the issue detail view now loads comments in the same query via `expand: ['details:comment']` instead of a second per-issue query

## 0.0.1 (2026-07-04)
