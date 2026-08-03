# Test Track / Queries — coverage report

**Scope:** all 14 manual scenarios under `Queries/` are now covered by some combination of an automated Playwright test and a manual `-ui.md` companion for the parts that automation can't reliably reach.

## Automated tests (11 files, all green on public)

| # | Scenario | Test file |
|---|---|---|
| 1–4 | adding / edit / browser / deleting | `postgres-query-lifecycle.test.ts` |
| 5 | ms-sql | `mssql-query-lifecycle.test.ts` |
| 6 | transformations | `postgres-query-transformations.test.ts` |
| 7 | browse-and-save-project | `chembl-parameterized-and-project.test.ts` |
| 8 | columns-inspect | `columns-inspect.test.ts` |
| 9 | get-all-get-top-100 | `get-all-get-top-100.test.ts` |
| 10 | new-sql-query | `new-sql-query.test.ts` |
| 11 | new-visual-query | `visual-query-and-params.test.ts` |
| 12 | query-layout | `query-layout.test.ts` |
| 13 | query-postprocessing | `query-postprocessing.test.ts` |
| 14 | visual-query-advanced | `visual-query-advanced.test.ts` |

The whole suite drives the real UI wherever it can — clicking menu items, typing into the SQL editor, right-clicking tree nodes, etc. — and only falls back to the JavaScript API for setup, cleanup, and the few cases where the UI is technically unreachable from the browser-automation layer.

## Manual companions (4 files)

Some manual steps depend on platform behaviour that browser automation cannot drive — drag-and-drop into canvas-rendered viewers, the visual-query builder's column pickers, Layout-tab viewer placement, and a couple of timing-related platform quirks specific to in-session edits. Each affected scenario has a `-ui.md` file that lists exactly the steps a human still needs to walk through:

- `new-visual-query-ui.md` — visual-builder fields and Debug-tab debug
- `query-layout-ui.md` — Layout-tab viewers and project save/reopen
- `query-postprocessing-ui.md` — Post-Process editor typing, Layout-tab viewers, balloon-on-Run end-to-end
- `visual-query-advanced-ui.md` — full visual-builder + share + layout + project save/reopen + Refresh-with-Enrich

Each `-ui.md` clearly states what the autotest already covers and what the human is being asked to verify, so manual passes stay short and focused.

## What's verified vs. what stays manual

- **Verified automatically every run:** query creation, SQL editing, running, saving, editing, deletion, browsing, parameterised-query parameter dialogs, toolbar parameter refresh, transformations, schema/column inspection, project save (with in-memory layout check), and the post-process body round-tripping through the server.
- **Stays manual (covered by `-ui.md`):** visual-query builder configuration, Layout-tab viewer setup, the green info-balloon on Run for freshly-edited queries, sharing dialogs, opening saved projects from the Browse tree, and Refresh-with-Enrich.

Net effect: every manual scenario is now either automated end-to-end or has a precisely scoped manual checklist for the parts that still need a human eye.
