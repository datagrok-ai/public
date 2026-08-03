# Sticky Meta — manual test docs

Manual test scenarios for the **Sticky Meta** feature (Browse > Platform > Sticky Meta),
plus the **Database meta** context-panel surface for database tables and columns.

Sticky Meta lets you attach custom metadata ("sticky meta") to any object in Datagrok.
Metadata is stored in Postgres and shared instantly across all contexts. The feature is
in **Beta**. Wiki: `public/help/govern/catalog/sticky-meta.md`.

These docs are written first as **manual** scenarios. After review they are rewritten as
Playwright autotests under this folder (`*.test.ts`), exercising the UI as much as possible
and falling back to the JS API only where the UI cannot do the step (setup, verification
reads, cleanup).

## Target & accounts

| Item | Value |
|------|-------|
| Server | `DATAGROK_URL` from `.env` (e.g. `https://dev.datagrok.ai`) |
| Account | `DATAGROK_LOGIN` / `DATAGROK_PASSWORD` — must be an **Administrator** (creating entity types / schemas requires admin) |
| Auth | Pre-authenticated via `globalSetup` → `e2e/.auth.json` |

## Scenario index

| # | File | Area | UI-only? |
|---|------|------|----------|
| 1 | `01-schema-and-type.md` | Create / edit / delete entity type + metadata schema | Yes |
| 2 | `02-add-and-edit.md` | Add metadata to a cell, sticky columns, batch edit | Yes |
| 3 | `03-persistence-copy-delete.md` | Metadata survives clone / project / refresh / relogin; delete values | Mostly (project/space/export are JS-API) |
| 4 | `04-database-meta.md` | Database meta on a DB table and column | Yes |

## Fixtures and naming

- All test-created schemas and entity types use the prefix **`PW_SM_`** plus a unique suffix
  (e.g. timestamp), so parallel runs and leftover state never collide.
- Scenarios that need a molecule-matching schema create their own (matching expression
  `semtype=molecule`) — they do **not** depend on a pre-seeded `TestSchema1`.
- Datasets are opened from `System:DemoFiles/...` (e.g. `SPGI.csv`, which has a `Structure`
  column detected as `Molecule`).

## Cleanup rule (mandatory)

**Every schema, entity type, or metadata value a test creates must be removed at the end of
that test** (`try` / `finally`). Leftover schemas pollute the Sticky meta context-panel pane
for every molecule column on the shared dev server.

Cleanup paths:

- **UI**: in the Schemas / Types view, search for the entity by name → right-click the card →
  **Delete** → confirm **DELETE**.
- **JS API** (used in autotest `finally` blocks): a schema is deleted **by id**, and the JS
  `Schema` wrapper does not expose `id`, so read it from the Dart handle:

  ```js
  const schema = (await grok.dapi.stickyMeta.getSchemas()).find((s) => s.name === name);
  const id = window.grok_Entity_Get_Id(schema.dart); // Schema wrapper has no .id getter
  await grok.dapi.stickyMeta.deleteSchema(id);        // deleteSchema(id: string)
  ```

## Key findings from recon (2026-06, dev)

- **Schema property types**: `string`, `int`, `bool`, `double`, `datetime`, `string_list`.
- **Database meta** context-panel section exists on a **DB table** and a **DB column**, reached
  via `Browse > Databases > Postgres > <connection> > Schemas > <schema> > <table> > <column>`.
  - Table fields: Domains, Row Count, Comment, LLM Comment + **Save**.
  - Column fields: Is Unique, Min, Max, Values, Sample Values, Unique Count, Quality, Comment,
    LLM Comment + **Save**.
- **Database meta is NOT shown on the connection (DbInfo) entity** — only on tables and columns.
  Connection-level metadata is therefore out of scope for `04-database-meta.md`.
- Known platform issues to watch on column save: **GROK-19427** (column meta save error
  `Value test doesn't match expected type string_list`) and **GROK-19429** (Row count deleting
  not working).

## Run (after autotests are written)

```bash
cd reddata/playwright-tests
npx playwright test e2e/stickymeta
```
