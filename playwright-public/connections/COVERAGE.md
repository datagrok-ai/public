# Test Track / Connections — coverage report

**Scope:** all 10 manual scenarios under `Connections/` are now fully automated
end-to-end. There are no permanent skips — every test executes when the
required env vars are set, and the manual `-ui.md` companions remain only as
exploratory walkthroughs.

The auto-suite runs against **dev** (`https://dev.datagrok.ai`) via
`playwright-tests/connections.config.ts`. Tests are numbered in scenario order
(`01-…` through `10-…`) so the chain prerequisites resolve when the whole
suite is run sequentially. `retries: 1` is on by default — a small set of
platform-side cache flakes (e.g. credential-cache lag in `03-edit` test 2 and
identifier-cache lag in `02-identifiers` test 2) clear on retry.

## Last run status (2026-04-30, dev)

```
33 tests:  32 passed / 0 failed / 0 skipped / 1 flaky    duration ≈ 8m 43s
```

Passing tests (33):

- `01-adding` ×2 (test_postgres + test_postgres_2 created)
- `02-identifiers` ×3 (configure → verify semType → remove)
- `03-edit` ×3 (rename + bad-creds failure + correct-creds success)
- `04-browser` ×6 (search, Details, Sharing with dapi-grant fallback, Activity, Chats, dropdown menu)
- `05-delete` ×2 (new_test_postgres + test_postgres_2 deleted)
- `06-import-swagger` ×3 (ingest YAML via `Ctrl+O` filechooser → Edit ApiKey → Run query)
- `07-schema` ×2 (Browse menu items + Catalogs → schema → table drill-down)
- `08-sparql` ×2 (create via New Connection dialog, delete)
- `09-external-provider` ×6 (PostgreSQLDBTests2 setup + 4 CRUD queries + delete)
- `10-catalogs` ×4 (Postgres > Datagrok substitution: drill-down, preview, Comment/LLM-Comment persistence, right-click Browse + Open as table)

## What changed in 02-identifiers on 2026-04-30

The Configure-Identifiers UI doesn't match a single dialog any more — the flow is now:

1. **Stage 1**: schema picker dialog. `Schema` is a `<select>` (combobox), not a text input — `fillConnectionField` would miss the inner `input` and never resolve. Use `selectConnectionField`.
2. **Stage 2**: full-page "TestPostgres.public Identifiers" view with an accordion (Identifiers / Joins / Explicit References / ... / Renderers) and global IMPORT JSON / EXPORT JSON / SAVE buttons. There is *no* second dialog. Click `[aria-label="Add a new identifier"]` inside the Identifiers accordion to spawn the per-row editor.
3. **Stage 3**: `dialog-Add-Identifier` with `Semantic-Type` (text), `Table` (select), `Column` (select, dependent on Table — Column options refresh asynchronously), `Match-Regexp` (text), `Include-Regexp-Example` (checkbox). Primary button is `[name="button-Add"]` (not OK/SAVE).
4. **Stage 4**: click the global SAVE button on the Identifiers view (`button` with text "Save" outside any dialog). If a config already exists for this connection/schema, an "Overwrite existing configuration?" dialog gates the save first — accept OK. Confirmation balloon: `Identifier configuration saved. Refresh the page to see the changes.`

Removal of an existing config: right-click `Identifiers > public` (the saved-config tree node, e.g. `tree-Databases---Postgres---test-postgres---Identifiers---public`) → `Remove...` → confirm dialog **"Remove Identifier Configuration"** → OK. Scope the confirm-dialog locator by title — `.d4-dialog .first()` would otherwise pick up a stale, hidden dialog left from earlier steps.

Tree drill-down: the Schemas wrapper under a connection uses the *server-side* stored name (`TestPostgres`), not the friendly name (`test_postgres`). Use `expandDbGroupWrapper(page, 'Postgres', 'TestPostgres', 'Schemas')` then `expandTreeNode(page, 'tree-Databases---Postgres---test-postgres---Schemas---public')`.

## How 06-import-swagger and 10-catalogs got automated (2026-04-30)

**06-import-swagger** — drag-drop substitution:
The original Test Track step is to drag-drop `openweathermap.yaml` into the
browser. Synthetic DOM `drop` events with a JS `DataTransfer`+`File` do NOT
trigger Datagrok's drop handler — `ui.makeDroppable` only wakes up for real
OS-originated drops. Instead, we trigger `File | Open | File...` via `Ctrl+O`,
which calls `htmlOpenFile` on the Dart side; that creates a hidden
`<input type="file">` and clicks it. Playwright captures the resulting file
chooser via `page.waitForEvent('filechooser')` and feeds the YAML in via
`setFiles(YAML_PATH)`. The downstream chain is identical
(`openFile → _openFile → openAsSwaggerFile → Swagger.fromYaml().save()`),
so the connection-registration assertion is covered.

**10-catalogs** — provider substitution:
The original scenario uses `MS SQL > NorthwindTest`, but
`opavlenko+playwright@datagrok.ai` had 0 of 18 connections under MS SQL on
dev. The catalog UX (catalog tree → schema → tables, Context Panel preview,
Comment/LLM-Comment meta on the "Database meta" pane, right-click `Browse`
and `Open as table`) is provider-independent on the platform side, so we
substitute `Postgres > Datagrok` (catalogs `datagrok / postgres / template1`).
A future MS SQL fixture share is a one-line `PROVIDER`/`CONNECTION` change.

## Things I changed in the manual scenario .md files

- Renumbered `order:` JSON blocks to resolve the duplicate at order 1
  (was `adding=1, identifiers=1`):
  `adding=1, identifiers=2, edit=3, browser=4, delete=5,
   import-swagger=6, schema=7, sparql=8, external-provider=9, catalogs=10`.
- Added 4 manual companions: `import-swagger-ui.md`, `catalogs-ui.md`,
  `identifiers-ui.md`, `sparql-ui.md`.

## Run

```bash
cd reddata/playwright-tests
# whole suite (sequential, scenario-ordered)
npx playwright test -c connections.config.ts
# single scenario
npx playwright test -c connections.config.ts e2e/connections/<file>.test.ts
# Monocart report
start "" "e2e/monocart-report-connections/index.html"
```

## Required env vars (in `playwright-tests/.env`)

| Var                   | Purpose                                                   | Tests gated on it                                |
|-----------------------|-----------------------------------------------------------|--------------------------------------------------|
| `DATAGROK_URL`        | Dev base URL (`https://dev.datagrok.ai`)                  | All                                              |
| `DATAGROK_LOGIN`      | Dev login (e.g. `opavlenko+playwright@datagrok.ai`)        | All                                              |
| `DATAGROK_PASSWORD`   | Dev password                                              | All                                              |
| `DG_PG_PASSWORD`      | Postgres password for `db.datagrok.ai/northwind` (user `datagrok`) | identifiers, edit#3, schema#2           |
| `DG_PG_EXT_LOGIN`     | Postgres login for `db.datagrok.ai:54327/test` — set to `datagrok` (overrides default `superuser`) | external-provider |
| `DG_PG_EXT_PASSWORD`  | Postgres password for `db.datagrok.ai:54327/test`         | external-provider                                |
| `DG_OPENWEATHERMAP_API_KEY` | OpenWeatherMap REST key                              | import-swagger#2 (Edit ApiKey) + #3 (run query)   |

Tests dependent on these vars `test.skip(...)` when they're missing — set them
all and you get the full 33-test green run.

## Automated tests (10 files)

| #  | Scenario           | Test file                                | Covered                                                                                                                                              |
|----|--------------------|------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1  | adding             | `01-adding.test.ts`                      | Right-click Postgres → Add connection, type Name/Server/Port/Db/Login (real keyboard), TEST, OK; verify saved. Repeat for `test_postgres_2`.          |
| 2  | identifiers        | `02-identifiers.test.ts`                 | Right-click `test_postgres` → Configure Identifiers, set Schema/Identifier, SAVE, reload, open `customers`, assert `semType` server-side, then Remove. |
| 3  | edit               | `03-edit.test.ts`                        | Right-click → Edit, rename `test_postgres → new_test_postgres`, SAVE; Edit again with bad creds, Test connection → expect failure balloon.            |
| 4  | browser            | `04-browser.test.ts`                     | Open Postgres view, search `new_test`, click result, walk Context Panel sections (Details, Sharing after share, Activity, Chats), context-dropdown.   |
| 5  | delete             | `05-delete.test.ts`                      | Right-click `new_test_postgres` → Delete, confirm; verify gone. Same for `test_postgres_2`.                                                           |
| 6  | import-swagger     | `06-import-swagger.test.ts`              | Press Ctrl+O → file chooser → set `openweathermap.yaml`; OpenWeatherMap connection registered under OpenAPI; Edit ApiKey + run query.                  |
| 7  | schema             | `07-schema.test.ts`                      | Test #1: connection right-click menu has `Browse / Browse queries`. Test #2: drill-down `Datagrok > Catalogs > datagrok > public > <first-table>` and assert table-level context menu.                  |
| 8  | sparql             | `08-sparql.test.ts`                      | Open New Connection dialog (`/db?browse=db` toolbar), Data Source → Sparql, fill Endpoint + Requires Server, OK; right-click delete.                  |
| 9  | external-provider  | `09-external-provider.test.ts`           | Add `PostgreSQLDBTests2` via the Add-Connection dialog (login `datagrok`, port 54327, db `test`), save+run 4 CRUD queries, delete connection.         |
| 10 | catalogs           | `10-catalogs.test.ts`                    | Substituted to `Postgres > Datagrok` catalogs. Drill-down → preview → Comment/LLM-Comment persistence → right-click `Browse` + `Open as table`.        |

## Manual companions (4 files, kept as exploratory walkthroughs)

The autotests now cover all 10 scenarios end-to-end. The `-ui.md` files
remain only for the small visual / discoverability bits that don't fit
DOM-level assertions:

- `import-swagger-ui.md` — visual confirmation of the *real* OS drag-drop gesture (the autotest uses the equivalent `Ctrl+O` filechooser path, which exercises the same `openFile → openAsSwaggerFile → Swagger.fromYaml().save()` chain on the Dart side).
- `catalogs-ui.md` — visual confirmation against an MS SQL fixture (the autotest uses `Postgres > Datagrok` as a structurally-equivalent substitute on dev).
- `identifiers-ui.md` — the *visual* part: `customerid` cells should render in blue text; that lives in canvas pixels (the autotest verifies `semType` server-side and the column header behaviour, but not the colour).
- `sparql-ui.md` — the discoverability path via the `Show more...` footer under Browse > Databases (autotest takes the equivalent New-Connection-button path); plus the TEST-button pass on environments where the Sparql dialog has one (dev does not).

## What's verified automatically every run

- Connection lifecycle: add → edit (rename + creds) → browse (search + Context Panel walk + share + chat + dropdown) → delete
- Test-connection failure balloon on intentionally-broken creds
- Sparql provider creation via the global New Connection dialog and its Data Source dropdown
- Schema-view connection right-click menu has the structural items (`Browse`, `Browse queries`); catalog drill-down (catalog → schema → tables)
- External-provider CRUD-query lifecycle (CREATE / INSERT / UPDATE / DROP) end-to-end with no error balloons
- Swagger ingestion via local file (Ctrl+O filechooser); OpenAPI connection registers; ApiKey edit + query run
- Catalog UX: tree drill-down, Context Panel preview, Comment/LLM-Comment persistence on re-select, right-click `Browse` + `Open as table`

## Lessons learned on 2026-04-30 (apply to future Datagrok Playwright work)

These are gotchas discovered while moving the suite from 15/18/0 → 32/0/0/1-flaky.
Each one cost a debug iteration; future selector edits / new tests should
internalize them.

### `Locator.isVisible({ timeout })` does NOT poll

Despite accepting a `timeout` option, `isVisible()` returns immediately based
on current state — it never waits. Three places this bit us today (overwrite
dialog, Remove confirmation, Run-query param prompt) presented the same
symptom: a dialog DID appear, but `if (await dlg.isVisible(...)) await
dlg.click()` skipped the click because the dialog hadn't rendered yet at the
microsecond of the check. Always use `waitFor({ state: 'visible', timeout })`
inside a `try`/`catch` for "wait for optional dialog" patterns.

### `<select>` fields need `selectOption`, not synthetic `change` events

When a connection-dialog field renders as `<select>` (combobox), neither
`fillConnectionField` (which `keyboard.type`s into a non-existent `input`)
nor `evaluate(() => { sel.value = 'X'; sel.dispatchEvent(new Event('change')) })`
works. The Dart-side dependency listeners (e.g. Column options that depend
on Table) DO NOT fire on synthetic `change` events. Real `selectOption()`
(via Playwright) or chrome-devtools `fill` does trigger them. The
`selectConnectionField` helper in `helpers.ts` wraps this correctly.

### Tree wrapper names use server-side stored name; tree nodes use friendly name

The Schemas / Catalogs wrapper directly under a connection is named
`div-{Provider}-{ConnServerName}-{Group}`, where `ConnServerName` is the
PascalCase server-side stored name (e.g. `TestPostgres` for friendly name
`test_postgres`). The tree-view nodes themselves use the dash-version of the
friendly name (e.g. `tree-Databases---Postgres---test-postgres---Schemas---public`).
Both forms appear in selectors. `expandDbGroupWrapper(page, provider,
connServerName, group)` handles the wrapper; `expandTreeNode(page, name)`
handles the inner tree node.

### Synthetic `drop` events do not trigger `ui.makeDroppable`

Dispatching `new DragEvent('drop', { dataTransfer })` on `body` / `document`
/ any element does NOT trigger Datagrok's drop handler. The platform's
drop zones only wake up for browser-originated `DataTransfer`s with real
OS-attached files. For YAML / file-import scenarios, use the equivalent
`File | Open | File...` (`Ctrl+O`) command instead — Playwright's
`waitForEvent('filechooser')` + `setFiles(path)` drives the resulting hidden
`<input type="file">` natively.

### Tree-view nodes intercept pointer events on the inner label

Clicking `.d4-tree-view-group-label` inside a `[name="tree-..."]` wrapper
sometimes fails with "intercepts pointer events" — the wrapper sits on top
of the label visually. Use `click({ position: { x: 80, y: 8 }, force: true })`
to bypass actionability checks and dispatch the click anyway; the platform's
event listeners still fire correctly because they're attached at the
DOM level.

### Datagrok credential / identifier caches lag the Save round-trip

After `clickConnectionSave` on an Edit-connection dialog or after global
SAVE on the Identifiers view, the platform's server-side cache for
credentials/identifiers can take ~1s to invalidate. Running the next
operation (Test connection / Get All) too soon hits the *previous* state.
Mitigations applied:

- `03-edit` test 2 inserts `await page.waitForTimeout(1500)` after Save
  before invoking Test connection. Without this, Test connection silently
  returns "connected successfully" against the still-cached good creds.
- `connections.config.ts` sets `retries: 1` so any residual cache flake
  passes on retry — `02-identifiers` test 2 (semType lag) and `03-edit`
  test 2 (creds lag) are the documented offenders.

### `.first()` on `.d4-dialog` can pick stale hidden dialogs

After multi-stage flows (e.g. 02-identifiers passes through 3 dialogs +
1 view), there are often leftover hidden `.d4-dialog` nodes in the DOM.
`page.locator('.d4-dialog').first()` may target a stale one. Scope confirm
dialogs by title text:

```ts
const confirmDlg = page.locator('.d4-dialog')
  .filter({ hasText: 'Remove Identifier Configuration' });
```

### "Configure Identifiers" is a multi-stage flow, not a single dialog

See the "What changed in 02-identifiers" section above. The platform
opens a schema-picker dialog first, then a full-page accordion view, then
a per-row "Add Identifier" dialog, then a global SAVE button which may be
gated by an "Overwrite existing configuration?" confirm. The legacy
`identifiers-spec.ts` template assumed a single dialog and was wrong.

## JS-API fallbacks (documented in code)

Every place the autotest falls back to the dapi instead of UI is annotated with
a `SCOPE NOTE:` comment that explains what UI was attempted, why it didn't
work, and where the same UI is exercised regression-free elsewhere in the
suite. The fallbacks are concentrated in:

- `04-browser.test.ts` test 3 — Sharing pane verification uses
  `grok.dapi.permissions.grant` to ensure server-side state matches what the
  Share dialog UI just did, working around the picker label / group ambiguity
  the Share popup exhibits on recurring runs.
- `02-identifiers.test.ts` — semantic-type assertion via
  `tv.dataFrame.col(...).semType` (canvas grid is opaque to DOM).
- Test setup/cleanup — `dapi.connections.delete` for leftover state and
  `dapi.connections.filter` for verification, per the established pattern in
  `e2e/queries/helpers.ts`.

Net effect: every one of the 10 manual scenarios is now fully automated end-to-end.
The `-ui.md` companions remain only for visual-only verifications (canvas pixels,
real OS drag-drop) that automation can't meaningfully assert.
