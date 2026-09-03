# Diff Studio — EMS (Entity Management System) Migration Plan

## Overview

The goal is to move Diff Studio model storage from files (`.ivp` in file storage) to an **EMS
domain schema**: the plugin declares PostgreSQL tables mapped to platform entities and gets,
out of the box, managed CRUD with permission checks, row/column-level security, an audit trail,
filtering, history, sharing, and a typed JS API — with no hand-written SQL or backend code.

Run history (`run`) as first-class entities was considered the headline addition, but is
**deferred to a future version** (see [Deferred: run history](#deferred-run-history)): it
structurally depends on models already living in the DB and cannot cover built-in library models
or unsaved drafts without extra machinery. **v1 is scoped to migrating model storage** from files
to the domain schema.

EMS status at time of writing: **Beta**, behind a feature flag, v1 merged into master
2026-08-13, targeted for the 1.27 release. For that reason the migration is planned to be
**additive and behind a feature flag**, not a hard cutover.

EMS documentation:
- dev: `public/help/develop/how-to/db/domain-schemas.md`
- user: `public/help/govern/catalog/domains.md`
- UI library: `public/libraries/domain-ui/README.md`
- reference packages: **Grit** (full app on domain-ui), **Inventory** (batch/upsert/aggregate),
  **PlatesFixture** (security modes)

## Current state (what we migrate)

| What | Where today | Mechanism in code |
|------|-------------|-------------------|
| Model = `.ivp` text | — | declarative format |
| User models | `<Login>:Home/<name>.ivp` | `app.ts` → `saveToMyFiles()` (`grok.dapi.files.writeAsText`) |
| Custom Library | `System:AppData/DiffStudio/library` + `external-models.json` | `hub.ts` → `buildExternalModelCards()`, `openCustomModelSettings()` |
| Built-in Templates/Library | `files/` inside the package | static resources |
| Recent (≤10) | `diff-studio-recent.d42` (binary DataFrame, columns `Info`, `IsCustom`) | `app.ts` → `saveModelToRecent()`; rendered in `hub.ts` → `buildRecent()` |
| Runs/results (**out of v1 scope**) | **nowhere** | ephemeral: solver → DataFrame |

Key files: `src/app.ts` (~2300 lines, UI and all path/file handling), `src/hub.ts`
(Templates/Library/Recent gallery), `src/utils.ts` (file/recent caching and reads),
`src/ui-constants.ts` (`PATH`, `MISC`), `src/scripting-tools.ts` (`.ivp` parser).

## Target data model

Schema `diffstudio`, file `databases/diffstudio/schema.json`. **Two tables in v1** (`model`,
`library_model`); `scenario` is optional/later; `run` is deferred (see
[Deferred: run history](#deferred-run-history)).

### Table `model` — user models (private, shared per row)

- `securityMode: "row"` + `defaultRowVisibility: "none"` — reproduces MyFiles semantics:
  a row is visible only to its author until they share it.
- **No `businessKey`**: model names may collide across users, whereas a business key / `unique`
  applies globally across live rows. Deduplication is not needed for user models. Overwrite
  semantics are reproduced by tracking the row id in the editor (see Stage 1), not by a natural key.
- **No `category` column**: the `.ivp` format has no model-level category (the `category:` seen in
  `.ivp` is a per-input UI annotation), and no current model declares `#tags`. A classification
  dimension for user models is left to a future version once there is a UI to set it.

### Table `library_model` — built-in library (public, read-only)

- `securityMode: "table"`, View granted to everyone on deploy; only the package author edits.
- `businessKey: ["name"]` — names are package-controlled and unique.
- `category` is an optional freeform string, filled by the package author on seeding.
- Populated from today's `files/library/*.ivp` (and `templates/*.ivp`) — see Stage 3.

### Table `scenario` (optional, later stage) — named input sets

- `securityMode: "master"`, `delegate: "model_id"` — reusable/shareable parameter sets for
  sensitivity/fitting.

### Draft `schema.json`

```json
{
  "name": "diffstudio",
  "version": "1.0.0",
  "tables": {
    "library_model": {
      "securityMode": "table",
      "businessKey": ["name"],
      "friendlyName": "Library model",
      "columns": {
        "name":        {"type": "string", "required": true, "isName": true},
        "description": {"type": "string"},
        "category":    {"type": "string"},
        "source":      {"type": "string", "required": true},
        "icon":        {"type": "string"}
      }
    },
    "model": {
      "securityMode": "row",
      "defaultRowVisibility": "none",
      "friendlyName": "Model",
      "columns": {
        "name":        {"type": "string", "required": true, "isName": true},
        "description": {"type": "string"},
        "source":      {"type": "string", "required": true},
        "method":      {"type": "string", "choices": ["ros34prw", "lsoda", "cvode", "mrt", "ros3prw", "rk3", "rk4", "rkdp", "ab4", "ab5"], "default": "ros34prw"}
      }
    }
  }
}
```

Reminder: system columns (`id`, `version`, `created_on`, `updated_on`, `author_id`,
`is_deleted`) are added automatically — do not declare them.

## Security (summary)

| Table | Mode | Who sees | Who edits |
|-------|------|----------|-----------|
| `library_model` | `table` | everyone (View on deploy) | package author |
| `model` | `row` + `defaultRowVisibility: none` | author; plus those it is shared with | author / anyone with Edit |

Sharing a model is standard (View/Edit/Delete/Share on the row); on first share the row is
promoted to a full entity (sharing dialog, favorites, comments, global search).

## Feature flag

The EMS path is gated by a **package property** `diffStudioEms` declared in `package.json`
(`propertyType: bool`, `defaultValue: false`). It is read synchronously through a helper
`isEmsEnabled()` in `utils.ts`:

```ts
export function isEmsEnabled(): boolean {
  return _package.settings?.['diffStudioEms'] === true;
}
```

An admin flips it per deployment (**Manage > Plugins > DiffStudio**) with no rebuild; tests toggle
it with `_package.setSettings({diffStudioEms: true/false}, group)`. Compare explicitly with
`=== true` (values may arrive as strings). The **database schema always deploys on publish** —
the flag gates only the code paths that choose DB vs files, never the creation of the tables
(empty tables are harmless).

## Migration stages

Each stage stands on its own and does not break the previous one. User-facing flows switch
**behind the feature flag** (`diffStudioEms`); the file layer stays as a fallback until the final
stage.

### Stage 0. Setup and prototype
- Create `databases/diffstudio/schema.json` with table `model` (`library_model` is added in
  Stage 3; the draft above shows the eventual v1 shape).
- Add the `diffStudioEms` package property to `package.json` and the `isEmsEnabled()` helper.
- Run `grok api` → generates `src/generated/db.ts` with the typed client `diffstudioDb`.
- Deploy to dev; verify the generated interfaces and generic access
  `grok.dapi.domains.table('diffstudio.model')`.
- **Exit criterion:** insert/query/get/update/delete on `diffstudio.model` work from the console.

### Stage 1. Save/load a model via the client (behind the flag)
- `saveToMyFiles()` (`app.ts:1228`): when the flag is on → the DB path instead of
  `files.writeAsText`. The editor tracks the current model's row id (`currentModelId`):
  - first save → `diffstudioDb.models.insert(...)`, remember the returned `id`;
  - subsequent saves of the same model → `diffstudioDb.models.update(currentModelId, ..., {version})`;
  - "Save As" → clear `currentModelId`, then `insert`.

  This reproduces file-overwrite semantics without a `businessKey` (names are not unique across
  users). `currentModelId` is a transient per-model field — **drain it in `resetForReuse()`**
  (see CLAUDE.md, or it leaks between consecutive previews).
- Field mapping: `name`/`description` from the parsed `.ivp` (`getIVP`); `method` from
  `#meta.solver`; `source` = the full model text. (No `category` — not present in the format.)
- Opening a model: wherever a `.ivp` is currently read by path, add a branch that reads the row
  by `id` (`models.get(id)`) and parses `source` with the same `getIVP`.
- Deep-link / `PATH`/URL params: introduce entity-`id` addressing (`?model=<uuid>`) alongside
  the file-based one; old links keep working via the fallback.
- **Exit criterion:** behind the flag a model saves (insert, then update-in-place), opens from the
  table, and re-saving updates rather than duplicates; with the flag off, behavior is unchanged.

### Stage 2. UI: Recent and gallery backed by the DB
- `hub.ts` `buildRecent()`: source becomes a query over `model` (e.g. sorted by `updated_on`)
  instead of `diff-studio-recent.d42`. Or keep Recent as-is initially.
- Consider replacing gallery sections with `@datagrok-libraries/domain-ui` components
  (`domains.table('diffstudio.model')` → `form/grid/list/app`) where it simplifies the code.
- **Exit criterion:** gallery and Recent read models from the table; UX no worse than today.

### Stage 3. Built-in Library as `library_model`
- Add `library_model` to the schema (additive change). Populate from `files/library/*.ivp` and
  `files/templates/*.ivp` (seed on deploy or via a one-off script); carry over descriptions and
  icons. `icon` = the relative path under `files/` (e.g. `icons/pk.png`) — no `file`-typed column
  needed.
- `hub.ts` `buildLibrary()`/`buildExternalModelCards()` and `external-models.json` → read from
  `library_model`. Decide the fate of `external-models.json` (can be deprecated).
- **Exit criterion:** Templates/Library render from the table, public and read-only.

### Stage 4. Migrating existing users' data
- One-off import of `.ivp` from `<Login>:Home/*.ivp` into `model` (by author), or lazy dual-read:
  if a model is not found in the DB, read the file and offer to import it.
- Keep deep-link backward compatibility (old file links resolve and import on open).
- **Exit criterion:** no user "loses" models; old links open.

### Stage 5. Removing the flag and cleanup
- Once EMS leaves Beta and data is migrated, make EMS the default path; keep the file layer only
  for upload/export.
- Remove the dead `.d42` recent and file-based model-save code.

## Deferred: run history

Storing solver runs (`run`) as first-class, searchable, shareable entities was the original
headline. It is **out of v1 scope** for three reasons:

- A `run` must reference a persisted `model` (`model_id`), so it can only cover **saved user
  models** — not built-in library models (a separate table) or unsaved drafts, which is the more
  common "just try it" flow.
- EMS has **no entity-typed column** yet: a `run` cannot natively reference its `FuncCall`,
  `Space`, or `Func` — inputs/metadata would have to be a JSON string until an entity column
  exists (echoes #pharm-sphere: "funccall as an entity first, then v2 on domain dbs").
- Run results can be large; the storage policy (opt-in `file` column vs none) needs its own design.

When revisited, the shape is a `run` table in `securityMode: "master"` delegating to `model_id`
(Edit on the model = Edit on its runs), plus a decision on library/draft coverage
(auto-save-to-`model` on first run, or a nullable second reference). Note also that an append-only
run log likely wants **no `businessKey`** (the system `id` already makes rows unique; a
`[model_id, started_on]` key risks dropping same-millisecond runs as duplicates).

## Code map (where changes land)

- `databases/diffstudio/schema.json` — **new**, source of truth.
- `src/generated/db.ts` — generated by `grok api`, never hand-edited.
- `package.json` — **new** `properties` entry for the `diffStudioEms` flag.
- `src/utils.ts`: `isEmsEnabled()` helper (**new**); recent/file caches
  (`getCachedRecentModelsTable`, `getCachedFileInfo`, etc.) — adapt to DB queries.
- `src/app.ts`: `saveToMyFiles()` (1228) + a `currentModelId` field (drained in
  `resetForReuse()`), `saveModelToRecent()` (2655), model opening/preview, deep-link and
  `PATH`/URL logic.
- `src/hub.ts`: `buildRecent()`, `buildLibrary()`, `buildExternalModelCards()`,
  `openCustomModelSettings()`.
- `src/ui-constants.ts`: `PATH`, `MISC` — new `id`-addressing constants.
- `src/scripting-tools.ts`: `getIVP()` — reused as the `source` parser (unchanged).

## Known EMS limitations and workarounds

- **Beta / feature flag** — hence the whole plan is additive and flagged.
- **Platform→EMS references:** sticky-meta/entity properties cannot reference an EMS table
  (EMS cannot be used as a controlled vocabulary) — no impact on Diff Studio in v1.
- **Destructive schema changes** require an explicit `migrations` section; a **package-managed
  schema** can be purged only by an admin.

## Open questions and resolved decisions

Resolved during review:
- **`model.category`** — dropped (format has no model-level category; no model declares `#tags`).
  `library_model.category` stays as an optional freeform string filled by the package author.
- **User model name uniqueness** — non-unique, no `businessKey`; the editor tracks `currentModelId`
  and updates in place (see Stage 1).
- **Feature flag mechanism** — package property `diffStudioEms`, read via `_package.settings`.
- **Run history** — deferred (see above).

Still open:
1. `source` as `string` vs `file`. Recommendation: `string` (compact `.ivp`, easy to
   search/parse). — **proposed: string**
2. Fate of `external-models.json` and `.d42` recent after Stages 2–3 (deprecate?).
3. Migration strategy for legacy `.ivp`: one-off import vs lazy dual-read.

## Testing

- Unit/integration: CRUD on `diffstudio.model` (insert/query/get/update/delete, versions,
  soft-delete), the insert-then-update-in-place save flow (re-saving does not duplicate),
  permissions (author sees, another user does not; after sharing — sees).
- Run the existing Playwright suites (`playwright/`, `diff-studio-hub-test-plan.md`) under both
  feature-flag states (toggle via `_package.setSettings({diffStudioEms: ...})`).
- Check the Grit reference as the model for expected UI/security behavior.

## Definition of Done (pilot)

- `schema.json` deploys, `grok api` generates the client, CRUD works.
- Behind the flag: a model saves/opens from the table; re-saving updates in place (no duplicates);
  privacy as in MyFiles; sharing works.
- With the flag off, behavior is unchanged.
- Existing tests green in both flag states.
- Legacy `.ivp` migration defined and verified on test data.
