# Diff Studio — EMS (Entity Management System) Migration Plan

## Overview

The goal is to move Diff Studio model storage from files (`.ivp` in file storage) to an **EMS
domain schema**: the plugin declares PostgreSQL tables mapped to platform entities and gets,
out of the box, managed CRUD with permission checks, row/column-level security, an audit trail,
filtering, history, sharing, and a typed JS API — with no hand-written SQL or backend code.

The additional (and headline new) value is **run history** (`run`) as first-class entities:
today runs are ephemeral; after migration they can be stored, searched, shared, and discussed.

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
| Runs/results | **nowhere** | ephemeral: solver → DataFrame |

Key files: `src/app.ts` (~2300 lines, UI and all path/file handling), `src/hub.ts`
(Templates/Library/Recent gallery), `src/utils.ts` (file/recent caching and reads),
`src/ui-constants.ts` (`PATH`, `MISC`), `src/scripting-tools.ts` (`.ivp` parser).

## Target data model

Schema `diffstudio`, file `databases/diffstudio/schema.json`. Three tables (a fourth is optional).

### Table `model` — user models (private, shared per row)

- `securityMode: "row"` + `defaultRowVisibility: "none"` — reproduces MyFiles semantics:
  a row is visible only to its author until they share it.
- **No `businessKey`**: model names may collide across users, whereas a business key / `unique`
  applies globally across live rows. Deduplication is not needed for user models.

### Table `library_model` — built-in library (public, read-only)

- `securityMode: "table"`, View granted to everyone on deploy; only the package author edits.
- `businessKey: ["name"]` — names are package-controlled and unique.
- Populated from today's `files/library/*.ivp` (and `templates/*.ivp`) — see Stage 4.

### Table `run` — run history (greenfield)

- `securityMode: "master"`, `delegate: "model_id"` — a run inherits its model's security
  (Edit on the model = Edit on its runs).
- `businessKey: ["model_id", "started_on"]`.

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
        "category":    {"type": "string", "choices": ["PK", "PK-PD", "chemistry", "bioreactor", "other"]},
        "source":      {"type": "string", "required": true},
        "method":      {"type": "string", "choices": ["ros34prw", "lsoda", "cvode", "mrt", "ros3prw", "rk3", "rk4", "rkdp", "ab4", "ab5"], "default": "ros34prw"}
      }
    },
    "run": {
      "securityMode": "master",
      "delegate": "model_id",
      "businessKey": ["model_id", "started_on"],
      "friendlyName": "Run",
      "columns": {
        "model_id":    {"type": "ref", "ref": "model", "required": true, "onDelete": "cascade"},
        "method":      {"type": "string"},
        "params":      {"type": "string"},
        "started_on":  {"type": "datetime", "required": true},
        "duration_ms": {"type": "int", "min": 0},
        "status":      {"type": "string", "choices": ["ok", "error", "timeout"], "default": "ok"},
        "result":      {"type": "file"}
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
| `run` | `master → model` | inherited from the model | inherited from the model |

Sharing a model is standard (View/Edit/Delete/Share on the row); on first share the row is
promoted to a full entity (sharing dialog, favorites, comments, global search).

## Migration stages

Each stage stands on its own and does not break the previous one. User-facing flows switch
**behind a feature flag** (`diffStudioEms`); the file layer stays as a fallback until the final
stage.

### Stage 0. Setup and prototype
- Create `databases/diffstudio/schema.json` (tables `model` + `run`; `library_model` and
  `scenario` come later).
- Run `grok api` → generates `src/generated/db.ts` with the typed client `diffstudioDb`.
- Deploy to dev (behind the flag); verify the generated interfaces and generic access
  `grok.dapi.domains.table('diffstudio.model')`.
- **Exit criterion:** insert/query/get/update/delete on `diffstudio.model` work from the console.

### Stage 1. Save/load a model via the client (behind the flag)
- `saveToMyFiles()` (`app.ts:1228`): when the flag is on → `diffstudioDb.models.insert/save(...)`
  instead of `files.writeAsText`. Name/description/category/method come from the parsed `.ivp`
  (`getIVP`); `source` = the full model text.
- Opening a model: wherever a `.ivp` is currently read by path, add a branch that reads the row
  by `id` (`models.get(id)`) and parses `source` with the same `getIVP`.
- Deep-link / `PATH`/URL params: introduce entity-`id` addressing (`?model=<uuid>`) alongside
  the file-based one; old links keep working via the fallback.
- **Exit criterion:** behind the flag a model saves and opens from the table; with the flag off,
  behavior is unchanged.

### Stage 2. Run history `run` (greenfield, fastest payoff)
- In the solver path (after `solve*` in `solver-tools.ts` / the launch point in `app.ts`), insert
  a `run` row: `model_id`, `method`, `params` (input JSON), `started_on`, `duration_ms`,
  `status`; `result` (CSV/DataFrame) — optional into the `file` column.
- A "Runs" pane/tab on the model:
  `diffstudioDb.runs.query().where({model_id}).orderBy('started_on', false)`.
- **Exit criterion:** runs are saved, listed under the model, filterable, and inherit the model's
  permissions.

### Stage 3. UI: Recent and gallery backed by the DB
- `hub.ts` `buildRecent()`: source becomes a query over `model` (e.g. sorted by `updated_on`)
  instead of `diff-studio-recent.d42`. Or keep Recent as-is initially.
- Consider replacing gallery sections with `@datagrok-libraries/domain-ui` components
  (`domains.table('diffstudio.model')` → `form/grid/list/app`) where it simplifies the code.
- **Exit criterion:** gallery and Recent read models from the table; UX no worse than today.

### Stage 4. Built-in Library as `library_model`
- Populate `library_model` from `files/library/*.ivp` and `files/templates/*.ivp` (seed on
  deploy or via a one-off script); carry over icons/descriptions.
- `hub.ts` `buildLibrary()`/`buildExternalModelCards()` and `external-models.json` → read from
  `library_model`. Decide the fate of `external-models.json` (can be deprecated).
- **Exit criterion:** Templates/Library render from the table, public and read-only.

### Stage 5. Migrating existing users' data
- One-off import of `.ivp` from `<Login>:Home/*.ivp` into `model` (by author), or lazy dual-read:
  if a model is not found in the DB, read the file and offer to import it.
- Keep deep-link backward compatibility (old file links resolve and import on open).
- **Exit criterion:** no user "loses" models; old links open.

### Stage 6. Removing the flag and cleanup
- Once EMS leaves Beta and data is migrated, make EMS the default path; keep the file layer only
  for upload/export.
- Remove the dead `.d42` recent and file-based model-save code.

## Code map (where changes land)

- `databases/diffstudio/schema.json` — **new**, source of truth.
- `src/generated/db.ts` — generated by `grok api`, never hand-edited.
- `src/app.ts`: `saveToMyFiles()` (1228), `saveModelToRecent()` (2654), model opening/preview,
  deep-link and `PATH`/URL logic.
- `src/hub.ts`: `buildRecent()`, `buildLibrary()`, `buildExternalModelCards()`,
  `openCustomModelSettings()`.
- `src/utils.ts`: recent/file caches (`getCachedRecentModelsTable`, `getCachedFileInfo`, etc.) —
  adapt to DB queries.
- `src/ui-constants.ts`: `PATH`, `MISC` — new `id`-addressing constants.
- `src/scripting-tools.ts`: `getIVP()` — reused as the `source` parser (unchanged).

## Known EMS limitations and workarounds

- **Beta / feature flag** — hence the whole plan is additive and flagged.
- **No entity-typed column:** a `run` cannot yet natively reference its `FuncCall`, `Space`, or
  `Func`. v1 workaround: store `params`/metadata as a JSON string; defer the FuncCall link until
  an entity column exists (echoes #pharm-sphere: "funccall as an entity first, then v2 on domain
  dbs").
- **Platform→EMS references:** sticky-meta/entity properties cannot reference an EMS table
  (EMS cannot be used as a controlled vocabulary) — no impact on Diff Studio in v1.
- **No write-once/immutable columns** — if a "frozen" snapshot of run inputs is needed, rely on
  the audit trail or store the snapshot in `params`.
- **Destructive schema changes** require an explicit `migrations` section; a **package-managed
  schema** can be purged only by an admin.

## Open questions (decide before coding)

1. `source` as `string` vs `file`. Recommendation: `string` (compact `.ivp`, easy to
   search/parse). — **proposed: string**
2. Uniqueness of user model names: keep non-unique (like files) or add a `key`. —
   **proposed: non-unique, no businessKey**
3. Run `params`: JSON string vs property schema. For v1 — JSON string. — **proposed: JSON**
4. Whether to store the run `result` (can be large): off by default, opt-in. —
   **proposed: optional**
5. Fate of `external-models.json` and `.d42` recent after Stages 3–4 (deprecate?).
6. Migration strategy for legacy `.ivp`: one-off import vs lazy dual-read.

## Testing

- Unit/integration: CRUD on `diffstudio.model`/`run` (insert/query/get/update/delete, versions,
  soft-delete), permissions (author sees, another user does not; after sharing — sees), `run`
  permission inheritance.
- Run the existing Playwright suites (`playwright/`, `diff-studio-hub-test-plan.md`) under both
  feature-flag states.
- Check the Grit reference as the model for expected UI/security behavior.

## Definition of Done (pilot)

- `schema.json` deploys, `grok api` generates the client, CRUD works.
- Behind the flag: a model saves/opens from the table; privacy as in MyFiles; sharing works.
- Run history is written and shown under the model, inheriting permissions.
- Existing tests green in both flag states.
- Legacy `.ivp` migration defined and verified on test data.
