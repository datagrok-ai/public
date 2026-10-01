# Forge architecture

How Forge works now. Forge is predictive modeling as a package: engines that live in other packages, models stored in
an EMS schema owned by Forge, and the UI on top.

## Key decisions

- **Own schema instead of the platform model entity.** Models, training runs and applications are rows of the EMS
  schema `forge`. That gives typed storage, row security and an audit trail without core changes.
- **Engines stay in their packages.** Forge calls them through the `mlname`/`mlrole` contract and copies no engine code.
- **Two tools coexist.** The built-in `ML | Models` tool keeps working until the switchover. Forge uses its own
  identifiers (`ML | Forge` menu, app, schema) and never touches the old tool's.
- **Training data never leaves the browser without consent.** Only the `copy` storage mode uploads it.
- **Application records inherit the model's security.** A record is visible to whoever can see its model. It can be
  written by anyone the `model` table lets insert, which today is every user through the table's Edit grant, so
  applications by users who can only view a model are recorded too.
- **Later:** MLflow models and the migration of existing models.

## Components and data flow

```
package.ts ──> ui/forge-app.ts ──> engines/           (DG.Func registry)
                               └─> generated/db.ts    (grok.dapi.domains -> EMS schema forge)
```

- `package.ts` holds the annotated functions only: the app `forgeApp` (shown as "Forge", `Browse > Apps`) and
  `forgeModels` (`ML | Forge | Models`). Both delegate to `ForgeApp`.
- `ForgeApp.create()` discovers the engines, loads the catalog frame and the column captions, then builds the view.
  `ForgeApp.open()` adds it to the shell and is the only error boundary: `ForgeError` becomes a warning, anything
  else becomes an error balloon plus `_package.logger.error`.
- Logic folders (`engines/` now; `preparation/`, `training/`, `metrics/`, `storage/`, `apply/` later) work without the
  DOM, throw instead of catching, and never import `ui/`. The same logic serves the UI, API functions and tests.

Next stages, not implemented yet:

- **Train:** prepare the table, call the engine's `train`, compute the metrics, record a training run.
- **Save:** write the model blob to file storage and the `model` row.
- **Apply:** match the columns, replay the preparation, call the engine's `apply`, record an application.

## EMS schema `forge`

Manifest: `databases/forge/schema.json`. `grok publish` deploys it; a debug publish applies destructive changes
without migration scripts, except a `promotion` change on a row table, which is always refused (`[promotion-change]`).
To change it, publish a manifest with only a placeholder table, then the real one, bumping `version` both times; the
first publish drops the tables and their data. `grok api` generates the typed client `src/generated/db.ts`
(`forgeDb.models`, `forgeDb.trainingRuns`, `forgeDb.applications`). Every table also has the system columns `id`,
`version`, `created_on`, `updated_on` and `author_id`; who and when always come from them.

| Table | Security | Why |
|---|---|---|
| `model` | `row`, `promotion: lazy`, `defaultRowVisibility: none`, grants `All users: view, edit` | A model is private until shared, like the old models. Lazy promotion makes a model a platform entity on its first share; sharing, favorites and comments work from then on. The table Edit grant lets any user insert models; with visibility `none` it reveals no foreign rows |
| `training_run` | `row`, lazy promotion, `defaultRowVisibility: none`, same grants | The trainer's experiment history, including failed and cancelled runs. `model_id` is optional with `onDelete: setnull`, so the history outlives the model |
| `application` | `master`, `delegate: model_id`, `audit: false`, no grants | Secured by the model (see Key decisions). Master tables refuse grants. No audit: high-churn records. `onDelete: cascade` from the model |

Deletes are soft: a deleted row stays with `is_deleted` and disappears from queries and counts. Current EMS limitation:
a soft delete does not clean up a promoted row's entity and permissions, so a model deleted after being shared leaves
a live entity until the platform fixes it.

### JSON column shapes

EMS `json` columns hold objects only; a list is wrapped in an object.

| Column | Shape |
|---|---|
| `model.target` | `{name, type, semType?, categories?}` |
| `model.features`, `training_run.features` | `{columns: [{name, type, semType?}]}`, in training order |
| `options` | Preparation replay, see below |
| `hyperparameters` | `{<train function input>: value}` |
| `metrics` | `{validation: {<metric id>: number}, train?: {<metric id>: number}}`; ids `mse`, `rmse`, `mae`, `r2`, `accuracy`, `f1`, `auc` |
| `splitting` | `{scheme: none \| kfold \| holdout, folds?, trainFraction?, isStratified?}` |
| `model.dataset_ref` | Reference mode: `{kind: file \| query \| script, id?, path?, params?, script?}` |
| `dataset_fingerprint` | `{rowCount, columnCount, hash, columns: [{name, type, min?, max?, mean?, categories?}]}` |

`options` keeps the keys of the platform's built-in models, so their preparation replays unchanged:
`preprocessingInfo`, `postprocessingInfo`, `positiveClass`, `negativeClass`, `binaryClassificationThreshold`,
`targetType`, `allowNulls`.

### Model blob

The `model.blob` column (type `file`) stores `file://System:DomainFiles/forge/model/<uuid>/model.bin`, the layout the
EMS row editor uses for file columns. The random path segment makes the file unguessable, and the row security of
`model` protects the pointer. The platform's file renderer and download work on the column as is.

### Storage modes

`model.storage_mode` (default `none`) records what a model keeps of its training data:

- `none`: only `dataset_fingerprint`, enough to check that a new table fits the model.
- `reference`: `dataset_ref` points to the source (file, query or script); the data is not uploaded.
- `copy`: `dataset_table_id` is the id of an uploaded copy of the training table (no foreign key; checked on read).

`has_training_rows` marks blobs that embed training rows (SVM, KNN). `legacy_id` is reserved for migrated models.

### Catalog captions

The catalog grid uses the schema's friendly names, so it matches `Browse > Platform > Domains`.
`ForgeApp.create()` reads them with `grok.dapi.domains.registry.rowProperties('forge.model')` and sets
`column.meta.friendlyName`; the registry derives captions for columns without one (`engine_kind` -> "Engine
kind"). The system column `created_on` is labelled "Created" explicitly. The grid shows only the catalog columns, in
order (`grid.columns.setOrder` and `setVisible`); `id`, `version`, `updated_on` and `author_id` stay hidden.

## Engine contract

An engine is the set of functions that share `meta.mlname` (the engine id). `meta.mlrole` gives each function's role:

| Role | Signature | Notes |
|---|---|---|
| `train` | `(df, predictColumn, ...hyperparameters) -> model blob` | The other inputs are the hyperparameters |
| `apply` | `(df, model) -> dataframe` | |
| `isApplicable` | `(df, predictColumn) -> bool` | |
| `isInteractive` | `(df, predictColumn) -> bool` | Optional. `meta.mlupdate: 'false'` turns off live retraining |
| `visualize` | `(df, targetColumn, predictColumn, model)` | Optional |

- An engine is complete, and usable, only with `train`, `apply` and `isApplicable`.
- A function engine receives the feature table and the target column. A script engine (`DG.Script`) receives one
  table with a copy of the target appended, and the target name.
- `EngineRegistry.discover()` calls `DG.Func.find({meta: {mlrole}})` once per role, groups by `mlname`, takes kind and
  namespace from the `train` function, and lists function engines first, then scripts, each by name.
  If a package is registered twice, the later function of a role wins.
- The namespace is the package name for function engines (EDA registers as `Eda`) and the `nqName` prefix for scripts.
- `isApplicable` on an engine without that role throws `ForgeError`; `isInteractive` without it returns `false`.
- Discovery reads the client function registry synchronously; script engines that are not loaded into it yet are
  missed. No script engines exist yet to verify this.
