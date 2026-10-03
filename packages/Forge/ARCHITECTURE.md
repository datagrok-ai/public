# Forge architecture

How Forge works now. Forge is predictive modeling as a package: engines that live in other packages, models stored in
an EMS schema owned by Forge, and the UI on top.

## Key decisions

- **Own schema instead of the platform model entity.** Models, training runs and applications are rows of the EMS
  schema `forge`. That gives typed storage, row security and an audit trail without core changes.
- **Engines stay in their packages.** Forge calls them through the `mlname`/`mlrole` contract and copies no engine code.
- **"Method" for users, "engine" in code.** Captions, messages and docs for users say method; code, the `engines/`
  folder, the `engine_*` columns and the `mlname` contract keep engine.
- **Two tools coexist.** The built-in `ML | Models` tool keeps working until the switchover. Forge uses its own
  identifiers (`ML | Forge` menu, app, schema) and never touches the old tool's.
- **Training data never leaves the browser without consent.** Only the `copy` storage mode uploads it; this version
  saves every model in `none` mode, with a fingerprint of the data instead of the data.
- **Every training attempt is recorded.** Completed, failed and cancelled runs are `training_run` rows; a model row
  appears only on Save.
- **Application records inherit the model's security.** A record is visible to whoever can see its model. It can be
  written by anyone the `model` table lets insert, which today is every user through the table's Edit grant, so
  applications by users who can only view a model are recorded too.
- **Later:** applying models, more engines, MLflow models and the migration of existing models.

## Components and data flow

```
package.ts ──> ui/forge-app.ts ───> engines/                     (DG.Func registry)
           │                    └─> storage/, generated/db.ts    (catalog query, deleteModel, modelsChanged)
           └─> ui/train-view.ts ──> training/ ──> engines/, metrics/
                                └─> storage/  ──> generated/db.ts  (grok.dapi.domains -> EMS schema forge)
                                              └─> grok.dapi.files  (System:DomainFiles/forge/model)
```

- `package.ts` holds the annotated functions only: the app `forgeApp` (shown as "Forge", `Browse > Apps`),
  `forgeModels` (`ML | Forge | Models`) and `forgeTrain` (`ML | Forge | Train...`). They delegate to `ForgeApp` and
  `TrainView`.
- `ForgeApp` lists the methods (section **Methods**, columns Method, Package, Method type, Roles, Hyperparameters)
  and the catalog. `ForgeApp.loadModels()` queries the catalog columns and sets the captions; the refresh icon reloads
  the grid. The trash icon (**Delete model**) is disabled while the catalog has no current row (re-evaluated on
  current-row change and on every reload; the subscription to the current frame is replaced on reload and dropped in
  `detach()`); it confirms and calls `deleteModel`. The two icons sit on the baseline of the **Models** title
  (`forge-pane-header` in `css/forge.css`).
- `modelsChanged` (an rxjs `Subject` in `storage/model-store.ts`) fires after every successful `saveModel` and
  `deleteModel`, whoever calls them; every open `ForgeApp` subscribes in `this.subs` and reloads itself.
- `TrainView` (tab **Predictive model** with the core's model icon, `ui.iconSvg('model')` returned by `getIcon()`)
  holds the form, the live validation and the results. The view keeps the latest training in `lastTraining`:
  `result`, `runId`, the `datasetName` and `fingerprint` taken at training time, and `isSaved`. **Save** (`ui.bigButton` in the view's
  ribbon) is disabled until a training result exists, while a training runs, and once that result is saved: OK in
  `ui/save-model-dialog.ts` (Name, Description) calls `saveModelAs`, which marks the training saved before it
  writes, so one training becomes at most one model; a failed write enables **Save** again.
- Table, Target, Features and Method carry short Forge tooltips (`tooltipText`, the caption's tooltip). Hyperparameter
  inputs come from `ui.input.forProperty`, whose caption tooltip is the `train` input's own `description`. For an
  invalid input the platform shows the tooltip text with the validator messages in red below it.
- UI entry points and handlers are the only error boundaries (`ui/report-error.ts`): `ForgeError` becomes a warning,
  anything else an error balloon plus `_package.logger.error`.
- Logic folders (`engines/`, `training/`, `metrics/`, `storage/`; `preparation/` and `apply/` later) work without the
  DOM, throw instead of catching, and never import `ui/`. The same logic serves the UI, API functions and tests.

Next stage: **Apply** - match the columns, call the engine's `apply`, record an application.

## Training pipeline

`training/train-model.ts`. A `TrainingRequest` is the engine, the feature frame (the checked columns), the target
column, the hyperparameters, a seed and the number of folds (5).

**Defaults in the view.** Target: the table's last column. Features (`defaultFeatures`): every numerical column
except the target and except integer columns without missing values whose values are all different (row numbers and
ids, such as the `col 1` row-number column of `iris.csv`, which is sorted by species). Float columns are never
excluded. Hyperparameters: `defaultHyperparameters(engine)`, read from the `train` inputs' initial values. Changing
the target leaves the features as they are.

**Problems.** `trainingProblems(request)` returns the messages per input; `checkTrainable` throws them as one
`ForgeError`. Rules 1-7 are cheap checks; rule 8 calls the engine and runs only when 1-7 pass:

| # | Condition | Input |
|---|---|---|
| 1 | no feature | Features |
| 2 | the target is also a feature | Features |
| 3 | fewer rows than `2 * folds` (10) | Target |
| 4 | the target has empty values | Target |
| 5 | a classification target with fewer than two values | Target |
| 6 | the target is neither numerical, text nor boolean | Target |
| 7 | a feature is not numerical (names listed) | Features |
| 8 | the engine's `isApplicable` says no | Features |

The view adds a validator to **Target** and **Features** that returns the view's current problems for that input.
Every change of Table, Target or Features numbers a new check, clears the problems, disables **Train** (a tooltip already
shown on it then reads "Checking the selection...") and schedules the check 300 ms later; a check that finishes after a newer change is discarded, so
**Train** stays disabled while a check is pending. After a check the inputs are validated again (red input, messages
in its tooltip) and **Train** stays disabled with the first message as its tooltip while any problem exists or a
training runs (one training per view at a time). The tooltip overlay of the disabled button is built once per
disabled period. A check that fails keeps **Train** disabled with the error as its tooltip and shows each distinct
error once. Platform validators do not run on empty values, so an empty Features shows the platform's "Can't be
empty". **Train** runs `checkTrainable` again before training. A rejected selection writes nothing.

**Task.** A numerical target is regression; text and boolean targets are classification.

**Validation.** `kFold(rowCount, folds, seed)` shuffles the rows with a mulberry32 generator (Fisher-Yates) and
deals them round-robin into the folds, so fold sizes differ by at most one and every row is validated exactly once.
For each fold the engine trains on the other folds and predicts the held-out rows; the predictions form one
out-of-fold column, and the **Validation** metrics are computed once on it. The **Train** metrics come from the final
model, trained on all rows and applied to them. No stratification: a tiny table can leave a class out of a fold, and
the engine then fails. The UI draws a random 31-bit seed per training and shows it; the seed and
`splitting = {scheme: 'kfold', folds: 5, isStratified: false}` are saved.

**Progress and cancel.** `trainModel(request, progress)` reports "Fold i of 5" after each fit and throws
`ForgeError('Training was cancelled.')` before a fit when the progress indicator was cancelled. A running fit cannot
be interrupted.

**Run records.** `TrainView` writes the `training_run` row right after `trainModel` settles: `completed` with the
metrics, `cancelled` when the indicator was cancelled, `failed` with the error message otherwise. `model_id` is set by
`linkTrainingRun` when the model is saved; deleting the model keeps the run (`setnull`).

The result (`TrainingResult`) carries the blob, the task, the metrics, the target and feature schemas, the
preparation options (empty for now), the splitting, the seed, the hyperparameters and the row count.

## Metrics

`metrics/metrics.ts`, ported from the built-in tool's definitions; labels and descriptions are the built-in ones,
values are shown with 3 decimals.

- Regression: `mse` = mean of squared errors, `rmse` = its root, `mae` = mean of absolute errors (new),
  `r2` ("R squared") = 1 - SSres / SStot; for a constant target 1 if the predictions are exact, else 0.
- Classification: labels compare as text (a boolean `true` equals the label `'true'`). `accuracy` = correct / n.
  When at most two labels occur, the shares are computed for the positive class, the target's first category:
  `sensitivity` = TP / (TP + FN), `specificity` = TN / (TN + FP), `precision` = TP / (TP + FP),
  `npv` ("Negative Predicted Value") = TN / (TN + FN); a share whose denominator is 0 is left out, as the built-in
  tool hid it. `f1` (new): of the positive class for two classes, the macro average over the classes otherwise.
  For a two-category target the positive class is saved as `metrics.positiveClass`.
- Rows where the actual or the predicted value is empty are skipped; no rows left is a `ForgeError`.
- AUC-ROC is not computed yet: the engines' `apply` returns labels, not scores.

## Saving

`storage/model-store.ts`. `saveModel(fields, blob)` writes the blob to
`System:DomainFiles/forge/model/<uuid>/model.bin` first and then inserts the `model` row with
`blob = file://<that path>`. A failed insert leaves an unreachable file, never a row without its blob.
`deleteModel(id)` soft-deletes the row; it deletes a folder only when `blob` matches
`file://System:DomainFiles/forge/model/<uuid>/<file>`, and then exactly `System:DomainFiles/forge/model/<uuid>`. Any
other value (a file picked in the generic row editor elsewhere, an empty segment, `..`, no blob) leaves the files
alone. In the UI, the ribbon **Save** asks for Name and Description, then calls `saveModel` and `linkTrainingRun`;
the catalog's delete calls `deleteModel`. Both storage functions fire `modelsChanged`.

`model-fields.ts` builds the row payloads: `modelFieldsOf` (`storage_mode: 'none'`, `has_training_rows: false`,
feature and row counts, dataset name) and `trainingRunOf`. The json-column types are declared in
`training/train-model.ts`.

**Dataset fingerprint** (`dataset-fingerprint.ts`), the only trace of the training data on the server, computed in
the browser over the features in training order and the target last: `rowCount`, `columnCount`, an 8-hex-digit hash
and per column `name`, `type`, `missingCount`, plus `min`/`max`/`mean` for numerical columns and `categories` for
text and boolean columns with up to 20 of them. The hash is FNV-1a (32 bit) over each column's name and type and
then its data: the raw buffer cut to the column length for numerical columns, the values as text otherwise (bigint
columns as text). It depends on the column order.

## What differs from the built-in tool

The form starts ready to train: the last column is the target and the numerical columns except row numbers and ids
are the features, where the built-in tool starts with nothing selected. Problems are shown live on the input they
concern and **Train** stays unavailable until they are fixed. Validation uses five folds that cover every row once,
from a seed saved with the model, instead of five overlapping random samples, and it runs from 10 rows instead of
being skipped below 100. MAE and F1 are new; AUC-ROC is not computed yet. Training and saving are separate steps, and
every attempt is recorded. The training table is not uploaded: a model keeps a fingerprint of it.

## EMS schema `forge`

Manifest: `databases/forge/schema.json`, version `0.1.4`. `grok publish` deploys it; a debug publish applies
destructive changes without migration scripts, except a `promotion` change on a row table, which is always refused
(`[promotion-change]`). To change it, publish a manifest with only a placeholder table, then the real one, bumping
`version` both times; the first publish drops the tables and their data. `grok api` generates the typed client
`src/generated/db.ts` (`forgeDb.models`, `forgeDb.trainingRuns`, `forgeDb.applications`). Every table also has the
system columns `id`, `version`, `created_on`, `updated_on` and `author_id`; who and when always come from them.

| Table | Security | Why |
|---|---|---|
| `model` | `row`, `promotion: lazy`, `defaultRowVisibility: none`, grants `All users: view, edit` | A model is private until shared, like the old models. Lazy promotion makes a model a platform entity on its first share; sharing, favorites and comments work from then on. The table Edit grant lets any user insert models; with visibility `none` it reveals no foreign rows |
| `training_run` | `row`, lazy promotion, `defaultRowVisibility: none`, same grants | The trainer's experiment history, including failed and cancelled runs. `model_id` is optional with `onDelete: setnull`, so the history outlives the model |
| `application` | `master`, `delegate: model_id`, `audit: false`, no grants | Secured by the model (see Key decisions). Master tables refuse grants. No audit: high-churn records. `onDelete: cascade` from the model |

Deletes are soft: a deleted row stays with `is_deleted` and disappears from queries and counts. Current EMS limitation:
a soft delete does not clean up a promoted row's entity and permissions, so a model deleted after being shared leaves
a live entity until the platform fixes it. The generic Domains view offers **New model...** to anyone with the
table's Edit grant; such a row has no engine blob, and Forge lists it and can delete it (the row only, unless its
file sits in Forge's own `forge/model/<uuid>/` layout).

### JSON column shapes

EMS `json` columns hold objects only; a list is wrapped in an object. TS types: `training/train-model.ts`
and `storage/dataset-fingerprint.ts`.

| Column | TS type | Shape |
|---|---|---|
| `model.target` (caption "Target details") | `TargetSchema` | `{name, type, semType?, categories?}`; `categories` for classification targets |
| `model.features`, `training_run.features` | `FeaturesSchema` | `{columns: [{name, type, semType?}]}`, in training order |
| `options` | `PreparationOptions` | Preparation replay, see below; `{preprocessingInfo: [], postprocessingInfo: []}` for now |
| `hyperparameters` | `Hyperparameters` | `{<train function input>: value}` |
| `metrics` | `MetricsRecord` | `{train: {<metric id>: number}, validation: {<metric id>: number}, positiveClass?}`; ids `mse`, `rmse`, `mae`, `r2`, `accuracy`, `f1`, `sensitivity`, `specificity`, `precision`, `npv` (`auc` reserved) |
| `splitting` | `Splitting` | `{scheme: none \| kfold \| holdout, folds?, trainFraction?, isStratified?}` |
| `model.dataset_ref` | - | Reference mode: `{kind: file \| query \| script, id?, path?, params?, script?}` |
| `dataset_fingerprint` | `DatasetFingerprint` | `{rowCount, columnCount, hash, columns: [{name, type, missingCount, min?, max?, mean?, categories?}]}` |

`options` keeps the keys of the platform's built-in models, so their preparation replays unchanged:
`preprocessingInfo`, `postprocessingInfo`, `positiveClass`, `negativeClass`, `binaryClassificationThreshold`,
`targetType`, `allowNulls`.

### Model blob

The `model.blob` column (type `file`) stores `file://System:DomainFiles/forge/model/<uuid>/model.bin`, the layout the
EMS row editor uses for file columns. The random path segment makes the file unguessable, and the row security of
`model` protects the pointer. The generic Domains UI shows the column but offers no download of the file.

### Storage modes

`model.storage_mode` (default `none`) records what a model keeps of its training data:

- `none`: only `dataset_fingerprint`, enough to check that a new table fits the model. The only mode used so far.
- `reference`: `dataset_ref` points to the source (file, query or script); the data is not uploaded.
- `copy`: `dataset_table_id` is the id of an uploaded copy of the training table (no foreign key; checked on read).

`has_training_rows` marks blobs that embed training rows (SVM, KNN); `false` for XGBoost. `legacy_id` is reserved for
migrated models.

### Captions

Every column of `model` and `training_run` declares a `friendlyName`, because the platform renders undeclared ones
differently in different places: the Domains **grid** (`Browse > Platform > Domains`) shows the raw name with
underscores ("Engine_kind") and no json columns at all, while the row's **context panel** shows every column with a
caption the client registry derives by splitting the name ("Engine kind") and json values as raw text. The server
drops a declared caption that equals the capitalized name (Name, Task, Features, ...); those render the same anyway.

| Column | Caption | Column | Caption |
|---|---|---|---|
| `name` | Name | `storage_mode` (model) | Data storage |
| `description` (model) | Description | `dataset_name` | Table |
| `engine_name` | Method | `row_count` | Training rows |
| `engine_namespace` | Package | `dataset_ref` (model) | Data source |
| `engine_kind` | Method type | `dataset_table_id` (model) | Data copy |
| `task` | Task | `dataset_fingerprint` | Data summary |
| `target_name` | Target | `blob` (model) | Model file |
| `target` (model) | Target details | `has_training_rows` (model) | Contains training rows |
| `features` | Features | `legacy_id` (model) | Legacy id |
| `feature_count` (model) | Feature count | `model_id` (run) | Model |
| `options` | Preparation | `status` (run) | Status |
| `hyperparameters` | Hyperparameters | `error` (run) | Error |
| `metrics` | Metrics | `started_on` (run) | Started |
| `seed` | Seed | `duration_ms` (run) | Duration (ms) |
| `splitting` | Data split | | |

The `model` filter on `engine_name` is labelled **Method**. The catalog grid reads the same captions with
`grok.dapi.domains.registry.rowProperties('forge.model')` in `ForgeApp.loadModels()` and sets
`column.meta.friendlyName`; the system column `created_on` is labelled "Created" explicitly. The grid shows only the
catalog columns, in order (`grid.columns.setOrder` and `setVisible`); `id`, `version`, `updated_on` and `author_id`
stay hidden. "No models yet." is shown while the catalog is empty.

## Engine contract (methods in the UI)

An engine is the set of functions that share `meta.mlname` (the engine id). `meta.mlrole` gives each function's role:

| Role | Signature | Notes |
|---|---|---|
| `train` | `(df, predictColumn, ...hyperparameters) -> model blob` | The other inputs are the hyperparameters |
| `apply` | `(df, model) -> dataframe` | Forge takes the first column as the prediction |
| `isApplicable` | `(df, predictColumn) -> bool` | |
| `isInteractive` | `(df, predictColumn) -> bool` | Optional. `meta.mlupdate: 'false'` turns off live retraining |
| `visualize` | `(df, targetColumn, predictColumn, model)` | Optional |

- An engine is complete, and usable, only with `train`, `apply` and `isApplicable`.
- A function engine receives the feature table and the target column. A script engine (`DG.Script`) receives one
  table with a copy of the target appended, and the target name.
- `train` returns the blob as a `Uint8Array` or a `DG.FileInfo` (its `data` is used); anything else is a
  `ForgeError`. `apply` must return a dataframe with at least one column.
- Hyperparameter defaults come from the `train` annotation (`//input: int iterations = 20`). The platform keeps them
  as text in `Property.initialValue` (`Property.defaultValue` stays null); `defaultHyperparameters` converts them by
  the input type.
- `EngineRegistry.discover()` calls `DG.Func.find({meta: {mlrole}})` once per role, groups by `mlname`, takes kind and
  namespace from the `train` function, and lists function engines first, then scripts, each by name.
  If a package is registered twice, the later function of a role wins.
- The namespace is the package name for function engines (EDA registers as `Eda`) and the `nqName` prefix for scripts.
- `isApplicable` on an engine without that role throws `ForgeError`; `isInteractive` without it returns `false`.
- Discovery reads the client function registry synchronously; script engines that are not loaded into it yet are
  missed. No script engines exist yet to verify this.
