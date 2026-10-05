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
- **Every application is recorded.** Completed, failed and cancelled applications are `application` rows; a mapping
  the model cannot use is refused before any work and not recorded.
- **Application records inherit the model's security.** A record is visible to whoever can see its model. It can be
  written by anyone the `model` table lets insert, which today is every user through the table's Edit grant, so
  applications by users who can only view a model are recorded too.
- **No copies of the user's data.** Training and application work on frames that share the user's columns
  (`preparation/shared-frame.ts`: `sharedFrame`); a column is copied only where Forge would otherwise change it (see
  "Missing values" and "Application pipeline"). A frame stays a parent of its columns and receives their events, so
  every such frame is given back with `releaseFrame` (its columns removed) as soon as the work is done, and the
  live checks work on column lists.
- **Missing values never reach a method.** Rows with gaps are skipped or the gaps are imputed before the engine is
  called; dates and very large whole numbers (bigint) are not numerical features.
- **Later:** more engines, MLflow models and the migration of existing models.

## Components and data flow

```
package.ts ──> ui/forge-app.ts ───> engines/                     (DG.Func registry)
           │                    ├─> storage/, generated/db.ts    (catalog query, deleteModel, modelsChanged)
           │                    └─> ui/apply-model-dialog.ts
           ├─> ui/train-view.ts ──> training/ ──> preparation/, engines/, metrics/
           │                    └─> storage/  ──> generated/db.ts  (grok.dapi.domains -> EMS schema forge)
           │                                  └─> grok.dapi.files  (System:DomainFiles/forge/model)
           ├─> ui/apply-model-dialog.ts ──> apply/, engines/, preparation/, generated/db.ts (model query)
           └─> apply/ (runApplyModel) ──> preparation/, engines/, storage/ (blob path), generated/db.ts,
                                          grok.dapi.files (blob read)
```

- `package.ts` holds the annotated functions only: the app `forgeApp` (shown as "Forge", `Browse > Apps`),
  `forgeModels` (`ML | Forge | Models`), `forgeTrain` (`ML | Forge | Train...`), `forgeApply` (`ML | Forge |
  Apply...`) and the API function `applyModel`. They delegate to `ForgeApp`, `TrainView`, `openApplyDialog` and
  `runApplyModel`.
- `ForgeApp` lists the methods (section **Methods**, columns Method, Package, Method type, Roles, Hyperparameters)
  and the catalog. `ForgeApp.loadModels()` queries the catalog columns and sets the captions; the refresh icon reloads
  the grid. The play icon (**Apply model**) and the trash icon (**Delete model**) are disabled while the catalog has
  no current row (re-evaluated on current-row change and on every reload; the subscription to the current frame is
  replaced on reload and dropped in `detach()`). The play icon opens the Apply dialog for the current row's model,
  with the current table (the catalog is not a table view, so in practice the first open table; without tables,
  "Open a table first."), and shows that table's view after applying. The trash icon confirms and calls
  `deleteModel`. The icons sit on the baseline of the **Models** title (`forge-pane-header` in `css/forge.css`); the
  grid takes the full width of the view: an inline `width: 100%`, the override the platform's 400px rule for boxes
  in a panel (`ui.css`) asks for.
- `modelsChanged` (an rxjs `Subject` in `storage/model-store.ts`) fires after every successful `saveModel` and
  `deleteModel`, whoever calls them; every open `ForgeApp` subscribes in `this.subs` and reloads itself.
- `TrainView` (tab **Predictive model** with the core's model icon, `ui.iconSvg('model')` returned by `getIcon()`)
  holds the form, the live validation and the results. The form has two groups, **Data** (Table, Target, Features,
  Missing values with Neighbors and Distance) and **Method** (Method and the hyperparameters), with **Train** below
  them, right-aligned with the 350px column inputs (`forge-train-row`); both start expanded. A Table change rebuilds
  the form; the groups' subscriptions (`formSubs`) are dropped then and in `detach()`. The view keeps the latest
  training in `lastTraining`: `result`, `runId`, the `datasetName` and `fingerprint` taken at training time, and
  `isSaved`. **Save** (`ui.bigButton` in the view's ribbon) is disabled until a training result exists, while a
  training runs, and once that result is saved: OK in `ui/save-model-dialog.ts` (Name, prefilled, and Description)
  calls `saveModelAs`, which marks the training saved before it writes, so one training becomes at most one model; a
  failed write enables **Save** again. **Name** is not nullable and its validator refuses a blank name ("Enter a name
  for the model."), so OK is disabled with that tooltip and the dialog's own validation blocks Enter too.
- `CollapsibleGroup` (`ui/collapsible-group.ts`, classes `forge-group*` in `css/forge.css`) is the folding block of
  the Train view groups and of the Apply dialog's **Columns** and **More options**: a header with a chevron, the
  caption and an optional summary (red with `forge-group-invalid`), in Diff Studio's style; the chevron sits in a
  fixed-width box (`forge-group-chevron`), so captions line up whether a group is open or folded. `expandOnError(inputs)`
  opens the group when one of its inputs is validated with an error; the state is not remembered.
- Table, Target, Features and Method carry short Forge tooltips (`tooltipText`, the caption's tooltip). Hyperparameter
  inputs come from `ui.input.forProperty`, whose caption tooltip is the `train` input's own `description`. For an
  invalid input the platform shows the tooltip text with the validator messages in red below it.
- **Missing values** (`ui/missing-values-inputs.ts`, `MissingValuesInputs`, used by the Train view and the Apply
  dialog) is a radio of **Skip rows** / **Impute** with its options on one line (`forge-inline-radio`), shown only
  when its columns have gaps, with the tooltip "Rows with missing values in: <col>: N, ...". **Impute** is offered
  only when `Eda:knnImpute` exists; choosing it shows the function's own inputs (**Neighbors**, **Distance**, from
  `imputeSettingsOf` with `defaultValuesOf`) right under it. Hiding them puts their defaults back, so a hidden invalid
  value cannot block the form.
- UI entry points and handlers are the only error boundaries (`ui/report-error.ts`): `ForgeError` becomes a warning,
  anything else an error balloon plus `_package.logger.error`. API functions do not catch: the error reaches the
  caller.
- Logic folders (`engines/`, `preparation/`, `training/`, `metrics/`, `storage/`, `apply/`) work without the DOM,
  throw instead of catching, and never import `ui/`. The same logic serves the UI, API functions and tests. The one
  catch in them is `applyAndRecord`, which records a failed or cancelled application and rethrows.

## Training pipeline

`training/train-model.ts`. A `TrainingSelection` is what the user chose: the engine, the checked feature columns
themselves (a list, no frame), the target column, the hyperparameters, a seed, the number of folds (5) and the
missing-value handling. `prepareTraining(selection)` builds a shared frame of the columns, handles the missing values
(see "Missing values") and returns the `TrainingRequest` the model is trained on: the prepared features and target,
and the `options` that record the preparation (`options.missingValues.skippedRows` is the one record of the skipped
rows); the frames it does not return are given back at once.
The request's `features` may still share the user's columns: the Train view releases them right after training and
recording the run. After preparation fewer than `2 * folds` rows is a `ForgeError` ("After skipping rows with missing
values, N rows remain; training needs at least 10.").

**Missing values in the view.** **Missing values** (see Components) sits under **Features** and follows the checked
features. The **Target** tooltip adds "<target>: N rows without a value are skipped." ("1 row without a value is
skipped.") when the target has gaps; it is information, not a problem. After a training with skipped rows,
**Results** shows "Rows: N used, M skipped (missing values)." above the validation line.

**Defaults in the view.** Target: the table's last column. Features (`defaultFeatures`): every column the methods
read as numbers (`isReadableNumber`: numerical, not dates, not bigint) except the target and integer columns without
missing values whose values are all different
(row numbers and ids, such as the `col 1` row-number column of `iris.csv`, which is sorted by species). Float columns
are never excluded. Hyperparameters: `defaultHyperparameters(engine)`, read from the `train` inputs' initial values.
Changing the target leaves the features as they are.

**Problems.** `trainingProblems(selection)` returns the messages per input (`target`, `features`, `missingValues`);
`checkTrainable` throws them as one `ForgeError`. "Numerical" means `numerical_no_datetime`: dates are not numbers
for Forge. Rules 1-7 and 9 are cheap checks on the column list; rule 8 calls the engine, whose check takes a table, on
a shared frame given back right after the call, and runs only when 1-7 pass:

| # | Condition | Input |
|---|---|---|
| 1 | no feature | Features |
| 2 | the target is also a feature | Features |
| 3 | fewer rows with a target value than `2 * folds` (10) | Target |
| 4 | a bigint target ("holds very large whole numbers", `bigIntProblem` in `training/default-features.ts`) | Target |
| 5 | a classification target with fewer than two values | Target |
| 6 | the target is neither numerical, text nor boolean (dates included) | Target |
| 7 | a feature is not numerical (names listed; dates included) | Features |
| 7b | a bigint feature, one message per column | Features |
| 8 | the engine's `isApplicable` says no | Features |
| 9 | Impute chosen but not available, a yes/no feature with gaps, or fewer than two features to impute from | Missing values |

Rows with a missing target are skipped by `prepareTraining`, never refused.

The view adds a validator to **Target**, **Features** and **Missing values** that returns the view's current problems
for that input. Every change of Table, Target, Features, Missing values or the imputation inputs numbers a new check,
clears the problems, disables **Train** (a tooltip already shown on it then reads "Checking the selection...") and
schedules the check 300 ms later; a check that finishes after a newer change is discarded, so **Train** stays disabled
while a check is pending. After a check the inputs are validated again (red input, messages in its tooltip) and
**Train** stays disabled with the first message as its tooltip while any problem exists (an invalid **Neighbors**
included) or a training runs (one training per view at a time). A `ButtonGate` (`ui/button-gate.ts`, shared with the
Save and Apply dialogs' OK) builds the tooltip overlay of the disabled button once per disabled period. A check that
fails keeps **Train** disabled with the error as its tooltip and shows each distinct
error once. Platform validators do not run on empty values, so an empty Features shows the platform's "Can't be
empty". **Train** runs `checkTrainable` again before training. A rejected selection writes nothing.

**Task.** A numerical target (not a date) is regression; text and boolean targets are classification.

**Validation.** `kFold(rowCount, folds, seed)` shuffles the rows with a mulberry32 generator (Fisher-Yates) and
deals them round-robin into the folds, so fold sizes differ by at most one and every row is validated exactly once.
For each fold the engine trains on the other folds and predicts the held-out rows; the predictions form one
out-of-fold column, and the **Validation** metrics are computed once on it. The **Train** metrics come from the final
model, trained on all rows and applied to them. No stratification: a tiny table can leave a class out of a fold, and
the engine then fails. The UI draws a random 31-bit seed per training and shows it; the seed and
`splitting = {scheme: 'kfold', folds: 5, isStratified: false}` are saved.

**Progress and cancel.** `trainModel(request, progress?)` (a `LoopProgress`, the cancel flag and `update` of a
progress indicator, in `engines/engine-calls.ts`) reports "Fold i of 5" after each fit and throws
`ForgeError('Training was cancelled.')` before a fit when the progress indicator was cancelled. Before that check it
gives the event loop a turn (`yieldToEventLoop` in `engines/engine-calls.ts`): an engine call that finishes without
I/O resumes as a microtask, so without the turn a click on the progress's cancel (and the progress repaint) would wait
for the whole training. The turn is a `MessageChannel` message, which a background tab does not throttle the way it
throttles timers (to about one a second). A running fit cannot be interrupted.

**Run records.** `TrainView` writes the `training_run` row right after `trainModel` settles: `completed` with the
metrics, `cancelled` when the indicator was cancelled, `failed` with the error message otherwise. `model_id` is set by
`linkTrainingRun` when the model is saved; deleting the model keeps the run (`setnull`).

The result (`TrainingResult`) carries the blob, the task, the metrics, the target and feature schemas, the
preparation options, the splitting, the seed, the hyperparameters and the row count (the prepared rows the model
learned from); the skipped rows are in `options.missingValues`. The run and model rows record that row count, and the
dataset fingerprint is computed on the prepared data.

## Missing values

`preparation/missing-values.ts`. `prepareMissingValues(features, target, settings)` returns `PreparedData`: the
features (and target) to use, `keptRows` (the input rows kept, or `null` when every row was kept), `skippedRows`,
`imputedColumns` and `failedRows`. Rows with a missing target are always skipped, and the target is never imputed.
`settings.mode` is one of `MISSING_VALUES_MODES`:

- **`skip`** (Skip rows): rows with a missing value in any feature (or the target) are skipped. Only then are the
  features and the target copied (`rowCopy` / `columnRowCopy` in `preparation/row-copy.ts`); with nothing to skip
  the input frame itself is returned. A masked clone of a text column keeps every category of the source, the empty
  one included, so the copies drop the categories no kept row has: otherwise the method would count a class for the
  skipped empty cells, and the saved target and data summary would list it. The imputed text copies are compacted the
  same way. The k-fold copies in `trainModel` are not: they keep the target's full category list, so every fold model
  has the classes of the final model.
- **`impute`** (Impute): the rows with a missing target are skipped first; then the feature gaps are filled by EDA's
  `Eda:knnImpute` (k nearest neighbors, in place). Because it writes in place, it gets a frame whose gapped columns are
  copies and whose other columns are still the user's (no extra copy when the target skip already made one). Cells it
  cannot fill (a row whose features are all missing) stay empty; those rows are then skipped and counted in
  `failedRows` and `skippedRows`. `settings.impute` holds `neighbors` and `distance`; `imputeSettingsOf(func)` lists
  the function's own inputs for them, and `defaultValuesOf` reads their defaults.

The returned `features` may share the user's columns (the input frame itself, or the imputation frame); the caller
releases both its input frame and the returned one. An imputation frame that is not returned is released inside.

`missingValuesProblems(columns, settings)` refuses Impute when the function is missing ("Imputation is not
available: the EDA package has no knnImpute function."), for a yes/no column with gaps ("Impute cannot fill the yes/no
column 'X'. Choose Skip rows.") and with fewer than two features the imputer can read ("Impute needs at least two
features. Choose Skip rows or check more features."); it says nothing for Skip rows or for features without gaps.
`missingColumnsOf(columns)` lists the columns with gaps and their counts. Both take column lists, so the UI checks
build no frame.

At training, `prepareTraining` records the handling in `options`: `preprocessingInfo` gains `impute-missing` when
imputation ran and `ignore-missing` when a feature has gaps under Skip rows, or imputation left rows unfilled (the ids
of the built-in tool), and
`missingValues` is `{mode, neighbors?, distance?, skippedRows}`. Both are information only: an application never
replays them and chooses its own handling. Imputing before cross-validation lets the held-out rows' values take part
in the imputation, as the built-in tool did; a later version may impute per fold.

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
feature and row counts, dataset name) and `trainingRunOf`. The json-column types are listed in "JSON column shapes"
below.

**Dataset fingerprint** (`dataset-fingerprint.ts`), the only trace of the training data on the server, computed in
the browser over the features in training order and the target last: `rowCount`, `columnCount`, an 8-hex-digit hash
and per column `name`, `type`, `missingCount`, plus `min`/`max`/`mean` for numerical columns and `categories` for
text and boolean columns with up to 20 of them. The hash is FNV-1a (32 bit) over each column's name and type and
then its data: the raw buffer cut to the column length for numerical columns, the values as text otherwise (bigint
columns as text). It depends on the column order.

## Application pipeline

`apply/apply-model.ts`. `loadModel(idOrName)` reads, in one query, the `APPLY_COLUMNS` of the model with that id, or
of the one visible model with that name ("No Forge model 'X' is available to you." / "Several models are named 'X'.
Use the model id."). `loadedModelOf(row, engines)` then refuses, in this order, a row without a blob ("has no model
file"), a blob outside Forge's own `forge/model/<uuid>/` layout (`ownBlob` in `storage/model-store.ts`), a row
without a `{name, type}` feature list, and an engine that is not discovered or incomplete ("The method 'X' is not
installed. Install the <package> package."). The result (`LoadedModel`) holds the row, the engine, the features, the
options (`preparationOptionsOf`) and the blob path; the target is the row's `target_name`.

An `ApplyRequest` is the model, the table, the column mapping, the batch size (`DEFAULT_BATCH_SIZE` 10000) and the
missing-value settings. `applyAndRecord(request, source, progress?)`:

1. Refuses the request when `mappingProblems` reports anything or the table has no rows; a refusal is not recorded.
2. Builds the feature frame (`featureFrame`, `apply/feature-frame.ts`): the mapped columns in training order, shared
   with the table. A column whose name differs from its feature is a renamed copy: renaming a shared column would
   rename it in the user's table, and the engines find features by name.
3. `prepareMissingValues` without a target (see "Missing values"); every row skipped is a `ForgeError` ("Every row
   has a missing value in the columns the model needs, so nothing can be predicted. Fill the missing values or
   choose Impute.").
4. `replayPreprocessing` (`preparation/preparation-steps.ts`) replays `options.preprocessingInfo` on a frame of the
   same columns: `one-hot` replaces every text and boolean column with `<column>=<category>` 0/1 columns (the
   categories of the applied data), `skip-unique-categories` removes text columns whose values are all different;
   `ignore-missing` and `impute-missing` are skipped; any other id is refused ("The model uses the preparation step
   'X', which Forge cannot replay yet."). Without steps the frame itself is used.
5. Reads the blob and calls the engine's `apply`: once on the frame itself when it fits one batch, otherwise once per
   batch of rows (a copy each, of a range mask built from a word array), the results written into one column
   (numerical predictions through their raw data). Before every call a cancelled progress
   indicator stops it ("Application was cancelled."); after every call the progress shows "Rows a-b of n". At most
   every 50 ms the loop gives the event loop a turn (`yieldToEventLoop`), so a click on the cancel and the progress
   repaint get through: fast engine calls resume as microtasks and would otherwise hold the page to the last batch.
6. `replayPostprocessing` replays `binary-classification`: scores at or above `binaryClassificationThreshold` become
   `positiveClass`, the others `negativeClass`, in a column of `targetType`.
7. With skipped rows, the predictions are spread over a full-length column; skipped rows stay empty.
8. Names the column `predictionName(table, target)`: `<target> (predicted)`, or `<target> (predicted 2)`, ... when a
   column of that name exists (ignoring case); tags it `forge.model` = the model id (`PREDICTION_TAG` in
   `constants.ts`) and adds it to the table.

No step writes into a shared column: the engines read their inputs and return new frames, the replay steps add and
remove columns of their own frame, and the imputer only gets copies of the columns it fills. Every frame built on the
way (feature frame, prepared frame, replay frame) is released in `finally`, whatever happens. The raw engine call
(`apply` in `engines/engine-calls.ts`) takes the prediction column out of the engine's result frame, so the user's
table is the column's only parent (`column.dataFrame` is the table). Only the first column of the engine's result is
used; other output columns are dropped.

After the pipeline, `applyAndRecord` writes one `application` row (`apply/application-store.ts`): `completed` with
`column_name` and `skipped_rows`, `cancelled` when the progress indicator was cancelled, `failed` with the error
otherwise; always `model_id`, `table_name`, `row_count`, `source` and `duration_ms`. When the pipeline failed and the
record cannot be written either, the pipeline's error is the one thrown. `applyWithProgress(request, source)` runs
`applyAndRecord` under the cancellable task-bar indicator "Predicting <target>", closed when it ends; the Apply dialog
and `Forge:applyModel` (with `showProgress`) use it.

## Column matching

`apply/column-matching.ts`. A `ColumnMapping` maps a feature name to a table column name.

**Compatibility** (`compatibility(feature, col, method)`), one rule for every path: a numerical feature accepts any
numerical column except bigint ("'X' holds very large whole numbers, which <method> cannot read. Convert the column
to a decimal type or choose another column.") and dates or other kinds ("'X' holds dates but 'F' needs numbers.
Choose a column with numbers."), by the same `isReadableNumber` rule and `bigIntProblem` sentence as training; a semantic type the column lacks is only a hint ("'X' is not marked as <SemType>;
check that it is the same kind of value."). A text, boolean or date feature needs a column of the same type ("'X'
holds numbers but 'F' needs text. Choose a column with text.") and, when the feature has a semantic type, the same one
("'X' is not a <SemType> column, which 'F' needs. Choose a <SemType> column."). Semantic types compare ignoring case;
a feature without one accepts any column of the right kind.

- `exactMapping`: the column named like the feature, ignoring case (`DataFrame.col`, whose name lookup ignores case;
  a table cannot hold two names that differ only in case); without one, every column whose name contains
  `<feature>=` maps to itself (a pre-encoded one-hot column). No compatibility check.
- `suggestMapping`: compatible exact names first; then the Hungarian assignment (ported from the built-in tool's
  `MatchingSolver`) over the remaining features and compatible columns, with `nameDistance` as the cost and no pair
  above `MAX_NAME_DISTANCE` (0.3); when it has no solution, each feature takes the closest free column within the
  threshold. Features without a close column stay unmapped. `nameDistance` is the Jaro-Winkler distance of the
  lower-cased names (Levenshtein for one-character feature names), from `DG.StringUtils`. `isSuggested` is true when
  every feature gets a column.
- `mappingProblems(features, mapping, table, method)`, per feature in order: "Choose a column for 'F'." (unmapped,
  `isUnmapped`), "The column 'X' is no longer in the table. Choose another column for 'F'.", the compatibility error,
  "'X' is also used for 'G'. Choose a different column." (on the later feature).

## Apply dialog

`ui/apply-model-dialog.ts`. `openApplyDialog(table, {modelId?, switchToTable?})` is the boundary of `ML | Forge |
Apply...` (the current table) and of the catalog's play icon (the current table, else the first open one); without a
table it says "Open a table first.". It calls `applyModelDialog({table, modelId?, switchToTable?})`, which reads up
to 10000 visible model rows once (`MAX_MODELS`), only the columns applying needs (`APPLY_COLUMNS`; no models: "No
Forge models yet. Train and save one with ML | Forge | Train... first.") and discovers the engines once. The dialog
**Apply predictive model** has:

- **Table** (`Table to add the prediction to.`): a change re-orders **Model**, keeps the chosen model and re-prefills
  the rows; the table's column add and remove events re-check the rows while the dialog is open, through a change
  event of **Batch size**, so the dialog's own Enter check (refreshed only on input changes) follows too.
- **Model** (`Saved model to apply.`): labels are unique (`modelLabels`): the model name, with ` (YYYY-MM-DD HH:mm)`
  when models share a name, to the second when they share the minute too, then ` #2`, ` #3`, ... in creation order;
  the dialog keys the models by these labels. The models that fit the table (`isSuggested` over the row's
  `featureSchemasOf`, computed once per distinct feature list) come first, newest first in each group; the default is
  the preset model, else the first. A
  model `loadedModelOf` refuses shows its message instead of the rows, and the same message marks **Model**.
- **Columns** (a `CollapsibleGroup`): the summary `N of M matched` (features without a mapping problem out of all) and,
  inside, one column input per feature, captioned with the feature name, in a block (`forge-apply-rows`) that scrolls
  past 40% of the window height (at least 160 px). It starts collapsed when every feature has a valid column, expanded
  otherwise; the summary turns red while any feature has a problem, and a row validated with an error expands it.
  The list offers only columns `compatibility` does not refuse; the prefill is `suggestMapping`; the tooltip is
  "Table column used as the feature 'F' (<kind>)." plus "'X' has N missing values." and a semantic-type hint when
  they apply. The validator is the row's `mappingProblems` message. The dialog keeps every input it was given, so one
  input per feature name is reused across models; an input of a feature the chosen model lacks is emptied and made
  nullable. Folded rows are still validated by the dialog and still gate **OK**. The rows block opens scrolled to the
  top: the dialog focuses the input it was given last on `show()`, so the rows are given to it before the other
  inputs, and the scroll is reset when the rows are rebuilt or **Columns** opens.
- **Missing values** with **Neighbors** and **Distance** (see Components), following the mapped columns; its validator
  is `missingValuesProblems` for Impute. It stays outside the groups.
- **More options** (a `CollapsibleGroup`, collapsed): **Batch size** (10000, 1..100000, `Rows predicted in one step.
  Lower it for heavy methods.`); an invalid batch size expands it.

After every change the dialog recomputes the problems, re-validates the visible inputs and disables **OK** with the
first problem in form order as its tooltip (the model refusal, the rows' messages in feature order, then any other
invalid input). The gaps of the mapped columns are counted once per change and feed both **Missing values** and the
rows' tooltips. **OK** closes the dialog at once (its handler only starts the work, as the built-in dialog did) and
runs `applyWithProgress(..., 'ui')` from the static `ApplyForm.run`, so the running application does not hold the
form; then it
shows "Added the column "C" to T." ("...; N rows skipped (missing values)." when rows were skipped) and, from the
catalog, switches to the table's view. Errors go to `reportError` (a cancel is the yellow "Application was
cancelled.").

## `Forge:applyModel`

`applyModel(model, table, columnNamesMap?, showProgress = true) -> table`, for scripts and other packages, never opens
a dialog: `runApplyModel` resolves the model by id or name, maps the features by exact names, overlays the given
`columnNamesMap` pairs, and refuses the mapping problems before any work, the unmapped features in one sentence first
("The table has no column 'a' the model needs. Map it in columnNamesMap." / "The table has no columns 'a', 'b' the
model needs. Map them in columnNamesMap."). Rows with missing values are skipped; batches of 10000 rows; with
`showProgress` a cancellable task-bar indicator "Predicting <target>". The application is recorded with `source: api`.

## What differs from the built-in tool

The form starts ready to train: the last column is the target and the numerical columns except row numbers and ids
are the features, where the built-in tool starts with nothing selected. Problems are shown live on the input they
concern and **Train** stays unavailable until they are fixed. Validation uses five folds that cover every row once,
from a seed saved with the model, instead of five overlapping random samples, and it runs from 10 rows instead of
being skipped below 100. MAE and F1 are new; AUC-ROC is not computed yet. Training and saving are separate steps, and
every attempt is recorded. The training table is not uploaded: a model keeps a fingerprint of it. Rows with a missing
target are skipped instead of blocking training, and dates and bigint columns are refused as features and targets.
Missing values are a **Missing values** choice (Skip rows / Impute) with the imputation settings inline, instead of
two checkboxes and EDA's own dialog, and the target is never imputed.

Applying lists every visible model, up to 10000 (those that fit first), instead of only the suggested ones, unfolds
the feature rows in place under a **Columns** summary instead of a sub-dialog, checks every mapped column's kind
before the engine (the API too, which used to check presence only), suggests close names only within a distance
threshold instead of always assigning a column, makes every check a validator with **OK** unavailable until all pass,
skips or imputes rows with missing values at application, checks for a cancel before every batch, records every
application, and names and tags the prediction column with Forge's own `<target> (predicted)` and `forge.model`.

## EMS schema `forge`

Manifest: `databases/forge/schema.json`, version `0.1.6`. `grok publish` deploys it; a debug publish applies
destructive changes without migration scripts, except a `promotion` change on a row table, which is always refused
(`[promotion-change]`). To change it, publish a manifest with only a placeholder table, then the real one, bumping
`version` both times; the first publish drops the tables and their data. `grok api` generates the typed client
`src/generated/db.ts` (`forgeDb.models`, `forgeDb.trainingRuns`, `forgeDb.applications`). Every table also has the
system columns `id`, `version`, `created_on`, `updated_on` and `author_id`; who and when always come from them.

| Table | Security | Why |
|---|---|---|
| `model` | `row`, `promotion: lazy`, `defaultRowVisibility: none`, grants `All users: view, edit` | A model is private until shared, like the old models. Lazy promotion makes a model a platform entity on its first share; sharing, favorites and comments work from then on. The table Edit grant lets any user insert models; with visibility `none` it reveals no foreign rows |
| `training_run` | `row`, lazy promotion, `defaultRowVisibility: none`, same grants | The trainer's experiment history, including failed and cancelled runs. `model_id` is optional with `onDelete: setnull`, so the history outlives the model |
| `application` | `master`, `delegate: model_id`, `audit: false`, no grants | Secured by the model (see Key decisions). Master tables refuse grants. No audit: high-churn records. `onDelete: cascade` from the model. `status` is `completed`, `failed` or `cancelled`; `skipped_rows` counts the rows left without a prediction because of missing values. Inserts pass `source` and `status` explicitly |

Deletes are soft: a deleted row stays with `is_deleted` and disappears from queries and counts. Current EMS limitation:
a soft delete does not clean up a promoted row's entity and permissions, so a model deleted after being shared leaves
a live entity until the platform fixes it. The generic Domains view offers **New model...** to anyone with the
table's Edit grant; such a row has no engine blob, and Forge lists it and can delete it (the row only, unless its
file sits in Forge's own `forge/model/<uuid>/` layout).

### JSON column shapes

EMS `json` columns hold objects only; a list is wrapped in an object. TS types: `training/train-model.ts`,
`preparation/preparation-options.ts` (`PreparationOptions`), `engines/engine.ts` (`Hyperparameters`) and
`storage/dataset-fingerprint.ts`.

| Column | TS type | Shape |
|---|---|---|
| `model.target` (caption "Target details") | `TargetSchema` | `{name, type, semType?, categories?}`; `categories` for classification targets |
| `model.features`, `training_run.features` | `FeaturesSchema` | `{columns: [{name, type, semType?}]}`, in training order |
| `options` | `PreparationOptions` | Preparation replay, see below; Forge writes `{preprocessingInfo, postprocessingInfo: [], missingValues: {mode, neighbors?, distance?, skippedRows}}` |
| `hyperparameters` | `Hyperparameters` | `{<train function input>: value}` |
| `metrics` | `MetricsRecord` | `{train: {<metric id>: number}, validation: {<metric id>: number}, positiveClass?}`; ids `mse`, `rmse`, `mae`, `r2`, `accuracy`, `f1`, `sensitivity`, `specificity`, `precision`, `npv` (`auc` reserved) |
| `splitting` | `Splitting` | `{scheme: none \| kfold \| holdout, folds?, trainFraction?, isStratified?}` |
| `model.dataset_ref` | - | Reference mode: `{kind: file \| query \| script, id?, path?, params?, script?}` |
| `dataset_fingerprint` | `DatasetFingerprint` | `{rowCount, columnCount, hash, columns: [{name, type, missingCount, min?, max?, mean?, categories?}]}` |

`options` keeps the keys of the platform's built-in models, so their preparation replays unchanged:
`preprocessingInfo`, `postprocessingInfo`, `positiveClass`, `negativeClass`, `binaryClassificationThreshold`,
`targetType`, `allowNulls` (`allowNulls` is read and ignored). `PreparationOptions` and `preparationOptionsOf`, which
narrows a stored value to it, are in `preparation/preparation-options.ts`.

### Model blob

The `model.blob` column (type `file`) stores `file://System:DomainFiles/forge/model/<uuid>/model.bin`, the layout the
EMS row editor uses for file columns. The random path segment makes the file unguessable, and the row security of
`model` protects the pointer. The generic Domains UI shows the column but offers no download of the file.

### Storage modes

`model.storage_mode` (default `none`) records what a model keeps of its training data:

- `none`: only `dataset_fingerprint`, enough to check that a new table fits the model. The only mode used so far.
- `reference`: `dataset_ref` points to the source (file, query or script); the data is not uploaded.
- `copy`: `dataset_table_id` is the id of an uploaded copy of the training table (no foreign key; not written yet).

`has_training_rows` marks blobs that embed training rows (SVM, KNN); `false` for XGBoost. `legacy_id` is reserved for
migrated models.

### Captions

Every column of the three tables declares a `friendlyName`, because the platform renders undeclared ones
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

`application`: `model_id` Model, `table_name` Table, `row_count` Rows, `skipped_rows` Skipped rows, `column_name`
Prediction column, `source` Source, `status` Status, `error` Error, `duration_ms` Duration (ms).

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
| `apply` | `(df, model) -> dataframe` | Forge takes the first column as the prediction, out of the returned frame |
| `isApplicable` | `(df, predictColumn) -> bool` | |
| `isInteractive` | `(df, predictColumn) -> bool` | Optional. `meta.mlupdate: 'false'` turns off live retraining |
| `visualize` | `(df, targetColumn, predictColumn, model)` | Optional |

- An engine is complete, and usable, only with `train`, `apply` and `isApplicable`.
- A function engine receives the feature table and the target column. A script engine (`DG.Script`) receives one
  table with a copy of the target appended, and the target name.
- `train` returns the blob as a `Uint8Array` or a `DG.FileInfo` (its `data` is used); anything else is a
  `ForgeError`. `apply` must return a dataframe with at least one column.
- Hyperparameter defaults come from the `train` annotation (`//input: int iterations = 20`). The platform keeps them
  as text in `Property.initialValue` (`Property.defaultValue` stays null), a string one in quotes (`'RBF'`);
  `defaultValuesOf` (used by `defaultHyperparameters` and for the imputation inputs) converts them by the input type
  and strips the quotes.
- `EngineRegistry.discover()` calls `DG.Func.find({meta: {mlrole}})` once per role, groups by `mlname`, takes kind and
  namespace from the `train` function, and lists function engines first, then scripts, each by name.
  If a package is registered twice, the later function of a role wins.
- The namespace is the package name for function engines (EDA registers as `Eda`) and the `nqName` prefix for scripts.
- `isApplicable` on an engine without that role throws `ForgeError`; `isInteractive` without it returns `false`.
- Discovery reads the client function registry synchronously; script engines that are not loaded into it yet are
  missed. No script engines exist yet to verify this.
