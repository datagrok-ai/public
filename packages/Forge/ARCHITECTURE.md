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
- **Training data never leaves the browser without consent.** Only the `copy` storage mode, chosen in the Save
  dialog, uploads it; `none` keeps a fingerprint of the data and `reference` a link to its source.
- **Every training attempt is recorded.** Completed, failed and cancelled runs are `training_run` rows, from the Train
  view (`source: ui`) and from `Forge:trainModel` (`source: api`, saved or not); a model row appears only on Save. A
  training superseded by a newer one (`TrainingQueue`) is not an attempt the user finished and is not recorded.
- **Every application is recorded.** Completed, failed and cancelled applications are `application` rows; a mapping
  the model cannot use is refused before any work and not recorded.
- **Application records inherit the model's security.** A record is visible to whoever can see its model. It can be
  written by anyone the `model` table lets insert, which today is every user through the table's Edit grant, so
  applications by users who can only view a model are recorded too.
- **No copies of the user's data.** Training and application pass lists of the user's columns themselves; a column
  is copied only where Forge would otherwise change it (see "Missing values" and "Application pipeline"), and the
  preparation steps add only Forge's own columns. A function that takes a table gets a frame of the columns
  (`preparation/shared-frame.ts`: `sharedFrame`) for the length of the call: every engine call
  (`engines/engine-calls.ts`), the imputer (`preparation/missing-values.ts`) and the training-copy upload
  (`storage/dataset-copy.ts`). A frame stays a parent of its columns and receives their events, so each runs its call
  through `onFrame`, which gives the frame back with `releaseFrame` (its columns removed) when the call ends, whatever
  happens.
- **Missing values reach only a method that handles them.** Rows with gaps are skipped or the gaps are imputed before
  the engine is called, unless the method has `meta.mlnativenans` (see "Missing values"); dates and very large whole
  numbers (bigint) are not numerical features.
- **Later:** more engines, MLflow models and the migration of existing models.

## Components and data flow

```
package.ts ──> ui/forge-app.ts ───> engines/                     (DG.Func registry, methodsFrame)
           │                    ├─> storage/, generated/db.ts    (catalog query, modelsChanged)
           │                    ├─> catalog/                     (Applicable to, Compare)
           │                    └─> ui/model-actions.ts, ui/model-handler.ts, ui/model-comparison.ts,
           │                        ui/apply-model-dialog.ts, ui/data-grid.ts
           ├─> ui/train-view.ts ──> training/ ──> preparation/, engines/, metrics/, apply/ (evaluation.ts)
           │                    ├─> ui/results-pane.ts ──> training/insights.ts, ui/metrics-table.ts, ui/results-layout.ts,
           │                    │                      │   the method's visualize (preparation/shared-frame.ts)
           │                    │                      └─> ui/results-viewers.ts ──> metrics/, preparation/shared-frame.ts,
           │                    │                                                    apply/ (scatteredColumn, copyInto), ui/data-grid.ts,
           │                    │                                                    ui/results-layout.ts (classesOf)
           │                    ├─> ui/runs-panel.ts ──> metrics/, ui/data-grid.ts, ui/results-viewers.ts (resultsTitle),
           │                    │                       engines/, preparation/, training/, generated/db.ts (types)
           │                    ├─> ui/save-model-dialog.ts, ui/apply-model-dialog.ts, ui/forge-app.ts (the balloon)
           │                    └─> storage/  ──> generated/db.ts  (grok.dapi.domains -> EMS schema forge)
           │                                  ├─> grok.dapi.files  (System:DomainFiles/forge/model)
           │                                  └─> grok.dapi.tables (training data copies)
           ├─> ui/apply-model-dialog.ts ──> apply/, engines/, preparation/, generated/db.ts (model query),
           │                                constants.ts (MAX_MODELS)
           ├─> apply/ (runApplyModel) ──> preparation/, engines/, storage/ (findModelRow, blob path),
           │                              generated/db.ts, grok.dapi.files (blob read)
           ├─> training/train-api.ts (runTrainModel) ──> training/, preparation/, engines/,
           │                                             storage/ (saveTrainedModel, the run record)
           ├─> storage/model-store.ts (runDeleteModel: findModelRow, then deleteModelRow)
           ├─> training/train-api.ts (runListMethods) ──> engines/method-list.ts (methodList), training/
           ├─> catalog/list-models.ts (runListModels) ──> generated/db.ts, catalog/ (applicableTables, metricColumns),
           │                                              apply/ (featureSchemasOf), grok.dapi (entities)
           ├─> catalog/get-model.ts (runGetModel) ──> storage/ (findModelRow, ownBlob, dataset-ref), catalog/
           │                                          (applicableTables, loginsOf, model-activity), apply/, metrics/,
           │                                          training/
           ├─> catalog/model-edit.ts (runUpdateModel) ──> storage/ (findModelRow, modelsChanged), generated/db.ts
           ├─> training/evaluate-api.ts (runEvaluateModel) ──> training/ (evaluation.ts, insights.ts, metricsRecordOf),
           │                                                  apply/ (loadedModelOf, apiMapping), engines/,
           │                                                  storage/ (findModelRow, dataset-ref)
           ├─> ui/model-handler.ts ──> ui/rich-tooltip.ts (the tooltip), generated/db.ts (features read once),
           │                       apply/ (featureSchemasOf), catalog/ (featureFit),
           │                       storage/model-fields.ts (tagsOf), storage/model-store.ts (modelsChanged, followed),
           │                       ui/save-model-dialog.ts (STORAGE_CAPTIONS)
           │                       └─> ui/model-panes.ts (modelAccordion, modelIcon, userOf, isShownVersion) ──>
           │                                                  apply/application-store.ts (applicationRecorded),
           │                                                  catalog/, generated/db.ts, grok.dapi (users,
           │                                                  entities, permissions), DG.DomainObjectHandler
           │                                                  (audit, sharing), ui/data-grid.ts,
           │                                                  ui/save-model-dialog.ts (writeModelInfo),
           │                                                  ui/metrics-table.ts (MetricsTable, storedSummary),
           │                                                  ui/evaluation-view.ts ──> training/evaluation.ts
           │                                                  (storedTrainingTable), apply/, catalog/ (applicableTables),
           │                                                  engines/, ui/metrics-table.ts, ui/results-pane.ts,
           │                                                  ui/button-gate.ts
           ├─> ui/model-actions.ts ──> catalog/, storage/, ui/model-handler.ts, ui/apply-model-dialog.ts,
           │                           ui/edit-model-dialog.ts, DG.DomainObjectHandler (openShareDialog: shareRow)
           ├─> ui/model-comparison.ts ──> catalog/, ui/data-grid.ts, ui/model-panes.ts (modelIcon),
           │                              PowerGrid's Forms viewer
           ├─> ui/prediction-column-panel.ts ──> ui/model-handler.ts (the card)
           └─> ui/table-models-panel.ts ──> catalog/ (modelsFor, model-activity), constants.ts (MAX_MODELS), generated/db.ts,
                                            storage/model-fields.ts (tagsOf), ui/model-handler.ts,
                                            ui/model-panes.ts (modelFacts), ui/apply-model-dialog.ts (openApplyDialog),
                                            ui/train-view.ts
```

- `package.ts` holds the annotated functions only: the app `forgeApp` (shown as "Forge", `Browse > Apps`),
  `forgeModels` (`ML | Forge | Models`), `forgeTrain` (`ML | Forge | Train...`), `forgeApply` (`ML | Forge |
  Apply...`), the API functions `applyModel`, `trainModel`, `deleteModel`, `listMethods`, `listModels`, `getModel`,
  `updateModel`, `evaluateModel` and `shareModel`, the autostart `_initForge` (registers
  the model handler with the model commands and the context panel's follower, and the comparison handler),
  `isPredictionColumn`, the column panel
  **Predicted by** (function `predictedByPanel`) and the table panel **Forge models** (function `tableModelsPanel`).
  They delegate to `ForgeApp`, `TrainView`, `openApplyDialog`, `runApplyModel`, `runTrainModel`, `runDeleteModel`,
  `runListMethods`, `runListModels`, `runGetModel` (with `grok.shell.tables`), `runUpdateModel`, `runEvaluateModel`,
  `openShareDialog`, `ForgeModelHandler.registerOnce`, `registerModelActions`, `ForgeModelHandler.followChanges`,
  `ModelComparisonHandler.registerOnce`, `isForgePrediction`, `predictedByPane` and `tableModelsPane`. `_initForge` reads the registered handlers' names once
  (`DG.ObjectHandler.list()`) and both `registerOnce` calls skip a handler already there, so a second `_initForge` (or
  the test bundle) registers nothing twice; the model commands and the follower are set up only with a newly
  registered handler.
- Every table of data in the UI is `readOnlyGrid` (`ui/data-grid.ts`): `DG.Viewer.grid` without editing, row header,
  current-row indicator and new-row icon, as wide as its host and as tall as its rows (up to 15 by default,
  then it scrolls: an inline height of `(rows + 1) * rowHeight + 2` px; the app's Methods grid shows every method), with tooltips from `{headers, cells}` (by column name; the cells
  as functions of the row) through `onCellTooltip` (`gridTooltip`; a cell without its own text shows its value).
  `numberColumn` builds the grids' float columns (format `0.000`); `textColumn` and `boolColumn` (`columns.ts`, not UI,
  also used by `catalog/`, `engines/method-list.ts` and `training/train-api.ts` for the API lists) their string and bool columns.
- `ForgeApp` lists the methods (section **Methods**, the grid `methodsGrid` over `methodsFrame` of every discovered
  method, `engines/method-list.ts`, shared with `Forge:listMethods`: Method, Package, Method type, Roles,
  Hyperparameters; the Method type cell's tooltip says what a function or a script method is, the Roles cell's gives
  each role's meaning, the Hyperparameters cell's each input's description) and the catalog. Its `aiDescription` (the
  view's briefing for the platform's AI assistant) describes both sections, **Applicable to**, the context panel and
  the icons, and names the API functions and `agents/skills.md` (see "Grokky knowledge (agents/)").
  `ForgeApp.loadModels()` queries the catalog columns and sets the captions; the refresh icon reloads the grid. The
  play icon (**Apply model**) and the trash icon (**Delete model**) are disabled while the catalog has no current row
  (re-evaluated on current-row change and on every reload; the subscriptions to the current frame are replaced on
  reload and dropped in `detach()`). The play icon opens the Apply dialog for the current row's model through
  `openModelApply` (`ui/model-actions.ts`, shared with the model's **Apply...**): the current table (the catalog is not
  a table view, so in practice the first open table; without tables, "Open a table first.") and **Applicable to** as
  `preferredTable` (see "Apply dialog"), and that table's view after applying. The trash icon opens
  `confirmDeleteModel` (`ui/model-actions.ts`), whose OK calls `deleteModel`. The
  Compare icon (`compareSelected`, tooltip **Compare in a new view**) is enabled while two or more rows are selected
  and opens `openComparisonView` (`ui/model-comparison.ts`) of their rows, read with `COMPARE_COLUMNS` and kept in
  catalog order: the table view **Compare models** of `compareModels` with the PowerGrid **Forms** viewer added
  (without PowerGrid, a balloon says to install it). **Applicable to** is a `ui.input.table` (empty = every model; the
  platform's input tracks the open tables and keeps its Open file icon) and filters the grid's rows in the browser with
  `requiredFeaturesOf` and `isSuggested`, as `applicableTables` does, once per distinct required feature list in a pass; the input empties itself when its table is closed but
  reports no change, so the catalog filters again on `onTableRemoved` when the table it filters by closed, and treats
  a closed table as none. A new current row becomes the current object
  (`grok.shell.setCurrentObject(ForgeModelHandler.rowOf(...), true, true)` of the row's values, the `features` text
  parsed; forced, since a plain set within a second of the previous one is ignored); no current row, two or more
  selected rows, or the model the panel already shows, leave the context panel as it is. A selection change, debounced
  300 ms, makes a `ModelComparison` of two or more selected rows the current object (read with `COMPARE_COLUMNS`; a
  sequence number drops a slower, older read); with fewer it shows the current row's model, without a current row the
  one selected row's, and with neither it empties a comparison (`setCurrentObject(null)`, the platform's empty panel,
  `property_panel.dart:116-118`); `detach()` outdates a read still running. A reload (the refresh icon or
  `modelsChanged`) replaces the grid's frame, so `refresh()` reads the current row and the selected rows when the new
  frame arrives and finds them again by id (a sequence number drops an older reload that ends later); a deleted row is
  simply gone. The panel keeps the object it shows, unless that is a comparison with a model the reload no longer has
  (`isStaleComparison`, a delete): then it shows the selection again by the rule above (once: a restored selection does it
  through its change event), and is emptied rather than left on the comparison. A single model the panel shows follows
  its changes and its delete by itself, catalog or not (`followChanges`, see "Model handler and panels"). Right-clicking a row adds `CATALOG_ACTIONS` (`addModelItems`, with
  **Applicable to** as the preferred table) and, with two or more selected rows, **Compare** to the grid's menu: the
  grid's menu is its own, so the platform's model commands are not in it. The header holds the **Models** title and
  the icons on the title's baseline (`forge-pane-header` in `css/forge.css`), the icons in `--blue-1`, the blue of an
  input's own icons (ui.css), and in the platform's `d4-disabled` grey while disabled; **Applicable to** is on the
  line under it. The grid takes the full width of the view: an inline `width: 100%`, the override the platform's
  400px rule for boxes in a panel (`ui.css`) asks for. The grid's `onCellTooltip` (through `modelTooltip(cell)`) shows the
  model's rich tooltip (`forgeModelHandler.renderTooltip` of the row's catalog values, as the Domains card and the
  panel title) for a **Name** cell; any other cell gets null and keeps the grid's own tooltip.
- `modelsChanged` (an rxjs `Subject` in `storage/model-store.ts`) fires after every successful `saveModel`,
  `updateModelInfo` and `deleteModel`, whoever calls them; every open `ForgeApp` subscribes in `this.subs` and reloads
  itself, 300 ms after the last of a burst (tag chips written one by one reload it once); the context panel follows it
  too (`ForgeModelHandler.followChanges`, subscribed once by `_initForge`).
- `TrainView` (tab **Predictive model** with the core's model icon, `ui.iconSvg('model')` returned by `getIcon()`)
  holds the form, the live validation, the live training and the results. `TrainView.openFor(table)` opens it on a
  table (the boundary of the table pane's **Train**); `TrainView.open()` (`ML | Forge | Train...`) calls it with the
  current table, or says "Open a table first.". The form is one `ui.form`, so the platform
  aligns every label in it: the group **Data** (Table, Target, Features), the group **Preparation**
  (`PreparationInputs`, below; hidden while none of its inputs applies), the group **Method** (Method and the
  hyperparameters), all expanded at the start, and the **Train** row
  (`ui.buttonsInput`, `trainButton`). The view is a resizable split (`ui.splitH(..., true)`): the form's panel starts 400 px wide
  (`FORM_WIDTH`) and scrolls by itself; the results take the rest (`forge-train-results`, scrolling), with the
  training loader (`resultsStatus`) in the header of the pane's **Performance** panel.
  The `trainRow` is right-aligned (`forge-train-row`). The ribbon holds only **Save** (`ui.bigButton`). The platform's column inputs keep a column removed from the
  table until their dialog is opened again, so the form listens to the table's `onColumnsRemoved` (in `formSubs`) and
  empties **Target** or unchecks the **Features** that are gone (`dropRemovedColumns`), each a change as any other. A
  Table change rebuilds the form; the groups'
  subscriptions (`formSubs`) are dropped then and in `detach()`, the hyperparameter inputs' (`methodSubs`) on every
  Method change too. The view keeps the latest training in `lastTraining`: `result`, `runId`, the method, the table and the user's feature and target
  columns, the `datasetName` and `fingerprint` taken at training time, and `isSaved`; any change drops it at once.
  **Save** is disabled until a training result exists, while a training runs, and once that result is saved: OK in
  `ui/save-model-dialog.ts` (Name, prefilled, Description, Tags and **Data storage**, see "Saving") calls
  `saveModelAs(info, choice)`, which marks the training saved before it writes (`saveTrainedModel`, see "Saving"), so
  one training becomes at most one model; a failed write enables **Save** again and shows the write's error. A
  failure of the run's link comes after the row is written, so it also enables **Save** again and a second Save
  makes a second model of the same training. **Name** is not nullable and its validator refuses a blank name ("Enter
  a name for the model."), so OK is disabled with that tooltip and the dialog's own validation blocks Enter too. The
  same form (`modelInfoDialog`; `saveModelDialog` adds the storage block to it) is **Edit model**. After a save the balloon `Model "<name>"
  saved.` has the links **Apply...** (`openApplyDialog` on the training table, the new model preset) and **Show in the
  catalog** (`ForgeApp.open()`, which focuses an open catalog view instead of opening a second one). The results are
  a `ResultsPane` (`resultsPane`, see "Results"); `Training.request` keeps the `TrainingRequest` the
  result was trained on (the prepared columns, the selection's target, `keptRows`) for it. The Performance panel's grid is a
  `MetricsTable` of `ui/metrics-table.ts`, shared
  with the model's **Performance** pane (it builds one there). The pane keeps one while results come (a training,
  a cutoff re-cut): `update(summary)` rewrites the values in place, refills the bullet list and repaints the marks when
  the metric rows and the value columns are the same, so the grid neither flickers nor loses its column widths; other
  rows (task, AUC-ROC) build a new one; a hint shown in between detaches the grid and keeps it. It is a frame of **Metric**, **Train**, **Validation**
  and, with `MetricsSummary.evaluation`, **Evaluation** (format `0.000`) in a
  `readOnlyGrid` with the header texts, the metric descriptions and the full-precision values as tooltips, then a
  bullet list (`ul`): `Missing values: kept (the method handles them)` (only under `keep` when the training features
  had gaps), `Rows: <used> used, <skipped> skipped (missing values)` (only when rows were skipped),
  `Validation: <folds>-fold cross-validation on <rows> rows`, `Seed: <seed>` with the seed in a selectable
  `forge-seed` span and the copy icon right after it (**Copy the seed**, balloon `Seed <seed> copied.`),
  `Positive class: <class>` (two classes only) and `Evaluation: <table> (<rows> rows[, <n> skipped (missing values)])`
  (with an evaluation). The `metricFindings` of the metrics mark their cells: `onCellPrepare` paints a marked value
  (the finding's `splits`) in `DG.Color.failure`, and its tooltip ends with the finding's text (`R2 is below 0.3.`).
- `CollapsibleGroup` (`ui/collapsible-group.ts`, classes `forge-group*` in `css/forge.css`) is the folding block of
  the Train view groups and of the Apply dialog's **Columns** and **More options**: a header with a chevron, the
  caption and an optional summary (red with `forge-group-invalid`), in Diff Studio's style; the chevron sits in a
  fixed-width box (`forge-group-chevron`), so captions line up whether a group is open or folded.
  `expandOnError(inputs)` opens the group when one of its inputs is validated with an error; the state is not
  remembered. The header carries `aria-expanded` and the summary `aria-invalid` (see "Automation surface").
- Table, Target, Features and Method carry short Forge tooltips (`tooltipText`, the caption's tooltip). Hyperparameter
  inputs come from `ui.input.forProperty`, whose caption tooltip is the `train` input's own `description`. For an
  invalid input the platform shows the tooltip text with the validator messages in red below it.
- **Missing values** (`ui/missing-values-inputs.ts`, `MissingValuesInputs`, used by the Train view and the Apply
  dialog) is a radio of **Skip rows** / **Impute** with its options on one line (`forge-inline-radio`), shown only
  when its columns have gaps, with the tooltip "Rows with missing values in: <col>: N, ...". **Impute** is offered
  only when `Eda:knnImpute` exists; choosing it shows the function's own inputs (**Neighbors**, **Distance**, from
  `imputeSettingsOf` with `defaultValuesOf`) right under it. Hiding them puts their defaults back, so a hidden invalid
  value cannot block the form. `setKeeps(true)` (the Train view, for a method with `handlesMissingValues`) and
  `restore` of a `keep` setting (the Apply dialog) hide it all and make `settings()` `keep`, with or without gaps.
- **Preparation** (`ui/preparation-inputs.ts`, `PreparationInputs`, Train view only) owns its `CollapsibleGroup`:
  **Missing values** (the `MissingValuesInputs` above), **One-hot encoding** (shown while a checked feature
  is text or yes/no; checked when it appears if every such feature, the all-unique ones aside while Skip unique
  categories is checked, has at most `ONE_HOT_MAX_CATEGORIES` (20) categories (`defaultOneHot` in
  `preparation/pipeline.ts`, also `Forge:trainModel`'s default), else unchecked: a deviation from the
  built-in tool, which always starts unchecked; it follows that default until the user toggles it, `isOneHotChosen`,
  and hiding resets it), **Skip unique categories** (checked; shown while a checked feature has `hasUniqueCategories`),
  **Predict probability** (unchecked; shown while `twoClasses(target)` is not null, its tooltip naming the positive
  class) and under it **Positive class cutoff** (`ui.input.float` 0..1, step 0.01, with a slider, `DEFAULT_CUTOFF`; shown while
  Predict probability is checked; a validator and `cutoffProblem` refuse an empty value or one outside 0..1). The show conditions are the pipeline's own predicates (`preparation/pipeline.ts`).
  `update(features, target)` runs on every Features and Target change; a hidden input gets its default back, so it
  reappears as it first appeared, and the group is hidden while nothing in it is shown. `steps()` is the
  `PreparationSteps` of the selection; a hidden input's step is off (Skip unique categories, checked while hidden, is
  sent only while it is shown).
- UI entry points and handlers are the only error boundaries (`ui/report-error.ts`): `ForgeError` becomes a warning,
  anything else an error balloon plus `_package.logger.error`. API functions do not catch: the error reaches the
  caller.
- Logic folders (`engines/`, `preparation/`, `training/`, `metrics/`, `storage/`, `apply/`, `catalog/`) work without
  the DOM, throw instead of catching, and never import `ui/`. The same logic serves the UI, API functions and tests. The
  catches in them are `applyAndRecord` and `runTrainModel`, which record a failed or cancelled application or training
  and rethrow, `saveTrainedModel`, which deletes the copy of a failed save and rethrows, `runEvaluateModel`, which
  turns an error after a cancel into `Evaluation was cancelled.`, and `TrainingQueue`, whose answer is the outcome
  with the error in it. `catalog/` (see
  "Catalog logic") builds on `apply/`, `storage/`, `metrics/` and `training/`.

## Training pipeline

`training/train-model.ts`. A `TrainingSelection` is what the user chose: the engine, the checked feature columns
themselves (a list, no frame), the target column, the hyperparameters, a seed, the number of folds (5), the
missing-value handling and the preparation steps (`steps`, a `PreparationSteps`, see "Preparation"; the Train view
takes them from `PreparationInputs.steps()`). `prepareTraining(selection)`
handles the missing values (see "Missing values"), then applies the steps (`prepareFeatures`), and returns the
`TrainingRequest`: `features` and `target`, the selection's columns with the missing values handled (the user's own
columns where nothing changed; the model's feature list, target schema and data summary describe them), `keptRows`
(the selection's rows they hold, a `BitSet` over the table, or `null` when no row was skipped), and
`prepared` (`PreparedFeatures`), what the method is trained on, with the `options` that record the whole preparation
(`options.missingValues.skippedRows` is the one record of the skipped rows). No frame is built: the engine calls
build their own. After preparation fewer than `2 * folds` rows is a `ForgeError` ("After skipping rows with missing
values, N rows remain; training needs at least 10.").

**Missing values in the view.** **Missing values** (see Components) is the first input of **Preparation** and follows
the checked features; `fillMethod`, which every method change goes through, hides it for a method that handles
missing values (`PreparationInputs.setKeepsMissingValues`), and every training of it records `keep`, with or without
gaps. The **Target** tooltip adds "<target>: N rows without a value are skipped." ("1 row without a value is
skipped.") when the target has gaps; it is information, not a problem. After a training with skipped rows,
**Performance** lists `Rows: N used, M skipped (missing values)` first in the bullet list under the metrics grid;
after a training under `keep` whose features had gaps, `Missing values: kept (the method handles them)` comes before
it (`keepsMissingValues` of the summary: the recorded mode and `featuresHadGaps` of the training's fingerprint, the
feature columns' `missingCount`; a saved model's Performance pane, `storedSummary`, reads the same from
`options.missingValues` and `dataset_fingerprint`). The Runs cell says `keep missing values` likewise only when the
features had gaps (`SessionRun.featureGaps`, from `missingColumnsOf` of the selection: a run that failed before
its data was prepared has no fingerprint).

**Defaults in the view.** Target: the table's last column. Features (`defaultFeatures`): every column the methods
read as numbers (`isReadableNumber`: numerical, not dates, not bigint) except the target and integer columns without
missing values whose values are all different
(row numbers and ids, such as the `col 1` row-number column of `iris.csv`, which is sorted by species). Float columns
are never excluded. Hyperparameters: `defaultHyperparameters(engine)`, read from the `train` inputs' initial values.
Changing the target leaves the features as they are.

**Problems.** `checkSelection` (below) returns the messages per input (`target`, `features`, `missingValues`,
`method`). "Numerical" means `numerical_no_datetime`: dates are not numbers for Forge. Rules 1-7 and 9 are cheap
checks on the column list (rule 7 on the columns `prepareFeatures` gives: the selection's preparation steps applied,
as a method would get them); rule 8 calls the engines on those prepared columns, and runs only when
1-7 pass (`hasDataProblems` is false):

| # | Condition | Input |
|---|---|---|
| 1 | no feature | Features |
| 1b | Skip unique categories leaves no prepared feature: `Skip unique categories leaves no feature. Check more features.` | Features |
| 2 | the target is also a feature | Features |
| 3 | fewer rows with a target value than `2 * folds` (10) | Target |
| 4 | a bigint target ("holds very large whole numbers", `bigIntProblem` in `training/default-features.ts`) | Target |
| 5 | a classification target with fewer than two values | Target |
| 6 | the target is neither numerical, text nor boolean (dates included) | Target |
| 7 | a prepared feature is text or yes/no, which happens only while One-hot encoding is off: `<Method> needs numerical features. Check One-hot encoding in Preparation, or uncheck: SEX, RACE.`; any other prepared feature that is not numerical (dates included): `<Method> needs numerical features. Uncheck: <names>.` | Features |
| 7b | a bigint feature, one message per column | Features |
| 8 | the engine's `isApplicable` says no: `<Method> cannot learn from this selection. It needs numerical features and a numerical, text or boolean target.` | Method |
| 9 | Impute chosen but not available, a yes/no feature with gaps, or fewer than two features to impute from | Missing values |

Rows with a missing target are skipped by `prepareTraining`, never refused.

**Choosing the method.** `checkSelection(selection, engines)` is the check for a view that offers the methods of
`engines` (the view passes its own list, to keep its `Engine` objects). It runs rules 1-7 and 9; while a target or
features rule fails it asks no method and returns no engines (`failedCheck(problems)`, the shape of every check that
lists no method). Otherwise it asks every
complete method at once (`applicableEngines`, see "Choosing a method" under the engine contract) about the prepared
columns and target (so with Predict probability the regressors are listed, as for any numerical target), and rule 8
becomes: the selection's method is not among them (the rule 8 text), or none is (`No method can learn from this
selection. Check the features and the target.`). The result (`SelectionCheck`) also has `best` (`selectBestEngine`
over the listed methods and the prepared columns), `failed` (the methods whose check threw, for the caller to log),
`isInteractive`: the selection's method is listed and `retrainsLive` (its `isLiveUpdate` is on and its
`isInteractive` says yes on the prepared columns), `prepared`, those columns and target (none from
`failedCheck`), and `insights`, the data insights of the prepared columns and the selection's target (see "Insights";
empty from `failedCheck`). A method's `isInteractive` that throws makes the whole check throw.

The view adds a validator to **Target**, **Features**, **Method** (the `method` messages) and **Missing values** that
returns the view's current problems for that input. Every change of Table, Target, Features, a **Preparation** input
but the cutoff, Method or a hyperparameter (`requestCheck`) numbers a new check, clears the problems, drops
`lastTraining`, supersedes the running training (`TrainingQueue.supersede()`), disables **Train** (a tooltip
already shown on it then reads "Checking the selection...") and schedules the check `CHECK_DELAY_MS` (200 ms, the
built-in tool's delay) after the last change; one timer serves the check and the live training. A check that finishes
after a newer change is discarded. The check (`revalidate`) calls `checkSelection` with the view's own list of
complete methods and logs a method whose check threw once per view (`_package.logger.error`). Then **Method**
(`methodOf`): while a target or features rule fails, or no method is listed, the method stays; otherwise the user's
choice (`userChoice`, set when the user picks a method) while it is listed, else the suggested one (`best`), with the
balloon `<Method> cannot be used with this selection; <suggested> is chosen.` when it replaces the user's choice, which
is then forgotten; for a changed method, which is listed, only `retrainsLive` is asked again, on the check's
`prepared` columns. **Method** then lists the methods by name (empty with the
`No method can learn ...` problem); a failing target or features rule leaves the list as it is. Hyperparameter inputs
are rebuilt on a Method change from the defaults overlaid with `hyperparameterValues`, the values per method name
kept for the view's session (filled on every hyperparameter change; a Table change keeps them). After the check the
inputs, the hyperparameters' own range validators included, are validated again (red input, messages in its tooltip).
An invalid **Neighbors** or hyperparameter (`<caption>: <message>`, e.g. `Iterations: Value must be less than 100`) is
a problem like the check's own: nothing trains, live or on **Train**, and nothing is recorded. **Train** is shown only
for a method known not to be interactive on this table (`interactivity`, kept per method from the last check that
could ask, that is one without a data problem; cleared on a Table change; an unknown method hides it), never for an
interactive one, a problem or not. For a shown **Train**, the tooltip is `Train the model on the current selection.`
while enabled, and it stays
disabled with the first problem as its tooltip, while a training runs ("Training is in progress.") and after a
completed or failed training of the selection until an input changes ("The model is trained on this selection. Change
an input to train again."). A `ButtonGate` (`ui/button-gate.ts`, shared with the Save and Apply dialogs' OK) builds
one tooltip overlay per disabled period (`ui.setDisabled` adds an overlay and a timer per call, both living until
the button is enabled or leaves the document) whose content is a reusable element holding the reason: a new reason
only changes its text (a function would show as its source), and a new overlay is built only when the last one is
gone. On enabling it binds the button's own tooltip again (the enabled tooltip, or none), so
the last disabled reason does not stay on it. A check that fails shows **Train** disabled with
the error as its tooltip and shows each distinct error once. Platform validators do not run on empty values, so an empty
Features shows the platform's "Can't be empty". A rejected selection writes nothing.

**Live training.** A check without problems starts a training at once when `isInteractive` is true (`startTraining`
live); otherwise **Performance** shows `<Method> takes a while on <N> rows, so it does not retrain on every change. Click
Train.` (N, the target's length). **Train** (`train()`) starts one when the selection has no problem and nothing
trains. Every training runs through the view's `TrainingQueue` with a task-bar indicator `Training <Method> model`,
created when the training starts (a request superseded while waiting shows none; cancellable: the indicator's cancel
stops it before its next fit); a loader without text
stands next to the **Performance** header, over the previous results (live) or an empty panel (**Train**). The
queued function prepares the data, builds the run record from the request and trains. Before any training, and after a problem, a failure or a cancel, **Performance** reads `Choose the
target and the features; the model trains as you change them.` (the not-interactive hint after a cancel of such a
method); when the only problems are invalid settings (hyperparameters, missing-values settings or the **Positive class
cutoff**), it reads `Fix the settings.`. `isTraining` is true from the start of a training until its run is recorded; **Save** waits for it. A result
that arrives after a newer change is recorded but not shown.

**Task.** A numerical target (not a date) is regression; text and boolean targets are classification. The task
follows the target the method learns: with Predict probability that is the float 0/1 target, so the task is
regression, while the target schema stays the selection's (with its two categories).

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

**Training queue.** `TrainingQueue` (`training/training-queue.ts`) runs one training at a time for a view that
retrains as the user changes the form. `run(train, indicator?)` waits until the trainings before it have settled,
then, unless superseded meanwhile, creates the indicator (`indicator` is a factory, so a request that never starts
has none) and calls `train` with a `LoopProgress` whose `canceled` is true once the request is superseded or the
indicator is cancelled, and whose `update` goes to the indicator. A newer `run` supersedes the running training,
which stops before its next fit (`trainModel`'s check), and every waiting one, which never starts; `supersede()` does
the same without a newer request (the view's change that may not lead to a training). It resolves, never rejects,
with a `QueuedTraining`, a union on `outcome`: `completed` with `result`; `superseded` whatever the training did;
`cancelled` (the indicator was cancelled) or `failed`, each with `error`. A superseded training is not recorded; the
other outcomes are the run statuses below.

**Run records.** `TrainView` writes the `training_run` row (`source: ui`) right after the queue settles: `completed`
with the metrics, `cancelled` (the yellow `Training was cancelled.`), `failed` with the error message; a `superseded`
training writes nothing, and neither does a failure before the data was prepared (too few rows left). `model_id` is set by
`linkTrainingRun` when the model is saved; deleting the model keeps the run (`setnull`).

**Runs.** Each settled training that is not superseded also becomes a `SessionRun` (`ui/runs-panel.ts`) put first in
`TrainView.runs`: its number in the view, the status, the method (the `Engine`), the target and feature names, the hyperparameters, the missing-values settings, the preparation steps, the
record's row count (0 when nothing was recorded), the seed, and the metrics or the error. A training of a table the view has left adds nothing
(`createForm` replaces the list, the training compares it). `RunsPanel` shows them under the results
(`ui.panel([resultsPane.root, runsPanel.root], 'forge-train-results')`; the panel's own root is `forge-runs`, 12 px under the rows): the core's title `Runs` with its tooltip
(`PANEL_TEXTS.runs`, `resultsTitle`), hidden while the list is empty (`detach()` empties it, which drops the grid's
click subscription), and a `readOnlyGrid` with the current-row
indicator on, built again on every `update` over a copy of the list (a new frame, a new grid of the new height; the
replaced grid goes with `ui.empty`, the platform's widget disposal): **#**, **Method**, **Target**,
**Features** (the count), **Preparation** (`impute (<n> neighbors, <distance>)`, `one-hot`, `skip unique
categories`, `probability (cutoff <c>)`, joined by commas), **Rows**, **Status**, **Seed**, then `<metric>
(validation)` for every metric of `METRIC_IDS` any run's validation has, with header tooltips (the metrics'
`METRIC_DESCRIPTIONS`) and cell tooltips (the hyperparameters, the feature names, the error or `Stopped from the task
bar.`, `Train <0.000>, validation <0.000>` with `-` for a missing value). The newest row is current; a click on a
row (`grid.onCellClick` on a table cell, the current row too) calls `restoreRun`, while code that moves the current
row restores nothing (the read-only grid moves nothing with the arrow keys): a run whose target or a feature the
table no longer has warns `The table no longer has the columns of this run: <names>.` and sets nothing; a run whose
inputs equal the form's (`isSameInputs` over `runInputsOf` of the current selection: method, target, features,
hyperparameters, missing values and preparation steps but the cutoff, the same values a run keeps) while the view
shows the completed result of those inputs (`lastTraining`, which any change but the cutoff drops) trains nothing:
it only sets the run's cutoff when that differs, which re-cuts the shown result; after a cancelled or failed
training there is no such result, so the click sets the same inputs and trains; otherwise it stores the run's hyperparameters for its
method, makes that method the user's choice and the current method (its inputs rebuilt with those values), sets
**Target** and **Features**, then `PreparationInputs.restore`
(Missing values with the impute settings, then Skip unique categories, One-hot encoding, Predict probability and the
cutoff, each only while shown; a hidden input keeps its default) and requests a check: the view retrains as after any
change, and the new training is a new first row.

The result (`TrainingResult`) carries the blob, the task, the metrics, the target and feature schemas, the
preparation options, the splitting, the seed, the hyperparameters, the row count (the prepared rows the model
learned from; the skipped rows are in `options.missingValues`) and `prediction`, the final model's prediction of
those rows as the method returned it (with Predict probability the scores, the same column as `scores.train`);
`trainPredictionOf(result)` gives its labels (the scores cut at the result's cutoff). The run and model rows record that row count, and the
dataset fingerprint is computed on the request's `features` and `target` (before the preparation steps). With
Predict probability it also keeps `scores` (`ProbabilityScores`: the method's train and out-of-fold predictions, the
selection's target, the positive class and both AUC-ROC values, measured once), so `recutResult(result, cutoff)` cuts
them again at another cutoff without retraining: it returns the result with `options.binaryClassificationThreshold`
and the metrics of the new labels; AUC-ROC does not change. A result without scores is returned as it is.

**The cutoff in the view.** A **Positive class cutoff** change requests no check: it fires the same `CHECK_DELAY_MS`
timer (`isCheckRequested` tells the two apart; a check wins), and `recut()` replaces `lastTraining.result` with
`recutResult` at the input's value and shows it, without a training or a run row; **Save** then writes the re-cut
`options` and metrics. A completed training is re-cut at the current value before it is shown (the cutoff may have
moved while it trained); its run row keeps the cutoff it trained with. An invalid cutoff (`cutoffProblem`) computes
nothing: `recut()` keeps the old result out of sight and shows `Fix the settings.`, **Save** is disabled, and
`problem()` blocks **Train** and the live training, as for any invalid setting. When the timer fires with no result
to re-cut, no training running and no finished training of the selection (another input changed while the cutoff
was invalid), it runs the check instead, so a valid cutoff then trains as any fixed setting does.

## Missing values

`preparation/missing-values.ts`. `prepareMissingValues(features, target, settings)` returns `PreparedData`: the
features (and target) to use, `keptRows` (the input rows kept, or `null` when every row was kept), `skippedRows`,
`imputedColumns` and `failedRows`. Rows with a missing target are always skipped, and the target is never imputed.
`settings.mode` is one of `MISSING_VALUES_MODES`:

- **`skip`** (Skip rows): rows with a missing value in any feature (or the target) are skipped. Only then are the
  features and the target copied (`columnRowCopy` in `preparation/row-copy.ts`); with nothing to skip the input
  columns themselves are returned. A masked clone of a text column keeps every category of the source, the empty
  one included, so the copies drop the categories no kept row has: otherwise the method would count a class for the
  skipped empty cells, and the saved target and data summary would list it. The imputed text copies are compacted the
  same way. The k-fold copies in `trainModel` are not: they keep the target's full category list, so every fold model
  has the classes of the final model.
- **`impute`** (Impute): the rows with a missing target are skipped first; then the feature gaps are filled by EDA's
  `Eda:knnImpute` (k nearest neighbors, in place). Because it writes in place, it gets a frame whose gapped columns are
  copies and whose other columns are still the user's (no extra copy when the target skip already made one), given
  back when the call ends. Cells it
  cannot fill (a row whose features are all missing) stay empty; those rows are then skipped and counted in
  `failedRows` and `skippedRows`. `settings.impute` holds `neighbors` and `distance`; `imputeSettingsOf(func)` lists
  the function's own inputs for them, and `defaultValuesOf` reads their defaults.
- **`keep`**: only the rows with a missing target are skipped; the features go to the method with their gaps. It is
  the default for a method with `meta.mlnativenans` (`handlesMissingValues(engine)`), with or without gaps: in the
  Train view and in `Forge:trainModel` (resolved after the suggestion may have replaced the chosen method). An
  explicit `skip` or `impute` (the API's `missingValues` option) is still honoured for such a method; `keep` for any
  other method is refused ("<Method> needs rows without missing values. Choose Skip rows or Impute.",
  `needsRowsWithoutGaps`). The Train view hides **Missing values** for such a method. Only when the training features
  had gaps do the Performance list (after training and of a saved model) say `Missing values: kept (the method
  handles them)` and the Runs cell `keep missing values` (see "Missing values in the view").

A model is applied (the Apply dialog, `Forge:applyModel`) and evaluated (Run Evaluation, `Forge:evaluateModel`) in the
mode it was trained with, `missingValuesOf(model)` in `apply/apply-model.ts`: the recorded `skip`, `keep`, or `impute`
with the recorded neighbors and distance; a model that records no mode gets `keep` when its method handles missing
values, else `skip`. The Apply dialog hides **Missing values** only when that mode is `keep`; otherwise it shows the
choice preset to the model's mode (`MissingValuesInputs.restore` when the model is chosen).

It takes and returns column lists; the returned `features` may be the user's own columns.

`missingValuesProblems(columns, settings)` refuses Impute when the function is missing ("Imputation is not
available: the EDA package has no knnImpute function."), for a yes/no column with gaps ("Impute cannot fill the yes/no
column 'X'. Choose Skip rows.") and with fewer than two features the imputer can read ("Impute needs at least two
features. Choose Skip rows or check more features."); it says nothing for Skip rows or for features without gaps.
`missingColumnsOf(columns)` lists the columns with gaps and their counts. Both take column lists, so the UI checks
build no frame.

At training, `prepareTraining` records the handling in `options`: `preprocessingInfo` gains `impute-missing` when
imputation ran and `ignore-missing` when a feature has gaps under Skip rows, or imputation left rows unfilled (the ids
of the built-in tool), and
`missingValues` is `{mode, neighbors?, distance?, skippedRows}`. `preprocessingInfo` is information only: an
application never replays it; `missingValues` sets the mode a model is applied and evaluated with (`missingValuesOf`,
see above). Imputing before cross-validation lets the held-out rows' values take part
in the imputation, as the built-in tool did; a later version may impute per fold.

The methods' own handling of gaps matters only for `mlnativenans` methods: the Training tests `<Method> trains and
applies with missing values under Skip rows and Impute` train and apply every other EDA method on iris with two gaps
in both modes and check that no table or column passed to `train` or `apply` has a missing value; `<Method> trains and
applies on the rows with missing values under keep` train XGBoost, Random Forest, Extra Trees, Hist Gradient
Boosting, CatBoost, Decision Tree and TabM on all 150 rows, check that the gaps reached `train` and that every row gets a
prediction; and an explicit Skip rows or Impute for XGBoost still keeps the gaps from it.

## Preparation

`preparation/pipeline.ts`. The steps after the missing values, in the built-in tool's order: skip unique categories,
one-hot, predict probability. `PreparationSteps` is `{oneHot, skipUniqueCategories, predictProbability, cutoff}`
(`DEFAULT_CUTOFF` 0.5).

`prepareFeatures(columns, target, steps, options)` (training) returns `PreparedFeatures` `{columns, target, options}`:
`options` is a copy of the given ones with each step that changed something recorded; a step that does not apply
changes nothing and is not recorded. Only Forge's own columns are created; the user's columns are never copied or
changed, and those no step touches are passed on as they are.

- **Skip unique categories**: drops the text and yes/no columns whose values are all different (`hasUniqueCategories`:
  `categories.length` equals the length), such as ids; `preprocessingInfo` gains `skip-unique-categories` and
  `options.skippedColumns` records the dropped names.
- **One-hot**: every remaining text or yes/no column becomes one int 0/1 column `<column>=<category>` per category, after
  the other columns, in the order of the columns and their categories; `preprocessingInfo` gains `one-hot`, and
  `options.oneHotCategories` records `{<column>: [categories]}` in that order.
- **Predict probability**, for a text or yes/no target with exactly two classes (`twoClasses`, the empty one aside): the method gets
  a float target named like the selection's, 1 for the first category (the positive class) and 0 for the other; float,
  because XGBoost rounds the predictions of an int target. `postprocessingInfo` gains `binary-classification`, and
  `positiveClass`, `negativeClass`, `binaryClassificationThreshold` (the cutoff) and `targetType` (the target's type)
  are recorded: the keys of the built-in tool's models. Training then measures the scores (see "Metrics"). As in the
  built-in tool, the score is the regression model's prediction on the 0/1 target, not a calibrated probability (it can
  fall outside 0..1); real class probabilities come when the EDA methods return them (phase 9).

`replayPreprocessing(columns, options)` (application) replays `options.preprocessingInfo` on the mapped columns and
returns the columns the method gets (the given list itself when there is nothing to replay): `ignore-missing` and
`impute-missing` are skipped, any other unknown id is refused ("The model uses the preparation step 'X', which Forge
cannot replay yet."). With `oneHotCategories`, `one-hot` builds exactly the recorded columns (a value training never
saw is 0 in each, a category the table lacks is a column of zeros); with `skippedColumns`, `skip-unique-categories`
drops exactly those columns by name, whatever their values (the application does not map them, see "Application
pipeline", so they are there only when passed directly); without it (the built-in tool's models) it drops the
all-unique columns (the core's rule); without `oneHotCategories` one-hot uses the applied data's categories. `replayPostprocessing`
(also in `pipeline.ts`) cuts the scores into the two classes (see "Application pipeline").

## Metrics

`metrics/metrics.ts`, ported from the built-in tool's definitions; labels and descriptions are the built-in ones (but
`R2` for the built-in "R squared"), values are shown with 3 decimals.

- Regression: `mse` = mean of squared errors, `rmse` = its root, `mae` = mean of absolute errors (new),
  `r2` ("R2") = 1 - SSres / SStot; for a constant target 1 if the predictions are exact, else 0.
- Classification: labels compare as text (a boolean `true` equals the label `'true'`). `accuracy` = correct / n.
  When at most two labels occur, the shares are computed for the positive class, the target's first category:
  `sensitivity` = TP / (TP + FN), `specificity` = TN / (TN + FP), `precision` = TP / (TP + FP),
  `npv` ("Negative Predicted Value") = TN / (TN + FN); a share whose denominator is 0 is left out, as the built-in
  tool hid it. `f1` (new): of the positive class for two classes, the macro average over the classes otherwise.
  For a two-category target the positive class is saved as `metrics.positiveClass`.
- Rows where the actual or the predicted value is empty are skipped; no rows left is a `ForgeError`.
- `auc` ("AUC-ROC"), Predict probability only: `aucOf(actual, score, positiveClass)`, the trapezoid over the rows
  sorted by score, highest first, the rows of one score taken as one step (one diagonal segment, as scikit-learn; the
  built-in ROC curve stepped row by row, so tied scores counted in their row order); undefined (left out) when the
  rows hold one class only. A
  Predict probability training measures its train and out-of-fold scores: the classification metrics above on the
  labels the cutoff gives (`replayPostprocessing` of the scores) against the selection's target, with its positive
  class, and `auc` on the scores themselves. The scores are regression outputs, not calibrated probabilities (real
  class probabilities come when the EDA methods return them, phase 9).
- `rocCurveOf(actual, score, positiveClass)` (the points of **ROC Curve**, ported from EDA's probabilistic MPO): at each
  of `ROC_THRESHOLDS` (0, 0.01, ..., 1) a row is positive when its score is at least the threshold (as the cutoff cuts
  the scores); `tpr` and `fpr` are float arrays of 101 values, a rate whose class has no row 0. The curve is drawn,
  not measured: `auc` stays the metric, so the AUC of the 101 points may differ slightly from it.

**Findings.** `metricFindings(metrics)` lists the values worth a word, per metric in `METRIC_IDS` order (a
`MetricFinding`: `id`, `kind`, the `splits` it marks, `text`, `detail` `Train 0.912, validation 0.405.`, 3 decimals, `-`
for a missing value):

- `threshold`: the splits whose value is below `METRIC_THRESHOLDS` (`r2` 0.3, `auc` 0.5, the built-in tool's):
  `R2 is below 0.3.` The built-in tool tested the train value only; Forge tests both splits.
- `mismatch`: a validation value far from the train one (`areValuesSimilar` false), marking the validation value:
  `The validation <label> is far from the training one: the model may be overfitting.`

`areValuesSimilar(a, b)` is the built-in tool's closeness: within 0.1 when the mean magnitude of the two is below 1,
otherwise a relative difference (to that mean) of at most 0.9, so only values several times apart differ.

**Statistics** (`metrics/statistics.ts`). `statisticsOf(actual, predicted)` gives a `StatisticRow` per id of
`STATISTIC_IDS` (`min, max, med, avg, stdev, variance, skew, kurt, q1, q2, q3`, the built-in tool's order; `q2` is the
median again) from `DG.Stats` of each column; the ids are shown as they are, `STATISTIC_DESCRIPTIONS` holds the
platform's texts of the statistics. `residualsOf(target, prediction)` is a float column `residuals`, prediction minus target per
row, empty where either is.

## Insights

`training/insights.ts`, the built-in tool's tips; an `Insight` is `{text, detail?, items?, panel?}`: the detail holds
the numbers (one line per `\n`), `items` a list of them, `panel` (`InsightPanel`: `classBalance`, `statistics`,
`correlation`, `residuals`, `confusions`, `metrics`) the results panel that shows what it says. No DOM: the UI shows
them.

`dataInsights(columns, target)` runs on the columns a method gets and the selection's target (`SelectionCheck.insights`;
with Predict probability that is the text target, not the 0/1 one, so its classes are checked too):

| Check | Columns | Text |
|---|---|---|
| Unbalanced classes | text and yes/no, the target first; a class whose rows are more than 20% off the average per class (empty cells and unused categories aside); the first four such columns | the target: `The classes are unbalanced.` (panel `classBalance`); a feature: `The classes of <column> are unbalanced.`; detail `<class>: 24.6%` per class, one decimal, rounded so the shares add up to exactly 100% (largest remainders); past ten classes the nine largest and `other <n> classes: <share>%` |
| Highly correlated | numerical (not dates), 2 to 99 of them; the pairs with Pearson correlation above 0.9 | one insight `Features are highly correlated.`; `items` `<a> and <b>: 0.96` per pair; panel `correlation` when a pair has no target |
| Almost unique | text and yes/no without a semantic type, more than 80% of the values different | `Almost every value of <column> is different; it will not help the model.`; detail `90% of its values are unique (above 80%).` |

The built-in tool's `<column> is text` tip is gone: rule 7 already reports such a feature, and the insights are
computed only when no target or features rule fails.

`resultInsights({values, target, prediction, positiveClass?})` checks a trained model's predictions of `target` (the
kept rows), `values` being its metrics:

- text or yes/no target: accuracy below 0.2 -> `Only 15% of the rows are predicted correctly.` (detail
  `Accuracy 0.150 (below 0.2).`, panel `metrics`); with sensitivity and specificity, one ratio of them below 0.6 ->
  `The model finds <negative> much better than <positive>.` or the other way round (detail `Sensitivity 0.900,
  specificity 0.300.`, panel `confusions`);
- numerical target: the residuals' quartiles: the top quarter more than four times as wide as the middle half ->
  `Large errors are mostly on one side: the model over-predicts.`, the bottom quarter -> `... under-predicts.` (detail
  the built-in arithmetic `Q4 - Q3 = a - b = c` and `Q3 - Q1 = ...`, panel `residuals`); then one insight for the
  statistics of `statisticsOf` whose predicted value is not similar to the actual one (`areValuesSimilar`; a value that
  is not finite is skipped): `The predictions' statistics deviate significantly from the actual ones.`, `items`
  `variance: predicted 0.214, actual 0.762` per such statistic, panel `statistics`.

The metric findings are not part of `resultInsights`: `evaluationInsights(result, stored)` (the Results pane and
`Forge:evaluateModel`) adds `metricInsights(metrics)` of the stored metrics after them, the findings as
insights with panel `metrics`: each `threshold` finding as it is (its text and detail), and all `mismatch` findings in
one `The model may be overfitting.` with `items` `<label>: train 0.912, validation 0.405` per metric (the grid's cell
tooltips keep the finding's own text per metric).

## Results

`ui/results-pane.ts`, `ResultsPane`: the Train view's results and the content of the `EvaluationView` (see "Run
Evaluation"). The pane holds the insights, the hint and the rows; the viewers and the data under them are built,
refilled and given back by `ui/results-viewers.ts` (`buildShown`, `refillShown`, `markRoc`, `updateWrongPredictions`,
`releaseShown`, `releaseLinkedFrame`; the panel titles `PANEL_TEXTS`, `resultsTitle`, `resultsPanel`). Its constructor takes the element shown next to the **Performance** title (the Train view's loader host;
the evaluation passes none). Two panels are built once and never rebuilt: the **Performance** panel
(`forge-results-panel forge-results-metrics`: its title with the status, then the metrics host `forge-results-body`
with a hint text or the `MetricsTable` root) and **Insights & Tips** (`forge-insights`: one `li` per insight, its text
only; its tooltip the detail lines, the `items` as a bulleted list and, last, `See <panel title>.` in grey
(`forge-tooltip-minor`) when the insight's `panel` is mounted; repeated texts once; hidden while empty). The insights of a viewer panel
(every `InsightPanel` but `metrics`, whose grid marks its own findings) are also its viewer's `description`
(`descriptionVisibilityMode: 'Always'`, the texts joined by a space, escaped as literal Markdown), refreshed with the
list. Every panel's title is the built-in tool's widget title, `ui.div([span], 'd4-dialog-header d4-dialog-title
forge-results-title')`; `PANEL_TEXTS` holds each title (`<name>` for the target's or the method's name) with what the
panel shows, the span's tooltip (the span takes the mouse, `d4-dialog-title` lets it through). The root holds the hint
host (the two panels while there is no result) and the rows host (`forge-results-rows`).

- `setDataInsights(insights)`: the check's `SelectionCheck.insights`, set after every check (also when a problem or a
  hint follows), shown first until the next call.
- `show(data: ResultsData)`: `engine`, `table`, `columns` (what the method saw), `target` (the selection's), `prediction`
  (labels or values; the Train view passes `trainPredictionOf(result)`), `score` (Predict probability), `keptRows`,
  `positiveClass`, `cutoff` (Predict probability: the Train view's `binaryClassificationThreshold` of the result, the
  evaluation's of the model), `blob`, `summary`. The grid is the `MetricsTable` in place at once (see "Components");
  the result insights are `evaluationInsights` of `summary.evaluation?.values ?? summary.metrics.train` and
  `summary.metrics`.
  Then the method's `visualize` is awaited, and the rows are mounted in one step: the previous rows stay until the
  whole result is ready, and a show that a newer show or hint outdated (a `version` counter) detaches the views it got
  and mounts nothing. A result of the same data keeps the viewers (see "Keeping the viewers" below); the method's own
  views are new every time.
- `showHint(text)`: the two panels go back to the hint host first (`append` moves without killing, `ui.empty` kills
  every widget under it), then the rows go, the frames are given back, the result insights go (the data insights stay)
  and the metrics host shows the text (`''` for nothing). The grid is detached, not disposed: the next result with the
  same metric rows shows it again. `clear()`: everything goes, the grid disposed. The view calls it on a Table change
  and in `detach()`.
- `frames`: the shared frames of the shown result (tests); kept until `clear()`, a hint or a result of other data.

`ui/results-layout.ts` (no DOM): `resultsLayout(kind, available, visualizeCount)` gives the rows, each a height, its
panel keys and their weights (the initial shares of the row's width). `kind` (`resultKindOf(target)`) is `regression` for a numerical target,
`binary` for exactly two non-empty classes (Predict probability too: the pane's target is the selection's text
column), `multiclass` otherwise. A key the result lacks is dropped and an empty row disappears; a row keeps its
template height. When row 2 is left with fewer than three panels, `classBalance` moves from row 3 to its right end
(weight 1.5). Weights in brackets, 1 where none is given:

| kind | row 1 (350 px) | row 2 (300 px) | row 3 (450 px) |
|---|---|---|---|
| `binary` | `predictedVsActual` (2), `rocCurve` (2), `metrics` | `insights`, `confusions` (1.5), `wrongPredictions` (2.5) | `features`, `correlation`, `classBalance` |
| `multiclass` | `predictedVsActual`, `confusions` | `insights`, `metrics` (1.5), `wrongPredictions` (2.5) | `features`, `correlation`, `classBalance` |
| `regression` | `predictedVsActual`, `residuals`, `distributions` | `insights`, `metrics` (1.5), `statistics` (1.5) | `features`, `correlation` |

Then `visualize:0`, `visualize:1`, ... in rows of at most `VISUALIZE_PER_ROW` (3), `VISUALIZE_ROW_HEIGHT` (350 px).

`mountRows(layout, panels)`: a row is a `ui.splitH(slots, null, true)` (`forge-results-row`, inline `height`) of slots,
`ui.box(panel)` each, so the splitter's dividers bind the slots themselves. A row whose signature (the keys and the
height) equals a mounted row's is reused: a slot that holds its panel already is left alone, another is emptied and
takes its panel (a panel the new layout uses elsewhere is moved out first, the two stable panels go back to the hint
host, the rest is disposed), so the widths and heights the user dragged survive; another row is built at its initial
sizes, and a mounted row that is not reused is emptied the same way. The rows go into a `ui.splitV(rows, null, true)`,
rebuilt only when the sequence of rows changes (its dividers bind the rows); a single row is mounted as is (a lone box
in a `ui-split-v` gets `flex-basis: 50%`). Platform facts this relies on (`js-api/ui.ts` `splitH` / `splitV`):

- once in the DOM, `splitH` sets each slot's `flexGrow` from the widths it measures (0 px in a hidden view), and on
  every size change of the row writes the rendered widths back as inline widths (`handleResize`), which can keep a
  passing layout; a drag sets both neighbours' widths and re-splits their `flexGrow` by the dragged widths. So the slots
  share the row by `flexGrow` alone (`flex-basis: 0`, `min-width: 0` on `.forge-results-row` slots): each slot is built
  with the growth of its weight and a width of `SLOT_WIDTH` (100 px) per weight unit, so whatever widths the splitter
  finds (the preset ones, or rendered ones in the same ratio) it derives the weights' ratio; a row whose slots all
  lost their growth (measured in a hidden view) gets the weights back;
- `splitV`'s `handleResize` (`ui.ts` ~1947) sets every row's height to `$(row).height()` plus a share of the size
  change, with no check that the root is in the document (`splitH` has one); `cash-dom` measures a detached or hidden
  element as 0 px, so switching to another view zeroed every row. The pane keeps each row's inline height while the
  rows are shown (a `MutationObserver` on the mounted rows' and slots' `style`, observed anew on every mount and
  disconnected when the rows go) and puts it back when it changes while they are hidden; a workaround to drop once
  `splitV` has `splitH`'s guard;
- the rows host is a `ui-div`, and the platform rule `div.ui-div > div.ui-box` (`ui.css`) makes a box in a div 400 x
  300 px: `div.ui-div.forge-results-rows > div.ui-box` sets the split root to the full width and its content's height
  (the sum of the rows), so the scrolling `forge-train-results` panel scrolls the rows and the splitters resize inside;
- `div.ui-box > div.ui-div` makes a panel (a `ui.divV`) fill its slot and scroll, so the title scrolled away with the
  content: Forge's `div.ui-box > div.ui-div.forge-results-panel` keeps the panel from scrolling (`overflow: hidden`,
  `!important` as the platform's); the viewer, the panel's last child (a `ui-box`, so 400 x 300 px by the same platform
  rule), fills the panel's width and the rest under the title (`width: 100%`, `flex: 1 1 0`, a grid's own height
  overridden), and in a row the Performance panel's body and the insights list do the same and scroll inside.

One frame, `sharedFrame([...columns, target, prediction, marks, score?])`, holds the user's columns where the method
saw them (no row skipped, no preparation step) and Forge's own columns: the prediction renamed `Predicted <target>`,
`residuals` (`residualsOf`, a regression) or `correctness` (a bool, right label or not; markers `true` circle border,
`false` cross; a classifier) and `score`. Forge's columns get a name neither the table nor the frame has (` 2`, ` 3`,
... appended). The panels (the built-in tool's titles and conditions):

| Target | Panel | Viewer |
|---|---|---|
| numerical | **Predicted <target> vs Actual** | `scatterPlot`: target by prediction, regression line |
| | **Residuals** | `scatterPlot`: residuals by prediction, regression line, y from `-m` to `m`, `m` the target's largest magnitude |
| | **Distributions of predicted and actual <target>** | `histogram` of Forge's own frame of 2N rows: `targets` (the target, then the prediction), split by `label` (`TARGET` / `PREDICTION`) |
| | **Statistics** | `readOnlyGrid` of **Statistic / Actual / Predicted** (`statisticsOf`, `0.000`), header tooltips, the statistic's description on **Statistic** |
| text or yes/no | **ROC Curve** (with Predict probability) | `scatterPlot` over Forge's own frame of `Threshold`, `FPR`, `TPR` (`rocCurveOf` of the score for the positive class, 101 rows), as EDA's probabilistic MPO: a line (`linesOrderColumnName` `FPR`, `linesWidth` 5, `markerType` `dot`), the dashed formula line `${TPR} = ${FPR}` `Non-informative baseline`, no column selectors; the frame's current row is the cutoff's threshold (`round(cutoff * 100)`); the panel title `ROC Curve (AUC = 0.912)`, the AUC-ROC of **Performance** (`summary.evaluation?.values.auc ?? summary.metrics.train.auc`) |
| | **Predicted <target> vs Actual** | `scatterPlot` of the first two numerical columns, colored by the target, markers by `correctness` (two or more numerical columns) |
| | **Confusions** | `confusionMatrix`: target by prediction (fewer than 10 classes) |
| | **Wrong predictions** | when a row is wrong: a `readOnlyGrid` of a second shared frame, every column of the table and the prediction spread over its rows (`scatteredColumn`; the same column when no row was skipped), filtered to the wrong rows, with the row header (the table's row numbers); header tooltip `The model's prediction for the row.` |
| | **Class balance** | `barChart` of the same frame: split and stacked by the target, the row count (`count`), so one bar per class in the class's own color (the target column's colors), the Stack, Category and Value selectors hidden; the kept rows only, as every viewer of the frame; no predicted bars |
| both | **Features** | `pcPlot` of the numerical columns (and, for a numerical target, the target and the prediction), colored by the target (3 to 20 numerical columns) |
| | **Features correlation** | `correlationPlot` of the numerical columns (3 to 20, as **Features**) |
| | the method's name | the method's `visualize(df, targetColumn, predictColumn, model)` (on a shared frame of `columns`), one panel per element it returns (EDA's PLS: its viewers' roots), rows 4+; a failure is reported once (`reportError`, repeating) and the other panels are mounted |

The hovered row is linked both ways between each frame and the table: through `keptRows` (the frame's row `i` is the
table's `kept[i]`, the inverse `-1` for a skipped row), directly for the wrong-predictions frame; a value is set only
when it differs, and the events a set fires are not passed back (a skipped table row would come back as none). A
frame's extra parent changes nothing for the user's columns (`Column.parent` is the first one, the table).

Keeping the viewers: the pane keeps the shown result's data (table, kind, `columns`, `target`, `keptRows`), its frame,
Forge's columns, Forge's own frames of **Distributions ...**, **Statistics** and **ROC Curve** (with its title), the
wrong-predictions frame and the viewer panels. A new result whose table and kind are the same, with or without a score as before, the same kept rows
and the same columns and target (the same objects, or copies with the same name, type and values: a skipped row or a
preparation step gives new copies on every training) only refills: the new prediction and score are written into the
old columns in place (`copyInto` of `apply/apply-model.ts`, or `init` for a text column so its categories are rebuilt, then `fireValuesChanged`; only when the type and the length match, otherwise everything
is rebuilt; a score equal to the shown one, as after a cutoff change, is not written; a new score refills the ROC
frame's `FPR` and `TPR` in place, its 101 rows kept), then the marks, the distribution values and the statistics are recomputed in place, so the viewers redraw with
their state (zoom, sizes) kept. A scatter plot recomputes its regression line only when the state of its inputs
changes (its columns, its filter's version, its settings: `getRegressionState`, `d4/.../scatterplot/features/
regression_line.dart`), not when the values of its columns do, so a regression refill fires the frame's filter
(`filter.fireChanged()`, the filter itself unchanged). **Wrong predictions** follows the new prediction: its scattered
column refilled and its filter re-set while a row is wrong, built when one becomes wrong, dropped (and its frame given
back after the mount) when none is. **ROC Curve** is on its own frame, so no refill of the preview frame reaches it: a
cutoff change only moves its current row, and the title and the current row are set only when they change (on the
platform's ROC curve viewer, which shared the preview frame, every refill of the prediction and the marks redrew the
curve). A hyperparameter, Method or cutoff change gives such a result; a change of the
data (Table, Target, Features, other skipped rows or imputed values, Predict probability, another task kind) builds
everything again. A panel that appears or goes (Wrong predictions, Insights & Tips, the method's own views) changes
only the rows it is in, and row 3 too when **Class balance** moves between rows 2 and 3. A mount disposes the viewers it no longer shows with `ui.empty` of their slots and rows (the platform's
`Widget.kill`; `Viewer.close` works only in a table view); after a rebuild the previous result's subscriptions are
dropped and its shared frames given back (`releaseFrame`); a hint and `clear` do the same for the shown result. A
viewer or the **Wrong predictions** grid that fails to build gives its frame and hover links back before the error
goes on (the shown result stays).

## Run Evaluation

`ui/evaluation-view.ts`. `runEvaluation(model)` (the Performance button; a UI boundary, `reportError`) loads the model
(`loadedModelOf`) and finds data with known targets with `storedTrainingTable` (`training/evaluation.ts`, shared with
`Forge:evaluateModel`): a stored reference (`storedDatasetRef` -> `openDatasetRef`), else the uploaded copy
(`grok.dapi.tables.getTable(dataset_table_id)`); neither adds the table to the workspace. Such a
table is evaluated on its exact column names (`exactMapping`) and its column `target_name` (none: `The data source has
no column '<target>' the model predicts. Choose another table.`). A model that keeps neither gets
`evaluateModelDialog(model, loaded, tables)` over the open tables it fits (`applicableTables`; none: `No open table has
the columns the model needs. Open its training data or a similar table.`): **Evaluate model** with **Table** (the
current table when listed, else the first) and **Target** (the columns `compatibility` accepts for the stored target
schema, `targetSchemaOf(model.target)`, else a `double` for a regression and a `string` otherwise; `target_name`
preset when it passes), refilled on a Table change; OK, kept disabled by a `ButtonGate` (`Choose a table.` / `Choose a
target.`), evaluates with `suggestMapping`.

`EvaluationView` (a `DG.ViewBase`, tab `Evaluation: <model>` with the model icon) holds a `ResultsPane` (`resultsPane`)
in a scrolling `forge-train-results` panel: the same rows as a training, the **Performance** header without a loader.
`EvaluationView.open(model, loaded, table, target, mapping)` adds the view, shows `Evaluating the model...` and runs
`evaluateModel` under the cancellable task-bar progress `Evaluating <model>` (closed when it ends). A cancel
(`progress.canceled`; the batch loop's message is not read) shows `Evaluation was cancelled.` and closes the view;
another failure is reported and closes it. The result goes to `ResultsPane.show` with the stored summary
(`storedSummary`, as Performance reads it) plus `evaluation` (the values, the table's name, the kept rows and the
skipped ones), so the grid has the **Evaluation** column, the bullets end with `Evaluation: <table> (<rows> rows)` and
the result insights read the evaluation's values; the positive class is the stored one, else the options'. The view's
`detach()` clears the pane (its frames given back); a result that comes after the view closed is dropped. The evaluated
table is held by the view only; no column is added and no application is recorded.

## Saving

`storage/model-store.ts`. `saveModel(fields, blob)` writes the blob to
`System:DomainFiles/forge/model/<uuid>/model.bin` first and then inserts the `model` row with
`blob = file://<that path>`. A failed insert leaves an unreachable file, never a row without its blob.
`deleteModel(id)` soft-deletes the row; it deletes a folder only when `blob` matches
`file://System:DomainFiles/forge/model/<uuid>/<file>`, and then exactly `System:DomainFiles/forge/model/<uuid>`. Any
other value (a file picked in the generic row editor elsewhere, an empty segment, `..`, no blob) leaves the files
alone. It also deletes the uploaded copy of a `copy` model (`deleteTrainingCopy` of `dataset_table_id`, see "Storage
modes"); both deletes are tried whatever the other does, the first failure is thrown after the row is gone, and
`modelsChanged` fires anyway. `saveTrainedModel(model)` is the save of a training, shared by the ribbon **Save** and
`Forge:trainModel`: `copy` uploads the user's feature and target columns first (`uploadTrainingCopy`) and, when the
save then fails, deletes the copy (`deleteTrainingCopy`; its own failure is dropped) before rethrowing the save's
error; `reference` stores the origin it is given; then `saveModel` and `linkTrainingRun`. Nothing is uploaded
without `copy`. A `linkTrainingRun` failure is thrown after the row is written and keeps the model (and its copy): the
view's **Save** is enabled again, so a second Save makes a second model, and `Forge:trainModel` throws with the model
saved and returns no `modelId`. `findModelRow(idOrName, columns)` resolves a model for the API functions: the visible model with that
id (a uuid) or that name, read in one query of at most two rows ("No Forge model 'X' is available to you." /
"Several models are named 'X'. Use the model id."). In the UI, the ribbon **Save** asks for Name, Description, Tags
and **Data storage**, then calls `saveTrainedModel`; the catalog's delete and **Delete model** call `deleteModel`.
Both storage functions fire `modelsChanged`.

**Data storage in the Save dialog** (`StorageOffer` of `saveModelDialog`): a radio (`forge-inline-radio`, tooltip `What
the model keeps of its training data.`) of `Reference` (only when `datasetRefOf(table)` gives a reference; then the
default), `None` (the default otherwise) and `Copy`, with one `forge-note` line under it that follows the choice:
`A link to <path | query name | a script> is saved; the data stays where it is.` / `Only a summary of the
data is saved.` / `The training columns (<N> rows) will be uploaded to the server.` (N: the table's rows); there are no warnings
about a model file with training rows or a method that runs on the server (a platform setting is planned).
OK passes the choice (`StorageChoice`; `reference` carries the reference the dialog offered, so `datasetRefOf` runs
once per Save) to `saveModelAs(info, choice)`, which saves through `saveTrainedModel` (above).

`model-fields.ts` builds the row payloads: `modelFieldsOf` (the storage mode and its `dataset_ref` or
`dataset_table_id` from a `ModelStorage`, `has_training_rows` from the method (`hasTrainingRows`), feature and row
counts, dataset name, the tags as their column text, see "Tags") and `trainingRunOf`. The json-column types are listed
in "JSON column shapes" below.

**Dataset fingerprint** (`dataset-fingerprint.ts`), the only trace of the training data on the server, computed in
the browser over the features in training order and the target last: `rowCount`, `columnCount`, an 8-hex-digit hash
and per column `name`, `type`, `missingCount`, plus `min`/`max`/`mean` for numerical columns and `categories` for
text and boolean columns with up to 20 of them. The hash is FNV-1a (32 bit) over each column's name and type and
then its data: the raw buffer cut to the column length for numerical columns, the values as text otherwise (bigint
columns as text). It depends on the column order.

## Application pipeline

`apply/apply-model.ts`. `loadModel(idOrName)` reads, in one query (`findModelRow`, see "Saving"), the
`APPLY_COLUMNS` of the model with that id, or of the one visible model with that name ("No Forge model 'X' is
available to you." / "Several models are named 'X'. Use the model id."). `loadedModelOf(row, engines)` then refuses, in this order, a row without a blob ("has no model
file"), a blob outside Forge's own `forge/model/<uuid>/` layout (`ownBlob` in `storage/model-store.ts`), a row
without a `{name, type}` feature list, and an engine that is not discovered or incomplete ("The method 'X' is not
installed. Install the <package> package."). The result (`LoadedModel`) holds the row, the engine, the features a
table must provide (`requiredFeaturesOf`: the row's feature list without `options.skippedColumns`, which the model
keeps as the record of the user's choice but never needs), the options (`preparationOptionsOf`) and the blob path;
the target is the row's `target_name`.

An `ApplyRequest` is the model, the table, the column mapping, the batch size (`DEFAULT_BATCH_SIZE` 10000) and the
missing-value settings. `applyAndRecord(request, source, progress?)`:

1. Refuses the request when `mappingProblems` reports anything or the table has no rows (`checkApplyRequest`); a
   refusal is not recorded.
2. Takes the feature columns (`featureColumns`, `apply/feature-columns.ts`): the mapped columns of the table in training
   order, the table's own. A column whose name differs from its feature is a renamed copy: renaming the table's
   column would rename it in the user's table, and the engines find features by name.
3. `prepareMissingValues` without a target (see "Missing values"); every row skipped is a `ForgeError` ("Every row
   has a missing value in the columns the model needs, so nothing can be predicted. Fill the missing values or
   choose Impute.").
4. `replayPreprocessing` (see "Preparation") replays `options.preprocessingInfo`: one-hot with the recorded training
   categories, or the applied data's for a model without the record; without steps the columns themselves are used.
5. Reads the blob and calls the engine's `apply`: once on the columns themselves when they fit one batch, otherwise
   once per batch of rows (a copy each, of a range mask built from a word array), the results written into one column
   (numerical predictions through their raw data). Before every call a cancelled progress
   indicator stops it ("Application was cancelled."); after every call the progress shows "Rows a-b of n". At most
   every 50 ms the loop gives the event loop a turn (`yieldToEventLoop`), so a click on the cancel and the progress
   repaint get through: fast engine calls resume as microtasks and would otherwise hold the page to the last batch.
6. `replayPostprocessing` replays `binary-classification`: scores at or above `binaryClassificationThreshold` become
   `positiveClass`, the others `negativeClass`, in a column of `targetType`.
7. With skipped rows, the predictions are spread over a full-length column (`scatteredColumn`); skipped rows stay
   empty.
8. Names the column `predictionName(table, target)`: `<target> (predicted)`, or `<target> (predicted 2)`, ... when a
   column of that name exists (ignoring case); tags it `forge.model` = the model id (`PREDICTION_TAG` in
   `constants.ts`) and adds it to the table.

Steps 2-6 are `predictRows(request, progress?)`, which adds nothing to the table: it returns `RowPredictions`, the
columns the method got, its own output (`raw`), the labels (`raw` with the postprocessing replayed, `raw` itself
without a step), `keptRows`, `skippedRows` and the blob read. `targetSchemaOf(model.target)` reads a stored target
schema with its categories.

**Evaluation** (`training/evaluation.ts`). `evaluateModel(request, progress?)` measures a saved model on a table whose
`target` column holds the actual values (an `EvaluationRequest`: the `LoadedModel`, the table, the target, the column
mapping, the positive class): `checkApplyRequest`, then `predictRows` with the default batch size and the
missing-values mode the model was trained with (`missingValuesOf`, see "Missing values"); the target is the given column, or a copy of its kept rows
(`columnRowCopy`). A Predict probability model is measured as at training: the classification metrics of the labels with the positive class (the request's, else the options') plus
`auc` of the scores (`probabilityValues`, shared with training), the scores coming back as `score`; any other model
uses `metricsOf(taskOf(target), ...)` (`taskOf`: regression for a numerical target). A
row without an actual value is left out by the metrics. The `EvaluationResult` holds the columns, the target and the
prediction of the kept rows, `keptRows`, `skippedRows`, `rowCount` (the kept rows), the metric `values`, the blob
and `positiveClass` (the request's, else the options'). Nothing is added to the table and no application is recorded.

No step writes into the table's columns: the engines read their inputs and return new frames, the replay builds new
lists with Forge's own one-hot columns, and the imputer only gets copies of the columns it fills. The pipeline builds
no frame of its own; the imputer's and the engine's frames are given back by the calls that build them. The raw
engine call
(`apply` in `engines/engine-calls.ts`) takes the prediction column out of the engine's result frame, so the user's
table is the column's only parent (`column.dataFrame` is the table). Only the first column of the engine's result is
used; other output columns are dropped.

After the pipeline, `applyAndRecord` writes one `application` row (`apply/application-store.ts`, which then fires
`applicationRecorded` with the model's id): `completed` with
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
Choose a column with numbers."), by the same `isReadableNumber` rule and `bigIntProblem` sentence as training; a
semantic type the column lacks is only a hint ("'X' is not marked as <SemType>; check that it is the same kind of
value."). A text, boolean or date feature needs a column of the same type ("'X'
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

`ui/apply-model-dialog.ts`. `openApplyDialog(table, {modelId?, switchToTable?, preferredTable?})` is the boundary of
`ML | Forge | Apply...` (the current table), of a model's **Apply...** and the catalog's play icon (the current
table, else the first open one) and of the table pane's run icon (the pane's table, also as `preferredTable`);
without a table it says "Open a table first.". It calls `applyModelDialog({table,
modelId?, switchToTable?, preferredTable?})`, which reads up to 10000 visible model rows once (`MAX_MODELS` of
`constants.ts`), only the
columns applying needs (`APPLY_COLUMNS`; no models: "No Forge models yet. Train and save one with ML | Forge |
Train... first.") and discovers the engines once. With a `modelId` among them, the dialog opens on `presetTable`:
`preferredTable` (the catalog's **Applicable to**) while it is open, else the current table if the model fits it
(`applicableTables`), else the first open table it fits, else `table`; without one (`ML | Forge | Apply...`) on `table`.
The dialog **Apply predictive model** has:

- **Table** (`Table to add the prediction to.`): a change re-orders **Model**, keeps the chosen model and re-prefills
  the rows; the table's column add and remove events re-check the rows while the dialog is open, through a change
  event of **Batch size**, so the dialog's own Enter check (refreshed only on input changes) follows too.
- **Model** (`Saved model to apply.`): labels are unique (`modelLabels`): the model name, with ` (YYYY-MM-DD HH:mm)`
  when models share a name, to the second when they share the minute too, then ` #2`, ` #3`, ... in creation order;
  the dialog keys the models by these labels. The models that fit the table (`isSuggested` over the row's
  `requiredFeaturesOf`, computed once per distinct list) come first, newest first in each group; the default is
  the preset model, else the first. A
  model `loadedModelOf` refuses shows its message instead of the rows, and the same message marks **Model**.
- **Columns** (a `CollapsibleGroup`): the summary `N of M matched` (features without a mapping problem out of the
  model's `LoadedModel.features`, so the columns Skip unique categories left out are neither asked for nor counted) and,
  inside, one column input per such feature, captioned with the feature name, in a block (`forge-apply-rows`) that
  scrolls past 40% of the window height (at least 160 px). It starts collapsed when every feature has a valid column, expanded
  otherwise; the summary turns red while any feature has a problem, and a row validated with an error expands it.
  The list offers only columns `compatibility` does not refuse; the prefill is `suggestMapping`; the tooltip is
  "Table column used as the feature 'F' (<kind>)." plus "'X' has N missing values." and a semantic-type hint when
  they apply. The validator is the row's `mappingProblems` message. The dialog keeps every input it was given, so one
  input per feature name is reused across models; an input of a feature the chosen model lacks is emptied and made
  nullable. Folded rows are still validated by the dialog and still gate **OK**. The rows block opens scrolled to the
  top: the dialog focuses the input it was given last on `show()`, so the rows are given to it before the other
  inputs, and the scroll is reset when the rows are rebuilt or **Columns** opens.
- **Missing values** with **Neighbors** and **Distance** (see Components), following the mapped columns; its validator
  is `missingValuesProblems` for Impute. It stays outside the groups. It is preset to the mode the model was trained
  with (`restore(missingValuesOf(model))` once the model's rows are filled), and hidden when that mode is `keep`: the
  rows with gaps are then predicted.
- **More options** (a `CollapsibleGroup`, collapsed): **Batch size** (10000, 1..100000, `Rows predicted in one step.
  Lower it for heavy methods.`); an invalid batch size expands it.

After every change the dialog recomputes the problems, re-validates the visible inputs and disables **OK** with the
first problem in form order as its tooltip (the model refusal, the rows' messages in feature order, then any other
invalid input). The gaps of the mapped columns are counted once per change and feed both **Missing values** and the
rows' tooltips. **OK** closes the dialog at once (its handler only starts the work, as the built-in dialog did) and
runs `applyWithProgress(..., 'ui')` from `applyAndReport`, outside the form, so the running application does not hold
the form; then it shows "Added the column "C" to T." ("...; N rows skipped (missing values)." when rows were
skipped) and, from the catalog, switches to the table's view. Errors go to `reportError` (a cancel is the yellow
"Application was cancelled.").

## `Forge:applyModel`

`applyModel(model, table, columnNamesMap?, showProgress = true) -> table`, for scripts and other packages, never opens
a dialog: `runApplyModel` resolves the model by id or name and, in `apiMapping` (shared with `Forge:evaluateModel`),
maps the features by exact names, overlays the given `columnNamesMap` pairs, and refuses the mapping problems before
any work, the unmapped features in one sentence first
("The table has no column 'a' the model needs. Map it in columnNamesMap." / "The table has no columns 'a', 'b' the
model needs. Map them in columnNamesMap."). Missing values are handled in the mode the model was trained with
(`missingValuesOf`, see "Missing values"); batches of 10000 rows; with
`showProgress` a cancellable task-bar indicator "Predicting <target>". The application is recorded with `source: api`.

## `Forge:trainModel`

`trainModel(table, target, features?, method?, hyperparameters?, options?, name?, description?, tags?, storage =
'none', showProgress = true) -> object`, for scripts, other packages and fixtures, is the Train view's pipeline
without the view (`runTrainModel` in `training/train-api.ts`; the view and the function share, from
`training/train-model.ts`, `trainingEngines`, `FOLDS`, `randomSeed`, `problemMessages`, `checkSelection`,
`prepareTraining` and `trainModel`, and `defaultFeatures`, `DEFAULT_STEPS`, `defaultOneHot`, `cutoffProblem`,
`modelNameOf` (`NO_NAME`; the name is written trimmed) and `saveTrainedModel`):

1. The target and the features by name ("The table has no column 'x'." / "The table has no columns 'a', 'b'."); without
   `features`, `defaultFeatures`.
2. The method by name among the complete methods ("No method 'x' is installed. The methods are: ..."); without one,
   the check starts with XGBoost (`PREFERRED_ENGINE` in `train-model.ts`) and takes the method it suggests (`best`), as the view does.
3. `options`, narrowed key by key (`OPTION_NAMES`; an unknown key is refused with the list) to the view's defaults:
   `missingValues` left unset when not given (the selection's `missingValues` is then undefined and
   `selectedMissingValues` gives the training method's `defaultMissingValues`: `keep` for a method that handles
   missing values, with or without gaps, else `skip`; `impute` with `neighbors` and `distance` checked against
   `Eda:knnImpute`'s inputs, its defaults otherwise; an explicit `keep` for a method that does not handle missing
   values is the check's Missing values problem, `needsRowsWithoutGaps`: when the suggestion replaces the chosen
   method, the selection is checked again against the suggested one), `oneHot` `defaultOneHot`, `skipUniqueCategories` `true`, `predictProbability`
   `false` (refused on a target without exactly two classes), `cutoff` `DEFAULT_CUTOFF` (with Predict probability,
   0..1 or the view's `Enter a number from 0 to 1.`), `seed` random, `folds` 5.
4. The save, checked before any work: a blank name ("Enter a name for the model."), an unknown storage, `reference`
   without `datasetRefOf(table)`.
5. `checkSelection`; every problem goes into one `ForgeError`, target and features first.
6. The hyperparameters: `defaultHyperparameters` overlaid with the given ones, each checked against its `train`
   input (an unknown name with the list of settings; the kind, a whole number for an int, `min`/`max`, `choices`).
7. With `showProgress`, the task-bar indicator `Training <Method> model` (cancellable) is the `LoopProgress` of
   `trainModel`; `prepareTraining`, then `trainModel`.
8. The `training_run` row (`source: api`): `completed` with the metrics; `cancelled` or `failed` (with the error)
   before the error is rethrown. A refused input or a failed preparation records nothing.
9. With a name, `saveTrainedModel` with the user's feature and target columns, which links the run; nothing is
   uploaded without `storage: copy`. When the link fails after the row is written, the error is thrown with the model
   saved, so the caller gets no `modelId` for it (the model is found by its name).

It returns `TrainModelResult` (`modelId` or null, `runId`, `name`, `method`, `task`, `target`, `features`, `rowCount`,
`skippedRows`, `seed`, `metrics`) and adds nothing to the table. Errors are rethrown unchanged (`trainAndSave` records
a failed or cancelled run first): `ForgeError` texts reach the caller as they are.

## `Forge:deleteModel`

`deleteModel(model) -> id` (`runDeleteModel` in `storage/model-store.ts`): `findModelRow` resolves the id or the
unique name with the row's `blob` and `dataset_table_id` (one read), then `deleteModelRow` deletes that row as the
catalog's delete does (`deleteModel(id)` reads the row itself first): the row, its own blob folder, its uploaded copy; the
application rows go with the row (`onDelete: cascade`), the training runs keep their rows without the model
(`setnull`). The server's refusal of a row the caller may not delete passes through unchanged.

## `Forge:listMethods`

`listMethods(table?, target?) -> dataframe` (`runListMethods` in `training/train-api.ts`, which adds the table's
columns to `methodList` of `engines/method-list.ts`) lists the complete methods
(the ones `trainModel` accepts) in discovery order as the frame **Methods**:

1. `methodColumns` gives the app grid's five columns (Method, Package, Method type, Roles, Hyperparameters); the app
   builds its grid from them (`methodsFrame`) over every discovered method. The frame is built at once
   (`DG.DataFrame.fromColumns`) in the order below.
2. **Runs on** (`server` when `isServerEngine`, else `browser`), **Keeps training rows** (`hasTrainingRows`) and
   **Retrains live** (`isLiveUpdate`) go after Method type; **Settings** last: per method the json text of
   `methodSettings(engine)`, one `{name, type, default?, min?, max?, choices?, description?}` per `train` input
   (`hyperparametersOf`; `default` from `defaultHyperparameters`, a key the input does not declare is left out).
3. With a `table`: the target column is `target` by name ("The table has no column 'x'.", as `trainModel`), by
   default the table's last column (the Train view's default; no columns: "The table has no columns."); the features are `defaultFeatures(table, target)`.
   Then `trainModel`'s steps without a method: `trainingEngines`, `optionsOf({})`, `checkSelection` and `engineOf`,
   which throws the check's problems as `trainModel` does. **Applicable** is the check's `engines` (the prepared
   columns, an engine whose check throws is not applicable), **Suggested** the one `engineOf` gives (`best`, else the
   first applicable). Nothing is copied: the checks run on a shared frame of the prepared columns.

`grok.functions.call('Forge:listMethods', {})` is accepted with the optional `dataframe` input omitted (no `nullable`
needed); the body treats a missing `table` or an empty `target` as not given.

## `Forge:listModels`

`listModels(name?, target?, task?, method?, tags?, table?, limit = 100) -> dataframe` (`runListModels` in
`catalog/list-models.ts`) lists the visible models, newest first, as the frame **Forge models**:

1. A `task` other than `regression` or `classification` fails with "The task must be regression or classification.",
   a `limit` that is not a whole number of at least 1 with "The limit must be a whole number of at least 1.". An empty
   text counts as not given.
2. One query of the `model` columns `LIST_COLUMNS`, newest first, at most `MAX_MODELS` (10000, the cap the Apply
   dialog and the table pane share; at most `limit` when no browser filter below is given): `task` and `method` (the
   `engine_name`, exact; a method no longer installed lists its saved models) filter on the server.
3. In the browser, one pass over the rows that stops at `limit` matches: `name` is a part of the name, `target` the
   whole target name, both ignoring case; `tags` (a comma-separated text, `tagsOf`) keeps a model that has every one
   of them, ignoring case; `table` keeps the models that fit it (`applicableTables`; nothing is uploaded).
4. The authors' logins come from `loginsOf`, one `grok.dapi.getEntities` read (a deleted user gives an empty cell). The columns:
   Id, Name, Description, Method, Task, Target, Features (the feature names, comma-separated), Training rows, Data
   storage (`none`, `reference`, `copy`), Tags (the normalized text, `tagsText`), Author, Created, then `<Metric> (validation)` for
   every metric a listed model has (`metricColumns` of `catalog/compare-models.ts`, which the Compare view uses for
   both splits).

## `Forge:getModel`

`getModel(model) -> object` (`runGetModel(idOrName, tables)` in `catalog/get-model.ts`; `package.ts` passes
`grok.shell.tables`) reads the model as the Details, Performance and Activity panes show it, in one object
(`ModelDetails`): `findModelRow` with every column but `legacy_id` (`GET_MODEL_COLUMNS`), then `modelActivity` and the
author (`loginsOf`; a deleted or unknown user gives an empty `author`) together. The keys: `id`, `version`, `name`,
`description`, `tags` (a list, `tagsOf`), `method`, `package` (the engine's namespace), `task`, `target`, `targetType`
and `classes` (`targetSchemaOf`; null without them), `features` (in training order) and `featureTypes` (`{name: type}`)
from one `featureSchemasOf`, `applicableTo` (the names of the open tables that fit, `applicableTables`), `hyperparameters`,
`preparation` (`preparationOptionsOf` of `options`), `metrics` (`metricsRecordOf`, null without them), `findings`
(`{text, detail}` of `metricFindings`, none without metrics), `seed`, `splitting`, `storage` (`mode`, with a reference
its `source`, `datasetRefCaption`, with a copy its `tableId`), `dataset` (`name`, `rows`, `summary`: the stored
fingerprint), `hasModelFile` (`ownBlob`), `keepsTrainingRows`, `author` (login), `created` and `updated` (ISO texts),
`lastRun` (`{when, status}` of the newest application, null when never applied) and `applications` (the count). The
model is resolved as for `deleteModel`: "No Forge model 'X' is available to you.", "Several models are named 'X'. Use
the model id."

## `Forge:updateModel`

`updateModel(model, name?, description?, tags?) -> object` (`runUpdateModel` in `catalog/model-edit.ts`) writes what
**Edit model...** writes, without a dialog: a `name` blank after trimming fails with "Enter a name for the model."
(`modelNameOf`, as `trainModel`: `NO_NAME`, the name is written trimmed), a call with none of the three with "Give a name, a description or tags to
change.", both before the model is read; `tags` is a comma-separated text that replaces the current tags (`tagsOf`;
an empty text removes them); a field not given stays. `findModelRow` resolves the model, `updateModelInfo` writes the
given fields without a version (the last write wins) and fires `modelsChanged`, so an open catalog reloads itself and
a context panel showing the model shows it again (`followChanges`). It returns `{id, name, description, tags,
version}` after the write. `Forge:deleteModel` fires it too: the panel showing the deleted model is emptied.

## `Forge:evaluateModel`

`evaluateModel(model, table?, target?, columnNamesMap?, showProgress = true) -> object` (`runEvaluateModel` in
`training/evaluate-api.ts`) is Run Evaluation without the view:

1. `findModelRow` reads the `EVALUATE_COLUMNS` (`APPLY_COLUMNS` plus `metrics`, `dataset_ref`, `dataset_table_id`),
   `loadedModelOf` loads it (its refusals as for `applyModel`).
2. The table: the given one, else `storedTrainingTable` (the reference, else the copy; not added to the workspace),
   else "The model keeps no training data. Pass a table with the actual values."
3. The target: the column named `target`, by default `target_name` (an empty text counts as not given); none: "The
   table has no column 'x' with the actual values. Pass the column's name in target."
4. The features: `apiMapping`, the mapping and the refusals of `Forge:applyModel`.
5. `evaluateModel` with the stored positive class (it falls back to the options' one and returns the class it used
   as `positiveClass`, which the Evaluation view reads too), under the cancellable task-bar progress
   `Evaluating <model>` unless `showProgress` is false (the batch loop relabels it `Rows i-j of n`); an error after a
   cancel becomes "Evaluation was cancelled.".

It returns `EvaluateModelResult`: `id`, `name`, `table` (the evaluated table's name), `target`, `rowCount`,
`skippedRows`, `positiveClass` (when there is one), `metrics` (`train` and `validation` from the stored record, `{}`
each without one, and `evaluation`) and `insights`: `evaluationInsights` of the evaluation and the
stored metrics, as the Evaluation view's Insights & Tips, as `{text, detail?, items?}` (no `panel`). No column is
added and no application is recorded.

## `Forge:shareModel`

`shareModel(model)` (`openShareDialog` in `ui/model-actions.ts`, a UI function) resolves the model with
`findModelRow` (a missing model rejects the call, no dialog) and opens the platform's Share dialog
(`DG.DomainObjectHandler.shareRow` of `ForgeModelHandler.rowOf` the row), as the Sharing pane's **Share...** does. It
does not wait for the dialog, so a caller such as the AI assistant is not held while the user fills it; a refusal of
`shareRow` before the dialog opens is shown with `reportError`. No headless share exists: a grant on a row from JS is
refused (see "Sharing").

## Grokky knowledge (agents/)

The Datagrok AI assistant (the Grokky package) runs in a container and reaches Forge only through its `datagrok_exec`
tool: JavaScript run in the user's browser, where `grok.functions.call('Forge:<function>', ...)` runs the nine API
functions above. Its server-side function call (`/public/v1/functions/<name>/call`) cannot run a JS package
function, the same limit as `grok s functions run`; its function search reads the `description` of a function and of
its inputs, which is why every API annotation describes its inputs and, for an `object` output, the result's keys (the
generated `package-api.ts` types `object`, `map` and `list` as `any`).

- `agents/package-knowledge.yaml` (the shape of the other packages' files: `packageName`, `description`, `keywords`,
  `overview`, `apiRef`, `docsRef`) says what Forge does, where its UI is and which functions it has.
- `agents/skills.md` is the assistant's skill: a routing table from the user's intent to a function, every signature
  with its result shape and error sentences, `datagrok_exec` recipes and the rules (no upload unless a copy is asked
  for, apply only a model that fits, confirm a delete, `ML | Models` is not Forge). It is self-contained: the runtime
  reads a package's `agents/` in its clone of the public repository's `master` (`workspace/packages/Forge/agents/`)
  and, for an installed package that is not in the clone, syncs the published flat `agents/` files from
  `packages/published/files/<name>/<version>/<buildHash>/<buildNumber>/agents/` as `Forge-package-knowledge.yaml` and
  `Forge-skills.md`, where `apiRef` resolves nowhere.
- `grok publish` uploads `agents/` with the package (`.npmignore` does not exclude it); nothing in the manifest names
  it. No package test reads the published files.
- `ForgeApp.aiDescription`, the catalog's briefing for the assistant, names the API functions and `agents/skills.md`.

`skills.md` repeats the annotations of `package.ts` and the README's "From scripts and Datagrok AI": a change of a
function's name, inputs, result keys or error sentences updates it in the same change.

## Catalog logic

`catalog/`, one function per file, for the model catalog and the model's panels:

- `applicableTables(row, tables)` (`applicable-tables.ts`): the tables in which `isSuggested` finds a close column for
  every feature the row needs (`requiredFeaturesOf` of its `features` and `options`); none for a row without a
  feature list. It runs in the browser on the open tables; nothing is uploaded. `featureFit(row, tables)` gives the
  names of every feature (`featureSchemasOf`) with them (the card and Details). `modelsFor(rows, table)` keeps the
  rows that fit one table, in their order (the table pane).
- `modelActivity(id)` (`model-activity.ts`): the model's `application` records, newest first, at most
  `MAX_ACTIVITY_ROWS` (100), the full `count` (the list's length below the cap, a count query at it), and `lastRun`,
  the newest application of any status (`when`, `status`), undefined for a model never applied.
- `compareModels(rows)` (`compare-models.ts`): the table "Compare models", one row per model in the given order, with
  the columns Name, Description, Method, Task, Target, Training rows, Created, then `<metric> (train)` and `<metric>
  (validation)` (labels of `METRIC_LABELS`) for every metric of `METRIC_IDS` at least one of the models has; a metric
  a model lacks is empty (stored metrics are read with `metricsRecordOf` of `training/train-model.ts`); `metricColumns(rows,
  splits)` builds those columns, shared with `Forge:listModels` (validation only).
  `COMPARE_COLUMNS` are the `model` columns it reads. Fewer than two models: "Select at least two models to compare."
- `readModelFile(row)` (`model-file.ts`): the bytes of a model's blob in Forge's own storage (`ownBlob`), named
  `<model name>.bin` with `\ / : * ? " < > |` replaced by `_`; for any other blob: "The model 'X' has no model file in
  Forge's storage, so there is nothing to download."
- `updateModelInfo(id, version, {name?, description?, tags?})` (`model-edit.ts`): writes the given fields with
  `forgeDb.models.update`; with a `version`, a row changed since is refused with the platform's
  `DomainVersionConflictError`, without one the write wins. The tags (`string[]`) are written as their column text
  (`tagsText`; an empty list clears the column). It fires `modelsChanged` and returns the row's new version. The
  same file holds `runUpdateModel` (`Forge:updateModel`); `list-models.ts` and `get-model.ts` hold the bodies of
  `Forge:listModels` and `Forge:getModel`.

**Tags.** `model.tags` is a `string` column with the tags joined by ", " (`tagsText` in `storage/model-fields.ts`);
`tagsOf` splits the text on commas. Both go through `normalizedTags`: trimmed, without empty ones and repeats, in
their order. A model without tags leaves the column empty (the server stores the empty text as null). A comma typed
inside a tag therefore splits it into two tags. The code works with `string[]` (`modelFieldsOf`, `updateModelInfo`);
only these two helpers see the text. Why not a `string_list` column: on Datagrok 1.28 EMS refuses a list in a
single-row insert or update ("Types differ. Expected: string_list, passed: list": the batch loader special-cases
`string_list`, the single-row validation does not) and reads one back as the PostgreSQL array text (`{a,b}`). Once
the platform fixes that, the column can become a `string_list`. In the UI, **Tags** is the platform's chips input
(`ui.input.tags` with new items allowed, `tagsInput` in `ui/tags-input.ts`, its text box's placeholder `Type a tag and
press Enter`): a tag is typed and confirmed with Enter; `tagsOfInput` reads the list. While the box holds text, the
input keeps that Enter to itself (the dialog's Enter OK never sees it) and adds the chip asynchronously; a test types a
chip into the box (`typeTag`), since assigning `value` skips that path. The platform builds the box without
`autocomplete="off"` (its own text inputs set it, `text_input.dart:10`), so the browser offered its autofill history
there; a pick fills the box with no keystroke, makes no chip and saves nothing. `tagsInput` sets `autocomplete = 'off'`
on the box, as the platform's text inputs do.

## Model handler and panels

`ForgeModelHandler` (`ui/model-handler.ts`, a `DG.DomainObjectHandler` of `forge.model`, named "Forge model handler")
renders model rows wherever the platform shows them: the catalog's context panel, the Domains view, the Browse tree
and the **Predicted by** pane (which reads only the card's columns). The autostart function `_initForge` registers it
once (`ForgeModelHandler.registerOnce`), together with the model commands and the context panel's follower
(`followChanges`, below). Its members read only the row's values; a
row made from catalog values (`ForgeModelHandler.rowOf`, through the platform's `rowFrom`, with the system dates as
text) carries the catalog columns only. The Domains view's rows come without the json columns: the card and the
tooltip take the count from `feature_count` and read `features` and `options` once per model and session
(`withFeatures`), filling the feature list and the **Applicable to** icon in when they come.

- Icon: the model icon. Markup: the icon and the name. Tooltip: Forge's rich tooltip (`ui/rich-tooltip.ts`,
  `richTooltip`: a bold title, the main lines, the minor lines in grey, at most 320 px wide; `shortList` cuts a list
  after about two lines with `+N more`; the same file's `bulletList` fills a `ul` for the insights' tooltips and the metrics
  grid's bullets): the name; `Predicts <target> from <N> features with <method> · <task> · <rows>
  rows` (the parts the row has: no `from` part without features, no `with` part without a method, no `Predicts` part
  without a target); the features; in grey `Created <YYYY-MM-DD> by <login> · Data storage: <mode> · Tags: <tags>` (the parts the
  row has; the author's login is read through `userOf` (`ui/model-panes.ts`: `grok.dapi.users.find` once per user and session, a failed read retried; shared with Details and the table pane's tooltip) and filled in when it comes).
  Grid header and cell tooltips stay single short lines.
- Card, the built-in tool's in five lines of one height (`forge-model-card`: 90 px of content, each line 18 px, cut
  with an ellipsis; `content-box`, because `ui.card` gives a gallery card 10 px of padding and `div.ui-div` counts it
  into the height): the name (`grok-gallery-grid-item-title`) with a table icon when an open table fits the model (its
  tooltip `Applicable to <tables>`), `Predict <target>`, `by <first five features>...` (`by <N> features` until the
  list is read), `using <method>`, `Created on <YYYY-MM-DD>`; a missing part leaves its line empty. A double-click is the platform's Open (the row's own view).
- `renderProperties` reads the full row once and returns `modelAccordion` (`ui/model-panes.ts`), whose title is the
  one the platform gives a row: the model icon, the favorites star (`ui.star`), the name and the row's commands
  (`ui.contextActions`), added before the panes (`Accordion.addTitle` appends); a row the user can no longer read shows
  "The model is no longer available.". It replaces the platform's whole panel of the row (generic
  Details, Shared with, Chats, History, Actions): the platform has no JS surface for its Chats, and its History is
  hosted as a pane. The accordion's panes are built when first opened; Details and Activity share one
  `modelActivity` query. The newest accordion follows `applicationRecorded` (`apply/application-store.ts`, an
  in-memory Subject fired by `recordApplication` in this session, after a successful insert, whatever the record's
  status; another user's or another tab's application is seen only by a panel built afresh): while the accordion is in
  the page, for its model it drops the query and shows again Details' **Last run** and **Applications** cells and the
  Activity pane if built (a pane not built yet reads afresh when opened). Details' other rows and its Tags input (with
  the version of its last write) stay, as do the current object, the expanded panes and the other panes; an evaluation
  records nothing. The roots of
  Performance and Activity carry `forge-pane-fill` (`width: 100%`, `css/forge.css`): the accordion's pane content sizes
  its children to their content, where a grid's own `width: 100%` resolves to the grid's width, so without it the
  grids stay narrower than the panel; with it they span the panel and follow its resizing:


| Pane | Content |
|---|---|
| **Details** (open) | The description, then **Author**, **Created**, **Updated**, **Table** (`<table> (<rows> rows)`), **Data storage** (`None` / `Reference` / `Copy`), **Data source** (a reference: its `path`, else its query `name`, else `a script`), **Data copy** (a copy: `ui.render` of the uploaded `TableInfo`, `missing` when `grok.dapi.tables.find` gives none), **Last run** (the newest application's time, with ` (failed)` / ` (cancelled)` when it did not complete; `Never` when none), **Applications** (the count), **Features**, **Target**, **Method**, **Task**, **Applicable to** (the open tables it fits; left out when none), then **Tags**: every change is written with `writeModelInfo` (`ui/save-model-dialog.ts`, shared with Edit model; tags only), one write at a time, each against the version the previous one returned, so quick changes are no conflict; a model changed elsewhere meanwhile asks with the platform's conflict dialog to reload (Details is rebuilt from a fresh read) or overwrite |
| **Performance** | `MetricsTable` of the stored metrics (`storedSummary`; the Train view's grid: Metric / Train / Validation, the findings marked in red, and its bullet list of rows, validation, seed with the copy icon, positive class); "No metrics were recorded for this model." without metrics; then the button **Run Evaluation** (`runEvaluation`, see "Run Evaluation"; tooltip `Apply the model to its training data, or to an open table, and show how well it predicts.`) when the model's blob is Forge's own (`ownBlob`), so a row made by hand has none |
| **Activity** | `N application(s)` and a grid of the newest applications (up to 100): **When**, **Who** (login, the authors read with one `grok.dapi.getEntities`, `loginsOf`; empty for a deleted user), **Table**, **Rows**, **Prediction column**, **Status**, **Source**, **Duration (ms)**, with a tooltip per header, the time with seconds on **When**, the meaning on **Status** (`Completed: the column was added.`, `Failed: <error>`, `Cancelled: stopped between batches, no column added.`) and on **Source** (`The Apply dialog or the catalog.`, `A script, through Forge:applyModel.`); "Not applied yet." without one |
| **Sharing** | The groups the model was shared with (**Can view**, **Can edit**) or "Not shared yet. Only its author and administrators can see this model.", the button **Share...** (`DG.DomainObjectHandler.shareRow`: the platform's sharing dialog; its one `try/catch` reports a refusal before the dialog opens), and "Sharing a model needs the Share permission on this model; ask an administrator." for a user without Share on the row (`~can_share`; only for such a user is a refused read of the shares expected, and it leaves the groups out) |
| **History** | The platform's audit pane (`DG.DomainObjectHandler.auditPane`): one line per change of the row (user, time, `insert` / `promote` / `update`), newest first |

The Sharing pane reads the shares with `grok.dapi.permissions.get`, which reads the grants of the entity's wrapping
project, where **Share...** writes them. The author's own View, Edit, Delete and Share grants, written on the row by
eager promotion, are not visible to it, so a model nobody shared it with reads "Not shared yet". `refreshSharing`
fills the pane, when it opens and again on every `grok.events.onEntityShared` for the model (the sharing dialog fires
it after OK; `shareRow` itself resolves when the dialog opens). The pane resolves the model's entity once; each fill
reads the shares and `~can_share` (a one-column query with access) side by side. There is one subscription: each new
Sharing pane replaces the previous one's (the context panel shows one model). The read is not stale within a
session: a share made elsewhere shows on the next read.

The context panel follows the model it shows (`ForgeModelHandler.followChanges`, subscribed once per page by
`_initForge` with the handler's registration): 300 ms after the last `modelsChanged` of a burst, whoever wrote
(**Edit model**, the catalog's delete, Details' Tags, `Forge:updateModel`, `Forge:deleteModel`), it reads that model
again (`CARD_COLUMNS`). A model gone (deleted, or no longer visible) empties the panel (`setCurrentObject(null)`, as
the catalog does); a model whose `version` differs from the one the newest accordion shows (`isShownVersion` in
`ui/model-panes.ts`; Details' own Tags writes move it, so typing chips rebuilds nothing) is set again, forced, from the
fresh row, and the title and every pane are built anew; a change of another model costs one read and rebuilds
nothing. A panel that moved on during the read is left alone. `modelsChanged` lives per bundle: the package's
follower does not hear the test bundle's writes, so a test of a UI path that writes from the test bundle starts its
own `followChanges` and unsubscribes it.

## Model commands

`ui/model-actions.ts` has two command sets. `MODEL_ACTIONS` (**Apply...**, **Download**) are registered once by
`registerModelActions` as param funcs of `forge.model`: they show wherever the platform shows a model (the Domains
view, the Browse tree, the panel title's `ui.contextActions`), next to the platform's own Open, Edit..., Clone,
Delete, Share..., History, Copy link and Watch (no JS API lists param funcs). `CATALOG_ACTIONS` (**Apply...**, **Edit
model...**, **Download**, **Delete model**) are the catalog grid's menu (`addModelItems`): the grid's menu carries a
`GridCell`, not the row, so the platform adds nothing there.

| Command | Set | Shown | Does |
|---|---|---|---|
| **Apply...** | both | always | The Apply dialog for the model (`modelId`), on the current table or the first open one; from the catalog with **Applicable to** as `preferredTable` |
| **Edit model...** | catalog | always | **Edit model** (`ui/edit-model-dialog.ts`): Name (not blank), Description and Tags, prefilled from a fresh read of those columns; OK writes them with `writeModelInfo` against the read version, a model changed meanwhile asks to reload (the dialog reopens) or overwrite; then the balloon `Model "<name>" updated.`; a context panel showing the model follows the write (`followChanges`, see "Model handler and panels") |
| **Download** | both | the model's file is Forge's own (`ownBlob`) | `readModelFile`, saved by the browser as `<model name>.bin` |
| **Delete model** | catalog | always | `confirmDeleteModel`: "Delete the model "<name>"?", then `deleteModel` (row, file, applications) |

The platform's generic **Edit...** changes every column and its **Delete** removes the row only, leaving the file; a
JS handler cannot remove them. Every command is a UI boundary (`reportError`).

## Model comparison

`ui/model-comparison.ts`. A `ModelComparison` holds the compared rows (`CompareModelRow`, read with
`COMPARE_COLUMNS`); `ModelComparisonHandler` (type `forge.model.comparison`, "Forge model comparison handler")
claims it by its `kind` (not its class: the test bundle has its own copy), captions it `Compare N models` and renders
`comparisonForms`. The context panel shows no caption of its own, so `comparisonForms` starts with the title
`Compare N models` and the model icon, as the platform's accordion title (`d4-accordion-title`, the look of the model
panel's title); then PowerGrid's viewer `DG.Viewer.fromType('Forms', ...)` over `compareModels` with every row selected,
the fields `compareFormFields` (Name, Method, Task, Target, Training rows, Created and the `(validation)` metric
columns; the viewer shows at most 20), at an inline width of 100% (the panel's width) and height of 160 px per model up
to four. The viewer lays its forms out once, when attached: the root follows a later panel resize, the forms inside do
not (the library viewer has no resize handling). Without PowerGrid (`hasFormsViewer`: no `PowerGrid:formsViewer`)
the panel shows `Install the PowerGrid package to see the comparison as forms.` and the grid. `openComparisonView` is
the Compare view: the table view of `compareModels` with the same Forms viewer added.

## Predicted by

`isPredictionColumn(col)` (`Forge:isPredictionColumn`) is true when the column carries the `forge.model` tag
(`PREDICTION_TAG`). The panel function `predictedByPanel` (friendly name **Predicted by**, `PANEL_PREDICTED_BY`;
`meta.role: panel`, condition `Forge:isPredictionColumn(col)`) shows the card (with its rich tooltip; it reads the
card's and the tooltip's columns) of the model whose id the tag holds, or
"The model that predicted this column is no longer available to you.". The raw tag stays in the column's tags.

## Table pane

The panel function `tableModelsPanel` (friendly name **Forge models**, `PANEL_TABLE_MODELS`; `meta.role: panel`, one
`dataframe` input, no condition, so the platform lists it for every table and builds it only when the pane is
expanded) returns `tableModelsPane(table)` (`ui/table-models-panel.ts`): under `ui.wait`, it reads up to `MAX_MODELS`
rows with the card's columns (`CARD_COLUMNS` of `ui/model-handler.ts`, shared with **Predicted by**) and
`dataset_name`, newest first,
keeps those `modelsFor` the table (in the browser; nothing is uploaded) and shows one line per model
(`forge-model-line`), the built-in pane's look (`table_meta.dart` `renderModelsSection`): the run icon
(`ui.iconSvg('model-run')` with `d4-entity-markup-icon`, 12 px in the built-in row `span.d4-entity-markup-row >
div[data-wrapper=true]` of `d4.css`, inside a `d4-flex-col` list, centred on the name's line by
`forge-model-line`; **Run prediction**: `openApplyDialog(table, {modelId,
preferredTable: table, switchToTable: true})`), the name as a link (`d4-link-label`, no model icon) and the menu arrow
(`ui.contextActions(row)`: the platform's row commands and Forge's **Apply...** / **Download**). The name is not
`ui.bind`-ed, since a bound JS row gets the generic menu (only **Properties...**) and Forge's rich tooltip; instead a
click sets the current object (forced: for 2 s after an accordion pane is expanded the platform drops an unforced
change), a right-click opens the arrow's menu, `ui.makeDraggable` drags the row (onto a Space),
and its tooltip is the built-in model tooltip (`modelSummary`: the name, then Created by, Created, Updated, Trained on,
Last run, Tags, Inputs, Outputs, Method and Applicable to when an open table fits; the built-in **Links** row has no
Forge data and is left out). The tooltip is built from the pane's row with `modelFacts` (`ui/model-panes.ts`, shared
with Details); its only reads are the author, through the session's cached user lookup (`userOf`, shared with the
model's rich tooltip), and the last run, read once per line (a failed read is not kept). It is built inside `ui.wait`
(the platform builds a tooltip in a timer and silently drops it when the builder throws synchronously; `ui.wait`
shows the error instead). Without a fitting model the text `No model fits this table.`; then
the big button **Train** (`Train a new model on this table.`, `TrainView.openFor(table)`). The pane holds no subscription: the platform never detaches a panel's widget, so a
subscription would leak; the pane is built again whenever the context panel shows the table again, as the built-in
**Models** pane is.

## Sharing

Every model and training run is a platform entity from its insert, and its author holds View, Edit, Delete and Share
on it, so the author sees, applies, edits and deletes his models whatever the table's grants. Sharing goes through
the platform: **Share...** opens the sharing dialog (the row is already promoted), which writes the grant on the
entity's wrapping project. A non-administrator author shares his own model too: eager promotion gives him Share, so
GROK-21041 no longer blocks own models (WO-3, `forge.tester` on the release build). A user a model is shared with
sees it in the catalog, applies it, and the application is recorded (`application` is secured by its model). He
trains, saves and deletes models of his own. Linking a model into a shared Space is another way to share it.

**Spaces.** A model reaches a Space only by the platform's drag-and-drop: a `forge.model` row dropped on a Space node
goes to `DomainSpaceLinks.drop` (`xamgle/.../meta/domain_ui.dart`), which confirms, promotes the row (already eager)
and saves the `ProjectRelation` itself, reusing an existing one, so a second drop links nothing twice; the Space's
members then gain View of the model. Forge's model elements are drag sources: the Domains card (`renderCard`) and
the context-panel title's name are bound with `ui.bind`, the table pane's name is a `ui.makeDraggable` source with the
row as the drag object; either way the drag object is the row's Dart `DomainRow`, the type the drop accepts. No JS surface links a row: `Project.addLink` and
`grok.dapi.spaces.id(x).addEntity(id, true)` both pass `Project.canLink`, which accepts packages, files, domain schemas,
tables, saved filters and func calls but no row, and the server refuses with "Cannot add this entity to the space"
(`spaces_service.dart`); a regular project cannot hold a row at all. So Forge has no **Add to** command.

Observed on the local stand (image 1.28.0): an administrator sees only the `model` and `training_run` rows granted to
him, not every user's (platform observation P3: this checkout's EMS has an admin bypass for row predicates,
`predicate_builder.dart:140-141`, which the stand does not apply), and training runs are not shared with their model.
From JS, `grok.dapi.permissions.grant` on a model row is refused in the package tests ("You don't have a permission
to share this object"), so the tests share nothing; `grok s shares add` (the public API) shares a row.

## What differs from the built-in tool

The form starts ready to train: the last column is the target and the numerical columns except row numbers and ids
are the features, where the built-in tool starts with nothing selected. Problems are shown live on the input they
concern and nothing trains until they are fixed. The live retraining and the 200 ms delay are the built-in tool's,
but a method chosen by the user stays chosen while it applies (the built-in tool reset it on every data change), its
hyperparameter values are kept per method, an empty method list is a red **Method** instead of a line of text, the
ribbon holds only **Save** (no TRAIN / SAVE switch of one button; a **Train** button under the inputs appears only
for a method that does not retrain on every change), and the task bar names the method.
Validation uses five folds that cover every row once,
from a seed saved with the model, instead of five overlapping random samples, and it runs from 10 rows instead of
being skipped below 100. MAE and F1 are new; AUC-ROC is measured on the Predict probability scores, whose target is a
float 0/1 column (the built-in tool's int target gets rounded 0/1 predictions from XGBoost). One-hot records the
training categories and replays them (the built-in tool encoded the applied table's categories), and text and yes/no
features are refused while One-hot is off (which starts checked when every such feature has at most 20 categories).
Skip unique categories records the columns it dropped, and applying neither asks for them nor passes them to the
method (the built-in tool re-ran its rule on the applied table). AUC-ROC takes tied scores as one diagonal
segment, as scikit-learn (the built-in ROC curve stepped row by row, so the order of tied rows changed it; the user's
decision). Training and saving are separate steps, and
every completed, failed or cancelled attempt is recorded. The training table is uploaded only when the user chooses
**Copy** (the built-in tool uploaded it on every save); otherwise a model keeps a fingerprint of it, and a
**Reference** to its source when the platform recorded one. The DONE link of the task bar became the **Apply...** link
of the balloon after Save. Scripts train, save and delete models with `Forge:trainModel` and `Forge:deleteModel`
(the built-in tool trained on the server and offered no JS function for it). Rows with a missing
target are skipped instead of blocking training, and dates and bigint columns are refused as features and targets.
Missing values are a **Missing values** choice (Skip rows / Impute) with the imputation settings inline, instead of
two checkboxes and EDA's own dialog, and the target is never imputed.

Applying lists every visible model, up to 10000 (those that fit first), instead of only the suggested ones, unfolds
the feature rows in place under a **Columns** summary instead of a sub-dialog, checks every mapped column's kind
before the engine (the API too, which used to check presence only), suggests close names only within a distance
threshold instead of always assigning a column, makes every check a validator with **OK** unavailable until all pass,
skips or imputes rows with missing values at application, checks for a cancel before every batch, records every
application, and names and tags the prediction column with Forge's own `<target> (predicted)` and `forge.model`.

The catalog stays a grid (with **Applicable to** filtering open tables in the browser instead of uploading one), the
card shows in the Domains view. The model's activity comes from its application records, not from log events, and
"Last run" from the newest of them. Performance shows the stored metrics, with the values below a threshold or far
from the training one in red as in the Train view (the built-in tool marked them only in its preview). Run Evaluation
opens the model's reference or copy, or asks for an open table, shows the stored metrics with an **Evaluation**
column and the plots in a view of its own, and records nothing. The card has five lines of one height, **Applicable to** an icon on its first line, and the
model's tooltip is a compact rich one with the full information (the built-in cards had three to six lines). The
results are rows of panels in the platform's splitters, with initial heights of 350 / 300 / 450 / 350 px, initial
widths by weight and resizable (the built-in previews were cards of 300 to 600 px wide and 350 px high that wrapped),
kept with their viewers while a retrain keeps the data (the built-in previews were built again), with the
**Performance** title over the metrics panel, a tooltip on every title and a **Class balance** bar chart of the
target's classes (the built-in tool's stacked target and prediction chart was never shown); an insight about a panel
points to it and is that viewer's description, the unbalanced classes are one line with the shares, and the deviating
prediction statistics are one line instead of one per statistic (the built-in tool listed each); the threshold is tested on both splits (the built-in tool tested the train value), and **Wrong predictions** shows the
table's own columns and the prediction (the built-in tool cloned the target and the prediction). The panel
has no Chats pane. Two or more selected models show in the panel as forms (the built-in tool compared them only as a
command). Tags are a column of the model, not entity tags. The model file downloads as `<name>.bin` instead
of a zip.

The table's pane is **Forge models** (the built-in tool's pane is **Models**): it lists the models that fit the table,
found in the browser by close column names and kinds, instead of the server's suggested models; each line's run
icon (`model-run`) is **Run prediction** as in the built-in pane, the model's menu is the platform's row menu with Forge's
commands (the built-in menu had **Run** and the model's commands), the tooltip has no **Links** row, and the empty
list says `No model fits this table.` above **Train**. The catalog view's `aiDescription` is Forge's own text.

The view's history is the **Runs** grid instead of a card per training: one row per training that was not superseded (a cancelled or
failed one too, also a failure before any run was recorded), every row side by side, so there is no pinning; a click
on a row whose inputs differ from the form's sets them again and trains (a click on a row the form already matches
and shows the completed result of only makes it current; a cutoff-only difference re-cuts), as a click on a pinned
card reloaded its parameters.

The model's menu has no **Add to** (the built-in tool added a model to an open project): a model reaches a Space by
the platform's drag (see Sharing). Forge's Browse node is **Apps > Forge**; **Browse > Platform > Predictive models**
and its **MLFlow** node stay the built-in tool's until the switchover. The Domains view of `model` can be filtered by
**Tags** (the built-in catalog's tags filter); its **Author** filter has no counterpart.

## Automation surface

The views publish their state in the DOM, so a test waits on what a view says instead of on a timer:

- The Train view's results (`.forge-train-results`, `TrainView.publishState`, run by every `updateControls` and when a
  `ResultsPane.show` settles): `aria-busy="true"` while a check is pending (`check` undefined), a training is queued,
  running or being recorded (`isTraining`) or results are being shown (`showing`, the `ResultsPane.show` calls not
  settled); `aria-invalid="true"` once checked while the selection has a problem (`problem()`, an invalid setting
  included) or its last training ended `failed` or `cancelled` (`isFailed`, cleared by any change but the cutoff and by
  a training of the current selection that completes, a **Train** click after a cancel too). Ready is
  `aria-busy="false"` without `aria-invalid="true"`; a slow method waiting for **Train** is ready.
- On the same element, `data-forge-run-id` (`SHOWN_RUN`) is the id of the `training_run` row whose result the results
  show. `showResults` sets it when a completed training of the current selection is shown, and a cutoff-only re-cut
  sets the same id again. `showHint` removes it: before the first training, for a selection problem or an invalid
  setting, for a method waiting for **Train**, after a failed or a cancelled training, when a **Train** click starts
  one, and on a Table change. A training that runs while a result is shown (a live retrain) leaves the old id until
  its own result is shown, and a superseded training leaves it as it is: it records no run and shows nothing. A new
  training is thus a ready view with an id other than the one before; a "Ready" alone is not one (a method waiting
  for **Train** is ready, and the metrics grid is updated in place).
- The Evaluation view's results (`.forge-train-results` of `EvaluationView`): `aria-busy="true"` from `open` until the
  results are shown or the evaluation ends; a failure or a cancel closes the view, so it has no `aria-invalid`. It has
  no `data-forge-run-id`: an evaluation records no `training_run` (`training/evaluation.ts`), and the view
  shows the one result it opened for.
- The Forge view (`ForgeApp`): its constructor adds the class `forge-app` to the view root, so a claim can be scoped to
  Forge's own view. Its models pane, `.forge-catalog`, publishes the catalog's loads (`ForgeApp.publishLoad`): the
  constructor sets `data-forge-catalog-load="0"` and `aria-busy="false"` once the first catalog is bound. A
  `modelsChanged` and a **Refresh** click each make a new load due at once (`markStale`: the next load number, per view,
  and `aria-busy="true"`), a change before the 300 ms debounce of its reload; a load still running is then superseded and
  publishes nothing, so a load read in between is never taken for the change's. The `refresh()` that follows loads
  under that number and, once the reloaded frame is bound and the current row and the selection are found again, sets
  `data-forge-catalog-load` to it and `aria-busy="false"`; a failed load sets only `aria-busy="false"`, and only while
  its number is the newest. Each change or click advances the number by one. A reload is thus a load number above the
  one remembered with `aria-busy` not true; every open catalog, hidden ones too, reloads on `modelsChanged`.
- `CollapsibleGroup`: `aria-expanded` on `.forge-group-header` follows `setExpanded`; `aria-invalid` on
  `.forge-group-summary` follows `setSummary` (the Apply dialog's red `N of M matched`).
- Every table of data is the platform's grid, whose widget status publishes the hit areas `cell <row> of <column>`,
  `header <column>` and `row header <row>` (table rows counted from 1) and the readings `rows`, `rows shown`, `columns
  shown`, `current row`, `current column` and `text of cell <row> of <column>`.

The user flows are Gherkin features on `@datagrok-libraries/bdd` in `bdd/` (`bdd/README.md`), run by Playwright
against a stand. Forge's vocabulary is in `bdd/bindings/`: the kinds `group` and `results panel` over the classes
above, `detail` (a caption-value row of a `ui.tableFromMap`, in the context panel or a tooltip) and `model line` (a line of the table
pane **Forge models**: its run icon, name and menu arrow), a named table set as the current object (the platform's
corner click does not make it current again), forced and not frozen, the catalog's loads over the
`.forge-catalog` of every open Forge view (a hidden view's root is read through `grok.shell.views`), the captions a grid
draws (the text of the header cell at each column's `header <column>` area), the demo-file opener (`OpenServerFile`
unprocessed, so the table keeps its creation script), the API fixtures (`Forge:trainModel`, a row inserted without a
model file), and the server claims, which run the `grok s` CLI from Node as the page's user (a guard compares the two
logins). A feature that saves data names it
`forge-test-<purpose>-{time}` and sweeps it at its start and end: its models through `Forge:deleteModel` in the
page (with their folders and copies), its tables' training runs through `grok s domains delete`, each read back gone
by the exact names and ids it deleted. A step that claims a training waits for a new `data-forge-run-id` and asks the
server whether that id is the table's newest run, with the method and the settings the step chose; a claim that
nothing trained reads once the view is settled (`aria-busy` not true).

## EMS schema `forge`

Manifest: `databases/forge/schema.json`, version `0.1.14`. `grok publish` deploys it; a debug publish applies
destructive changes without migration scripts, except a `promotion` change on a row table, which is always refused
(`[promotion-change]`). To change it, publish a manifest with only a placeholder table, then the real one, bumping
`version` both times; the first publish drops the tables with their data and their rows' entities and permissions,
but no files (the blob folders of the dropped models stay). `grok api` generates the typed client
`src/generated/db.ts` (`forgeDb.models`, `forgeDb.trainingRuns`, `forgeDb.applications`). Every table also has the
system columns `id`, `version`, `created_on`, `updated_on` and `author_id`; who and when always come from them.

| Table | Security | Why |
|---|---|---|
| `model` | `row`, `promotion: eager`, `defaultRowVisibility: none`, grants `All users: view, edit` | A model is private until shared, like the old models. Eager promotion makes every inserted model a platform entity at once, with View, Edit, Delete and Share for its author's personal group: under visibility `none` a user sees only rows he holds a grant on (administrators see all by design; the 1.28.0 stand does not apply that, see Sharing), so without it a user would not see his own models. Sharing, favorites and comments work on every model. The table Edit grant lets any user insert models; with visibility `none` it reveals no foreign rows |
| `training_run` | `row`, eager promotion, `defaultRowVisibility: none`, same grants | The trainer's experiment history, including failed and cancelled runs, visible to the trainer for the same reason. `model_id` is optional with `onDelete: setnull`, so the history outlives the model. `source` (0.1.14) is `ui` (the Train view, the default for the rows recorded before it) or `api` (`Forge:trainModel`); inserts pass it explicitly |
| `application` | `master`, `delegate: model_id`, `audit: false`, no grants | Secured by the model (see Key decisions). Master tables refuse grants. No audit: high-churn records. `onDelete: cascade` from the model. `status` is `completed`, `failed` or `cancelled`; `skipped_rows` counts the rows left without a prediction because of missing values. Inserts pass `source` and `status` explicitly |

Deletes are soft: a deleted row stays with `is_deleted` and disappears from queries and counts. Current EMS limitation:
a soft delete does not clean up a promoted row's entity and permissions, so every deleted model (and every model a
test run saves and deletes) leaves its entity until the platform fixes it or the next schema reset drops it; no
client API removes it. Since the stand's 2026-10-20 build (1.28.0) the entities API no longer returns a newly deleted
model's entity, while `grok s shares list` still lists its grants. The generic Domains view offers **New model...** to anyone with the
table's Edit grant; such a row has no engine blob, and Forge lists it and can delete it (the row only, unless its
file sits in Forge's own `forge/model/<uuid>/` layout).

### JSON column shapes

EMS `json` columns hold objects only; a list is wrapped in an object. TS types: `training/train-model.ts`,
`preparation/preparation-options.ts` (`PreparationOptions`), `engines/engine.ts` (`Hyperparameters`),
`storage/dataset-ref.ts` (`DatasetRef`) and `storage/dataset-fingerprint.ts`.

| Column | TS type | Shape |
|---|---|---|
| `model.target` (caption "Target details") | `TargetSchema` | `{name, type, semType?, categories?}`; `categories` for classification targets |
| `model.features`, `training_run.features` | `FeaturesSchema` | `{columns: [{name, type, semType?}]}`, in training order |
| `options` | `PreparationOptions` | Preparation replay, see below and "Preparation"; Forge writes `{preprocessingInfo, postprocessingInfo, missingValues: {mode, neighbors?, distance?, skippedRows}, skippedColumns?: [names], oneHotCategories?: {<column>: [categories]}}`, plus `positiveClass`, `negativeClass`, `binaryClassificationThreshold`, `targetType` with Predict probability |
| `hyperparameters` | `Hyperparameters` | `{<train function input>: value}` |
| `metrics` | `MetricsRecord` | `{train: {<metric id>: number}, validation: {<metric id>: number}, positiveClass?}`; ids `mse`, `rmse`, `mae`, `r2`, `accuracy`, `f1`, `sensitivity`, `specificity`, `precision`, `npv`, `auc` (Predict probability only) |
| `splitting` | `Splitting` | `{scheme: none \| kfold \| holdout, folds?, trainFraction?, isStratified?}` |
| `model.dataset_ref` | `DatasetRef` | Reference mode: `{kind: file \| query \| script, script, path?, id?, name?}`, see "Storage modes" |
| `dataset_fingerprint` | `DatasetFingerprint` | `{rowCount, columnCount, hash, columns: [{name, type, missingCount, min?, max?, mean?, categories?}]}` |

`options` keeps the keys of the platform's built-in models, so their preparation replays unchanged:
`preprocessingInfo`, `postprocessingInfo`, `positiveClass`, `negativeClass`, `binaryClassificationThreshold`,
`targetType`, `allowNulls` (`allowNulls` is read and ignored). Forge adds `missingValues`, `skippedColumns` and `oneHotCategories`; a
model without `oneHotCategories` replays one-hot, and one without `skippedColumns` skip-unique, as the built-in tool did. `PreparationOptions` and
`preparationOptionsOf`, which narrows a stored value to it, are in `preparation/preparation-options.ts`.

### Model blob

The `model.blob` column (type `file`) stores `file://System:DomainFiles/forge/model/<uuid>/model.bin`, the layout the
EMS row editor uses for file columns. The random path segment makes the file unguessable, and the row security of
`model` protects the pointer. The generic Domains UI shows the column but offers no download of the file.

### Storage modes

`model.storage_mode` (default `none`) records what a model keeps of its training data (`ModelStorage` in
`storage/model-fields.ts`), chosen in the Save dialog (see "Saving"):

- `none`: only `dataset_fingerprint`, enough to check that a new table fits the model.
- `reference`: `dataset_ref` (`DatasetRef`) points to the source; the data is not uploaded.
- `copy`: `dataset_table_id` is the id of an uploaded copy of the feature and target columns (no foreign key).

**Reference** (`storage/dataset-ref.ts`). `datasetRefOf(table)` reads the origin the platform recorded in the
table's tags: the creation script (`.script`, `DG.Tags.CreationScript`; one call per line, the first assigning the
table's variable, each line ending in a `//{"timestamp": ...}` comment, which is dropped). With `.DataQuery.id` it is a
`query` reference (`id`, `name` from `DataQuery.name`); a first line `OpenFile("<path>")` (or `OpenServerFile`) is a
`file` reference with `path`; any other script is a `script` reference. Without a creation script, a `source.file`
tag holding a server path (it has a `:`; a local file's tag holds only the file name) gives a `file` reference with
the script `data = OpenFile("<path>")` (the path escaped as a string literal). Anything else gives null: no reference
to offer. Which tables carry the
script: the platform records it after a function run, with data history on (the default), that is unprocessed or
runs without default result handling (`shell.dart`). **Browse > Files** opens a file that way (an `OpenServerFile`
call without default result handling, `file_editors.dart`), recorded as `<Name> = OpenFile("<path>")`; the tests
reproduce it with an unprocessed `OpenServerFile` call (`openIrisFromFile`). `grok.data.files.openTable`,
`grok.functions.call` and `grok.functions.eval` record nothing, and neither does a table built in code.
`openDatasetRef(ref)` opens a `file` reference with `grok.data.files.openTable(path)`. A `query` reference needs its
query (`grok.dapi.queries.find(id)`; a missing one is a `ForgeError`) and is then replayed like a `script` reference,
since its parameter values live only in the script: as the platform's data sync does, each line is evaluated in a
fresh `DG.Context` (`grok.functions.eval(line, context)`), and the table is the first line's variable. A stored
reference can be edited in the generic row editor, so a script is run only when every line is `<variable> =
<call>(<arguments>)` (an output accessor allowed) whose arguments, quoted strings aside, hold no call and no `;`;
anything else is refused with a `ForgeError`. The table is not added to the workspace. `storedDatasetRef(value)` narrows a stored `dataset_ref` to a `DatasetRef` (null when it is not one).

**Copy** (`storage/dataset-copy.ts`). `uploadTrainingCopy(columns, modelName)` uploads a shared frame of the given
columns (every row) with `grok.dapi.tables.uploadDataFrame` and returns the table id; the frame is named
`trainingCopyName(modelName)`, `<model name> (training data)`. The server keeps that name as the table's
`friendlyName` and makes `name` an identifier (`ForgeTestCopy...TrainingData`), so `grok s tables list --filter
"training data"` finds the copies. `deleteTrainingCopy(id)` deletes the table only when its `friendlyName` ends with
` (training data)`, so a `dataset_table_id` edited in the generic row editor cannot delete another table; a table
already gone is no error (`grok.dapi.tables.find` resolves undefined).

`has_training_rows` marks model files that embed training rows: `hasTrainingRows(engine)`, the method's
`meta.mlhasrows: true` on `train` (EDA's SVM keeps its support vectors). `legacy_id` is reserved for migrated models.

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
| `splitting` | Data split | `tags` (model) | Tags |

`application`: `model_id` Model, `table_name` Table, `row_count` Rows, `skipped_rows` Skipped rows, `column_name`
Prediction column, `source` Source, `status` Status, `error` Error, `duration_ms` Duration (ms).

The Domains view of `model` offers the filters **Method** (`engine_name`), **Task**, **Data storage**, **Tags** (the
comma-separated text) and **Created** (a range). A `training_run` json column carries the same `description` as the
`model` column of the same name (the schema test compares them). The catalog grid reads the same captions with
`grok.dapi.domains.registry.rowProperties('forge.model')` in `ForgeApp.loadModels()` and sets
`column.meta.friendlyName`; the system column `created_on` is labelled "Created" explicitly. The grid shows Name,
Method, Task, Target, Data storage, Training rows, Tags and Created, in this order (`grid.columns.setOrder` and
`setVisible`); `id`, `version`, `updated_on`, `author_id` and the loaded `features`, `options` and `blob` (for the
card, **Applicable to** and **Download**) stay hidden. "No models yet." is shown while the catalog is empty.

## Engine contract (methods in the UI)

An engine is the set of functions that share `meta.mlname` (the engine id). `meta.mlrole` gives each function's role:

| Role | Signature | Notes |
|---|---|---|
| `train` | `(df, predictColumn, ...hyperparameters) -> model blob` | The other inputs are the hyperparameters |
| `apply` | `(df, model) -> dataframe` | Forge takes the first column as the prediction, out of the returned frame |
| `isApplicable` | `(df, predictColumn) -> bool` | |
| `isInteractive` | `(df, predictColumn) -> bool` | Optional. `meta.mlupdate: 'false'` turns off live retraining |
| `visualize` | `(df, targetColumn, predictColumn, model)` | Optional |

Three more `meta` keys on `train`: `mlhasrows: 'true'` (the model file embeds training rows; EDA's SVM and
k-Nearest Neighbors), `mlserver: 'true'` (the method runs on the server) and `mlnativenans: 'true'` (the method takes
rows with missing feature values as they are; EDA's XGBoost, Random Forest, Extra Trees, Hist Gradient Boosting,
CatBoost, Decision Tree and TabM; see "Missing values").

- An engine is complete, and usable, only with `train`, `apply` and `isApplicable`.
- A function engine receives the feature table and the target column. A script engine (`DG.Script`) receives one
  table with a copy of the target appended, and the target name.
- The calls (`isApplicable`, `isInteractive`, `train`, `apply` in `engines/engine-calls.ts`) take the feature columns
  and build the table themselves, a frame of the columns given back by `onFrame` when the call ends.
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

### Choosing a method

- `applicableEngines(engines, features, target)` (`engines/applicable-engines.ts`) asks every complete engine's
  `isApplicable` at once (`Promise.allSettled`) on the feature columns (each call with its own frame). `applicable` keeps the given order (discovery order: functions first); an engine whose check throws is
  left out and returned in `failed` with its error, for the caller to log (nothing is logged here). On iris, Species
  gets every EDA method but Linear Regression and PLS Regression (14); Petal.Length every one but Softmax, Naive Bayes
  and LDA (13).
- `selectBestEngine(engines, features, target)` (`engines/best-engine.ts`) is the built-in tool's suggestion, ported
  as is: a feature with the `Molecule` semantic type gives `Chemprop`; a classification (the target is not numerical)
  with at least one categorical and one numerical feature gives `XGBoost`; a regression with five or more numerical
  features gives `PLS Regression`, any other regression `Linear Regression`; everything else `XGBoost`. "Numerical"
  and "categorical" are the platform's `isNumerical` / `isCategorical` (dates count as numerical, as in the built-in
  tool); a feature named like the target is not counted. A suggestion not in the list gives the first engine, an
  empty list undefined. The built-in tool matched `<package>: <name>`; Forge's registry keeps one engine per name, so
  it matches the name.
- `isServerEngine(engine)` (`engines/engine.ts`): the training data leaves the browser (no UI uses it yet; `Forge:listMethods` reports it as **Runs on**). True for a script `train` in
  a language other than JavaScript, a `train` with `meta.mlserver: 'true'`, or a name in `SERVER_ENGINES`
  (`['Chemprop']`, whose functions do not say so yet).
- `hasTrainingRows(engine)`: `meta.mlhasrows: 'true'` on `train`.
- `handlesMissingValues(engine)`: `meta.mlnativenans: 'true'` on `train`; it makes `keep` the default missing-values
  mode (see "Missing values").
