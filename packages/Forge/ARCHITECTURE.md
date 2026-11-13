# Forge architecture

How Forge works now. Forge is predictive modeling as a package: engines that live in other packages, models stored in
an EMS schema owned by Forge, and the UI on top.

## Key decisions

- **Own schema instead of the platform model entity.** Models, training runs and applications are rows of the EMS
  schema `forge`. That gives typed storage, row security and an audit trail without core changes.
- **Engines stay in their packages.** Forge calls them through the `mlname`/`mlrole` contract and copies no engine code.
- **"Method" for users, "engine" in code.** Captions, messages and docs for users say method; code, the `engines/`
  folder, the `engine_*` columns and the `mlname` contract keep engine.
- **The platform's entry points open Forge.** On a server with Forge, the core's `ML | Models` menu, **Browse >
  Platform > Predictive models** and `grok.ml.applyModel` call Forge's functions by name and hide the earlier tool's
  panes (see "Core entry points"). Forge registers no menu of its own; its prediction columns carry the platform's
  `predictive.model` tag. Its storage stays its own (app, schema).
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
           │                    │                                                    ui/results-layout.ts
           │                    ├─> ui/runs-panel.ts ──> metrics/, ui/data-grid.ts, ui/results-viewers.ts (resultsTitle),
           │                    │                       engines/, preparation/, training/, generated/db.ts (types),
           │                    │                       ui/tune-settings.ts (TuneStatus)
           │                    ├─> ui/model-cards.ts ──> ui/metrics-table.ts (the exported pieces), ui/runs-panel.ts,
           │                    │                        ui/tune-settings.ts (STATUS_TOOLTIPS), ui/model-handler.ts
           │                    │                        (featureCount), metrics/, engines/, training/k-fold.ts (types)
           │                    ├─> ui/tune-settings.ts ──> training/search.ts (searchSpaceOf, rangeProblem), engines/, metrics/
           │                    ├─> ui/search-panel.ts ──> training/search.ts, metrics/, ui/data-grid.ts,
           │                    │                         ui/results-viewers.ts (resultsTitle)
           │                    ├─> ui/validation-inputs.ts ──> training/k-fold.ts, ui/collapsible-group.ts
           │                    ├─> ui/button-gate.ts (Train, Tune)
           │                    ├─> training/search.ts (tuneSearch) ──> training/k-fold.ts, training/train-model.ts, engines/
           │                    ├─> training/browser-search.ts (browserSearch: the trials) ──> training/search.ts,
           │                    │                                                            training/train-model.ts, metrics/
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
           ├─> ui/prediction-column-panel.ts ──> apply/prediction-tag.ts, ui/model-handler.ts (the card)
           └─> ui/table-models-panel.ts ──> catalog/ (modelsFor, model-activity), constants.ts (MAX_MODELS), generated/db.ts,
                                            storage/model-fields.ts (tagsOf), ui/model-handler.ts,
                                            ui/model-panes.ts (modelFacts), ui/apply-model-dialog.ts (openApplyDialog),
                                            ui/train-view.ts
```

- `package.ts` holds the annotated functions only: the app `forgeApp` (shown as "Forge", `Browse > Apps`),
  `forgeModels` (opens the app; no menu), `forgeTrain` and `forgeApply` (no menu of their own: the core's `ML | Models
  | Train Model...` and `Apply Model...` call them, see "Core entry points"), the API functions `applyModel`, `trainModel`, `deleteModel`, `listMethods`, `listModels`, `getModel`,
  `updateModel`, `evaluateModel` and `shareModel`, the autostart `_initForge` (registers
  the model handler with the model commands and the context panel's follower, and the comparison handler),
  `isPredictionColumn`, the column panel
  **Predicted by** (function `predictedByPanel`) and the table panel **Models** (function `tableModelsPanel`).
  They delegate to `ForgeApp`, `TrainView`, `openApplyDialog`, `runApplyModel`, `runTrainModel`, `runDeleteModel`,
  `runListMethods`, `runListModels`, `runGetModel` (with `grok.shell.tables`), `runUpdateModel`, `runEvaluateModel`,
  `openShareDialog`, `ForgeModelHandler.registerOnce`, `registerModelActions`, `ForgeModelHandler.followChanges`,
  `ModelComparisonHandler.registerOnce`, `isForgePrediction`, `predictedByPane` and `tableModelsPane`. `_initForge` reads the registered handlers' names once
  (`DG.ObjectHandler.list()`) and both `registerOnce` calls skip a handler already there, so a second `_initForge` (or
  the test bundle) registers nothing twice; the model commands and the follower are set up only with a newly
  registered handler.
- Every table of data in the UI but the model cards' (plain `ui.table` by the user's decision of
  2026-11-05, see "Model cards") is `readOnlyGrid` (`ui/data-grid.ts`): `DG.Viewer.grid` without editing, row header,
  current-row indicator and new-row icon, as wide as its host and as tall as its rows (up to 15 by default,
  then it scrolls: an inline height of the header plus `rows * rowHeight + 2` px, the header one row high until the
  grid draws, then `colHeaderHeight` after each `onAfterDrawContent` (a header whose names wrap, as the Runs metric
  columns', is taller), never less than a row; the app's Methods grid shows every method), with tooltips from `{headers, cells}` (by column name; the cells
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
  table (the boundary of the table pane's **Train**); `TrainView.open()` (`forgeTrain`, the core's `ML | Models | Train Model...`) calls it with the
  current table, or says "Open a table first.". The form is one `ui.form`, so the platform
  aligns every label in it: the group **Data** (Table, Target, Features), the group **Preparation**
  (`PreparationInputs`, below; hidden while none of its inputs applies), the group **Method** (Method, whose first item
  is `Auto`, and the method's `TuneSettings`, see "Search": **Settings**, its status, then the hyperparameters in
  Manual or **Search options** and the Tune bar with `tuneButton` in Auto-tune; under Auto the form's `FindSettings`,
  `form.find`, instead: the status `Not searched`, **Search options** with Metric, Strategy, Limit, Max time or
  Max trials and Candidates, and the bar
  with `findButton`) and the group **Validation** (`ValidationInputs`, below; collapsed
  at the start), the others expanded at the start, and last the row `trainRow` (`forge-train-row`, a
  `ui.buttonsInput` with `trainButton` at its right; `position: sticky; bottom: 0` on the `--white` background, so it
  sits right under Validation on a short form and stays at the bottom edge while a taller form scrolls). The view is a
  resizable split (`ui.splitH(..., true)`): the left box starts 400 px wide (`FORM_WIDTH`) and holds the form's panel
  (`forge-train-form`), which scrolls by itself; the results take the rest (`forge-train-results`, scrolling), with the
  training loader (`resultsStatus`) in the header of the pane's **Performance** panel.
  The ribbon holds only **Save** (`ui.bigButton`). The platform's column inputs keep a column removed from the
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
  a threshold re-cut): `update(summary)` rewrites the values in place, refills the bullet list and repaints the marks when
  the metric rows and the value columns are the same, so the grid neither flickers nor loses its column widths; other
  rows (task, AUC-ROC, the Crux ids) build a new one (the rows are the ids with a validation value or an evaluation
  value, in `METRIC_IDS` order: `balancedAccuracy`, `mcc`, `averagePrecision`, `logLoss` and `brier` after `auc`; a
  model saved before AUC-ROC, the Crux ids or the probability ids existed gets those rows with Train and Validation empty when evaluated; a change of the
  Test or Evaluation column builds a new one too); a hint shown in between detaches the grid and keeps it. It is a frame of **Metric**, **Train**, **Validation**,
  with `metrics.test` **Test** and, with `MetricsSummary.evaluation`, **Evaluation** (format `0.000`) in a
  `readOnlyGrid` with the header texts, the metric descriptions and the full-precision values as tooltips, then a
  bullet list (`ul`): `Missing values: kept (the method handles them)` (only under `keep` when the training features
  had gaps), `Rows: <used> used, <skipped> skipped (missing values)` (only when rows were skipped),
  `Validation: <folds>-fold cross-validation on <rows> rows` or `Validation: hold out of <p>% (<n> rows)` (from
  `MetricsSummary.splitting`, `splittingOf(model.splitting)` for a stored model; `, stratified` after it when
  `isStratified`; `heldOutCount` of the stored `trainFraction`), `Split seed: <seed>` (absent when `isShuffled` is
  false; an older splitting without it counts as shuffled) with the seed in a selectable
  `forge-seed` span and the copy icon right after it (**Copy the seed**, balloon `Split seed <seed> copied.`),
  then, with `MetricsSummary.search` (the training's `search`, or `storedSummary`'s `searchRecordOf(model.search)`),
  `Search: <Tune> by <metric>, <limit>: <N> trials of <M> method(s), chosen by the one-standard-error rule` (Find
  best model: `chosen by the highest <metric> − SE`, `the lowest <metric> + SE` for a minimized one)
  (`limitCaption`: `limit 5 min` (`<n> s` for seconds that are not whole minutes) or `limit 100 trials`, both joined
  by ` + ` for a continued search of both kinds, then `, Focus on leaders` for successive halving; `<N> trials of ` is
  left out without a time in the limit, as the limit gives the count: `limit 100 trials: 3 methods`) and either
  `Test: <rows> rows held out before the search` or `Validation after a search is optimistic: the best of the trials
  was chosen on it`, `Positive class: <class>` (two classes only), `Classification threshold: <t>`
  (`MetricsSummary.threshold`: the result's `binaryClassificationThreshold` in the view, `storedSummary`'s for a
  model whose options record `class-probabilities`; with `thresholdMetric`, the view's result's
  `options.thresholdMetric`, the optimized variant, see "The threshold in the view") and `Evaluation: <table> (<rows> rows[, <n> skipped
  (missing values)])` (with an evaluation). The `metricFindings` of the metrics mark their cells: `onCellPrepare` paints a marked value
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
  **Predict probability** (unchecked; shown while the target is categorical with two or more `classesOf` classes;
  `setProbabilityReason(reason)` from the view's `probabilityReason(form)` (`<Method> gives no class probabilities`
  when the method gives none; under Auto the Find block's reason, `No checked method gives class
  probabilities`, see "Search"): with a reason it is unchecked, `enabled = false` and the reason is its tooltip,
  otherwise `Adds the probability of each
  class when the model is applied` plus, for two classes, `; for two classes, lets you choose the classification
  threshold of <positive>.`; while the form is locked `fillMethod` keeps that `enabled` for the unlock and disables
  it again; its own callback `onProbabilityChanged`), under it **Classification threshold**
  (`ui.input.float` 0..1, step 0.01, with a slider, `DEFAULT_THRESHOLD`, `format: '3 significant digits'`: the
  platform's `FloatInput` keeps the full value and only shows it formatted, `0.520` for 0.5201..., `0.0423`, but `0.9`
  and `0.5` as they are, since the formatter drops the digits a value does not have (kept by the user's decision); a
  typed value is parsed in full; its tooltip `Rows whose probability of <b>positive</b> is at least this value are
  predicted as <b>positive</b>, the others as <b>negative</b>.` is the input's own (`setTooltip` with HTML, set by
  `update` from the target's classes): the platform binds its caption and editor to the input's tooltip, and those
  bindings win over one on the root, which therefore never showed while hovering them; a
  validator and `thresholdProblem` refuse an empty value or one outside 0..1; its own callback `onThresholdChanged`),
  and under that **Best threshold by** (`thresholdMetric`, a `ui.input.choice` of the `BEST_THRESHOLD_METRICS` labels,
  `DEFAULT_THRESHOLD_METRIC` F1; its `tooltipText` is `The metric the classification threshold maximizes on the
  validation rows.`, each `option`'s `title` the metric's description, Balanced accuracy's naming Youden's J; `metric`
  maps the label back to the id) with the button **Set** (`setThresholdButton`, a `ui.button` as Tune is, appended to the
  input's root after the select, not to its options, whose box would underline it and shrink its font; the root has
  `forge-threshold-metric`, and the Tune bar's rule drops the button's right margin and padding, so its letters end at
  the row's right edge where Tune's do). The order is `stepInputs`'. Both are shown
  while Predict probability is checked and the target has exactly two classes, and hide together. The show conditions
  count the classes `update` takes from a categorical target with `classesOf` (`metrics/metrics.ts`): Predict
  probability two or more, the threshold and its metric exactly two. The threshold's hand-edit state `isThresholdEdited`: a value change not made
  by `setThreshold` (which sets the value with the input's `notify` off) sets it and calls `onThresholdChanged`; `setThreshold(value,
  metric?)` sets the value (and the metric) without telling the view and clears it, or sets it without a metric.
  Hiding keeps the metric, the threshold and the hand edit (the user's decision 3b): checked again, the view's re-cut
  optimizes an unedited threshold for the kept metric, and a hand edit stays. `update` resets them
  (`setThreshold(DEFAULT_THRESHOLD, DEFAULT_THRESHOLD_METRIC)`) when the target's `classesOf` classes change (another
  target column, or none). `steps()` carries `thresholdMetric` (the chosen metric) unless the threshold was edited
  by hand, and `predictProbability` as the checkbox's value (a locked form, disabled, keeps it; a hidden input or a
  method without probabilities has it unchecked). The view's `bestThreshold(metric)`
  callback gives the best value, undefined without one; `updateSet` (called on a metric change, by `setThreshold`,
  on a visibility change and from the view's `updateControls`) shows **Set** only while a best exists and the input
  differs from it (hidden, not disabled, by the user's decision; its tooltip, bound with `ui.tooltip.bind`, is `Set the
  threshold that maximizes <metric> on the validation rows.`), and the status `thresholdStatus` under the threshold
  (`TuneSettings.statusLine()`, the italic line Settings has under its radios: `forge-tune-status-text` in a
  `ui-input-root forge-tune-status` row) only while a best exists: `Best for <chosen metric>` while the value is its
  best, else `Customized` (the Settings status's word) after a hand edit, else `Best for <optimizedFor>`, the metric
  `setThreshold` last set the best of. A click
  (ignored while the form is locked) calls `setThreshold(best, metric)` and `onThresholdChanged`, which re-cuts.
  `restoreThreshold(steps)` sets a shown threshold from a run's steps (with their metric, else as a hand edit).
  `update(features, target)` runs on every Features and Target change; a hidden input but the threshold pair gets its
  default back, so it reappears as it first appeared, and the group is hidden while nothing in it is shown. `steps()` is the
  `PreparationSteps` of the selection; a hidden input's step is off (Skip unique categories, checked while hidden, is
  sent only while it is shown), and Predict probability only while enabled.
- **Validation** (`ui/validation-inputs.ts`, `ValidationInputs`, Train view only) owns its `CollapsibleGroup`
  (collapsed): **Scheme** (`K-fold, stratified` | `K-fold` | `Hold-out`; `update(target)` on every Target change lists
  the stratified one only for a classification target, turns it into `K-fold` for a regression and gives it back on
  the next classification target unless the user picked K-fold, `isSchemeChosen`), **Folds** (`ui.input.int` 2..20,
  shown for the K-fold schemes), **Test fraction** (`ui.input.float` 0.05..0.5, shown for Hold-out), **Shuffle** and
  **Split seed** (`ui.input.int` 0..`MAX_SEED`, 42, `showPlusMinus: false`, enabled while Shuffle is checked, no icons
  after it: the user removed the dice and the copy icon). Every input has a tooltip of at most five words (`How rows
  are split`, `Seed for random splitting`, ...).
  `settings()` is the `ValidationSettings` (`training/k-fold.ts`; Hold-out is stratified, which a regression ignores),
  `seedValue` the split seed, `usedInputs` the inputs `problem()` checks (Folds or Test fraction, the seed while
  shuffled), `restore(validation, seed)` sets them as a user would (a **Runs** row). A Table change builds a new form
  whose `ValidationInputs` takes the old one's values (`createForm(table, previous)`).
- **Settings the method ignores** (`meta.mlenabled`, see "Engine contract"). `applyRules(form)` runs at the end of
  `fillMethod`, on every hyperparameter change and on every Target change (`formSubs`): it
  evaluates `isEnabled(enabledRules(engine), name, values, task)` per hyperparameter, with the task the method trains
  (`taskOf(target)`, unknown without a target), sets `input.enabled` (never
  hides) and passes the ignored names to `form.tune.setInactive` (read back by `isActive`; in Auto-tune with the tuned
  settings free and the target's task, see "Search"). `problem()` and `revalidate()` skip the inactive inputs
  (`activeHyperparameters()`), so a disabled invalid value blocks nothing; `selection()` still sends every value. A
  malformed `mlenabled` logs one warning per method (`loggedFailures`) and enables everything.
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
themselves (a list, no frame), the target column, the hyperparameters, the split seed (`seed`), the validation
(`validation`, a `ValidationSettings` of `training/k-fold.ts`: `scheme` `kfold` | `holdout`, `folds`, `testFraction`,
`stratified`, `shuffle`; the view's **Validation** group, `DEFAULT_VALIDATION` for the API), the
missing-value handling and the preparation steps (`steps`, a `PreparationSteps`, see "Preparation"; the Train view
takes them from `PreparationInputs.steps()`). `prepareTraining(selection)`
handles the missing values (see "Missing values"), then applies the steps (`prepareFeatures`), and returns the
`TrainingRequest`: `features` and `target`, the selection's columns with the missing values handled (the user's own
columns where nothing changed; the model's feature list, target schema and data summary describe them), `keptRows`
(the selection's rows they hold, a `BitSet` over the table, or `null` when no row was skipped), and
`prepared` (`PreparedFeatures`), what the method is trained on, with the `options` that record the whole preparation
(`options.missingValues.skippedRows` is the one record of the skipped rows). No frame is built: the engine calls
build their own. After preparation fewer than `minimumRows(validation)` rows (`2 * folds`, or two held out rows:
`ceil(2 / testFraction)`) is a `ForgeError` ("After skipping rows with missing values, N rows remain; training needs
at least 10.").

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
| 3 | fewer rows with a target value than `minimumRows(validation)` (`2 * folds`, 10 by default; `ceil(2 / testFraction)` for a hold out) | Target |
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
columns and target (the selection's target: Predict probability lists the same classifiers), and rule 8
becomes: the selection's method is not among them (the rule 8 text), or none is (`No method can learn from this
selection. Check the features and the target.`). The result (`SelectionCheck`) also has `best` (`selectBestEngine`
over the listed methods, the prepared columns (one-hot copies count as numerical features) and the selection's target;
the molecule rule reads the selection's features, as
preparation drops or one-hot encodes a molecule column), `failed` (the methods whose check threw, for the caller to log),
`isInteractive`: the selection's method is listed and `retrainsLive` (its `isLiveUpdate` is on and its
`isInteractive` says yes on the prepared columns), `prepared`, those columns and target (none from
`failedCheck`), and `insights`, the data insights of the prepared columns and the selection's target (see "Insights";
empty from `failedCheck`). A method's `isInteractive` that throws makes the whole check throw.

The view adds a validator to **Target**, **Features**, **Method** (the `method` messages) and **Missing values** that
returns the view's current problems for that input. Every change of Table, Target, Features, a **Preparation** input
but Predict probability (of a method without a probability input, see "The threshold in the view") and the
threshold, Method or a hyperparameter (`requestCheck`) numbers a new check, clears the problems, drops
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
gone. The disabled reason is also bound as the button's own tooltip (the last binding wins), so on enabling the gate
writes the enabled tooltip into that element instead of binding another listener (Train and Tune), and the last
disabled reason does not stay on it; a gate without an enabled tooltip (the dialogs' OK) binds none on enabling, since
an emptied element would still show as an empty box. `update(reason)` takes a reason computed by the caller (the view computes `problem()` once
for both gates); `release()` enables a hidden button and forgets the reason, so no overlay timer polls for it. A check that fails shows **Train** disabled with
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
method); when the only problems are invalid settings (hyperparameters, missing-values settings or the
**Classification threshold**), it reads `Fix the settings.`. `isTraining` is true from the start of a training until its run is recorded; **Save**
waits for it. A result
that arrives after a newer change is recorded but not shown.

**Task.** A numerical target (not a date) is regression; text and boolean targets are classification. The task and
the target schema follow the selection's target, which is also the prepared one (Predict probability changes no
target since 0.1.20).

**Validation.** `validationSplit(request)` gives every row's fold by the request's `validation` and split `seed`
(`splitRows` in `training/k-fold.ts`), stratified by the selection's target when `validation.stratified` and that
target is a classification (a regression never is).
`kFold(rowCount, folds, seed, stratify?, shuffle = true)` shuffles the rows with a mulberry32 generator
(Fisher-Yates) and deals them round-robin into the folds, so fold sizes differ by at most one and every row is
validated exactly once. With `stratify` the shuffled rows
are grouped by class and each class is dealt in turn with one running counter, so every fold holds each class's share
(a class's counts differ by at most one across the folds). Without `shuffle` the rows keep their order and the seed is
unused: contiguous folds (`floor(j * folds / n)`), or each class dealt round-robin in row order. `holdOutFold(rowCount,
testFraction, seed, stratify?, shuffle)` marks `heldOutCount` rows (`round(testFraction * n)`, at least one, one row
left) 0 and the others -1: the last ones of the (shuffled) order, with `stratify` each class's share by `groupQuotas`
(the largest remainders round up; shared with the search's own hold out) from its own rows. For each predicted fold
(`validationFolds`: k, or the one hold out) the engine trains on the other rows and
predicts the fold's rows (`foldPrediction`); the predictions form one out-of-fold column (`outOfFoldColumn`; a row of
fold -1 is null there, and `FLOAT_NULL` in `outOfFoldProbabilities`), and the **Validation** metrics are computed once,
on the held out rows' clones for a hold out (so the EDA report never sees a null). `outOfFold` keeps that column and
the out-of-fold probabilities with their nulls; `recutResult` measures them with `metricsOf`, which skips rows
without a prediction, and the threshold optimization reads only the rows with a probability. The **Train** metrics come from the final model, trained on all rows
and applied to them. A class with fewer rows than folds still misses some folds. The split seed is the view's **Split
seed** (42 by default, the API's `options.seed` too); the seed and `splitting` (`trainingSetupOf`: `{scheme: 'kfold',
folds, isStratified, isShuffled}` or `{scheme: 'holdout', trainFraction: 1 - testFraction, isStratified, isShuffled}`;
the schema's column description since 0.1.18) are saved; `splittingOf` reads a stored one. A request with `holdOut`
(a search's, see "Search") also applies the final model to the held out rows: `metrics.test` (`measure`, without
per-class values or curves).

**Progress and cancel.** `trainModel(request, progress?)` (a `LoopProgress`, the cancel flag and `update` of a
progress indicator, in `engines/engine-calls.ts`) reports "Fold i of k" (a hold out: "Hold out") after each fit and throws
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
record's row count (0 when nothing was recorded), the validation settings and the split seed, the settings status
(`settingsStatus`: `form.tune.status` read when `startTraining` starts, so a search's winner is `Auto-tuned`), and the
metrics or the error. A training of a table the view has left adds nothing
(`createForm` replaces the list, the training compares it). `RunsPanel` shows them under the results
(`ui.panel([resultsPane.root, searchPanel.root, runsPanel.root], 'forge-train-results')`, the **Search** panel of
"Search" between them; the panel's own root is `forge-runs`, 12 px under the rows): the core's title `Runs` with its tooltip
(`PANEL_TEXTS.runs`, `resultsTitle`), hidden while the list is empty (`detach()` empties it, which drops the grid's
click subscription), and a `readOnlyGrid` with the current-row
indicator on, built again on every `update` over a copy of the list (a new frame, a new grid of the new height; the
replaced grid goes with `ui.empty`, the platform's widget disposal): **#**, **Method**, **Target**,
**Features** (the count), **Preparation** (`impute (<n> neighbors, <distance>)`, `one-hot`, `skip unique
categories`, `probability (threshold <t>, best <metric>)` for two classes (`SessionRun.twoClasses`, from
`twoClasses(target)` in `runInputsOf`) of a completed run whose steps carry `thresholdMetric`, `probability (threshold
<t>)` without it, `probability` for more, joined by commas; the model cards show the same items), **Rows**, **Status**, **Split seed** (null when the run
was not shuffled), then `<metric>
(validation)` for every metric of `METRIC_IDS` any run's validation has, with header tooltips (the metrics'
`METRIC_DESCRIPTIONS`) and cell tooltips (the hyperparameters without the method's `probabilityInput`, the feature
names, the error or `Stopped from the task
bar.`, `Train <0.000>, validation <0.000>` with `-` for a missing value). The newest row is current; a click on a
row (`grid.onCellClick` on a table cell, the current row too) calls `restoreRun`, while code that moves the current
row restores nothing (the read-only grid moves nothing with the arrow keys): a run whose target or a feature the
table no longer has warns `The table no longer has the columns of this run: <names>.` and sets nothing; a run whose
inputs equal the form's (`isSameInputs` of `ui/runs-panel.ts` over `runInputsOf` of the current selection: method, target, features,
hyperparameters, missing values, preparation steps but Predict probability, the threshold and its metric, the
validation and the split seed, the same values a run keeps; SVM's runs differ by their `probability` hyperparameter)
while the view shows the completed result of those inputs (`lastTraining`, which any change but those drops) trains
nothing: it sets the run's Predict probability when it differs and its threshold with `restoreThreshold` (with the
run's metric, else as a hand edit), then re-cuts the shown result; after a
cancelled or failed
training there is no such result, so the click sets the same inputs and trains; otherwise it stores the run's hyperparameters for its
method, makes that method the user's choice and the current method (its inputs rebuilt with those values), sets
**Target** and **Features**, then `PreparationInputs.restore`
(Missing values with the impute settings, then Skip unique categories, One-hot encoding, Predict probability and the
threshold with its metric, each only while shown; a hidden input keeps its default), then `ValidationInputs.restore` (the run's
validation and split seed) and requests a check: the view retrains as after any change, and the new training is a new
first row. While the form is locked (`locked`: a **Train** training or a search, see "Search") `restoreRun` and `restoreLeader`
set nothing.

**Model cards.** `ModelCards` (`ui/model-cards.ts`, `TrainView.modelCards`) is a column of cards at the end of the
`ui.form`, after the **Train** row (`formRoot.append(...groups, trainRow, modelCards.root)`): the cards take the
form's width (`ui.css` caps a `.ui-form` at 700 px), so their right edge is the inputs', **Tune**'s and the Settings
reset's in any pane width. The sticky **Train** row is still held at the pane's bottom edge only while its place is
below it, then scrolls away above the cards. The colors are the core's (`predictive_modeling.css`) in tokens: a pinned
card `color-mix` of `--grey-1` and `--white` (#f9f9f9), the current one `--blue-1` at 8% over transparent
(#b8d3ed3b), the border `--grey-2` (#ccc). A card is a plain `div` (`ui.element`), not a `ui.div`: a tab control in
a `ui-div` gets the platform's 400 x 300 px box (`ui.css div.ui-div > div.ui-box`), so a card is as tall as its shown
pane.
`addRun` of a completed run calls `setCurrent(run)`: the current card (`forge-card-current`) is
built from the same `SessionRun` as the Runs row, and the previous current card is disposed; a failed or cancelled run
leaves it. A card is a header (the method in a bold `forge-card-title`, then its icon: the current card's "+", a
pinned card's trash), for a run with a `search` its caption (`CARD_CAPTIONS` by `search.mode`: `Auto-tuned model`
after Find best model, `Auto-tuned settings` after Tune, in `forge-card-caption`, the italic look of
`forge-tune-status-text`; nothing for a run without a search; the user's decision of 2026-11-10 dropped the line of
the features and the preparation, which Runs shows) and a `ui.tabControl` with two panes:
**Performance** (`performance`, built once with the card), a `ui.table` of the
run's metrics: the first row the run's `search.metric`, else `defaultMetric` of the task (Balanced accuracy, R2), then
the rest of `metricIdsOf` in its order, the rows after the first `SHOWN_METRICS` (3) hidden until **Show more**
(`forge-card-more`, a `ui.link` with the tooltip `Show all metrics`; `Show less`, `Show fewer metrics`, while expanded;
none with 3 metrics or fewer) shows them; the state is the card's own (the user's decision: unlike the tab, not
shared), collapsed on a new card, and a click only shows or hides those rows and relabels the link; a pinned card's
click on it restores nothing. The value columns (`valueColumns`),
the findings (`findingsOf` over `metricFindings`) and the tooltips (`METRIC_HEADERS`, `metricTooltip`;
`METRIC_DESCRIPTIONS` of `metrics/`) that `ui/metrics-table.ts` shares with `MetricsTable`, values `toFixed(3)`, uncolored (the
user's decision: no failure color on a card; a value's findings stay in its tooltip); **Settings** (lazy), a
`ui.table` (**Setting**, **Value**; `No settings` without one) of every setting of `visibleHyperparametersOf(run.engine)` (captioned as
`ui.input.forProperty` captions it: `camelCaseToSentence` of the caption or name; the description as the row
tooltip). Both tables are `forge-card-table`: one row height and padding, the names in `--text-color-light`, the
values in `--grey-5` (as the platform's key-value tables) and right-aligned with their headers, the header row in
`--text-color` at weight 500 (darker than both, lighter than the bold title). The plain tables are the one exception to "tables of data are grids" (the
user's decision: the core's look, no Dart grid per card). The current card's "+" (`forge-card-pin`) calls `pin(run)`:
a pinned card of the same run goes right under the current card, newest first, with a trash (`forge-card-unpin`)
calling `unpin(run)` (the pinned card of that `SessionRun`); a pin whose inputs a pinned card has (`isSameInputs`)
adds nothing. `refreshPin` hides the
current card's "+" (`visibility: hidden`, so the header keeps its height) while a pinned card has its inputs, as the
core did; the current card has no trash, and `unpin` of that pinned card shows "+" again. A click on a pinned
card outside its tab headers (`TabControl.header`) and its icon calls `restoreRun` (its guards apply). Each card has a
key (`ModelCards.keyOf`: the target, the scheme, the folds or the test fraction, stratified, shuffle and the seed when
shuffled); `updateControls` calls `update(key)` with the form's key (undefined while there is no form or target) and
`show` displays only the cards of that key, the current one too, at once on any change. Tabs switch together:
`onTabChanged` of a card calls `syncTabs`, which keeps the shared tab and sets `currentPane` of every card on another
tab (the echo of that set finds every card on the shared tab and sets nothing); a new card opens on the shared tab.
`dispose` unsubscribes the card's `onTabChanged` and `ui.remove`s the card (the TabControl's widget disposal);
`clear()` runs in `createForm`
before `ui.empty(formPane)` and in `detach`. Nothing is written: the cards and the pins live with the view.

The result (`TrainingResult`) carries the blob, the task, the metrics, the target and feature schemas, the
preparation options, the splitting, the seed, the hyperparameters, the row count (the prepared rows the model
learned from; the skipped rows are in `options.missingValues`) and `prediction`, the final model's prediction of
those rows as the method returned it; `trainPredictionOf(result)` gives its labels (`replayPostprocessing(prediction,
options, probabilities)`: with Predict probability those of the class probabilities; for a result `recutResult` built,
the train labels it built, kept in a `WeakMap` by result). The run and model rows record
that row count, and the dataset fingerprint is computed on the request's `features` and `target` (before the
preparation steps). It also keeps `actual` (the selection's target the predictions are measured against), `outOfFold`
(`{prediction, probabilities?}`: the out-of-fold column, null for a row of no fold, and the out-of-fold
probabilities) and, after a search with a hold out, `test` (`{actual, prediction, probabilities?}` of the held out
rows). When the method returns class probabilities (`proba_<class>` columns after its prediction, see "Engine
contract") the result keeps the final model's ones of the training rows as `probabilities` (`ClassProbabilities`
`{classes, columns}`); the out-of-fold ones are gathered by class name (`outOfFoldProbabilities`, named as the final
model's columns; `positionsOf(fold)`, each row's index in its fold, is computed first and passed to both it and
`outOfFoldColumn`). `trainModel` returns the
plain classifier's result, put through `recutResult` only when the prepared options record `class-probabilities`.
`recutResult(result, {predictProbability, threshold, thresholdMetric?})` (synchronous, no Crux call) applies or
removes the step without retraining: the options lose the id, the three two-class keys and `thresholdMetric`, then
with `predictProbability` get the id back and, for two classes (`twoClasses(actual)`), `positiveClass`,
`negativeClass` and `binaryClassificationThreshold`: with `thresholdMetric`, `bestThresholdOf(result, metric)`
(`bestThreshold` of the out-of-fold probabilities of `metrics.positiveClass`, see "Metrics"; computed once per
training and metric, cached in a `WeakMap` by the result's `outOfFold`, which every re-cut keeps, so the re-cuts and
**Set** read the same value) and `thresholdMetric`
itself, or `threshold` without a metric or when there is no best (no out-of-fold probabilities, one class on the
validation rows); `trainModel`'s own re-cut passes the prepared options' `thresholdMetric` (recorded by
`prepareFeatures` from the steps), so a training through the view, the API or a search's winner returns the optimized
threshold. The labels of train, `outOfFold` and `test` are
`replayPostprocessing(prediction, options, probabilities)`, and the result gets the TS metrics of those labels
(`metricsOf` per split; per class `classLabelMetrics`: precision, sensitivity, F1), the probability values
(AUC-ROC, average precision, log loss, Brier: `SCORE_METRIC_IDS`, overall and per class) kept from the trained
result; the result's `curves` stay. A result without `probabilities` gets no step: it is returned as it is, or
without the step when its options had it (a method that declares `mlprobabilities` but returned none for this data,
such as EDA on classes that differ only in case), so Apply then adds nothing. Each split is one `measure(task, actual, predicted, positiveClass?,
probabilities?, curves)` call (see "Metrics"; the two run together): the train split with `curves`, whose `curvesOf`
becomes `TrainingResult.curves` (two classes with a positive class), both splits' reports giving `metrics.classes`
(`classMetricsOf`; both in `recordOf`).

**The threshold in the view.** A **Classification threshold** change requests no check: it fires the same
`CHECK_DELAY_MS` timer (`isCheckRequested` tells the two apart; a check wins), and `recut()` replaces
`lastTraining.result` with `recutResult(result, preparation.steps())` and shows it, without a training or a run row;
**Save** then writes the re-cut `options` and metrics. A **Predict probability** toggle (`onProbabilityToggled`)
does the same for a method without a probability input; for one with it (`probabilityInput(engine)`, SVM) it
requests a check, so the method retrains with the input set (`selection()` sets it through `withProbabilityInput`
from `steps().predictProbability`, so runs, searches, cards and the API carry it in `hyperparameters`; the form never
shows the input: `fillMethod` and `TuneSettings` read `visibleHyperparametersOf`). A completed training is re-cut at
the current steps before it is shown (they may have changed while it trained); its run row keeps the steps it trained
with, the threshold and its metric as the training resolved them (below). An invalid threshold (`thresholdProblem`) computes nothing: `recut()` keeps the old result out of sight and
shows `Fix the settings.`, **Save** is disabled, and `problem()` blocks **Train** and the live training (`Classification
threshold: <problem>`), as for any invalid setting. When the timer fires with no result to re-cut, no training running
and no finished training of the selection (another input changed while the threshold was invalid), it runs the check
instead, so a valid threshold then trains as any fixed setting does.

**The optimized threshold.** While the threshold is not edited by hand, `preparation.steps()` carries the
**Best threshold by** metric as `thresholdMetric`, so every `recut()` (Predict probability checked, a completed
training, **Set**, a restored run) and every training (its prepared options) resolve the threshold to that metric's best on
the result; `recut()` then writes the resolved value back into the input with `setThreshold(value, metric)` (no
`onThresholdChanged`, no hand edit). A hand edit drops the metric from the steps, so re-cuts and trainings keep the
typed value until **Set**; the state lives in `PreparationInputs` and survives `requestCheck` and the hiding of the
inputs (unchecking Predict probability), and goes only when the target's classes change. The input shows the value
with 3 significant digits (its `format`), while `value`, the steps, the options, the texts and **Set**'s comparison
use the full value. A metric change only
updates **Set** and the status: the threshold moves at the next re-cut or training. The view's
`bestThreshold(metric)` gives the best value of `lastTraining` (`bestThresholdOf`), undefined without a training,
validation probabilities or two classes in the validation rows; `updateControls` calls `updateSet()`. `showResults` passes `threshold` and
`thresholdMetric` from the result's options, so **Performance** reads `Classification threshold: <t>, the best <metric>
on the validation rows: its validation values are optimistic` while the result's threshold is the optimized one, and
`Classification threshold: <t>` after a hand edit; `storedSummary` passes the saved options' `thresholdMetric`, so a
saved model's **Performance** pane reads the same line through the same `MetricsTable.fillList`. A completed
training's run row and its `training_run` record get
the result's options and the resolved `threshold` / `thresholdMetric` in their steps (a failed or cancelled run keeps
the steps it was started with; its Preparation text shows no metric).

## Search

**Tune** searches the settings of the current method, **Find best model** (mode `best`: **Method** Auto, or
`Forge:trainModel` with `method: 'auto'`) several methods; the search is Forge's own loop and every trial is a regular
fold training (`foldPrediction` through the method's `train` / `apply`), so a trial sees the same preparation, split
and metrics as a training. The logic is `training/search.ts` (types, spaces, sampling, hold out, candidates, plan,
leaders, record) and `training/browser-search.ts` (the executor); the UI is `ui/tune-settings.ts` (Settings and Search
options), `ui/find-settings.ts` (the Find best model block), `ui/search-panel.ts` and `TrainView`.

**Spaces.** `searchSpaceOf(engine)` reads `meta.mlspace` of the method's `train` function (`searchSpaceText`; see
"Engine contract"): the `ml.tune` `space` JSON (`{"<input>": {"type": "int" | "float" | "categorical", "low", "high",
"log"?, "step"?, "choices"?}}`, shorthands `[low, high]` and a list of choices) plus an optional `simpler: "low" |
"high" | "first"`. It checks every name against the train inputs (`The search space of <Method> names no setting
<name>.`), numeric bounds by `rangeProblem` (below; `The search space of <Method> gives <name> an invalid range:
<problem>.`), the choices against the input's, and `simpler` (`first` only for a categorical); a malformed space throws.
`tuneCandidate(engine, start, tune?, space?)` is the one candidate of Tune (`<Method> declares no search space, so it
cannot be tuned.` without a space; `space` is the parsed one the view already has, read from the engine when absent). With the user's `tune` (`{space: TuneSpace, task}`, from Search options;
`SearchChoice.space`, which the API never sets) its space is `userSpace(space, tune, rules, start, task)`: the pinned
dimensions (`TuneSpace.pinned`) left out, and those whose `meta.mlenabled` rule fails with the fixed values and the
task whatever the tuned ones are (`isEnabled(..., free)`: a clause on a free, tuned setting may hold); the others with
the user's bounds or choice subset (`TuneSpace.ranges`: `TunedRange {low?, high?, choices?}`), their type, `log`,
`step` and `simpler` kept, so `simplicityOf` measures within the user's bounds. The candidate keeps `pinned`.
`rangeProblem(dim, property, low, high)` is the UI's check of a range (two numbers, whole for an int, `low < high`
(`FROM_BELOW_TO`), inside the input's `min` / `max`, `low > 0` on a log scale, `>= 1` for a log int).

**Sampling.** `drawnPoints(space, start, seed, name)` yields the candidate's draws one after another from its own
stream `mulberry32(seed ^ fnv1a(FNV_OFFSET, utf8(name)))` (the FNV-1a of `storage/dataset-fingerprint.ts`), every
dimension drawn in space order, one `u` each (`sampled`: int `low + floor(u (high - low + 1))`, log int rounded from
log-uniform, float uniform or log-uniform rounded to `FLOAT_DIGITS` (3) significant digits (`toPrecision`) and clamped
to `[low, high]`, so the model trains on the value the form shows, categorical uniform; `step` rounds ints; the start
is never rounded), the other settings kept from `start` (the form's values, or
the method's defaults with the API's `hyperparameters`), so the i-th point of a method does not depend on the
interleaving. `trialPoints(candidate, seed, task)`, built once per search, yields the trials' points: each passes
through `activePoint(space, point, rules, defaults, task)`: each dimension whose `meta.mlenabled` rule
(`enabledRules`, `isEnabled`) fails for the drawn point, with the task the trials train (the search request's prepared
target), takes the method's default (`defaultHyperparameters`; not the form's value, so the leaderboard shows the
value the method documents; a dimension without a default keeps its draw), so points that differ only in ignored
settings collapse. The rules, the defaults and the start's place are read once. The first point is `start` only when
`isInsideSpace(space, start, rules, task)`: each value the method uses (by its rule) lies inside its dimension (a
choice of a categorical, within the bounds and whole for an int); otherwise the draws start at once, so every trial
lies inside the space (XGBoost's default `alpha` 0 is outside its range). Without a space the start is the only point.
A point already tried is skipped; `MAX_DUPLICATES` (20) in a row exhaust the candidate (the rounding makes every
space finite: a float range narrower than its rounding, `[0.1, 0.1001]`, draws one value).

**Hold out.** `holdOutSplit(request, seed)`: with at least `HOLD_OUT_MIN_ROWS` (500) prepared rows,
exactly `testRows = clamp(round(0.2 n), 100, 2000)` rows become `HoldOut {columns, target, rows, fraction}`: the first
ones of the rows shuffled by `seed + 1` (`shuffledRows`); for a classification each class (`classGroups`) gives its
first `share * testRows` rows, the quotas floored and the largest remainders rounded up so they sum to `testRows`
(each class within one row of its share) (clones of the prepared columns and target); the
search request holds masked clones of the rest with the whole table's preparation (`oneHotCategories`, options) and
`keptRows` mapped to the table's rows (`maskedKeptRows`), so the results' viewers still link to the table. Below 500
rows the request is returned as it is.

**Plan and executor.** `searchPlanOf(mode, request, candidates, choice, seed)` checks the metric against the task
(`SEARCH_METRICS`: regression `r2`, `rmse`, `mae`; classification `balancedAccuracy`, `accuracy`, `f1`, `mcc`, `auc`, `logLoss`, `brier`;
the first is the default; `The search metric '<id>' does not apply to a <task>. Choose one of: ...`), sets the
direction (`METRIC_DIRECTIONS`: -1 for errors), the ONE limit (`budgetMs` from `budgetSeconds`, or `maxTrials`;
both throw `COMBINED_LIMITS`, `The search options 'budgetSeconds' and 'maxTrials' cannot be combined: choose one
limit.`; neither gives `DEFAULT_BUDGET_SECONDS`, 60), for `best` the `strategy` (`equal` by default, `halving`; an
equal split by fewer trials than candidates throws `At least one trial per method: <N>.`) and the split (`validationSplit` of the
request with the plan's seed: the **Validation** group's folds, stratification and shuffle, or its one hold out; the
same split the winner's training gets), and `parallelism = PARALLEL_FITS` (`max(1, min(4,
floor(hardwareConcurrency / 2)))`, EDA's pool size) when every candidate's `train` declares `meta.mlparallel: 'true'`
(`isParallel`), else 1. `runSearch(plan, progress, executor)` runs `SearchExecutor.run(plan, progress)` and then the
leaders and the choice: the executor is the seam where a Crux-side search (or another scheduler) can replace the
browser one without touching the plan, the record or the UI. `BrowserSearchExecutor` (`browserSearch`) launches
trials while slots are free (a `Semaphore(parallelism)` gates fold fits; each fit gets its own masked clones, so
concurrent fits never share a column): a trial runs folds 0-1 together, then the checks, then the other folds
together; each wave waits for all its fits (`Promise.allSettled`) before its first error is rethrown, so `run()` never
returns with a fit in flight; once a fold of a trial fails, its folds still waiting for a slot skip their fits. The
fold masks (`foldMasks`) and the held out targets are built once per search. Abandon (a time limit only): once a
trial has completed (an equal split: once a trial of its own candidate has), a trial whose first two folds forecast
`elapsed + forecast` past the deadline (`budgetMs`; an equal split: its candidate's share end) stops (`abandoned`), where `forecast = fitMs / 2 * rest / min(parallelism, rest)` (`fitMs`: the
sum of the first two fits' times, `rest`: the other folds, fitted up to the parallelism at once; other trials' fits
sharing the slots are not counted). Pruning (as `ml.tune`'s median pruner), per candidate: trial `n` of a candidate
(`Trial.number`, the candidate's own count; 0 is its first trial, never pruned) is judged only against the
candidate's trials numbered below `n` that passed their first folds: each trial publishes `early[candidate][n]`, a promise of its first
`PRUNE_AFTER_FOLDS` (2) fold scores once it goes on to its other folds, or undefined when it was pruned, abandoned or
failed there (resolved in `finally`, so it always settles), and keeps the passed scores in `passedScores`; from
trial 5 on a trial waits only for the `early` of the earlier trials still in flight (`inFlight`; not for them to
complete), then reads the passed ones' scores from `passedScores`, and with at least `PRUNE_MIN_TRIALS` (5) of them passed, a trial whose
mean of its first two fold scores (signed by the direction) is below their median is `pruned`. So the prune decisions
do not depend on the finish order, and a parallel search gives the same trials, states and leaderboard as a sequential
one; the 30-trial iris test runs in about 1160 ms one
fit at a time and 480 ms with four. With no folds left after the first ones (a hold out, or 2 folds) nothing is
pruned. A complete trial's `score` is the metric of the pooled out-of-fold predictions (as
`trainModel`'s Validation, with the out-of-fold probabilities for `auc`, `logLoss`, `brier`, passed to `metricsOf`
only for a metric of `SCORE_METRIC_IDS`; with a hold out, the held
out rows only, the TS metrics skipping the nulls); `foldScores` are the per-fold values (one with a hold out). Both
score the labels `replayPostprocessing(prediction, request.prepared.options, probabilities)` gives, so with Predict
probability a trial is scored at the form's threshold (or the most probable class); the winner's training then
optimizes the threshold again through its prepared options' `thresholdMetric` (absent after a hand edit or an API
`threshold`, see "Training"). A
metric the method cannot give (`auc`, `logLoss` or `brier` without probabilities) fails the trial (`The metric <label> needs class
probabilities; <Method> gives none.`), and so does any error of a fit (`failed`, with `error`). Stops: `cancelled`
(`progress.canceled`, read before every fit: fits in flight finish, nothing new starts, a cancel between fits marks
the trial `abandoned`), `budget` (checked before every launch), `trials` (`maxTrials` trials launched, whatever their
state: complete, pruned, failed or abandoned, so the candidates' trials never add up to more and parallel and
sequential runs launch the same trials), `failed` (`MAX_FAILED_WITHOUT_COMPLETE` (10)
failed trials while none has completed; launches are capped meanwhile, so such a search has exactly 10 trials; not
for an equal split, where the same count of one candidate's own failures ends its share, and a split where every
share ended without a complete trial stops `failed`) and `exhausted`. Trials are kept in launch order;
`progress.update(percent of the time, or of the trials settled in this run, text)` follows every launch and every
settled trial (counts, not a pass over the trials), at most once per `REPORT_MS` (250 ms; one within it is sent when
it ends), and once more when the last trial settled. The text is short so the task bar never cuts it: `Finding best
model` (`mode: 'best'`) or `Tuning <Method>` (the candidate's `engine.name`), then in parentheses, with a trial limit,
`settled/maxTrials` (this run's finished trials, as the bar's share, so it moves until the last running trial
settles), with a time limit `N s left` (whole seconds, under a minute) or `N min left` (rounded up), and once the time
is up while trials still run `finishing N trials` (`finishing 1 trial`; launched minus settled, falling as they
settle); no score and no mark of a continued search. The first update replaces the indicator's title (`searchTitle`: `Tuning <Method>` or
`Finding best model`, the executor's action too; `Continuing: …` set by the view).

**Continuing.** `run(plan, progress, previous?)` continues `previous` (an outcome of a plan `canContinue` accepts):
its trials come first in `trials`, `elapsedMs` adds up, and the limit and the failed-trials count
apply to the new trials only, while `completed` (for the abandon rule) counts all. The
outcome's `resume` (`ExecutorState`) is what the executor needs for that, per candidate: `drawn` (the points taken
from its stream, duplicates too, so a continued stream goes on where it stopped), `exhausted`, and `passed` (each
trial's first fold scores when it passed them, else null: the pruning's `passedScores`); with
successive halving also `halving` (an equal split keeps nothing more: a continued one splits its new limit equally
again over the candidates not exhausted). The tried keys and each candidate's trial count come from the trials. So a
halving or Tune search
continued after N trials with a cap of M launches the same trials, states and points as one run with a cap of N + M
(a trial abandoned at the end leaves its point tried and its halving slot free).

**Leaders and the choice.** `leadersOf` applies the one-standard-error rule within each candidate: among its complete
trials whose score is within one SE of its best trial's (`se` = the sample standard deviation of `foldScores` over
`sqrt(folds)`), the smallest `simplicityOf` its point (the distance from each `simpler` end, 0..1 per dimension, in log
space for a log dimension, index / (n - 1) for a categorical `first`; dimensions without one count 0), then the
earlier trial. That trial is the leader (its `hyperparameters`, `score`, `se`, `foldScores`; `trials` counts all of
the candidate's, pruned and failed ones too); leaders best first. `chooseLeader` takes, across the leaders, the
first by `isRankedBefore`: the best `cautiousScore` (`score − direction · se`: the highest score − SE, the lowest
score + SE for a minimized metric) as a frame stores it (`Math.fround`), a tie by the method's name in code-unit order
(capitals first), exactly the order the Search bar chart sorts its bars in (measured 2026-11-10; the user's decision of
that day, replacing the one-SE rule between methods, so the chosen one is always the top bar and the grid's first row)
and marks it `chosen`; with one candidate (Tune)
its leader, the one-SE rule's within the method, is chosen. A continued search chooses again over all its trials. `searchRecordFrom` makes the `SearchRecord` saved with the run and the model (`search`, see "JSON column
shapes"); for Tune it keeps `space: {searched, pinned}` (the candidate's dimensions as searched and its pinned
settings). `searchRecordOf` narrows a stored value (undefined for anything malformed; a malformed `space` alone is
dropped, as is a missing one). `tuneSearchPlan(request, engine, start, choice, seed, space?)` gives Tune's
`PlannedSearch {plan, holdOut?}` (`holdOutSplit`, `searchPlanOf('tune', ...)` with `tuneCandidate`) and
`bestSearchPlan(request, candidates, choice, seed)` Find best model's; `tuneSearch(request, engine, start, choice,
seed, progress, executor, space?)` and `bestSearch(request, candidates, choice, seed, progress, executor)` run them
(`runSearch`) for the view and the API and make the record; `failedSearchMessage` gives `No trial completed: <the
first failed trial's error>` when a trial failed. Both return a `SearchResult {plan, holdOut?, outcome, leaders,
chosen?, record, state}`; `state` (`SearchState {plan, holdOut?, outcome, fingerprint}`, in memory: the plan holds
the columns) is what continuing needs. `canContinue(previous, plan)` (pure) accepts a plan of the same mode, strategy, seed,
metric and validation, the same candidates in order (each its engine name, space, pinned names and its start's values
outside the space: Tune's searched values may differ, as the form shows the last winner) and the same search rows
(`requestFingerprint`: the selection's feature and target names, the prepared columns' and target's
`datasetFingerprint` hash, the row count and the prepared options); the limit (its kind too), `topK` and the
parallelism may differ. `continueSearch(previous, planned, progress, executor)` refuses another plan (`The search
cannot be continued: the data, the validation or the search settings changed.`) and runs `planned` on top of
`previous.outcome`: the leaders, the choice and the record over all trials, the record's `budgetSeconds` and
`maxTrials` each the sum of the runs' of that kind (so are the state's `plan.budgetMs` and `plan.maxTrials`; a record
of runs of both kinds holds both), `trials`, `completed`, `pruned`, `failed` and `elapsedMs` cumulative.
`Forge:trainModel` does not continue searches.

**Find best model.** The candidates: `searchableEngines(applicable, discovery)` keeps the browser methods of a check
(`!isServerEngine`) sorted by `costRankOf` (the place in `COST_RANK`, from Linear Regression to RealMLP, the cost order
of fitting; an unknown method after it in discovery order); `quickEngines(engines, columns, target,
predictProbability)` asks each `retrainsLive` at once and drops a method whose probabilities are a train input
(`probabilityInput`: SVM) while Predict probability is on, as its fits are several times slower; at a time limit of
`QUICK_BUDGET_SECONDS` (60) or less, or with a trial limit, the default methods are the quick ones, above it every
searchable one.
`hasProbabilities(engine, on)` (`engines/engine.ts`) tells whether a method's trials carry class probabilities with
Predict probability `on` (always for `mlprobabilities: 'true'`, only when on for a probability input).
`bestCandidates(engines, discovery, predictProbability)` gives each method its full space (`searchSpaceOf`; a
malformed one becomes `note` and the method is tried with its defaults only), `start` = its defaults with the
probability input set from Predict probability (`withProbabilityInput`; the input is outside the visible space, so every
drawn point keeps it) and `costRank`, cheapest first; no method throws `No method can take part in the search. Check
the features and the target.` `bestSearch(request, candidates, choice, seed, progress, executor)` is `tuneSearch`'s
sibling with `searchPlanOf('best', ...)` (`topK` = `choice.topK ?? TOP_K`, 3). With `strategy: 'equal'` the
executor runs the candidates not exhausted one after another in plan (cost) order, each from its first point (its
defaults when inside its space): with `maxTrials`, a quota set when the candidate starts, `ceil((cap - launched) /
(candidates left))` (so the remainder goes to the first ones, and the sum reaches the cap while any candidate has
points: `trials` stops it); with `budgetMs`, a share that starts when the candidate does and lasts `(budgetMs -
elapsed) / (candidates left)`; so trials or time a candidate leaves unused (exhausted, finished, failed) go to the
ones after it, and a trial of it forecast past its share end is abandoned (never its candidate's first
complete one). A share also ends after `MAX_FAILED_WITHOUT_COMPLETE` own failures without a complete trial; the next
share starts only once the current candidate's trials settled (`nextTrial` returns null meanwhile, a wait). With
`strategy: 'halving'` the executor runs successive halving
(`Halving {round, pool, given, eliminated, finalists?, next}`, part of `resume`): round `r` gives each candidate of
the pool `2^r` trials, round 0 (the pool: every candidate in the plan's cost order) one each, its first point;
`nextTrial` serves the pool member with the fewest trials of the round below the quota (ties in pool order, so a
round goes round-robin), each its next untried point from its own stream (`drawOf` over its state: `trialPoints`, the
tried keys, its own counter); `given` counts a member's trials of the round (an abandoned one gives its slot back).
Once every member has its quota or is exhausted, nothing (`null`, the loop awaits a running trial) until the round's
trials have settled (`pending`); then `closeRound` ranks the pool members that have a non-empty space, are not
exhausted and have a complete trial by their best complete score so far (ties by index), keeps `min(n, max(topK,
ceil(n / 2)))` of them as the next round's pool, best first, records the round in `eliminated` for the others and
starts the next round. A member without a space (it gets its one start), exhausted, or without a complete trial (all
failed) leaves without being ranked or eliminated. A pool of `topK` or fewer becomes `finalists`, tuned in turn from
`next` as before (a candidate whose stream is exhausted leaving the cycle). A null while a round's trial is in flight
is a wait, never exhaustion (the loop stops with `exhausted` only when nothing runs). One candidate (Tune) draws in
turn and never waits. The budget, the abandon rule (every trial, round 0's too, once any trial completed), `maxTrials`
(the trials of all methods launched) and `MAX_FAILED_WITHOUT_COMPLETE` are global and stop the search in any round;
the choice is then made over the leaders as always. The trials are scored at the prepared request's threshold as
Tune's; the needs-probabilities failure names the trial's method. The order inside a round is fixed, a round is
closed only after all its trials settled and each candidate draws and prunes on its own, so a parallel search gives
the trials and leaderboard of a sequential one. The record's `methods` (`best` only, `candidateOutcomes(plan,
outcome)`): per candidate in cost order `{method, trials, outcome, note?, eliminatedRound?}`, `complete` when one of
its trials completed, `failed` when none did and one failed (`note` its first error), `abandoned` when its trials were
abandoned or pruned only, `unstarted` without a trial; a malformed-space note is kept otherwise; `eliminatedRound` the
round it was eliminated after (from `resume.halving.eliminated`), absent for a finalist or a candidate that left
otherwise; `searchRecordOf` drops a malformed `methods` list alone.

**Settings and Search options.** `fillMethod` builds a `TuneSettings` (`ui/tune-settings.ts`, `form.tune`) per method
over the method's own hyperparameter inputs. **Settings** is a `ui.input.radio` (`forge-inline-radio`; Manual / Auto-tune,
Manual by default; a change runs no check and no training: it shows the mode's layout, saves, and calls the view's
`onChanged(true)`: `applyRules`, `showTrainRow`, which shows **Train** in Manual only, and `updateControls`); the reset
button (`resetButton`, a native `button` with `fa-undo`) sits in its options, the status row (`forge-tune-status`,
30 px, a label column and the italic text, `aria-live`) under it. The status is computed from the values (`status`):
the view's `lastTuned` of the method (the values a tune set with its `TuneContext`, the table, the target's and the
features' names, and `isFound` for Find best model's, in memory like `hyperparameterValues`) when the host's
`context()` is the same (`isSameContext`: the table by `dart`, the target by name, the features as a set of names; the
preparation and the validation are not part of it) -> `Auto-tuned`, unless the values are the method's defaults
(`isDefault`, their text kept from the constructor) and a Tune set them; else the defaults -> `Default`, else
`Customized`, the values compared as text over every setting. A Find winner thus reads `Auto-tuned` even on its
defaults; the reset button follows `isDefault`, not the status. The status's tooltip (`statusTooltip`) is its
`STATUS_TOOLTIPS` entry (`Default`: `The method's recommended settings`, `Auto-tuned`: `Found by the search`,
`Customized`: `Changed by you`, `Not searched`: `Click Find to search`; no search lines are appended, those stay in
the Runs and Search panels), except for the tuned values on other data: `Tuned for <target> on <table name>; tune again for
this data`. Each hyperparameter change calls `updateStatus` from the settings' own `onChanged` subscription, a Target or
Features change from the view's form subscription (a Table change builds new settings)
(the status read once: text, the reset button's `aria-disabled` and label, the same as its tooltip: `<Method> settings
are already default` / `Reset <Method> settings to defaults`). Reset sets each own input to its default at once under
`isResetting` (the settings' own error and status pass runs once at the end), so the usual change flow (check, live
retrain) follows. A method without a space (`noSpaceReason`: `<Method> declares no search space.`, or the malformed
space's error) has its Auto-tune option (`autoItem`, found once) disabled with that tooltip, at construction and on
every `setLocked` (the lock owns the mode: `mode.enabled` re-enables its options); the Settings label and the other options share `How the settings
are chosen: by hand or by a search on this data`. A method without settings (no hyperparameter inputs: Naive Bayes)
lays out no Settings row at all (no radio, reset or status) and is always Manual (`load` refuses a stored Auto-tune
for it), so its Train shows as for any method. Both layouts are built once in the constructor (the radio, the status
row, the own inputs, the Search options block `forge-tune-auto` with one `forge-tune-cell` per row: Metric, Limit,
Max time and Max trials first, then the settings' rows in the method's order; and the Tune bar
`forge-tune-bar`, `tuneButton` at the right); a mode switch only shows Manual's own inputs or Auto-tune's block and
Tune bar (`showMode`, `ui.setDisplay`). The Tune bar stays in the document in Manual, so a disabled Tune's tooltip
overlay lives on in Auto-tune; in Manual the view releases Tune (`ButtonGate.release`), so a hidden Tune keeps no
overlay timer. With the labels above the inputs the
status row drops its empty label and the toggle has `margin-top: 10.5px`, so the status sits right under the radios and
the gaps radios -> status -> Search options stay equal. The radios' editor drops the padding the platform gives it
for its options (the reset button follows the radios in the row), so Auto-tune stays on the radios' line.

**The Search options block is plain DOM, by the user's decision** (2026-11-04): platform inputs there were laid out
wrongly by the platform form (it measures every `ui-input-root` for its label width and its switch to labels above the
inputs, and nested inputs broke both), so the block follows the reference prototype
with native controls in its own CSS grid (`forge-tune-grid`: the label column, then `minmax(0, 1fr)`, so nothing
overflows at any width). The label column is the form's: `TuneSettings.alignLabels` copies the width the platform
gave the Settings label (`captionLabel`, followed by `ui.onSizeChanged`; Find best model's block follows **Method**'s)
into `--forge-tune-label`, with the labels' 8 px gap. The toggle (`optionsToggle`) and the Metric row (`metricRow`,
`setMetrics`) are static builders that `FindSettings` shares, and so are the limit rows (`LimitInputs`,
`ui/limit-inputs.ts`): **Limit** (`ui.input.radio`, `Time` / `Trials`, `forge-inline-radio`), then **Max time**
(`ui.input.int` in whole minutes, `min: 1`, `max: 60`, `nullable: false`, `step: 1`, `showPlusMinus`,
`addPostfix('min')`, the platform's unit) or **Max trials** (`ui.input.int`, `min: 1`, `max: 10000`, `nullable:
false`, `showPlusMinus: false`, since a range would add the clicker); the platform validates both (an integer, the range, not
empty: `'abc' is not an integer.`, `Value must be greater than 1`, `Value must be less than 60`, `Can't be empty`), marks
an invalid editor `d4-invalid` (a red underline) and shows its messages in the input's tooltip under its description
(`tooltipText`; a `ui.tooltip.bind` on the editor would hide them). Find's own rule is one `addValidator` on Max
trials, run only on a count in range; these platform inputs (and Find's Strategy and Candidates) sit in grid cells with
`forge-tune-input`: their own labels hidden, the caption a `forge-tune-label` in the grid's label column
(`LimitInputs.row`), an int editor as wide as its cell (ui.css caps it at 150 px) and at the bottom of the 30 px row
(`align-self: flex-end`): ui.css puts the options (clicker and unit) there, each with its own underline, so the
editor at the top drew a second line 2 px above theirs. The platform form sets its label width a while after such an input shows, so a test
that measures the labels waits for them to settle; when the
platform switches the form to `ui-form-condensed` (a form narrower than its labels and widest input: the Train form at
400 px sits near that line, so a resize such as opening Search options can switch it), the grid puts its labels above
the fields too. The toggle (`toggle`, `forge-tune-toggle`, a `button`: an empty `forge-group-chevron` box and the
groups' header gap first, so its chevron, `chevron-right` / `chevron-down` swapped as in `CollapsibleGroup`, starts
where the Method header's text does; `aria-expanded` / `aria-controls`, collapsed by default, saved, as far under the
status line as the status is under Settings) shows or hides the grid (`setExpanded`). Every input of a row shows the
setting's Manual tooltip exactly, its `Property.description` (nothing without one): the label, the range fields, the
pinned / fixed value, the choices button and its checkboxes. A pin's tooltip is `Pin to exclude from tuning. <its
aria-label>` (a fixed one `<Caption>: Fixed, can’t be tuned`); the Metric label and select `Metric the search
optimizes`, Limit (label and options) `What ends the search: a time or a trial limit`, Max time (label and field)
`Minutes the search may run`, Max trials `Maximum number of trials, pruned and failed ones included`, the toggle `Show / Hide search options`. What
the user sees ends on one vertical line, the block's right edge: the reset and pin glyphs (their 30 px buttons put the
glyph at the right, `justify-content: flex-end`, the hit area extending left), the selects' underlines with their
chevrons, and Tune's letters (no right margin or padding). The toggle (`margin-top: -3px`) puts the Search options text
as far under the status text as the status text is under the Settings label; the layout test measures the text
through DOM Ranges and the glyph boxes, beside and above the inputs. The fields, selects and choices button share
`forge-tune-control` (the underline, focus, `aria-invalid` and disabled rules, scoped to the grid so they beat the
platform's `select.ui-input-editor`). The selects are `ui-input-editor`s with the platform choice input's chevron and
no focus frame; focused, they underline in `--blue-1`. The chevron rule (`url('/images/chevron.svg')`) is in
`css/forge-platform-images.css`, listed in `package.json` `sources`: the platform loads it with the package, unbundled,
so the path stays the platform's (the bundled `forge.css` would have the build resolve it). `togglePin` swaps only its
row's slot and pin state, and hides the tooltip first: the pin's tooltip text changes with the click, and a hovered
element that leaves the document gets no `mouseleave` (the stale `Rate` over Max Depth). Per setting in
the train order (`TuneRow`, `row(name)`): a label, then on one line (`forge-tune-row`) the value and the pin, and the
error line (`forge-tune-err`, `role="alert"`) under them. A tuned numeric setting is its range: two text fields
(`forge-tune-field`, `from` / `to`, `aria-label` `<Caption> from` / `to`) around a dash; a tuned categorical one a
`forge-tune-select`-looking button listing the checked choices that opens a list of checkboxes (`choices`; closed by a
click outside or Escape); a pinned one, or one outside the space, its own value (`value`: a native `select` for a
choice setting, a checkbox for a yes/no, a field otherwise) kept in step with the setting's input both ways
(`own.stringValue` as typed, then `own.validate()`; `own.onChanged` writes it back), so its value is the setting's.
Each line ends with its pin, a native `button` (`forge-tune-icon forge-tune-pin`, `fa-thumbtack`: `fal` `--grey-4`
when tuned, `fas` `--blue-1` with `forge-tune-pinned` when pinned, `fas` `--grey-3` and `aria-disabled` outside the
space; `aria-pressed`, `aria-label` `<Caption>: will be tuned` / `<Caption>: pinned,
excluded from tuning` / `<Caption>: Fixed, can’t be tuned`). Before them, first in the grid, **Metric** (a native `select` of the
task's `SEARCH_METRICS`, values the metric ids, set by `setTask` from `applyRules`; the chosen one kept while it
applies, else the stored one, else the default) and the limit rows (`limit`; only the input Limit chooses is
shown; both are validated, only the shown one counts). The colours are the platform tokens (`--text-color`, `--text-color-light`,
`--grey-1..5`, `--blue-1`, `--red-3`), 13 / 12 px. `updateErrors` (after every change, in Auto-tune only; entering
Auto-tune runs it) writes each row's problem (`rowProblem`: the tuned range by `rangeProblem`, `Choose at least one`,
or a pinned / fixed value's own `validity`; none for a setting the method ignores), keeps it on the row, validates
the limit (`LimitInputs.validate`: `validate()` of both inputs, the shown one's result; `problem()` is the shown
input's first `validationResults` message captioned, `Max trials: Value must be greater than 1`), marks the fields
`aria-invalid`, and opens the grid on a problem. `setInactive` (one set, from
`applyRules`, which the view's `activeHyperparameters` reads too through `isActive`) disables the value controls of a
setting the method ignores, so it is not checked; `setLocked` disables every control and the mode. In Auto-tune
`applyRules` evaluates the rules with the tuned settings free and the target's task. `TuneSettings` reads
`visibleHyperparametersOf(engine)`, so SVM's hidden `probability` neither shows nor makes it read `Customized`.
The view's `applyRules` runs only when the values, the task, the mode or the pins change (`onChanged(affectsRules)`),
not on a range, the choices, Metric or the limit, and sets an input's `enabled` only when it differs; the rules are read
once per method (`fillMethod`). `problem()` is the Tune part of `tuneBlocker`, from the rows' problems as
`updateErrors` left them: `Fix the highlighted range`, `All settings are fixed` (no tuned setting the method uses),
the limit's (`Max time: <problem>` / `Max trials: <problem>`). `choice()` is the `SearchChoice` (metric, the one limit:
`budgetSeconds` = the minutes · 60 or `maxTrials`, `space` = pins and ranges). Everything but the values is saved after
every change, debounced (`TuneSettings.saver`, see the Find block), in `grok.userSettings` (storage `TuneSettings.storage`, `Forge Tune`; the UI tests set their own, `Forge
Tune tests`, and empty it after each category, so the user's settings stay as they are; key = method name, JSON
`{mode, expanded, pinned, ranges, metric, limit, minutes, trials}`, an invalid value saved as its default);
`TuneSettings.load` checks each item (a pin or a range of a setting no longer in the space, a range `rangeProblem`
refuses, a choice subset that is empty, a mode the method cannot have) and falls back to its default;
`LimitInputs.load` reads the limit: the kind (`Time` or `Trials` as stored; an old record's `seconds` is `Time`; else
the default, `Trials`), the minutes (1..60) or else an old record's `seconds` as whole minutes rounded up and clamped to
1..60, the trials (1..10000), else the defaults (Trials, 1 min, 50 trials; Find best model's 100); a metric of the other task falls back in `setTask`.

**In the view.** The **Tune** button (`startSearch('tune', form.tune.choice())`, gated by `tuneBlocker`: the
selection's problem, a training, `noSpaceReason`, or `TuneSettings.problem()`; a trained selection can be tuned; its
tooltip `searchTooltip('tune')`: `Find the best settings`, or `Continue the search: <N> trials so far` while a kept
search would be continued). `startSearch(mode, choice)` runs through the view's `TrainingQueue` like a training (it counts in `activeTrainings`, so **Train**,
**Save** and `aria-busy` behave as during one; `SHOWN_RUN` goes with the hint), on the selection as the form has it
(`searchSelection`: Predict probability and its threshold included) and the form's missing-values setting: `prepareTraining`, the `datasetFingerprint` of all its rows,
the search with the queue's progress (`Tuning <Method>`, cancelable, closed in a `finally` before the winner's
training starts). **Continuing.** The view keeps the last `SearchState` of each mode (`searches`, keyed by
`searchKey`: `best`, or `tune <Method>`; in memory, as the plan holds the columns; a Table change clears it) with the
`searchInputs` it ran on (JSON of the target, the feature names, the preparation steps with the threshold and its
metric, as `requestFingerprint` keeps them in the options, `selectedMissingValues`, the validation, the seed, the
metric, the strategy, the checked methods, Tune's space and the values of the settings not tuned; not the limit), so
the tooltip promises a continuation only for a search `canContinue` accepts (a Predict-probability threshold moved by
the winner's training starts the next search afresh). `startSearch` starts nothing while `tuneBlocker` / `findBlocker`
holds (an invalid limit included), whatever `choice` it is given. `keptSearch` gives the state while the inputs are the
same; the run plans the search once (`tuneSearchPlan` / `bestSearchPlan`) and calls `continueSearch` when
`canContinue(state, plan)`, the progress then reading `Continuing: <title>`, else `runPlanned` starts it afresh. Either way the result's `state` replaces the kept one for its key. The leaders, the trials
and the record are cumulative (the record's limits the sums), so **Search** shows them all. The split seed is the form's **Split seed**. Until the search and the winner's training end, the
form is locked: `lockForm()` keeps each form input's `enabled` (Table, Target, Features, the Preparation and
Validation inputs, Method and the hyperparameters) in `locked` and disables it, and `TuneSettings.setLocked(true)`
disables Settings and every Search options control (`controls`: the reset button, the fields,
choices, values and pins, Metric; and the limit's inputs; not the toggle, which only shows and hides, as the group headers) until
`setLocked(false)` (`fillMethod` under the lock locks the new settings too), **Train** and **Tune** wait as during a
training, `restoreRun` / `restoreLeader` (the **Runs** and **Search** row clicks) set nothing, and `applyRules` under
the lock records a rebuilt input's state (`setWinner`) instead of enabling it; `unlockForm()` in the `finally` gives
every input its state back and applies the rules again, whatever the end (completed, failed, cancelled, an error).
Cancel stays in the task bar. `train()` (the **Train** button) takes the same lock around its `startTraining(false)`
and gives it back in its `finally`; a live retrain (`startTraining(true)`) does not lock the form. Locks never nest:
both callers start only while nothing trains (`isTraining`, `tuneBlocker`). Then
`SearchPanel.update({plan, leaders, outcome})`; without a chosen leader the warning `No trial completed: <error>` when a trial
failed, else `No trial completed.`, and the hint (`idleHint`); else `setWinner` (`setMethod`: the method's values set to the chosen ones,
the inputs rebuilt by `fillMethod` and **Settings** to Manual; **Method** set under `isSettingMethod`: no check is
requested), `lastTuned` of the method = the form's values with the search's table, target and features (the status
then reads `Auto-tuned`), the balloon (`Optimal settings for <Method> found.` or
`Search cancelled; the best so far is trained.`) and `startTraining(false, {selection, request, fingerprint})`: the
override skips `prepareTraining`, records the fingerprint of all the selection's rows (as the API) and trains the search request with the chosen hyperparameters, `holdOut` and `search`, with the
search's seed and split, so the model's Validation equals the leader's score. `trainModel` measures `metrics.test` on
the hold out and returns `search`; the run record, the `SessionRun` (`Found by Tune: <N> trials` in the Runs Method
tooltip), the Performance bullets and the saved model carry it. `SearchPanel` (`forge-results-panel forge-search`,
a bottom block of the results pane, see below; its title's tooltip `Each method's best settings, ranked by <Metric> −
SE.` once a search is shown, through `resultsTitle`'s `tooltip`) shows only a search of two rows or more, so never
after Tune or a Find of one method (`update` hides it; Performance and the Runs tooltip describe such a search). It is
a `readOnlyGrid`, its leaders sorted by `isRankedBefore`, as
the chart sorts its bars (the bar chart ignores a column's category order, measured 2026-11-10), so the chosen one is
the first row and the top bar; its own copy:
`result.search.leaders` keeps its order for the record, Runs and the cards), of
`#` (header `Rank by <Metric> − SE.`), `Method`, `<Metric> (CV)` (`<Metric>
(hold-out)` when the search request's validation is a hold out), `SE`, `<Metric> − SE` (`cautiousScore`, so
`<Metric> + SE` for a minimized metric, empty without a score; header `<Metric> minus its standard error over the
folds: a cautious estimate of <Metric>.` (`plus` for a minimized one; `(0 on a hold out: one score)` for a hold out),
then ` The method with the highest value is chosen.` (`lowest`)), `Trials` (header: `Every trial
of the method, pruned and failed ones included.`), `Chosen` with header and cell tooltips
(the settings, then `Searched: <name> <low>–<high>, <name> (<choices>), ...` and `Pinned: ...` from the plan's
candidate; `Folds: ...`; Chosen: `The method with the highest <Metric> − SE.` / `Chosen: the highest
<Metric> − SE.`), the chosen row current; with the `outcome`, the
candidates without a leader (`candidateOutcomes`) follow in candidate order with an empty value and SE, their trials
and the Method tooltip `No completed trial: <note>` / `Abandoned: the budget ran out` / `Not started: the budget ran
out` above `spaceLines`; a method the halving dropped (`CandidateRecord.eliminatedRound` r) reads `Eliminated in round
<r + 1>` as its Trials tooltip and in its Method tooltip, one without a space (which leaves the pool after round 0,
`closeRound`) `No settings to tune` as its Trials tooltip; under an equal split the Trials tooltip explains the count
(`Its space was exhausted` from `outcome.resume.exhausted`, else `Its share of the limit ran out`; nothing for a
failed method). With a leader, a `DG.Viewer.barChart` over the grid's own frame (so the current row and the selection
are shared) sits left of the grid in a `ui.splitH` (`forge-results-row`; each viewer in a `forge-panel-fill` box, as
the splitter wraps anything but a box and a grid root as the pane would keep its own height; preset widths 200 / 100
px, then growth 2 / 1: two thirds for the chart). The chart: split `Method`, value `<Metric> − SE`, `valueAggrType`
sum (one row per bar), `filter` `IsNotEmpty(${<Metric> − SE})` (the methods without a score have no bar), sorted by
value, best first (`desc`, `asc` for a minimized metric), the value labels on, no title, the
`description` `<Metric> minus its standard error: a cautious estimate` (`plus` for a minimized one;
`descriptionVisibilityMode: 'Always'`), no formula line, stack, color column, legend, value, category or stack selector, `showSelectedRows` and `showFilteredRows` off (one bar color).
`onClick: 'Select'`: a bar click selects its method's row, and the panel restores the method of a one-row selection
(`frame.onSelectionChanged`), as a row click does (the bar chart fires `onCategoryClicked` only for a stack segment).
A click calls `restoreLeader` (`onRestore` takes `Pick<Leader, 'method' | 'hyperparameters'>`:
a candidate row passes its `start`; `setMethod`: the method and its values set, **Settings** to Manual, then a check,
as after any change; a `restoreRun` switches to Manual first, whatever it then sets). `createForm` (a Table change)
and `detach` hide it. **Search** and **Runs** (`forge-results-panel`, their grid hosts `forge-panel-fill`) are the
`ResultsPane`'s bottom blocks (its constructor's `bottom`, each in its own box, 320 / 260 px at first): `layout()`
mounts the result rows and the shown bottom blocks in ONE `ui.splitV` (a single box as is), rebuilt only when the
sequence of boxes changes; the view calls it after every `searchPanel.update` / `runsPanel.update`. A hidden block's
box is taken out before the old splitter is emptied, so it leaves no pane or divider and keeps its widgets; the size
keeper notes and restores the bottom boxes' heights as the rows'.

**Find best model in the view.** `AUTO` (`Auto`) is the first **Method** item (`createForm`, `showMethods`);
`userChoice` holds it while chosen (`isAutoMethod`), and `this.engine` stays a real method, the suggestion (`methodOf`
returns `best ?? engines[0]` under Auto, with no balloon), so the hidden settings belong to it. `chooseMethod(AUTO)`
refills the Method group on entering Auto and requests a check with `requestCheck(true)`, which keeps `lastTraining`
(and `isTrained` / `isFailed`): the selection is the same, so the results pane keeps the model's viewers, metrics and
insights, and nothing trains; a method name leaves it. `fillMethod` appends `form.method`, `form.tune` and
`form.find` and shows the tune root or, under Auto, the find root; under Auto **Missing values** never hides
(`setKeepsMissingValues(false)`), Predict probability follows `probabilityReason(form)` (`form.find.probabilityReason`
under Auto, else whether the method gives probabilities, through `PreparationInputs.setProbabilityReason`) and the Method tooltip is `AUTO_TOOLTIP`. `form.find` is one `FindSettings` per form
(`ui/find-settings.ts`), laid out as Auto-tune's block from the same builders (requirement: switching between a method
and Auto changes nothing visible but the rows): `TuneSettings.statusLine()` reading `Not searched` (`NOT_SEARCHED`,
tooltip from `STATUS_TOOLTIPS`) right under **Method**, then `forge-tune-auto` with the `optionsToggle` (collapsed at
first; `setExpanded` saves the state, `expanded` in the stored settings, restored; an invalid limit opens it;
the Find bar stays shown while it is collapsed) and the `forge-tune-grid` (label column aligned to **Method**'s caption through `alignLabels`, its
subscription in `find.subs`, which the view keeps in `formSubs`): `metricRow`, **Strategy** (`ui.input.radio` of
`STRATEGY_LABELS`, `Equal` / `Leaders` (short, so the row fits the 400 px form on one line; the printed name of
halving stays `limitCaption`'s `Focus on leaders`), Equal by default; option tooltips `Every method gets
an equal share of the limit: a fair comparison` / `Each round the worse half drops out; the rest get more trials`,
label `How the checked methods share the limit`), the limit rows (`LimitInputs`: 100 trials by default, Max time 1
min; its `trialsProblem` validator of Max trials, under Equal, `At least one trial per checked method: <N>` while the
trials are fewer than the checked methods, re-validated through `updateErrors` on every change of the methods
(`rebuild` ends with it) or the strategy) and **Candidates**,
a `forge-tune-label` and the platform `ui.input.multiChoice<string>('Candidates')` (`forge-tune-input
forge-find-candidates`: its own label hidden, every check shown, the first level with the caption), then the
`forge-tune-bar` with **Find** (`findButton`, gated by `findGate` / `findBlocker`: the selection's problem, a training,
or `FindSettings.problem()`: `Check at least one method`, else the limit's (`Max time: ...`, `Max trials: ...`); tooltip
`searchTooltip('best')`: `Find the best model`, or `Continue the search: <N> trials so far`).
`setMethods(entries)` takes `MethodEntry {name, isLive, isServer, probabilities: always | input | none}`; `rebuild`
lists only the entries without a reason (`reasonOf`: a server method, `<Method>: Runs on the server; not searched in
the browser`; one with `none` while Predict probability is on, `<Method> gives no class probabilities`; one without
probabilities now under a score metric, `<Method> gives no class probabilities, which <Metric> needs`: `exclusionOf`,
shared with `Forge:trainModel`), the
reasons going to the Candidates caption's tooltip as `Not listed: <reason>` lines after `METHODS_TOOLTIP`, and checks
them by the user's toggle (`overrides`, kept only against the default) or the default (not a quick limit,
`isQuickLimit`, or a quick method, `isQuickMethod`: live and not a probability input while Predict probability is on);
it clears the value first (the platform keeps a checked item it no longer lists) and ignores its own value changes
(`isRebuilding`). A change of the metric rebuilds, validates, saves and tells the view; one of the limit rebuilds only
when the limit becomes quick or stops being quick (`isQuick`), one of the strategy only validates, saves and tells it.
Stored in `grok.userSettings` (`TuneSettings.storage`, key `Find best` = `FIND_KEY`; Tune's and Find's saves go
through `TuneSettings.saver`, written 300 ms after the last change into the storage of that change, and at once by
`flushSave`, which the settings' `subs` call when the view drops them): `{expanded, metric, strategy, limit,
minutes, trials, methods}`, each item checked by `FindSettings.load` over `TuneSettings.storedRecord` (the limit by
`LimitInputs.load`, an old `seconds` migrated as Tune's). The Candidates caption's `METHODS_TOOLTIP`: `The methods to try, each with its
defaults first`. `inputs` (Strategy, Limit, Max time, Max trials, Candidates) join `formInputs()` (the lock) and
`expandOnError`; `setLocked` disables them and Metric. `revalidate` under Auto clears `problems.method` while
a method can learn, asks the browser methods `quickEngines(..., false)` at once (into `interactivity`), lists them
(cost order) then the server ones, sets Predict probability's reason, validates the find inputs and, instead of
training, shows `AUTO_HINT` only without a `lastTraining` (`idleHint` under Auto shows the kept training again after a
search that trained nothing); `showTrainRow` hides **Train**, `problem()` skips the hidden settings, `updateControls` gates
**Find** and releases **Tune** (and the other way round outside Auto); `onProbabilityToggled` only tells `form.find`
(`setPredictProbability`). `startSearch('best', form.find.choice())` (`choice()`: metric, strategy, the one limit,
`methods` the checked names; the plan's `topK` is `TOP_K`) runs the same flow as Tune with the selection's missing
values by `searchMissingValues` (`impute` kept, else `skip`: every method gets the same folds), `bestSearchPlan` over
`bestCandidates(candidateEngines(choice.methods), ...)` (the last check's browser methods so named) and the progress
`Finding best model`; the winner (`plan.candidates`' engine of the chosen leader) goes to `setWinner` (Method
leaves Auto), `lastTuned`, the balloon `Best model: <Method> (<Metric> 0.000).` and `startTraining` with `engine:
winner` in both the selection and the request (what `trainModel` trains and `trainingRunOf` records).

**API.** `Forge:trainModel` with `options.search` (`searchChoiceOf`: keys `metric`, ONE limit `budgetSeconds` (whole,
>= 1) or `maxTrials` (whole, >= 1), both refused with `COMBINED_LIMITS`, neither: `budgetSeconds` 60; with `method:
'auto'` only `strategy` (`equal` or `halving`, else `The search option 'strategy' needs equal or halving.`), `methods`
(a list of names) and `topK` (whole, >= 1); any of the three with a method: `The search options 'strategy', 'methods'
and 'topK' need method auto.`; `searchPlanOf` refuses an equal split by fewer trials than methods) tunes the chosen or suggested method in `trainAndSave` before its training,
with Predict probability as the call sets it, as in the view (the trials score the labels its threshold gives, and
the winner's threshold is optimized for `thresholdMetric` unless `threshold` was given; the selection's `hyperparameters` carry the method's probability input): `searched` runs the same `tuneSearch` (or, with `method: 'auto'`,
`bestSearch` over the candidates of `autoEngines`, see "`Forge:trainModel`") with `browserSearch` and the call's progress
(`searchTitle`: `Tuning <Method>` or `Finding best model`, then a new `Training <Method> model`), and trains the request with the
winner's method (`request.engine`, which also names the saved model and the result's `method`), the chosen settings,
the hold out and the record; no completed trial throws `No trial completed: <the first failed trial's error>` when a trial
failed, else `No trial completed within the budget. Give the search more time.` (or `The search was cancelled before a
trial completed.` first), and records nothing. `options.seed` (42 by default) is the one seed (the split, the
search's hold out, sampling); `options.scheme`, `folds`, `testFraction`, `stratified` (ignored for a regression) and
`shuffle` make the `ValidationSettings`, used by the search as by the training. The run's `duration_ms` includes the
search.

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
one-hot, predict probability. `PreparationSteps` is `{oneHot, skipUniqueCategories, predictProbability, threshold,
thresholdMetric?}` (`DEFAULT_THRESHOLD` 0.5; `thresholdProblem` refuses a value outside 0..1; `thresholdMetric`, a
`BestThresholdMetric` with `DEFAULT_THRESHOLD_METRIC` `f1`, asks for the threshold that maximizes it, resolved by
`recutResult`; absent, `threshold` is used as it is).

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
- **Predict probability**, for a categorical target with two or more classes (`classesOf`, the empty one aside): the
  target is not changed (the method trains the classifier it is); `postprocessingInfo` gains `class-probabilities`
  (`CLASS_PROBABILITIES`), and for two classes `positiveClass` (the first class), `negativeClass` and
  `binaryClassificationThreshold` (the threshold) are recorded, the keys of the built-in tool's models, so
  `preparationOptionsOf` reads them as before, plus Forge's `thresholdMetric` when the steps carry one (`recutResult`
  keeps it only when it found the best; `preparationOptionsOf` keeps a known metric id only). Since 0.1.20 Forge never records `binary-classification`, the former
  regression mode (a float 0/1 target, `targetType`); its models still replay. The labels come from the method's class
  probabilities in `replayPostprocessing` (below); a result without them drops the step (`recutResult`, see
  "Training").

`replayPreprocessing(columns, options)` (application) replays `options.preprocessingInfo` on the mapped columns and
returns the columns the method gets (the given list itself when there is nothing to replay): `ignore-missing` and
`impute-missing` are skipped, any other unknown id is refused ("The model uses the preparation step 'X', which Forge
cannot replay yet."). With `oneHotCategories`, `one-hot` builds exactly the recorded columns (a value training never
saw is 0 in each, a category the table lacks is a column of zeros); with `skippedColumns`, `skip-unique-categories`
drops exactly those columns by name, whatever their values (the application does not map them, see "Application
pipeline", so they are there only when passed directly); without it (the built-in tool's models) it drops the
all-unique columns (the core's rule); without `oneHotCategories` one-hot uses the applied data's categories.
`replayPostprocessing(prediction, options, probabilities?)` (also in `pipeline.ts`) is the one place that turns a
method's output into labels, for training (`recutResult`, `trainPredictionOf`), the search trials, Apply and
Evaluation: `binary-classification` cuts the raw score into the two classes (`binaryClasses`, see "Application
pipeline"); `class-probabilities` takes the labels from `probabilities` (`probabilityLabels`): without them the
prediction itself; with the two-class keys and the positive class's column, `positiveClass` where its probability is
at least `binaryClassificationThreshold`, else `negativeClass`; otherwise the most probable class (ties: the first
column). A row without a prediction stays empty, a row without probabilities keeps the prediction's label; the column
keeps the prediction's type and name and is compacted (no stray empty category).

## Metrics

`metrics/metrics.ts`, ported from the built-in tool's definitions; labels and descriptions are the built-in ones (but
`R2` for the built-in "R squared", `Recall (Sensitivity)` for "Sensitivity" and `Negative Predictive Value` for
"Negative Predicted Value"), values are shown with 3 decimals.

- Regression: `mse` = mean of squared errors, `rmse` = its root, `mae` = mean of absolute errors (new),
  `r2` ("R2") = 1 - SSres / SStot; for a constant target 1 if the predictions are exact, else 0.
- Classification: labels compare as text (a boolean `true` equals the label `'true'`). `accuracy` = correct / n.
  When at most two labels occur, the shares are computed for the positive class, the target's first category:
  `sensitivity` = TP / (TP + FN), `specificity` = TN / (TN + FP), `precision` = TP / (TP + FP),
  `npv` ("Negative Predictive Value") = TN / (TN + FN); a share whose denominator is 0 is left out, as the built-in
  tool hid it. `f1` (new): of the positive class for two classes, the macro average over the classes otherwise.
  For a two-class target (`classesOf`: its empty category aside) the positive class, the first class, is saved as
  `metrics.positiveClass`.
- Rows where the actual or the predicted value is empty are skipped; no rows left is a `ForgeError`.
- `auc` ("AUC-ROC"): `aucOf(actual, score, positiveClass)`, the trapezoid over the rows
  sorted by score, highest first, the rows of one score taken as one step (one diagonal segment, as scikit-learn; the
  built-in ROC curve stepped row by row, so tied scores counted in their row order); undefined (left out) when the
  rows hold one class only. For any classifier with probabilities, `metricsOf(task, actual, predicted,
  positiveClass?, probabilities?)` adds `auc` from `aucOfProbabilities(actual, probabilities, positiveClass?)`: the
  positive class's column for two classes (`probabilityOf`), the macro average of the
  one-vs-rest AUCs otherwise (scikit-learn's `multi_class='ovr', average='macro'`) over the classes whose rows hold
  both classes; the train value on the final model's probabilities, the validation value on the out-of-fold columns
  gathered by class name (a class a fold model never saw is 0 in its rows). With probabilities `metricsOf` also adds
  `probabilityMetrics(actual, probabilities, positiveClass?)` (TS, over the rows where the actual value and every
  probability column have a value, classes compared as text): `logLoss` ("Log loss") = the mean of -ln of the actual
  class's probability, clipped to [1e-15, 1 - 1e-15] (0 for a class without a column, so about 34.5 for that row);
  `brier` ("Brier score") = for two classes with a positive class the mean of (p_positive - y)^2, otherwise the mean
  over the rows of the sum over the columns of (p_k - y_k)^2. A Predict probability training measures the labels of
  the probabilities (`recutResult`), the probability metrics as without the step.
- `bestThreshold(actual, probability, positiveClass, metric)` (the optimized threshold, through `bestThresholdOf` in
  `training/train-model.ts`): the paired rows sorted by probability,
  highest first; each distinct probability is a candidate (`p >= t` positive) whose TP / FP accumulate (TN / FN from
  the totals) and whose metric of `BEST_THRESHOLD_METRICS` (`f1`, `balancedAccuracy`, `mcc`, `accuracy`;
  `BestThresholdMetric`) comes from the counts (`f1Of` of the positive class; the others from `countMetrics`, the
  per-class counts core `classificationMetrics` uses too, on the two classes' counts); the best candidate wins, a tie
  the one nearest 0.5; the value is the midpoint between it and the next lower distinct probability (the candidate
  itself when it is the lowest), rounded to 3 decimals when the rounding stays in that gap; the gap is checked with
  `isAtThreshold(p, t)` (`p >= Math.fround(t)`: single precision), the predicate the `class-probabilities` replay cuts
  the float32 probabilities with (so a typed 0.7 keeps a row stored as 0.699999988 positive); undefined when the rows
  hold one class only.
- `balancedAccuracy` ("Balanced accuracy": the mean recall over the actual classes) and `mcc` ("MCC": the multiclass
  formula of scikit-learn, which is the usual one for two classes) are `classificationMetrics`' too, from the same
  label counts, for any number of classes; `averagePrecision` ("Average precision") comes from Crux (below). They
  follow `auc` in `METRIC_IDS`, then `logLoss` and `brier`. `classLabelMetrics` gives each label's precision, sensitivity and F1 (the re-cut's
  per-class values); `classesOf(target)` the target's classes, its empty category aside. `CLASS_METRIC_IDS` (`precision`, `sensitivity`, `f1`, `auc`,
  `averagePrecision`) are the ids of the per-class values: `ClassMetrics` `{label, rows, train, validation}` of a
  training (`MetricsRecord.classes`), `ClassValues` `{label, rows, values}` of an evaluation. `Curves` holds the
  positive class's ROC (`fpr`, `tpr`, `threshold`) and PR (`precision`, `recall`, `threshold`) curves as float arrays.

**Crux metrics** (`metrics/report.ts`). `Eda:classificationReport` (found by `DG.Func.find`, `reportFunction()`) runs
the Crux (scikit-learn) metric functions; `classificationReport(actual, predicted, positiveClass?, probabilities?,
curves)` calls it through `func.apply` with the standalone `actual` and `predicted` columns and a shared frame of the
probability columns (`onFrame`; `null` without them), and resolves to the report (`ClassificationReport`, mirroring
EDA's type) or undefined when EDA has no such function or the call throws (Crux not loaded, a Crux error; caught
there, not logged). `measure(task, actual, predicted, positiveClass?,
probabilities?, curves, score?)` = `metricsOf(...)` for a regression; for a classifier, `metricsOf` of the labels plus
`probabilityMetrics` (with probabilities), then `auc` and `averagePrecision` (with probabilities, or with `score`, the
positive class's score of a model of the former regression step, which the report gets as that class's only
probability) are the report's, and the report itself is returned (the positive class defaults to the target's first
class, `classesOf`); without the report `auc` is `aucOfProbabilities`', as in `metricsOf` (`aucOf` of the `score`). `SCORE_METRIC_IDS` (`auc`, `averagePrecision`, `logLoss`,
`brier`) are the values no threshold changes, which a re-cut keeps. The report gives the per-class values (`classMetricsOf(train, validation)` by label from both
splits' reports, a label of one report only with empty values for the other; `classValuesOf(report)` for an
evaluation; the report's recall is `sensitivity`) and the curves (`curvesOf`: float arrays, the PR's null last
threshold becomes Infinity; undefined without curves: three classes, a missing positive class, one class only).
Parity rule: the TS ids other than `auc` stay the source of their values (balanced accuracy and MCC included);
`auc` is the report's when there is one; the `Training` and
`Metrics` tests compare them with the report's at 1e-9. Without EDA's function (version skew), or when it fails, Forge
measures the TS ids only (no average precision), with no per-class values and no curve panels.

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
  `The model finds <negative> much better than <positive>.` or the other way round (detail `Recall (Sensitivity)
  0.900, Specificity 0.300.`, panel `confusions`);
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

`summaryInsight(stored, evaluation?)` is the one line the Results pane shows when there is no other: Balanced accuracy
for a classifier, R2 for a regressor (`METRIC_LABELS`), read from `evaluation` (`on the evaluated rows`, Run
Evaluation), else from the stored validation values when they have either metric (`on the validation rows`), else
from the train ones (`on the training rows`): `Balanced accuracy 0.850 on the validation rows: on average 85% of each
class is predicted correctly.`, `R2 0.800 on the training rows: the model explains 80% of the target's variation.`
(a negative R2: `the model predicts worse than the target's mean.`); the metric's `METRIC_DESCRIPTIONS` entry is its
detail, panel `metrics`. Without either metric (only an old model's stored record) it is `undefined`.

## Results

`ui/results-pane.ts`, `ResultsPane`: the Train view's results and the content of the `EvaluationView` (see "Run
Evaluation"). The pane holds the insights, the hint and the rows; the viewers and the data under them are built,
refilled and given back by `ui/results-viewers.ts` (`buildShown` with `perClassPanel`, `refillShown`, `markCurves`,
`updateWrongPredictions`,
`releaseShown`, `releaseLinkedFrame`; the panel titles `PANEL_TEXTS`, `resultsTitle`, `resultsPanel`). Its constructor takes the element shown next to the **Performance** title (the Train view's loader host;
the evaluation passes none). Two panels are built once and never rebuilt: the **Performance** panel
(`forge-results-panel forge-results-metrics`: its title with the status, then the metrics host `forge-results-body`
with a hint text or the `MetricsTable` root) and **Insights & Tips** (`forge-insights`: one `li` per insight, its text
only; its tooltip the detail lines, the `items` as a bulleted list and, last, `See <panel title>.` in grey
(`forge-tooltip-minor`) when the insight's `panel` is mounted; repeated texts once; hidden only while empty in the
hint state; with a result never (the `summaryInsight` line of `training/insights.ts` when there is no other, the
metric's description as its detail; nothing when it is `undefined`)). The insights of a viewer panel
(every `InsightPanel` but `metrics`, whose grid marks its own findings) are also its viewer's `description`
(`descriptionVisibilityMode: 'Always'`, the texts joined by a space, escaped as literal Markdown), refreshed with the
list; **Confusions** shows it under the matrix (`descriptionPosition: 'Bottom'`), the others above the viewer. Every panel's title is the built-in tool's widget title, `ui.div([span], 'd4-dialog-header d4-dialog-title
forge-results-title')`; `PANEL_TEXTS` holds each title (`<name>` for the target's or the method's name) with what the
panel shows, the span's tooltip (the span takes the mouse, `d4-dialog-title` lets it through). The root holds the hint
host (the two panels while there is no result) and the rows host (`forge-results-rows`).

- `setDataInsights(insights)`: the check's `SelectionCheck.insights`, set after every check (also when a problem or a
  hint follows), shown first until the next call.
- `show(data: ResultsData)`: `engine`, `table`, `columns` (what the method saw), `target` (the selection's), `prediction`
  (labels or values; the Train view passes `trainPredictionOf(result)`), `score` (only the Run Evaluation of a model
  of the former regression step, `binary-classification`; the Train view passes none),
  `curves` (the result's `Curves` from EDA's report: the positive class's ROC and PR curves for two classes with a
  probability or score), `keptRows`,
  `positiveClass`, `threshold` (Predict probability: the Train view's `binaryClassificationThreshold` of the result, the
  evaluation's of the model, circled as the current **Classification threshold**), `blob`, `summary`. The grid is the `MetricsTable` in place at once (see "Components");
  the result insights are `evaluationInsights` of `summary.evaluation?.values ?? summary.metrics.train` and
  `summary.metrics`, the fallback line `summaryInsight(summary.metrics, summary.evaluation?.values)` (cleared by
  `showHint`); **Insights & Tips** is always one of the mounted panels. A `multiclass` result gets the **Per-class** parallel coordinates plot, built with the other viewers (`buildShown`)
  and refilled with them (see "Keeping the viewers").
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
`binary` for exactly two non-empty classes, `multiclass` otherwise. `metrics` and `insights` are always available; any
other key the result lacks is dropped and an empty row disappears; a row keeps its
template height (`ROW_HEIGHTS` 350 / 350 / 450). Weights in brackets, 1
where none is given:

| kind | row 1 (350 px) | row 2 (350 px) | row 3 (450 px) |
|---|---|---|---|
| `binary` | `predictedVsActual`, `rocCurve`, `prCurve` | `insights`, `confusions` (1.5), `metrics` (1.5) | `features`, `correlation`, `classBalance` |
| `multiclass` | `predictedVsActual`, `confusions`, `perClass` | `insights`, `metrics` (1.5), `classBalance` (1.5) | `features`, `correlation` |
| `regression` | `predictedVsActual`, `residuals`, `distributions` | `insights`, `metrics` (1.5), `statistics` (1.5) | `features`, `correlation` |

Then `visualize:0`, `visualize:1`, ... in rows of at most `VISUALIZE_PER_ROW` (3), `VISUALIZE_ROW_HEIGHT` (350 px),
then `wrongPredictions` alone as the last row (`WRONG_ROW_HEIGHT`, 300 px; in no template), right above the Train
view's **Runs**.

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
| text or yes/no | **ROC Curve** (two classes: with Predict probability, or a method that returns class probabilities; `ResultsData.curves` present) | `scatterPlot` over Forge's own frame of `Threshold`, `FPR (1 − Specificity)`, `Recall (Sensitivity)`, `Mark` (`curveFrame` of the `rocCurve` entry of `CURVES`, from Crux's exact `ml.roc_curve`: one row per distinct probability or score plus the first `(0, 0)` at threshold Infinity, thresholds decreasing), as EDA's probabilistic MPO: a line (`linesOrderColumnName` `Threshold`, `linesWidth` 3) with the markers column `Mark` (`markersColumnName`; `meta.markers` `point` -> `dot`, one pixel whatever the size, `cutoff` -> `circle`; `markerDefaultSize` 8, explicit so a new frame does not resize them; `showCurrentPoint` draws the circle's outline; `legendVisibility: 'Never'`), automatic axis ranges (no `xMin` / `xMax` / `yMin` / `yMax`, so no range sliders; `autoLayout: false`), the dashed formula line `${Recall (Sensitivity)} = ${FPR (1 − Specificity)}` `Non-informative baseline`, the X and Y column selectors, no size or color selector; the frame's current row is the cutoff's point and its only `cutoff` mark (the last row whose threshold is at least the cutoff, as `replayPostprocessing` cuts `score >= cutoff`; set when the frame is built, `markCurves` then moves it: the previous row back to `point`, the new one to `cutoff`, `CurvePanel.marked`), none (-1, every mark `point`) without a cutoff (probabilities only); more classes: no panel; the panel title `ROC Curve (AUC = 0.912)`, the AUC-ROC of **Performance** (`summary.evaluation?.values.auc ?? summary.metrics.train.auc`), the area under the drawn curve |
| | **PR Curve** (the same cases) | the same scatter plot over `Threshold`, `Recall (Sensitivity)`, `Precision`, `Mark` (`curveFrame` of the `prCurve` entry, from `ml.precision_recall_curve`: thresholds increasing, the last point `(recall 0, precision 1)` at Infinity), automatic axis ranges, with the dashed horizontal `Non-informative baseline` `${Precision} = <pr.precision[0]>` (the lowest threshold predicts every row positive: its precision is the positive class's share of the rows); the current row and the circle the first whose threshold is at least the cutoff; the title `PR Curve (AP = 0.812)` from `averagePrecision` the same way |
| | **Per-class** (more than two classes, with `metrics.classes` or an evaluation's `classes`) | `perClassPanel(shown, summary)` in `buildShown`: `DG.Viewer.pcPlot` over Forge's own frame (`Shown.perClass`) of **Class** and a column per `PER_CLASS_AXES` id with a value (`precision`, `sensitivity`, `f1`, `auc`, named by `PER_CLASS_NAMES`: `Precision`, `Recall`, `F1`, `AUC-ROC`; an all-missing axis left out, so SVM gets three) of the shown split (`perClassValues`: the evaluation's when there is one, else the validation's; labels: the union, the evaluation's first); **Class** colored with the target's row colors (`classColors`: `DG.Color.getRowColor` of each class's first row, a user color coding kept, `setCategorical`; a class the rows lack gets `getCategoricalColor`), `colorColumnName` **Class**, the plot's own per-axis ranges (`normalizeEachColumn` default), `showCurrentLine: false` (a click on a line draws nothing extra), `showFilters: false`, `showColorSelector: false`, `legendPosition: 'Top'`, `legendVisibility` `Always` when the class names plus `LEGEND_ITEM_CHARS` (3) per class come to at most `MAX_TOP_LEGEND_CHARS` (40) characters - one line in a third of the pane - else `Never`; in `viewers` like **Features** |
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
Forge's columns, Forge's own frames of **Distributions ...**, **Statistics**, **ROC Curve** and **PR Curve** (each
curve with its viewer and title), the wrong-predictions frame and the viewer panels. A new result whose table and kind
are the same, with or without a score as before, with or without curves as before, the same kept rows
and the same columns and target (the same objects, or copies with the same name, type and values: a skipped row or a
preparation step gives new copies on every training) only refills: the new prediction and score are written into the
old columns in place (`copyInto` of `apply/apply-model.ts`, or `init` for a text column so its categories are rebuilt, then `fireValuesChanged`; only when the type and the length match, otherwise everything
is rebuilt; a score equal to the shown one, as after a threshold change, is not written; curves that are not the shown
object (a retrain: a new report) give the two curve viewers new frames, `viewer.dataFrame = frame` with the axes and the markers column set
again (each new frame has its own `Mark` column, the threshold already circled; `replaceCurves`; the viewers and their panels stay), while a re-cut keeps the result's `curves` object and
the frames), then the marks, the per-class values (each axis column `init` + `fireValuesChanged`, which the plot
redraws on; different classes or axes rebuild everything), the distribution values and the statistics are recomputed
in place, so the viewers redraw with their state (zoom, sizes) kept. A scatter plot recomputes its regression line only when the state of its inputs
changes (its columns, its filter's version, its settings: `getRegressionState`, `d4/.../scatterplot/features/
regression_line.dart`), not when the values of its columns do, so a regression refill fires the frame's filter
(`filter.fireChanged()`, the filter itself unchanged). **Wrong predictions** follows the new prediction: its scattered
column refilled and its filter re-set while a row is wrong, built when one becomes wrong, dropped (and its frame given
back after the mount) when none is. **ROC Curve** and **PR Curve** are on their own frames, so no refill of the
preview frame reaches them: a threshold change or a Predict probability toggle only moves their current rows and the
`cutoff` mark (`markCurves`: `set` of the two rows, then `fireValuesChanged`), and the title, the current row and the marks are set only when they change (on the
platform's ROC curve viewer, which shared the preview frame, every refill of the prediction and the marks redrew the
curve). A hyperparameter, Method, threshold or Predict probability change gives such a result; a change of the
data (Table, Target, Features, other skipped rows or imputed values, another task kind) builds
everything again. A panel that appears or goes (Wrong predictions, the method's own views) changes
only the rows it is in. A mount disposes the viewers it no longer shows with `ui.empty` of their slots and rows (the platform's
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
(`storedSummary`, as Performance reads it) plus `evaluation` (the values, the table's name, the kept rows, the
skipped ones and the per-class `classes`) and the result's `curves`, so the grid has the **Evaluation** column,
**Per-class** shows the evaluation's values, the bullets end with `Evaluation: <table> (<rows> rows)` and
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
6. `replayPostprocessing(raw, options, probabilities)` (see "Preparation"): `binary-classification` (a model of the
   former regression mode): scores at or above `binaryClassificationThreshold` become `positiveClass`, the others
   `negativeClass`, in a column of `targetType`; `compact()` then drops the unused `''` category the new text column
   starts with; `class-probabilities`: the labels of the class probabilities.
7. With skipped rows, the predictions are spread over a full-length column (`scatteredColumn`); skipped rows stay
   empty.
8. Names the column `predictionName(table, target)`: `<target> (predicted)`, or `<target> (predicted 2)`, ... when a
   column of that name exists (ignoring case); tags it `predictive.model` = the model id (`PREDICTION_TAG` in
   `constants.ts`, the tag the platform's earlier tool wrote too) and adds it to the table. A model whose options record `class-probabilities` and whose method
   returned probabilities then adds, after it, the positive class's column (`positiveClass` recorded and among the
   classes) or one per class, each spread over the kept rows as the labels are, named `probabilityName(table, target,
   class)` (`P(<target> = <class>)`, `P(<target> = <class>) 2`, ...; both names come from `freeName`), tagged the same way and described (`meta.
   description`, the grid's header tooltip) `Probability that <target> is <class>, predicted by <model name>`:
   `ApplyResult.probabilityColumns` (empty otherwise). `recordApplication` keeps `column_name` = the label column.

Steps 2-6 are `predictRows(request, progress?)`, which adds nothing to the table: it returns `RowPredictions`, the
columns the method got, its own output (`raw`), the labels (`raw` with the postprocessing replayed, `raw` itself
without a step), `keptRows`, `skippedRows`, the blob read and the class probabilities of the kept rows. `targetSchemaOf(model.target)` reads a stored target
schema with its categories.

**Evaluation** (`training/evaluation.ts`). `evaluateModel(request, progress?)` measures a saved model on a table whose
`target` column holds the actual values (an `EvaluationRequest`: the `LoadedModel`, the table, the target, the column
mapping, the positive class): `checkApplyRequest`, then `predictRows` with the default batch size and the
missing-values mode the model was trained with (`missingValuesOf`, see "Missing values"); the target is the given column, or a copy of its kept rows
(`columnRowCopy`). A model of the former regression step (`binary-classification`) is measured on its scores: the
classification metrics of the labels with the positive class (the request's, else the options') plus `auc` and
`averagePrecision` of the scores (`measure` with no probabilities and the scores as `score`), but no `logLoss` or
`brier` (scores are not probabilities, so they are not passed as probabilities), the scores coming back as `score`; any
other model, one saved with `class-probabilities` too (its labels already cut by the replay), uses `measure(taskOf(target), ..., rows.probabilities,
true)` (`taskOf`: regression for a numerical target; `auc` of the class probabilities when the method returned them:
the positive class's column for two classes, else the macro average). A
row without an actual value is left out by the metrics. The `EvaluationResult` holds the columns, the target and the
prediction of the kept rows, `keptRows`, `skippedRows`, `rowCount` (the kept rows), the metric `values`, the blob,
`positiveClass` (the request's, else the options'), `probabilities` (the class probabilities of the kept rows, when
the method returns them), `classes` (`classValuesOf` of the report) and `curves` (two classes with a probability or
score). Nothing is added to the table and no application is recorded.

No step writes into the table's columns: the engines read their inputs and return new frames, the replay builds new
lists with Forge's own one-hot columns, and the imputer only gets copies of the columns it fills. The pipeline builds
no frame of its own; the imputer's and the engine's frames are given back by the calls that build them. The raw
engine call
(`apply` in `engines/engine-calls.ts`) takes the prediction column out of the engine's result frame, so the user's
table is the column's only parent (`column.dataFrame` is the table). The first column is the prediction; the
`proba_<class>` columns after it are the class probabilities (`EnginePrediction` `{prediction, probabilities?}`, taken
out of the frame too); other columns are dropped. `predictInBatches` concatenates the probabilities of the batches
(the first batch's classes: a model gives the same ones in every batch) and `predictRows` returns them as
`RowPredictions.probabilities`; the application adds the prediction column only.

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
`forgeApply` (the core's `ML | Models | Apply Model...`; the current table), of a model's **Apply...** and the catalog's play icon (the current
table, else the first open one) and of the table pane's run icon (the pane's table, also as `preferredTable`);
without a table it says "Open a table first.". It calls `applyModelDialog({table,
modelId?, switchToTable?, preferredTable?})`, which reads up to 10000 visible model rows once (`MAX_MODELS` of
`constants.ts`), only the
columns applying needs (`APPLY_COLUMNS`; no models: "No Forge models yet. Train and save one with ML | Models |
Train Model... first.") and discovers the engines once. With a `modelId` among them, the dialog opens on `presetTable`:
`preferredTable` (the catalog's **Applicable to**) while it is open, else the current table if the model fits it
(`applicableTables`), else the first open table it fits, else `table`; without one (`forgeApply`) on `table`.
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
`training/train-model.ts`, `trainingEngines`, `problemMessages`, `checkSelection`,
`prepareTraining` and `trainModel`, from `training/k-fold.ts` `DEFAULT_VALIDATION` and `DEFAULT_SEED`, and
`defaultFeatures`, `DEFAULT_STEPS`, `defaultOneHot`, `thresholdProblem`,
`modelNameOf` (`NO_NAME`; the name is written trimmed) and `saveTrainedModel`):

1. The target and the features by name ("The table has no column 'x'." / "The table has no columns 'a', 'b'."); without
   `features`, `defaultFeatures`.
2. The method by name among the complete methods ("No method 'x' is installed. The methods are: ..."); without one,
   the check starts with XGBoost (`PREFERRED_ENGINE` in `train-model.ts`) and takes the method it suggests (`best`), as the view does.
   `auto` (`isAuto`) does the same for the selection's checks; the suggestion is not trained: Find best model
   chooses the method (steps 6-7). `hyperparameters` with `auto` are refused (`The hyperparameters cannot be given
   with method auto: every method starts from its defaults.`).
3. `options`, narrowed key by key (`OPTION_NAMES`; an unknown key is refused with the list) to the view's defaults:
   `missingValues` left unset when not given (the selection's `missingValues` is then undefined and
   `selectedMissingValues` gives the training method's `defaultMissingValues`: `keep` for a method that handles
   missing values, with or without gaps, else `skip`; `impute` with `neighbors` and `distance` checked against
   `Eda:knnImpute`'s inputs, its defaults otherwise; an explicit `keep` for a method that does not handle missing
   values is the check's Missing values problem, `needsRowsWithoutGaps`: when the suggestion replaces the chosen
   method, the selection is checked again against the suggested one), `oneHot` `defaultOneHot`, `skipUniqueCategories` `true`, `predictProbability`
   `false` (kept with `search`; refused on a numerical target: `Predict probability needs a text or boolean target.`),
   `threshold` `DEFAULT_THRESHOLD` (`cutoff` accepted as its old name: the unknown-key check takes it on top of
   `OPTION_NAMES`, which the "Unknown option" list shows, then becomes `threshold` when that is absent, so `threshold`
   wins when both are given; with Predict probability on two classes, 0..1 or the
   view's `Enter a number from 0 to 1.`), `thresholdMetric` `f1` (`DEFAULT_THRESHOLD_METRIC`; one of
   `BEST_THRESHOLD_METRICS`, else `The option 'thresholdMetric' needs one of: f1, balancedAccuracy, mcc, accuracy.`;
   the steps carry it only when neither `threshold` nor `cutoff` is given, so the training optimizes the threshold and
   the run and the model store the metric; a completed run's record takes the result's options), `seed` 42 (`DEFAULT_SEED`), `scheme` `kfold`, `folds` 5, `testFraction` 0.2 (0.05..0.5),
   `stratified` `true`, `shuffle` `true` (`DEFAULT_VALIDATION`), `search` none (`searchChoiceOf`:
   an object of `metric` (a search metric of either task; the task's check is `searchPlanOf`'s), one limit:
   `budgetSeconds` (whole, >= 1) or `maxTrials` (whole, >= 1), never both (`The search options 'budgetSeconds' and
   'maxTrials' cannot be combined: choose one limit.`; neither: 60 s), `strategy` (`equal` or `halving`), `methods`
   (`The search option 'methods' needs a list of method names.`) and `topK` (whole, >= 1), anything else refused with
   that key list; `strategy`, `methods` or `topK` without `auto`: `The search options 'strategy', 'methods' and 'topK'
   need method auto.`). With `auto`: `search` defaults to `{budgetSeconds: 60}`, `missingValues` to `skip`, and `keep` is refused (`The option 'missingValues'
   keep cannot be used with method auto. Use skip or impute.`).
4. The save, checked before any work: a blank name ("Enter a name for the model."), an unknown storage, `reference`
   without `datasetRefOf(table)`.
5. `checkSelection`; every problem goes into one `ForgeError`, target and features first.
6. With `predictProbability`, a final method without `givesProbabilities` is refused (`<Method> gives no class
   probabilities, so Predict probability cannot be used with it.`; not under `auto`). Under `auto`, `autoEngines`
   picks the methods: the pool is `searchableEngines(check.engines)`; each name of `search.methods` must be installed
   (`engineNamed`), in the pool (`<Method> cannot take part in the search: it runs on the server.` / `...: it cannot
   learn from this selection.`), and pass `exclusionOf` (`search.ts`, the rule the Find block's Candidates use, each
   caller with its own words): give probabilities while `predictProbability` is on (the sentence above) and have
   them under a score metric (`isScoreMetric`: `<Method> cannot take part in the search: it gives no class
   probabilities, which <Metric> needs.`); without `methods`, the pool minus those two kinds, and under a quick limit
   (`isQuickLimit`: a trial limit or `QUICK_BUDGET_SECONDS` or less) `quickEngines` of it (`isQuickMethod`);
   `bestCandidates` of the result. The hyperparameters: `defaultHyperparameters`
   overlaid with the given ones, each checked against its `train` input of `visibleHyperparametersOf` (an unknown
   name with the list of settings, so SVM's `probability` is refused; the kind, a whole number for an int,
   `min`/`max`, `choices`), then `withProbabilityInput(engine, ..., predictProbability)`, so the training, `searched`
   and the run record carry the method's probability input.
7. With `showProgress`, the task-bar indicator `Training <Method> model` (cancellable) is the `LoopProgress` of
   `trainModel`; `prepareTraining`, then `trainModel`. With `search`, `searched` first runs the search on the prepared
   request (see "Search": the indicator `Tuning <Method>` or, with candidates, `Finding best model`, then a new
   `Training <Method> model` of the winner) and `trainModel` gets the search rows with the winner's method
   (`request.engine`: also the saved model's and the result's `method`), its settings, `holdOut` and `search`.
8. The `training_run` row (`source: api`): `completed` with the metrics; `cancelled` or `failed` (with the error)
   before the error is rethrown. A refused input or a failed preparation records nothing.
9. With a name, `saveTrainedModel` with the user's feature and target columns, which links the run; nothing is
   uploaded without `storage: copy`. When the link fails after the row is written, the error is thrown with the model
   saved, so the caller gets no `modelId` for it (the model is found by its name).

It returns `TrainModelResult` (`modelId` or null, `runId`, `name`, `method`, `task`, `target`, `features`, `rowCount`,
`skippedRows`, `seed`, `metrics` with, for a classifier, `classes`, and `test` after a search with a hold out, and
`search` after a search) and adds nothing to the table. Errors are rethrown unchanged (`trainAndSave` records
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
   (`visibleHyperparametersOf`: a method's probability input is left out, as from the Hyperparameters column; `default`
   from `defaultHyperparameters`, a key the input does not declare is left out).
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
`preparation` (`preparationOptionsOf` of `options`), `metrics` (`metricsRecordOf` with `classes` and `test`, null without them), `search` (`searchRecordOf`, null without one), `findings`
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
each without one, and `evaluation`), `classes` (the evaluation's per-class values, a classifier with EDA's report)
and `insights`: `evaluationInsights` of the evaluation and the
stored metrics, as the Evaluation view's Insights & Tips without its `summaryInsight` fallback line, as
`{text, detail?, items?}` (no `panel`; empty when there is nothing to say). No column is
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
  for, apply only a model that fits, confirm a delete, where Forge is in the UI). It is self-contained: the runtime
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

`isPredictionColumn(col)` (`Forge:isPredictionColumn`) is true when the column carries a non-empty
`predictive.model` tag (`PREDICTION_TAG`). `predictionTagOf` (`apply/prediction-tag.ts`) reads it: a value with the
earlier tool's model markup (`#{x.<uuid>."<name>"}`, as that tool wrote it) is `legacy` with the name (null when
empty), any other value `forge` with the value as the model id. The panel function `predictedByPanel` (friendly name
**Predicted by**, `PANEL_PREDICTED_BY`; `meta.role: panel`, condition `Forge:isPredictionColumn(col)`) shows, for a
`forge` tag, the card (with its rich tooltip; it reads the card's and the tooltip's columns) of the model whose id the
tag holds, or "The model that predicted this column is no longer available to you."; for a `legacy` tag the text
`Predicted by "<name>", a model of the earlier predictive modeling tool. Forge cannot show it yet.` (without a name:
`Predicted by a model of the earlier predictive modeling tool. ...`), as a text node. The raw tag stays in the
column's tags.

## Table pane

The panel function `tableModelsPanel` (friendly name **Models**, `PANEL_TABLE_MODELS`; `meta.role: panel`, one
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

## Core entry points

The core detects Forge by its presence: when `Forge:applyModel` is registered in the user's browser, the platform's
predictive modeling entry points call Forge's functions by name, and each falls back to the earlier tool when its
function is missing:

| Function | Inputs | The core calls it from |
|---|---|---|
| `forgeTrain` | none (the current table) | `ML \| Models \| Train Model...` |
| `forgeApply` | none (the current table) | `ML \| Models \| Apply Model...` |
| `forgeApp` | none; output `view` | **Browse > Platform > Predictive models**, the URL `/models`, a restored view of the old catalog |
| `applyModel` | `model`, `table`, `columnNamesMap`, `showProgress`, passed by name | `grok.ml.applyModel(name, table, map, showProgress)` for any name that is not a model of the earlier tool |

With Forge installed the core also hides the earlier tool's table pane **Models**, its **Predicted By** row in a
column's **Details**, its **MLFlow** node and its other panes; models of the earlier tool are shown without commands.
Renaming one of these functions or its inputs silently brings the earlier tool back, or breaks `grok.ml.applyModel`;
`ui-tests.ts` ("the core's entry points exist without a menu of their own") pins them. Publishing or removing Forge
takes effect in a session after a page reload.

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
ribbon holds only **Save** (no TRAIN / SAVE switch of one button; a **Train** button right under the inputs appears only
for a method that does not retrain on every change), and the task bar names the method.
Validation uses five folds that cover every row once,
from a seed saved with the model, instead of five overlapping random samples, and it runs from 10 rows instead of
being skipped below 100. MAE and F1 are new; AUC-ROC is measured on the class probabilities of any classifier whose
method returns them (new; the macro one-vs-rest average for more than two classes). Predict probability is not the
built-in tool's regression on a 0/1 target: it keeps the classifier and takes the labels from its class
probabilities, for any number of classes, and the applied model adds the probability columns (new; the built-in
tool's models with `binary-classification` still replay on their scores); the threshold optimized for a metric
(**Best threshold by**, **Set**), log loss and the Brier score are new. Balanced accuracy and MCC are new; average precision, the per-class values (**Per-class**) and the exact ROC
and PR curves come from Crux (new; the built-in ROC curve was the platform's viewer, there was no PR curve). One-hot records the
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
application, and names the prediction column with Forge's own `<target> (predicted)`, tagged `predictive.model` with the model id
(the built-in tool wrote the model's markup there).

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

The table's pane is **Models**, as the built-in tool's (which the core hides when Forge is installed); it is a panel
function's pane, so the platform places it after its own panes instead of right after **Actions**. It lists the models that fit the table,
found in the browser by close column names and kinds, instead of the server's suggested models; each line's run
icon (`model-run`) is **Run prediction** as in the built-in pane, the model's menu is the platform's row menu with Forge's
commands (the built-in menu had **Run** and the model's commands), the tooltip has no **Links** row, and the empty
list says `No model fits this table.` above **Train**. The catalog view's `aiDescription` is Forge's own text.

The view's history is the **Runs** grid and the model cards. **Runs** has one row per training that was not
superseded (a cancelled or failed one too, also a failure before any run was recorded), every row side by side; a
click on a row whose inputs differ from the form's sets them again and trains (a click on a row the form already
matches and shows the completed result of only makes it current; a difference of only Predict probability and the
threshold re-cuts). The model cards
are the built-in tool's compare cards (the current model, "+" to pin, hidden while the current model is pinned, a
trash to remove, a click on a pinned card setting its inputs again, and their colors), extended with the **Performance** and **Settings** tabs switched together, the
settings status as trained, and the
same-data rule: a card is shown only on the target and validation it was trained with (the core cleared every card on
a data change).

The model's menu has no **Add to** (the built-in tool added a model to an open project): a model reaches a Space by
the platform's drag (see Sharing). Forge's Browse node is **Apps > Forge**; with Forge installed, the core opens Forge's
catalog from **Browse > Platform > Predictive models** too and hides its **MLFlow** node (see "Core entry points").
There is no catalog menu item (the built-in tool had none either). The Domains view of `model` can be filtered by
**Tags** (the built-in catalog's tags filter); its **Author** filter has no counterpart.

**Tune** is new (the built-in tool had no tuning or method search): a seeded random search of the method's
declared settings, set in **Settings** (Manual / Auto-tune) with per-setting ranges, choices and pins, a
one limit (time or trials), the **Search** leaderboard and a held out **Test** measurement. **Find best model** (Method
**Auto**, `Forge:trainModel` with `method: 'auto'`) is new too: where the built-in tool guessed a method from the data
(`selectBestModelEngine`, ported as Forge's suggestion), it measures the checked methods on the data, trying each
with its defaults, then an equal share of the limit each (or successive halving down to the best three, tuning
those), and trains the one with the best score less its standard error. The **Validation** group is new
too (the built-in tool always validated on five overlapping random samples): the scheme (stratified K-fold by default for a classifier, K-fold, Hold-out), the folds or the test fraction,
the shuffle and a split seed (42 by default); settings a method ignores are disabled (`meta.mlenabled`).

## Automation surface

The views publish their state in the DOM, so a test waits on what a view says instead of on a timer:

- The Train view's results (`.forge-train-results`, `TrainView.publishState`, run by every `updateControls` and when a
  `ResultsPane.show` settles): `aria-busy="true"` while a check is pending (`check` undefined), a training is queued,
  running or being recorded (`isTraining`) or results are being shown (`showing`, the `ResultsPane.show` calls not
  settled); `aria-invalid="true"` once checked while the selection has a problem (`problem()`, an invalid setting
  included) or its last training ended `failed` or `cancelled` (`isFailed`, cleared by any change but the threshold and a re-cutting Predict probability toggle, and by
  a training of the current selection that completes, a **Train** click after a cancel too). Ready is
  `aria-busy="false"` without `aria-invalid="true"`; a slow method waiting for **Train** is ready.
- On the same element, `data-forge-run-id` (`SHOWN_RUN`) is the id of the `training_run` row whose result the results
  show. `showResults` sets it when a completed training of the current selection is shown, and a re-cut (a threshold
  change or a Predict probability toggle of a method without a probability input) sets the same id again. `showHint` removes it: before the first training, for a selection problem or an invalid
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
- The model cards (`.forge-cards` at the end of the form): each card is a `.forge-card` with `data-forge-run="<number>"` (the
  run's number in the view; a pinned card has its run's number), the current one also `forge-card-current`; a card
  of other data has `style.display` `none`. The current card has `.forge-card-pin` ("+") only, with `style.visibility`
  `hidden` while a pinned card has its inputs; a pinned card has `.forge-card-unpin` (its trash) only; each icon carries its tooltip as
  `aria-label`. A search's winner has `.forge-card-caption` under the header (`Auto-tuned model` / `Auto-tuned
  settings`). The tabs' headers are the platform's `.d4-tab-header[data-source="tab-pane-<name>"]` (`Performance`,
  `Settings`); both panes' tables are `.forge-card-table` (`tr.header` the header row), the Settings pane's root
  `.forge-card-settings` (built on its first show).
- **Search** (`.forge-search`): `SearchPanel.grid` (the leaderboard) and `SearchPanel.chart` (the bar chart, null
  without a leader), each in a `.forge-panel-fill` box of one `.ui-split-h`; the chart's frame is the grid's
  (`chart.dataFrame === grid.dataFrame`), its widget status has a hit area `bar <Method>` per bar. Search and Runs
  sit with the result rows in the results pane's `.ui-split-v`.
- Preparation's probability inputs: the captions **Predict probability** (disabled, `enabled` false, for a method
  without class probabilities), **Classification threshold** (its editor's text has 3 significant digits, its `value`
  the full threshold), its status (`PreparationInputs.thresholdStatus`, a `.forge-tune-status-text` span with
  `aria-live="polite"` in the `.forge-tune-status` row right after the threshold's root: `Best for <metric>` or
  `Customized`; the row has `style.display` `none` while no best exists), then **Best threshold by**
  (`PreparationInputs.thresholdMetric`, its root the next sibling of the status row; a choice of `F1`, `Balanced
  accuracy`, `MCC`, `Accuracy`) with the button `Set` (`PreparationInputs.setThresholdButton`, a child of that input's root
  `.forge-threshold-metric`, after the select and outside `.ui-input-options`; `style.display` `none` while the
  threshold is the chosen metric's best or no best exists); the inputs hide together (`style.display` `none` on their
  roots), the status with them.
- Find best model: **Method**'s first item is `Auto` (`TrainView.isAutoMethod` while it is chosen; `findSettings` is
  the form's `FindSettings`, `findButton` its **Find**). The block is `.forge-find` (`style.display` `none` outside
  Auto, as the method's `.forge-tune` root is under Auto): the `.forge-tune-status` row (`statusText`, `Not searched`),
  the `.forge-tune-toggle` (`toggle`, `aria-expanded`, `false` at first), the `.forge-tune-grid` (`grid`) with the `.forge-tune-label`s
  `Metric`, `Strategy`, `Limit`, `Max time`, `Max trials`, `Candidates` (`candidatesLabel`, the reasons of the
  unlisted methods in its tooltip; `Max time` or `Max trials` hidden by Limit), the select `metric`, the platform radios
  `strategy` and `limit.kind`, the int inputs `limit.minutes` and `limit.trials` (an invalid one has a `validity` and
  its editor the platform's `d4-invalid`) and the platform multi-choice `candidates`
  (`.forge-find-candidates`, one `.ui-input-multi-choice-checks > div` with its checkbox per listed method), then the
  `.forge-tune-bar` with **Find**. A card's Performance has `.forge-card-more` (Show more / Show less) when the run has
  more than three metrics.
- `CollapsibleGroup`: `aria-expanded` on `.forge-group-header` follows `setExpanded`; `aria-invalid` on
  `.forge-group-summary` follows `setSummary` (the Apply dialog's red `N of M matched`).
- Every table of data but the model cards' is the platform's grid, whose widget status publishes the hit areas `cell <row> of <column>`,
  `header <column>` and `row header <row>` (table rows counted from 1) and the readings `rows`, `rows shown`, `columns
  shown`, `current row`, `current column` and `text of cell <row> of <column>`.

The user flows are Gherkin features on `@datagrok-libraries/bdd` in `bdd/` (`bdd/README.md`), run by Playwright
against a stand. Forge's vocabulary is in `bdd/bindings/`: the kinds `group` and `results panel` over the classes
above, `detail` (a caption-value row of a `ui.tableFromMap`, in the context panel or a tooltip) and `model line` (a line of the table
pane **Models**: its run icon, name and menu arrow), a named table set as the current object (the platform's
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

Manifest: `databases/forge/schema.json`, version `0.1.24` (0.1.17 added the `search` json column to `model` and
`training_run`; 0.1.18 describes the current `splitting` shape, 0.1.19 the `space` of a Tune `search`, 0.1.20 the
`postprocessingInfo` ids of `options` (`class-probabilities`), the `metrics` ids `logLoss` and `brier` and the
`search` score "by the chosen validation", and 0.1.21 the `methods` of a Find best model `search` and its `test`
(listed ahead of the Test setting that will write it), 0.1.22 the `eliminatedRound` of those `methods`, 0.1.23 the `strategy`, 0.1.24 the `chosen` of Find best model by
the best score - se, description changes that need only the bump; a description stays under
the 1024 characters of its column). `grok publish` deploys it; a debug publish applies
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
| `options` | `PreparationOptions` | Preparation replay, see below and "Preparation"; Forge writes `{preprocessingInfo, postprocessingInfo, missingValues: {mode, neighbors?, distance?, skippedRows}, skippedColumns?: [names], oneHotCategories?: {<column>: [categories]}}`; with Predict probability `postprocessingInfo` holds `class-probabilities` and, for two classes, `positiveClass`, `negativeClass`, `binaryClassificationThreshold` are set; models of the former regression mode (before 0.1.20) hold `binary-classification` with the same keys and `targetType` |
| `hyperparameters` | `Hyperparameters` | `{<train function input>: value}` |
| `metrics` | `MetricsRecord` | `{train: {<metric id>: number}, validation: {<metric id>: number}, test?: {<metric id>: number} (the final model on the rows held out before a search), positiveClass?, classes?: [{label, rows, train, validation}]}`; ids `mse`, `rmse`, `mae`, `r2`, `accuracy`, `f1`, `sensitivity`, `specificity`, `precision`, `npv`, `auc`, `balancedAccuracy`, `mcc`, `averagePrecision`, `logLoss`, `brier` (`auc`, `averagePrecision`, `logLoss` and `brier` for classifiers whose method returns probabilities; a model of the former regression mode has `auc` and `averagePrecision` of its scores; `balancedAccuracy` and `mcc` for every classifier; `classes` from Crux through `Eda:classificationReport`); class ids `precision`, `sensitivity`, `f1`, `auc`, `averagePrecision` |
| `splitting` | `Splitting` | `{scheme: none \| kfold \| holdout, folds?, trainFraction?, isStratified?, isShuffled?}` (`folds` for kfold, `trainFraction` for holdout; `isStratified` when the split kept each class's share; `isShuffled` false for a split in row order, absent in older models: shuffled; read by `splittingOf`); `none` only in older models |
| `search` (`model`, `training_run`) | `SearchRecord` | `{mode: tune \| best, metric, budgetSeconds?, maxTrials? (the limit; a continued one sums its runs'), strategy?: equal \| halving (best), seed, trials, completed, pruned, failed, elapsedMs, stopped: budget \| cancelled \| exhausted \| trials \| failed, holdOut?: {rows, fraction}, space?: {searched: [{name, type, low?, high?, log?, step?, choices?, simpler?}], pinned: [names]}, methods?: [{method, trials, outcome: complete \| failed \| abandoned \| unstarted, note?, eliminatedRound?}], leaders: [{method, hyperparameters, score, se, foldScores, trials, chosen, note?}]}`; absent without a search (`mode` `best` is Find best model, `methods` its candidates in cost order, `eliminatedRound` the round of successive halving one was dropped after); `space` (Tune) is in the schema's description since 0.1.19, `methods` and `test?: auto \| off \| <fraction>` (not written yet) since 0.1.21, `eliminatedRound` since 0.1.22, strategy and the one limit (budgetSeconds or maxTrials) since 0.1.23, chosen by the best score − se (score + se if minimized) for Find since 0.1.24 |
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
| `apply` | `(df, model) -> dataframe` | Forge takes the first column as the prediction and the `proba_<class>` columns after it as the class probabilities, out of the returned frame |
| `isApplicable` | `(df, predictColumn) -> bool` | |
| `isInteractive` | `(df, predictColumn) -> bool` | Optional. `meta.mlupdate: 'false'` turns off live retraining |
| `visualize` | `(df, targetColumn, predictColumn, model)` | Optional |

Three more `meta` keys on `train`: `mlhasrows: 'true'` (the model file embeds training rows; EDA's SVM and
k-Nearest Neighbors), `mlserver: 'true'` (the method runs on the server: `isServerEngine`, which also keeps a script
in another language than JavaScript and Chemprop out of Find best model) and `mlnativenans: 'true'` (the method takes
rows with missing feature values as they are; EDA's XGBoost, Random Forest, Extra Trees, Hist Gradient Boosting,
CatBoost, Decision Tree and TabM; see "Missing values"). Two optional keys serve the search (see "Search"):
`mlspace` (the JSON search space over the `train` inputs, read by `searchSpaceText` / `searchSpaceOf`; without it the
method cannot be tuned; EDA declares one for every method but Naive Bayes) and `mlparallel: 'true'` (the method takes
concurrent `train` / `apply` calls, read by `isParallel`; EDA's methods, run on its pool of Crux workers). The
optional `mlenabled` says when a `train` input is used (`enabledRules` / `isEnabled` in `engines/engine.ts`, next to
`searchSpaceText`): JSON `{"<input>": <condition> | [<condition>, ...]}`, a condition an object whose clauses all hold,
a list any of them; a clause is `"<other input>": [<values>]` (compared as text), `"<other input>": {"gt": n}` / `{"lt":
n}` or `"task": ["classification" | "regression"]` (the task the method trains; an unknown task holds). `enabledRules`
gives `{}` without the key and undefined for a malformed one (the logic logs nothing; the view warns once and treats it
as no rules). The view disables an input whose rule fails (see "Components"), the search fixes such a dimension to the
default (see "Search"). EDA declares it for SVM (gamma, degree, coef0 by kernel, epsilon for a regression), Linear Regression (iterations while alpha > 0) and XGBoost (seed while subsample or colsample < 1).
The optional `mlprobabilities` says the method returns class probabilities (`proba_<class>` from `apply`): `'true'`,
or the name of a bool `train` input that turns them on (EDA: `'true'` on every classifier but SVM, whose `probability`
input it names; none on Linear and PLS Regression). In `engines/engine.ts`: `givesProbabilities(engine)` (the meta is
`'true'` or names such an input), `hasProbabilities(engine, on)` (whether its trials carry them with Predict
probability `on`: Find best model's metric and probability rules), `probabilityInput(engine)` (that input's name when the meta names a `train` input of
type `DG.TYPE.BOOL`, else undefined), `withProbabilityInput(engine, hyperparameters, on)` (the hyperparameters with
the input set to `on`; unchanged without one) and `visibleHyperparametersOf(engine)` (`hyperparametersOf` without the
input: the form, `TuneSettings`, the cards' Settings tab, the API's `hyperparametersFor`, `methodColumns`,
`Forge:listMethods`, the catalog's Methods tooltip and `searchSpaceOf`'s default, so a space naming the input is refused
by the view and the API alike) and `visibleValuesOf(engine, values)` (the values without the input: the **Runs** Method
tooltip). `hyperparametersOf` and `defaultHyperparameters` keep the input (the contract).

- An engine is complete, and usable, only with `train`, `apply` and `isApplicable`.
- A function engine receives the feature table and the target column. A script engine (`DG.Script`) receives one
  table with a copy of the target appended, and the target name.
- The calls (`isApplicable`, `isInteractive`, `train`, `apply` in `engines/engine-calls.ts`) take the feature columns
  and build the table themselves, a frame of the columns given back by `onFrame` when the call ends.
- `train` returns the blob as a `Uint8Array` or a `DG.FileInfo` (its `data` is used); anything else is a
  `ForgeError`. `apply` must return a dataframe with at least one column; columns named `proba_<class>` after the
  first are read as the class probabilities (AUC-ROC, log loss, Brier, the two-class curves, Predict probability), any
  other column is dropped.
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
- `selectBestEngine(engines, features, target, selected = features)` (`engines/best-engine.ts`) is the built-in tool's
  suggestion, ported as is: a `selected` column with the `Molecule` semantic type gives `Chemprop`; a classification (the target is not numerical)
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
