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
- **Every training attempt is recorded.** Completed, failed and cancelled runs are `training_run` rows; a model row
  appears only on Save. A training superseded by a newer one (`TrainingQueue`) is not an attempt the user finished
  and is not recorded.
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
- **Missing values never reach a method.** Rows with gaps are skipped or the gaps are imputed before the engine is
  called; dates and very large whole numbers (bigint) are not numerical features.
- **Later:** more engines, MLflow models and the migration of existing models.

## Components and data flow

```
package.ts ──> ui/forge-app.ts ───> engines/                     (DG.Func registry)
           │                    ├─> storage/, generated/db.ts    (catalog query, modelsChanged)
           │                    ├─> catalog/                     (Applicable to, Compare)
           │                    └─> ui/model-actions.ts, ui/model-handler.ts, ui/model-comparison.ts,
           │                        ui/apply-model-dialog.ts, ui/data-grid.ts
           ├─> ui/train-view.ts ──> training/ ──> preparation/, engines/, metrics/
           │                    ├─> ui/save-model-dialog.ts, ui/apply-model-dialog.ts, ui/forge-app.ts (the balloon)
           │                    └─> storage/  ──> generated/db.ts  (grok.dapi.domains -> EMS schema forge)
           │                                  ├─> grok.dapi.files  (System:DomainFiles/forge/model)
           │                                  └─> grok.dapi.tables (training data copies)
           ├─> ui/apply-model-dialog.ts ──> apply/, engines/, preparation/, generated/db.ts (model query)
           ├─> apply/ (runApplyModel) ──> preparation/, engines/, storage/ (blob path), generated/db.ts,
           │                              grok.dapi.files (blob read)
           ├─> ui/model-handler.ts ──> ui/model-panes.ts ──> catalog/, generated/db.ts, grok.dapi (users,
           │                                                  entities, permissions), DG.DomainObjectHandler
           │                                                  (audit, sharing), ui/data-grid.ts,
           │                                                  ui/save-model-dialog.ts (writeModelInfo)
           ├─> ui/model-actions.ts ──> catalog/, storage/, ui/model-handler.ts, ui/apply-model-dialog.ts,
           │                           ui/edit-model-dialog.ts
           ├─> ui/model-comparison.ts ──> catalog/, ui/data-grid.ts, ui/model-panes.ts (modelIcon),
           │                              PowerGrid's Forms viewer
           └─> ui/prediction-column-panel.ts ──> ui/model-handler.ts (the card)
```

- `package.ts` holds the annotated functions only: the app `forgeApp` (shown as "Forge", `Browse > Apps`),
  `forgeModels` (`ML | Forge | Models`), `forgeTrain` (`ML | Forge | Train...`), `forgeApply` (`ML | Forge |
  Apply...`), the API function `applyModel`, the autostart `_initForge` (registers the model handler with the model
  commands, and the comparison handler), `isPredictionColumn` and the column panel **Predicted by** (function
  `predictedByPanel`). They delegate to `ForgeApp`, `TrainView`, `openApplyDialog`, `runApplyModel`,
  `ForgeModelHandler.registerOnce`, `registerModelActions`, `ModelComparisonHandler.registerOnce`,
  `isForgePrediction` and `predictedByPane`. `_initForge` reads the registered handlers' names once
  (`DG.ObjectHandler.list()`) and both `registerOnce` calls skip a handler already there, so a second `_initForge` (or
  the test bundle) registers nothing twice; the model commands are registered only with a newly registered handler.
- Every table of data in the UI is `readOnlyGrid` (`ui/data-grid.ts`): `DG.Viewer.grid` without editing, row header,
  current-row indicator and new-row icon, as wide as its host and as tall as its rows (up to 15, then it scrolls: an
  inline height of `(rows + 1) * rowHeight + 2` px), with tooltips from `{headers, cells}` (by column name; the cells
  as functions of the row) through `onCellTooltip` (`gridTooltip`; a cell without its own text shows its value).
  `textColumn` builds the grids' string columns.
- `ForgeApp` lists the methods (section **Methods**, the grid `methodsGrid`: Method, Package, Method type, Roles,
  Hyperparameters; the Method type cell's tooltip says what a function or a script method is, the Roles cell's gives
  each role's meaning, the Hyperparameters cell's each input's description) and the catalog.
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
  simply gone. The panel keeps the object it shows, unless that is a model or a comparison with a model the reload no
  longer has (`isGone`, a delete): then it shows the selection again by the rule above (once: a restored selection
  does it through its change event), and is emptied rather than left on the deleted model; **Edit model** sets it
  again itself (see "Model commands"). Right-clicking a row adds `CATALOG_ACTIONS` (`addModelItems`, with
  **Applicable to** as the preferred table) and, with two or more selected rows, **Compare** to the grid's menu: the
  grid's menu is its own, so the platform's model commands are not in it. The header holds the **Models** title and
  the icons on the title's baseline (`forge-pane-header` in `css/forge.css`), the icons in `--blue-1`, the blue of an
  input's own icons (ui.css), and in the platform's `d4-disabled` grey while disabled; **Applicable to** is on the
  line under it. The grid takes the full width of the view: an inline `width: 100%`, the override the platform's
  400px rule for boxes in a panel (`ui.css`) asks for.
- `modelsChanged` (an rxjs `Subject` in `storage/model-store.ts`) fires after every successful `saveModel`,
  `updateModelInfo` and `deleteModel`, whoever calls them; every open `ForgeApp` subscribes in `this.subs` and reloads
  itself, 300 ms after the last of a burst (tag chips written one by one reload it once).
- `TrainView` (tab **Predictive model** with the core's model icon, `ui.iconSvg('model')` returned by `getIcon()`)
  holds the form, the live validation, the live training and the results. The form is one `ui.form`, so the platform
  aligns every label in it: the group **Data** (Table, Target, Features), the group **Preparation**
  (`PreparationInputs`, below; hidden while none of its inputs applies), the group **Method** (Method and the
  hyperparameters), all expanded at the start, and the **Train** row
  (`ui.buttonsInput`, `trainButton`). The view is a resizable split (`ui.splitH(..., true)`): the form's panel starts 400 px wide
  (`FORM_WIDTH`) and scrolls by itself; **Results** takes the rest, its header line holding the training loader.
  The `trainRow` is right-aligned (`forge-train-row`). The ribbon holds only **Save** (`ui.bigButton`). A Table change rebuilds the form; the groups'
  subscriptions (`formSubs`) are dropped then and in `detach()`, the hyperparameter inputs' (`methodSubs`) on every
  Method change too. The view keeps the latest training in `lastTraining`: `result`, `runId`, the method, the table and the user's feature and target
  columns, the `datasetName` and `fingerprint` taken at training time, and `isSaved`; any change drops it at once.
  **Save** is disabled until a training result exists, while a training runs, and once that result is saved: OK in
  `ui/save-model-dialog.ts` (Name, prefilled, Description, Tags and **Data storage**, see "Saving") calls
  `saveModelAs(info, choice)`, which marks the training saved before it writes, so one training becomes at most one
  model; a failed write enables **Save** again and deletes a copy it uploaded (a failure of that delete is logged, the
  write's error is shown). **Name** is not nullable and its validator refuses a blank name ("Enter
  a name for the model."), so OK is disabled with that tooltip and the dialog's own validation blocks Enter too. The
  same form (`modelInfoDialog`; `saveModelDialog` adds the storage block to it) is **Edit model**. After a save the balloon `Model "<name>"
  saved.` has the links **Apply...** (`openApplyDialog` on the training table, the new model preset) and **Show in the
  catalog** (`ForgeApp.open()`, which focuses an open catalog view instead of opening a second one). The Results grid is a `MetricsTable` of `ui/model-panes.ts`, shared
  with the model's **Performance** pane (it builds one there). The Train view keeps one while results come (a training,
  a cutoff re-cut): `update(summary)` rewrites the values in place and refills the bullet list when the metric rows are
  the same, so the grid neither flickers nor loses its column widths; other rows (task, AUC-ROC, a hint shown in between)
  build a new one. It is a frame of **Metric**, **Train** and **Validation** (format `0.000`) in a
  `readOnlyGrid` with the header texts, the metric descriptions and the full-precision values as tooltips, then a
  bullet list (`ul`): `Rows: <used> used, <skipped> skipped (missing values)` (only when rows were skipped),
  `Validation: <folds>-fold cross-validation on <rows> rows`, `Seed: <seed>` with the seed in a selectable
  `forge-seed` span and the copy icon right after it (**Copy the seed**, balloon `Seed <seed> copied.`), and
  `Positive class: <class>` (two classes only).
- `CollapsibleGroup` (`ui/collapsible-group.ts`, classes `forge-group*` in `css/forge.css`) is the folding block of
  the Train view groups and of the Apply dialog's **Columns** and **More options**: a header with a chevron, the
  caption and an optional summary (red with `forge-group-invalid`), in Diff Studio's style; the chevron sits in a
  fixed-width box (`forge-group-chevron`), so captions line up whether a group is open or folded.
  `expandOnError(inputs)` opens the group when one of its inputs is validated with an error; the state is not
  remembered.
- Table, Target, Features and Method carry short Forge tooltips (`tooltipText`, the caption's tooltip). Hyperparameter
  inputs come from `ui.input.forProperty`, whose caption tooltip is the `train` input's own `description`. For an
  invalid input the platform shows the tooltip text with the validator messages in red below it.
- **Missing values** (`ui/missing-values-inputs.ts`, `MissingValuesInputs`, used by the Train view and the Apply
  dialog) is a radio of **Skip rows** / **Impute** with its options on one line (`forge-inline-radio`), shown only
  when its columns have gaps, with the tooltip "Rows with missing values in: <col>: N, ...". **Impute** is offered
  only when `Eda:knnImpute` exists; choosing it shows the function's own inputs (**Neighbors**, **Distance**, from
  `imputeSettingsOf` with `defaultValuesOf`) right under it. Hiding them puts their defaults back, so a hidden invalid
  value cannot block the form.
- **Preparation** (`ui/preparation-inputs.ts`, `PreparationInputs`, Train view only) owns its `CollapsibleGroup`:
  **Missing values** (the `MissingValuesInputs` above), **One-hot encoding** (shown while a checked feature
  is text or yes/no; checked when it appears if every such feature, the all-unique ones aside while Skip unique
  categories is checked, has at most `ONE_HOT_MAX_CATEGORIES` (20) categories, else unchecked: a deviation from the
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
  catches in them are `applyAndRecord`, which records a failed or cancelled application and rethrows, and
  `TrainingQueue`, whose answer is the outcome with the error in it. `catalog/` (see
  "Catalog logic") builds on `apply/`, `storage/`, `metrics/` and `training/`.

## Training pipeline

`training/train-model.ts`. A `TrainingSelection` is what the user chose: the engine, the checked feature columns
themselves (a list, no frame), the target column, the hyperparameters, a seed, the number of folds (5), the
missing-value handling and the preparation steps (`steps`, a `PreparationSteps`, see "Preparation"; the Train view
takes them from `PreparationInputs.steps()`). `prepareTraining(selection)`
handles the missing values (see "Missing values"), then applies the steps (`prepareFeatures`), and returns the
`TrainingRequest`: `features` and `target`, the selection's columns with the missing values handled (the user's own
columns where nothing changed; the model's feature list, target schema and data summary describe them), and
`prepared` (`PreparedFeatures`), what the method is trained on, with the `options` that record the whole preparation
(`options.missingValues.skippedRows` is the one record of the skipped rows). No frame is built: the engine calls
build their own. After preparation fewer than `2 * folds` rows is a `ForgeError` ("After skipping rows with missing
values, N rows remain; training needs at least 10.").

**Missing values in the view.** **Missing values** (see Components) is the first input of **Preparation** and follows
the checked features. The **Target** tooltip adds "<target>: N rows without a value are skipped." ("1 row without a value is
skipped.") when the target has gaps; it is information, not a problem. After a training with skipped rows,
**Results** lists `Rows: N used, M skipped (missing values)` first in the bullet list under the metrics grid.

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
`isInteractive` says yes on the prepared columns), and `prepared`, those columns and target (none from
`failedCheck`). A method's `isInteractive` that throws makes the whole check throw.

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
live); otherwise **Results** shows `<Method> takes a while on <N> rows, so it does not retrain on every change. Click
Train.` (N, the target's length). **Train** (`train()`) starts one when the selection has no problem and nothing
trains. Every training runs through the view's `TrainingQueue` with a task-bar indicator `Training <Method> model`,
created when the training starts (a request superseded while waiting shows none; cancellable: the indicator's cancel
stops it before its next fit); a loader without text
stands on the **Results** header line, above the previous results (live) or an empty pane (**Train**). The
queued function prepares the data, builds the run record from the request and trains. Before any training, and after a problem, a failure or a cancel, **Results** reads `Choose the
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

**Run records.** `TrainView` writes the `training_run` row right after the queue settles: `completed` with the
metrics, `cancelled` (the yellow `Training was cancelled.`), `failed` with the error message; a `superseded` training
writes nothing, and neither does a failure before the data was prepared (too few rows left). `model_id` is set by
`linkTrainingRun` when the model is saved; deleting the model keeps the run (`setnull`).

The result (`TrainingResult`) carries the blob, the task, the metrics, the target and feature schemas, the
preparation options, the splitting, the seed, the hyperparameters and the row count (the prepared rows the model
learned from); the skipped rows are in `options.missingValues`. The run and model rows record that row count, and the
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
`missingValues` is `{mode, neighbors?, distance?, skippedRows}`. Both are information only: an application never
replays them and chooses its own handling. Imputing before cross-validation lets the held-out rows' values take part
in the imputation, as the built-in tool did; a later version may impute per fold.

The methods' own handling of gaps never matters (SVM refuses them at training; Linear Regression, PLS Regression and
Softmax have none): the Training tests `<Method> trains and applies with missing values under Skip rows and Impute`
train and apply every EDA method on iris with two gaps in both modes and check that no table or column passed to
`train` or `apply` has a missing value.

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

## Saving

`storage/model-store.ts`. `saveModel(fields, blob)` writes the blob to
`System:DomainFiles/forge/model/<uuid>/model.bin` first and then inserts the `model` row with
`blob = file://<that path>`. A failed insert leaves an unreachable file, never a row without its blob.
`deleteModel(id)` soft-deletes the row; it deletes a folder only when `blob` matches
`file://System:DomainFiles/forge/model/<uuid>/<file>`, and then exactly `System:DomainFiles/forge/model/<uuid>`. Any
other value (a file picked in the generic row editor elsewhere, an empty segment, `..`, no blob) leaves the files
alone. It also deletes the uploaded copy of a `copy` model (`deleteTrainingCopy` of `dataset_table_id`, see "Storage
modes"); both deletes are tried whatever the other does, the first failure is thrown after the row is gone, and
`modelsChanged` fires anyway. In the UI, the ribbon
**Save** asks for Name, Description, Tags and **Data storage**, then calls `saveModel` and `linkTrainingRun`; the
catalog's delete and **Delete model** call `deleteModel`. Both storage functions fire `modelsChanged`.

**Data storage in the Save dialog** (`StorageOffer` of `saveModelDialog`): a radio (`forge-inline-radio`, tooltip `What
the model keeps of its training data.`) of `Reference` (only when `datasetRefOf(table)` gives a reference; then the
default), `None` (the default otherwise) and `Copy`, with one `forge-note` line under it that follows the choice:
`A link to <path | query name | a script> is saved; the data stays where it is.` / `Only a summary of the
data is saved.` / `The training columns (<N> rows) will be uploaded to the server.` (N: the table's rows); there are no warnings
about a model file with training rows or a method that runs on the server (a platform setting is planned).
OK passes the choice (`StorageChoice`; `reference` carries the reference the dialog offered, so `datasetRefOf` runs
once per Save) to `saveModelAs(info, choice)`: `copy` uploads the user's feature and target columns first
(`uploadTrainingCopy`) and, when `saveModel` then fails, deletes the copy (`deleteTrainingCopy`; its own failure is
logged) before rethrowing the save's error.

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

`apply/apply-model.ts`. `loadModel(idOrName)` reads, in one query, the `APPLY_COLUMNS` of the model with that id, or
of the one visible model with that name ("No Forge model 'X' is available to you." / "Several models are named 'X'.
Use the model id."). `loadedModelOf(row, engines)` then refuses, in this order, a row without a blob ("has no model
file"), a blob outside Forge's own `forge/model/<uuid>/` layout (`ownBlob` in `storage/model-store.ts`), a row
without a `{name, type}` feature list, and an engine that is not discovered or incomplete ("The method 'X' is not
installed. Install the <package> package."). The result (`LoadedModel`) holds the row, the engine, the features a
table must provide (`requiredFeaturesOf`: the row's feature list without `options.skippedColumns`, which the model
keeps as the record of the user's choice but never needs), the options (`preparationOptionsOf`) and the blob path;
the target is the row's `target_name`.

An `ApplyRequest` is the model, the table, the column mapping, the batch size (`DEFAULT_BATCH_SIZE` 10000) and the
missing-value settings. `applyAndRecord(request, source, progress?)`:

1. Refuses the request when `mappingProblems` reports anything or the table has no rows; a refusal is not recorded.
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
7. With skipped rows, the predictions are spread over a full-length column; skipped rows stay empty.
8. Names the column `predictionName(table, target)`: `<target> (predicted)`, or `<target> (predicted 2)`, ... when a
   column of that name exists (ignoring case); tags it `forge.model` = the model id (`PREDICTION_TAG` in
   `constants.ts`) and adds it to the table.

No step writes into the table's columns: the engines read their inputs and return new frames, the replay builds new
lists with Forge's own one-hot columns, and the imputer only gets copies of the columns it fills. The pipeline builds
no frame of its own; the imputer's and the engine's frames are given back by the calls that build them. The raw
engine call
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
`ML | Forge | Apply...` (the current table) and of a model's **Apply...** and the catalog's play icon (the current
table, else the first open one); without a table it says "Open a table first.". It calls `applyModelDialog({table,
modelId?, switchToTable?, preferredTable?})`, which reads up to 10000 visible model rows once (`MAX_MODELS`), only the
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
  is `missingValuesProblems` for Impute. It stays outside the groups.
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
a dialog: `runApplyModel` resolves the model by id or name, maps the features by exact names, overlays the given
`columnNamesMap` pairs, and refuses the mapping problems before any work, the unmapped features in one sentence first
("The table has no column 'a' the model needs. Map it in columnNamesMap." / "The table has no columns 'a', 'b' the
model needs. Map them in columnNamesMap."). Rows with missing values are skipped; batches of 10000 rows; with
`showProgress` a cancellable task-bar indicator "Predicting <target>". The application is recorded with `source: api`.

## Catalog logic

`catalog/`, one function per file, for the model catalog and the model's panels:

- `applicableTables(row, tables)` (`applicable-tables.ts`): the tables in which `isSuggested` finds a close column for
  every feature the row needs (`requiredFeaturesOf` of its `features` and `options`); none for a row without a
  feature list. It runs in the browser on the open tables; nothing is uploaded. `featureFit(row, tables)` gives the
  names of every feature (`featureSchemasOf`) with them (the card and Details).
- `modelActivity(id)` (`model-activity.ts`): the model's `application` records, newest first, at most
  `MAX_ACTIVITY_ROWS` (100), the full `count` (the list's length below the cap, a count query at it), and `lastRun`,
  the newest application of any status (`when`, `status`), undefined for a model never applied.
- `compareModels(rows)` (`compare-models.ts`): the table "Compare models", one row per model in the given order, with
  the columns Name, Description, Method, Task, Target, Training rows, Created, then `<metric> (train)` and `<metric>
  (validation)` (labels of `METRIC_LABELS`) for every metric of `METRIC_IDS` at least one of the models has; a metric
  a model lacks is empty (stored metrics are read with `metricsRecordOf` of `training/train-model.ts`).
  `COMPARE_COLUMNS` are the `model` columns it reads. Fewer than two models: "Select at least two models to compare."
- `readModelFile(row)` (`model-file.ts`): the bytes of a model's blob in Forge's own storage (`ownBlob`), named
  `<model name>.bin` with `\ / : * ? " < > |` replaced by `_`; for any other blob: "The model 'X' has no model file in
  Forge's storage, so there is nothing to download."
- `updateModelInfo(id, version, {name?, description?, tags?})` (`model-edit.ts`): writes the given fields with
  `forgeDb.models.update`; with a `version`, a row changed since is refused with the platform's
  `DomainVersionConflictError`, without one the write wins. The tags (`string[]`) are written as their column text
  (`tagsText`; an empty list clears the column). It fires `modelsChanged` and returns the row's new version.

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
once (`ForgeModelHandler.registerOnce`), together with the model commands. Its members read only the row's values; a
row made from catalog values (`ForgeModelHandler.rowOf`, through the platform's `rowFrom`, with the system dates as
text) carries the catalog columns only.

- Icon: the model icon. Markup: the icon and the name. Tooltip: the name and **Method**, **Task**, **Target**,
  **Training rows**.
- Card, the built-in tool's (`grok-gallery-grid-item-title`): the name, `Applicable to <open tables it fits>` (only
  when one fits), `Predict <target>`, `by <first five features>...`, `using <method>`, `Created on <YYYY-MM-DD>`. A
  double-click is the platform's Open (the row's own view).
- `renderProperties` reads the full row once and returns `modelAccordion` (`ui/model-panes.ts`), whose title is the
  one the platform gives a row: the model icon, the favorites star (`ui.star`), the name and the row's commands
  (`ui.contextActions`), added before the panes (`Accordion.addTitle` appends); a row the user can no longer read shows
  "The model is no longer available.". It replaces the platform's whole panel of the row (generic
  Details, Shared with, Chats, History, Actions): the platform has no JS surface for its Chats, and its History is
  hosted as a pane. The accordion's panes are built when first opened; Details and Activity share one
  `modelActivity` query:

| Pane | Content |
|---|---|
| **Details** (open) | The description, then **Author**, **Created**, **Updated**, **Table** (`<table> (<rows> rows)`), **Data storage** (`None` / `Reference` / `Copy`), **Data source** (a reference: its `path`, else its query `name`, else `a script`), **Data copy** (a copy: `ui.render` of the uploaded `TableInfo`, `missing` when `grok.dapi.tables.find` gives none), **Last run** (the newest application's time, with ` (failed)` / ` (cancelled)` when it did not complete; `Never` when none), **Applications** (the count), **Features**, **Target**, **Method**, **Task**, **Applicable to** (the open tables it fits; left out when none), then **Tags**: every change is written with `writeModelInfo` (`ui/save-model-dialog.ts`, shared with Edit model; tags only), one write at a time, each against the version the previous one returned, so quick changes are no conflict; a model changed elsewhere meanwhile asks with the platform's conflict dialog to reload (Details is rebuilt from a fresh read) or overwrite |
| **Performance** | `MetricsTable` of the stored metrics (the Train view's grid: Metric / Train / Validation and its bullet list of rows, validation, seed with the copy icon, positive class); "No metrics were recorded for this model." without metrics |
| **Activity** | `N application(s)` and a grid of the newest applications (up to 100): **When**, **Who** (login, the authors read with one `grok.dapi.getEntities`; empty for a deleted user), **Table**, **Rows**, **Prediction column**, **Status**, **Source**, **Duration (ms)**, with a tooltip per header, the time with seconds on **When**, the meaning on **Status** (`Completed: the column was added.`, `Failed: <error>`, `Cancelled: stopped between batches, no column added.`) and on **Source** (`The Apply dialog or the catalog.`, `A script, through Forge:applyModel.`); "Not applied yet." without one |
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
| **Edit model...** | catalog | always | **Edit model** (`ui/edit-model-dialog.ts`): Name (not blank), Description and Tags, prefilled from a fresh read of those columns and `blob`; OK writes them with `writeModelInfo` against the read version, a model changed meanwhile asks to reload (the dialog reopens) or overwrite; then the balloon `Model "<name>" updated.`, and when the current object is that model, it is set again (a row of the values in hand; the panel reads the model itself), so the context panel shows the new name, description and tags |
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
`meta.role: panel`, condition `Forge:isPredictionColumn(col)`) shows the card of the model whose id the tag holds, or
"The model that predicted this column is no longer available to you.". The raw tag stays in the column's tags.

## Sharing

Every model and training run is a platform entity from its insert, and its author holds View, Edit, Delete and Share
on it, so the author sees, applies, edits and deletes his models whatever the table's grants. Sharing goes through
the platform: **Share...** opens the sharing dialog (the row is already promoted), which writes the grant on the
entity's wrapping project. A non-administrator author shares his own model too: eager promotion gives him Share, so
GROK-21041 no longer blocks own models (WO-3, `forge.tester` on the release build). A user a model is shared with
sees it in the catalog, applies it, and the application is recorded (`application` is secured by its model). He
trains, saves and deletes models of his own. Linking a model into a shared Space is another way to share it.

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
of the balloon after Save. Rows with a missing
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
"Last run" from the newest of them. Performance shows the stored metrics; there is no Run Evaluation yet. The panel
has no Chats pane. Two or more selected models show in the panel as forms (the built-in tool compared them only as a
command). Tags are a column of the model, not entity tags. The model file downloads as `<name>.bin` instead
of a zip.

## EMS schema `forge`

Manifest: `databases/forge/schema.json`, version `0.1.12`. `grok publish` deploys it; a debug publish applies
destructive changes without migration scripts, except a `promotion` change on a row table, which is always refused
(`[promotion-change]`). To change it, publish a manifest with only a placeholder table, then the real one, bumping
`version` both times; the first publish drops the tables with their data and their rows' entities and permissions,
but no files (the blob folders of the dropped models stay). `grok api` generates the typed client
`src/generated/db.ts` (`forgeDb.models`, `forgeDb.trainingRuns`, `forgeDb.applications`). Every table also has the
system columns `id`, `version`, `created_on`, `updated_on` and `author_id`; who and when always come from them.

| Table | Security | Why |
|---|---|---|
| `model` | `row`, `promotion: eager`, `defaultRowVisibility: none`, grants `All users: view, edit` | A model is private until shared, like the old models. Eager promotion makes every inserted model a platform entity at once, with View, Edit, Delete and Share for its author's personal group: under visibility `none` a user sees only rows he holds a grant on (administrators see all by design; the 1.28.0 stand does not apply that, see Sharing), so without it a user would not see his own models. Sharing, favorites and comments work on every model. The table Edit grant lets any user insert models; with visibility `none` it reveals no foreign rows |
| `training_run` | `row`, eager promotion, `defaultRowVisibility: none`, same grants | The trainer's experiment history, including failed and cancelled runs, visible to the trainer for the same reason. `model_id` is optional with `onDelete: setnull`, so the history outlives the model |
| `application` | `master`, `delegate: model_id`, `audit: false`, no grants | Secured by the model (see Key decisions). Master tables refuse grants. No audit: high-churn records. `onDelete: cascade` from the model. `status` is `completed`, `failed` or `cancelled`; `skipped_rows` counts the rows left without a prediction because of missing values. Inserts pass `source` and `status` explicitly |

Deletes are soft: a deleted row stays with `is_deleted` and disappears from queries and counts. Current EMS limitation:
a soft delete does not clean up a promoted row's entity and permissions, so every deleted model (and every model a
test run saves and deletes) leaves a live entity until the platform fixes it or the next schema reset drops it; no
client API removes it. The generic Domains view offers **New model...** to anyone with the
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

The `model` filter on `engine_name` is labelled **Method**. The catalog grid reads the same captions with
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

Two more `meta` keys on `train`: `mlhasrows: 'true'` (the model file embeds training rows; EDA's SVM) and
`mlserver: 'true'` (the method runs on the server).

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
  gets XGBoost, SVM and Softmax; Petal.Length gets XGBoost, SVM, Linear Regression and PLS Regression.
- `selectBestEngine(engines, features, target)` (`engines/best-engine.ts`) is the built-in tool's suggestion, ported
  as is: a feature with the `Molecule` semantic type gives `Chemprop`; a classification (the target is not numerical)
  with at least one categorical and one numerical feature gives `XGBoost`; a regression with five or more numerical
  features gives `PLS Regression`, any other regression `Linear Regression`; everything else `XGBoost`. "Numerical"
  and "categorical" are the platform's `isNumerical` / `isCategorical` (dates count as numerical, as in the built-in
  tool); a feature named like the target is not counted. A suggestion not in the list gives the first engine, an
  empty list undefined. The built-in tool matched `<package>: <name>`; Forge's registry keeps one engine per name, so
  it matches the name.
- `isServerEngine(engine)` (`engines/engine.ts`): the training data leaves the browser (data only; no UI uses it yet). True for a script `train` in
  a language other than JavaScript, a `train` with `meta.mlserver: 'true'`, or a name in `SERVER_ENGINES`
  (`['Chemprop']`, whose functions do not say so yet).
- `hasTrainingRows(engine)`: `meta.mlhasrows: 'true'` on `train`.
