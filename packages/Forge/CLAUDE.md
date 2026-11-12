# Forge

Predictive modeling (train, apply, catalog) behind the platform's `ML | Models` when installed. Design: [ARCHITECTURE.md](ARCHITECTURE.md).

## Entry points

- `src/package.ts` — annotated functions only, each delegating: `forgeApp`, `forgeModels`, `forgeTrain`, `forgeApply`
  (no menus: the core's `ML | Models` and Browse node call them), the API `applyModel`, `trainModel`, `deleteModel`, `listMethods` (`training/train-api.ts` over
  `engines/method-list.ts`, whose `methodColumns` the app's Methods grid, `methodsFrame`, shares), `listModels`, `getModel`, `updateModel` (`catalog/list-models.ts`,
  `get-model.ts`, `model-edit.ts`), `evaluateModel` (`training/evaluate-api.ts`), `shareModel` (`openShareDialog` in
  `ui/model-actions.ts`), `_initForge` (autostart), the column panel
  **Predicted by** and the table panel **Models** (`ui/table-models-panel.ts`).
- `src/ui/`: `train-view.ts` (+ `runs-panel.ts`, `tune-settings.ts`: Settings and Search options, `find-settings.ts`:
  the Find best model block of Method Auto, built from Tune's static builders (status line, Search options toggle,
  Metric row) plus the platform radio Strategy and multi-choice Candidates, `limit-inputs.ts`: both blocks' Limit
  and Max time / Max trials rows (platform int inputs with `min`/`max`, the platform's validation) and the
  stored-limit migration, `search-panel.ts`: the Search bar chart
  and leaderboard in a splitter, `validation-inputs.ts`: the Validation group, `model-cards.ts`: the cards under the form), `evaluation-view.ts` (on `results-*.ts`), `forge-app.ts` (catalog),
  `apply-model-dialog.ts`, `model-*.ts`; their readiness attributes: ARCHITECTURE.md, "Automation surface".
- `src/training/` (`train-api.ts`: `Forge:trainModel`; `training-queue.ts`: one at a time; `search.ts`: the Tune and
  Find best model (equal split or successive halving) searches, their record and continuation (`canContinue`,
  `continueSearch`);
  `browser-search.ts`: its trials; `evaluation.ts`: records
  nothing, `storedTrainingTable` is the data Run Evaluation and `evaluateModel` default to), `src/apply/`, `src/preparation/` (`pipeline.ts`: steps and their replay; `Eda:knnImpute`), `src/catalog/`.
- `databases/forge/schema.json` — the EMS schema; `src/generated/db.ts` is its generated typed client `forgeDb`.
- `agents/` — Datagrok AI (Grokky) knowledge, published with the package: `package-knowledge.yaml` and the
  self-contained `skills.md` (ARCHITECTURE.md, "Grokky knowledge (agents/)").

## Glossary

| Term | Code | Meaning |
|---|---|---|
| Model | `forge.model`, `ModelRow` | A trained model: what is needed to apply it, plus the training summary |
| Engine (Method in the UI) | `Engine`, `engine_*` columns | Functions sharing `meta.mlname`, roles by `meta.mlrole` |
| TrainingRun | `forge.training_run` | One recorded training: completed, failed or cancelled (a superseded one is not) |
| Application | `forge.application` | One application of a model to a table; secured by the model |
| DatasetFingerprint | `DatasetFingerprint` | Row count, column stats and a hash, kept instead of the training data |
| StorageMode | `model.storage_mode` | `none` / `reference` / `copy`: what a model keeps of its training data |

## Conventions and traps

- Comment annotations only; a `// Words: text` comment right above any function is read as an annotation tag.
- Logic folders (all but `ui/`) never import `ui/`, log or call `grok.shell`; expected failures throw `ForgeError`.
  Catches: `applyAndRecord`, `runTrainModel`, `saveTrainedModel` (record or clean up, rethrow), `runEvaluateModel` (a
  cancel becomes `Evaluation was cancelled.`), `TrainingQueue.run`.
- Pass column lists, never frames: the engine calls, the imputer and the copy upload run on a `sharedFrame` through
  `onFrame`, which gives it back with `releaseFrame` (a frame gets its columns' events). Never write into a user's column.
  `ResultsPane` keeps frames and viewers until `clear()`, a hint or other data; `ui.empty` kills widgets (`close` needs a view).
- A masked `clone` of a text column keeps every source category (the empty one too): copy rows with `columnRowCopy`,
  never compact a user's column; the k-fold copies keep the full list on purpose.
- The method gets `TrainingRequest.prepared` (its target is always the selection's); `features`/`target` are the
  selection (feature list, fingerprint). Replay: `oneHotCategories`; `skippedColumns` are not needed at apply.
- Predict probability is post-processing (`class-probabilities`): every label comes from `replayPostprocessing(prediction,
  options, probabilities)`; a toggle or a threshold change re-cuts (`recutResult`) without retraining. A method's
  probability input (`meta.mlprobabilities` naming a `train` input: SVM's `probability`) is hidden from the form
  (`visibleHyperparametersOf`) and set from the checkbox (`withProbabilityInput`), so toggling it retrains.
  `binary-classification` is the former regression mode: replayed for old models, never recorded.
- A plan has ONE limit, `budgetMs` or `maxTrials` (`searchPlanOf` refuses both). The trial cap counts every trial
  launched. An equal split runs the methods one after another and waits for a method's trials before the next share
  starts. A search is continued only through its `SearchResult.state`
  (in memory: its plan holds the columns); the executor's `resume` keeps each candidate's stream position, so a
  continued Tune or halving search draws the points one longer run would (an equal split splits the new limit
  equally again).
- A loop of awaited engine calls must `yieldToEventLoop` (microtasks hold back a cancel click and the repaint).
- `ignore-missing` / `impute-missing` are recorded in `options`, never replayed. Prediction columns carry `predictive.model`
  (`PREDICTION_TAG`; `predictRows` adds none; an old model's markup is recognized by `apply/prediction-tag.ts`). Inserts of applications and runs need explicit `source`, `status`.
- Header defaults (`= 20`) are in `Property.initialValue`, strings in quotes (`'RBF'`): use `defaultValuesOf`.
- API lists are frames, one-model results objects; `object` outputs describe their keys in the output description
  (the generated API types them `any`). `agents/skills.md` repeats the API annotations and the README's errors:
  a change of an API function's name, inputs, result keys or error sentences updates it too.
- `deleteModel` deletes only Forge's own: a blob's `forge/model/<uuid>/` folder, a table named `<x> (training data)`.
- Validators read the form's problems, not their argument, and run only on non-empty values. A dialog validates every
  input it was given (hidden and folded too), focuses the last one on `show()` and re-checks Enter only on an input
  change: the Apply dialog reuses one row input per feature name, resets hidden inputs and adds its rows first.
- User-facing text says "method", never "engine". Every schema column declares a `friendlyName`. Never edit generated
  code: `src/package.g.ts`, `src/package-api.ts`, `src/generated/`.
- Schema change: bump `version`, `grok publish local` (destructive diffs applied; a `promotion` change is refused:
  publish a placeholder-only manifest first, which drops the data). `json` columns hold objects only; deletes are soft.
- `model` and `training_run` are eager: every insert (tests too) creates an entity, which a delete leaves with its
  grants (EMS bug); since the 2026-10-20 build the entities API no longer returns a newly deleted model's entity.
- `model.tags` is comma-separated text (EMS refuses a list in a single-row write): convert through `tagsText`/`tagsOf`.
- `renderProperties` of the handler replaces the platform's whole panel of a `forge.model` row (History is
  `auditPane`; no Chats); param funcs only add commands. A catalog row's `DomainRow` carries the catalog columns only.
- The context panel follows `modelsChanged` only through `ForgeModelHandler.followChanges` (subscribed once by
  `_initForge`; a gone model empties the panel, a new `version` rebuilds it). `modelsChanged` is per bundle: a test
  that writes from the test bundle and checks the panel starts its own `followChanges` and unsubscribes it.
- `grok.dapi.permissions.get` reads the wrapper project; a `grant` on a row is refused in tests (`shareRow`, `grok s`).
- Tables of data are grids (`readOnlyGrid`), never `ui.table`; the one exception is the model cards (`ui.table`
  by the user's decision of 2026-11-05). A plain `grok.shell.o =` freezes the panel for a
  second and ignores the next plain set: the catalog sets with `setCurrentObject(x, true, true)`.
- Tests type a tag chip into `input.d4-tags-selector-input` (`typeTag`); assigning the Tags `value` skips Enter.
- The core calls `forgeTrain`, `forgeApply`, `forgeApp` and `applyModel` (inputs by name) by name: renaming one or its
  inputs breaks `ML | Models` and `grok.ml.applyModel` (core#1518).
- A panel function's widget is never detached: a pane holds no subscription. `grok s functions list --filter
  'name = "X"'` matches the friendly name. A JS `ui.bind` row gets only **Properties...**: use `ui.contextActions(row)`.
- A row reaches a Space only by the platform's drag of a model element (`ui.bind` / `ui.makeDraggable` with the row).
  Column inputs keep a removed column: the Train view drops it on `onColumnsRemoved`.
- UI tests switch `TuneSettings.storage` to `Forge Tune tests` in the hooks loop at the end of `ui-tests.ts`: a new UI
  category that opens a Train view must be added there. The storage holds one key per method and `Find best` (the
  Find best model block); a test that toggles its methods clears the storage in `finally`.
- Method `Auto` is a view state (`userChoice`), not an engine: `this.engine` stays the suggested method, and
  `Forge:trainModel` takes `method: 'auto'`; the winner's engine must reach both the selection and the request.
  Choosing it keeps `lastTraining` (`requestCheck(true)`): nothing retrains, the results stay.
- The view continues a search only through its kept `SearchState` (`searches`, per mode, Tune per method) while
  `searchInputs` match; a repeated `startSearch` in a test continues the earlier one (the limit adds up). It starts
  nothing while Tune's or Find's problem holds (an invalid limit included), whatever `choice` it is given.
- The bar chart fires `onCategoryClicked` only for a stack segment: the Search chart (no stack) selects its row
  (`onClick: 'Select'`) and the panel restores on a one-row selection of the shared frame. `ui.box` returns a viewer's
  root itself (it is a box), and ui.css `div.ui-div > div.ui-box` is 400 x 300 px: a fill rule must be more specific
  (`div.ui-div.forge-panel-fill > div.ui-box`); a grid keeps its own inline height unless a flex column sizes it.
- Search and Runs are bottom blocks of `ResultsPane` (one `splitV` with the result rows): after showing or hiding one,
  call `resultsPane.layout()`; never `ui.empty` the rows host without taking those boxes out first.
- Test rows are `forge-test-*`, deleted in `finally`; `expect(x, undefined)` checks against `true`; close a Train view
  with `closeTrainView` (waits for `isTraining`); `openTable` records no creation script (use `openIrisFromFile`).

## Commands

```bash
grok build --typecheck   # runs `grok api` first: regenerates src/generated/db.ts and package-api.ts
pnpm run lint            # excludes src/generated/
grok publish local       # debug publish; deploys the schema
grok test --host local   # --category <Area> --skip-build (publishes the last build: `grok publish local` first)
npx grok-bdd run         # bdd/: the user flows (login env: bdd/README.md); bdd/generated/ comes from `grok-bdd compile`
```
