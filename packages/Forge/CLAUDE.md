# Forge

Predictive modeling (train, apply, catalog), replacing `ML | Models` at the switchover. Design: [ARCHITECTURE.md](ARCHITECTURE.md).

## Entry points

- `src/package.ts` — annotated functions only, each delegating: `forgeApp`, `forgeModels`, `forgeTrain`, `forgeApply`
  (`ML | Forge | ...`), the API `applyModel`, `_initForge` (autostart), the column panel **Predicted by** and the table
  panel **Forge models** (`ui/table-models-panel.ts`).
- `src/ui/`: `train-view.ts` (+ `runs-panel.ts`), `evaluation-view.ts` (both on `results-*.ts`, `metrics-table.ts`),
  `apply-model-dialog.ts`, `forge-app.ts` (catalog), `model-{handler,panes,actions,comparison}.ts`.
- `src/training/` (`training-queue.ts`: one at a time; `insights.ts`; `evaluation.ts`: records nothing), `src/apply/`,
  `src/preparation/` (`pipeline.ts`: steps and their replay; `Eda:knnImpute`), `src/catalog/` (one function per file).
- `databases/forge/schema.json` — the EMS schema; `src/generated/db.ts` is its generated typed client `forgeDb`.

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
- Logic folders (`engines`, `preparation`, `training`, `metrics`, `storage`, `apply`, `catalog`) never import `ui/`,
  log or call `grok.shell`; expected failures throw `ForgeError`. Catches: `applyAndRecord`, `TrainingQueue.run`.
- Pass column lists, never frames: the engine calls, the imputer and the copy upload run on a `sharedFrame` through
  `onFrame`, which gives it back with `releaseFrame` (a frame gets its columns' events). Never write into a user's column.
  `ResultsPane` keeps frames and viewers until `clear()`, a hint or other data; `ui.empty` kills widgets (`close` needs a view).
- A masked `clone` of a text column keeps every source category (the empty one too): copy rows with `columnRowCopy`,
  never compact a user's column; the k-fold copies keep the full list on purpose.
- The method gets `TrainingRequest.prepared` (Predict probability: a float 0/1 target); `features`/`target` are the
  selection (feature list, fingerprint). Replay: `oneHotCategories`; `skippedColumns` are not needed at apply.
- A loop of awaited engine calls must `yieldToEventLoop` (microtasks hold back a cancel click and the repaint).
- `ignore-missing` / `impute-missing` are recorded in `options`, never replayed. Prediction columns carry
  `PREDICTION_TAG` (`forge.model`); `predictRows` adds none. `ApplicationInsert` needs explicit `source`, `status`.
- Header defaults (`= 20`) are in `Property.initialValue`, strings in quotes (`'RBF'`): use `defaultValuesOf`.
- `deleteModel` deletes only Forge's own: a blob's `forge/model/<uuid>/` folder, a table named `<x> (training data)`.
- Validators read the form's problems, not their argument, and run only on non-empty values. A dialog validates every
  input it was given (hidden and folded too), focuses the last one on `show()` and re-checks Enter only on an input
  change: the Apply dialog reuses one row input per feature name, resets hidden inputs and adds its rows first.
- User-facing text says "method", never "engine". Every schema column declares a `friendlyName`. Never edit generated
  code: `src/package.g.ts`, `src/package-api.ts`, `src/generated/`.
- Schema change: bump `version`, `grok publish local` (destructive diffs applied; a `promotion` change is refused:
  publish a placeholder-only manifest first, which drops the data). `json` columns hold objects only; deletes are soft.
- `model` and `training_run` are eager: every insert (tests too) creates an entity, which a delete leaves (EMS bug).
- `model.tags` is comma-separated text, not `string_list` (EMS refuses a list in a single-row write): convert only
  through `tagsText` / `tagsOf` in `storage/model-fields.ts`.
- `renderProperties` of the handler replaces the platform's whole panel of a `forge.model` row (History is
  `auditPane`; no Chats); param funcs only add commands. A catalog row's `DomainRow` carries the catalog columns only.
- `grok.dapi.permissions.get` reads the wrapper project; a `grant` on a row is refused in tests (`shareRow`, `grok s`).
- Tables of data are grids (`readOnlyGrid`), never `ui.table`. A plain `grok.shell.o =` freezes the panel for a
  second and ignores the next plain set: the catalog sets with `setCurrentObject(x, true, true)`.
- Tests type a tag chip into `input.d4-tags-selector-input` (`typeTag`); assigning the Tags `value` skips Enter.
- Never reuse the old tool's identifiers: no `ML | Models` items, `Models` node or pane, `predictive.model` tag.
- A panel function's widget is never detached: a pane holds no subscription. `grok s functions list --filter
  'name = "X"'` matches the friendly name. A JS `ui.bind` row gets only **Properties...**: use `ui.contextActions(row)`.
- A row reaches a Space only by the platform's drag of a model element (`ui.bind` / `ui.makeDraggable` with the row).
  Column inputs keep a removed column: the Train view drops it on `onColumnsRemoved`.
- Test rows are `forge-test-*`, deleted in `finally`; `expect(x, undefined)` checks against `true`; close a Train view
  with `closeTrainView` (waits for `isTraining`); `openTable` records no creation script (use `openIrisFromFile`).

## Commands

```bash
grok build --typecheck   # runs `grok api` first: regenerates src/generated/db.ts and package-api.ts
pnpm run lint            # excludes src/generated/
grok publish local       # debug publish; deploys the schema
grok test --host local   # --category <Area> --skip-build
```
