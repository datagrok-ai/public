# Forge

Predictive modeling (train, apply, model catalog) as a package; it replaces the built-in `ML | Models` tool at the
switchover. Design, data shapes and what works so far: [ARCHITECTURE.md](ARCHITECTURE.md).

## Entry points

- `src/package.ts` — annotated functions only, each delegating: `forgeApp`, `forgeModels`, `forgeTrain`, `forgeApply`
  (`ML | Forge | ...`), the API `applyModel`, `_initForge` (autostart) and the column panel **Predicted by**.
- `src/ui/`: the screens `train-view.ts`, `apply-model-dialog.ts`, `forge-app.ts` (the catalog); the model wherever it
  shows: `model-handler.ts`, `model-panes.ts`, `model-actions.ts`, `model-comparison.ts`; `data-grid.ts` for tables.
- `src/training/` (`training-queue.ts`: one training at a time), `src/apply/`, `src/preparation/` (`pipeline.ts`: skip
  unique, one-hot, predict probability and their replay; `Eda:knnImpute`), `src/catalog/` (one function per file).
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
- A masked `clone` of a text column keeps every source category (the empty one too): copy rows with `columnRowCopy`,
  never compact a user's column; the k-fold copies keep the full list on purpose.
- The method gets `TrainingRequest.prepared` (Predict probability: a float 0/1 target); `features`/`target` are the
  selection (feature list, fingerprint). Replay: `oneHotCategories`; `skippedColumns` are not needed at apply.
- Awaited engine calls without I/O resume as microtasks: a loop of them must `yieldToEventLoop`, or a cancel click
  and the progress repaint wait for the end (batches: every 50 ms; training: before each fit).
- `ignore-missing` / `impute-missing` are recorded in `options`, never replayed. Prediction columns carry
  `PREDICTION_TAG` (`forge.model`). `ApplicationInsert` needs explicit `source` and `status`.
- Header defaults (`= 20`) are in `Property.initialValue`, strings in quotes (`'RBF'`): use `defaultValuesOf`.
- `deleteModel` deletes only Forge's own: a blob's `forge/model/<uuid>/` folder, a table named `<x> (training data)`.
- Validators read the form's problems, not their argument, and run only on non-empty values. A dialog validates every
  input it was given (hidden and folded too), focuses the last one on `show()` and re-checks Enter only on an input
  change: the Apply dialog reuses one row input per feature name, resets hidden inputs and adds its rows first.
- User-facing text says "method", never "engine". Every schema column declares a `friendlyName`. Never edit generated
  code: `src/package.g.ts`, `src/package-api.ts`, `src/generated/`.
- Schema change: bump `version`, `grok publish local` (destructive diffs applied; a `promotion` change is refused:
  publish a placeholder-only manifest first, which drops the data). `json` columns hold objects only; deletes are soft.
- `model` and `training_run` are eager: every insert, a test's too, creates a platform entity with its author's
  grants, and a deleted row leaves that entity behind (EMS bug; only a schema reset clears them).
- `model.tags` is comma-separated text, not `string_list` (EMS refuses a list in a single-row write): convert only
  through `tagsText` / `tagsOf` in `storage/model-fields.ts`.
- `renderProperties` of the handler replaces the platform's whole panel of a `forge.model` row (History is
  `auditPane`; no Chats); param funcs only add commands. A catalog row's `DomainRow` carries the catalog columns only.
- `grok.dapi.permissions.get` reads an entity's wrapper project, not the row's own grants; `grant` on it is refused in
  tests: share with the platform's dialog (`shareRow`, it fires `onEntityShared`) or `grok s shares add`.
- Tables of data are grids (`readOnlyGrid`), never `ui.table`. A plain `grok.shell.o =` freezes the panel for a
  second and ignores the next plain set: the catalog sets with `setCurrentObject(x, true, true)`.
- The Tags input keeps Enter while its box holds text and adds the chip asynchronously; tests type a chip into
  `input.d4-tags-selector-input` (`typeTag`), since assigning `value` skips that path.
- Never reuse the old tool's identifiers: no `ML | Models` items, no `Models` Browse node, no `predictive.model` tag.
- Test rows are `forge-test-*`, deleted in `finally`; `expect(x, undefined)` checks against `true`; `openTable` records
  no creation script (a reference test opens iris with `openIrisFromFile`). The Train view trains after every change
  (a superseded training writes no run): close it with `closeTrainView`, which waits for `isTraining` before cleanup.

## Commands

```bash
grok build --typecheck   # runs `grok api` first: regenerates src/generated/db.ts and package-api.ts
pnpm run lint            # excludes src/generated/
grok publish local       # debug publish; deploys the schema
grok test --host local   # --category <Area> --skip-build
```
