# Forge

Predictive modeling (train, apply, model catalog) as a package; it replaces the built-in `ML | Models` tool at the
switchover. Design, data shapes and what works so far: [ARCHITECTURE.md](ARCHITECTURE.md).

## Architecture

- `src/package.ts` — annotated functions only, each delegating to `ui/`: `forgeApp` (app "Forge"), `forgeModels`
  (`ML | Forge | Models`), `forgeTrain` (`ML | Forge | Train...`).
- `src/ui/forge-app.ts` (`ForgeApp`: methods, catalog, refresh and delete icons),
  `src/ui/train-view.ts` (`TrainView`, tab "Predictive model": form, live validation, results, ribbon Save),
  `src/ui/save-model-dialog.ts`, `src/ui/report-error.ts` (the error boundary); `css/forge.css`.
- `src/engines/` — the `Engine` contract, `EngineRegistry.discover()`, `engine-calls.ts` (`isApplicable`, `train`,
  `apply` for function and script engines), `defaultHyperparameters`.
- `src/training/train-model.ts` — `trainModel` (k-fold + final fit), `trainingProblems` / `checkTrainable`, and the
  json-shape types (`TargetSchema`, `MetricsRecord`, ...); `default-features.ts`, `k-fold.ts`.
- `src/metrics/metrics.ts` — metrics ported from the core's `predictive_modeling_metrics.dart` and
  `confusion_matrix_core.dart`, plus MAE and F1.
- `src/storage/` — `model-store.ts` (`saveModel`, `deleteModel`, `modelsChanged`), `model-fields.ts` (row payloads;
  the json types live in `training/train-model.ts`), `dataset-fingerprint.ts` (`DatasetFingerprint`),
  `training-run-store.ts`.
- `databases/forge/schema.json` — the EMS schema; `src/generated/db.ts` — its typed client `forgeDb`, generated.
- `src/tests/<area>-tests.ts` — categories `Engines`, `Metrics`, `Training`, `Storage`, `UI`.

## Glossary

| Term | Code | Meaning |
|---|---|---|
| Model | `forge.model`, `ModelRow` | A trained model: what is needed to apply it, plus the training summary |
| Engine (shown as Method in the UI) | `Engine`, `engine_*` columns | Functions sharing `meta.mlname`; `meta.mlrole` gives each one's role |
| Hyperparameters | `Hyperparameters` | The `train` inputs except `df` and `predictColumn` |
| TrainingRun | `forge.training_run` | One training attempt (completed, failed, cancelled), saved as a model or not |
| TrainingProblems | `TrainingProblems` | Why a selection cannot be trained, per input (`target`, `features`) |
| DatasetFingerprint | `DatasetFingerprint` | Row count, column stats and a hash of the training data, kept instead of it |
| Application | `forge.application` | One application of a model to a table; secured by the model |
| StorageMode | `model.storage_mode` | `none` / `reference` / `copy`: what a model keeps of its training data |

## Conventions and traps

- Comment annotations only; a `// Word: text` comment right above an exported function is read as an annotation tag.
- Logic folders (`engines`, `training`, `metrics`, `storage`) never import `ui/`, catch, log or call `grok.shell`;
  expected failures throw `ForgeError`. Run records are written by `TrainView` (the UI boundary), not `training/`.
- `apply` in `engine-calls.ts` is the raw engine call, not the application pipeline (a later `src/apply/`).
- Header defaults (`= 20`) are in `Property.initialValue`, not `defaultValue`: use `defaultHyperparameters(engine)`.
- Blobs are written by `saveModel` and removed only by `deleteModel`, and only when `blob` matches
  `file://System:DomainFiles/forge/model/<uuid>/<file>` (then the whole `<uuid>/` folder); otherwise the row only.
- `TrainView` validators read its `problems`, not their argument, and run only on non-empty values.
- Binary metrics use the target's first category as the positive class, stored as `metrics.positiveClass`.
- User-facing text says "method", never "engine". Every `model` / `training_run` column declares a `friendlyName`:
  without one the Domains grid shows the raw name ("Engine_kind"), while the client registry derives "Engine kind".
- `saveModel` and `deleteModel` emit `modelsChanged` (`storage/model-store.ts`); every open catalog reloads.
- Never edit generated code: `src/package.g.ts`, `src/package-api.ts`, `src/generated/`.
- Schema change: edit `schema.json`, bump `version`, `grok publish local` (a debug publish applies destructive
  diffs, refuses `promotion` changes). `json` columns hold objects only (`{columns: [...]}`); deletes are soft.
- Never reuse the old tool's identifiers: no `ML | Models` items, no `Models` Browse node, no `predictive.model` tag.
- Test rows are named `forge-test-*` and deleted in `finally`. `expect(x, undefined)` checks against `true`.

## Commands

```bash
grok build --typecheck   # runs `grok api` first: regenerates src/generated/db.ts and package-api.ts
pnpm run lint            # excludes src/generated/
grok publish local       # debug publish; deploys the schema
grok test --host local   # --category <Area> --skip-build
```
