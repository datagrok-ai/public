# Forge

Predictive modeling (train, apply, model catalog) as a package. It replaces the built-in `ML | Models` tool, which
stays in the platform until the switchover. So far: engine discovery, the EMS schema `forge` and an app that lists the
engines and the model catalog. Design and data shapes: [ARCHITECTURE.md](ARCHITECTURE.md).

## Architecture

- `src/package.ts` — annotated functions only, each delegating to a module: the app `forgeApp` (`//tags: app`,
  shown as "Forge") and `forgeModels` (`//top-menu: ML | Forge | Models`).
- `src/ui/forge-app.ts` — `ForgeApp` view: the engines table and the model catalog grid. `ForgeApp.open()` is the
  only place that catches errors.
- `src/engines/` — `engine.ts` (the `Engine` contract, roles, `isComplete`, `hyperparametersOf`), `engine-registry.ts`
  (`EngineRegistry.discover()`), `engine-calls.ts` (`isApplicable`, `isInteractive` for function and script engines).
- `databases/forge/schema.json` — the EMS schema; `src/generated/db.ts` — its typed client `forgeDb`, generated.
- `src/constants.ts` — identifiers that change at switchover (`APP_NAME`, `MENU_PATH`).
- `src/tests/<area>-tests.ts` — one file per area; the category is the area name (`Engines`, `Storage`, `UI`).

## Glossary

| Term | Code | Meaning |
|---|---|---|
| Model | `forge.model`, `ModelRow` | A trained model: what is needed to apply it, plus the training summary |
| Engine | `Engine` | Functions sharing `meta.mlname`; `meta.mlrole` gives each one's role |
| TrainingRun | `forge.training_run` | One training attempt, saved as a model or not |
| Application | `forge.application` | One application of a model to a table; secured by the model |
| StorageMode | `model.storage_mode` | `none` / `reference` / `copy`: what a model keeps of its training data |

## Conventions and traps

- Comment annotations only, no decorators. The registered name is the TS identifier; `//name:` becomes the friendly
  name.
- Logic folders (`engines/` and later areas) never import `ui/`, never catch, never log, never call `grok.shell`.
  Expected failures throw `ForgeError`; only `ui/` handles them.
- Never edit generated code: `src/package.g.ts`, `src/package-api.ts`, `src/generated/`.
- Schema change: edit `schema.json`, bump `version`, `grok publish local`; a debug publish applies destructive diffs,
  but always refuses `promotion` changes on row tables. The manifest has no root `description`.
- EMS `json` columns hold objects only: `features` is `{columns: [...]}`, never a bare array.
- EMS deletes are soft: deleted rows stay with `is_deleted`. A model becomes a platform entity on its first share.
- EDA is registered on servers as `Eda`, not `EDA`. An engine appears only when its package is installed.
- Never reuse the old tool's identifiers: no `ML | Models` items, no `Models` Browse node, no `predictive.model` tag.
- Test rows are named `forge-test-*` and deleted in `finally`.

## Commands

```bash
grok build --typecheck   # the package build script runs `grok api` first: regenerates src/generated/db.ts
pnpm exec grok api       # codegen only; a bare global `grok api` may be older and reject the manifest
pnpm run lint            # excludes src/generated/
grok publish local       # debug publish; deploys the schema
grok test --host local
```
