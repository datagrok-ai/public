# Retrosynthesis

Cheminformatics plugin that computes retrosynthetic routes for a target molecule using
[AiZynthFinder](https://github.com/MolecularAI/aizynthfinder). The TS side renders an
SVG tree of precursors and reaction steps inside the Datagrok context panel and a demo
view; the heavy lifting runs server-side in a Celery container.

## Architecture

| File | Role |
|---|---|
| [src/package.ts](src/package.ts) | Registers `Chemistry \| Retrosynthesis` panel and `Retrosynthesis Demo`. Debounces molecule changes through an `rxjs` `Subject` (2 s) before calling `updateRetrosynthesisWidget`. |
| [src/utils.ts](src/utils.ts) | `updateRetrosynthesisWidget`: validates the input, calls the backend `run_aizynthfind` function, builds the tabbed tree UI, and wires the "add paths to workspace" / settings icons. |
| [src/tree-creation-utils.ts](src/tree-creation-utils.ts) | Builds tab control + SVG layout from a `Tree[]`. `TAB_ID` attribute on each pane content links a pane back to its precursors object. |
| [src/config-utils.ts](src/config-utils.ts) | User-config storage (`grok.userSettings`) + dialog. Lists config folders under `System:AppData/Retrosynthesis/configs/`, reads `config.yml`, reconciles stored selections against currently available expansion/stock/filter policies. |
| [src/aizynth-api.ts](src/aizynth-api.ts) | TypeScript shape of the Python backend response (`Tree`, `TreeNode`, `Scores`, `ReactionMetadata`, `DataEntry`). |
| [src/const.ts](src/const.ts) | Layout constants, `BASE_PATH`, `CONFIGS_PATH`, and the `DEMO_MOLECULE` sentinel. |
| [src/mock-data.ts](src/mock-data.ts) | `DEMO_DATA` powers the demo when `DEMO_MOLECULE` is the target. `SAMPLE_TREE` is currently unused (kept for reference). Large file — do not load casually. |
| [dockerfiles/aizynthfinder/app.py](dockerfiles/aizynthfinder/app.py) | Celery tasks `run_aizynthfind`, `sync_config`, `check_health`. Function metadata is declared in `#name:`/`#input:` headers above each task. |

## Backend call pattern (cancellation-aware)

`updateRetrosynthesisWidget` runs the backend via `DG.Func.find(...).prepare(...).call(true)`
and tracks in-flight calls in an `activeFuncCalls` map keyed by `FuncCall.id`. Before
issuing a new call it cancels the previous one (`fc.cancel()`); after the await it checks
whether its own entry still exists in the map — if not, the call was superseded and the
result is discarded silently. Preserve this pattern when adding new backend invocations
so a fast-typing user does not see stale results.

`funcId` is initialised to `null` so the catch block can distinguish three states:
**pre-call setup error** (`funcId === null` → always show the message),
**in-call error of the current call** (`funcId` still in `activeFuncCalls` → show),
**superseded call** (`funcId` already deleted → swallow). Preserve this tri-state when
adding new error paths.

## Config resolution

1. `getStoredUserConfig()` reads `grok.userSettings(STORAGE_NAME='retrosynthesis', KEY='config')`.
   Tolerates a legacy plain-string form (older builds stored only the config name) by
   wrapping it into `{configName, expansion: '', filter: '', stock: ''}`.
2. `getConfigFilesFromAppData()` lists subfolders of `System:AppData/Retrosynthesis/configs/`.
3. `updateConfigConsideringConfigYml()` reconciles the stored selection against the
   folder's current `config.yml`. Missing folder/policy → reset that field and surface
   a warning via `grok.shell.warning`.
4. The reconciled config is cached in the module-scoped `currentUserConfig`. Always call
   `await setValidUserConfig()` before reading it from a fresh entry point.

The `DEFAULT_CONFIG_NAME = 'default'` sentinel is sent to the backend as an empty
`config: ''` string — that branch in `app.py` uses the bundled `DEFAULT_CONFIG`.

## Function registration: decorators (not `//name:` comments)

This package uses `@grok.decorators.panel` / `@grok.decorators.func` on a
`PackageFunctions` class instead of the `//name:` comment headers used by most packages.
`grok api` understands both — keep the decorator style when adding new functions to this
package for consistency. Python tasks in `app.py` still use the `#name:` header convention.

## DEMO_MOLECULE sentinel

`DEMO_MOLECULE = 'demo_molecule'` short-circuits the backend call and uses bundled
`DEMO_DATA` from `mock-data.ts`. The `Retrosynthesis Demo` view passes this on first
render so the demo loads instantly without hitting the container.

## Glossary (AiZynthFinder terms)

- **Target** — molecule to retrosynthesize (`Tree.smiles` at the root).
- **Tree / Route / Path** — one synthesis plan; a `Tree` whose nodes alternate
  `type: 'mol'` and `type: 'reaction'`.
- **TreeNode** — `is_chemical` (mol) or `is_reaction`; reactions carry
  `ReactionMetadata` (template, classification, policy info).
- **Precursor** — leaf `mol` node; `in_stock: true` means it's purchasable per the
  selected stock.
- **Scores** — `state score` (combined heuristic, used as the tab label),
  `number of reactions`, `number of pre-cursors`, `number of pre-cursors in stock`,
  `average template occurrence`.
- **Expansion / Stock / Filter policy** — three orthogonal selections from the
  config's `config.yml`. Empty string means "use the config's default" (expansion
  → `select_all()`, stock → all `stock.items`, filter → none).
- **Config folder** — a subfolder of `System:AppData/Retrosynthesis/configs/<name>/`
  containing a `config.yml` plus model/stock files. Synced into the container on
  first use via `_sync_user_config`.

## Backend container

`run_aizynthfind` is the only function the TS side calls in production. The container
syncs the user's config folder lazily via `DatagrokClient.files.sync_dir` guarded by a
`FileLock` (240 s timeout). When changing the Python signature, update both the
`#input:` headers and the `func.prepare({...})` keyword args inside
`updateRetrosynthesisWidget` (utils.ts).
