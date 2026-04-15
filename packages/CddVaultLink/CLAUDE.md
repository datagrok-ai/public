# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Datagrok plugin that integrates the platform with the **CDD Vault** registration system (https://www.collaborativedrug.com). Exposes a `CDD Vault` app under `Browse | Apps | Chem` with tabs for Protocols, Collections, Saved Searches, Molecules, and Search, plus a context panel for sketched/selected molecules.

Auth: the CDD Vault API key is read from the package credentials manager under the `apiKey` key (lazy-loaded in `cdd-vault-api.ts` on the first request).

## Architecture

Small plugin, flat `src/` layout — no deep hierarchy to navigate:

- **`package.ts`** — entry point. `PackageFunctions` class registers the `cddVaultApp` (via `@grok.decorators.app`), the app tree browser, context panel, and search function. Wires up the tabbed UI (Protocols / Collections / Saved Searches / Molecules / Search) against a `Vault` tree node.
- **`cdd-vault-api.ts`** — thin wrapper over the CDD Vault REST API. Defines all response types (`Vault`, `Protocol`, `Collection`, `SavedSearch`, `MoleculesQueryResult`, `ApiResponse<T>`, etc.) and query functions (`queryVaults`, `queryMolecules`, `queryMoleculesAsync`, `queryReadoutRowsAsync`, `queryProtocolsAsync`, `queryCollectionsAsync`, `queryBatchesAsync`, `querySavedSearches`, `querySavedSearchById`, `queryProjects`). The `*Async` variants use CDD's async job endpoints and are polled via `getAsyncResults` / `getAsyncResultsAsDf` in `utils.ts`. All HTTP goes through `grok.dapi.fetchProxy` (never raw `fetch`).
- **`utils.ts`** — the bulk of the UI glue. Builds tree nodes (`createVaultNode`, `createNestedCDDNode`), the single tab-loading entry point `createCDDTableView` (see *Tab data loading* below), molecule dataframe construction (`createMoleculesDfFromObjects`, `prepareDataForDf`, `reorderColumns`), context panel (`createCDDContextPanel`), ID→link rendering (`createLinks`, `createLinksFromIds`), URL routing (`handleInitialURL`, `createPath`, `getExportId`), and vault stats.
- **`search-function-editor.ts`** — custom function editor (`SeachEditor`) for similarity / substructure / exact / identity search. Exposed as a docked filters side-panel inside the Molecules tab, wired up by `initializeFilters` in `utils.ts` (not a standalone tab). Two API quirks worth knowing: `getParams()` returns resolved values for the search endpoint (e.g. protocol *id*, not name) and is **not** symmetric — to persist, pre-fill, or share form state, work with the input objects directly. `init()` is async and builds the protocol/run choice lists at the end, so any code that targets those inputs must run after init resolves.
- **`constants.ts`** — tab names (`PROTOCOLS_TAB`, `COLLECTIONS_TAB`, `MOLECULES_TAB`, `SAVED_SEARCHES_TAB`, `SEARCH_TAB`) and `CDDVaultSearchType`.
- **`detectors.js`** — semantic type detectors (loaded separately by platform; do not bundle).
- **`tests/`** + `package-test.ts` — grok-test suite.

### Async query pattern

Large CDD endpoints return an `id` for an async job. Always go through the helpers in `utils.ts`:

- `runAsyncExport(vaultId, () => queryXxxAsync(...), timeoutMinutes, text?)` — returns the raw `ApiResponse`.
- `runAsyncExportAsDf(vaultId, () => queryXxxAsync(...), timeoutMinutes, sdf?)` — returns a `DG.DataFrame`.

Both wrap the three-step *start query → `getExportId` → poll via `getAsyncResults` / `getAsyncResultsAsDf`* sequence. Do not call those three primitives directly from new code unless you need fine-grained control.

### Tab data loading (preview + Load all)

Every CDD endpoint that accepts `?async=true` **must have both a sync and an async wrapper** in `cdd-vault-api.ts` (e.g. `queryMolecules` / `queryMoleculesAsync`, `queryBatches` / `queryBatchesAsync`). The sync wrapper is used for a fast `PREVIEW_ROW_NUM`-row preview; the async wrapper fetches the full result.

All tabs load through a single helper — **`createCDDTableView` in `utils.ts`**:

```ts
createCDDTableView(
  viewName, progressMessage,
  syncFuncName, syncFuncParams,        // sync call with page_size: PREVIEW_ROW_NUM
  asyncFuncName, asyncFuncParams,      // full fetch; pass null to opt out
  vault, treeNode, addFilters?
);
```

Behavior:
1. The sync call runs and its df is shown immediately.
2. If it returned `< PREVIEW_ROW_NUM` rows, the full dataset is already loaded — nothing more happens.
3. Otherwise a ribbon row appears on the `DG.TableView`: *"Showing first N rows"* + **Load all** button.
4. Clicking **Load all** runs the async call in the background with a `DG.TaskBarProgressIndicator`. When it resolves, the view's `dataFrame` is swapped and the ribbon updates to *"Showing all N rows"*.
5. If the user navigates to another tab before the background async resolves, the result is **dropped silently** (guarded by the module-level `openedView` token inside `utils.ts`). No dataframe swap happens on a stale view.
6. Errors: sync failure → empty df + toast, no ribbon. Async failure → preview stays visible, toast, ribbon/button return to initial state so the user can retry.

**Opt-out** (single-call, no preview): pass `null` for `asyncFuncName` and `asyncFuncParams`. Used by Saved Searches, whose flow is an SDF export rather than a paged query.

**When adding a new tab**: add sync+async wrappers in `cdd-vault-api.ts`, add corresponding package-level functions (e.g. `getFoos` / `getFoosAsync`), and call `createCDDTableView` with both names. Do not invent a different loading flow.

### URL / deep-linking

The app is URL-addressable: `cddVaultApp` receives a `path` param (meta.url), and `handleInitialURL` in `utils.ts` walks the `Apps | Chem | CDD Vault` tree to the correct vault/tab/entity. `createPath` builds the reverse. Keep new tabs consistent with this scheme so deep links keep working.

## Glossary — CDD Vault concepts → code

All types live in `src/cdd-vault-api.ts` unless noted. Fields below are what the code actually reads; CDD returns more.

| Concept | Code type(s) | One-line meaning | Key relationships |
|---|---|---|---|
| **Vault** | `Vault` | Top-level tenant / workspace in CDD. Every API call is scoped by `vaultId`. | Contains everything below. |
| **Project** | `Project` | Access-control + organizational grouping inside a vault. | Molecules, Batches, Runs each carry a `projects` list. |
| **Molecule** | `Molecule` | The abstract compound — structure (SMILES/InChI/molfile), calculated props, synonyms, registry #. | Has many `Batch`es; appears in `Collection`s; referenced by `ReadoutRow.molecule`. |
| **Batch** | `Batch` | A specific physical lot of a molecule (salt form, solvent, formula weight, stoichiometry). | Belongs to one `Molecule`; referenced by `ReadoutRow.batch`. |
| **Collection** / **VaultCollection** | `Collection`, `VaultCollection` | User-curated set of molecules. Intended split: `Collection = {id, name}` stub used inside `Molecule.collections`; `VaultCollection` = the full object with the molecule-ID list and project. | `VaultCollection.molecules: number[]` → `Molecule.id`s. |
| **Protocol** | `Protocol` | An assay definition — what is measured, how, and its readout schema. | Has `readout_definitions`, `calculations`, `runs`, `protocol_statistics`. |
| **Protocol.category** | `string` (on `Protocol`) | Vault-configured controlled vocabulary (e.g. `ADME`, `Binding`, `Cellular`). **Values vary per vault** — do not hardcode an enum. Wire format is plain string. | If filtering by category in the UI, derive the list from the vault's own protocol responses. |
| **ReadoutDefinition** | `ReadoutDefinition` | Schema for a single measurement within a Protocol (name, unit, data type, precision, aggregation). | Owned by `Protocol.readout_definitions`. |
| **Calculation** | `Calculation` | Derived-value rule inside a Protocol: takes input readout definitions, emits an output one via `formula`. | Owned by `Protocol.calculations`. |
| **Run** | `Run` | One execution of a Protocol (date, person, place, conditions, source files). | Owned by `Protocol.runs`; referenced by `ReadoutRow.run`. |
| **Readout** | `Readout` (internal) | A single measured value `{value, outlier, modifier?}`. | Appears as values inside `ReadoutRow.readouts`. |
| **ReadoutRow** | `ReadoutRow` | One row of measurements: `(protocol, molecule, batch, run) → {readoutName: Readout}`. The `type` field selects the aggregation level (see below). | Fan-in of Protocol + Molecule + Batch + Run. |
| **SavedSearch** | `SavedSearch` | Server-side *saved query criteria* (not a frozen result set). Fetching one re-runs it; the API returns an async export (`?format=sdf` → `ExportStatus`). | Scoped to a vault. |
| **ExportStatus** | `ExportStatus` | Handle + status of an async CDD job (the `*Async` queries and saved-search fetches return this). | Polled in `utils.ts` via `getAsyncResults` / `getAsyncResultsAsDf`. |

### `ReadoutRow.type` — aggregation levels (fine → coarse)

| Value | Grain |
|---|---|
| `detail_row` | One individual measurement (e.g., a single well / replicate) for one (molecule, batch, run). |
| `batch_run_aggregate_row` | Aggregate across replicates of **one batch in one run** (e.g., mean of 3 replicates on one plate). |
| `batch_protocol_aggregate_row` | Aggregate across **all runs** of a protocol for **one batch** (rolls up multiple days/plates). |
| `molecule_protocol_aggregate_row` | Aggregate across **all batches** of a molecule for **one protocol** — the molecule-level summary. |

### Relationships

```
Vault ─┬─ Project (access grouping)
       ├─ Collection / VaultCollection ── Molecule[]
       ├─ Molecule ── Batch[]                           ← physical lot of a molecule
       ├─ Protocol ─┬─ ReadoutDefinition[]              ← what is measured
       │            ├─ Calculation[]                    ← derived readouts
       │            └─ Run[]                            ← one execution
       ├─ ReadoutRow  ── (protocol, molecule, batch, run) → {name: Readout}
       └─ SavedSearch (saved query criteria; fetching → ExportStatus → SDF)
```

## Conventions specific to this package

- Do not edit `package.g.ts` or `package-api.ts` — regenerated by `grok api`.
- Functions are registered via `@grok.decorators.*` on static methods of `PackageFunctions` (this package uses the decorator form, not `//name:` comment metadata).
- Never call CDD endpoints with raw `fetch` — always go through `grok.dapi.fetchProxy` (handled centrally in `cdd-vault-api.ts`).
- Plugin changelog: add a `## v.next` bullet to `CHANGELOG.md` for any sizable change (see repo-level rule).
- **Per-user persistence** — use `grok.userSettings.add(key, subKey, value)` / `grok.userSettings.getValue(key, subKey)` for state that should survive reloads (e.g. remembered form inputs). Values are strings — `JSON.stringify` / `JSON.parse` for objects. Convention in this package: top-level key `CDDVaultLink.<purpose>`, subkey scoped to whatever makes the stored state unambiguous (typically `vaultId`). Do not use `localStorage` or files for this.
