# Compute2 Run Comparison — Code Overview

The tool lives in `packages/Compute2/src/components/RunComparison/` and is split into pipeline-stage files with a deliberate separation:

| Stage | File | Dependencies | Role |
|---|---|---|---|
| Shared types | `types.ts` | `datagrok-api` (type-only) | node infos, bindings, targets, `ComparisonEntry`, `isNumericType` |
| Matching | `matching.ts` | none (unit-testable) | name similarity, units compatibility, clustering into targets |
| Selection | `selection.ts` | none (unit-testable) | index-row computation, entry statuses, filters, multi-value compatibility |
| Entry extraction | `entry-extraction.ts` | `datagrok-api`, `compute-utils` (incl. RTD) | FuncCall/DataFrame → `ComparisonEntry` (the only file with server calls) |
| Chart building | `comparison-builders.ts` | `datagrok-api` (no server) | entries + targets → long-format chart DataFrames |
| UI | `RunComparison.tsx` | Vue 3, `webcomponents-vue` | all state, reactivity, rendering |

It's launched from the Model Hub ribbon (`Tools > Compare Runs...` → `Compute2:CompareRuns` in `package.ts:292`), which mounts `RunComparisonApp` — a trivial wrapper around the `RunComparison` component (forwarding `roleOnlyFilter`).

## The overall pipeline

```
FuncCall runs / open tables
   → ComparisonEntry            (entry-extraction: metadata + DataFrames)
   → ComparisonEntryNodes       (plain-data description: scalars + tables)
   → ComparisonTarget[]         (matching: name-clustered candidates)
   → user picks a target        (selection: pickers, statuses, compatibility)
   → chartDf                    (comparison-builders: long-format DataFrame)
   → DG.Viewer                  (bar chart or line chart)
```

## `entry-extraction.ts` — turning runs into entries

**`entryFromFuncCall(call)`** builds a `ComparisonEntry` from a history run. The key trick: if the FuncCall has a `PIPELINE_CONFIG` option (`CONFIG_PATH`, imported from RTD's `funccall-utils`), it's a workflow meta-call, and `collectWorkflowSteps` walks the deserialized pipeline state to find all leaf `funccall` steps; each step's run is loaded via `historyUtils.loadRun` and flattened into the entry. Otherwise the call itself is extracted directly.

`collectWorkflowSteps` is a typed walk built on RTD primitives: `buildTraverseD` (the same DFS traversal RTD's `StateTreeFactory` uses) over `PipelineSerializedState`, with `isFuncCallSerializedState` as the leaf guard — so the serialized format pinned by the LibTests persistence tests is enforced at compile time here. Each leaf gets a stable slash `path` from the `configId` chain (root included) and a ` · `-joined `friendlyPath` from friendly names (root name elided, since it's shared by every step). The walk is wrapped in a try/catch: a legacy or malformed serialized config degrades to an entry with no steps instead of breaking the dialog.

**`extractCallNodes`** does the per-call work: it iterates all input/output params, collecting numeric scalars (int/float/bigint) into `ScalarNodeInfo[]` and DataFrames into `TableNodeInfo[]` + a `Map<path, DG.DataFrame>`. Each node gets a stable `path` (config-id based, the identity) and a `friendlyPath` (step friendly names + captions, for display). It also reads two function-annotation options off DataFrame outputs — `{comparisonIndex: ...}` and `{comparisonSplit: ...}` — which become default picker values later.

**`entryFromDataFrame(df)`** wraps an open workspace table as a "raw" entry with one table and no scalars, so ad-hoc data can join the comparison.

## `comparison-builders.ts` — turning entries into charts

All builders produce a *long/concatenated* DataFrame — there is deliberately no row alignment or delta computation (dropped in commit `2b73683a`); the chart's split-by-Run does the visual comparison:

- **`buildScalarComparison`** — one row per run: `Run | Path | <value>` columns. Run names are forced to string type so numeric-looking names don't break legends.
- **`buildMultiColumnComparison(targets, entries)`** — the workhorse. Concatenates raw rows of every participating run into one frame: index column (typed float/datetime/string via `getIndexKind`), optional split column, a `Run` column, and one float value column per target (nulls where a target has no binding for that run). Duplicate display names get `(2)`, `(3)` suffixes. Returns `null` when fewer than two runs participate.
- **`buildColumnComparison`** — the single-target special case, delegating to the multi version and stripping `valueColumnNames` (whose presence is the downstream marker for "multi-value result").

`getIndexKind` decides how the index behaves: all-datetime → datetime index, any non-numeric → "key index" (string, drives a bar chart instead of a line chart).

## `matching.ts` — matching values across runs

The central problem: different runs (possibly of different models) name the same quantity slightly differently. The solution is name clustering.

- **`nameSimilarity`** — Sørensen–Dice bigram similarity on normalized names (lowercased, separators collapsed). **`nameMatchConfidence`** grades a pair as `exact` / `normalized` / `fuzzy` (≥ 0.7 similarity, `FUZZY_NAME_THRESHOLD`) / no match.
- **`clusterByName`** — greedy clustering: each item joins the best existing cluster that (a) doesn't already contain an item from the same run (at most one binding per entry per cluster), (b) name-matches, and (c) has no hard units mismatch. A higher-confidence match beats an earlier lower-confidence one; ties are broken by similarity score, with the table name as a secondary tiebreaker for columns. A cluster's confidence is the *weakest* member match; a units `warn` (one side has units, the other doesn't) sets `unitsWarning`.
- **`matchScalarTargets(entries)`** / **`matchColumnTargets(entries, indexColumns, splitColumns)`** — feed numeric scalars / numeric table columns into the clusterer and keep clusters covering ≥ 2 runs. These become `ComparisonTarget`s (`kind: 'scalar' | 'column'`) with `bindings` (per-run pointers to the actual value), `confidence`, `coverage/total`, `unitsWarning`. Crucially, **columns only participate if the user has picked an index column for their table**; the index and split columns themselves are excluded. `dedupeTargetKeys` suffixes colliding keys since distinct clusters can share a canonical name.

## `selection.ts` — pickers, statuses, and multi-value compatibility

- **`computeIndexRows(...)`** — builds the rows of the "Index columns" picker UI. One row per table, or — with "Merge same functions" on — one row per group of tables that come from the same function (`nqName` + output name), so picking an index once applies to all runs of that model. Merged rows only offer columns present with the same type in *every* group member; stored selections that no longer match candidates are treated as unset rather than kept stale.
- **`compatibleTargetsFor(anchor, targets, getColumnType)`** — for multi-value mode: returns column targets whose `bindingSignature` (sorted `entryId|tablePath|indexCol|splitCol` tuples) is identical to the anchor's, provided the index is line-chartable (numeric/datetime). Only such targets can share one chart.
- **`getEntryStatuses`** — per-run participation status for the selected target, so the UI can flag excluded runs with a reason (`'no similar data'` vs `'index not set'`).
- Small helpers: `matchesFilter` (substring + fuzzy-token filter for the list search boxes), `isSplitCandidate` (string column ≠ index), `selectionToMap` (validated `Record` → `Map` conversion).

## `RunComparison.tsx` — Vue state and reactivity

### Important state (all in `setup()`)

| State | Purpose |
|---|---|
| `models`, `selectedModel` | Model Hub function list and the model whose History is shown |
| `historySelection` | runs checked in the embedded `History` component |
| `entries` (`shallowRef<ComparisonEntry[]>`) | **the comparison set** — the root of everything downstream |
| `indexSelection`, `splitSelection` | `Record<entryId, Record<tablePath, columnName>>` — user's index/split picks |
| `allowFloatIndex`, `allowDatetimeIndex`, `mergeSameFuncs` | toggles feeding the index pickers |
| `selectedTargetKey` | which target is being charted |
| `multiMode`, `multiKeys` | multi-value mode flag and its selected target keys |
| `indexFilter`, `targetFilter` | list search boxes |
| `chartViewer` | last-created `DG.Viewer`, kept for the workspace-snapshot export |
| `historyHeight`, `sidebarWidth`, `chartHeight` | resizable layout dims |

### The reactive chain of `computed`s

Everything re-derives from `entries` + the selections:

1. `indexColumnsMap` / `splitColumnsMap` — validated Map form of the selections (a pick becomes invisible if its type toggle is off, unless it's the annotated default).
2. `targets` — `matchScalarTargets` + `matchColumnTargets`, sorted by coverage. Changing an index selection immediately re-runs matching.
3. `selectedTarget`, `compatibleTargets` (multi-value candidates for the current anchor), `filteredTargets`.
4. `comparison` — the final result object: calls the right builder based on `target.kind` and `multiMode`, wrapped in `markRaw` (DataFrames must not be made reactive).
5. `indexRows` / `filteredIndexRows` — via `computeIndexRows` for the picker UI.
6. `entryStatuses` — exclusion chips.

### Notable mechanics

- **`addEntry`** dedupes by id and calls `applyAnnotatedDefaults`, which pre-fills index/split pickers from the `comparisonIndex`/`comparisonSplit` annotations — but only where the user hasn't already chosen.
- **`addSelectedRuns`** reloads each checked history run fully (`historyUtils.loadRun`) before extraction, guarded by `isAddingRuns`.
- The **table input** (`ui.input.table`) resets its value in a `setTimeout` because doing it synchronously re-enters the Dart stream controller mid-dispatch (fix from `7ddbfc3b`).
- **Shift+click** on a target row enters multi-value mode grid-style: if the clicked target is signature-compatible with the current selection, both are kept; otherwise the click re-anchors.
- **Chart choice** in `renderComparison`:

  | Situation | Chart |
  |---|---|
  | scalar target | bar chart, split by `Run` |
  | column target, key (string) index | bar chart, split by index, stacked by `Run` |
  | column target, numeric/datetime index | line chart: x = index, y = value column(s), split by `Run` (+ inner split column), `multiAxis: false` |

  Multi-value mode scales the minimum chart height by the number of value columns (`250px` per value).
- **`openInWorkspace`** snapshots: clones `chartDf` into a new table view and re-adds a viewer using `chartViewer.getOptions()`, so user tweaks to the chart carry over.

### UI color coding

| Element | Value | Color |
|---|---|---|
| Source badge | `workflow` | `#7b6fb3` (purple) |
| Source badge | `function` | `#4a90d9` (blue) |
| Source badge | `raw` | `#8a8a8a` (grey) |
| Confidence chip | `exact` | `#3cb173` (green) |
| Confidence chip | `normalized` | `#d9a544` (yellow) |
| Confidence chip | `fuzzy` | `#d97b44` (orange) |
| Exclusion chip | error | `#fbeaea` bg / `#a94442` text (red) |

## Entry points and tests

Besides the Model Hub menu, `History.tsx` has an `allowCompare` prop with a `compareSelected` handler — but that path opens the *old* compute-utils `RunComparisonView`; the new tool explicitly passes `allowCompare={false}` to its embedded History.

Tests live in `src/test/`, mirroring the source split: `run-comparison-matching.ts` (names/units, clustering), `run-comparison-selection.ts` (statuses, filters, index rows), and `run-comparison-builders.ts` (DataFrame building), with shared fixtures in `run-comparison-fixtures.ts`. User docs are at `help/compute/run-comparison.md`.
