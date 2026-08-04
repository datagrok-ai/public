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
   → chartDf                    (comparison-builders: long- or wide-format DataFrame)
   → DG.Viewer                  (bar/line/scatter chart; radar or PC plot for multi-scalar)
```

## `entry-extraction.ts` — turning runs into entries

**`entryFromFuncCall(call)`** builds a `ComparisonEntry` from a history run. The key trick: if the FuncCall has a `PIPELINE_CONFIG` option (`CONFIG_PATH`, imported from RTD's `funccall-utils`), it's a workflow meta-call, and `collectWorkflowSteps` walks the deserialized pipeline state to find all leaf `funccall` steps; each step's run is loaded via `historyUtils.loadRun` and flattened into the entry. Otherwise the call itself is extracted directly.

`collectWorkflowSteps` is a typed walk built on RTD primitives: `buildTraverseD` (the same DFS traversal RTD's `StateTreeFactory` uses) over `PipelineSerializedState`, with `isFuncCallSerializedState` as the leaf guard — so the serialized format pinned by the LibTests persistence tests is enforced at compile time here. Each leaf gets a stable slash `path` from the `configId` chain (root included) and a ` · `-joined `friendlyPath` from friendly names (root name elided, since it's shared by every step). The walk is wrapped in a try/catch: a legacy or malformed serialized config degrades to an entry with no steps instead of breaking the dialog.

**`extractCallNodes`** does the per-call work: it iterates all input/output params, collecting numeric scalars (int/float/bigint) into `ScalarNodeInfo[]` and DataFrames into `TableNodeInfo[]` + a `Map<path, DG.DataFrame>`. Each node gets a stable `path` (config-id based, the identity) and a `friendlyPath` (step friendly names + captions, for display). It also reads the `{comparison: {...}}` function-annotation option off DataFrame outputs (index/split/mode/units, with `{comparisonIndex}`/`{comparisonSplit}` as legacy aliases) — these become default picker values later.

**`entryFromDataFrame(df)`** wraps an open workspace table as a "raw" entry with one table and no scalars, so ad-hoc data can join the comparison.

## `comparison-builders.ts` — turning entries into charts

All builders produce a *long/concatenated* DataFrame — there is deliberately no row alignment or delta computation (dropped in commit `2b73683a`); the chart's split-by-Run does the visual comparison:

- **`buildScalarComparison`** — one row per run: `Run | Path | <value>` columns. Run names are forced to string type so numeric-looking names don't break legends.
- **`buildMultiScalarComparison(targets, entries)`** — the one *wide*-format exception: one row per run, `Run` plus one float column per selected scalar target (null where a run lacks that scalar), because radar and PC plot consume axes as columns. A run participates if it has a binding in at least one selected target; returns `null` when fewer than two runs participate. Shares the display-name dedupe (`(2)`, `(3)` suffixes) with the column builder.
- **`buildMultiColumnComparison(targets, entries, axisModes?)`** — the workhorse. Consumes each target's *enabled* bindings (matching guarantees at most one per run). Concatenates raw rows of every participating run into one frame: index column (typed via `getIndexAxis`), optional split column, a `Run` column, and one float value column per target (nulls where a target has no binding for that run). Duplicate display names get `(2)`, `(3)` suffixes. Returns `null` when fewer than two runs participate. In multi-value scatter mode the wide value columns are replaced by a *melted* pair — a `Data` column holding the (deduped) target name and a `Value` column with the number, rows repeated per target — because a scatterplot charts a single y column; the melted names are deduped against the index/split/`Run` columns and returned as `melted: {seriesColumnName, valueColumnName}`.
- **`buildColumnComparison`** — the single-target special case, delegating to the multi version and stripping `valueColumnNames` (whose presence is the downstream marker for "multi-value result").

`getIndexAxis` decides how the index behaves. Legacy rules (no `AxisConfigMap` passed, or the require-all gate fails): all-datetime → datetime index, any non-numeric → "key index" (string, drives a bar chart instead of a line chart), else float. With the `timeseries` mode configured on every participating table (and every index actually numeric/datetime), the axis is `elapsed`: each run's index converts to ms (datetime via epoch ms, numeric scaled by its declared units), always shifts by its own table's min non-null value (auto-alignment — first point is 0), and lands in one display unit as a 64-bit float column named `` `<index> (<unit>)` `` (`Column.fromFloat64Array` — `fromList(FLOAT)` can store 32-bit, which would quantize large values; label deduped against value columns). The result carries `timeSeriesUnit` on the elapsed path. The `points` mode does not touch the axis at all — the same require-all gate just sets `isScatter` on the result (a key index disables it), which flips the chart from a line to a scatterplot. The multi-value require-all gates are evaluated on `targets[0]`'s participating bindings — secondary targets' cross-table picks are padded rows and can't break the mode independently.

## `matching.ts` — matching values across runs

The central problem: different runs (possibly of different models) name the same quantity slightly differently. The solution is name clustering.

Two design rules keep the model simple:

1. **A target maps the same quantity across runs — at most one item per run.** Matching is a cross-source mapping with a preview, nothing more.
2. **Several series within one run are the split column's job**, the native pattern for model results — never several mapped columns of one run. (Suffix-style wide data, e.g. `temp_A`/`temp_B`, should be reshaped into a split column, not multi-mapped.)

- **`nameSimilarity`** — Sørensen–Dice bigram similarity on normalized names (lowercased, separators collapsed). **`nameMatchConfidence`** grades a pair as `exact` / `normalized` / `fuzzy` (≥ 0.7 similarity, `FUZZY_NAME_THRESHOLD`) / no match.
- **`tablesCompatible`** — the table-level gate for columns: two tables are comparable only if their index columns name-match and their split columns are either both absent or name-match (a split table charts per-category series, so pairing it with an unsplit one would be misleading). Column-vs-column compatibility stays per-column (name + units) — there is deliberately no whole-column-set matching.
- **`clusterByName`** — greedy clustering: each item joins the best existing cluster that (a) doesn't already contain an item from the same run (at most one binding per entry per cluster), (b) name-matches, (c) has no hard units mismatch, and (d) passes `tablesCompatible` against the cluster seed when both carry a table key (scalars don't). A higher-confidence match beats an earlier lower-confidence one; ties are broken by similarity score, with the table name as a secondary tiebreaker for columns. A cluster's confidence is the *weakest* member match; a units `warn` (one side has units, the other doesn't) sets `unitsWarning`.
- **`matchScalarTargets(entries)`** — feeds numeric scalars into the clusterer and keeps clusters covering ≥ 2 runs; a scalar target's `bindings` are exactly the greedy cluster members (no candidate model).
- **`matchColumnTargets(entries, indexColumns, splitColumns, overrides)`** — the candidate model for columns, implementing rule 1 above. The greedy pass seeds clusters as before; every other item compatible with a cluster's seed (name + units + table gate) is listed as a `ColumnCandidate` (`auto` = greedy member; `confidence`/`unitsWarn` vs the seed; `enabled`). Enablement is radio-style — at most one enabled candidate per run; an explicit user pick beats a default one, ties resolve to the first in candidate order. Defaults: non-raw candidates are enabled iff `auto`; raw (standalone) candidates are enabled in every cluster they match non-fuzzily (exact/normalized name) — fuzzy matches stay listed but unchecked. Cluster survival requires ≥ 2 distinct runs among *default-enabled* candidates, so user toggles never remove a target and a raw attachment can resurrect a single-run cluster. `overrides` (targetKey → `candidateId` → enabled) apply after `dedupeTargetKeys`; the derived fields — `bindings`, `coverage`, `confidence`, `unitsWarning` — reflect the enabled subset only. `candidateId` (`entryId|tablePath|columnName`, index/split excluded) keeps toggles stable across picker changes. Crucially, **columns only participate if the user has picked an index column for their table**; the index and split columns themselves are excluded.

## `selection.ts` — pickers, statuses, and multi-value compatibility

- **`computeIndexRows(...)`** — builds the rows of the "Index columns" picker UI. One row per table, or — with "Merge same functions" on — one row per group of tables that come from the same function (`nqName` + output name), so picking an index once applies to all runs of that model. Merged rows only offer columns present with the same type in *every* group member; stored selections that no longer match candidates are treated as unset rather than kept stale. Takes the raw `AxisModeSelection` and sets `IndexRow.axis` (`{mode, units?}`) when the current index is numeric/datetime — `units` for numeric only, prefilled from column metadata; merged rows agree on a non-default mode like `current` does (any disagreement displays as `series`).
- **`multiValueOverlap(anchor, other)`** / **`compatibleTargetsFor(anchor, targets, getColumnType)`** — multi-value mode is split into a *suggestion* predicate and builder-enforced validity. Per anchor run, the other target's pick is `aligned` (same table — rows can be shared), `missing` (no pick), or `conflicting` (picked from another table); scalar bindings carry no table, so scalars are aligned or missing, never conflicting. `compatibleTargetsFor` never mixes kinds: for a scalar anchor every scalar target is compatible (all scalars are numeric by construction); for a column anchor it *suggests* column targets with no conflicts and ≥ 1 aligned run, provided the anchor's index is line-chartable (numeric/datetime) everywhere. It never decides validity of an already selected combination: the builder pads missing/conflicting runs with nulls, so editing picks can never eject a selected target. Index/split agreement is automatic for aligned picks — the pickers are per (run, table). Note the layering: cross-source reconciliation (which run/step/raw table holds the same quantity) is already settled by the clusterer and recorded in each target's bindings — so a single target freely mixes a workflow step's table with a plain CSV; the overlap check only asks whether two *targets* can read their columns from the same physical rows.
- **`getEntryStatuses`** — per-run participation status for the selected target, so the UI can flag excluded runs with a reason (`'no similar data'`, `'index not set'`, or `'disabled'` — the run has candidates but the user toggled them all off). With an `AxisConfigMap` passed, matched runs on a *partially* configured target additionally get `warning: 'relative timeseries not set'` or `'independent points not set'` (amber chip — the run still charts, just as a plain series).
- **Axis mode helpers** — `parseTimeUnit` (units-metadata prefill; bare `m` deliberately unrecognized), `resolveAxisModes` (stored picks → non-`series`-only `AxisConfigMap`, gated on a numeric/datetime index; units for `timeseries` on a numeric index only), `targetAxisMode` (require-all gate per mode: `full`/`partial`/`none` over a target's bindings), `resolveDisplayUnit` (first numeric `timeseries` binding's units, datetime-only → `'auto'`), `pickAutoUnit` (largest of days/h/min/s spanning ≥ 2 units, else ms).
- **`isTargetEqualAcrossRuns(target, getColumnValues)`** — backs the "Hide equal values" toggle: a target with ≥ 2 bound runs is equal when all scalar values match, or when every run's value, index, and split column contents are element-wise the same (NaN-safe). Numbers compare via the PEP 485 `isclose` formula with `EQUALITY_REL_TOL = 1e-3` (0.1% relative, zero absolute tolerance — near-zero is never conflated with zero). Unfetchable column data counts as differing, so nothing is hidden on missing data.
- Small helpers: `matchesFilter` (substring + fuzzy-token filter for the list search boxes), `isSplitCandidate` (string column ≠ index), `selectionToMap` (validated `Record` → `Map` conversion).

## `RunComparison.tsx` — Vue state and reactivity

### Important state (all in `setup()`)

| State | Purpose |
|---|---|
| `models`, `selectedModel` | Model Hub function list and the model whose History is shown |
| `historySelection` | runs checked in the embedded `History` component |
| `entries` (`shallowRef<ComparisonEntry[]>`) | **the comparison set** — the root of everything downstream |
| `indexSelection`, `splitSelection` | `Record<entryId, Record<tablePath, columnName>>` — user's index/split picks |
| `axisModeSelection` | per-table `{mode, units}` row-semantics picks — `series` (default), `timeseries` (relative), or `points` (independent) — stale picks kept but ignored; `resolveAxisModes` derives the non-default-only `AxisConfigMap` consumed by builders and statuses |
| `mergeSameFuncs` | merge toggle feeding the index pickers (numeric/datetime/string indexes are always offered — the old float/datetime gating toggles are gone) |
| `hideEqual` | "Hide equal values" toggle (on by default): targets whose data is identical across all bound runs (`isTargetEqualAcrossRuns` in selection.ts — scalar values, or value/index/split column contents) are dropped from `filteredTargets`; multi-mode selections stay visible |
| `selectedTargetKey` | which target is being charted |
| `candidateOverrides` | `Record<targetKey, Record<candidateId, boolean>>` — manual enable/disable toggles, fed into `matchColumnTargets`; a watcher on `targets` prunes entries whose target/candidate no longer resolves |
| `expandedTargetKeys` | which target rows have their candidate checklist expanded |
| `multiMode`, `multiKeys` | multi-value mode flag and its selected target keys; the mode never auto-exits — keys are pruned only when their target vanishes structurally, and the chart pads conflicting/missing runs (`partial` chip on the row) |
| `scalarChartType` | radar/PC-plot switch for multi-scalar charts; `radarAvailable` (Charts package deployed?) forces PC plot and hides the switch when radar can't render |
| `indexFilter`, `targetFilter` | list search boxes |
| `chartViewer` | last-created `DG.Viewer`, kept for the workspace-snapshot export |
| `historyHeight`, `sidebarWidth`, `chartHeight` | resizable layout dims |

### The reactive chain of `computed`s

Everything re-derives from `entries` + the selections:

1. `indexColumnsMap` / `splitColumnsMap` — validated Map form of the selections (a pick becomes invisible if its type toggle is off, unless it's the annotated default).
2. `targets` — `matchScalarTargets` + `matchColumnTargets`, sorted scalars first, then columns, alphabetically by display name — stable, so toggling candidates never reorders the list; the row displays the live enabled `coverage`. Changing an index selection immediately re-runs matching.
3. `selectedTarget`, `compatibleTargets` (multi-value candidates for the current anchor), `filteredTargets`.
4. `comparison` — the final result object: calls the right builder based on `target.kind` and `multiMode`, wrapped in `markRaw` (DataFrames must not be made reactive).
5. `indexRows` / `filteredIndexRows` — via `computeIndexRows` for the picker UI.
6. `entryStatuses` — exclusion chips.

### Notable mechanics

- **`addEntry`** dedupes by id and calls `applyAnnotatedDefaults`, which pre-fills the index/split pickers and the axis mode/units from the `comparison` annotation — but only where the user hasn't already chosen.
- **`addSelectedRuns`** reloads each checked history run fully (`historyUtils.loadRun`) before extraction, guarded by `isAddingRuns`.
- The **table input** (`ui.input.table`) resets its value in a `setTimeout` because doing it synchronously re-enters the Dart stream controller mid-dispatch (fix from `7ddbfc3b`).
- **Shift+click** on a target row enters multi-value mode grid-style: if the clicked target is suggestion-compatible with the current selection, both are kept; otherwise the click re-anchors.
- **Chart choice** in `renderComparison`:

  | Situation | Chart |
  |---|---|
  | scalar target | bar chart, split by `Run` |
  | multiple scalar targets | radar (default, Charts package `Radar` viewer — axes = targets, one polygon per run) or PC plot, via a small switch above the chart; PC plot is the silent fallback when Charts isn't deployed |
  | column target, key (string) index | bar chart, split by index, stacked by `Run` |
  | column target, numeric/datetime index | line chart: x = index, y = value column(s), split by `Run` (+ inner split column), `multiAxis: false` |
  | column target, relative timeseries (all tables configured) | same line chart over the converted `<index> (<unit>)` elapsed axis — no dedicated chart branch |
  | column target, independent points (all tables configured) | scatterplot: x = index pick, y = value, color = split ?? `Run` (runs move to marker shapes when split takes color); multi-value melts all values into one y column colored by target name, markers = `Run`, split unencoded |

  Column multi-value mode scales the minimum chart height by the number of value columns (`250px` per value); multi-scalar uses a flat minimum (one chart, not stacked series).
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
| Axis-mode warning chip | `relative timeseries not set` / `independent points not set` | `#fdf3e3` bg / `#8a6d3b` text (amber) |

## Entry points and tests

Besides the Model Hub menu, `History.tsx` has an `allowCompare` prop with a `compareSelected` handler — but that path opens the *old* compute-utils `RunComparisonView`; the new tool explicitly passes `allowCompare={false}` to its embedded History.

Tests live in `src/test/`, mirroring the source split: `run-comparison-matching.ts` (names/units, clustering), `run-comparison-selection.ts` (statuses, filters, index rows), and `run-comparison-builders.ts` (DataFrame building), with shared fixtures in `run-comparison-fixtures.ts`. User docs are at `help/compute/run-comparison.md`.
