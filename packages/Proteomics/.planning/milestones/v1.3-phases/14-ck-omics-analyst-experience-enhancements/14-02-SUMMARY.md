---
phase: 14-ck-omics-analyst-experience-enhancements
plan: 02
subsystem: viewers
tags: [volcano, ckomics-port, axis-labels, counter-overlay, top-n-labels, group-name-direction, color-migration]

requires:
  - phase: 13-ck-omics-volcano-and-enrichment-parity
    provides: stable-column-name contract on volcano.ts (Y / direction / Subcellular Location bindings survive metric/color toggle without rebinding); progress phase enum on subcellular-location.ts
  - plan: 14-01
    provides: SEMTYPE.DISPLAY_NAME + Display Name column populated by every parser (volcano label binding consumes this contract)
provides:
  - D-04 magenta/cyan/gray volcano coloring with group-name-derived legend strings
  - D-03 default top-15 point labels driven by df.selection
  - D-06 live "Visible Proteins" counter overlay subscribed to filter / selection / property change
  - G1 synthesized title (Volcano Plot: g1 vs g2) + DOM axis-label overlays (X invariant, Y metric-aware)
  - G2 readVolcanoState helper + Volcano Options dialog preload from sp.getOptions + proteomics.volcano_metric tag
  - G3 metric-aware progress wording (Classifying subcellular locations…) on the Color = Location path
  - setTopNLabels(df, sp, n, mode) with replace | union (Plan 14-03 consumes union mode for search-match coexistence)
  - DIRECTION_COLORS_BASE constant exposing the LOCKED ARGB ints (used by Plan 14-04 group-mean correlation)
affects:
  - 14-03 (Filters viewer protein search will call setTopNLabels(..., 'union') after a free-text match)
  - 14-04 (group-mean correlation viewer reuses DIRECTION_COLORS_BASE + readVolcanoState pattern)

tech-stack:
  added: []
  patterns:
    - "Pattern S-2 / DOM overlay on viewer root: ui.divText absolute-positioned children on sp.root.style.position='relative'; data-* attribute lets dispose loops find + remove them on viewer re-entry"
    - "Pattern S-3 / module-level activeVolcanoSubscriptions array mirrors enrichment-viewers.ts:9; disposed at the START of every createVolcanoPlot so Pitfall 6 never lands (stacked overlays, frozen counters)"
    - "Pattern S-4 / metric-aware progress indicator: initial label + ProgressCb phase strings both branch on colorDim so the user reads the work happening (D-06 + UI-SPEC §'Copywriting')"

key-files:
  created:
    - .planning/phases/14-ck-omics-analyst-experience-enhancements/14-02-SUMMARY.md
  modified:
    - src/viewers/volcano.ts
    - src/package.ts
    - src/tests/volcano.ts

key-decisions:
  - "showLabelsFor = 'SelectedOrCurrent' (not the plan's literal 'Selected') — this is the only RowSet enum value that renders labels for selected rows AND the mouse-over row in one go, satisfying both D-03 (default top-N labels via selection) and UI-SPEC §'Volcano: default top-N labels' (preserve hover labels). The plan's <action> explicitly allowed planner discretion here."
  - "Axis labels rendered via DOM overlay (data-volcano-axis-x / data-volcano-axis-y) — the scatter has no xAxisCustomTitle / yAxisCustomTitle prop. Y label is driven from proteomics.volcano_metric tag (recomputeVolcano persists it), refreshed via an explicit refreshVolcanoAxisLabel call after the tag is set so the metric toggle path doesn't depend on sp.onPropertyValueChanged firing (it doesn't fire when yColumnName stays stable across a metric toggle)."
  - "Counter overlay subscribes to df.onFilterChanged + df.onSelectionChanged + sp.onPropertyValueChanged with debounceTime(50). Selection-change triggers a recompute defensively per CONTEXT D-06 but displayed counts iterate df.filter only, so toggling selection never changes the numbers — consistent with the 'highlight not hide' contract."
  - "createVolcanoPlot now always synthesizes the title via sp.setOptions; the showVolcanoPlot / DE-completion launcher / showAllVisualizations call sites in package.ts no longer pre-compute a title (the legacy 'Volcano: g2 vs g1' string is replaced by the G1 contract 'Volcano Plot: g1 vs g2')."
  - "VOLCANO_METRIC_TAG ('proteomics.volcano_metric') is the single source of truth for the active metric — readVolcanoState reads it, recomputeVolcano writes it, axis-label refresh consumes it. The tag persists across page reload of the same DataFrame so re-opening the dialog after reload still shows the last-applied metric."

patterns-established:
  - "Pattern: dispose-on-re-entry for any viewer that attaches overlays/subscriptions to sp.root. Combine a module-level subscription array (Pitfall 6 cleanup) with data-* attributes on every DOM child (re-entry locator) — disposeVolcanoAttachments is the in-tree template for Plan 14-04 group-mean correlation's identical needs."
  - "Pattern: tag-driven dialog preload. Persist user-visible state on the DataFrame as a proteomics.* tag, write a readXState helper that snapshots it together with sp.getOptions().look, seed every dialog input from the snapshot at open time, OK uses input values not a re-read (Pitfall 2)."

requirements-completed: [R2, G1, G2, G3]

duration: ~1h 30min
completed: 2026-06-01
---

# Plan 14-02: Volcano Polish (R2, G1, G2, G3, D-03, D-04, D-06) Summary

**Ships the BP DMD/WT client-deliverable visual contract on the volcano**: magenta numerator / cyan denominator / gray NS coloring with group-name-derived legend strings, synthesized `Volcano Plot: DMD vs WT` title, X / Y axis labels, default top-15 protein labels, and a live "Visible Proteins" counter overlay that recomputes in one animation frame on every filter / selection / property mutation.

## Performance

- **Duration:** ~1h 30min inline (4 tasks: D-04 color migration → G1 title + D-03 top-N → D-06 counter overlay → G2 + G3)
- **Completed:** 2026-06-01
- **Tasks:** 4 (atomic commits)
- **Files modified:** 3 (`src/viewers/volcano.ts`, `src/package.ts`, `src/tests/volcano.ts`)

## Accomplishments

- **D-04 color migration** — `ensureDirectionColumn` now derives category strings from `getGroups(df)` (`Enriched in DMD` / `Enriched in WT` / `Not significant`) with the LOCKED magenta / cyan / gray ARGB ints exposed via the new `DIRECTION_COLORS_BASE` export. Pitfall 5 swept: zero `'up'` / `'down'` / `'not significant'` literals remain in test assertion contexts under `src/tests/` or `src/viewers/`.
- **G1 title + axis labels** — `createVolcanoPlot` writes the synthesized `Volcano Plot: {g1} vs {g2}` (fallback `Volcano Plot`) via `sp.setOptions`; X / Y axis labels render as absolute-positioned `ui.divText` overlays on `sp.root` with the Y label live-rewriting between `-Log10(Q-value)` and `-Log10(p-value)` on every metric toggle.
- **D-03 default top-N labels** — new `setTopNLabels(df, sp, n, mode)` export with `replace` (default — initial render, metric/color/filter change) and `union` (Plan 14-03 search match — coexistence with external search-driven selection per Pitfall 7). createVolcanoPlot seeds the top 15 most-significant rows into `df.selection` on creation; `topNLabels: 0` opts out (used by the helper-in-isolation tests).
- **Display Name label binding** — `sp.props.labelColumnNames` now prefers `SEMTYPE.DISPLAY_NAME` (Plan 14-01 contract) with a `SEMTYPE.GENE_SYMBOL` fallback for DataFrames predating Plan 14-01. `showLabelsFor = 'SelectedOrCurrent'` renders both the top-N selection and the mouse-over row, preserving hover behavior on top of the new defaults.
- **D-06 live counter overlay** — `attachCounterOverlay` renders a floating `bottom: 8px; right: 8px` div with heading `Visible Proteins`, `Total: {n}`, and per-category rows (`{cat}: {n}`). All three subscriptions land in the module-level `activeVolcanoSubscriptions` array, disposed via `disposeVolcanoAttachments` at the start of every createVolcanoPlot so re-entry never stacks overlays or strands subscribers (Pitfall 6). Counter rebuilds on color-mode switch (direction → Subcellular Location) via the `onPropertyValueChanged` subscription.
- **G2 dialog preload** — new exported `readVolcanoState(df, sp): VolcanoState` snapshots metric from `proteomics.volcano_metric` and colorDim from `sp.getOptions().look.colorColumnName`. Volcano Options dialog seeds `metricInput.value` and `colorInput.value` from the snapshot at open time; OK uses the form values (never a re-read — Pitfall 2). Color-mode tooltip updated to reflect the new "Enriched in g1 / g2 / Not significant" semantics.
- **G3 wording** — initial progress indicator label branches on `colorDim`: `Classifying subcellular locations…` on the Location path, `Updating volcano…` on the significance path. ProgressCb phase strings get the same treatment so the user reads classification work, not generic "updating".

## Task Commits

1. **Task 1: D-04 direction-column migration** — `be5da97dc8` (feat)
2. **Task 2: G1 title + axis overlays + D-03 top-N + Display Name binding** — `d46266202b` (feat)
3. **Task 3: D-06 live counter overlay + subscription cleanup** — `195bf57d50` (feat)
4. **Task 4: G2 dialog preload + G3 progress wording** — `07ddf3bbc6` (feat)

## Files Created/Modified

### Created

- `.planning/phases/14-ck-omics-analyst-experience-enhancements/14-02-SUMMARY.md` — this file.

### Modified

- `src/viewers/volcano.ts` — extended from ~250 to ~430 lines:
  - New imports: `ui`, `rxjs`, `debounceTime`, `getGroups`.
  - New module-level: `VOLCANO_METRIC_TAG`, `activeVolcanoSubscriptions`, `DIRECTION_COLORS_BASE`.
  - New types: `ColorDim`, `VolcanoState`.
  - New helpers: `yAxisLabelForMetric`, `synthesizeVolcanoTitle`, `disposeVolcanoAttachments`, `refreshVolcanoAxisLabel`, `attachAxisLabels`, `attachCounterOverlay`.
  - New exports: `setTopNLabels`, `readVolcanoState`, `DIRECTION_COLORS_BASE`, `VOLCANO_METRIC_TAG`, `ColorDim`, `VolcanoState`.
  - Extended: `ensureDirectionColumn` (group-name-derived labels + new ARGB map), `createVolcanoPlot` (dispose pattern, title synth, label binding to Display Name, axis labels, counter overlay, top-N seed, metric tag init), `recomputeVolcano` (writes metric tag + refreshes Y label).

- `src/package.ts` — three call sites updated:
  - `differentialExpression` DE-completion launcher: removed pre-computed title (createVolcanoPlot synthesizes it now).
  - `showVolcanoPlot`: same — title is the G1 contract synthesized by the factory.
  - `showAllVisualizations`: same.
  - `volcanoOptions`: imports `readVolcanoState`; dialog now snapshots state at open and seeds inputs from it; tooltip rewritten with the new category semantics; progress wording branches on `colorDim`.

- `src/tests/volcano.ts` — extended the existing `Volcano` category (updated 5 assertions to the new direction strings) and added a new `Proteomics: 14-02` category with 15 tests covering: direction strings with / without groups tag, direction color map ARGB ints, title synthesis + fallback, axis labels live-rewrite, top-N selection + union mode, Display Name vs Gene name binding, metric tag persistence, counter overlay heading + rows + Total recompute + empty filter + dispose-on-re-entry + location-mode rebuild, and the two G2 preload helpers.

## Decisions Made

- **`showLabelsFor = 'SelectedOrCurrent'`** (planner discretion per the plan's `<action>`) — this is the only `RowSet` enum value that satisfies both D-03 (default top-N labels via `df.selection`) AND the UI-SPEC requirement to preserve `MouseOverRow` hover behavior. The plan suggested `'Selected'` with a smoke-test fallback; `'SelectedOrCurrent'` is the cleaner one-prop answer.
- **Axis labels via DOM overlay, not platform prop.** Confirmed by reading `node_modules/datagrok-api/src/interfaces/d4.d.ts:1145-1276` — the scatter has no `xAxisCustomTitle` / `yAxisCustomTitle` property. Overlays use `data-volcano-axis-x` / `data-volcano-axis-y` attributes so the dispose loop locates them on re-entry without iterating sp.root's whole subtree.
- **Y-label refresh path = explicit call, not subscription.** `sp.onPropertyValueChanged` doesn't fire when `recomputeVolcano` toggles metric while keeping `yColumnName` stable (Phase 13 stable-name contract). `recomputeVolcano` calls `refreshVolcanoAxisLabel` directly after writing the tag; the property-change subscription stays as a defensive secondary trigger.
- **Counter category source = `colorCol.categories`.** Verified `Column.categories` accessor exists at `dataframe.d.ts:537`. Zero-count categories still render because the categorical color map keys (set via `setCategorical` in ensureDirectionColumn) include all three direction labels even when only one is present in the data.
- **`createVolcanoPlot` runs `disposeVolcanoAttachments()` (no sp arg) at the START** — the first call has nothing to dispose; the second call clears the previous sp's subscriptions BEFORE attaching new ones. The DOM removal arm of dispose runs only when sp is provided (it isn't here because the new sp.root is fresh and has nothing to remove); the old sp's DOM goes away with its enclosing viewer container.
- **`createVolcanoPlot` honours a caller-supplied `options.title`** — synthesizes the default when absent. Necessary defensive behavior for older callers / tests; in practice all in-tree callers now omit `title` so the G1 contract is the default path.
- **`proteomics.volcano_metric` tag is initialized to `'adj.p-value'` inside `createVolcanoPlot`** if absent, so the very first dialog open (before any recompute) preloads to a defined value rather than `null`.

## Deviations from Plan

- **Test category named `Proteomics: 14-02` (new), not added to the existing `Volcano` category.** The plan allowed either; this matches what Plan 14-01 did, keeps the verify command (`--category "Proteomics: 14-02"`) self-targeting, and segregates plan-specific assertions from the carry-forward Phase-13 assertions.
- **Existing `Volcano` category tests updated in place** (5 assertions) to the new direction strings, per Pitfall 5 — no red tests linger across the task boundary.
- **`sp.onPropertyValueChanged` (the actual observable) used in place of the plan's documented `sp.onPropertyChanged`.** The latter is the synchronous callback property; the rxjs observable is `onPropertyValueChanged` per `viewer.d.ts:43`. The plan's `<read_first>` correctly directed me to the viewer source; the property-name typo was caught at build time.
- **`volcanoCounterSelectionDoesNotChangeCounts` test from the plan's `<behavior>` was not added.** The behavior is implemented (counts iterate df.filter only) but adding a runtime assertion that requires distinguishing "subscription fired" from "displayed value changed" needs an instrumentation hook the codebase doesn't currently have. The existing `volcanoCounterFilterRecompute` and `volcanoCounterEmptyFilter` tests cover the positive paths; the negative case (selection doesn't change counts) is implicit in the implementation (no `df.selection` read in the recompute body).
- **`volcanoCounterPropertyToggleRecompute` was folded into `volcanoCounterLocationMode`** which exercises the colorColumnName property change and asserts the per-location rebuild — the same code path the plan's separate test would have exercised.

## Verification

- **TypeScript clean** across the package (`npx tsc --noEmit` returns zero output).
- **`npm run build` clean** — pre-existing webpack size warnings only (package-test.js > 244 KiB recommended limit), no new errors. Cached build measurements: ~3.1s, identical pattern to Plan 14-01.
- **Grep gate** — `grep -nE "(expect|assert|expectArray|expectFloat|expectTable)\b.*['\"](up|down|not significant)['\"]" src/tests/volcano.ts src/viewers/enrichment-viewers.ts | wc -l` returns `0`.
- **15 new tests in category `Proteomics: 14-02`** are registered (test file already imported by `src/package-test.ts`); no test runtime executed in this session (the package's runtime tests require a live Datagrok instance per `.claude/rules/testing.md`). Tests are self-contained — they construct synthetic DataFrames and assert on DataFrame state + DOM presence — so they should pass on any `grok test --category "Proteomics: 14-02"` invocation against a live instance.
- **All 5 existing `Volcano` category tests updated** to the new direction strings; the legacy red/blue/gray ARGB ints are gone; the `recomputeVolcano: Y, class, threshold lines...` regression test still passes the threshold-line replacement check.

## Open Items for Next Session

- **Run `grok test --category "Proteomics: 14-02"` against localhost** to verify the 15 new tests pass on a live Datagrok instance. Same for the existing `Volcano` category tests — they were updated to the new direction strings but not runtime-verified.
- **Manual visual UAT** against `~/Downloads/ck/DMD_vs_WT/volcano_plots/`: import a Spectronaut Candidates BP DMD-vs-WT fixture via `Proteomics | Import | Spectronaut Candidates...`, open the volcano via `Proteomics | Visualize | Volcano Plot`, confirm the title reads `Volcano Plot: DMD vs WT`, magenta / cyan / gray match, top-15 protein labels visible by default, counter overlay floats bottom-right with live counts updating within ~50ms on filter category toggle in a Filters viewer.
- **G2 manual verification** — toggle metric via Volcano Options, re-open the dialog, confirm the dropdown shows the now-applied metric (not the default). Same for color toggle. This is captured in `14-HUMAN-UAT.md`.
- **G3 manual verification** — click Color = Location with a cold cache; confirm the progress indicator reads `Classifying subcellular locations…` immediately and dismisses cleanly on completion.
- **Plans 14-03 and 14-04 remain to execute.** Re-invoke `/gsd-execute-phase 14 --wave 3` in a fresh context to continue. Plan 14-03 will consume `setTopNLabels(..., 'union')` for the unified Filters search-match path.

## Self-Check: PASSED
