# Charts — confirmed defects

Every entry below was raised by reading and then **adversarially re-verified** against the working tree
(2026-09-09) by an independent pass; refuted claims are in `findings.md` Appendix A, not here.
Severity is user-facing: **high** = fails/misleads in realistic use, **medium** = degraded quality or a
maintenance trap, **low** = polish/hygiene. Line numbers refer to the working tree (includes the
uncommitted radar changes).

## High

- **CH-01 [resolved] — Globe leaks a live WebGL render loop per closed viewer.**
  `src/viewers/globe/globe-viewer.ts:86-90` starts a self-scheduling `requestAnimationFrame` loop with no
  stored handle; `detach()` (:158-160) only unsubscribes subs. No `cancelAnimationFrame`,
  `renderer.dispose()`, or scene teardown anywhere. Every closed Globe keeps its WebGL context, scene,
  camera, orbit controls, and the viewer instance alive forever.

- **CH-02 [resolved] — Multiplot duplicates every data row; category split admits wrong rows.**
  `src/viewers/multiplot/utils.ts:145-154` (`getUniversalData`): two consecutive conditional pushes both
  fire when `condition` is undefined (every plot without `splitByColumnName`) → every filtered row pushed
  twice. On the split path `condition.value` is a single string, so `row[f] === value` and
  `value.includes(row[f])` both match → duplicates again, and `.includes` is a *substring* test, so
  category `'Alpha'` admits rows of `'Alphabet'`. Affects ClinicalCase (the main Multiplot consumer).

- **CH-03 [resolved] — Radar plots missing values as zeros and skews percentile bands.**
  `src/viewers/radar/radar-viewer.ts:445-448` maps the int-null sentinel to `0`; nothing filters nulls
  upstream (:313-342, :490). Float nulls (`FLOAT_NULL`/`null` → `Number(...)` ≈ 0) are not handled at
  all, and `getQuantile` (:561-575) filters only the int sentinel — the min/max percentile areas are
  computed over null-polluted data. With the default axis floor of 0, missing data renders as
  legitimate-looking values.

- **CH-04 [resolved] — Group Analysis runs and leaks after close; cached cell charts go stale.**
  `src/viewers/group-analysis/group-analysis-viewer.ts:145-149`: the `onRowsFiltered` subscription is
  never pushed to `this.subs`, so `updateGrid()` (grid + viewer rebuilds) keeps firing after the viewer
  is closed. `viewersStorage` (:108) is keyed by grid column + **viewport row** (`gc.gridRow`, not
  `tableRowIndex`) and never reset in `updateGrid()` (:301); each cached viewer also subscribes
  `onDartPropertyChanged` (:342-345) without ever unsubscribing → stale per-row charts after any
  filter/regroup, and unbounded subscription growth.

## Medium

- **CH-05 [resolved] — Tree/Sunburst call `this.detach()` on every render; the u2 scope is disposed after the first.**
  `src/viewers/tree/tree-viewer.ts:920-929`, `src/viewers/sunburst/sunburst-viewer.ts:440-449`: each
  render disposes the chart, calls `this.detach()`, and re-subscribes. Platform contract
  (js-api `viewer.ts:536-544`, `u2core/scope.ts:63-73`): `detach()` sets `isDetached`, disposes the u2
  effect scope **irreversibly** (property-signal refresh silently stops), and does **not** clear
  `this.subs` — the array regrows ~6 entries per render for the viewer's lifetime. Concretely lost after
  the first render: the base `onDataChanged → render` subscription (`echart-viewer.ts:59-71` is never
  re-added by `addSubs()`), so Tree stops reacting to in-place data edits; the constructor's
  `chartDiv` size sub is orphaned (re-init targets `this.root`, not the chartDiv). The correct pattern is
  `subs.forEach(unsubscribe); subs.length = 0` — never `detach()`.

- **CH-06 [resolved] — Word cloud ignores the filter.**
  `src/viewers/word-cloud/word-cloud-viewer.ts:128-139`: word sizes come from `strColumn.categories` and
  full-column counts; `this.filter` is never consulted, and there is no `onSourceRowsChanged` override —
  the `filter.onChanged` sub (:76) just re-renders identical unfiltered counts. Also O(categories)
  full-column `toList()` materializations per render.

- **CH-07 [resolved] — Word cloud clicks lose modifier keys.**
  `word-cloud-viewer.ts:183-188` passes the echarts param object as the event to
  `selection.handleClick` (hidden by `@ts-ignore`); `ctrlKey/shiftKey/metaKey` read as undefined, so
  every click is a plain replace-selection. Peer viewers pass `params.event.event`.

- **CH-08 [resolved] — Sankey renders ALL rows when the filter excludes everything.**
  `src/viewers/sankey/sankey.ts:143-147`: `selectedIndexes.length > 0 ? selectedIndexes : all rows` —
  the empty-filter case silently inverts into the unfiltered graph.

- **CH-09 [resolved] — Sankey drag breaks on node names with CSS-special characters.**
  `sankey.ts:336` writes the node name (spaces→underscores only) as a class; `:350` builds
  `select('text.' + name)` — digit-leading names (`"2020"`) throw `SyntaxError` mid-drag; names with
  `.`/`:`/`#` silently desync the label from the dragged rect.

- **CH-10 [resolved] — Radar: the first row is unclickable.**
  `radar-viewer.ts:180-181`: `idx = parseInt(name.replace(/\D/g,'')) - 1; if (idx)` — `"row 1"` → 0 →
  falsy → `currentRowIdx` never set for row 1, while the tooltip path (:189-196) works, making the
  inconsistency visible.

- **CH-11 [resolved] — Error messages stack into duplicate banners (radar, tree, sunburst, word cloud, globe).**
  `_showMessage` appends unconditionally and the failing-render paths return before the corresponding
  remove: `radar-viewer.ts:536-541`, `tree-viewer.ts:896/900`, `sunburst-viewer.ts:427/431`,
  `word-cloud-viewer.ts:110-120`, `globe-viewer.ts:207-210`. Renders fire from ~8 subscriptions plus
  resize, so a viewer in the error state accumulates one banner per event. Bonus (tree/sunburst): the
  early error-return also skips the dispose/detach block, leaving that render's duplicate chart handlers
  and subs alive (the one path where handler double-binding persists).

- **CH-12 [resolved] — Tree color scale drops zero-valued aggregates.**
  `src/utils/tree-utils.ts:133-136` (`if (!value) continue`): 0 is excluded from min/max, so zero-valued
  nodes clamp to a scale extreme and the whole color range shifts for aggregations that legitimately
  produce 0 (sum, avg, variance, #selected).

- **CH-13 — Category colors come from a stale cross-table cache.**
  `tree-utils.ts:8, 82-97`: static LRU keyed by bare column *name*, shared across all dataframes, never
  invalidated — and it caches the *aggregated* column, whose category set changes with every filter. On
  the default Tree path (`inherit` undefined, :746) same-named columns across tables collide and changed
  categories yield `undefined` → echarts palette fallback.

- **CH-14 [resolved] — One failed molecule image permanently kills molecule labels.**
  `tree-viewer.ts:766-771` / `sunburst-viewer.ts:297-302`: `moleculeRenderQueue = queue.then(...)` with
  no `.catch` — the first rejection (e.g. Chem missing, bad SMILES) poisons the chain; every later
  molecule label is skipped with an unhandled rejection each. `TreeUtils.getMoleculeImage`
  (`tree-utils.ts:223-229`) has no availability check and no caller catches (verified all 4 sites).

- **CH-15 [resolved] — Table rebind keeps the old dataframe's subscriptions.**
  `radar-viewer.ts:210-250`, `chord-viewer.ts:143-146`, `sankey.ts:112-114`,
  `word-cloud-viewer.ts:76-77`: `onTableAttached` pushes subs without dropping prior ones; a rebind
  (`table` property change / project rebind re-fires attach) leaves the old dataframe subscribed —
  duplicate renders plus a retained reference to the old dataframe.

- **CH-16 [resolved] — Multiplot ships dev scaffolding as product.** *(timeline red/green colors deliberately kept — restoring the commented-out color logic is a ClinicalCase-visible behavior change, out of scope per plan)*
  `multiplot.ts`: 15 live `console.*` (:84,85,104,134,390,418,488,496,589,594,610,614,618,622,890 —
  verified present in `dist/package.js`) firing on every property change/click/zoom/brush; dev
  properties visible in the property panel — `paramA` default `'string inside'` (:42), `paramOptions`
  default `'none22'` (:60); timeline items hardcode `fill:'red'` (:939) / `'green'` (:1002) with the
  real color logic commented out (:941-945, :998-999); the constructor's `ui.onSizeChanged` subscription
  is unmanaged and there is no `detach()` at all — the echarts instance is never disposed (:92-100).

- **CH-17 [resolved] — Surface plot: wrong background default; dataShape diverges from filtered data.** *(the `grok.shell.error/warning` balloons from `onTableAttached` on 1-2-column tables are left as-is)*
  `surface-plot.ts:66`: `this.int('backgroundColor', 0xFFF)` = 4095, assigned raw (unparseable color) to
  `option.backgroundColor` on every render (:257) until the property is first changed (:221). `:272-273`
  `dataShape = [√n, √n]` uses the **unfiltered** length while `series[0].data` is filtered (:277) —
  shape and data diverge under any filter (garbled surface). Also `grok.shell.error/warning` fire from
  `onTableAttached` (:177, :180) on small tables.

- **CH-18 [resolved] — Timelines: inside x-zoom is inert; type-mismatch warning fires per row.**
  `src/viewers/timelines/echarts-options.ts:30-34`: `{type:'inside', xAxisIndex:[1,2]}` targets axes
  that don't exist (there is one xAxis) → wheel/pinch x-zoom does nothing, and the dead entry still
  occupies a `zoomState` slot (`timelines-viewer.ts:465-469, 102-107`). `isSameDate`
  (`timelines-viewer.ts:488-495`) calls `grok.shell.warning` from per-row click/tooltip predicates
  (:126-138, :217-220) — one balloon per row on column-type mismatch.

- **CH-19 [resolved] — `lodash` is imported but not a dependency.**
  `radar-viewer.ts:11`, `tree-viewer.ts:12`, `sunburst-viewer.ts:11` import `lodash`;
  `package.json` declares only `@types/lodash` (dev). Resolves through hoisting today; a clean isolated
  install or hoisting change breaks the build, and the version is unpinned.

## Low

- **CH-20 [resolved]** — Tree `onPropertyChanged`: missing `break` after `case 'initialTreeDepth'`
  (`tree-viewer.ts:510-517`) falls into `'symbolSize'`. Impact currently neutralized by the
  unconditional reassign in `_render` (:913-916) — net effect is a wasted second full render; still a
  live trap.
- **CH-21** — Dead null-category styling: `tree-utils.ts:171-175` grey `itemStyle` immediately
  overwritten at :177-179.
- **CH-22** — `aggToStat` uses `eval` over a fixed lookup (`src/utils/utils.ts:26-47`). No injection
  (closed key set) but it already disables terser mangling of the enclosing scope (verified in
  `dist/package.js`); a minifier change silently breaks Tree size/color aggregation.
- **CH-23** — Chord silently and persistently overwrites the user's `sortBy` with `'alphabet'`
  mid-render via `@ts-ignore` (`chord-viewer.ts:286-287`); no recursion/spam (guard at :157), but the
  layout-persisted property changes without the user's intent.
- **CH-24** — `super.detach()` skipped in sankey (:221-223), globe (:158-160), word-cloud (:100-102),
  chord (:162-164), group-analysis (:126-128): `isDetached` never set, u2 scope never disposed,
  re-entrant detach possible. Subs are drained manually, hence low.
- **CH-25** — Two source `debugger` statements (`multiplot.ts:344, 971`) — stripped by terser in the
  production bundle (verified absent in `dist/`), so source-hygiene only. Unreachable second return
  :400; undeclared `this.count` (:991-993, dead behind `const overlap = false`); double-push of
  show/hide toggles into `typeComboElements` (:754) — harmlessly unreachable extra entries.
- **CH-26** — Dead code: `src/viewers/multiplot/timeLinesRender.ts` (123 lines, 100% commented, never
  imported); `layout.ts:1-3` `MPlotLayout2` stub; `utils.ts:190-191` empty `splitToMultipleSeries`;
  `utils.ts:165-169` `getBitByIndex32` re-implements `BitSet.get`; `utils.ts` `normalize100` unused;
  the whole `src/deprecated/` tree (~650 loc) is imported by nothing.
- **CH-27** — Group Analysis polish: user-visible column named `pValue(AGE` — unbalanced parenthesis
  (`group-analysis-viewer.ts:269`); first chart cell renders empty until a later cellPrepare
  (:334-349 inverted cache logic); no numeric-type gate on the T-test column choice.
- **CH-28** — README drift: Timelines documents `colorByColumnName` (README.md:24) — the property is
  `colorColumnName`; `autoSize` undocumented; Tree section documents nonexistent `edgeShape` and
  `expandAndCollapse` (README.md:196-197); `left`/`right` documented as common properties but removed
  (`echart-viewer.ts:48-52`).
- **CH-29** — Flag cell renderer is a stub shipped as a registered renderer
  (`src/renderers/flag-cell-renderer.ts`): renders the literal text "flag" in hardcoded black, ignoring
  cellStyle/theme; paired detector fires on any string column literally named `flag`.
- **CH-30** — Radar `showCurrentRow` description reads "Hides max and min values" (copy-paste;
  `radar-viewer.ts:61`); word-cloud's persisted property is named `columnColumnName`
  (`word-cloud-viewer.ts:40`) — frozen by saved layouts.
- **CH-31** — `test-report.csv` (test artifact) sits untracked at the package root — should be
  gitignored; `detectors.js` `detectMagnitude` mutates `col.semType` inside the detector.
