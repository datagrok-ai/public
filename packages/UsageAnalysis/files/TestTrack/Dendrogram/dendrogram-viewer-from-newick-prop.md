---
feature: dendrogram
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [dendrogram.cp.viewer-from-newick-prop]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - dendrogram-viewer-from-newick-prop-spec.ts
gate_verdicts:
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T16:57:43Z
    spec_runs:
      - spec: dendrogram-viewer-from-newick-prop-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 68
        failure_keys: []
---

# Dendrogram — Viewer from literal newick property + full property-panel surface

Opens the Dendrogram viewer through `addViewer('Dendrogram', {newick: '...'})`
on a small synthetic 4-leaf DataFrame, then exercises every property-panel
surface on the viewer so each property is observed to either restyle the
canvas (style props) or rebuild the view (data props) — verifying the
viewer's property-change dispatch and its table-attached lifecycle
(re-populating the `newickTag` choice list when a new table is attached).

Scenario 1 is the smoke load (verifies the viewer mounts and the canvas
renders four leaves from the literal newick). Scenario 2 is the regression
sweep over the 17 property surfaces grouped by dispatch behavior (the four
data props that trigger a full rebuild, then the 13 style/behavior props that
trigger a restyle only) plus the `newickTag` repopulation that proves the
`onTableAttached` lifecycle.

## Setup

- A clean Datagrok view (no preloaded tables, no other dendrograms attached).
- The Dendrogram package is installed and registered (the
  `[name="viewer-Dendrogram"]` viewer type is registered via
  `_package.registerViewer`).
- A small synthetic 4-leaf DataFrame named `leaves4` is built with two
  columns:
  - `leaf` (string) — the four leaf names `['a', 'b', 'c', 'd']` matching
    the literal newick used in Scenario 1.
  - `value` (int) — numeric column `[1, 2, 3, 4]` used to drive the
    `colorColumnName` + `colorAggrType` rebuilds in Scenario 2.
- The literal newick used throughout is `((a:1,b:1):1,(c:1,d:1):1);` — a
  binary 4-leaf tree whose leaf names match the `leaf` column values.
- DataFrame tags are populated with two distinct newick strings so the
  `newickTag` choices are non-empty after `onTableAttached`:
  - tag `.newick` = `((a:1,b:1):1,(c:1,d:1):1);`
  - tag `.newick-alt` = `(((a:1,b:1):1,c:1):1,d:1);` (a left-leaning
    4-leaf newick — same leaves, different topology so the rebuild is
    observable.)

## Scenarios

### Scenario 1 — Viewer mounts and renders from literal newick property (smoke)

Steps:

1. Build the `leaves4` DataFrame programmatically (4 rows, columns `leaf`
   and `value` per Setup); set DataFrame tags `.newick` and `.newick-alt`
   to the two newick strings from Setup.
2. Open the DataFrame in a `TableView` via `grok.shell.addTableView(df)`.
3. Add the Dendrogram viewer via
   `tv.addViewer('Dendrogram', {newick: '((a:1,b:1):1,(c:1,d:1):1);',
   nodeColumnName: 'leaf'})`. Wait for the viewer host to mount
   (`[name="viewer-Dendrogram"]` present in DOM).

   * Expected result: a `[name="viewer-Dendrogram"]` element appears in
     the active view. Its inner `canvas` is rendered and non-empty
     (`canvas.width > 0`, `canvas.height > 0`). No console error is
     emitted during mount.
4. Assert the viewer's `treeNewick` resolves to the literal property per
   the priority order (property > `newickTag` > `.newick` tag >
   `NEWICK_EMPTY`) per
   `public/packages/Dendrogram/src/viewers/dendrogram.ts#L309`. Read it
   off the viewer instance via the viewer's `props.newick` (the literal
   string set in step 3).

   * Expected result: `viewer.props.newick === '((a:1,b:1):1,(c:1,d:1):1);'`
     — the literal property won the priority race against `.newick` and
     `.newick-alt` tags.
5. Assert that the four leaf names from the newick are present in the
   rendered tree's leaf set. The viewer exposes leaf names via its
   internal tree model; assert through `TreeHelper.getLeafList(rootNode)`
   on the resolved tree (or via the canvas-text observable if
   `showLabels` is on — Scenario 2 covers that path explicitly).

   * Expected result: the leaf set is `{a, b, c, d}` (set semantics, no
     duplicates). No `Non unique key tree leaf name` exception is raised.

### Scenario 2 — Property panel exercises all 17 viewer properties + newickTag lifecycle

Continuing from the Scenario 1 state (viewer mounted, DataFrame attached,
literal newick property in effect).

Steps:

1. Open the property panel for the Dendrogram viewer (right-click viewer
   header → `Properties`, or
   `grok.shell.o = viewer`). Wait for the property panel pane to render.

   * Expected result: the property pane lists the Dendrogram viewer's
     property categories — `Data` (newick, newickTag, nodeColumnName,
     colorColumnName, colorAggrType), `Style` (lineWidth, nodeSize,
     showGrid, mainColor, lightColor, currentColor, mouseOverColor,
     selectionsColor, showLabels, font, step), `Behavior` (stepZoom,
     showTooltip).

   Selector pinning note: per `.claude/skills/grok-browser/references/dendrogram.md`
   `## dendrogram-viewer-property-panel`, the property-row anchors use
   the framework's name-mangling rule — `nodeColumnName` mangles to
   `prop-node` (not `prop-node-column-name`) and `colorColumnName`
   mangles to `prop-color` (not `prop-color-column-name`). Tests that
   pin `prop-node-column-name` / `prop-color-column-name` will fail to
   locate the row. The full mangling matrix lives in the ref doc's
   `## Selector validation matrix`.

2. **newickTag rebuild path (data prop; exercises the table-attached
   choice-repopulation and the property-changed rebuild branch).** Set
   `newickTag` = `.newick-alt` via the property panel input. Wait for the
   canvas redraw.

   * Expected result: the `newickTag` choice list is non-empty and
     contains at least `.newick` and `.newick-alt` (the two DataFrame
     tags whose name starts with `.`). Selecting `.newick-alt` causes the
     view to rebuild — the new tree's topology differs from Scenario 1
     (left-leaning vs balanced). Leaf set still `{a, b, c, d}` (set
     semantics preserved). No console error.

   Then reset `newickTag` to empty (or to `.newick`) so subsequent steps
   operate against the original literal-newick view — the literal
   `newick` property still wins by priority order, so the canvas snaps
   back to the Scenario 1 topology.

3. **nodeColumnName rebuild path.** Set `nodeColumnName` via property
   panel — pick the `leaf` column (it is already set; toggle to a
   different column and back to verify the rebuild fires both ways).

   * Expected result: the property change is reflected in
     `viewer.props.nodeColumnName`; the view rebuilds; the
     current/mouseOver/selection sync now binds against the chosen
     column. No console error.

4. **colorColumnName + colorAggrType rebuild path.** Set
   `colorColumnName` to `value`. Set `colorAggrType` to `avg` (default),
   then cycle through the remaining choices `min`, `max`, `med`,
   `total-count` per
   `public/packages/Dendrogram/src/viewers/dendrogram.ts#L142`. After
   each change wait for the canvas redraw.

   * Expected result: each `colorAggrType` choice is accepted (no
     `Unsupported aggregation` exception); the view rebuilds per the
     `colorColumnName / colorAggrType` rebuild branch of
     `onPropertyChanged`; the canvas re-renders. No console error.

5. **lineWidth restyle path (style-prop sweep, exercises restyle branch
   of `onPropertyChanged`).** Set `lineWidth` to `0` (lower bound), then
   to `16` (upper bound per slider range `0..16, step 0.1`), then to
   `2.5` (mid-range).

   * Expected result: each value is accepted by the slider; the canvas
     restyles after each change (no full rebuild); no console error. The
     viewer instance reflects the new value on `viewer.props.lineWidth`.

6. **nodeSize restyle path.** Set `nodeSize` to `0`, `16`, then `4.0` —
   same sweep shape as `lineWidth`.

   * Expected result: each value accepted; restyle observed; no console
     error.

7. **showGrid restyle path.** Toggle `showGrid` off then on (default is
   tracked per the property's declared default).

   * Expected result: the canvas restyles after each toggle (the
     background grid disappears / reappears, observed via the canvas
     image-data fingerprint or simply via no-console-error on the
     restyle event). `viewer.props.showGrid` reflects the toggle.

8. **Color sweep (mainColor / lightColor / currentColor / mouseOverColor
   / selectionsColor).** For each of the five color props, set the value
   to `#FF0000` (red), then back to its declared default. After each
   change wait for the canvas to restyle.

   * Expected result: every color prop accepts the hex value via the
     property panel's color input; the canvas restyles after each change
     (no rebuild); the viewer's `props.<colorName>` reflects the new
     value; no console error.

9. **showLabels + font restyle path.** Toggle `showLabels` on. Then set
   `font` to `12pt monospace` (the default is `monospace 10pt` per
   `dendrogram.ts#L169`). Toggle `showLabels` off.

   * Expected result: with `showLabels` on, leaf labels (`a`, `b`, `c`,
     `d`) are rendered onto the canvas (verifiable via canvas text or
     simply via no-console-error + viewer prop reflection). The `font`
     change restyles the rendered labels. Toggling `showLabels` off
     removes labels on the next restyle. No console error.

10. **stepZoom behavior prop.** Set `stepZoom` to `-4`, `4`, `0.5` (range
    `-4..4, step 0.1`). After each change perform a small wheel-zoom
    interaction on the viewer canvas (Ctrl+wheel up by one notch) and
    verify the zoom step magnitude scales with the new `stepZoom` value.

    * Expected result: each value accepted; the wheel-zoom step size
      reflects the new `stepZoom`; no console error.

11. **showTooltip behavior prop.** Toggle `showTooltip` off. Hover a leaf
    on the canvas (move mouse to the leaf node's screen coordinate).
    Toggle `showTooltip` on. Hover again.

    * Expected result: with `showTooltip` off, no tooltip is shown on
      hover (no `[id^="tooltip"]` or equivalent tooltip element appears).
      With `showTooltip` on, hovering a node produces a tooltip
      containing the node name plus index/range (per
      `dendrogram.ts#L627`). `viewer.props.showTooltip` reflects the
      toggle. No console error.

12. **step style prop.** Set `step` to `0`, `64`, `1.5` (range
    `0..64, step 0.1`). After each change wait for restyle.

    * Expected result: the vertical leaf-row spacing reflects the new
      `step` value; the canvas restyles; no console error.

13. **No-console-error invariant (across the whole property sweep).** At
    the end of the scenario assert that no error was logged by the
    `onPropertyChanged` dispatcher during steps 2–12 (no
    `Unsupported property` / `Unsupported aggregation` /
    `Unsupported semType` / `Non unique key tree` thrown). The viewer
    remained attached, the canvas remained non-empty
    (`canvas.width > 0 && canvas.height > 0`), and
    `viewer.dataFrame === df` throughout.

    * Expected result: zero console errors across the full property
      sweep.

## Notes

- This file combines the smoke load and the full property sweep in one
  scenario (rather than splitting into a separate smoke file and a
  property-panel file) to preserve continuity — the smoke step seeds
  the viewer that the property sweep then mutates. It can be split
  later without changing what's covered.
- Selector pinning: per
  `.claude/skills/grok-browser/references/dendrogram.md`
  `## dendrogram-viewer-property-panel`, property-row anchors use a
  name-mangling rule — `nodeColumnName` mangles to `prop-node` (not
  `prop-node-column-name`), and `colorColumnName` mangles to
  `prop-color` (not `prop-color-column-name`). Tests that pin the
  un-mangled forms will not locate the row; see the ref doc's
  `## Selector validation matrix` for the full mapping before writing
  the Playwright spec.
- **Deferred.** Pixel-level canvas-output assertions (tree-layout
  geometry, exact leaf-row spacing under `step`, stroke-width clamping
  under `lineWidth`, label rendering under `showLabels`+`font`) are
  out of scope — pixel-diff infrastructure isn't available in the
  automated UI layer. The scenario instead asserts property-dispatch
  semantics (no console error, property reflected on the viewer
  instance, no-throw rebuild/restyle) without reaching for canvas
  pixel-equality.
