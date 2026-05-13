---
name: show-formula-lines
version: 0.1.0
description: |
  Overlay reference lines, fitted curves, and shaded acceptance bands on
  a Datagrok viewer by writing JSON into the dataframe's `.formula-lines`
  tag or the viewer's `formulaLines` property. For package developers
  who need persistent annotations on a scatter plot, line chart, density
  plot, box plot, histogram, or bar chart without subclassing the viewer.
  Use when asked to "overlay a fitted curve on a chart", "draw a y=x
  reference on a scatter plot", or "shade an acceptance band on a viewer".
triggers:
  - overlay a fitted curve
  - draw a reference line on a chart
  - shade an acceptance band
  - add a y=x diagonal
  - highlight a range on a plot
  - annotate a viewer with a curve
allowed-tools:
  - Read
  - Edit
  - Bash
harness-authored: true
---

# show-formula-lines

## When to use

You need to overlay non-data graphics on a chart — a reference line
(`y = 2 * x`), a fitted model curve, a band of acceptable range, an
identity diagonal — and have it persist with the dataframe or scoped
to a single viewer. Common phrasings: "overlay a 4PL fit on the
titration scatter plot", "draw a y=x diagonal on observed-vs-predicted",
"shade the in-spec range on the histogram".

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the
  package root, code lives under `src/`.
- `datagrok-api` imports — the article omits these:
  ```typescript
  import * as grok from 'datagrok-api/grok';
  import * as DG   from 'datagrok-api/dg';
  ```
- A `DG.DataFrame` (`grok.data.demo.demog()` works for demos) and a
  supported viewer: scatter plot, line chart, density plot, box plot,
  histogram, or bar chart (`DG-FACT-211`). Other viewers persist the
  JSON but silently drop it at render time.
- Familiarity with `manipulate-viewers` (sibling skill): the viewer
  is attached via `view.addViewer(DG.VIEWER.SCATTER_PLOT, opts)`,
  not the `view.scatterPlot(...)` wrapper which is documented as
  deprecated in its JSDoc (`js-api/src/views/view.ts:559-563`).

## Steps

1. **Pick the storage scope: dataframe vs. viewer.**
   Both expose the same `meta.formulaLines` helper (`DG-FACT-207`).
   Dataframe storage = JSON in tag `.formula-lines`
   (`DG.TAGS.FORMULA_LINES`, `js-api/src/const.ts:309`); renders on
   every viewer of that frame. Viewer storage =
   `viewer.props['formulaLines']`; renders only on that one viewer.
   ```typescript
   const df = grok.data.demo.demog();
   df.meta.formulaLines.addLine({                 // every viewer of df
     title: 'Identity', formula: '${height} = ${weight}',
   });
   ```

2. **Add a single line — `formula` is the only required field.**
   Left of `=` must be a single column; right side is any
   [Add New Column](../../../transform/add-new-column.md) expression
   over a second column or constant. Defaults: color `'#838383'`,
   `width: 1`, `style: 'solid'`, `zIndex: 100`, `opacity: 100`,
   `visible: true` (`DG-FACT-212`).
   ```typescript
   df.meta.formulaLines.addLine({
     title: 'Parabola',
     formula: '${height} = 180 + 0.01 * ${weight} * ${weight} - 1.5 * ${weight}',
     color: '#FFA500', width: 2, style: 'dashed', zIndex: -30,
   });
   ```
   Expected: line tooltip reads "Parabola"; canonical sample at
   `packages/ApiSamples/scripts/data-frame/metadata/formula-lines.js`.

3. **Add a band — both `formula` AND `column2` are required.**
   Band formula uses comparison/range syntax — `${col} < C`,
   `${col} > avg`, `${col} in(a, b)`, `${col} in(q1, q3)`. `column2`
   names the OTHER axis the band spans (`DG-FACT-209`). Default
   band color is `'#F0F0F0'`.
   ```typescript
   df.meta.formulaLines.addBand({
     formula: '${height} in(150, 180)',
     column2: 'weight',
     color: '#7FFFD4', opacity: 30,
   });
   ```

4. **Bulk-add via `addAll(items[])` — one re-serialization, not N.**
   Each `add*` call rewrites the entire JSON; for >2 items, batch
   through `addAll` (`DG-FACT-208`). Set `type: 'line' | 'band'`
   explicitly when using `add` / `addAll`.
   ```typescript
   df.meta.formulaLines.addAll([
     {type: 'line', formula: '${height} = 200', color: '#0000ff'},
     {type: 'band', formula: '${age} > 18', column2: 'sex'},
   ]);
   ```
   Real-world pattern: `packages/EDA/src/pls/pls-ml.ts:380`
   (`scatter.meta.formulaLines.addAll(getLines(names))`).

5. **Attach the viewer through `addViewer`, then add viewer-scoped lines.**
   Two viewer-level toggles control formula-line VISIBILITY without
   touching storage (`DG-FACT-210`): `showDataframeFormulaLines` and
   `showViewerFormulaLines`, both default `true`. Pass them at
   construction or via `setOptions(...)` — clearing storage is
   destructive; toggling these is reversible.
   ```typescript
   const view = grok.shell.addTableView(df);
   const plot = view.addViewer(DG.VIEWER.SCATTER_PLOT, {
     xColumnName: 'weight', yColumnName: 'height',
     showDataframeFormulaLines: true,
     showViewerFormulaLines: true,
   }) as DG.ScatterPlotViewer;

   plot.meta.formulaLines.addLine({                 // viewer-only line
     formula: '${weight} = 150', color: '#ff0000', width: 10,
   });
   ```
   Pattern: `packages/Chem/src/analysis/molecular-matched-pairs/mmp-viewer/mmp-viewer.ts:810-816`
   adds an Identity line on a paired grid's dataframe so every viewer
   of that frame renders it.

6. **Remove lines — `removeWhere`, NOT `removeAt`.**
   `removeAt(idx)` is broken: with the default `count=1` its slice
   `slice(idx, idx + count - 1) === slice(idx, idx)` returns `[]`,
   wiping the entire list (`js-api/src/helpers.ts:125-127`,
   `DG-FACT-208`). Use `removeWhere` until upstream is fixed.
   ```typescript
   df.meta.formulaLines.removeWhere((_, i) => i === 2);   // works
   df.meta.formulaLines.clear();                          // wipes all
   ```

## Common failure modes

- **Line added but never appears on the chart.** The viewer is in the
  unsupported set or you hit an axis it doesn't allow: box plot has
  no X-axis line (categorical), histogram has no Y-axis line (count
  is derived), bar chart accepts X only when `orientation:'horizontal'`,
  Y only when `'vertical'` (`DG-FACT-211`). The JS API persists the
  JSON either way — there is no `addLine`-time error.
- **Band silently missing — only `formula` was set.** `column2` is
  required for bands; without it the helper still serializes the
  item but rendering drops it. Always pass both (`DG-FACT-209`).
- **Opacity at `0.7` makes the line invisible.** `opacity` is on the
  `[0..100]` scale, NOT the CSS `[0..1]` convention; `0.7` rounds to
  nearly transparent. Use `70` (`DG-FACT-212`).
- **`removeAt(idx)` empties the whole list.** Upstream slice bug at
  `js-api/src/helpers.ts:125-127`; workaround is
  `removeWhere((_, i) => i === idx)` (`DG-FACT-208`).
- **Looped `addLine` for a large set causes UI hitch.** Each call
  re-stringifies the whole list (`DG-FACT-208`); collect and `addAll`.
- **Trellis Plot lines don't stick on the outer viewer.** Trellis is
  special-cased — its storage lives in `viewer.props['innerViewerLook']`,
  not on the outer viewer's own `formulaLines` (`DG-FACT-207`,
  `js-api/src/viewer.ts:838-851`).

## Verification

- `npm run build` (or `grok check`) exits `0`.
- `df.getTag('.formula-lines')` returns a JSON string whose parsed
  array length matches the count of `addLine`/`addBand` calls
  (`DG-FACT-207`).
- `df.meta.formulaLines.items.length` equals that same count, and
  every entry has the expected `type: 'line' | 'band'` and `formula`.
- Open the table view: lines/bands render on the scatter plot;
  `plot.setOptions({showDataframeFormulaLines: false})` hides them
  all without clearing storage; setting back to `true` restores them.
- For an unsupported viewer (e.g., bar chart with `orientation:
  'vertical'` and an X-axis formula), `items` still contains the
  entry but the chart shows nothing — confirms `DG-FACT-211`.

## See also

- Source articles:
  - `help/develop/how-to/viewers/show-formula-lines.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-207` through `DG-FACT-212`.
- Reference packages:
  - `packages/ApiSamples/scripts/data-frame/metadata/formula-lines.js` —
    canonical full-parameter example for both lines and bands.
  - `packages/EDA/src/pls/pls-ml.ts:371-383` — viewer constructed with
    `showViewerFormulaLines: true` followed by `addAll(...)` on the
    new viewer's `meta.formulaLines`.
  - `packages/Chem/src/analysis/molecular-matched-pairs/mmp-viewer/mmp-viewer.ts:810-816` —
    Identity line (`${Observed} = ${Predicted}`) added to a grid's
    dataframe so every paired viewer renders it.
- Related skills:
  - `manipulate-viewers` (sibling — defines the `addViewer` /
    `setOptions` surface this skill builds on).
