---
name: show-formula-lines
description: Draw reference lines and bands on a Datagrok viewer via dataframe or viewer formula-line storage
---

# show-formula-lines

## When to use

You need to overlay non-data graphics on a chart — a reference line
(`y = 2 * x`), a band of acceptable range, an identity diagonal — and
have it persist with the dataframe or scoped to one viewer. Triggers:
"draw a y=x line on the scatter plot", "shade an acceptance band",
"highlight a region of the histogram".

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
  histogram, or bar chart (`DG-FACT-211`). Other viewers silently drop
  formula-line JSON at render time.
- Familiarity with `manipulate-viewers` (sibling skill) — formula-line
  attach uses `view.addViewer(...)`, not the deprecated
  `view.scatterPlot(...)` wrapper (`DG-FACT-DRIFT-083`).

## Steps

1. **Pick the storage scope: dataframe vs. viewer.**
   Both expose the same `meta.formulaLines` helper (`DG-FACT-207`).
   Dataframe storage = JSON in tag `.formula-lines`
   (`DG.TAGS.FORMULA_LINES`), renders on every viewer of that frame.
   Viewer storage = JSON in `viewer.props['formulaLines']`, renders
   only on that one viewer:
   ```typescript
   const df = grok.data.demo.demog();
   df.meta.formulaLines.addLine({                 // every viewer of df
     title: 'Identity', formula: '${height} = ${weight}',
   });
   ```

2. **Add a single line — `formula` is the only required field.**
   Left of `=` must be a single column; right side is any
   [Add New Column](../../../transform/add-new-column.md) expression
   over a second column or constant. Defaults are line color
   `#838383`, `width: 1`, `style: 'solid'`, `zIndex: 100`,
   `opacity: 100`, `visible: true` (`DG-FACT-212`):
   ```typescript
   df.meta.formulaLines.addLine({
     title: 'Parabola',
     formula: '${height} = 180 + 0.01 * ${weight} * ${weight} - 1.5 * ${weight}',
     color: '#FFA500', width: 2, style: 'dashed', zIndex: -30,
   });
   ```
   Expected: line tooltip reads "Parabola"; sample at
   `packages/ApiSamples/scripts/data-frame/metadata/formula-lines.js:13-35`.

3. **Add a band — both `formula` AND `column2` are required.**
   Band formula uses the comparison/range syntax
   `${col} < C`, `${col} > avg`, `${col} in(a, b)`, `${col} in(q1, q3)`.
   `column2` names the OTHER axis the band spans (`DG-FACT-209`).
   Default band color is `#F0F0F0`:
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
   explicitly when going through `add` / `addAll`:
   ```typescript
   df.meta.formulaLines.addAll([
     {type: 'line', formula: '${height} = 200', color: '#0000ff'},
     {type: 'band', formula: '${age} > 18', column2: 'sex'},
   ]);
   ```
   Pattern: `packages/EDA/src/pls/pls-ml.ts:380`
   (`scatter.meta.formulaLines.addAll(getLines(names))`).

5. **Attach the viewer through `addViewer`, not `view.scatterPlot`.**
   Two viewer-level toggles control formula-line VISIBILITY without
   touching storage (`DG-FACT-210`): `showDataframeFormulaLines`
   and `showViewerFormulaLines`, both default `true`. The article's
   `view.scatterPlot(...)` example uses a deprecated wrapper
   (`DG-FACT-DRIFT-083`) — use `addViewer` instead:
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
   adds an Identity line to a grid's dataframe before the scatter renders.

6. **Remove lines — `removeWhere`, NOT `removeAt`.**
   `removeAt(idx)` is broken upstream — its slice math collapses the
   list to `[]` for the documented `count = 1` (`DG-FACT-DRIFT-081`).
   Use `removeWhere` until fixed; `clear()` empties the storage:
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
  `[0..100]` scale, NOT CSS `[0..1]`; `0.7` rounds to "nearly
  transparent". Use `70` (`DG-FACT-212`).
- **`removeAt(2)` empties the whole list.** Upstream slice bug;
  workaround is `removeWhere((_, i) => i === idx)`
  (`DG-FACT-DRIFT-081`).
- **Lines were `view.scatterPlot(...)` style and TypeScript flags
  deprecation.** `View.scatterPlot` carries `@deprecated`
  (`DG-FACT-DRIFT-083`); rewrite as
  `view.addViewer(DG.VIEWER.SCATTER_PLOT, opts)`.
- **Built large set with looped `addLine` and the page hitches.**
  Each call re-stringifies the whole list (`DG-FACT-208`); collect
  items first and call `addAll(items)` once.

## Verification

- TypeScript build (`npm run build` or `grok check`) exits `0` and
  emits no `deprecated` warnings on viewer-attach calls.
- `df.getTag('.formula-lines')` returns a JSON string whose parsed
  array length matches the count of `addLine`/`addBand` calls
  (`DG-FACT-207`).
- `df.meta.formulaLines.items.length` equals that same count, and
  every entry has the expected `type: 'line' | 'band'` and `formula`.
- Open the table view: lines/bands render on the scatter plot;
  toggling `viewer.setOptions({showDataframeFormulaLines: false})`
  hides them all without clearing storage; setting it back to `true`
  restores them.
- For an unsupported viewer (e.g., bar chart with `orientation:
  'vertical'` and an X-axis formula), `items` still contains the
  entry but the chart shows nothing — confirms behavior matches
  `DG-FACT-211`.

## See also

- Source articles:
  - `help/develop/how-to/viewers/show-formula-lines.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-207` … `DG-FACT-212` and drifts `DG-FACT-DRIFT-081`,
  `DG-FACT-DRIFT-082`, `DG-FACT-DRIFT-083`.
- Reference packages:
  - `packages/ApiSamples/scripts/data-frame/metadata/formula-lines.js` —
    canonical full-parameter example for both lines and bands.
  - `packages/EDA/src/pls/pls-ml.ts:371-383` — `DG.Viewer.scatterPlot`
    with `showViewerFormulaLines: true` followed by `addAll(...)` on
    the new viewer's `meta.formulaLines`.
  - `packages/Chem/src/analysis/molecular-matched-pairs/mmp-viewer/mmp-viewer.ts:810-816` —
    Identity line (`${Observed} = ${Predicted}`) added to a grid's
    dataframe so every paired viewer renders it.
- Related skills:
  - `manipulate-viewers` (sibling — defines the `addViewer` /
    `setOptions` surface this skill builds on; covers the
    `DG-FACT-DRIFT-080`/`-083` deprecation drift).
