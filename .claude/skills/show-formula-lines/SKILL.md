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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# show-formula-lines

## When to use

You need to overlay non-data graphics on a chart — a reference line
(`y = 2 * x`), a fitted model curve, a band of acceptable range, an
identity diagonal — and have it persist with the dataframe or scoped
to a single viewer. Common phrasings: "overlay a 4PL fit on the
titration scatter plot", "draw a y=x diagonal on observed-vs-predicted",
"shade the in-spec range on the histogram".

## Prerequisites

- A `DG.DataFrame` (`grok.data.demo.demog()` works for demos) and a
  supported viewer: scatter plot, line chart, density plot, box plot,
  histogram, or bar chart (`DG-FACT-211`). Other viewers persist the
  JSON but silently drop it at render time.
- Familiarity with `manipulate-viewers` (sibling skill): the viewer
  is attached via `view.addViewer(DG.VIEWER.SCATTER_PLOT, opts)`,
  not the `view.scatterPlot(...)` wrapper which is documented as
  deprecated in its JSDoc (`js-api/src/views/view.ts:559-563`).

## Steps

1. **Pick the storage scope: dataframe vs. viewer.** Both expose the
   same `meta.formulaLines` helper (`DG-FACT-207`). Dataframe storage
   renders on every viewer of that frame; viewer storage renders only
   on that one viewer.
   ```typescript
   const df = grok.data.demo.demog();
   df.meta.formulaLines.addLine({                 // every viewer of df
     title: 'Identity', formula: '${height} = ${weight}',
   });
   ```

2. **Add a single line — `formula` is the only required field.**
   Left of `=` must be a single column; right side is any *Add New
   Column* expression. Defaults at `DG-FACT-212` (color `#838383`,
   `width: 1`, `style: 'solid'`, `zIndex: 100`, `opacity: 100`).
   ```typescript
   df.meta.formulaLines.addLine({
     title: 'Parabola',
     formula: '${height} = 180 + 0.01 * ${weight} * ${weight} - 1.5 * ${weight}',
     color: '#FFA500', width: 2, style: 'dashed', zIndex: -30,
   });
   ```

3. **Add a band — both `formula` AND `column2` are required**
   (`DG-FACT-209`). Band formula uses comparison/range syntax:
   `${col} < C`, `${col} in(a, b)`.
   ```typescript
   df.meta.formulaLines.addBand({
     formula: '${height} in(150, 180)',
     column2: 'weight',
     color: '#7FFFD4', opacity: 30,
   });
   ```

4. **Bulk-add via `addAll(items[])` — one re-serialization, not N**
   (`DG-FACT-208`). Set `type: 'line'|'band'` explicitly when using
   `add`/`addAll`.
   ```typescript
   df.meta.formulaLines.addAll([
     {type: 'line', formula: '${height} = 200', color: '#0000ff'},
     {type: 'band', formula: '${age} > 18', column2: 'sex'},
   ]);
   ```

5. **Attach the viewer through `addViewer`, then add viewer-scoped lines.**
   `showDataframeFormulaLines`/`showViewerFormulaLines` (both default
   `true`) toggle visibility without destroying storage (`DG-FACT-210`).
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

6. **Remove lines — `removeWhere`, NOT `removeAt`.** `removeAt` is
   broken — wipes the entire list (`DG-FACT-208`).
   ```typescript
   df.meta.formulaLines.removeWhere((_, i) => i === 2);   // works
   df.meta.formulaLines.clear();                          // wipes all
   ```

## Common failure modes

- Line added but never appears — unsupported viewer/axis combo; JS API doesn't validate (`DG-FACT-211`).
- Band silently missing — `column2` required for bands (`DG-FACT-209`).
- `opacity: 0.7` near-invisible — scale is `[0..100]`, not `[0..1]`; use `70` (`DG-FACT-212`).
- `removeAt(idx)` empties list — upstream bug; use `removeWhere((_, i) => i === idx)` (`DG-FACT-208`).
- Looped `addLine` causes UI hitch — re-serializes per call; use `addAll` (`DG-FACT-208`).
- Trellis Plot lines don't stick — storage lives in `innerViewerLook`, not outer viewer (`DG-FACT-207`).

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
