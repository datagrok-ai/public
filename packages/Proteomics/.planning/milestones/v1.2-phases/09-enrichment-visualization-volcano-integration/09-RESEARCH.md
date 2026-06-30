# Phase 9: Enrichment Visualization & Volcano Integration - Research

**Researched:** 2026-03-06
**Domain:** Datagrok viewer configuration, cross-DataFrame event wiring, enrichment data visualization
**Confidence:** HIGH

## Summary

This phase adds three capabilities: (1) a dot plot of enriched terms using Datagrok's scatter plot viewer with size and color mapping, (2) a bar chart of enriched terms using Datagrok's native bar chart viewer, and (3) cross-DataFrame selection wiring that highlights member proteins on the volcano plot when a user selects an enrichment term.

The key insight from researching the Datagrok API is that the scatter plot viewer has `sizeColumnName` and `colorColumnName` properties, which means a dot plot can be built entirely from the existing scatter plot viewer -- no custom viewer is needed. The enrichment DataFrame from Phase 8 already has the required columns (`Gene Count`, `Gene Ratio`, `FDR`). The cross-DataFrame linkage follows an established pattern seen in Bio and other packages: subscribe to `onCurrentRowChanged` on the enrichment DataFrame and programmatically set `selection` on the protein DataFrame using gene membership lookup.

**Primary recommendation:** Build the dot plot as a configured `DG.Viewer.scatterPlot()` on a top-N filtered clone of the enrichment DataFrame. Build the bar chart as a `DG.Viewer.barChart()` on the same clone. Wire cross-DataFrame selection via `enrichmentDf.onCurrentRowChanged` -> parse Intersection column -> set `proteinDf.selection` on matching gene rows. No custom viewers or external libraries needed.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| VIZ-02 | User can view a dot plot of top enriched terms (gene ratio, gene count, adjusted p-value) | Scatter plot viewer supports categorical Y axis (Term Name), numerical X axis (Gene Ratio), `sizeColumnName` (Gene Count), `colorColumnName` (FDR). Top-N achieved by cloning enrichment DataFrame with BitSet filter (same pattern as heatmap top-N). |
| VIZ-03 | User can view a bar chart of top enriched terms ranked by significance or gene count | Bar chart viewer supports `splitColumnName` (Term Name), `valueColumnName` (Gene Count or -log10FDR), `barSortType`/`barSortOrder` for ranking. |
| ENRICH-04 | User can select a GO/pathway term in enrichment results and see its member proteins highlighted on the volcano plot | Cross-DataFrame pattern: `enrichmentDf.onCurrentRowChanged` -> read Intersection column -> find matching rows in protein DataFrame by gene symbol -> set `proteinDf.selection`. Established pattern from Bio package (sequence-similarity-viewer.ts line 88-93). |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| datagrok-api | workspace version | ScatterPlotViewer (dot plot), BarChartViewer, DataFrame events, BitSet selection | Platform API; all viewer capabilities needed are built-in |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Existing `enrichment.ts` exports | Phase 8 | `buildEnrichmentDf()`, `GostResult` type, enrichment DataFrame schema | Source of enrichment data; schema is locked |
| Existing `volcano.ts` | Phase 3 | `createVolcanoPlot()`, volcano viewer reference | Target for protein highlighting |
| Existing `column-detection.ts` | Phase 1 | `findColumn()`, `findProteomicsColumns()` | Finding gene symbol column in protein DataFrame |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| ScatterPlot as dot plot | Custom JsViewer canvas rendering | Massive effort; ScatterPlot already has size/color mapping, tooltips, selection |
| DataFrame clone for top-N | DataFrame filter (BitSet) on original | Filter would affect enrichment table grid; clone isolates views (same pattern as heatmap) |
| `onCurrentRowChanged` for cross-DF | Shared global event bus | No standard global bus in Datagrok; direct DataFrame events are the established pattern |

**No additional npm packages needed.** This phase uses only existing Datagrok API viewers and DataFrame events.

## Architecture Patterns

### Recommended Project Structure
```
src/
  viewers/
    enrichment-viewers.ts    # NEW: dot plot, bar chart, and cross-DF wiring
  analysis/
    enrichment.ts            # EXISTING: enrichment pipeline (Phase 8) -- NO changes needed
  viewers/
    volcano.ts               # EXISTING: volcano plot -- NO changes needed
  package.ts                 # MODIFIED: add menu entries for enrichment viz
```

### Pattern 1: Dot Plot via Configured ScatterPlot
**What:** Create a scatter plot with categorical Y axis (term names), numerical X axis (gene ratio), size mapped to gene count, color mapped to -log10(FDR).
**When to use:** VIZ-02 -- enrichment dot plot visualization
**Example:**
```typescript
// Source: Datagrok js-api IScatterPlotSettings interface (js-api/src/interfaces/d4.ts)
function createEnrichmentDotPlot(enrichDf: DG.DataFrame, topN: number = 20): DG.ScatterPlotViewer {
  // 1. Create top-N clone sorted by FDR ascending
  const fdrCol = enrichDf.col('FDR')!;
  const indexed: {idx: number; fdr: number}[] = [];
  for (let i = 0; i < enrichDf.rowCount; i++) {
    if (!fdrCol.isNone(i))
      indexed.push({idx: i, fdr: fdrCol.get(i) as number});
  }
  indexed.sort((a, b) => a.fdr - b.fdr);

  const mask = DG.BitSet.create(enrichDf.rowCount);
  for (let i = 0; i < Math.min(topN, indexed.length); i++)
    mask.set(indexed[i].idx, true);
  const topDf = enrichDf.clone(mask);

  // 2. Add -log10(FDR) column for color mapping
  const negLogCol = topDf.columns.addNewFloat('-log10(FDR)');
  const topFdr = topDf.col('FDR')!;
  for (let i = 0; i < topDf.rowCount; i++) {
    const fdr = topFdr.get(i) as number;
    negLogCol.set(i, fdr > 0 ? -Math.log10(fdr) : 10);
  }

  // 3. Create scatter plot with size and color mapping
  const sp = DG.Viewer.scatterPlot(topDf, {
    x: 'Gene Ratio',
    y: 'Term Name',                // categorical Y axis
    sizeColumnName: 'Gene Count',  // bubble size = gene count
    colorColumnName: '-log10(FDR)', // color = significance
    markerMinSize: 5,
    markerMaxSize: 25,
    showLabelsFor: 'None' as any,
    title: 'Enrichment Dot Plot',
  });

  return sp;
}
```

### Pattern 2: Bar Chart of Enrichment Terms
**What:** Create a horizontal bar chart showing top terms ranked by -log10(FDR) or gene count.
**When to use:** VIZ-03 -- enrichment bar chart visualization
**Example:**
```typescript
// Source: Datagrok js-api IBarChartSettings interface (js-api/src/interfaces/d4.ts)
function createEnrichmentBarChart(
  topDf: DG.DataFrame,  // reuse same top-N clone from dot plot
): DG.Viewer {
  const bar = DG.Viewer.barChart(topDf, {
    splitColumnName: 'Term Name',
    valueColumnName: '-log10(FDR)',
    valueAggrType: 'avg',
    barSortType: 'by value',
    barSortOrder: 'desc',
    orientation: 'Horizontal',
    title: 'Top Enriched Terms',
  } as any);

  return bar;
}
```

### Pattern 3: Cross-DataFrame Selection Wiring (Enrichment -> Volcano)
**What:** When user clicks a row in enrichment table/dot plot, highlight member proteins on the volcano plot by setting selection on the protein DataFrame.
**When to use:** ENRICH-04 -- linking enrichment terms to volcano plot
**Example:**
```typescript
// Source: Bio package pattern (sequence-similarity-viewer.ts:88-93)
// and Datagrok DataFrame API (data-frame.ts:161, bit-set.ts)
function wireEnrichmentToVolcano(
  enrichDf: DG.DataFrame,
  proteinDf: DG.DataFrame,
): rxjs.Subscription {
  // Build gene-to-row index for the protein DataFrame
  const geneCol = findColumn(proteinDf, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']);
  if (!geneCol) return rxjs.Subscription.EMPTY;

  const geneToRows = new Map<string, number[]>();
  for (let i = 0; i < proteinDf.rowCount; i++) {
    if (!geneCol.isNone(i)) {
      const gene = (geneCol.get(i) as string).trim();
      if (!geneToRows.has(gene)) geneToRows.set(gene, []);
      geneToRows.get(gene)!.push(i);
    }
  }

  // Subscribe to enrichment row changes
  return enrichDf.onCurrentRowChanged.subscribe(() => {
    const rowIdx = enrichDf.currentRowIdx;
    if (rowIdx < 0) return;

    const intersectionCol = enrichDf.col('Intersection');
    if (!intersectionCol) return;

    const memberGenesStr = intersectionCol.get(rowIdx) as string;
    if (!memberGenesStr) return;

    const memberGenes = memberGenesStr.split(',').map((g) => g.trim()).filter(Boolean);

    // Set selection on protein DataFrame
    proteinDf.selection.setAll(false, false);
    for (const gene of memberGenes) {
      const rows = geneToRows.get(gene);
      if (rows) {
        for (const row of rows)
          proteinDf.selection.set(row, true, false);
      }
    }
    proteinDf.selection.fireChanged();
  });
}
```

### Pattern 4: Top-N Filtered Clone (Shared by Dot Plot and Bar Chart)
**What:** Clone enrichment DataFrame with only top-N rows by FDR, isolating viewer state from the enrichment table grid.
**When to use:** Both dot plot and bar chart should use the same filtered clone.
**Example:**
```typescript
// Source: heatmap.ts pattern (lines 47-63) -- existing in codebase
function createTopNEnrichmentDf(enrichDf: DG.DataFrame, topN: number = 20): DG.DataFrame {
  const fdrCol = enrichDf.col('FDR')!;
  const indexed: {idx: number; fdr: number}[] = [];
  for (let i = 0; i < enrichDf.rowCount; i++) {
    if (!fdrCol.isNone(i))
      indexed.push({idx: i, fdr: fdrCol.get(i) as number});
  }
  indexed.sort((a, b) => a.fdr - b.fdr);

  const mask = DG.BitSet.create(enrichDf.rowCount);
  for (let i = 0; i < Math.min(topN, indexed.length); i++)
    mask.set(indexed[i].idx, true);

  const topDf = enrichDf.clone(mask);
  topDf.name = `Top ${Math.min(topN, indexed.length)} Enriched Terms`;

  // Add -log10(FDR) for visualization
  const negLogCol = topDf.columns.addNewFloat('-log10(FDR)');
  const topFdr = topDf.col('FDR')!;
  for (let i = 0; i < topDf.rowCount; i++) {
    const fdr = topFdr.get(i) as number;
    negLogCol.set(i, fdr > 0 ? -Math.log10(fdr) : 10);
  }

  return topDf;
}
```

### Pattern 5: Enrichment Visualization Dashboard
**What:** Open dot plot + bar chart in a docked layout alongside the enrichment table, with cross-DF wiring to the protein table.
**When to use:** Main entry point from menu or auto-opened after enrichment analysis.
**Example:**
```typescript
// Source: qc-dashboard.ts docking pattern (lines 136-142)
function openEnrichmentVisualization(
  enrichDf: DG.DataFrame,
  proteinDf: DG.DataFrame,
): void {
  const topDf = createTopNEnrichmentDf(enrichDf, 20);

  // Find or create table view for enrichment results
  const tv = grok.shell.getTableView(enrichDf.name)
    ?? grok.shell.addTableView(enrichDf);

  const dotPlot = createEnrichmentDotPlot(topDf);
  const barChart = createEnrichmentBarChart(topDf);

  // Dock viewers
  const dotNode = tv.dockManager.dock(dotPlot, DG.DOCK_TYPE.RIGHT, null, 'Dot Plot', 0.5);
  tv.dockManager.dock(barChart, DG.DOCK_TYPE.DOWN, dotNode, 'Bar Chart', 0.5);

  // Wire cross-DF selection (enrichment -> volcano/protein table)
  wireEnrichmentToVolcano(enrichDf, proteinDf);
}
```

### Anti-Patterns to Avoid
- **Building a custom canvas viewer for the dot plot:** The Datagrok scatter plot already supports categorical axes, size mapping, and color mapping. A custom viewer would be thousands of lines of canvas code for no benefit.
- **Modifying the enrichment DataFrame for visualization:** Adding visualization-only columns (like -log10(FDR)) to the original enrichment DataFrame would pollute the data. Use a cloned DataFrame.
- **Using DataFrame.filter for top-N in shared views:** Setting `filter` on the enrichment DataFrame would hide rows from the grid/table view too. Use `clone()` for isolation (established heatmap pattern).
- **Polling for changes across DataFrames:** Use `onCurrentRowChanged` subscription, not polling. Datagrok's event system is reactive.
- **Forgetting to call `fireChanged()` after programmatic selection:** BitSet operations with `notify: false` (third arg) require an explicit `fireChanged()` call to update the UI.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Dot plot visualization | Custom JsViewer with canvas drawing | `DG.Viewer.scatterPlot()` with `sizeColumnName` + `colorColumnName` | Scatter plot already renders categorical Y axis with size/color mapping; gets free tooltips, selection, zoom |
| Bar chart visualization | Custom bar rendering | `DG.Viewer.barChart()` with sort options | Native viewer with built-in sort, selection highlighting, formatting |
| Cross-DataFrame linking | Global event bus or shared state | `DataFrame.onCurrentRowChanged` + `selection.set()` | Established pattern in Bio/Peptides packages; direct, debuggable |
| Top-N filtering | DataFrame filter manipulation | `DataFrame.clone(BitSet)` | Isolates viewer from grid state; same pattern as heatmap (viewers/heatmap.ts) |
| Gene membership index | Linear scan per term selection | `Map<string, number[]>` built once | O(1) lookup per gene symbol instead of O(n) per selection event |

**Key insight:** Every visualization in this phase can be built from existing Datagrok viewers with configuration. No custom rendering, no external charting libraries, no npm dependencies.

## Common Pitfalls

### Pitfall 1: Scatter Plot with Categorical Axis Does Not Show Expected Layout
**What goes wrong:** Term names are too long and overlap, or the scatter plot treats them as numerical values and shows a blank plot.
**Why it happens:** Datagrok's scatter plot supports categorical axes, but very long term names may clip or overlap.
**How to avoid:** Truncate Term Name to ~50 characters when building the top-N clone. Set adequate `yAxisWidth` if needed via `sp.props.autoAxisSize = false` and `sp.props.yAxisWidth = 200`. Consider using `Term Name` only for the top-N clone where there are at most 20 terms.
**Warning signs:** Labels overlapping or plot area squeezed to a thin strip.

### Pitfall 2: Cross-DF Selection Not Updating Volcano
**What goes wrong:** Clicking a term in the enrichment table does nothing to the volcano plot.
**Why it happens:** (a) The gene symbol column is not found because the column was added by Phase 8 mapping but has a different name than expected. (b) The subscription was created but the enrichment DataFrame reference is stale (e.g., enrichment was re-run).
**How to avoid:** Use `findColumn()` with semantic type `SEMTYPE.GENE_SYMBOL` first, then name hints. Store the subscription and clean it up on re-run. Always check that `enrichDf.currentRowIdx >= 0` before reading column values.
**Warning signs:** No visual change on the volcano plot; check browser console for errors.

### Pitfall 3: Intersection Column Parse Failure
**What goes wrong:** Gene lookup fails because gene symbols in the Intersection column have extra whitespace or unexpected formatting.
**Why it happens:** `buildEnrichmentDf()` uses `memberGenes.join(', ')` -- the space after comma must be trimmed when parsing back.
**How to avoid:** Always `split(',').map(s => s.trim()).filter(Boolean)` when reading the Intersection column.
**Warning signs:** Selection highlights zero rows despite valid member genes.

### Pitfall 4: fireChanged() Not Called After Batch Selection
**What goes wrong:** Selection bits are set correctly but the volcano plot does not visually update.
**Why it happens:** `BitSet.set(idx, true, false)` suppresses notification to avoid firing per-bit. But `fireChanged()` must be called afterward to propagate the change.
**How to avoid:** Always call `proteinDf.selection.fireChanged()` after the selection loop completes.
**Warning signs:** Selection count shows 0 in UI even though `selection.trueCount` is non-zero when checked programmatically.

### Pitfall 5: Top-N Clone Loses Tag/Column Metadata
**What goes wrong:** The cloned DataFrame loses the `proteomics.enrichment` tag or column semantic types.
**Why it happens:** `DataFrame.clone()` copies tags and column metadata by default, but derived columns added after cloning (like `-log10(FDR)`) need their own metadata.
**How to avoid:** Verify that the tag `proteomics.enrichment` persists on the clone. If adding derived columns, set their tags explicitly.
**Warning signs:** Other code checking for `getTag('proteomics.enrichment')` does not recognize the clone.

### Pitfall 6: Memory Leak from Unmanaged Subscriptions
**What goes wrong:** Each time enrichment visualization is opened, a new `onCurrentRowChanged` subscription is created without cleaning up the old one.
**Why it happens:** No subscription management -- the old subscription stays active.
**How to avoid:** Return the subscription from the wiring function. Store it in a module-level variable and call `.unsubscribe()` before creating a new one. Alternatively, use `rxjs.take(1)` if the wiring is one-shot per session.
**Warning signs:** Multiple selection events firing for a single click; stale DataFrame references.

## Code Examples

### Enrichment DataFrame Schema (from Phase 8 -- locked)
```typescript
// Source: analysis/enrichment.ts buildEnrichmentDf()
// Columns available in enrichment DataFrame:
//   Source: string       (GO:BP, GO:MF, GO:CC, KEGG, REAC)
//   Term ID: string      (GO:0006915, KEGG:04210, etc.)
//   Term Name: string    (apoptotic process, Apoptosis, etc.)
//   P-value: float       (adjusted p-value from g:GOSt)
//   FDR: float           (same as P-value for FDR method)
//   Gene Count: int      (intersection_size from API)
//   Gene Ratio: float    (precision = intersection_size / query_size)
//   Intersection: string (comma-separated gene symbols: "TP53, BRCA1, EGFR")
//   Significant: bool    (FDR < 0.05)
// Tag: proteomics.enrichment = "true"
```

### ScatterPlot Settings for Dot Plot (Verified API)
```typescript
// Source: js-api/src/interfaces/d4.ts IScatterPlotSettings
// Key properties used:
//   x / xColumnName: string         - X axis column
//   y / yColumnName: string         - Y axis column (supports categorical)
//   sizeColumnName: string          - numerical column for marker size
//   colorColumnName: string         - numerical/categorical column for color
//   markerMinSize: number           - minimum marker size in pixels
//   markerMaxSize: number           - maximum marker size in pixels
//   markerOpacity: number           - 0.0 to 1.0
//   showSizeSelector: boolean       - show size column picker
//   showColorSelector: boolean      - show color column picker
//   title: string                   - viewer title
//   invertColorScheme: boolean      - invert the color ramp
```

### BarChart Settings (Verified API)
```typescript
// Source: js-api/src/interfaces/d4.ts IBarChartSettings
// Key properties used:
//   splitColumnName: string         - categorical column (each bar = one category)
//   valueColumnName: string         - value column (what determines bar height)
//   valueAggrType: string           - aggregation: 'avg', 'sum', 'count', etc.
//   barSortType: string             - 'by value' or 'by category'
//   barSortOrder: string            - 'asc' or 'desc'
//   orientation: string             - 'Horizontal' for horizontal bars
//   colorColumnName: string         - optional color coding
//   title: string                   - viewer title
```

### BitSet Selection Pattern (Verified API)
```typescript
// Source: js-api/src/dataframe/bit-set.ts
// BitSet methods used for selection:
//   selection.setAll(false, false)   - clear all; false = no notify
//   selection.set(idx, true, false)  - set bit; false = no notify
//   selection.fireChanged()          - fire event after batch update
//   selection.trueCount              - number of selected rows
```

### DataFrame Clone Pattern (Verified -- used in heatmap.ts)
```typescript
// Source: js-api/src/dataframe/data-frame.ts line 288
// DataFrame.clone(rowMask?: BitSet, columnIds?: string[], saveSelection?: boolean): DataFrame
// Creates a new DataFrame with specified rows. Tags are preserved.
const mask = DG.BitSet.create(df.rowCount);
mask.set(0, true);
mask.set(5, true);
const cloned = df.clone(mask);  // only rows 0 and 5
```

### DockManager Pattern (Verified -- used in qc-dashboard.ts)
```typescript
// Source: viewers/qc-dashboard.ts lines 136-142
const dm = tv.dockManager;
const nodeA = dm.dock(viewerA, DG.DOCK_TYPE.RIGHT, null, 'Title A', 0.5);
dm.dock(viewerB, DG.DOCK_TYPE.DOWN, nodeA, 'Title B', 0.5);
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| External charting (D3, Plotly) | Datagrok built-in viewers | Always (for Datagrok packages) | Native viewers integrate with selection, filtering, tooltip system |
| Custom dot plot viewer | ScatterPlot with size/color mapping | Datagrok supports categorical axes | No custom rendering needed |
| Manual DOM event wiring | DataFrame.onCurrentRowChanged | Datagrok reactive events | Cleaner, fewer leaks |

**GO redundancy reduction (REVIGO-like):** Out of scope per requirements. The requirements ask for "compact top-N view" which is handled by simply sorting by FDR and taking top 20 terms. Full semantic similarity clustering (rrvgo/REVIGO) would require a GO DAG and semantic similarity matrix -- massive scope that is explicitly out of scope per REQUIREMENTS.md ("GO DAG visualization: Large DAG rendering slow; rarely useful in practice").

## Open Questions

1. **Exact behavior of scatter plot with many categories on Y axis**
   - What we know: ScatterPlot supports categorical columns on both axes. With 20 terms, each gets a horizontal band.
   - What's unclear: Whether very long term names (>60 chars, common in GO) will cause layout issues.
   - Recommendation: Truncate term names to ~50 chars when building the top-N clone. Test with real data during implementation. If layout is poor, fall back to shorter term name + tooltip.

2. **Whether to auto-open enrichment viz after enrichment analysis or add separate menu entry**
   - What we know: DE auto-opens volcano plot (package.ts line 105). PCA opens in separate table view.
   - What's unclear: Whether user prefers auto-open or on-demand.
   - Recommendation: Add a menu entry `Proteomics | Visualize | Enrichment Plots...` AND auto-open after enrichment analysis completes (consistent with DE -> volcano pattern). Both share the same function.

3. **Whether cross-DF selection should also work with cloned top-N DataFrame**
   - What we know: The dot plot and bar chart operate on a cloned top-N DataFrame. The enrichment table grid shows the full DataFrame.
   - What's unclear: Whether clicking a dot in the dot plot (which is on the clone) should also highlight proteins.
   - Recommendation: Wire `onCurrentRowChanged` on BOTH the full enrichment DataFrame AND the top-N clone to the protein DataFrame. The user may interact with either view.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | @datagrok-libraries/test (Puppeteer-based) |
| Config file | packages/Proteomics/src/package-test.ts |
| Quick run command | `grok test --host localhost --category "Enrichment Visualization"` |
| Full suite command | `grok test --host localhost` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| VIZ-02 | createEnrichmentDotPlot returns a ScatterPlotViewer with correct axis/size/color mappings | unit | `grok test --host localhost --test "dot plot"` | No -- Wave 0 |
| VIZ-03 | createEnrichmentBarChart returns a BarChart viewer with correct split/value/sort | unit | `grok test --host localhost --test "bar chart"` | No -- Wave 0 |
| ENRICH-04 | wireEnrichmentToVolcano sets selection on protein DataFrame when enrichment row changes | unit | `grok test --host localhost --test "cross-DF selection"` | No -- Wave 0 |
| (shared) | createTopNEnrichmentDf filters to correct number of rows sorted by FDR | unit | `grok test --host localhost --test "top-N filter"` | No -- Wave 0 |

### Sampling Rate
- **Per task commit:** `grok test --host localhost --category "Enrichment Visualization"`
- **Per wave merge:** `grok test --host localhost`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `src/tests/enrichment-visualization.ts` -- new test file covering VIZ-02, VIZ-03, ENRICH-04, top-N filtering
- [ ] Add `import './tests/enrichment-visualization';` to `package-test.ts`
- [ ] No new framework install needed -- @datagrok-libraries/test already configured

### Test Strategy Notes
- **Unit tests (no network):** Test `createTopNEnrichmentDf()` with mock enrichment DataFrame (verify row count, sort order, -log10(FDR) column). Test `wireEnrichmentToVolcano()` by creating mock enrichment and protein DataFrames, setting currentRowIdx, and verifying selection bits. Test viewer creation functions return non-null viewers with expected property values.
- **Manual verification:** Visual layout of dot plot (bubble sizes, color gradient), bar chart (sort order, orientation), cross-DF highlighting on actual volcano plot. These are UI-level and cannot be fully automated.

## Sources

### Primary (HIGH confidence)
- Datagrok js-api `IScatterPlotSettings` interface (js-api/src/interfaces/d4.ts lines 2745-3117) -- all scatter plot properties including `sizeColumnName`, `colorColumnName`, categorical axis support
- Datagrok js-api `IBarChartSettings` interface (js-api/src/interfaces/d4.ts lines 4-173) -- bar chart split, value, sort properties
- Datagrok js-api `DataFrame` class (js-api/src/dataframe/data-frame.ts) -- `clone()`, `selection`, `onCurrentRowChanged`, `currentRowIdx`
- Datagrok js-api `BitSet` class (js-api/src/dataframe/bit-set.ts) -- `set()`, `setAll()`, `fireChanged()`, `create()`
- Existing Proteomics package code -- `enrichment.ts` (schema), `volcano.ts` (viewer), `heatmap.ts` (top-N clone pattern), `qc-dashboard.ts` (dock pattern)

### Secondary (MEDIUM confidence)
- Bio package `sequence-similarity-viewer.ts` (lines 88-93) -- cross-DataFrame row change subscription pattern
- Peptides package `mutation-cliffs-viewer.ts` (lines 172-220) -- additional cross-DF selection pattern confirmation

### Tertiary (LOW confidence)
- Scatter plot categorical axis behavior with long labels -- needs runtime verification with real enrichment data

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all APIs verified from source code in js-api
- Architecture: HIGH -- all patterns verified from existing Proteomics codebase (heatmap.ts, qc-dashboard.ts, enrichment.ts)
- Pitfalls: HIGH -- identified from actual API signatures and established patterns in Bio/Peptides packages
- Cross-DF wiring: HIGH -- verified from multiple package implementations in codebase
- Categorical scatter plot: MEDIUM -- API supports it (categorical column type on Y axis) but visual behavior with long GO term names needs runtime testing

**Research date:** 2026-03-06
**Valid until:** 2026-04-06 (30 days -- Datagrok API is stable)
