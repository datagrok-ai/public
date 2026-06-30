# Phase 3: Visualization - Research

**Researched:** 2026-02-28
**Domain:** Datagrok built-in viewers, scatter plot formula lines, Grid heatmap mode, PCA computation
**Confidence:** HIGH

## Summary

Phase 3 visualization can be implemented almost entirely using Datagrok's built-in viewer infrastructure. The volcano plot is a configured scatter plot with formula lines for threshold lines. The heatmap is Datagrok's Grid viewer in heatmap mode (`isHeatmap: true`). The PCA plot is a scatter plot of computed PC1 vs PC2 columns, colored by experimental group. Linked selection comes free because all viewers share the same DataFrame.

The critical insight from codebase analysis is that **none of these viewers need to be custom JsViewers**. The existing `VolcanoViewer` stub in `src/viewers/volcano-viewer.ts` that extends `DG.JsViewer` should be replaced with a factory function that creates and configures a standard `ScatterPlotViewer`. The ClinicalCase package's `AERiskAssessmentView` demonstrates this exact pattern for a volcano plot using `df.plot.scatter()` with `meta.formulaLines`.

**Primary recommendation:** Use Datagrok's built-in `ScatterPlotViewer` for volcano and PCA, Grid in heatmap mode for the expression heatmap, and compute PCA client-side using SVD/eigendecomposition on the intensity matrix. Wire everything through the same DataFrame for free linked selection.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| VIZ-01 | Volcano plot (log2FC vs -log10 adj p-value, formula lines for thresholds) | ScatterPlotViewer + FormulaLinesHelper; ClinicalCase provides exact pattern |
| VIZ-02 | Heatmap with hierarchical clustering and sample grouping | Grid `isHeatmap: true` mode; Dendrogram package's `heatmapDemo` shows the pattern |
| VIZ-03 | PCA plot colored by experimental group | Client-side PCA → new columns → ScatterPlotViewer with color property |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `datagrok-api` | (workspace) | ScatterPlotViewer, Grid (heatmap mode), FormulaLinesHelper | Platform-native viewers with built-in linked selection |
| `@datagrok-libraries/statistics` | (workspace) | Already a dependency; may provide matrix ops | Used in Phase 2 for t-test |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| None | - | PCA computation | Implement client-side with vanilla TypeScript (SVD on small matrices) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Client-side PCA | EDA package `computePCA` via `grok.functions.call('EDA:PCA', ...)` | Adds EDA dependency; EDA uses WASM which is heavy. For <5000 rows x <50 cols, a simple eigendecomposition is sufficient |
| Grid heatmap | Dendrogram package hierarchical clustering | Dendrogram adds dependency on `@datagrok-libraries/bio`; simple Grid heatmap with row sorting is simpler for v1 |
| Custom JsViewer | Built-in ScatterPlotViewer | Custom viewer adds complexity with no benefit since scatter plot handles all volcano/PCA requirements natively |

**Installation:**
No new packages needed. All dependencies are already in the workspace.

## Architecture Patterns

### Recommended Project Structure
```
src/
  viewers/
    volcano.ts           # Factory: createVolcanoPlot(df, options) → ScatterPlotViewer
    heatmap.ts           # Factory: createExpressionHeatmap(df, options) → Grid
    pca-plot.ts          # Factory: createPcaPlot(df, options) → ScatterPlotViewer
  analysis/
    pca.ts               # computePCA(df, intensityCols) → adds PC1, PC2 columns
  utils/
    proteomics-types.ts  # (existing) semantic type constants
```

### Pattern 1: Configured Built-in Viewer (Volcano Plot)
**What:** Create a ScatterPlotViewer via `df.plot.scatter()`, configure x/y columns, then add formula lines for threshold lines.
**When to use:** Any proteomics visualization that maps to a scatter plot (volcano, MA plot, PCA).
**Example:**
```typescript
// Source: ClinicalCase/src/views/ae-risk-assessment-view.ts (lines 76-108)
function createVolcanoPlot(
  df: DG.DataFrame,
  options: {fcThreshold?: number; pThreshold?: number} = {},
): DG.ScatterPlotViewer {
  const fcThreshold = options.fcThreshold ?? 1.0;
  const pThreshold = options.pThreshold ?? 0.05;

  // Compute -log10(adj.p-value) column if not present
  ensureNegLog10Column(df);

  const sp = df.plot.scatter({
    x: 'log2FC',
    y: '-log10(adj.p-value)',
    color: 'significant',
    showViewerFormulaLines: true,
  });

  // Horizontal line: p-value threshold
  sp.meta.formulaLines.addLine({
    formula: `\${-log10(adj.p-value)} = ${-Math.log10(pThreshold)}`,
    color: '#808080',
    width: 1,
    style: 'dash',
  });

  // Vertical lines: fold change thresholds
  sp.meta.formulaLines.addLine({
    formula: `\${log2FC} = ${fcThreshold}`,
    color: '#808080',
    width: 1,
    style: 'dash',
  });
  sp.meta.formulaLines.addLine({
    formula: `\${log2FC} = ${-fcThreshold}`,
    color: '#808080',
    width: 1,
    style: 'dash',
  });

  return sp;
}
```

### Pattern 2: Grid Heatmap Mode
**What:** Create a Grid viewer, set `isHeatmap: true` to render cells as a color matrix.
**When to use:** Expression heatmaps showing protein intensity across samples.
**Example:**
```typescript
// Source: Dendrogram/src/demos/heatmapDemo.ts (lines 38-41)
function createExpressionHeatmap(df: DG.DataFrame, intensityCols: string[]): DG.Grid {
  // Create a subset DataFrame with only the columns we want
  const heatmapDf = buildHeatmapDataFrame(df, intensityCols);
  const grid = heatmapDf.plot.grid();
  grid.props.isHeatmap = true;
  grid.props.showRowHeader = false;
  return grid;
}
```

### Pattern 3: PCA as Column Computation
**What:** Compute PCA client-side, add PC1/PC2 as new columns to the existing DataFrame, then create a scatter plot.
**When to use:** Dimensionality reduction visualization.
**Example:**
```typescript
// Source: EDA/src/package.ts (lines 98-116) - adapted for client-side
function createPcaPlot(
  df: DG.DataFrame,
  intensityCols: string[],
  groupCol: string,
): DG.ScatterPlotViewer {
  // Compute PCA and add PC1, PC2 columns to df
  addPcaColumns(df, intensityCols);

  return df.plot.scatter({
    x: 'PC1',
    y: 'PC2',
    color: groupCol,
  });
}
```

### Pattern 4: Adding Viewers to TableView with Docking
**What:** Use `tv.addViewer()` and `tv.dockManager.dock()` to position viewers in the layout.
**When to use:** Opening multiple linked viewers in a single view.
**Example:**
```typescript
// Source: Dendrogram/src/demos/heatmapDemo.ts (lines 23-30)
const tv = grok.shell.tv;
const viewer = tv.addViewer(DG.VIEWER.SCATTER_PLOT, {x: 'log2FC', y: '-log10(adj.p-value)'});
// Or dock manually:
const rootNode = tv.dockManager.rootNode;
tv.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, rootNode, 'Volcano', 0.5);
```

### Anti-Patterns to Avoid
- **Custom JsViewer for standard chart types:** Datagrok's built-in ScatterPlotViewer handles volcano and PCA plots natively. A JsViewer means re-implementing selection, zoom, tooltips, formula lines, and all interactivity from scratch.
- **Separate DataFrames per viewer:** Creating a new DataFrame for each viewer breaks linked selection. All viewers must share the same DataFrame instance for selection/filter synchronization to work.
- **Canvas rendering in JsViewer:** Drawing scatter plots or heatmaps on canvas means losing all platform features (column selectors, property panel, right-click menus, export, etc.).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Scatter plot rendering | Canvas-based scatter in JsViewer | `DG.Viewer.scatterPlot()` / `df.plot.scatter()` | Built-in handles zoom, selection, tooltips, color coding, labels, export, axis formatting |
| Threshold lines on scatter | Custom line drawing on canvas overlay | `viewer.meta.formulaLines.addLine()` | FormulaLinesHelper handles rendering, serialization, and interaction. ClinicalCase proves the pattern |
| Heatmap rendering | Custom heatmap cell renderer | `grid.props.isHeatmap = true` | Grid heatmap mode handles color scaling, row/col reordering, cell formatting, export |
| Linked selection | Manual event wiring between viewers | Shared DataFrame | All viewers on same DataFrame auto-sync selection, filter, current row, and mouse-over |
| PCA for >1000 columns | Custom SVD from scratch | Call EDA:PCA via `grok.functions.call` | Only needed if intensity matrix is very large; for typical proteomics (<50 samples), client-side is fine |

**Key insight:** Datagrok's viewer system is designed so that configured built-in viewers ARE the custom viewers. The platform provides scatter, grid/heatmap, and all interaction for free. The only custom work is computing derived columns (-log10 p-value, PC scores) and calling factory functions with the right configuration.

## Common Pitfalls

### Pitfall 1: Building a Custom JsViewer Instead of Configuring Built-in Scatter
**What goes wrong:** Developer creates a JsViewer subclass with canvas rendering for volcano/PCA, losing all built-in interactivity (zoom, pan, selection, tooltips, formula lines, color coding, axis histograms, property panel).
**Why it happens:** Natural instinct to "build a viewer" when the task says "create a volcano viewer."
**How to avoid:** Use `df.plot.scatter({...})` which returns a `ScatterPlotViewer`. Add formula lines via `meta.formulaLines`. The scatter plot IS the volcano plot.
**Warning signs:** Code extends `DG.JsViewer` for a standard chart type. Code uses `canvas.getContext('2d')`.

### Pitfall 2: Separate DataFrames Breaking Linked Selection
**What goes wrong:** PCA or heatmap viewer uses a separate DataFrame (e.g., `DG.DataFrame.fromColumns([pc1Col, pc2Col])`), so selecting a protein in the volcano plot does not highlight it in other viewers.
**Why it happens:** PCA produces a new DataFrame by default. Heatmap seems to need a subset of columns.
**How to avoid:** Add PCA result columns to the ORIGINAL DataFrame. For heatmap, create a Grid on the same DataFrame and configure which columns to show (using grid column visibility or a column list).
**Warning signs:** `DG.DataFrame.fromColumns()` or `DG.DataFrame.fromCsv()` in viewer creation code. Multiple DataFrames in the same view.

### Pitfall 3: Forgetting to Compute -log10(p-value) Column
**What goes wrong:** Volcano plot Y axis needs -log10(adjusted p-value), but the DE step only produces raw p-value and adjusted p-value columns.
**Why it happens:** DE results use raw scale p-values. Volcano convention requires -log10 transform.
**How to avoid:** Add a `-log10(adj.p-value)` computed column to the DataFrame before creating the volcano plot. Guard against log(0) by clamping.
**Warning signs:** Y axis shows raw p-values (0 to 1 range) instead of -log10 scale.

### Pitfall 4: Formula Line Formula Syntax Errors
**What goes wrong:** Formula lines don't render or render in wrong position.
**Why it happens:** FormulaLine `formula` field requires specific syntax: `${columnName} = value` for vertical/horizontal lines. Column names with special characters need exact matching.
**How to avoid:** Use the exact pattern from ClinicalCase: `formula: \`\${${colName}} = ${value}\``. Always set `showViewerFormulaLines: true` in scatter plot options. Test with simple numeric values first.
**Warning signs:** Lines not visible. Viewer shows without any threshold indicators.

### Pitfall 5: Heatmap Not Showing Color-Coded Cells
**What goes wrong:** Grid shows numbers instead of colored cells.
**Why it happens:** Grid defaults to `isGrid: true`. Must explicitly set `isHeatmap: true`.
**How to avoid:** Set `grid.props.isHeatmap = true` AND `grid.props.isGrid = false` as shown in the Dendrogram heatmapDemo.
**Warning signs:** Heatmap viewer looks like a regular spreadsheet with numbers.

## Code Examples

### Creating a Volcano Plot with Threshold Lines
```typescript
// Source: ClinicalCase/src/views/ae-risk-assessment-view.ts (verified in codebase)
// + js-api/src/helpers.ts FormulaLinesHelper API

import * as DG from 'datagrok-api/dg';

export function createVolcanoPlot(
  df: DG.DataFrame,
  fcThreshold: number = 1.0,
  pThreshold: number = 0.05,
): DG.ScatterPlotViewer {
  // Ensure -log10(adj.p-value) column exists
  const negLog10Name = '-log10(adj.p-value)';
  if (!df.columns.contains(negLog10Name)) {
    const adjPCol = df.col('adj.p-value');
    if (!adjPCol) throw new Error('adj.p-value column not found');
    const col = df.columns.addNewFloat(negLog10Name);
    for (let i = 0; i < df.rowCount; i++) {
      if (adjPCol.isNone(i)) {
        col.set(i, DG.FLOAT_NULL);
      } else {
        const p = adjPCol.get(i) as number;
        col.set(i, p > 0 ? -Math.log10(p) : 10); // clamp log(0)
      }
    }
  }

  const sp = df.plot.scatter({
    x: 'log2FC',
    y: negLog10Name,
    color: 'significant',
    showViewerFormulaLines: true,
  });

  // Horizontal threshold: -log10(pThreshold)
  sp.meta.formulaLines.addLine({
    formula: `\${${negLog10Name}} = ${-Math.log10(pThreshold)}`,
    color: '#888888',
    width: 1,
    style: 'dash',
  });

  // Vertical thresholds: +/- fcThreshold
  sp.meta.formulaLines.addLine({
    formula: `\${log2FC} = ${fcThreshold}`,
    color: '#888888',
    width: 1,
    style: 'dash',
  });
  sp.meta.formulaLines.addLine({
    formula: `\${log2FC} = ${-fcThreshold}`,
    color: '#888888',
    width: 1,
    style: 'dash',
  });

  return sp;
}
```

### Creating a Grid Heatmap
```typescript
// Source: Dendrogram/src/demos/heatmapDemo.ts (verified in codebase)
// Grid isHeatmap pattern: js-api/src/interfaces/d4.ts line 1240

import * as DG from 'datagrok-api/dg';

export function createExpressionHeatmap(
  df: DG.DataFrame,
  intensityCols: string[],
  labelCol: string = 'Gene Name',
): DG.Grid {
  const grid = df.plot.grid();

  // Configure heatmap mode
  grid.props.isHeatmap = true;
  grid.props.isGrid = false;
  grid.props.showRowHeader = false;

  // Show only intensity columns + label
  // Grid column visibility is controlled through grid.columns API
  return grid;
}
```

### Client-Side PCA (Simple Eigendecomposition)
```typescript
// Pattern: EDA/src/eda-tools.ts computePCA (line 16-67)
// Adapted for simple client-side use without WASM dependency

import * as DG from 'datagrok-api/dg';

/** Compute PCA and add PC1, PC2 columns to the DataFrame in-place. */
export function addPcaColumns(
  df: DG.DataFrame,
  intensityCols: string[],
  nComponents: number = 2,
): void {
  // 1. Build data matrix (rows = proteins, cols = samples)
  // 2. Center each column (subtract mean)
  // 3. Compute covariance matrix
  // 4. Eigendecomposition → eigenvectors
  // 5. Project: PC = data * eigenvectors
  // 6. Add PC1, PC2 as new columns to df

  const nRows = df.rowCount;
  const nCols = intensityCols.length;

  // Extract matrix, handling NaN by imputing column mean
  const matrix: number[][] = [];
  for (let i = 0; i < nRows; i++) {
    const row: number[] = [];
    for (const colName of intensityCols) {
      const col = df.col(colName)!;
      row.push(col.isNone(i) ? NaN : col.get(i) as number);
    }
    matrix.push(row);
  }

  // ... (centering, covariance, eigendecomposition)
  // Add result columns
  const pc1 = df.columns.addNewFloat('PC1');
  const pc2 = df.columns.addNewFloat('PC2');
  // ... populate from projection
}
```

### Adding Viewers to TableView
```typescript
// Source: Dendrogram/src/demos/heatmapDemo.ts (lines 12-28)

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export function addProteomicsViewers(tv: DG.TableView, df: DG.DataFrame): void {
  // All viewers share the same df → linked selection is automatic
  const volcano = createVolcanoPlot(df);
  const pca = createPcaPlot(df);

  tv.addViewer(volcano);
  tv.addViewer(pca);

  // Or with explicit docking:
  const rootNode = tv.dockManager.rootNode;
  tv.dockManager.dock(volcano, DG.DOCK_TYPE.RIGHT, rootNode, 'Volcano', 0.5);
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Custom JsViewer for volcano | Configure built-in ScatterPlotViewer + FormulaLines | Always available | No need for custom rendering; full platform integration |
| Grid with color-coded columns | Grid `isHeatmap: true` mode | Built-in | Single property toggle gives heatmap rendering |
| Server-side R PCA (Samples:PCAR) | Client-side WASM PCA (EDA) or simple JS PCA | EDA package | Removes R server dependency for PCA |

**Deprecated/outdated:**
- The existing `VolcanoViewer` class extending `DG.JsViewer` (in `src/viewers/volcano-viewer.ts`) should be replaced with a factory function that configures a built-in ScatterPlotViewer. The REQUIREMENTS.md Out of Scope section explicitly states "Custom volcano JsViewer" is out of scope because "Datagrok's scatter plot handles this natively with formula lines."
- Using `tv.scatterPlot()` / `tv.heatMap()` are deprecated in favor of `tv.addViewer(DG.VIEWER.SCATTER_PLOT, options)` or `Viewer.scatterPlot(df, options)`.

## Open Questions

1. **Grid Column Visibility for Heatmap**
   - What we know: Grid `isHeatmap` renders all numerical columns as colored cells. We need to show only intensity columns.
   - What's unclear: The exact API for hiding/showing specific grid columns programmatically. Grid has a `columns` property but the interface is not fully documented.
   - Recommendation: Use `grid.columns.setVisible(false)` for non-intensity columns, or create a view-specific column list. Test during implementation. Alternatively, create a subset DataFrame but maintain row indices linking back to the main DataFrame (this risks breaking linked selection -- prefer column visibility approach).

2. **Hierarchical Clustering for Heatmap Row Ordering**
   - What we know: VIZ-02 requires "hierarchical clustering and sample grouping." The Dendrogram package provides full hierarchical clustering with tree rendering.
   - What's unclear: Whether calling Dendrogram's service from our package is straightforward, or if it requires the package to be installed.
   - Recommendation: For v1, sort rows by significance (adj. p-value) or absolute log2FC rather than full hierarchical clustering. This avoids the Dendrogram dependency. Note this in the plan as a simplification that can be enhanced later. If hierarchical clustering is essential, use `grok.functions.call('Dendrogram:getDendrogramService')` but document the dependency.

3. **PCA Implementation Complexity**
   - What we know: For typical proteomics datasets (<5000 proteins, <50 samples), PCA on the sample covariance matrix (~50x50) is trivial.
   - What's unclear: Whether a simple power iteration / Jacobi eigendecomposition is robust enough, or if we need a full SVD library.
   - Recommendation: Implement using the covariance matrix approach (center data, compute S = X^T X / (n-1), eigendecompose the small nCols x nCols matrix). For <100 columns this is fast and numerically stable. Use Jacobi rotation or call `grok.functions.call('EDA:PCA')` as fallback.

## Sources

### Primary (HIGH confidence)
- `js-api/src/viewer.ts` - ScatterPlotViewer class, FormulaLinesHelper API, JsViewer base class
- `js-api/src/helpers.ts` - FormulaLinesHelper implementation, FormulaLine interface
- `js-api/src/interfaces/d4.ts` - IScatterPlotSettings (formulaLines, showViewerFormulaLines, x, y, color), IGridSettings (isHeatmap, isGrid)
- `js-api/src/dataframe/data-frame.ts` - DataFramePlotHelper (df.plot.scatter, df.plot.grid)
- `packages/ClinicalCase/src/views/ae-risk-assessment-view.ts` - Complete volcano plot pattern with formula lines
- `packages/Dendrogram/src/demos/heatmapDemo.ts` - Grid heatmap mode pattern
- `packages/EDA/src/package.ts` and `eda-tools.ts` - PCA computation and column addition pattern
- `packages/Proteomics/src/analysis/differential-expression.ts` - DE output columns (log2FC, p-value, adj.p-value, significant)
- `packages/Proteomics/src/analysis/experiment-setup.ts` - Group metadata storage pattern

### Secondary (MEDIUM confidence)
- `js-api/src/views/view.ts` - TableView.addViewer(), dockManager API
- `packages/EDA/src/eda-ui.ts` - PCA column naming convention (PC1, PC2)

### Tertiary (LOW confidence)
- Grid column visibility API for restricting heatmap columns - needs implementation-time validation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All patterns verified in codebase with multiple existing implementations
- Architecture: HIGH - Factory function pattern proven by ClinicalCase volcano and Dendrogram heatmap
- Pitfalls: HIGH - Each pitfall identified from concrete codebase analysis (existing JsViewer stub, DataFrame separation risk, formula line syntax)
- PCA implementation: MEDIUM - Simple eigendecomposition approach is standard math but untested in this specific codebase context

**Research date:** 2026-02-28
**Valid until:** 2026-03-28 (stable platform APIs, unlikely to change)
