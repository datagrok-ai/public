# Phase 7: QC Dashboard - Research

**Researched:** 2026-03-06
**Domain:** Proteomics QC visualization + Datagrok multi-viewer dashboard layout
**Confidence:** HIGH

## Summary

The QC dashboard requires computing proteomics-specific quality metrics (MA values, coefficient of variation, pairwise sample correlations, missing value patterns) and presenting them in a coordinated multi-viewer layout. The Datagrok platform provides all the necessary infrastructure: built-in viewers (scatter plot, box plot, correlation plot, grid/heatmap), a DockManager for programmatic layout, and automatic cross-viewer linking via shared DataFrame selection/filter/hover.

The core challenge is data preparation -- the QC metrics (MA values, CV, correlation matrix, missingness binary matrix) must be computed as new columns or separate DataFrames, then fed to standard Datagrok viewers. The MA plot and CV plot are standard scatter plots on computed columns. The missing values heatmap is a Grid in heatmap mode on a binary presence/absence matrix. The sample correlation matrix uses Datagrok's built-in correlation plot viewer. Per-sample intensity distributions use the built-in box plot viewer.

**Primary recommendation:** Compute QC metric columns on the existing protein-level DataFrame where possible (MA, CV). Use `tv.addViewer()` + `tv.dockManager.dock()` for layout. All viewers sharing one DataFrame get automatic cross-selection for free.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| QC-01 | User can open a QC dashboard showing all QC viewers in a linked layout | DockManager API for programmatic layout; `tv.addViewer()` + `tv.dockManager.dock()` with DockType (LEFT/RIGHT/TOP/DOWN) and ratio parameters |
| QC-02 | User can view an MA plot (M vs A) showing intensity-dependent bias between conditions | Scatter plot on computed `M` (log2 ratio) and `A` (mean average) columns; requires group annotations |
| QC-03 | User can view a missing values heatmap showing presence/absence pattern across samples | Grid in heatmap mode on binary 0/1 columns computed from intensity column nullness |
| QC-04 | User can view a sample correlation matrix as heatmap of pairwise Pearson correlations | Datagrok's built-in `CORR_PLOT` viewer with `xColumnNames`/`yColumnNames` set to intensity columns |
| QC-05 | User can view per-sample intensity distributions (box plots) | Unpivot intensity columns to long format, use Datagrok BOX_PLOT with category = sample name |
| QC-06 | User can view a CV plot showing measurement variability per protein within replicates | Scatter plot on computed `mean_intensity` vs `CV` columns per protein within each group |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| datagrok-api | (platform) | All viewers, DockManager, DataFrame, BitSet selection | Platform-provided; must not be bundled |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| (none) | - | - | All computation is pure TypeScript on DataFrame columns; no external libraries needed |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Built-in CORR_PLOT | Manual Pearson + custom heatmap | CORR_PLOT is built-in and handles rendering; no need to hand-roll |
| Built-in BOX_PLOT | Custom distribution viewer | BOX_PLOT is built-in with all needed features |
| Unpivot for box plot | Trellis plot of histograms | Box plot on unpivoted data is the standard proteomics approach |

**Installation:**
```bash
# No new dependencies needed -- all built-in Datagrok viewers
```

## Architecture Patterns

### Recommended Project Structure
```
src/
  viewers/
    qc-dashboard.ts      # Dashboard layout orchestration + menu entry
    qc-computations.ts   # MA, CV, missingness, correlation computations
  viewers/
    volcano.ts           # (existing)
    heatmap.ts           # (existing)
    pca-plot.ts          # (existing)
```

### Pattern 1: Dashboard Layout via DockManager
**What:** Create a TableView, add multiple viewers via `addViewer()`, then arrange them using `tv.dockManager.dock()` with DockType and ratio parameters.
**When to use:** Any time you need multiple viewers in a coordinated layout.
**Example:**
```typescript
// Source: js-api/src/docking.ts + packages/UITests/src/views/docking-nested.ts
const tv = grok.shell.addTableView(df);

// Add viewers -- they auto-link via shared DataFrame
const maPlot = tv.addViewer(DG.VIEWER.SCATTER_PLOT, {x: 'A', y: 'M'});
const boxPlot = tv.addViewer(DG.VIEWER.BOX_PLOT, {value: 'intensity', category1ColumnName: 'sample'});

// Arrange using dockManager
const rightNode = tv.dockManager.dock(maPlot, DG.DOCK_TYPE.RIGHT, null, 'MA Plot', 0.5);
tv.dockManager.dock(boxPlot, DG.DOCK_TYPE.DOWN, rightNode, 'Distributions', 0.5);
```

### Pattern 2: Computed Columns on Existing DataFrame
**What:** Add QC metric columns (M, A, CV, mean_intensity) directly to the protein-level DataFrame so all protein-row viewers share selection.
**When to use:** When the computed metric has one value per protein row (MA, CV).
**Example:**
```typescript
// MA plot: M = mean(group2) - mean(group1), A = (mean(group1) + mean(group2)) / 2
// Both are per-protein, so they belong as columns on the protein DataFrame
const mCol = df.columns.addNewFloat('M');
const aCol = df.columns.addNewFloat('A');
for (let i = 0; i < df.rowCount; i++) {
  const g1Mean = meanOfGroup(df, groups.group1.columns, i);
  const g2Mean = meanOfGroup(df, groups.group2.columns, i);
  mCol.set(i, g2Mean - g1Mean);
  aCol.set(i, (g1Mean + g2Mean) / 2);
}
```

### Pattern 3: Separate DataFrame for Different Row Granularity
**What:** When QC data has different row count than protein table (e.g., sample-level correlation matrix, unpivoted box plot data), create a separate DataFrame. Viewers on different DataFrames do NOT auto-link.
**When to use:** Box plot (unpivoted: rows = protein x sample), missing values heatmap (rows = proteins, but binary columns), correlation matrix (sample x sample).
**Example:**
```typescript
// Unpivot intensity columns for box plot
// Original: rows=proteins, cols=[log2(sample1), log2(sample2), ...]
// Unpivoted: rows=protein*sample, cols=[ProteinId, Sample, Intensity]
const longDf = unpivotIntensities(df, intensityCols);
const boxTv = grok.shell.addTableView(longDf);
boxTv.addViewer(DG.VIEWER.BOX_PLOT, {value: 'Intensity', category1ColumnName: 'Sample'});
```

### Pattern 4: Existing Viewer Pattern (follow volcano.ts)
**What:** Each viewer creation function returns a configured viewer. The dashboard function calls them and arranges the results.
**When to use:** All QC viewers.

### Anti-Patterns to Avoid
- **Creating custom JsViewers for standard plots:** Datagrok has scatter plot, box plot, correlation plot, and grid/heatmap built in. Use them.
- **Modifying the original DataFrame's intensity columns:** QC computations should ADD new columns, never modify existing log2 intensity columns.
- **Assuming cross-DataFrame selection:** Viewers on DIFFERENT DataFrames do NOT automatically share selection. Only viewers sharing the SAME DataFrame instance get free linking.
- **Adding too many columns to the main DataFrame:** If a QC metric requires many helper columns (e.g., one binary column per sample for missingness), consider a cloned DataFrame to avoid polluting the main one (same pattern as heatmap.ts).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Correlation matrix visualization | Custom heatmap renderer | `DG.VIEWER.CORR_PLOT` with `xColumnNames`/`yColumnNames` | Built-in viewer handles rendering, color scales, tooltips |
| Box plot for distributions | Custom distribution renderer | `DG.VIEWER.BOX_PLOT` with category column | Built-in handles quartiles, whiskers, outliers |
| Scatter plot for MA/CV | Custom canvas renderer | `DG.VIEWER.SCATTER_PLOT` with computed columns | Built-in handles zoom, selection, hover, labels |
| Multi-viewer layout | Custom HTML/CSS grid | `tv.dockManager.dock()` | DockManager handles resize, drag, split |
| Pearson correlation computation | External stats library | ~15 lines of TypeScript | Simple formula: sum of (x-mx)(y-my) / sqrt(sum(x-mx)^2 * sum(y-my)^2) |
| MA value computation | External biostats library | ~20 lines of TypeScript | M = mean(g2) - mean(g1), A = (mean(g1) + mean(g2))/2, per protein |
| CV computation | External stats library | ~15 lines of TypeScript | CV = sd / mean, per protein within group |

**Key insight:** All QC computations are straightforward per-row or per-column statistics. The Datagrok built-in viewers handle all the visualization complexity. The only real work is data preparation.

## Common Pitfalls

### Pitfall 1: Cross-DataFrame Selection Not Working
**What goes wrong:** User selects proteins in MA plot but box plot doesn't highlight.
**Why it happens:** Box plot uses unpivoted DataFrame (different row count), so selection doesn't propagate.
**How to avoid:** Keep MA plot and CV plot on the main protein DataFrame. For viewers that need different DataFrames (box plot, missing values), accept that cross-linking is limited or use event-based manual synchronization.
**Warning signs:** Selection in one viewer has no effect on another.

### Pitfall 2: Missing Group Annotations
**What goes wrong:** MA plot and CV plot require experimental group annotations, but user hasn't run "Annotate Experiment" yet.
**How to avoid:** Check `getGroups(df)` at dashboard entry; show clear warning message and return early if null. Same pattern as existing `showHeatmap()` and `showPcaPlot()`.
**Warning signs:** Null reference errors when accessing groups.

### Pitfall 3: All-Null Columns Breaking Statistics
**What goes wrong:** Protein with all-null intensities in a group causes NaN in mean/CV/MA calculations.
**Why it happens:** Division by zero when count is 0.
**How to avoid:** Guard every mean/sd/cv computation with a count > 0 check. Set result to DG.FLOAT_NULL when insufficient data.
**Warning signs:** NaN appearing in scatter plots, breaking axis scaling.

### Pitfall 4: Dock Layout Proportions
**What goes wrong:** Viewers are too small or overlap when docked.
**Why it happens:** The `ratio` parameter in `dock()` is relative to the reference node, and nesting order matters.
**How to avoid:** Plan the dock tree carefully. Typical 2x3 grid: split right first (0.5), then split each half vertically (0.5 each). Test with real data.
**Warning signs:** Viewers collapsed to tiny size, or one viewer taking up 90% of space.

### Pitfall 5: Performance with Many Proteins
**What goes wrong:** Computing pairwise sample correlations or unpivoting with 10,000+ proteins is slow.
**Why it happens:** Correlation matrix is O(n_samples^2 * n_proteins); unpivot creates n_proteins * n_samples rows.
**How to avoid:** Use typed arrays (Float32Array) for bulk computation. For unpivot, consider limiting to non-null values or sampling. Show a progress indicator.
**Warning signs:** UI freeze > 2 seconds on dashboard open.

### Pitfall 6: MA Plot Requires log2-Transformed Data
**What goes wrong:** MA plot shows nonsensical values because raw (non-log2) intensities were used.
**Why it happens:** M = difference of means only makes sense on log2 scale (where difference = log ratio).
**How to avoid:** Use `log2(...)` prefixed columns (same ones used throughout the pipeline). Check that columns start with `log2(` or have SEMTYPE.INTENSITY.
**Warning signs:** M values in thousands instead of -5 to +5 range.

## Code Examples

### Computing MA Values (per protein)
```typescript
// Source: Standard proteomics QC computation
// M = mean(group2) - mean(group1) on log2 intensity columns
// A = (mean(group1) + mean(group2)) / 2
function computeMA(df: DG.DataFrame, groups: GroupAssignment): void {
  const g1Cols = groups.group1.columns.map(n => df.col(n)!);
  const g2Cols = groups.group2.columns.map(n => df.col(n)!);
  const mCol = df.columns.addNewFloat('M');
  const aCol = df.columns.addNewFloat('A');

  for (let i = 0; i < df.rowCount; i++) {
    const g1Mean = groupMean(g1Cols, i);
    const g2Mean = groupMean(g2Cols, i);
    if (isNaN(g1Mean) || isNaN(g2Mean)) {
      mCol.set(i, DG.FLOAT_NULL);
      aCol.set(i, DG.FLOAT_NULL);
    } else {
      mCol.set(i, g2Mean - g1Mean);
      aCol.set(i, (g1Mean + g2Mean) / 2);
    }
  }
}

function groupMean(cols: DG.Column[], rowIdx: number): number {
  let sum = 0, count = 0;
  for (const col of cols) {
    if (!col.isNone(rowIdx)) {
      sum += col.get(rowIdx) as number;
      count++;
    }
  }
  return count > 0 ? sum / count : NaN;
}
```

### Computing CV (per protein within group)
```typescript
// CV = sd / mean for each protein within replicate group
function computeCV(df: DG.DataFrame, groupCols: string[], cvColName: string, meanColName: string): void {
  const cols = groupCols.map(n => df.col(n)!);
  const cvCol = df.columns.addNewFloat(cvColName);
  const meanCol = df.columns.addNewFloat(meanColName);

  for (let i = 0; i < df.rowCount; i++) {
    let sum = 0, sumSq = 0, count = 0;
    for (const col of cols) {
      if (!col.isNone(i)) {
        // CV on raw (non-log) scale: exponentiate log2 values
        const val = Math.pow(2, col.get(i) as number);
        sum += val;
        sumSq += val * val;
        count++;
      }
    }
    if (count < 2) {
      cvCol.set(i, DG.FLOAT_NULL);
      meanCol.set(i, count > 0 ? sum / count : DG.FLOAT_NULL);
    } else {
      const mean = sum / count;
      const variance = (sumSq - sum * sum / count) / (count - 1);
      const sd = Math.sqrt(variance);
      cvCol.set(i, mean > 0 ? sd / mean : DG.FLOAT_NULL);
      meanCol.set(i, mean);
    }
  }
}
```

### Missing Values Binary Matrix (clone approach)
```typescript
// Create a binary presence/absence DataFrame for missingness heatmap
function createMissingnessMatrix(df: DG.DataFrame, intensityCols: string[]): DG.DataFrame {
  const columns: DG.Column[] = [];

  // Label column (protein ID)
  const labelCol = findColumn(df, SEMTYPE.PROTEIN_ID, ['protein id', 'accession']);
  if (labelCol) columns.push(labelCol.clone());

  // Binary columns: 1 = present, 0 = missing
  for (const colName of intensityCols) {
    const srcCol = df.col(colName)!;
    const binCol = DG.Column.fromInt32Array(colName, new Int32Array(df.rowCount));
    for (let i = 0; i < df.rowCount; i++)
      binCol.set(i, srcCol.isNone(i) ? 0 : 1);
    columns.push(binCol);
  }

  return DG.DataFrame.fromColumns(columns);
}
```

### Dashboard Layout with DockManager
```typescript
// Source: js-api/src/docking.ts, packages/UITests/src/views/docking-nested.ts
function openQcDashboard(df: DG.DataFrame): void {
  const tv = grok.shell.tv!;  // or grok.shell.addTableView(df)

  // Compute QC columns (adds M, A, CV, meanIntensity to df)
  computeMA(df, groups);
  computeCV(df, groups.group1.columns, 'CV_group1', 'meanInt_group1');

  // Create viewers
  const maPlot = tv.addViewer(DG.VIEWER.SCATTER_PLOT, {x: 'A', y: 'M'});
  const cvPlot = tv.addViewer(DG.VIEWER.SCATTER_PLOT, {x: 'meanInt_group1', y: 'CV_group1'});
  const corrPlot = tv.addViewer(DG.VIEWER.CORR_PLOT, {
    xColumnNames: intensityCols,
    yColumnNames: intensityCols,
    showPearsonR: true,
  });

  // Layout: 2x2 grid
  // Grid occupies left half by default
  const rightNode = tv.dockManager.dock(maPlot, DG.DOCK_TYPE.RIGHT, null, 'MA Plot', 0.5);
  tv.dockManager.dock(cvPlot, DG.DOCK_TYPE.DOWN, rightNode, 'CV Plot', 0.5);
  // Dock correlation below the grid
  tv.dockManager.dock(corrPlot, DG.DOCK_TYPE.DOWN, null, 'Sample Correlation', 0.3);
}
```

### Getting Intensity Columns (reusable pattern)
```typescript
// Intensity columns are log2-prefixed and have SEMTYPE.INTENSITY
function getIntensityColumns(df: DG.DataFrame): string[] {
  return df.columns.toList()
    .filter(c => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('))
    .map(c => c.name);
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Standalone QC tools (Perseus, MSstats) | Integrated QC in analysis platforms | 2020+ | Scientists expect QC within the same tool |
| Static QC plots (PDF report) | Interactive linked viewers | 2018+ | Users want to click and explore, not just view |
| MA plot only for microarrays | MA plot standard for any -omics comparison | 2015+ | Applies to proteomics intensity data identically |

**Standard proteomics QC viewers (per Perseus/ProteomeDiscoverer/MSstats):**
- MA plot (aka Bland-Altman): intensity-dependent bias detection
- Missing values pattern: systematic vs random missingness
- Sample correlation heatmap: batch effects, outlier samples
- Per-sample distribution (box/violin): normalization assessment
- CV plot: technical reproducibility within replicates

## Open Questions

1. **Box plot data shape**
   - What we know: Datagrok BOX_PLOT expects a single value column + category column. Proteomics data has one column per sample.
   - What's unclear: Whether Datagrok has a built-in unpivot/melt function, or if we must manually create the long-format DataFrame.
   - Recommendation: Manually create long-format DataFrame (proteinId, sample, intensity) -- straightforward loop. This creates a SEPARATE DataFrame, so cross-selection with protein-level MA/CV plots will NOT work automatically. This is acceptable -- box plot shows sample-level distributions, not protein-level.

2. **Correlation plot column limit**
   - What we know: CORR_PLOT takes `xColumnNames` and `yColumnNames` arrays.
   - What's unclear: Performance with >50 samples (2,500+ cells in the matrix).
   - Recommendation: Proteomics experiments rarely exceed 50 samples. If they do, show a warning. Not a blocker.

3. **Missing values heatmap interactivity**
   - What we know: Grid in heatmap mode works (used in existing heatmap.ts). Binary 0/1 columns would render as two-color heatmap.
   - What's unclear: Whether binary Grid heatmap looks good visually with default color scale.
   - Recommendation: Use Grid heatmap mode with custom color scheme (white=missing, dark=present). If visual quality is poor, fall back to a simple bar chart showing % missing per sample.

## Validation Architecture

> `workflow.nyquist_validation` not set -- including validation section.

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Datagrok test framework (@datagrok-libraries/test) |
| Config file | webpack.config.js (test entry: package-test.ts) |
| Quick run command | `grok test --test "QC" --host localhost` |
| Full suite command | `grok test --host localhost` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| QC-01 | Dashboard opens with multiple viewers | integration | `grok test --test "QC Dashboard" --host localhost` | No -- Wave 0 |
| QC-02 | MA plot computes M and A correctly | unit | `grok test --test "MA computation" --host localhost` | No -- Wave 0 |
| QC-03 | Missing values binary matrix computed correctly | unit | `grok test --test "Missing values" --host localhost` | No -- Wave 0 |
| QC-04 | Correlation plot configured with correct columns | integration | `grok test --test "Correlation" --host localhost` | No -- Wave 0 |
| QC-05 | Box plot receives unpivoted data | unit | `grok test --test "Box plot" --host localhost` | No -- Wave 0 |
| QC-06 | CV computed correctly per group | unit | `grok test --test "CV computation" --host localhost` | No -- Wave 0 |

### Sampling Rate
- **Per task commit:** `grok test --test "QC" --host localhost`
- **Per wave merge:** `grok test --host localhost`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `src/tests/qc-dashboard.ts` -- covers QC-01 through QC-06
- [ ] Register test file in `src/package-test.ts`

## Sources

### Primary (HIGH confidence)
- `js-api/src/docking.ts` -- DockManager, DockNode, DockContainer API (read directly)
- `js-api/src/views/view.ts` -- TableView.addViewer(), TableView.dockManager (read directly)
- `js-api/src/const.ts` -- DOCK_TYPE enum (LEFT/RIGHT/TOP/DOWN/FILL), VIEWER enum (read directly)
- `js-api/src/interfaces/d4.ts` -- IBoxPlotSettings, ICorrelationPlotSettings, IScatterPlotSettings (read directly)
- `packages/UITests/src/views/docking-nested.ts` -- Working examples of nested dock layouts (read directly)
- `packages/Proteomics/src/viewers/` -- Existing viewer patterns: volcano.ts, heatmap.ts, pca-plot.ts (read directly)
- `packages/Proteomics/src/analysis/experiment-setup.ts` -- GroupAssignment interface, getGroups/setGroups (read directly)

### Secondary (MEDIUM confidence)
- Proteomics QC standards (MA plot, CV, missingness, correlation) -- well-established in Perseus, ProteomeDiscoverer, MSstats literature

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All Datagrok API verified by reading source
- Architecture: HIGH - DockManager pattern confirmed in UITests and multiple packages
- Pitfalls: HIGH - Based on understanding of DataFrame sharing model and existing code patterns
- QC computations: HIGH - Standard proteomics statistics (mean, sd, CV, Pearson correlation)

**Research date:** 2026-03-06
**Valid until:** 2026-04-06 (stable -- Datagrok API and proteomics statistics are well-established)
