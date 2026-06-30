# Phase 7: QC Dashboard - Context

**Gathered:** 2026-03-06
**Status:** Ready for planning

<domain>
## Phase Boundary

Scientists can assess data quality through a unified dashboard of linked QC viewers (MA plot, missing value heatmap + bar chart, sample correlation, intensity distributions, CV plot) displayed in a tiled layout within the current TableView. The dashboard works at any pipeline stage (post-import, post-normalization, post-imputation) and requires group annotations.

</domain>

<decisions>
## Implementation Decisions

### Dashboard Entry & Layout
- Single menu item: "Proteomics | Visualize | QC Dashboard..." — nested under Visualize alongside existing viewers
- All QC viewers docked into current TableView using `tv.dockManager.dock()` in a tiled grid layout
- Requires group annotations (shows warning if not annotated, same pattern as PCA)
- Re-running QC replaces existing QC viewers (close old, create fresh) — no viewer accumulation

### Pipeline Stage
- QC dashboard works at any stage — always uses current log2() intensity columns
- No side-by-side pre/post comparison — shows current data state only
- Scientists can re-open QC after each pipeline step to compare visually

### MA Plot
- Uses annotated groups automatically: M = mean(group2) - mean(group1), A = 0.5 * (mean(g1) + mean(g2)) per protein
- Include loess/lowess trend line showing systematic intensity-dependent bias
- Color by up/down/not-significant direction if DE has been run (reuse `ensureDirectionColumn()` from volcano.ts); single color if no DE yet
- M=0 horizontal reference line

### Missing Value Display
- Two visualizations: binary heatmap (protein x sample) AND per-sample % missing bar chart
- Binary heatmap in separate DataFrame (0/1 matrix) — no cross-linking with other QC viewers (acceptable)
- Bar chart: one bar per sample, colored by experimental group

### Sample Correlation
- Pairwise Pearson correlation matrix displayed as heatmap
- Uses log2() intensity columns from all samples

### Intensity Distributions
- Per-sample box plots showing intensity distribution
- Requires unpivoted (long-format) data in separate DataFrame

### CV Plot
- CV = sd/mean per protein within each replicate group
- Scatter or density visualization

### Claude's Discretion
- Individual QC viewers as separate menu items vs dashboard-only (lean toward dashboard-only to keep menu clean)
- Exact tiled grid arrangement (2x3, 2x2+bottom, etc.) based on what works with dockManager
- CV plot visualization type (scatter vs density vs box plot per group)
- Sample correlation heatmap implementation details
- Loess implementation approach (TypeScript approximation vs moving average)
- Box plot data unpivot approach and whether it links to protein-level viewers

</decisions>

<specifics>
## Specific Ideas

- Loess trend line on MA plot is standard in proteomics QC tools (Perseus, limma plotMA) — scientists expect it
- Direction coloring on MA plot reuses existing volcano.ts logic — conditional on DE completion
- Missing value bar chart colored by group helps spot group-specific dropout patterns (MNAR)

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `getGroups()` in experiment-setup.ts: retrieves group assignments from DataFrame tag — needed for MA, CV, bar chart group coloring
- `ensureDirectionColumn()` in volcano.ts: computes up/down/not-significant — reuse for conditional MA plot coloring
- `findColumn()` / `findProteomicsColumns()` in column-detection.ts: locate intensity, protein ID columns
- `SEMTYPE` constants in proteomics-types.ts: filter intensity columns by semantic type
- Existing `showAllVisualizations()` pattern in package.ts: model for multi-viewer creation

### Established Patterns
- Viewers added via `tv.addViewer()` — auto-dock into current TableView
- Intensity columns identified by `semType === SEMTYPE.INTENSITY && name.startsWith('log2(')`
- Group data stored as JSON in `proteomics.groups` DataFrame tag
- PCA uses separate TableView for different-shape data (sample-level) — same pattern for box plots and missing value heatmap
- Column computation pattern: add computed columns to DataFrame, then create viewer referencing them (MA: M and A columns, CV: cv column)

### Integration Points
- New menu entry in `PackageFunctions` class alongside existing `showAllVisualizations()`
- New file: `src/viewers/qc-dashboard.ts` (or `src/viewers/qc/` directory if multiple files)
- QC computation utilities could go in `src/analysis/qc-metrics.ts`
- Imports `getGroups` from experiment-setup, `ensureDirectionColumn` from volcano, column detection utils

</code_context>

<deferred>
## Deferred Ideas

- Automated QC pass/fail scoring — explicitly out of scope per REQUIREMENTS.md
- QC report PDF export — explicitly out of scope per REQUIREMENTS.md
- Pre/post normalization side-by-side comparison — could be a future enhancement

</deferred>

---

*Phase: 07-qc-dashboard*
*Context gathered: 2026-03-06*
