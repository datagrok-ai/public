# Phase 3: Visualization - Context

**Gathered:** 2026-02-28
**Status:** Ready for planning

<domain>
## Phase Boundary

Interactive proteomics viewers for differential expression results: volcano plot, heatmap, and PCA plot. Volcano and heatmap share linked selection via a common protein-level DataFrame. PCA is a sample-level visualization on a separate DataFrame (nSamples rows vs nProteins rows), so selection is not linked between PCA and the protein-level viewers -- this is architecturally correct since PCA shows samples, not proteins. This phase does NOT add new analysis capabilities — it visualizes results from Phase 2.

</domain>

<decisions>
## Implementation Decisions

### Volcano Plot
- Adjustable significance thresholds (FC and p-value cutoffs) via viewer property panel — users can explore different cutoffs without re-running DE
- Three-color point scheme: red = significantly up, blue = significantly down, gray = not significant (standard proteomics convention)
- Auto-label top N most significant proteins by gene name; user can toggle labels on/off
- Click selects the protein row in the DataFrame — linked viewers update automatically (standard Datagrok behavior)

### Heatmap
- Default to top N significant DE proteins (e.g., top 50 by adjusted p-value); user can change N
- Hierarchical clustering on rows (proteins) with dendrogram
- Blue-white-red diverging color scale (z-score normalized across rows)
- Columns grouped by experimental group (control left, treatment right) using Phase 2 group annotations

### PCA Plot
- Axis labels show variance explained: "PC1 (45.2%)" format
- Label all sample points with sample names (typical proteomics datasets have <50 samples)
- 95% confidence ellipses drawn around each experimental group
- Default to PC1 vs PC2; user can switch to other components via property panel
- Compute all components upfront

### Menu & Workflow
- Menu items under Proteomics | Visualize | {Volcano Plot, Heatmap, PCA}
- "Show All Visualizations" menu item opens all three in a docked/split layout
- Volcano and heatmap check for DE prerequisite; show helpful message if DE hasn't been run
- PCA available anytime after import/normalization (useful as QC step before DE)
- Auto-open volcano plot after DE completes, plus keep menu items for reopening

### Claude's Discretion
- Exact number of top proteins to label on volcano plot
- Default N for heatmap protein count
- PCA computation approach (eigendecomposition method)
- Exact docked layout arrangement for "Show All" mode
- Loading/progress indicators during PCA computation
- Confidence ellipse visual style (fill opacity, border)

</decisions>

<specifics>
## Specific Ideas

- Volcano plot should follow Perseus/EnhancedVolcano conventions familiar to proteomics scientists
- Heatmap should look publication-ready with standard diverging color scale
- PCA with confidence ellipses is standard in proteomics QC and publications

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `VolcanoViewer` JsViewer stub exists at `src/viewers/volcano-viewer.ts` — research recommends replacing with configured ScatterPlotViewer factory function
- `findColumn()` utility in `src/utils/column-detection.ts` — detects columns by semantic type, useful for auto-detecting log2FC, p-value columns
- `SEMTYPE` constants in `src/utils/proteomics-types.ts` — LOG2FC, P_VALUE, PROTEIN_ID semantic types already defined
- `getGroups()` from `src/analysis/experiment-setup.ts` — retrieves group annotations for heatmap column grouping and PCA coloring

### Established Patterns
- Viewer registration via `@grok.decorators.func()` with `tags: ['viewer']` — see existing Volcano registration in package.ts
- Menu items use `//tags: menuItem` and `//top-menu:` metadata comments
- Prerequisite checks pattern: check DataFrame tag (e.g., `df.getTag('proteomics.de_complete')`) before proceeding
- Analysis functions operate on single DataFrame in-place, adding columns — viewers should read from same DataFrame

### Integration Points
- DE step sets `df.setTag('proteomics.de_complete', 'true')` — viewers can check this tag
- DE adds columns: `log2FC` (SEMTYPE.LOG2FC), `p-value` (SEMTYPE.P_VALUE), `adj.p-value` (SEMTYPE.P_VALUE), `significant` (bool)
- Group annotations stored as DataFrame metadata via `getGroups()` — needed for heatmap grouping and PCA coloring
- Intensity columns from import have semantic types for identification

</code_context>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 03-visualization*
*Context gathered: 2026-02-28*
