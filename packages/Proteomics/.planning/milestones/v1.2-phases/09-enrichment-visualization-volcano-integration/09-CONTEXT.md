# Phase 9: Enrichment Visualization & Volcano Integration - Context

**Gathered:** 2026-03-07
**Status:** Ready for planning

<domain>
## Phase Boundary

Scientists can visually explore enrichment results through dot plot and bar chart viewers, and link enrichment terms back to proteins on the volcano plot via selection highlighting. This phase adds visualization and cross-DataFrame interaction on top of the enrichment analysis and table from Phase 8.

</domain>

<decisions>
## Implementation Decisions

### Visualization Placement
- Dot plot and bar chart dock into the enrichment table view (same view as the enrichment results grid)
- Layout: dot plot and bar chart side-by-side above the enrichment table grid
- Full cross-linking within the enrichment view — clicking a dot or bar highlights the corresponding row in the grid (native Datagrok viewer-DataFrame selection sync)

### Dot Plot Design
- Scatter plot with categorical Y axis (Term Name), X axis (Gene Ratio), size mapped to Gene Count, color mapped to -log10(FDR) gradient
- -log10(FDR) color gradient — standard in proteomics tools (clusterProfiler, g:Profiler)

### Bar Chart Design
- Bar chart of top enriched terms ranked by significance (FDR) or gene count
- Uses Datagrok's built-in bar chart viewer

### Top-N and Filtering
- Default: top 15 terms by FDR
- No redundancy reduction — just top-N by FDR, simple and transparent
- Users can change top-N via Datagrok's viewer property panel (right sidebar)
- Source filtering (GO:BP, GO:MF, KEGG, etc.) via table grid column filter — dot plot and bar chart auto-update through native filter-viewer sync

### Volcano Highlighting UX
- Click a row in enrichment table → automatically highlights member proteins on the volcano plot
- Trigger: onCurrentRowChanged event on enrichment DataFrame
- Visual feedback: Datagrok selection (blue dots) + gene name labels on selected points
- No auto-zoom — keep current volcano viewport, just highlight in place
- Clear: click empty area on volcano (standard Datagrok behavior), or click a different enrichment row to replace

### Entry Point & Workflow
- Auto-open: when enrichment analysis completes, enrichment table view opens WITH dot plot and bar chart already docked (modify existing enrichment dialog's OK handler)
- Re-open menu: add 'Proteomics | Visualize | Enrichment Charts...' for existing enrichment tables (detected by proteomics.enrichment tag)
- Auto-detect protein table: find source protein DataFrame by proteomics.de_complete tag, subscribe to onCurrentRowChanged immediately — zero setup

### Claude's Discretion
- Exact DockManager arrangement and sizing for side-by-side dot plot + bar chart above grid
- -log10(FDR) column computation approach (add computed column vs inline)
- Gene name label rendering on selected volcano points (labelColumnNames property vs custom approach)
- How to find and reference the volcano ScatterPlotViewer from the protein table view
- Subscription lifecycle management (cleanup on view close)

</decisions>

<specifics>
## Specific Ideas

- Dot plot follows clusterProfiler convention: Y axis = term names, X axis = gene ratio, dot size = gene count, dot color = -log10(FDR)
- Enrichment DataFrame already has proteomics.enrichment tag (set in Phase 8's buildEnrichmentDf) — use this for re-open detection
- Intersection column contains comma-separated gene symbols — parse these to find matching rows in protein DataFrame
- Phase 8's enrichment dialog OK handler at enrichment.ts:398 is the integration point for auto-opening viewers

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `buildEnrichmentDf()` in analysis/enrichment.ts: creates enrichment DataFrame with Gene Count, Gene Ratio, FDR, Intersection columns — all needed for dot plot/bar chart
- `showEnrichmentDialog()` OK handler (enrichment.ts:366-405): integration point for auto-opening viewers after enrichment completes
- `createVolcanoPlot()` in viewers/volcano.ts: creates ScatterPlotViewer — need to find existing volcano instances on protein table
- `findColumn()` in utils/column-detection.ts: find columns by semantic type
- QC Dashboard docking pattern in viewers/qc-dashboard.ts: model for docking multiple viewers into a table view

### Established Patterns
- Viewers docked via `tv.dockManager.dock()` with DockNode references (Phase 7 QC Dashboard)
- DataFrame clone with BitSet for top-N filtering without affecting source data (heatmap.ts pattern)
- Separate TableView for different-shaped data (PCA, enrichment table)
- proteomics.enrichment tag on enrichment DataFrame (set in buildEnrichmentDf line 182)
- proteomics.de_complete tag on protein DataFrame for detecting DE completion

### Integration Points
- Modify `showEnrichmentDialog()` OK handler to call new enrichment visualization function after `grok.shell.addTableView()`
- New file: `src/viewers/enrichment-viz.ts` for dot plot, bar chart creation and volcano linking
- New menu entry in PackageFunctions class for 'Enrichment Charts...' re-open
- Cross-DataFrame event subscription: enrichment onCurrentRowChanged → protein DataFrame selection

</code_context>

<deferred>
## Deferred Ideas

- GO semantic similarity clustering (REVIGO-like redundancy reduction) — future enhancement
- Auto-zoom volcano to fit selected proteins — could add as option later
- Source-specific color coding in dot plot (different colors per GO:BP, KEGG, etc.) — future enhancement

</deferred>

---

*Phase: 09-enrichment-visualization-volcano-integration*
*Context gathered: 2026-03-07*
