# Phase 5: Gap Closure & Hardening - Context

**Gathered:** 2026-03-03
**Status:** Ready for planning

<domain>
## Phase Boundary

Fix four discrete integration issues from the v1.0 milestone audit and render the dendrogram tree visual alongside the heatmap. The four issues: R script function name mismatch, detector regex gap for log2-transformed columns, heatmap filter leakage to volcano plot, and missing dendrogram tree rendering.

</domain>

<decisions>
## Implementation Decisions

### R Script Name Resolution
- Verify empirically at runtime which name the platform registers (camelCase `limmaDE` from `#name:` vs PascalCase `LimmaDE` from `grok api`)
- Fix whichever side is wrong to align both to the canonical registered name
- When R environment is unavailable and user chose limma/DEqMS, show a warning notification ("R environment unavailable — using client-side t-test") rather than silently falling back

### Heatmap Filter Isolation
- Clone the DataFrame for the heatmap display (matches PCA pattern already in the package)
- Selection linkage loss from cloning is acceptable — the heatmap is a derived view (z-scores of top-50)
- Z-score columns stay in the cloned DataFrame only — main grid stays clean

### Dendrogram Rendering
- Graceful fallback if Dendrogram package not installed: show heatmap with clustered row order but no tree visual, log a console warning

### Detector Robustness
- Fix the log2-transformed column name regression so semTypes survive project save/reload

### Claude's Discretion
- Whether to use `injectTreeForGridUI2` (lower-level, any Grid) or adapt `hierarchicalClusteringUI` (turnkey, TableView-oriented) — use whichever fits the embedded grid architecture better
- Whether to replace hand-rolled Euclidean distance matrix with `TreeHelper.calcDistanceMatrix()` — refactor if it simplifies the dendrogram integration
- Dendrogram panel width (100-150px range, pick a reasonable default)
- Whether to use package-api.ts generated wrappers vs manual `grok.functions.call()` strings — do whatever is more maintainable
- Detector regex pattern for matching log2-transformed column names — be robust enough to survive save/reload without false positives
- Whether log2-transformed columns get same semType ('Proteomics-Intensity') or a new one — minimize downstream code changes
- Whether to make heatmap top-N count configurable or keep fixed at 50
- Whether to add detector patterns for non-MaxQuant formats (Spectronaut, DIANN) — fix the log2 regression first, add others only if trivial

</decisions>

<specifics>
## Specific Ideas

- The Dendrogram package provides `hierarchicalClusteringUI` as a turnkey solution (handles null filtering, clustering, tree injection, loader UI, cleanup, persistence) — but it's designed for TableView grids. The lower-level `injectTreeForGridUI2` works on any `DG.Grid` and is what `hierarchicalClusteringUI` delegates to internally (line 188 of hierarchical-clustering.ts)
- The heatmap already calls `getDendrogramService()` and hand-rolls a distance matrix — the gap is calling the tree injection API after clustering

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `injectTreeForGridUI2` (Dendrogram package): Injects tree visual as GridNeighbor to any DG.Grid — handles rendering, scrolling, zooming, context menu, resize
- `hierarchicalClusteringUI` (Dendrogram package): Full turnkey solution for TableView grids — handles null filtering, loader, cleanup, persistence via onInitializedScript
- `TreeHelper` (@datagrok-libraries/bio): `calcDistanceMatrix()` for standard distance computation
- `getClusterMatrixWorker` (@datagrok-libraries/math): WASM-accelerated hierarchical clustering
- `package-api.ts`: Auto-generated typed wrappers for R script functions (`DeqmsDE`, `LimmaDE`)

### Established Patterns
- PCA viewer already uses a cloned DataFrame for filter isolation — heatmap should follow this pattern
- Dendrogram package uses `getDendrogramService()` → `TreeHelper` → `getClusterMatrixWorker` → `injectTreeForGridUI2` pipeline
- `detectors.js` uses column name prefix matching + value range guards for semType detection

### Integration Points
- `heatmap.ts:111` — filter mutation point to fix (clone DataFrame instead)
- `heatmap.ts:138-185` — existing clustering code to extend with tree injection
- `heatmap.ts:159-170` — hand-rolled distance matrix to potentially replace
- `differential-expression.ts:135,192` — R script call sites to align with registered names
- `detectors.js:69-79` — `detectIntensity` function to extend with log2 pattern

</code_context>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 05-gap-closure-and-hardening*
*Context gathered: 2026-03-03*
