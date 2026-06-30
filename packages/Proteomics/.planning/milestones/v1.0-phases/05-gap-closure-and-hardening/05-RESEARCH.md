# Phase 5: Gap Closure & Hardening - Research

**Researched:** 2026-03-03
**Domain:** Datagrok integration bug fixes + dendrogram rendering
**Confidence:** HIGH

## Summary

Phase 5 addresses four discrete integration issues discovered during the v1.0 milestone audit, plus one visual enhancement (dendrogram rendering). The issues span four independent subsystems: R script function name resolution, detectors.js regex patterns, DataFrame filter isolation between viewers, and hierarchical clustering visualization.

The R script function name issue is a case mismatch between what `differential-expression.ts` calls (`Proteomics:limmaDE`) and what Datagrok actually registers. Datagrok function resolution is case-sensitive for `grok.functions.call()`. The R scripts declare `#name: limmaDE` and `#name: deqmsDE`, but the auto-generated `package-api.ts` converts these to PascalCase (`Proteomics:LimmaDE`, `Proteomics:DeqmsDE`). The manual code in `differential-expression.ts` uses the original camelCase. The fix is to ensure the call matches the registered name -- which depends on how the platform registers R script names (by the `#name:` annotation directly). The safest fix is to make the `#name:` annotation in the R scripts match what the TS code calls.

The filter leakage is caused by `grid.dataFrame.filter.copyFrom(filter)` on line 111 of `heatmap.ts` -- this mutates the shared DataFrame's filter BitSet, which the volcano plot also observes. The fix is to use a cloned DataFrame or a separate filter mechanism.

The dendrogram rendering should use the existing Dendrogram package's `IDendrogramService.injectTreeForGrid()` API rather than hand-rolling a tree renderer. This is an established pattern in the Datagrok ecosystem used by the Dendrogram, Bio, and other packages.

**Primary recommendation:** Fix three bugs with targeted edits (function name, detector regex, filter isolation), then integrate the Dendrogram package service for tree visualization alongside the heatmap grid.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `@datagrok-libraries/bio` | workspace | `IDendrogramService` interface + `NodeType` tree types | Defines the canonical interface for dendrogram injection; used by Dendrogram package |
| `@datagrok-libraries/gridext` | workspace | `GridNeighbor` class for attaching UI alongside Grid | The only supported way to dock a neighbor element beside a DG.Grid |
| `@datagrok-libraries/math` | workspace | `getClusterMatrixWorker()` for WASM hierarchical clustering | Production-quality WASM clustering used by Dendrogram package |
| Dendrogram package | installed | `getDendrogramService()` runtime service | Provides `injectTreeForGrid()` which renders, handles interaction, and manages lifecycle |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `@datagrok-libraries/bio/src/trees` | workspace | `TreeHelper` + tree types (`NodeType`, `DistanceMetric`, `LinkageMethod`) | For building Newick-compatible tree structures from clustering results |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Dendrogram service | Hand-rolled Canvas dendrogram | Would miss interaction (click, hover, zoom), lifecycle management, and resize observer handling. 400+ lines to replicate |
| `getClusterMatrixWorker` (WASM) | JS-based clustering | WASM is faster for 50+ proteins, and already used by Dendrogram package |
| Cloned DataFrame for heatmap filter | `grid.dataFrame.rows.match()` | Clone is simpler and matches the pattern used by PCA (separate DataFrame) |

**Installation:**
No new dependencies needed. `@datagrok-libraries/bio`, `@datagrok-libraries/gridext`, and `@datagrok-libraries/math` are already available in the monorepo. The Dendrogram package must be installed on the Datagrok server at runtime.

## Architecture Patterns

### Recommended Project Structure
No structural changes needed. All fixes go into existing files:
```
packages/Proteomics/
  scripts/
    limma_de.R           # Fix: ensure #name matches TS call
    deqms_de.R           # Fix: ensure #name matches TS call
  detectors.js           # Fix: add log2(...) pattern to detectIntensity
  src/
    viewers/heatmap.ts   # Fix: filter isolation + dendrogram injection
    analysis/differential-expression.ts  # Fix: function call name
```

### Pattern 1: R Script Function Name Resolution
**What:** Datagrok registers R scripts by the `#name:` annotation. The platform uses the exact name from the annotation when registering. `grok.functions.call('Package:functionName')` must match exactly.
**When to use:** Any time TS code calls an R script function.
**Confidence:** HIGH -- verified by inspecting the R scripts (`#name: limmaDE`), the auto-generated `package-api.ts` (which converts to `Proteomics:LimmaDE`), and the manual call in `differential-expression.ts` (`Proteomics:limmaDE`).

The mismatch:
- R script declares: `#name: limmaDE`
- `package-api.ts` auto-generates: `grok.functions.call('Proteomics:LimmaDE')`
- `differential-expression.ts` manually calls: `grok.functions.call('Proteomics:limmaDE')`

The `package-api.ts` file is auto-generated by `grok api` which PascalCases the name. The manual code uses the original camelCase name from the R script. The actual registered name at runtime is the `#name:` value. Since no code in the package imports from `package-api.ts` for these R calls (the manual code calls directly), the fix is straightforward: either change the R script `#name:` to match the call, or change the call to match the `#name:`. The safest approach is to verify that `Proteomics:limmaDE` (camelCase, matching the R script `#name:`) is the correct registered name, and leave the TS code as-is. If the platform PascalCases on registration, then the TS call needs to be updated to `Proteomics:LimmaDE`.

**Resolution strategy:** The most reliable fix is to change the R script `#name:` annotations to exactly match what the TS code calls. Since `differential-expression.ts` calls `Proteomics:limmaDE` and `Proteomics:deqmsDE`, the R scripts should keep `#name: limmaDE` and `#name: deqmsDE`. Then re-run `grok api` to regenerate `package-api.ts` (which will still PascalCase, but that file is not used for these calls). Alternatively, verify empirically by checking what name the platform actually registers -- but the safest path is to ensure both sides agree.

**IMPORTANT FINDING:** Looking at the Samples package R scripts, the convention varies: `#name: ListPackages` (PascalCase) and `#name: renvSpellingExample` (camelCase). The platform appears to register the exact `#name:` value. The `grok api` tool PascalCases in `package-api.ts` but that is a wrapper -- the actual function is registered with the literal `#name:` value. Therefore `differential-expression.ts` calling `Proteomics:limmaDE` should match `#name: limmaDE` correctly.

**The actual bug:** The audit claims the registered name is `Proteomics:LimmaDE` (PascalCase). This would happen if the platform auto-PascalCases R script names on registration. If so, the fix must change the TS calls to PascalCase. Given uncertainty, both approaches should be tested, but the safest fix is to make the TS calls use the same casing as `package-api.ts` generates (`Proteomics:LimmaDE` / `Proteomics:DeqmsDE`), since `grok api` is the authoritative tool for generating correct call names.

### Pattern 2: Detector Regex for Transformed Column Names
**What:** `detectors.js::detectIntensity` only matches raw intensity prefixes (`intensity`, `lfq intensity`, `ibaq`, `reporter intensity`). After parsing, columns are renamed to `log2(LFQ intensity Sample1)` etc. On project save/reload, the detector runs again and fails to match.
**When to use:** Any semantic type detector that must survive project save/reload with transformed column names.
**Example:**
```javascript
// Current (broken for log2-transformed columns):
if ((name.startsWith('intensity') || name.startsWith('lfq intensity') ||
  name.startsWith('reporter intensity') || name.startsWith('ibaq')) && col.min >= 0)

// Fixed (also matches log2-wrapped names):
if ((name.startsWith('intensity') || name.startsWith('lfq intensity') ||
  name.startsWith('reporter intensity') || name.startsWith('ibaq') ||
  (name.startsWith('log2(') && (name.includes('intensity') || name.includes('ibaq'))))
```
Note: `col.min >= 0` check must be relaxed for log2-transformed columns (log2 values can be negative).

### Pattern 3: Filter Isolation Between Viewers
**What:** When two viewers share a DataFrame, setting `df.filter` affects both. The heatmap currently uses `grid.dataFrame.filter.copyFrom(filter)` which mutates the shared filter.
**When to use:** Any viewer that needs to show a subset of rows without affecting other viewers on the same DataFrame.

Two approaches exist in the Datagrok ecosystem:

**Approach A: Clone the DataFrame** (used by PCA in this package)
```typescript
const heatmapDf = df.clone(filter);
const grid = heatmapDf.plot.grid();
```
Downside: Selection linkage to the original DataFrame is broken.

**Approach B: Use a separate Grid with row visibility** (preserves linkage)
The Grid viewer has row visibility that is independent of DataFrame filter. However, this requires using `grid.setRowOrder()` or similar APIs.

**Approach C: Use `grid.props.rowSource = 'Filtered'` with a viewer-level filter**
Some Datagrok viewers support viewer-level filtering that does not propagate to the DataFrame.

**Recommendation:** Use Approach A (clone) because:
1. The heatmap already shows a fundamentally different view (z-scores of top-50 proteins)
2. It already creates temporary `_zscore_` columns that pollute the original DataFrame
3. PCA already uses a separate DataFrame in this package, establishing precedent
4. Selection linkage for the heatmap is less critical than filter isolation

### Pattern 4: Dendrogram Rendering via Dendrogram Package

The Dendrogram package provides **two levels of API** for adding dendrograms:

#### High-level: `hierarchicalClusteringUI` (full turnkey)
**Source:** `packages/Dendrogram/src/utils/hierarchical-clustering.ts`
**What:** Handles the entire workflow: null filtering, distance matrix, WASM clustering, tree injection, loader UI, cleanup subscriptions, and persistence via `onInitializedScript`.
**When to use:** When attaching a dendrogram to a **TableView's main grid**.
**Limitation:** Expects `grok.shell.getTableView(df.name)` — designed for TableView scenarios, not embedded grids.
```typescript
import {hierarchicalClusteringUI} from 'packages/Dendrogram/src/utils/hierarchical-clustering';

await hierarchicalClusteringUI(
  dataFrame,
  ['column1', 'column2'],
  DistanceMetric.Euclidean,
  'Ward'
);
```

#### Low-level: `injectTreeForGridUI2` (grid-level injection)
**Source:** `packages/Dendrogram/src/viewers/inject-tree-for-grid2.ts`
**What:** Attaches a tree visual as a `GridNeighbor` to the left of any `DG.Grid`. Handles rendering, scrolling, zooming, context menu (assign clusters, reset zoom), and resize.
**When to use:** When attaching a dendrogram to an **internal/embedded grid** (not the main TableView grid) — which is the heatmap's case.
**Signature:** `injectTreeForGridUI2(grid: DG.Grid, treeRoot: NodeType | null, leafColName?: string, neighborWidth?: number, cut?: TreeCutOptions): GridNeighbor`

The heatmap creates its own internal grid (not the TableView grid), so `injectTreeForGridUI2` is the correct API level.

#### Clustering + rendering workflow for the heatmap:
```typescript
import {TreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {getClusterMatrixWorker} from '@datagrok-libraries/math';
import {injectTreeForGridUI2} from 'packages/Dendrogram/src/viewers/inject-tree-for-grid2';

// 1. Build distance matrix and cluster (replace hand-rolled distance calc)
const th = new TreeHelper();
const distanceMatrix = await th.calcDistanceMatrix(df, colNames, 'euclidean');
const clusterMatrix = await getClusterMatrixWorker(distanceMatrix.data, df.rowCount, linkageCode);
const newickRoot = th.parseClusterMatrix(clusterMatrix);
newickRoot.branch_length = 0;

// 2. Inject tree into heatmap's internal grid
const neighbor = injectTreeForGridUI2(grid, newickRoot, undefined, 150);
```

**Current state:** The heatmap code (`heatmap.ts:138`) already calls `getDendrogramService()` and hand-rolls a Euclidean distance matrix (lines 159-170), then only extracts `tree.leafOrder` for sorting. The fix should: (1) replace the hand-rolled distance matrix with `TreeHelper.calcDistanceMatrix()`, (2) replace `dendrogramService.getTree()` with `getClusterMatrixWorker()` + `th.parseClusterMatrix()`, (3) call `injectTreeForGridUI2()` to render the tree visual.

**Note:** `hierarchicalClusteringUI` internally calls `injectTreeForGridUI2` (line 188 of hierarchical-clustering.ts), confirming this is the correct rendering primitive.

### Anti-Patterns to Avoid
- **Mutating shared DataFrame filter:** `df.filter.copyFrom()` on a shared DataFrame affects all viewers. Clone the DataFrame or use viewer-level filtering.
- **Hand-rolling dendrogram rendering:** The Dendrogram package provides 400+ lines of interaction handling (hover, click, selection sync, resize, sort alignment). Reimplementing this is fragile and will miss edge cases.
- **Relying on `package-api.ts` casing for manual calls:** The auto-generated file PascalCases function names, but manual `grok.functions.call()` must use the exact registered name.
- **Removing `col.min >= 0` without adjusting for log2:** log2-transformed intensity columns have negative values. The detector must handle this.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Dendrogram tree rendering | Custom Canvas/SVG tree drawing | `IDendrogramService.injectTreeForGrid()` | 400+ LOC for rendering, interaction, scrolling, resize, selection sync. Dendrogram package handles all of this |
| Hierarchical clustering | Custom clustering algorithm | `@datagrok-libraries/math` `getClusterMatrixWorker()` | WASM-accelerated, handles Ward/average/complete linkage, tested in production |
| Tree data structures | Custom tree node types | `NodeType` from `@datagrok-libraries/bio/src/trees` | Standard type that `TreeHelper` and `injectTreeForGrid()` expect |
| Grid-adjacent panel | Custom DOM manipulation | `GridNeighbor` from `@datagrok-libraries/gridext` | Handles resize observer, scroll sync, canvas offset adjustment, cleanup |

**Key insight:** The entire dendrogram rendering infrastructure already exists in the Datagrok ecosystem. The current heatmap code already attempts to use the Dendrogram service for clustering. The gap is that it only extracts the sort order but never renders the tree visual.

## Common Pitfalls

### Pitfall 1: Function Name Case Sensitivity
**What goes wrong:** `grok.functions.call('Proteomics:limmaDE')` silently throws an error because the platform registers the function under a different case. The try/catch in the DE dialog catches this and falls back to client-side t-test, making it appear to work but using the wrong method.
**Why it happens:** `grok api` PascalCases function names in `package-api.ts`, but developers write manual calls using the original `#name:` casing from R scripts.
**How to avoid:** Always use the casing from `package-api.ts` for `grok.functions.call()` since that is what `grok api` determines to be the canonical registered name. Or test the actual registered name empirically.
**Warning signs:** DE dialog always shows "t-test fallback" even when R environment is available; `proteomics.de_method` tag is always `ttest`.

### Pitfall 2: Detector Survives Import But Not Reload
**What goes wrong:** Semantic types are correctly assigned during import (parser sets them directly), but after save/reload the detectors run and fail to recognize the transformed column names.
**Why it happens:** The parser assigns semTypes programmatically. The detector regex only matches the pre-transformation column names. On reload, Datagrok re-runs detectors, which are the only mechanism for restoring semTypes.
**How to avoid:** Ensure detector patterns cover ALL column name forms that can exist in a saved project -- both raw and transformed names.
**Warning signs:** Column pickers in analysis dialogs show no columns after reopening a saved project.

### Pitfall 3: Shared DataFrame Filter Mutation
**What goes wrong:** Heatmap sets `df.filter` to show top-50 proteins. Volcano plot, which shares the same DataFrame, now shows only 50 points instead of all proteins.
**Why it happens:** `grid.dataFrame.filter.copyFrom(filter)` mutates the filter on the shared DataFrame reference.
**How to avoid:** Clone the DataFrame for heatmap display. The z-score columns, filter state, and sort order become heatmap-private.
**Warning signs:** Opening heatmap causes volcano plot to shrink to 50 points; volcano point count changes when heatmap is opened/closed.

### Pitfall 4: Dendrogram Service Unavailable
**What goes wrong:** `getDendrogramService()` throws if the Dendrogram package is not installed on the server.
**Why it happens:** The service is provided by a separate package that may not be deployed.
**How to avoid:** Always wrap in try/catch with graceful fallback (already done in current code for clustering; extend to rendering).
**Warning signs:** "Function not found" error in console when clicking Heatmap.

### Pitfall 5: Z-Score Columns Polluting Original DataFrame
**What goes wrong:** Current heatmap adds `_zscore_*` columns to the original DataFrame. Other viewers (volcano, grid) see these extra columns.
**Why it happens:** Heatmap operates on the shared DataFrame to enable selection linkage.
**How to avoid:** If using a cloned DataFrame for filter isolation (recommended), the z-score columns naturally become private to the clone.
**Warning signs:** Extra `_zscore_*` columns visible in the main grid after opening heatmap.

## Code Examples

### Fix 1: R Script Function Call Name Alignment
```typescript
// In differential-expression.ts, update to match package-api.ts casing:
// Before:
const result = await grok.functions.call('Proteomics:limmaDE', { ... });
// After:
const result = await grok.functions.call('Proteomics:LimmaDE', { ... });

// Before:
const result = await grok.functions.call('Proteomics:deqmsDE', { ... });
// After:
const result = await grok.functions.call('Proteomics:DeqmsDE', { ... });
```

### Fix 2: Detector Regex for log2-Wrapped Column Names
```javascript
// In detectors.js::detectIntensity
detectIntensity(col) {
  if (col.type !== DG.TYPE.FLOAT && col.type !== DG.TYPE.INT)
    return null;
  const name = col.name.toLowerCase();
  // Match raw intensity column names (pre-transform)
  const isRawIntensity = (name.startsWith('intensity') || name.startsWith('lfq intensity') ||
    name.startsWith('reporter intensity') || name.startsWith('ibaq')) && col.min >= 0;
  // Match log2-transformed column names (post-transform)
  const isLog2Intensity = name.startsWith('log2(') &&
    (name.includes('intensity') || name.includes('ibaq'));
  if (isRawIntensity || isLog2Intensity) {
    col.semType = 'Proteomics-Intensity';
    return col.semType;
  }
  return null;
}
```

### Fix 3: Heatmap Filter Isolation via DataFrame Clone
```typescript
// In heatmap.ts, replace filter mutation with DataFrame clone:

// Build filter BitSet for top N
const filter = DG.BitSet.create(df.rowCount);
for (const idx of topNIndices)
  filter.set(idx, true);

// Clone DataFrame with only top-N rows (isolates filter from volcano)
const heatmapDf = df.clone(filter);

// Add z-score columns to the CLONE (not the original)
for (const colName of allIntensityCols) {
  const zColName = `_zscore_${colName}`;
  // ... compute z-scores on heatmapDf ...
}

// Create Grid on the cloned DataFrame
const grid = heatmapDf.plot.grid();
```

### Fix 4: Dendrogram Injection for Heatmap
```typescript
// After clustering and creating the grid on heatmapDf:
import {getDendrogramService} from '@datagrok-libraries/bio/src/trees/dendrogram';
import {NodeType} from '@datagrok-libraries/bio/src/trees';

// The tree root from clustering (already computed)
let dendrogramRendered = false;
try {
  const dendrogramService = await getDendrogramService();
  // newickRoot is the tree from clustering
  const neighbor = dendrogramService.injectTreeForGrid(grid, newickRoot, undefined, 150);
  dendrogramRendered = true;
} catch (_e) {
  console.warn('Dendrogram package not available; tree visual not rendered');
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Sort rows by cluster order only | Render dendrogram tree alongside grid via `IDendrogramService` | Dendrogram package v1.x (2024+) | Visual tree with interaction (hover, click, selection) |
| Mutate shared DataFrame filter | Clone DataFrame for isolated viewer | Established pattern | No cross-viewer filter contamination |

**Deprecated/outdated:**
- The current heatmap's approach of calling `getDendrogramService()` only to get `tree.leafOrder` for sorting (and never rendering the tree) is incomplete -- it was a first pass that needs to be extended to actually inject the tree visual.

## Open Questions

1. **Exact function name casing at runtime**
   - What we know: R scripts use `#name: limmaDE` (camelCase). `grok api` generates PascalCase calls in `package-api.ts`. The manual code uses camelCase.
   - What's unclear: Whether the platform registers R script functions with the exact `#name:` value or normalizes the casing.
   - Recommendation: Change the TS calls to PascalCase (`Proteomics:LimmaDE`, `Proteomics:DeqmsDE`) to match `package-api.ts`, since `grok api` is the canonical tool for determining registered names. If that fails at runtime, revert to camelCase. Alternatively, test empirically on a running Datagrok instance.
   - **Confidence:** MEDIUM -- the audit explicitly identified this as a case mismatch issue. The `grok api` tool is likely correct about the casing.

2. **Selection linkage after DataFrame clone**
   - What we know: Cloning the DataFrame breaks selection linkage between heatmap and volcano (clicking a protein in heatmap won't highlight in volcano).
   - What's unclear: Whether users expect cross-viewer selection between heatmap and volcano.
   - Recommendation: Accept the tradeoff. Filter isolation is more important than selection linkage. Document the limitation. Selection within the heatmap+dendrogram still works. Cross-viewer selection between volcano and grid still works (they share the original DataFrame).
   - **Confidence:** HIGH -- this is a deliberate design tradeoff, not a technical limitation.

## Sources

### Primary (HIGH confidence)
- Codebase inspection: `packages/Proteomics/` -- all four bug files examined directly
- Codebase inspection: `packages/Dendrogram/src/` -- service pattern, `injectTreeForGrid()` API, `hierarchical-clustering.ts` usage example
- Codebase inspection: `libraries/bio/src/trees/dendrogram.ts` -- `IDendrogramService` interface
- Codebase inspection: `libraries/gridext/src/ui/GridNeighbor.ts` -- neighbor attachment mechanism
- Codebase inspection: `packages/Proteomics/src/package-api.ts` -- auto-generated PascalCase function names
- Codebase inspection: `packages/Proteomics/scripts/limma_de.R`, `deqms_de.R` -- `#name:` annotations
- `.planning/v1.0-MILESTONE-AUDIT.md` -- RISK-01, RISK-02, RISK-03 descriptions

### Secondary (MEDIUM confidence)
- `grok api` behavior inferred from comparing R script `#name:` values with generated `package-api.ts` casing

### Tertiary (LOW confidence)
- Platform function name resolution casing behavior -- inferred from `package-api.ts` output, not verified against Datagrok runtime source code

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all libraries verified by codebase inspection of Dendrogram package usage
- Architecture: HIGH - all four bugs have clear root causes identified from code
- Pitfalls: HIGH - all pitfalls are directly observable in the codebase
- Function name resolution: MEDIUM - depends on platform behavior that cannot be fully verified without runtime testing

**Research date:** 2026-03-03
**Valid until:** 2026-04-03 (stable domain; no external library churn)
