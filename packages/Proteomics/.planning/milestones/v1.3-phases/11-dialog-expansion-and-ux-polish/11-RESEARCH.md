# Phase 11: Dialog Expansion and UX Polish - Research

**Researched:** 2026-03-07
**Domain:** Datagrok UI dialogs, embedded viewers, reactive inputs, DataFrame naming
**Confidence:** HIGH

## Summary

Phase 11 is a pure UI/UX phase that wires Phase 10's algorithms (quantile normalization, VSN, kNN imputation, zero/mean/median imputation) into polished dialogs with method selectors, conditional parameter visibility, and visual feedback. It also adds descriptive viewer titles and standardizes DataFrame naming from filenames. No new analysis algorithms are needed.

The codebase already contains every pattern required: conditional visibility (`differential-expression.ts:276-279`), live count updates (`enrichment.ts:344-355`), viewer-in-dialog embedding (`generic-parser.ts:86-115`), viewer title configuration (`qc-dashboard.ts:62-66`), box plot factory (`DG.Viewer.boxPlot()` used in `qc-dashboard.ts:127`), and long-format data preparation (`unpivotIntensities()` in `qc-computations.ts`). The work is assembly of proven patterns, not invention.

**Primary recommendation:** Structure implementation as three dialog-focused waves (normalization, imputation, DE) plus one cross-cutting UX wave (viewer titles + DataFrame naming). Each dialog follows the same formula: add method selector via `ui.input.choice()`, wire conditional visibility via `onChanged.subscribe()` + `style.display`, add visual feedback element.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Normalization dialog: Box plot (one box per sample) for intensity distribution; reactive preview updates on method change via DataFrame clone; before-state only (no side-by-side); inline yellow/orange banner warning for Spectronaut pre-normalized data (`proteomics.preNormalized` tag); warning text: "This data may be pre-normalized (Spectronaut). Additional normalization may distort results."; non-blocking warning
- Imputation dialog: Method selector dropdown with conditional inline parameters (MinProb: downshift + width; kNN: k neighbors; Zero/Mean/Median: no extra params); per-group minimum valid values filter (proteins failing filter are removed, not hidden); live count updates: "Will keep X/Y proteins (Z removed)"; same reactive pattern as enrichment dialog
- DE dialog: Single dropdown for comparison direction with auto-generated pairs ("Treatment vs Control", "Control vs Treatment"); dynamic hint text below dropdown: "Positive log2FC = higher in [first group], Negative log2FC = higher in [second group]"; comparison picker inserted at top of dialog; consistent conditional visibility (DEqMS shows peptide column, Limma and t-test show no extra params)
- Viewer titles: Format "Viewer: Context" (e.g., "Volcano: Treatment vs Control", "PCA: All Groups", "Heatmap: Top 50 DE Proteins"); set via `setOptions({title: '...'})` or viewer creation options
- DataFrame naming: All import paths use original filename minus extension; auxiliary DataFrames include source context (e.g., "PCA: proteinGroups", "Enrichment Results: proteinGroups")

### Claude's Discretion
- Exact box plot sizing and positioning within dialog
- Clone strategy for reactive normalization preview (full clone vs column subset)
- Default minimum valid values threshold (likely 2 or 3)
- Exact wording of info messages and tooltips

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| NORM-03 | User can select normalization method from dialog with before/after distribution plots | Box plot viewer embedding pattern from QC dashboard; `DG.Viewer.boxPlot()` factory; `unpivotIntensities()` for long-format conversion; reactive clone+update pattern |
| NORM-04 | Normalization dialog warns if Spectronaut pre-normalized data detected | `proteomics.preNormalized` tag already set by Spectronaut parser; inline styled div pattern |
| IMP-03 | User can select imputation method from dialog with conditional parameters per method | Conditional visibility pattern from DE dialog (`peptideRow.style.display` toggle); `ui.input.choice()` for method selector |
| IMP-04 | User can filter proteins by minimum valid values before imputation | Live count pattern from enrichment dialog (`onChanged.subscribe` -> update div); per-group valid counting via `getGroups()` |
| DE-01 | User can select comparison direction from auto-generated group pairs | `getGroups()` returns group1/group2; generate two directional pairs; parse selected string to determine column ordering |
| DE-02 | DE dialog shows/hides method-specific parameters based on selected method | Existing conditional visibility already works for DEqMS peptide column; extend to t-test (third method, no extra params) |
| UX-01 | Volcano, PCA, and heatmap viewers display descriptive titles | QC dashboard proves `title` option works in viewer creation options (`title: 'MA Plot'`) |
| UX-02 | DataFrame retains imported filename as its name | MaxQuant import in `package.ts:66` already does `df.name = file.name.replace(/\.[^.]+$/, '')`; Spectronaut parser needs fix (currently hardcodes `'spectronaut'`) |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `datagrok-api` | Current | All UI: dialogs, inputs, viewers, DataFrame | Platform API -- only option |

### Supporting
No additional libraries needed. All Phase 11 work uses existing Datagrok API components already imported in the codebase.

### Alternatives Considered
None -- all decisions are locked.

**Installation:**
```bash
# No new dependencies required
```

## Architecture Patterns

### Files to Modify
```
packages/Proteomics/src/
  analysis/
    normalization.ts          # MODIFY: expand showNormalizationDialog() with method selector + box plot preview + warning
    imputation.ts             # MODIFY: expand showImputationDialog() with method selector + conditional params + valid-values filter
    differential-expression.ts # MODIFY: add comparison direction picker to showDEDialog(), add t-test to method list
  viewers/
    volcano.ts                # MODIFY: add title param to createVolcanoPlot()
    heatmap.ts                # MODIFY: add title param to createExpressionHeatmap()
    pca-plot.ts               # MODIFY: add title param to createPcaPlot()
  parsers/
    maxquant-parser.ts        # MODIFY: remove hardcoded df.name = 'proteinGroups' (line 134)
    spectronaut-parser.ts     # MODIFY: remove hardcoded df.name = 'spectronaut' (line 162)
  package.ts                  # MODIFY: pass titles to viewer creation calls, ensure all import paths set df.name from filename
```

### Pattern 1: Method Selector with Conditional Visibility
**What:** A `ui.input.choice()` dropdown that shows/hides related inputs based on selection.
**When to use:** All three dialogs need this pattern.
**Example (proven in existing DE dialog):**
```typescript
// Source: differential-expression.ts:255-279
const methodInput = ui.input.choice('Method', {
  value: 'limma', items: ['limma', 'DEqMS'], nullable: false,
});
const peptideRow = peptideColInput.root;
peptideRow.style.display = 'none';
methodInput.onChanged.subscribe(() => {
  peptideRow.style.display = methodInput.value === 'DEqMS' ? '' : 'none';
});
```

### Pattern 2: Reactive Live Count Update
**What:** A `ui.divText()` element updated via `onChanged.subscribe()` chains.
**When to use:** Imputation dialog's "Will keep X/Y proteins" counter.
**Example (proven in enrichment dialog):**
```typescript
// Source: enrichment.ts:344-355
const countDiv = ui.divText('');
const updateCount = () => {
  const result = countSignificantProteins(df, fc, p, cols.log2fc!, cols.pValue!);
  countDiv.textContent = `${result.significant} of ${totalFmt} proteins are significant`;
};
updateCount();
fcInput.onChanged.subscribe(updateCount);
pInput.onChanged.subscribe(updateCount);
```

### Pattern 3: Viewer Embedded in Dialog
**What:** Create a viewer via `DG.Viewer.boxPlot()` factory, append `.root` to dialog container.
**When to use:** Normalization dialog's distribution preview.
**Example (proven in QC dashboard + generic parser):**
```typescript
// Source: qc-dashboard.ts:127-130 (box plot creation)
const longDf = unpivotIntensities(df, intensityCols);
const boxPlot = DG.Viewer.boxPlot(longDf, {
  valueColumnName: 'Intensity',
  categoryColumnName: 'Sample',
} as any);

// Source: generic-parser.ts:105-109 (viewer root in dialog)
const previewContainer = document.createElement('div');
previewContainer.style.cssText = 'border:1px solid #ddd;margin-top:8px;width:100%;height:250px;';
boxPlot.root.style.cssText = 'width:100%;height:220px;';
previewContainer.appendChild(boxPlot.root);
```

### Pattern 4: Viewer Title Configuration
**What:** Pass `title` property in viewer creation options or via `setOptions()`.
**When to use:** Volcano, PCA, heatmap viewers.
**Example (proven in QC dashboard):**
```typescript
// Source: qc-dashboard.ts:62-66
const maPlot = tv.addViewer(DG.VIEWER.SCATTER_PLOT, {
  xColumnName: 'A', yColumnName: 'M',
  title: 'MA Plot',
  markerDefaultSize: 2,
});
```

### Pattern 5: Comparison Direction Pair Generation
**What:** Generate directional pairs from `getGroups()` for DE comparison dropdown.
**When to use:** DE dialog comparison direction selector.
**Example:**
```typescript
const groups = getGroups(df)!;
const g1 = groups.group1;
const g2 = groups.group2;
const pairs = [
  `${g2.name} vs ${g1.name}`,  // Default: treatment vs control
  `${g1.name} vs ${g2.name}`,  // Reversed
];
const comparisonInput = ui.input.choice('Comparison', {
  value: pairs[0], items: pairs, nullable: false,
});
// Parse on OK to determine column ordering
const selected = comparisonInput.value!;
const reversed = selected.startsWith(g1.name);
const numeratorCols = reversed ? g1.columns : g2.columns;
const denominatorCols = reversed ? g2.columns : g1.columns;
```

### Anti-Patterns to Avoid
- **Modifying original DataFrame during preview:** Always clone before running normalization preview. The original must not change until OK is pressed.
- **Destroying/recreating inputs on method switch:** Create all method-specific inputs upfront, toggle visibility with `style.display`. Destroying inputs breaks `onChanged` subscriptions.
- **Hardcoding group names:** Always derive from `getGroups()`. Never assume "Control" and "Treatment".
- **Creating new viewers without removing old ones:** When updating the box plot preview, remove the old viewer's `.root` from the container before adding the new one.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Box plot visualization | Custom SVG/canvas charts | `DG.Viewer.boxPlot()` | Platform viewer with tooltips, styling, interactivity |
| Long-format data conversion | Manual column melting | `unpivotIntensities()` from `qc-computations.ts` | Already handles intensity column reshaping correctly |
| Conditional visibility | Custom show/hide framework | `element.style.display = condition ? '' : 'none'` | Proven one-liner, no framework needed |
| Live reactive updates | Debounced polling or timers | `input.onChanged.subscribe()` | Fires synchronously on value change |
| Dialog layout | Custom HTML forms | `ui.dialog().add()` chain | Handles buttons, close, keyboard shortcuts |
| Viewer titles | Custom title div above viewer | `setOptions({title: ...})` | Built into all Datagrok viewers |

**Key insight:** Every UI pattern needed in this phase already exists in the codebase. The work is composing existing patterns into new dialog layouts.

## Common Pitfalls

### Pitfall 1: Viewer Sizing Inside Dialogs
**What goes wrong:** Embedded viewers render with zero height or overflow the dialog.
**Why it happens:** Viewers default to filling their container, but dialog content areas lack explicit dimensions.
**How to avoid:** Set explicit `width` and `height` on both the container div and the viewer's `.root` element. Use `width:100%;height:250px;` as baseline. The generic parser already demonstrates this at line 109.
**Warning signs:** Viewer appears as thin line or dialog becomes extremely tall.

### Pitfall 2: Box Plot Requires Long-Format Data
**What goes wrong:** Passing wide-format DataFrame (one column per sample) to `DG.Viewer.boxPlot()` shows a single box.
**Why it happens:** Box plots expect long-format with a value column and category column.
**How to avoid:** Use `unpivotIntensities()` from `qc-computations.ts` to convert wide-format intensity columns to `{Sample, Intensity}` long-format DataFrame. This is exactly what the QC dashboard does.
**Warning signs:** Box plot shows one box instead of per-sample boxes.

### Pitfall 3: VSN Preview Latency
**What goes wrong:** Selecting VSN in the method dropdown triggers a server round-trip for preview, freezing the dialog.
**Why it happens:** VSN requires R script execution; median and quantile are synchronous TypeScript.
**How to avoid:** Per CONTEXT.md decision, the dialog shows "before-state only" -- a single reactive box plot showing current distributions. The preview updates by cloning and running the selected method. For VSN, either show a loading state in the preview area or show the original distributions with a note "VSN preview requires server computation."
**Warning signs:** Dialog becomes unresponsive when VSN is selected.

### Pitfall 4: Valid-Values Filter Logic (ANY vs ALL Groups)
**What goes wrong:** Filter threshold removes proteins that have valid values in some groups but not others.
**Why it happens:** Implementing "at least N valid values in ALL groups" instead of "at least N valid values in at least ONE group."
**How to avoid:** The CONTEXT.md specifies "per-group minimum valid values" and "protein must have at least N valid values in at least one annotated group." Use `getGroups()` to iterate group columns separately; a protein passes if ANY group meets the threshold.
**Warning signs:** Excessive protein removal at low thresholds (e.g., threshold=2 removes 80% of proteins).

### Pitfall 5: DataFrame Name Set in Two Places
**What goes wrong:** Parser sets `df.name = 'proteinGroups'` (maxquant-parser.ts:134) or `df.name = 'spectronaut'` (spectronaut-parser.ts:162), then package.ts overwrites with filename. Or the parser name survives.
**Why it happens:** Naming happens in both the parser function AND the import handler in package.ts.
**How to avoid:** Remove hardcoded names from parser functions (they should not set `df.name` at all). Set `df.name` exclusively in the import handler in package.ts where `file.name` is available. For `parseMaxQuantText()` called from the demo function (no filename), set a default name after the call.
**Warning signs:** DataFrame shows "proteinGroups" or "spectronaut" instead of actual filename.

### Pitfall 6: Comparison Direction Mismatch with DE Computation
**What goes wrong:** User selects "Control vs Treatment" but DE still computes Treatment-vs-Control (group2 vs group1).
**Why it happens:** `runLimmaDE()` and `runDifferentialExpression()` take `group1Cols` and `group2Cols` positionally. The existing code hardcodes group1 as denominator and group2 as numerator.
**How to avoid:** Parse the comparison string to determine ordering. When "Control vs Treatment" is selected, pass Control columns as group2 (numerator) and Treatment as group1 (denominator). The hint text must match: "Positive log2FC = higher in [first group listed]."
**Warning signs:** Hint text says one thing but FC values are in opposite direction.

## Code Examples

### Normalization Dialog with Method Selector and Box Plot Preview
```typescript
// Verified patterns from: qc-dashboard.ts:127-130, enrichment.ts:344-355, differential-expression.ts:255-279
import {unpivotIntensities} from '../viewers/qc-computations';
import {medianNormalize, quantileNormalize, vsnNormalize} from './normalization';

export function showNormalizationDialog(df: DG.DataFrame): void {
  // NORM-04: Pre-normalized warning
  const warningDiv = document.createElement('div');
  warningDiv.style.cssText = 'background:#FFF3CD;border:1px solid #FFEEBA;' +
    'color:#856404;padding:8px;margin:4px 0;border-radius:4px;display:none;';
  warningDiv.textContent = 'This data may be pre-normalized (Spectronaut). ' +
    'Additional normalization may distort results.';
  if (df.getTag('proteomics.preNormalized') === 'true')
    warningDiv.style.display = '';

  // NORM-03: Method selector
  const methodInput = ui.input.choice('Method', {
    value: 'Median Centering',
    items: ['Median Centering', 'Quantile', 'VSN'],
    nullable: false,
  });

  // Columns input
  const log2ColNames = df.columns.toList()
    .filter((c) => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('))
    .map((c) => c.name);
  const log2Cols = log2ColNames.map((n) => df.col(n)!);
  const colsInput = ui.input.columns('Intensity columns', {
    table: df, available: log2ColNames, value: log2Cols,
  });

  // Box plot preview container
  const plotContainer = document.createElement('div');
  plotContainer.style.cssText = 'width:100%;height:250px;border:1px solid #ddd;margin-top:8px;';

  const updatePreview = () => {
    plotContainer.innerHTML = '';
    const selected: DG.Column[] = colsInput.value ?? [];
    if (selected.length === 0) return;
    const colNames = selected.map((c: DG.Column) => c.name);
    // Clone + run selected method for preview
    const cloneDf = df.clone();
    const method = methodInput.value;
    if (method === 'Median Centering')
      medianNormalize(cloneDf, colNames);
    else if (method === 'Quantile')
      quantileNormalize(cloneDf, colNames);
    // VSN: show original (async not suitable for preview)
    const longDf = unpivotIntensities(cloneDf, colNames);
    const boxPlot = DG.Viewer.boxPlot(longDf, {
      valueColumnName: 'Intensity', categoryColumnName: 'Sample',
    } as any);
    boxPlot.root.style.cssText = 'width:100%;height:220px;';
    plotContainer.appendChild(boxPlot.root);
  };
  methodInput.onChanged.subscribe(updatePreview);
  colsInput.onChanged.subscribe(updatePreview);

  ui.dialog('Normalize')
    .add(warningDiv)
    .add(methodInput)
    .add(colsInput)
    .add(plotContainer)
    .onOK(async () => {
      const selected = colsInput.value.map((c: DG.Column) => c.name);
      const method = methodInput.value;
      if (method === 'Median Centering') medianNormalize(df, selected);
      else if (method === 'Quantile') quantileNormalize(df, selected);
      else if (method === 'VSN') await vsnNormalize(df, selected);
      grok.shell.info(`Normalized ${selected.length} columns (${method})`);
    })
    .show();

  updatePreview(); // initial render
}
```

### Imputation Dialog with Conditional Parameters and Valid-Values Filter
```typescript
// Verified patterns from: differential-expression.ts:275-279, enrichment.ts:344-355
const methodInput = ui.input.choice('Method', {
  value: 'MinProb', items: ['MinProb', 'kNN', 'Zero', 'Mean', 'Median'],
  nullable: false,
});

// MinProb params (created upfront, visibility toggled)
const downshiftInput = ui.input.float('Downshift', {value: 1.8});
const widthInput = ui.input.float('Width', {value: 0.3});
// kNN params
const kInput = ui.input.int('k (neighbors)', {value: 10});

// Toggle visibility
const minProbContainer = ui.div([downshiftInput.root, widthInput.root]);
const knnContainer = ui.div([kInput.root]);
knnContainer.style.display = 'none';

methodInput.onChanged.subscribe(() => {
  minProbContainer.style.display = methodInput.value === 'MinProb' ? '' : 'none';
  knnContainer.style.display = methodInput.value === 'kNN' ? '' : 'none';
});

// IMP-04: Valid-values filter with live count
const minValidInput = ui.input.int('Min valid values per group', {value: 2});
const filterCountDiv = ui.divText('');
const updateFilterCount = () => {
  const threshold = minValidInput.value ?? 0;
  const groups = getGroups(df);
  if (!groups || threshold === 0) {
    filterCountDiv.textContent = `Will keep all ${df.rowCount} proteins`;
    return;
  }
  let kept = 0;
  for (let i = 0; i < df.rowCount; i++) {
    // Protein passes if ANY group has >= threshold valid values
    const g1Valid = groups.group1.columns.filter((n) => {
      const c = df.col(n); return c && !c.isNone(i);
    }).length;
    const g2Valid = groups.group2.columns.filter((n) => {
      const c = df.col(n); return c && !c.isNone(i);
    }).length;
    if (g1Valid >= threshold || g2Valid >= threshold) kept++;
  }
  const removed = df.rowCount - kept;
  filterCountDiv.textContent = `Will keep ${kept}/${df.rowCount} proteins (${removed} removed)`;
};
updateFilterCount();
minValidInput.onChanged.subscribe(updateFilterCount);
```

### DE Comparison Direction with Hint Text
```typescript
// Verified patterns from: experiment-setup.ts:19-27
const groups = getGroups(df)!;
const g1 = groups.group1;
const g2 = groups.group2;

const pairs = [
  `${g2.name} vs ${g1.name}`,
  `${g1.name} vs ${g2.name}`,
];
const comparisonInput = ui.input.choice('Comparison', {
  value: pairs[0], items: pairs, nullable: false,
});

const hintDiv = ui.divText('');
hintDiv.style.cssText = 'font-style:italic;color:#888;font-size:12px;margin-bottom:8px;';
const updateHint = () => {
  const parts = comparisonInput.value!.split(' vs ');
  hintDiv.textContent = `Positive log2FC = higher in ${parts[0]}, ` +
    `Negative log2FC = higher in ${parts[1]}`;
};
updateHint();
comparisonInput.onChanged.subscribe(updateHint);
```

### Viewer Title Setting
```typescript
// Verified pattern from: qc-dashboard.ts:62-66
// Volcano:
const sp = createVolcanoPlot(df, options);
sp.setOptions({title: `Volcano: ${g2Name} vs ${g1Name}`});

// PCA:
const {viewer: pcaSp, pcaDf} = createPcaPlot(df, allCols, groups);
pcaSp.setOptions({title: 'PCA: All Groups'});

// Heatmap:
const grid = await createExpressionHeatmap(df, {topN: 50});
grid.setOptions({title: `Heatmap: Top 50 DE Proteins`});
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Single-method dialogs (median only, MinProb only) | Method selector with conditional params | Phase 11 | Users choose from all available methods |
| No distribution preview in normalization | Reactive box plot in dialog | Phase 11 | Visual confirmation of normalization effect |
| No comparison direction picker in DE | Dropdown with directional pairs + hint | Phase 11 | Users control FC sign interpretation |
| Untitled viewers | "Viewer: Context" format titles | Phase 11 | Users understand what each viewer shows |
| Hardcoded df.name in parsers | Filename-based naming in import handlers | Phase 11 | Users see meaningful table names |

**Current state of existing dialogs (before Phase 11):**
- Normalization (`normalization.ts:156-180`): Only median centering, column picker, no method selector, no preview, no warning
- Imputation (`imputation.ts:242-272`): Only MinProb with downshift/width params, no method selector, no valid-values filter
- DE (`differential-expression.ts:235-360`): Has limma/DEqMS method selector with conditional peptide column visibility, but no comparison direction picker, no t-test in dropdown

## Open Questions

1. **Box Plot Reactive Update Strategy**
   - What we know: `DG.Viewer.boxPlot()` creates a viewer bound to a specific DataFrame
   - What's unclear: Whether replacing the underlying DataFrame data triggers automatic re-render, or if a new viewer must be created each time
   - Recommendation: Create a new cloned DataFrame and new box plot viewer on each method change. Clear the container div before adding new viewer root. This is simpler and avoids potential viewer lifecycle issues.

2. **VSN Preview Feasibility**
   - What we know: VSN requires async R server call; median and quantile are synchronous
   - What's unclear: Acceptable latency for in-dialog preview
   - Recommendation: For VSN, show the current (un-normalized) distributions with a note, or skip the clone-and-normalize step. The CONTEXT.md says "before-state only" which simplifies this -- the box plot can show current distributions regardless of method.

3. **Grid Title Support**
   - What we know: Grid extends Viewer, which has `setOptions()`. QC dashboard uses `title` in ScatterPlot options.
   - What's unclear: Whether Grid's `setOptions({title: ...})` actually renders a visible title
   - Recommendation: Try it first. If Grid does not support `title`, add a `ui.h3()` element above the grid root as fallback.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Datagrok Puppeteer-based test framework |
| Config file | packages/Proteomics/webpack.config.js (test entry point) |
| Quick run command | `grok test --test "Dialog" --host localhost` |
| Full suite command | `grok test --host localhost` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| NORM-03 | Normalization dialog shows method selector and box plot preview | manual-only | N/A -- requires visual dialog interaction | N/A |
| NORM-04 | Warning banner appears for Spectronaut-tagged data | unit | Test `proteomics.preNormalized` tag detection logic | Wave 0 |
| IMP-03 | Imputation dialog shows conditional parameters per method | manual-only | N/A -- requires DOM interaction | N/A |
| IMP-04 | Valid values filter removes proteins below threshold | unit | Test filter counting logic with mock DataFrame | Wave 0 |
| DE-01 | Comparison pairs auto-generated from groups | unit | Test pair generation from GroupAssignment | Wave 0 |
| DE-02 | Method-specific params shown/hidden | manual-only | N/A -- requires DOM interaction | N/A |
| UX-01 | Viewer titles include descriptive context | smoke | Verify viewer creation includes `title` option | Wave 0 |
| UX-02 | DataFrame name matches filename | unit | Test all parser output df.name values | Wave 0 |

### Sampling Rate
- **Per task commit:** `grok test --test "Norm" --host localhost` (or relevant subset)
- **Per wave merge:** `grok test --host localhost`
- **Phase gate:** Full suite green + manual dialog walkthrough

### Wave 0 Gaps
- [ ] Unit test for valid-values filter counting logic (IMP-04)
- [ ] Unit test for comparison pair generation from GroupAssignment (DE-01)
- [ ] Unit test for df.name derivation from filename across all parsers (UX-02)
- [ ] Smoke test for viewer title option passthrough (UX-01)

**Note:** Dialog interaction tests (NORM-03, IMP-03, DE-02) are manual-only because they require DOM rendering in a running Datagrok instance. The underlying logic (normalization methods, imputation methods, DE computation) is already tested through existing pipeline tests.

### Observable Behaviors for Manual Validation
1. **Normalization dialog:** Open dialog -> see method dropdown (Median Centering/Quantile/VSN) -> box plot shows current distributions -> change method -> box plot updates -> for Spectronaut data, yellow/orange banner appears -> click OK -> data normalizes correctly
2. **Imputation dialog:** Open dialog -> see method dropdown (MinProb/kNN/Zero/Mean/Median) -> select MinProb -> downshift/width inputs visible -> select kNN -> k input visible, MinProb params hidden -> select Zero -> no extra params -> adjust min valid values -> live count updates ("Will keep X/Y proteins (Z removed)") -> click OK -> proteins filtered then imputed
3. **DE dialog:** Open dialog -> comparison dropdown at top ("Treatment vs Control", "Control vs Treatment") -> hint text shows FC direction interpretation -> change selection -> hint updates -> select DEqMS -> peptide column input appears -> select limma -> peptide column hides -> click OK -> DE runs with correct direction
4. **Viewer titles:** Open volcano -> title reads "Volcano: Treatment vs Control" -> open PCA -> title reads "PCA: All Groups" -> open heatmap -> title reads "Heatmap: Top 50 DE Proteins"
5. **DataFrame naming:** Import proteinGroups.txt via MaxQuant -> table named "proteinGroups" -> Import HYE_mix.tsv via Spectronaut -> table named "HYE_mix" -> Import custom.csv via Generic -> table named "custom"

### Edge Cases to Validate
- Normalization dialog with only 1 intensity column (quantile requires 2+ -- should warn or disable)
- Imputation dialog with 0 missing values (should work, imputes nothing)
- Valid values filter threshold higher than max samples in any group (removes all proteins -- should show "Will keep 0/Y proteins" as warning)
- VSN selected but R server unavailable (fallback to quantile with warning message)
- Spectronaut data re-imported after normalization (`preNormalized` tag persists)
- Demo function call to `parseMaxQuantText()` (no filename available -- needs default name)
- Auxiliary DataFrame naming: PCA, Enrichment Results, Heatmap DataFrames include source table name

## Sources

### Primary (HIGH confidence)
- Existing codebase: `differential-expression.ts` (conditional visibility pattern, lines 255-279)
- Existing codebase: `enrichment.ts` (live count update pattern, lines 344-355)
- Existing codebase: `generic-parser.ts` (viewer-in-dialog pattern, lines 86-115)
- Existing codebase: `qc-dashboard.ts` (box plot factory + viewer titles, lines 62-142)
- Existing codebase: `qc-computations.ts` (`unpivotIntensities()` function)
- Datagrok js-api: `viewer.ts:226` (`Viewer.boxPlot()` factory)
- Datagrok js-api: `interfaces/d4.ts:252` (`IBoxPlotSettings` interface)
- Datagrok js-api: `const.ts:666` (`VIEWER.BOX_PLOT`)

### Secondary (MEDIUM confidence)
- QC dashboard viewer title pattern: confirmed working via `title: 'MA Plot'` in ScatterPlot options
- Grid `setOptions({title: ...})` -- Grid extends Viewer but title rendering not verified for Grid specifically

### Tertiary (LOW confidence)
- Box plot viewer refresh behavior when underlying DataFrame values change (needs runtime verification -- recommended to create new viewer on each update)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all Datagrok API, no external libraries needed
- Architecture: HIGH - every pattern proven in existing codebase with exact line references
- Pitfalls: HIGH - based on direct code analysis of existing implementations and API interfaces
- Viewer titles for Grid: MEDIUM - ScatterPlot title confirmed, Grid title needs runtime test

**Research date:** 2026-03-07
**Valid until:** 2026-04-07 (stable -- Datagrok API changes slowly, all patterns from current codebase)
