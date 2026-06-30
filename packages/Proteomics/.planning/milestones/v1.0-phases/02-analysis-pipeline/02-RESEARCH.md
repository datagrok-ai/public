# Phase 2: Analysis Pipeline - Research

**Researched:** 2026-02-28
**Domain:** Proteomics analysis pipeline: experiment annotation, normalization, imputation, differential expression
**Confidence:** HIGH

## Summary

Phase 2 implements the core analysis pipeline for proteomics data within the existing Datagrok package. The pipeline has four stages: (1) experiment annotation -- assigning sample columns to experimental groups via a dialog, persisting as DataFrame tags; (2) median centering normalization -- subtracting per-column median from log2 intensity values; (3) MinProb missing value imputation -- replacing NaN with random draws from a downshifted normal distribution; (4) differential expression -- computing per-protein fold changes and p-values between two groups.

The critical architectural decision for this phase concerns differential expression. The REQUIREMENTS.md specifies "via R script" for ANLY-03 (limma), but the phase context states "All computation happens client-side in TypeScript/JavaScript (no R/Python server)." The project STATE.md records the roadmap decision as "R for DE (limma/DEqMS), TypeScript for normalization/imputation." Since limma has no JavaScript implementation and requires R's empirical Bayes infrastructure, a pure client-side implementation must use Welch's t-test + Benjamini-Hochberg FDR correction instead. The `@datagrok-libraries/statistics` package already provides both `tTest()` (using jStat) and `fdrcorrection()` (BH/BY methods). This is a scientifically sound approach: for proteomics with 3+ replicates per group, Welch's t-test with BH correction produces valid differential expression results. Limma's advantage (variance shrinkage for small sample sizes) can be added later as an R-script enhancement.

**Primary recommendation:** Implement all four pipeline stages in pure TypeScript using Datagrok's DataFrame API and `@datagrok-libraries/statistics`. Store group assignments as DataFrame tags. Use `ui.dialog()` with `ui.input.columns()` for all user-facing dialogs. Apply Welch's t-test per protein with BH-corrected p-values for differential expression.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| SETUP-01 | User can assign sample columns to experimental groups via annotation dialog | `ui.dialog()` with two `ui.input.columns()` inputs for Group A and Group B; filter to show only log2 intensity columns |
| SETUP-02 | Group assignments persist as DataFrame metadata for downstream analysis | `df.setTag('proteomics.groups', JSON.stringify({group1: [...names], group2: [...names]}))` -- tags survive serialization |
| ANLY-01 | User can normalize intensity columns using median centering | Pure TypeScript: for each column, compute median via `col.stats.med`, subtract from each value |
| ANLY-02 | User can impute missing values using MinProb with configurable parameters | Pure TypeScript: per-column, compute mean/stdev of non-missing values, draw from N(mean - downshift*stdev, width*stdev) |
| ANLY-03 | User can run differential expression between two groups | Client-side Welch's t-test via `@datagrok-libraries/statistics` `tTest()` + `fdrcorrection()` for BH adjustment; log2FC as mean difference |
| ANLY-05 | DE results include log2FC, p-value, and adjusted p-value columns with correct semantic types | Add columns `log2FC` (semType `Proteomics-Log2FC`), `p-value` (semType `Proteomics-PValue`), `adj.p-value` (semType `Proteomics-PValue`) |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| datagrok-api | ^1.25.0 | DataFrame API, UI dialogs, column inputs, tags | Platform API -- the only way to interact with Datagrok |
| @datagrok-libraries/statistics | ^1.2.12 | `tTest()` for Welch's t-test, `fdrcorrection()` for BH p-value adjustment | Already in package.json; provides validated statistical implementations with jStat |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| @datagrok-libraries/test | ^1.1.0 | `category`, `test`, `expect` for unit tests | Test files in `src/tests/` |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Client-side Welch's t-test | R limma via `grok.functions.call()` | Limma provides variance shrinkage (better with n<3 replicates); requires R server infrastructure, conda environment setup, and server-side execution. Add as v2 enhancement |
| `@datagrok-libraries/statistics` tTest | Manual t-test implementation | Library is already validated and handles edge cases (jStat for t-distribution CDF) |
| DataFrame tags for group storage | DataFrame `temp` property | Tags persist across serialization/save; `temp` is ephemeral and lost on close |

**Installation:**
No additional packages needed. All dependencies already in `package.json`.

## Architecture Patterns

### Recommended Project Structure
```
src/
  package.ts              # Wire annotateExperiment(), normalizeProteomics(), imputeMissingValues(), differentialExpression()
  analysis/
    normalization.ts      # medianNormalize() -- EXISTS, needs implementation
    imputation.ts         # imputeMinProb() -- EXISTS, needs implementation
    differential-expression.ts  # runDifferentialExpression() -- EXISTS, needs rewrite (remove DEResult interface, use DataFrame)
    experiment-setup.ts   # NEW: showAnnotationDialog(), getGroups(), setGroups()
  utils/
    proteomics-types.ts   # SEMTYPE constants (exists)
    column-detection.ts   # findColumn, findProteomicsColumns (exists)
  tests/
    parsers.ts            # Parser tests (exists)
    analysis.ts           # NEW: normalization, imputation, DE tests
```

### Pattern 1: Experiment Annotation via DataFrame Tags
**What:** Store experimental group assignments (which columns belong to which group) as serializable DataFrame tags.
**When to use:** After user assigns groups via dialog, before normalization/imputation/DE.
**Example:**
```typescript
// Source: Verified from js-api/src/dataframe/data-frame.ts (setTag/getTag methods)

interface GroupAssignment {
  group1: {name: string; columns: string[]};
  group2: {name: string; columns: string[]};
}

/** Store group assignments as a DataFrame tag. */
function setGroups(df: DG.DataFrame, groups: GroupAssignment): void {
  df.setTag('proteomics.groups', JSON.stringify(groups));
}

/** Retrieve group assignments from DataFrame tag. Returns null if not set. */
function getGroups(df: DG.DataFrame): GroupAssignment | null {
  const tag = df.getTag('proteomics.groups');
  return tag ? JSON.parse(tag) : null;
}
```

### Pattern 2: Annotation Dialog with Column Selection
**What:** Dialog with two column-list inputs for assigning samples to groups.
**When to use:** The annotateExperiment() menu handler.
**Example:**
```typescript
// Source: Verified from js-api/ui.ts (ui.input.columns, ui.input.string, ui.dialog)

function showAnnotationDialog(df: DG.DataFrame): void {
  const intensityCols = df.columns.toList()
    .filter((c) => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('));

  const group1Input = ui.input.columns('Group 1', {table: df, available: intensityCols});
  const group1Name = ui.input.string('Group 1 Name', {value: 'Control'});
  const group2Input = ui.input.columns('Group 2', {table: df, available: intensityCols});
  const group2Name = ui.input.string('Group 2 Name', {value: 'Treatment'});

  ui.dialog('Annotate Experiment')
    .add(group1Name)
    .add(group1Input)
    .add(group2Name)
    .add(group2Input)
    .onOK(() => {
      const groups: GroupAssignment = {
        group1: {name: group1Name.value, columns: group1Input.value.map((c: DG.Column) => c.name)},
        group2: {name: group2Name.value, columns: group2Input.value.map((c: DG.Column) => c.name)},
      };
      setGroups(df, groups);
      grok.shell.info(`Groups assigned: ${groups.group1.columns.length} + ${groups.group2.columns.length} samples`);
    })
    .show();
}
```

### Pattern 3: Analysis Function Modifying DataFrame In-Place
**What:** Analysis functions (normalize, impute, DE) modify the current DataFrame by adding or modifying columns, rather than creating new DataFrames.
**When to use:** All analysis steps -- normalization modifies values in existing columns, imputation fills NaN in existing columns, DE adds new columns.
**Example:**
```typescript
// Normalization modifies columns in-place
function medianNormalize(df: DG.DataFrame, colNames: string[]): void {
  for (const name of colNames) {
    const col = df.col(name)!;
    const median = col.stats.med;
    for (let i = 0; i < df.rowCount; i++) {
      if (!col.isNone(i))
        col.set(i, col.get(i) - median);
    }
  }
}

// DE adds new columns
function runDE(df: DG.DataFrame, group1Cols: string[], group2Cols: string[]): void {
  const log2fcCol = df.columns.addNewFloat('log2FC');
  const pValCol = df.columns.addNewFloat('p-value');
  // ... compute per protein ...
  log2fcCol.semType = SEMTYPE.LOG2FC;
  pValCol.semType = SEMTYPE.P_VALUE;
}
```

### Pattern 4: Dialog Reads Groups from Tags
**What:** Analysis dialogs (normalize, impute, DE) auto-populate from persisted group assignments.
**When to use:** When the user has already annotated the experiment.
**Example:**
```typescript
function showDEDialog(df: DG.DataFrame): void {
  const groups = getGroups(df);
  if (!groups) {
    grok.shell.warning('Please annotate experimental groups first (Proteomics | Annotate Experiment)');
    return;
  }
  // Pre-populate dialog with group info from tags
  // ...
}
```

### Anti-Patterns to Avoid
- **Creating new DataFrames for each analysis step:** Modify the current DataFrame in-place (add columns, modify values). Creating new DataFrames loses the user's view state, filters, and selections.
- **Storing groups in a separate data structure:** Use DataFrame tags so they persist with the table across save/load.
- **Implementing custom t-distribution CDF:** Use `@datagrok-libraries/statistics` which wraps jStat for validated implementations.
- **Hand-rolling FDR correction:** Use `fdrcorrection()` from `@datagrok-libraries/statistics`.
- **Using the DEResult interface:** The existing `differential-expression.ts` defines a `DEResult` interface. Delete it -- results should be added as DataFrame columns directly, matching the Datagrok pattern.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Welch's t-test | Manual t-statistic + p-value calculation | `tTest()` from `@datagrok-libraries/statistics` | Handles degrees of freedom, t-distribution CDF via jStat correctly |
| FDR correction (BH) | Manual p-value ranking and adjustment | `fdrcorrection()` from `@datagrok-libraries/statistics` | Validated implementation of Benjamini-Hochberg and Benjamini-Yekutieli |
| Column median | Manual sorted array + midpoint | `col.stats.med` from Datagrok Stats API | Built-in, handles nulls, uses optimized Dart implementation |
| Column mean/stdev | Manual sum/variance loop | `col.stats.avg` and `col.stats.stdev` | Built-in Stats class provides all descriptive statistics |
| Random normal draws | Custom Box-Muller transform | `Math.sqrt(-2*Math.log(u1))*Math.cos(2*Math.PI*u2)` | Box-Muller is simple enough for MinProb, but if needed, jStat provides `jStat.normal.sample()` via statistics lib |
| Dialog with column selection | Custom HTML form | `ui.dialog()` + `ui.input.columns()` | Standard Datagrok pattern, handles table binding and column filtering |

**Key insight:** The combination of Datagrok's Stats API (median, mean, stdev per column) and `@datagrok-libraries/statistics` (t-test, FDR correction) covers the entire statistical foundation needed. The analysis code is mostly orchestration -- iterate proteins, extract group values, call statistical functions, write results to columns.

## Common Pitfalls

### Pitfall 1: Operating on Raw Instead of Log2 Columns
**What goes wrong:** Normalization or DE applied to raw intensity columns instead of log2-transformed ones produces incorrect results.
**Why it happens:** The DataFrame contains both raw (`LFQ intensity Sample1`) and log2 (`log2(LFQ intensity Sample1)`) columns. Analysis must operate on log2 columns.
**How to avoid:** Filter columns by checking for `log2(` prefix in name, or filter by both `SEMTYPE.INTENSITY` and name pattern. The annotation dialog should only show log2 columns.
**Warning signs:** Normalization shifts of thousands instead of single digits; non-symmetric fold change distributions.

### Pitfall 2: NaN Handling in Statistical Functions
**What goes wrong:** `tTest()` receives arrays containing `NaN`/null values and produces `NaN` p-values.
**Why it happens:** Datagrok columns use `DG.FLOAT_NULL` for missing values; when extracted to arrays, these become special float values. `col.isNone(i)` must be checked.
**How to avoid:** When extracting column values for t-test, explicitly filter out null values: `const values = []; for (let i = 0; i < col.length; i++) if (!col.isNone(i)) values.push(col.get(i));`
**Warning signs:** All p-values are NaN; mean differences are NaN.

### Pitfall 3: Insufficient Replicates for T-Test
**What goes wrong:** T-test requires at least 2 values per group. Proteins with only 0-1 valid values in a group cannot be tested.
**Why it happens:** Missing values are common in proteomics (30-50% missingness). After filtering nulls, some proteins may have too few values.
**How to avoid:** Before running t-test per protein, check that both groups have >= 2 non-null values. For proteins below threshold, set p-value to NaN (cannot test). Report the count of untestable proteins.
**Warning signs:** Error from `tTest()` about wrong sample size; large fraction of NaN p-values.

### Pitfall 4: MinProb Imputation Width Parameter
**What goes wrong:** Imputed values are unrealistically tight (all identical) or unrealistically spread.
**Why it happens:** The `width` parameter controls the standard deviation of the imputation distribution relative to the observed standard deviation. Default 0.3 means 30% of observed SD.
**How to avoid:** Use established defaults: `downshift = 1.8`, `width = 0.3`. These are the Perseus defaults used by the proteomics community. Allow user to adjust via dialog but clearly label defaults.
**Warning signs:** Imputed value distribution is either a spike (width too small) or overlaps with real values (downshift too small).

### Pitfall 5: fdrcorrection Expects Float32Array
**What goes wrong:** `fdrcorrection()` from `@datagrok-libraries/statistics` expects `Float32Array`, not regular arrays.
**Why it happens:** The function signature requires `Float32Array` for efficient typed array operations.
**How to avoid:** Convert p-value arrays to `Float32Array` before calling: `new Float32Array(pValues)`. Also handle proteins where p-value is NaN by excluding them from FDR correction and setting their adjusted p-value to NaN.
**Warning signs:** Type errors at runtime; incorrect FDR-adjusted p-values.

### Pitfall 6: Tag Serialization of Column Names
**What goes wrong:** Column names containing special characters (parentheses in `log2(LFQ intensity Sample1)`) break JSON serialization or retrieval.
**Why it happens:** Column names are stored as JSON strings in tags; but later lookup via `df.col(name)` must use the exact name.
**How to avoid:** Store column names as-is in JSON (they serialize fine). When retrieving, validate that stored column names still exist in the DataFrame: `const col = df.col(name); if (!col) throw ...`.
**Warning signs:** Group column lookup returns null after save/load.

### Pitfall 7: In-Place Column Modification Triggers Viewer Updates
**What goes wrong:** Modifying column values (normalization) causes visual glitches if viewers are actively rendering.
**Why it happens:** Datagrok viewers observe column changes in real-time.
**How to avoid:** Call `col.set(i, value)` in a tight loop -- Datagrok batches updates. After all modifications, call `df.fireValuesChanged()` once to notify all viewers of the bulk change.
**Warning signs:** Grid flickering during normalization; viewers showing intermediate states.

## Code Examples

### Median Centering Normalization
```typescript
// Source: Verified against Datagrok Stats API (js-api/src/dataframe/stats.ts)

import * as DG from 'datagrok-api/dg';

/** Centers each log2 intensity column by subtracting its median.
 * Operates in-place on the DataFrame columns. */
export function medianNormalize(df: DG.DataFrame, colNames: string[]): void {
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;

    // col.stats.med computes median of non-null values
    const median = col.stats.med;
    if (isNaN(median)) continue; // skip all-null columns

    for (let i = 0; i < df.rowCount; i++) {
      if (!col.isNone(i))
        col.set(i, (col.get(i) as number) - median);
    }
  }
  df.fireValuesChanged();
}
```

### MinProb Imputation
```typescript
// Source: Perseus algorithm (Tyanova et al., 2016)

import * as DG from 'datagrok-api/dg';

/** Box-Muller transform for generating normal random variates. */
function randomNormal(mean: number, stdev: number): number {
  const u1 = Math.random();
  const u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return mean + stdev * z;
}

/** Imputes missing values using the MinProb method.
 * For each column, draws random values from N(mean - downshift*sd, width*sd).
 * Operates in-place on the DataFrame columns. */
export function imputeMinProb(
  df: DG.DataFrame,
  colNames: string[],
  downshift: number = 1.8,
  width: number = 0.3,
): void {
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;

    const stats = col.stats;
    const mean = stats.avg;
    const sd = stats.stdev;

    if (isNaN(mean) || isNaN(sd) || sd === 0) continue;

    const imputeMean = mean - downshift * sd;
    const imputeSd = width * sd;

    for (let i = 0; i < df.rowCount; i++) {
      if (col.isNone(i))
        col.set(i, randomNormal(imputeMean, imputeSd));
    }
  }
  df.fireValuesChanged();
}
```

### Welch's T-Test Differential Expression
```typescript
// Source: @datagrok-libraries/statistics (tests.ts tTest, multiple-tests.ts fdrcorrection)

import * as DG from 'datagrok-api/dg';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {fdrcorrection} from '@datagrok-libraries/statistics/src/multiple-tests';
import {SEMTYPE} from '../utils/proteomics-types';

/** Extracts non-null float values from specified columns for a given row. */
function getGroupValues(df: DG.DataFrame, colNames: string[], rowIdx: number): number[] {
  // NOTE: For proteomics DE, each ROW is a protein, each COLUMN is a sample.
  // But tTest expects arrays of values across samples for a single protein.
  // So we need to extract values across columns for a single row.
  const values: number[] = [];
  for (const name of colNames) {
    const col = df.col(name);
    if (col && !col.isNone(rowIdx))
      values.push(col.get(rowIdx) as number);
  }
  return values;
}

/** Runs Welch's t-test differential expression between two groups.
 * Adds log2FC, p-value, and adj.p-value columns to the DataFrame. */
export function runDifferentialExpression(
  df: DG.DataFrame,
  group1Cols: string[],
  group2Cols: string[],
): void {
  const n = df.rowCount;
  const log2fcCol = df.columns.addNewFloat('log2FC');
  const pValCol = df.columns.addNewFloat('p-value');

  const pValues: number[] = [];
  const testableIndices: number[] = [];

  for (let i = 0; i < n; i++) {
    const vals1 = getGroupValues(df, group1Cols, i);
    const vals2 = getGroupValues(df, group2Cols, i);

    if (vals1.length < 2 || vals2.length < 2) {
      // Cannot test -- insufficient data
      log2fcCol.set(i, DG.FLOAT_NULL);
      pValCol.set(i, DG.FLOAT_NULL);
      continue;
    }

    // log2FC = mean(group2) - mean(group1)
    const mean1 = vals1.reduce((a, b) => a + b, 0) / vals1.length;
    const mean2 = vals2.reduce((a, b) => a + b, 0) / vals2.length;
    log2fcCol.set(i, mean2 - mean1);

    // Welch's t-test
    const result = tTest(vals1, vals2);
    const pVal = result['p-value'];
    pValCol.set(i, pVal);

    pValues.push(pVal);
    testableIndices.push(i);
  }

  // BH FDR correction on testable proteins only
  const adjPValCol = df.columns.addNewFloat('adj.p-value');
  if (pValues.length > 0) {
    const [, corrected] = fdrcorrection(
      new Float32Array(pValues), 0.05, 'i', // 'i' = Benjamini-Hochberg
    );
    for (let j = 0; j < testableIndices.length; j++)
      adjPValCol.set(testableIndices[j], corrected[j]);
  }

  // Set untestable proteins to null
  for (let i = 0; i < n; i++) {
    if (pValCol.isNone(i))
      adjPValCol.set(i, DG.FLOAT_NULL);
  }

  // Assign semantic types
  log2fcCol.semType = SEMTYPE.LOG2FC;
  pValCol.semType = SEMTYPE.P_VALUE;
  adjPValCol.semType = SEMTYPE.P_VALUE;

  df.fireValuesChanged();
}
```

### Normalization Dialog
```typescript
// Source: Verified from js-api/ui.ts (ui.dialog, ui.input.columns, ui.input.choice)

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';
import {medianNormalize} from '../analysis/normalization';

function showNormalizationDialog(df: DG.DataFrame): void {
  const log2Cols = df.columns.toList()
    .filter((c) => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('));

  const colsInput = ui.input.columns('Intensity columns', {
    table: df,
    available: log2Cols,
    value: log2Cols,
  });

  ui.dialog('Normalize Proteomics Data')
    .add(colsInput)
    .onOK(() => {
      const colNames = colsInput.value.map((c: DG.Column) => c.name);
      medianNormalize(df, colNames);
      grok.shell.info(`Normalized ${colNames.length} columns (median centering)`);
    })
    .show();
}
```

### Getting Intensity Column Names from Groups
```typescript
// Helper to get log2 intensity column names that match group assignments
function getGroupColumnNames(df: DG.DataFrame): {group1Cols: string[]; group2Cols: string[]} | null {
  const groups = getGroups(df);
  if (!groups) return null;

  // Validate columns still exist
  const group1Cols = groups.group1.columns.filter((n) => df.col(n) !== null);
  const group2Cols = groups.group2.columns.filter((n) => df.col(n) !== null);

  if (group1Cols.length === 0 || group2Cols.length === 0) return null;
  return {group1Cols, group2Cols};
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| DEResult interface for DE output | Add columns directly to DataFrame | Phase 2 decision | Delete DEResult interface; columns are the standard Datagrok pattern |
| R limma via grok.functions.call | Client-side Welch's t-test + BH correction | Phase 2 constraint (client-side only) | Simpler deployment, no R server needed; add R limma as optional enhancement later |
| Manual median/stdev calculation | `col.stats.med` and `col.stats.stdev` | Always available in Datagrok | Use built-in Stats API, no manual computation needed |

**Deprecated/outdated:**
- The `DEResult` interface in `differential-expression.ts` must be deleted. DE results are columns in the DataFrame, not a separate data structure.
- The `runLimmaDE` function signature taking `DG.Column[]` arrays should change to take column name strings for simpler serialization and tag storage.

## Open Questions

1. **R Script Support Conflict**
   - What we know: REQUIREMENTS.md says ANLY-03 uses "R script" for limma. The phase context says "all computation client-side."
   - What's unclear: Whether the user wants to defer R/limma entirely or implement a client-side t-test fallback now with R later.
   - Recommendation: Implement client-side Welch's t-test + BH correction for this phase. Flag the R limma script (`scripts/limma_de.R`) as a future enhancement. The existing R script stub can remain for later implementation.

2. **Column Modification vs New Columns for Normalization**
   - What we know: The product spec mentions `_norm` suffix columns or in-place replacement.
   - What's unclear: Whether normalization should modify log2 columns in-place or create new `_norm` columns.
   - Recommendation: Modify log2 columns in-place. Proteomics normalization is always applied (not optional to undo), and creating duplicate columns doubles the already-wide table. The product spec also mentions in-place as an option.


3. **Imputation: Per-Sample vs Per-Column**
   - What we know: MinProb typically operates per-sample (per-column in protein x sample matrix).
   - What's unclear: Whether to impute per-column (each sample independently) or use a global distribution.
   - Recommendation: Per-column (per-sample) imputation, matching Perseus behavior and the ANLY-02 requirement "(per-sample)".

## Sources

### Primary (HIGH confidence)
- Datagrok js-api source: `js-api/src/dataframe/data-frame.ts` (DataFrame.setTag/getTag, tags property, fireValuesChanged)
- Datagrok js-api source: `js-api/src/dataframe/stats.ts` (Stats class: med, avg, stdev, min, max)
- Datagrok js-api source: `js-api/ui.ts` (ui.dialog, ui.input.columns, ui.input.string, ui.input.choice)
- `@datagrok-libraries/statistics` source: `libraries/statistics/src/tests.ts` (tTest with Welch's implementation using jStat)
- `@datagrok-libraries/statistics` source: `libraries/statistics/src/multiple-tests.ts` (fdrcorrection with BH and BY methods)
- Existing Proteomics scaffold: `packages/Proteomics/src/analysis/` (stub files with TODO comments)
- Existing Proteomics parser: `packages/Proteomics/src/parsers/maxquant-parser.ts` (verified working pattern)
- Datagrok R script pattern: `packages/Samples/scripts/r/pca.R` (metadata comments, column_list input, action:join output)

### Secondary (MEDIUM confidence)
- [limma powers differential expression analyses](https://pmc.ncbi.nlm.nih.gov/articles/PMC4402510/) - Limma methodology reference confirming R-only implementation
- [Perseus proteomics platform](https://maxquant.net/perseus/) - MinProb imputation defaults (downshift=1.8, width=0.3)
- [Proteomics Data Analysis in R/Bioconductor](https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/DEA.html) - Standard DE workflow confirmation

### Tertiary (LOW confidence)
- No JavaScript/TypeScript implementation of limma's empirical Bayes variance shrinkage exists. This was confirmed via web search but should be validated against npm registry.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries already in package.json; APIs verified from source code
- Architecture: HIGH - Patterns follow established Datagrok conventions verified from existing packages
- Pitfalls: HIGH - Common proteomics data issues well-documented; API behavior verified from source
- DE approach: MEDIUM - Welch's t-test is scientifically valid but user may prefer limma; flagged as open question

**Research date:** 2026-02-28
**Valid until:** 2026-03-28 (stable domain, no fast-moving dependencies)
