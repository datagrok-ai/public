# Phase 10: Spectronaut Parser and Core Algorithms - Research

**Researched:** 2026-03-07
**Domain:** Spectronaut long-format TSV parsing, quantile normalization, VSN normalization, kNN/simple imputation
**Confidence:** HIGH

## Summary

Phase 10 adds a Spectronaut parser (long-to-wide pivot) and four new algorithmic functions (quantile normalization, VSN normalization, kNN imputation, zero/mean/median imputation) to the existing Proteomics package. The codebase has well-established patterns for all three concerns: the MaxQuant parser provides the filter-clone-semtype-transform template, `normalization.ts` provides the in-place-normalize-and-tag template, and `imputation.ts` provides the iterate-columns-and-tag template.

The core technical challenge is the Spectronaut pivot: converting peptide-level long-format rows (8760 rows in demo = 93 proteins x ~94 peptides/run x 8 runs) into a wide protein-by-sample matrix using PG.IBAQ as the quantity column. PG.IBAQ is a protein-level value that is constant across all peptide rows for the same protein+sample, so the pivot is a deduplication (not an aggregation). The demo file has 2 conditions ("HYE mix A", "HYE mix B") x 4 replicates = 8 samples, producing an output DataFrame of 93 rows x 8 intensity columns.

The normalization and imputation algorithms are well-defined mathematically and straightforward to implement in TypeScript. VSN is the only one requiring R (via the `vsn` Bioconductor package). Zero new npm dependencies are needed per the project decision.

**Primary recommendation:** Mirror the MaxQuant parser pattern exactly for Spectronaut (parse text, filter, pivot, clone, assign semtypes, transform intensities). Implement normalization and imputation functions as pure exported functions following existing module patterns. Register all in package.ts via decorator pattern.

<user_constraints>

## User Constraints (from CONTEXT.md)

### Locked Decisions
- Use PG.IBAQ as the quantity column for protein-level intensities (no aggregation needed)
- Name pivoted sample columns as `R.Condition + "_" + R.Replicate` (e.g., "HYE mix A_1")
- Store R.FileName as a column tag on each pivoted intensity column for traceability
- Q-value filtering: configurable threshold parameter with default 0.01; non-numeric values ('Profiled', 'NaN') treated as passing
- Contaminant/decoy filtering: filter by PG.ProteinGroups prefix (CON__, REV__), matching MaxQuant parser pattern
- Detect pre-normalization using existing `detectLog2Status()` heuristic from shared-utils.ts
- If data detected as pre-normalized, skip log2 transform and use `copyAsLog2Columns()` instead
- Tag as `proteomics.preNormalized = 'true'` (new tag, separate from `proteomics.normalized`)
- Retain both raw and log2 intensity columns in output DataFrame (VSN needs raw intensities)
- kNN: configurable k parameter with default k=10, row-wise (protein-wise) imputation, Euclidean distance on shared non-missing values
- When fewer than k neighbors available, use available neighbors; fall back to column mean if zero neighbors
- Progress indicator via `grok.shell.startProgress()` with percentage updates per protein row
- Export pure TypeScript functions: `quantileNormalize()`, `vsnNormalize()`, `imputeKnn()`, `imputeZero()`, `imputeMean()`, `imputeMedian()`
- Register all functions via `//name:` metadata so they appear in Datagrok's function browser
- Existing normalization/imputation dialogs remain unchanged -- Phase 11 adds method selectors
- VSN auto-falls back to quantile normalization if R is unavailable (matches limma/DEqMS fallback pattern)

### Claude's Discretion
- Exact distance weighting scheme for kNN (uniform vs inverse-distance)
- VSN R script parameter choices (lts.quantile, etc.)
- Spectronaut protein ID column extraction and primary column logic

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope

</user_constraints>

<phase_requirements>

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| SPEC-01 | User can import Spectronaut long-format TSV with auto-detection and long-to-wide pivot | Pivot algorithm using Map<protein, Map<sample, value>> with PG.IBAQ column; mirrors MaxQuant parser flow |
| SPEC-02 | Spectronaut parser filters decoy proteins and applies q-value threshold | `filterByIdPrefix()` pattern from maxquant-parser.ts for CON__/REV__; q-value threshold with configurable default 0.01 |
| SPEC-03 | Spectronaut parser extracts R.Condition/R.Replicate for auto-group annotation | Extract unique Condition values, map to column groups, call `setGroups()` from experiment-setup.ts |
| SPEC-04 | Spectronaut data is tagged as potentially pre-normalized for downstream warning | `detectLog2Status()` from shared-utils.ts, new `proteomics.preNormalized` tag |
| NORM-01 | User can normalize with quantile normalization (client-side TypeScript) | Pure TS implementation: sort columns, compute rank means, replace values; follows `medianNormalize()` pattern |
| NORM-02 | User can normalize with VSN (R script, reads raw intensity columns, fallback to quantile) | R script in `scripts/vsn_normalize.R` following `limma_de.R` pattern; TypeScript wrapper with try/catch fallback |
| IMP-01 | User can impute with kNN (client-side TypeScript, k=10 default, progress indicator) | Pure TS kNN with Euclidean distance, `DG.TaskBarProgressIndicator.create()` for progress |
| IMP-02 | User can impute with zero, mean, or median methods | Three trivial per-column imputation functions following `imputeMinProb()` pattern |

</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| datagrok-api | ^1.25.0 | DataFrame, Column, UI, grok.functions | Already in package.json; all parser/analysis code uses it |
| @datagrok-libraries/statistics | ^1.2.12 | Statistical utilities | Already in package.json; used by DE module |

### Supporting (R scripts)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| vsn (Bioconductor) | 3.72+ | Variance-stabilizing normalization | VSN R script; server-side only |

### No New Dependencies
Per project decision: zero new npm dependencies. Quantile normalization, kNN imputation, and simple imputation are all pure TypeScript. VSN uses an R script following the existing limma/DEqMS pattern.

**Installation:** No changes to package.json needed.

## Architecture Patterns

### Recommended Project Structure
```
src/
  parsers/
    spectronaut-parser.ts    # NEW: Spectronaut long-to-wide parser
    maxquant-parser.ts       # Existing: reference pattern
    shared-utils.ts          # Existing: reusable transform functions
  analysis/
    normalization.ts         # EXTEND: add quantileNormalize(), vsnNormalize()
    imputation.ts            # EXTEND: add imputeKnn(), imputeZero(), imputeMean(), imputeMedian()
  utils/
    proteomics-types.ts      # Existing: SEMTYPE constants
  package.ts                 # EXTEND: add Spectronaut import menu entry, register new functions
scripts/
  vsn_normalize.R            # NEW: VSN R script
```

### Pattern 1: Spectronaut Long-to-Wide Pivot
**What:** Convert peptide-level long-format TSV (one row per peptide per sample) into protein-by-sample wide matrix using PG.IBAQ as quantity.
**When to use:** When Spectronaut TSV is detected (has R.FileName, R.Condition, R.Replicate, PG.ProteinGroups, PG.IBAQ columns).

**Verified demo data structure:**
- 35 columns: R.FileName, R.Condition, R.Replicate, PG.ProteinGroups, PG.Organisms, PG.IBAQ, EG.Qvalue, EG.Identified, plus many fragment-level columns
- 8760 peptide rows -> 93 unique proteins
- 2 conditions: "HYE mix A", "HYE mix B"
- 4 replicates: 1, 2, 3, 4
- 8 samples total
- PG.IBAQ values like 34.8, 212.2 -- protein-level, constant per protein+sample across all peptide rows
- No CON__/REV__ prefixed proteins in demo file (but parser must handle them)
- EG.Qvalue: all numeric in demo file (range 10^-18 to 0.5), but real exports can contain 'Profiled' strings

**Algorithm:**
```
1. Parse TSV into long DataFrame via DG.DataFrame.fromCsv(text, {delimiter: '\t'})
2. Build sample key: `${R.Condition}_${R.Replicate}`
3. Q-value filter: skip rows where numeric EG.Qvalue > threshold; non-numeric treated as passing
4. CON__/REV__ filter on PG.ProteinGroups
5. Build Map<proteinGroup, Map<sampleKey, ibaqValue>> -- first-encountered value wins (all same)
6. Also collect: PG.Organisms per protein, R.FileName per sample
7. Build wide DataFrame using DG.DataFrame.fromColumns()
8. Detect log2 status, transform or copy intensities
9. Assign semtypes (PROTEIN_ID, INTENSITY)
10. Apply addPrimaryColumnIfNeeded for semicolon-delimited protein groups
11. Auto-populate groups from R.Condition using setGroups()
12. Tag: proteomics.source='spectronaut', proteomics.preNormalized if detected
```

**Critical detail:** PG.IBAQ is a protein-level value repeated across all peptide rows for the same protein+sample. The pivot is deduplication, NOT aggregation.

### Pattern 2: MaxQuant Parser Flow (template to mirror)
**What:** The established parser flow from `maxquant-parser.ts`.
**Source:** Full source in `src/parsers/maxquant-parser.ts` (147 lines).

Key flow:
1. `buildIntensityColumnOptions(text)` -- scan headers, force double type
2. `DG.DataFrame.fromCsv(text, {delimiter: '\t', columnImportOptions})` -- parse
3. Filter: `filterByMarker()`, `filterByIdPrefix()` -- set filter bits, `fireChanged()`
4. `raw.clone(raw.filter)` -- clone filtered rows
5. `assignSemanticTypes(df)` -- PROTEIN_ID, GENE_SYMBOL
6. `addPrimaryColumnFromVariants()` -- semicolon splitting
7. `processIntensityColumns(df)` -- detect + log2 transform

Spectronaut parser should follow the same conceptual flow but with the pivot step between parse and semtype assignment.

### Pattern 3: R Script Integration
**What:** Server-side R scripts called via `grok.functions.call()`.
**Source:** `scripts/limma_de.R` and `src/analysis/differential-expression.ts`.

R script metadata pattern:
```r
#name: VsnNormalize
#description: Variance-stabilizing normalization via vsn
#language: r
#environment: channels: [conda-forge, bioconda], dependencies: [bioconductor-vsn]
#input: dataframe exprDf
#output: dataframe result
```

TypeScript caller pattern (from `runLimmaDE`):
1. Build clean DataFrame with simple column names (s1, s2, ...) to avoid encoding issues
2. `await grok.functions.call('Proteomics:VsnNormalize', {exprDf})`
3. Copy results back to original DataFrame columns
4. Wrap in try/catch with fallback to client-side alternative

**VSN critical detail:** VSN operates on RAW (non-log2) intensities and produces variance-stabilized values on a log2-like (glog2) scale. The CONTEXT.md decision to retain both raw and log2 columns enables this.

### Pattern 4: Normalization Function Signature
**Source:** `medianNormalize()` in `src/analysis/normalization.ts`.
```typescript
// Signature: (df: DG.DataFrame, colNames: string[]): void
// In-place modification using col.getRawData() for bulk access
// df.fireValuesChanged() after
// df.setTag('proteomics.normalized', 'true')
```

### Pattern 5: Imputation Function Signature
**Source:** `imputeMinProb()` in `src/analysis/imputation.ts`.
```typescript
// Signature: (df: DG.DataFrame, colNames: string[], ...params): number
// Returns count of imputed values
// Uses col.isNone(i) for missing detection
// df.fireValuesChanged() after
// df.setTag('proteomics.imputed', 'true')
```

### Pattern 6: Function Registration in package.ts
**Source:** `src/package.ts` -- decorator-based registration.
```typescript
@grok.decorators.func({'top-menu': 'Proteomics | Import | Spectronaut...'})
static async importSpectronaut(): Promise<void> {
  DG.Utils.openFile({
    accept: '.tsv,.txt,.csv',
    open: async (file: File) => {
      const text = await file.text();
      const df = parseSpectronautText(text);
      df.name = file.name.replace(/\.[^.]+$/, '');
      grok.shell.addTableView(df);
      grok.shell.info(`Imported ${df.rowCount} protein groups`);
    },
  });
}
```

### Anti-Patterns to Avoid
- **Building a generic pivot library:** The Spectronaut pivot is specific. Don't over-generalize.
- **Modifying existing dialogs:** Phase 10 adds functions only. Dialogs are Phase 11's scope.
- **Using DG.DataFrame.fromCsv() for pivoted data:** Build the wide DataFrame using `DG.DataFrame.fromColumns()` with pre-constructed typed arrays for better control.
- **Sorting in place for quantile normalization:** Must preserve original order; use index arrays.
- **Applying VSN to log2 columns:** VSN must receive raw (pre-log2) intensity columns.
- **Using col.get()/col.set() in tight loops:** Use `col.getRawData()` for Float32/64Array access.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| CSV/TSV parsing | Custom line splitter | `DG.DataFrame.fromCsv()` with delimiter option | Handles quoting, encoding, type detection |
| Column statistics (mean, median) | Manual computation | `col.stats.avg`, `col.stats.med` | Built into Datagrok Column API |
| Progress indicators | Custom DOM elements | `DG.TaskBarProgressIndicator.create()` | Platform-native UI |
| Group annotation storage | Custom tags | `setGroups()` / `getGroups()` from experiment-setup.ts | Already established pattern |
| Log2 transformation | Manual log2 loop | `log2TransformColumns()` from shared-utils.ts | Handles nulls, zeros, semtype assignment |
| Pre-normalized detection | Custom heuristic | `detectLog2Status()` from shared-utils.ts | Already tuned and tested |
| Primary column extraction | Semicolon splitting | `addPrimaryColumnIfNeeded()` from shared-utils.ts | Handles edge cases |
| VSN normalization | TypeScript port | R's `vsn::justvsn()` via R script | Complex maximum-likelihood estimation |

**Key insight:** The existing codebase in `shared-utils.ts` and `experiment-setup.ts` already has robust utilities. The parser reuses these. The normalization/imputation functions follow identical patterns to `medianNormalize` and `imputeMinProb`.

## Common Pitfalls

### Pitfall 1: PG.IBAQ Values in Ambiguous Range
**What goes wrong:** PG.IBAQ values in the demo (34.8, 212.2, 0) partially overlap the log2 range [0,30] but some exceed it (212.2). `detectLog2Status()` checks whether >50% of values are >=1000 (raw) or >80% are in [0,30] (log2).
**Why it happens:** iBAQ = summed peptide intensities / theoretical peptides, spanning a wide range.
**How to avoid:** Run `detectLog2Status()` on the pivoted intensity columns. Demo values like 212.2 exceed log2 range, so heuristic will correctly identify them as raw. If ambiguous, default to raw (apply log2 transform).
**Warning signs:** `detectLog2Status()` returns "Could not determine data scale".

### Pitfall 2: Duplicate Protein+Sample Entries in Long Format
**What goes wrong:** Multiple peptide rows per protein per sample all carry the same PG.IBAQ value.
**Why it happens:** Long format has one row per peptide evidence, not per protein.
**How to avoid:** Use first-encountered value during pivot (Map.set). Do NOT sum PG.IBAQ -- it's already a protein-level quantity. Verify final protein count matches `unique(PG.ProteinGroups)`.
**Warning signs:** Pivoted protein count doesn't match unique protein count after filtering.

### Pitfall 3: Q-value Column Non-Numeric Values
**What goes wrong:** Spectronaut's EG.Qvalue can contain 'Profiled' or 'NaN' strings.
**Why it happens:** Spectronaut marks some identifications with text-based status.
**How to avoid:** Check `isNaN(parseFloat(val))` -- if true, treat as passing per CONTEXT.md decision. Note: when `DG.DataFrame.fromCsv` parses the column, non-numeric strings in a mixed column may cause the column to be imported as string type, so check column type.
**Warning signs:** Filter removes too many rows or all rows.

### Pitfall 4: kNN Distance with All-Missing Shared Columns
**What goes wrong:** Two proteins may have no shared non-missing values, making distance undefined.
**Why it happens:** DIA data can be 20-50% missing with complementary patterns.
**How to avoid:** Skip protein pairs with zero shared non-missing values. If a target protein has zero valid neighbors, fall back to column mean per CONTEXT.md decision.
**Warning signs:** NaN distances, imputed values that are column means for all proteins.

### Pitfall 5: Quantile Normalization with Missing Values
**What goes wrong:** Standard quantile normalization assumes complete data.
**Why it happens:** Proteomics data always has missing values.
**How to avoid:** Sort only non-missing values per column. Compute rank means from non-null values. Leave null cells null. Each column may have a different count of non-null values, requiring scaled rank alignment.
**Warning signs:** Nulls becoming zero; distributions not aligned; all values becoming NaN.

### Pitfall 6: VSN Input Must Be Raw Intensities
**What goes wrong:** Passing log2 values to VSN produces double-transformation.
**Why it happens:** VSN performs its own glog2-like transformation internally.
**How to avoid:** Extract the raw (non-log2-prefixed) intensity columns for the VSN R script. Replace log2 columns with VSN output. The CONTEXT.md decision to retain both raw and log2 columns enables this.
**Warning signs:** VSN output values in wrong range.

## Code Examples

### Spectronaut Pivot Core Logic
```typescript
// Source: Derived from demo file analysis and CONTEXT.md decisions

interface PivotResult {
  proteinMap: Map<string, Map<string, number>>;  // protein -> sample -> IBAQ
  sampleKeys: string[];   // ordered unique "Condition_Replicate" keys
  sampleFileNames: Map<string, string>;  // sampleKey -> R.FileName
  organisms: Map<string, string>;  // protein -> PG.Organisms
}

function pivotSpectronaut(
  longDf: DG.DataFrame,
  qValueThreshold: number = 0.01,
): PivotResult {
  const condCol = longDf.col('R.Condition')!;
  const replCol = longDf.col('R.Replicate')!;
  const protCol = longDf.col('PG.ProteinGroups')!;
  const ibaqCol = longDf.col('PG.IBAQ')!;
  const qvalCol = longDf.col('EG.Qvalue');
  const fileCol = longDf.col('R.FileName');
  const orgCol = longDf.col('PG.Organisms');

  const proteinMap = new Map<string, Map<string, number>>();
  const sampleSet = new Set<string>();
  const sampleFileNames = new Map<string, string>();
  const organisms = new Map<string, string>();

  for (let i = 0; i < longDf.rowCount; i++) {
    // Q-value filter: numeric values > threshold are excluded; non-numeric pass
    if (qvalCol && !qvalCol.isNone(i)) {
      const qval = Number(qvalCol.get(i));
      if (!isNaN(qval) && qval > qValueThreshold) continue;
    }

    const protein = protCol.get(i) as string;
    if (!protein) continue;
    if (protein.startsWith('CON__') || protein.startsWith('REV__')) continue;

    const condition = condCol.get(i) as string;
    const replicate = String(replCol.get(i));
    const sampleKey = `${condition}_${replicate}`;
    sampleSet.add(sampleKey);

    if (!proteinMap.has(protein))
      proteinMap.set(protein, new Map<string, number>());

    // First-encountered value wins (PG.IBAQ is constant per protein+sample)
    if (!proteinMap.get(protein)!.has(sampleKey) && !ibaqCol.isNone(i))
      proteinMap.get(protein)!.set(sampleKey, Number(ibaqCol.get(i)));

    if (fileCol && !sampleFileNames.has(sampleKey))
      sampleFileNames.set(sampleKey, fileCol.get(i) as string);
    if (orgCol && !organisms.has(protein))
      organisms.set(protein, orgCol.get(i) as string);
  }

  return {
    proteinMap,
    sampleKeys: Array.from(sampleSet).sort(),
    sampleFileNames,
    organisms,
  };
}
```

### Building Wide DataFrame from Pivot
```typescript
// Source: Pattern from DG.DataFrame.fromColumns() usage in existing tests

function buildWideDataFrame(result: PivotResult): DG.DataFrame {
  const proteins = Array.from(result.proteinMap.keys());
  const n = proteins.length;

  const proteinCol = DG.Column.fromStrings('PG.ProteinGroups', proteins);
  proteinCol.semType = SEMTYPE.PROTEIN_ID;

  const cols: DG.Column[] = [proteinCol];

  for (const sampleKey of result.sampleKeys) {
    const values = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      const sampleMap = result.proteinMap.get(proteins[i])!;
      values[i] = sampleMap.has(sampleKey) ? sampleMap.get(sampleKey)! : DG.FLOAT_NULL;
    }
    const col = DG.Column.fromFloat32Array(sampleKey, values);
    col.semType = SEMTYPE.INTENSITY;
    if (result.sampleFileNames.has(sampleKey))
      col.setTag('spectronaut.fileName', result.sampleFileNames.get(sampleKey)!);
    cols.push(col);
  }

  return DG.DataFrame.fromColumns(cols);
}
```

### Auto-Group Population from R.Condition
```typescript
// Source: experiment-setup.ts setGroups() pattern

function autoPopulateGroups(df: DG.DataFrame, sampleKeys: string[]): void {
  const conditionMap = new Map<string, string[]>();
  for (const key of sampleKeys) {
    const lastUnderscore = key.lastIndexOf('_');
    const condition = key.substring(0, lastUnderscore);
    if (!conditionMap.has(condition))
      conditionMap.set(condition, []);
    conditionMap.get(condition)!.push(`log2(${key})`);
  }

  const conditions = Array.from(conditionMap.keys());
  if (conditions.length === 2) {
    setGroups(df, {
      group1: { name: conditions[0], columns: conditionMap.get(conditions[0])! },
      group2: { name: conditions[1], columns: conditionMap.get(conditions[1])! },
    });
  }
  // >2 conditions: don't auto-assign, user picks in annotation dialog
}
```

### Quantile Normalization (Pure TypeScript)
```typescript
// Operates on log2 intensity columns in-place
export function quantileNormalize(df: DG.DataFrame, colNames: string[]): void {
  const nCols = colNames.length;
  if (nCols < 2) return;
  const nRows = df.rowCount;
  const cols = colNames.map(n => df.col(n)!);

  // For each column: build sorted indices of non-missing values
  const sortedIndices: number[][] = [];
  for (let c = 0; c < nCols; c++) {
    const col = cols[c];
    const valid: number[] = [];
    for (let r = 0; r < nRows; r++) {
      if (!col.isNone(r)) valid.push(r);
    }
    valid.sort((a, b) => (col.get(a) as number) - (col.get(b) as number));
    sortedIndices.push(valid);
  }

  // Compute rank means across columns (handling different valid counts)
  const maxLen = Math.max(...sortedIndices.map(s => s.length));
  for (let rank = 0; rank < maxLen; rank++) {
    let sum = 0;
    let count = 0;
    for (let c = 0; c < nCols; c++) {
      const indices = sortedIndices[c];
      const scaledIdx = Math.round(rank * (indices.length - 1) / Math.max(maxLen - 1, 1));
      if (scaledIdx < indices.length) {
        sum += cols[c].get(indices[scaledIdx]) as number;
        count++;
      }
    }
    const avg = sum / count;
    for (let c = 0; c < nCols; c++) {
      const indices = sortedIndices[c];
      const scaledIdx = Math.round(rank * (indices.length - 1) / Math.max(maxLen - 1, 1));
      if (scaledIdx < indices.length)
        cols[c].set(indices[scaledIdx], avg);
    }
  }

  df.fireValuesChanged();
  df.setTag('proteomics.normalized', 'true');
}
```

### VSN R Script
```r
#name: VsnNormalize
#description: Variance-stabilizing normalization via vsn
#language: r
#environment: channels: [conda-forge, bioconda], dependencies: [bioconductor-vsn]
#input: dataframe exprDf
#output: dataframe result

# exprDf columns are named s1, s2, ... containing RAW intensities (NOT log2)
exprMat <- as.matrix(exprDf)

hasVsn <- suppressWarnings(require(vsn, quietly = TRUE))

if (hasVsn) {
  # justvsn returns glog2-transformed, variance-stabilized values
  normalizedMat <- justvsn(exprMat)
  result <- as.data.frame(normalizedMat)
  colnames(result) <- colnames(exprDf)
} else {
  stop("vsn package is not available")
}
```

### VSN TypeScript Wrapper
```typescript
// Source: existing runLimmaDE pattern in differential-expression.ts
export async function vsnNormalize(
  df: DG.DataFrame, colNames: string[],
): Promise<void> {
  // Find raw (non-log2) intensity columns
  const rawColNames = colNames.map(n => n.startsWith('log2(') ?
    n.substring(5, n.length - 1) : n);

  try {
    // Build expression DataFrame from RAW intensity columns
    const exprDf = DG.DataFrame.create(df.rowCount);
    for (let i = 0; i < rawColNames.length; i++) {
      const src = df.columns.byName(rawColNames[i]);
      const dst = exprDf.columns.addNewFloat(`s${i + 1}`);
      for (let r = 0; r < df.rowCount; r++)
        dst.set(r, src.isNone(r) ? DG.FLOAT_NULL : src.get(r));
    }

    const result: DG.DataFrame = await grok.functions.call(
      'Proteomics:VsnNormalize', {exprDf},
    );

    // Copy VSN output into log2 columns (VSN output IS log2-scale)
    for (let i = 0; i < colNames.length; i++) {
      const dst = df.columns.byName(colNames[i]);
      const src = result.columns.byIndex(i);
      for (let r = 0; r < df.rowCount; r++)
        dst.set(r, src.isNone(r) ? DG.FLOAT_NULL : src.get(r));
    }

    df.fireValuesChanged();
    df.setTag('proteomics.normalized', 'true');
  } catch (e: any) {
    console.warn('VSN normalization failed, falling back to quantile:', e);
    grok.shell.warning('R environment unavailable -- using quantile normalization');
    quantileNormalize(df, colNames);
  }
}
```

### kNN Imputation
```typescript
export function imputeKnn(
  df: DG.DataFrame, colNames: string[], k: number = 10,
): number {
  const nRows = df.rowCount;
  const nCols = colNames.length;
  const cols = colNames.map(n => df.col(n)!);
  const pi = DG.TaskBarProgressIndicator.create('kNN imputation...');

  // Build matrix
  const matrix = new Float64Array(nRows * nCols);
  const missing: boolean[][] = [];
  for (let r = 0; r < nRows; r++) {
    missing.push([]);
    for (let c = 0; c < nCols; c++) {
      if (cols[c].isNone(r)) {
        matrix[r * nCols + c] = NaN;
        missing[r].push(true);
      } else {
        matrix[r * nCols + c] = Number(cols[c].get(r));
        missing[r].push(false);
      }
    }
  }

  // Column means for fallback
  const colMeans = new Float64Array(nCols);
  for (let c = 0; c < nCols; c++)
    colMeans[c] = isNaN(cols[c].stats.avg) ? 0 : cols[c].stats.avg;

  let totalImputed = 0;

  for (let r = 0; r < nRows; r++) {
    const missingCols = missing[r].reduce(
      (acc, v, i) => v ? [...acc, i] : acc, [] as number[]);
    if (missingCols.length === 0 || missingCols.length === nCols) continue;

    // Compute distances to all other rows
    const distances: {idx: number; dist: number}[] = [];
    for (let other = 0; other < nRows; other++) {
      if (other === r) continue;
      let sumSq = 0; let shared = 0;
      for (let c = 0; c < nCols; c++) {
        if (!missing[r][c] && !missing[other][c]) {
          const diff = matrix[r * nCols + c] - matrix[other * nCols + c];
          sumSq += diff * diff;
          shared++;
        }
      }
      if (shared > 0)
        distances.push({idx: other, dist: Math.sqrt(sumSq / shared)});
    }

    distances.sort((a, b) => a.dist - b.dist);
    const neighbors = distances.slice(0, k);

    for (const colIdx of missingCols) {
      let sum = 0; let count = 0;
      for (const nb of neighbors) {
        if (!missing[nb.idx][colIdx]) {
          sum += matrix[nb.idx * nCols + colIdx];
          count++;
        }
      }
      const imputed = count > 0 ? sum / count : colMeans[colIdx];
      cols[colIdx].set(r, imputed);
      totalImputed++;
    }

    pi.update(Math.round((r / nRows) * 100));
  }

  pi.close();
  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return totalImputed;
}
```

### Simple Imputation Methods
```typescript
export function imputeZero(df: DG.DataFrame, colNames: string[]): number {
  let count = 0;
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    for (let i = 0; i < col.length; i++) {
      if (col.isNone(i)) { col.set(i, 0); count++; }
    }
  }
  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return count;
}

export function imputeMean(df: DG.DataFrame, colNames: string[]): number {
  let count = 0;
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    const avg = col.stats.avg;
    if (isNaN(avg)) continue;
    for (let i = 0; i < col.length; i++) {
      if (col.isNone(i)) { col.set(i, avg); count++; }
    }
  }
  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return count;
}

export function imputeMedian(df: DG.DataFrame, colNames: string[]): number {
  let count = 0;
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    const med = col.stats.med;
    if (isNaN(med)) continue;
    for (let i = 0; i < col.length; i++) {
      if (col.isNone(i)) { col.set(i, med); count++; }
    }
  }
  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return count;
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| DDA-only proteomics tools | DIA support essential | 2022-2024 | Spectronaut parser addresses DIA workflows |
| Perseus-only analysis | Multiple platforms (FragPipe-Analyst, alphastats) | 2023-2025 | Validates need for integrated Datagrok approach |
| Single normalization method | Multiple methods with user selection | Established | Quantile + VSN + median provides proper coverage |
| Single imputation method | Multiple methods for different missingness | Established | kNN for MAR, MinProb for MNAR, simple for convenience |

**Deprecated/outdated:**
- None relevant -- all algorithms are well-established and current

## Open Questions

1. **PG.IBAQ scale in real-world Spectronaut exports**
   - What we know: Demo file values like 34.8, 212.2 (moderate range, likely raw not log2)
   - What's unclear: Whether real exports use very different scales based on Spectronaut normalization settings
   - Recommendation: Rely on `detectLog2Status()` heuristic. Values > 30 trigger "raw" classification, which is correct for the demo. Default to raw if ambiguous.

2. **kNN performance on large datasets**
   - What we know: O(n^2) where n = protein count. Demo has 93 proteins (trivial). Real DIA can have 5,000-15,000.
   - What's unclear: Whether 15,000^2 = 225M distance computations will be acceptable in main thread.
   - Recommendation: Use Float64Array typed arrays. Progress indicator keeps user informed. Web Worker explicitly out of scope per REQUIREMENTS.md.

3. **Spectronaut protein group separators**
   - What we know: Demo uses single UniProt accessions (e.g., "O15439"). MaxQuant uses semicolons.
   - What's unclear: Whether Spectronaut exports can have semicolon-delimited protein groups.
   - Recommendation: Apply `addPrimaryColumnIfNeeded()` -- it already handles the "no semicolons" case by skipping.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | @datagrok-libraries/test (Puppeteer-based) |
| Config file | webpack.config.js (package-test.ts entry point) |
| Quick run command | `grok test --host localhost --category "Spectronaut"` |
| Full suite command | `grok test --host localhost` |

### Phase Requirements to Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SPEC-01 | Spectronaut TSV pivot produces correct protein x sample matrix | unit | `grok test --host localhost --test "Spectronaut: pivot"` | No - Wave 0 |
| SPEC-02 | CON__/REV__ filtered, q-value threshold applied | unit | `grok test --host localhost --test "Spectronaut: filter"` | No - Wave 0 |
| SPEC-03 | R.Condition/R.Replicate auto-populate groups | unit | `grok test --host localhost --test "Spectronaut: groups"` | No - Wave 0 |
| SPEC-04 | Pre-normalized tag set when detected | unit | `grok test --host localhost --test "Spectronaut: preNormalized"` | No - Wave 0 |
| NORM-01 | Quantile normalization aligns distributions | unit | `grok test --host localhost --test "quantile"` | No - Wave 0 |
| NORM-02 | VSN calls R script, falls back to quantile | integration | `grok test --host localhost --test "VSN"` | No - Wave 0 |
| IMP-01 | kNN fills missing values using neighbors | unit | `grok test --host localhost --test "kNN"` | No - Wave 0 |
| IMP-02 | Zero/mean/median imputation fills correctly | unit | `grok test --host localhost --test "impute"` | No - Wave 0 |

### Sampling Rate
- **Per task commit:** `grok test --host localhost --category "Spectronaut"` or relevant category
- **Per wave merge:** `grok test --host localhost` (full suite)
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `src/tests/spectronaut-parser.ts` -- covers SPEC-01 through SPEC-04 with synthetic long-format data
- [ ] Tests for NORM-01, NORM-02, IMP-01, IMP-02 in `src/tests/analysis.ts` (extend existing file)
- [ ] Test helper: `makeLongFormatTsv()` utility for building Spectronaut-format test data
- [ ] Register new test file in `src/package-test.ts`

## Sources

### Primary (HIGH confidence)
- Existing codebase: `src/parsers/maxquant-parser.ts` (147 lines), `src/parsers/shared-utils.ts` (166 lines), `src/analysis/normalization.ts` (53 lines), `src/analysis/imputation.ts` (83 lines), `src/analysis/differential-expression.ts` (360 lines), `src/analysis/experiment-setup.ts` (63 lines), `src/package.ts` (246 lines) -- all fully read and analyzed
- Demo file: `files/demo/spectronaut-hye-mix.tsv` -- 35 columns, 8760 rows, 93 proteins, 2 conditions x 4 replicates verified
- R scripts: `scripts/limma_de.R` (38 lines), `scripts/deqms_de.R` -- R integration pattern verified
- Existing tests: `tests/parsers.ts` (179 lines), `tests/analysis.ts` (294 lines) -- test pattern verified

### Secondary (MEDIUM confidence)
- Quantile normalization algorithm: well-established (Bolstad et al. 2003); sort-rank-mean is standard
- VSN: Bioconductor `vsn` package; `justvsn()` function is the standard one-call interface
- kNN imputation: standard algorithm; uniform weighting recommended for MNAR proteomics data

### Tertiary (LOW confidence)
- PG.IBAQ scale variability across different Spectronaut versions -- based on domain knowledge only

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - zero new dependencies, all patterns established in codebase
- Architecture: HIGH - mirrors existing parser and analysis module patterns exactly
- Pitfalls: HIGH - verified against actual demo data structure
- Algorithm correctness: MEDIUM - standard algorithms but null-handling details need careful testing

**Research date:** 2026-03-07
**Valid until:** 2026-04-07 (stable -- no moving targets, well-established algorithms)
