# Phase 6: Generic Matrix Parser - Research

**Researched:** 2026-03-05
**Domain:** CSV/TSV proteomics matrix import with column mapping UI
**Confidence:** HIGH

## Summary

Phase 6 adds a generic CSV/TSV matrix importer to the Proteomics package. Unlike the MaxQuant parser, which auto-detects columns by known naming conventions, the generic parser presents a column mapping dialog where the user selects protein ID, gene name, and intensity columns from any arbitrary matrix file. The result must produce the exact same semantic type annotations (`SEMTYPE.INTENSITY`, `SEMTYPE.PROTEIN_ID`, `SEMTYPE.GENE_SYMBOL`) and `log2()` naming convention that the existing pipeline expects.

The codebase already contains most building blocks: `processIntensityColumns()` for log2 transformation, `addPrimaryColumn()` for semicolon splitting, `DG.DataFrame.fromCsv()` for parsing, and `ui.input.column()`/`ui.input.columns()` for column pickers. The main new work is (1) the column mapping dialog with auto-suggestion and live preview, (2) a shared log2 transform utility extracted from maxquant-parser.ts, (3) log2 auto-detection heuristic, and (4) backporting the dataset input widget and source tags to the MaxQuant importer.

**Primary recommendation:** Extract shared utilities from maxquant-parser.ts into a new `src/parsers/shared-utils.ts`, build the generic parser dialog in `src/parsers/generic-parser.ts`, and backport dataset input + source tags to the MaxQuant importer.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Top menu entry only: "Proteomics | Import | Generic Matrix..." -- no file handler sniffing
- Accept CSV and TSV files (.csv, .tsv, .txt) with auto-detected delimiter
- Use Datagrok's built-in dataset input widget (table selector + file picker + folder browser + database icons) -- NOT a bare OS file dialog
- Backport the same dataset input widget to the MaxQuant importer for consistency
- Single dialog with all fields visible at once
- Layout: dataset input (top), protein ID column dropdown, gene name column dropdown (optional), intensity columns multi-select, log2 transform toggle, data preview grid (bottom)
- Auto-suggest intensity columns by keyword matching: intensity, lfq, ibaq, tmt, reporter, abundance
- Auto-suggest protein ID column by keyword matching: protein, accession, uniprot
- Live-updating preview: 5 rows, updates as column selections change to show only mapped columns
- Import and Cancel buttons at bottom
- Auto-detect whether data is already log2-transformed by checking value ranges (0-30 = log2, 1e3+ = raw)
- Pre-set the log2 toggle based on detection, with a hint message
- User can override the auto-detection
- Keep both original intensity columns AND log2-prefixed copies (mirror MaxQuant behavior)
- Always rename with log2() prefix for downstream compatibility, even if data was already log2-transformed
- Assign SEMTYPE.INTENSITY to both original and log2 columns
- Assign SEMTYPE.PROTEIN_ID to selected protein ID column
- Assign SEMTYPE.GENE_SYMBOL to selected gene name column (if provided)
- Create 'Primary Protein ID' / 'Primary Gene Name' columns only when semicolons detected in source data
- No row filtering -- generic import keeps all rows as-is
- Tag source format: `proteomics.source = 'generic'` -- backport `proteomics.source = 'maxquant'` to MaxQuant importer
- DataFrame name reflects source dataset/file name -- backport this to MaxQuant importer too

### Claude's Discretion
- Exact auto-detection threshold for log2 vs raw intensities
- Column picker filter configuration details
- Preview grid styling and column ordering
- Error handling for malformed files

### Deferred Ideas (OUT OF SCOPE)
- Tool-specific parsers (FragPipe, DIA-NN, Spectronaut, Proteome Discoverer) -- tracked as IMPORT-07 through IMPORT-10
- Optional row filter column in import dialog
- Excel (.xlsx) file support
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| IMPORT-01 | User can import a generic CSV/TSV matrix file via menu entry | Top menu "Proteomics \| Import \| Generic Matrix..."; dataset input widget for file selection; `DG.DataFrame.fromCsv()` with auto-detected delimiter |
| IMPORT-02 | User sees a column mapping dialog to select protein ID column and intensity columns | `ui.input.column()` for single column picker, `ui.input.columns()` for multi-select; `ui.dialog()` for dialog construction |
| IMPORT-03 | System auto-suggests likely protein ID and intensity columns based on naming patterns | Column name keyword matching against known patterns; `filter` option on column/columns inputs |
| IMPORT-04 | User sees a data preview (first rows) during column mapping | Grid viewer via `DG.Viewer.grid(previewDf).root` embedded in dialog; clone first 5 rows for preview |
| IMPORT-05 | User can optionally log2-transform intensities on import | `ui.input.bool()` for toggle; auto-detect log2 vs raw via value range heuristic; extracted shared log2 transform utility |
| IMPORT-06 | Imported generic matrix receives same semantic types as MaxQuant import | Reuse `SEMTYPE` constants; assign via `col.semType`; tag with `proteomics.source = 'generic'`; `log2()` prefix naming |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| datagrok-api | (platform) | DataFrame, Column, UI inputs, dialog, grid | Platform API -- all UI and data operations |
| proteomics-types.ts | local | SEMTYPE constants | Already defines all 5 semantic types used by pipeline |
| column-detection.ts | local | findColumn / findProteomicsColumns | Used by all downstream viewers |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| @datagrok-libraries/test | (repo) | Test framework (`category`, `test`, `expect`) | Unit tests for parser logic |

No additional npm dependencies needed. Everything uses the Datagrok platform API.

**Installation:** No new packages to install.

## Architecture Patterns

### Recommended Project Structure
```
src/
  parsers/
    shared-utils.ts       # NEW: extracted log2 transform, addPrimaryColumnIfNeeded
    maxquant-parser.ts    # MODIFIED: import from shared-utils, add source tag, accept file name
    generic-parser.ts     # NEW: column mapping dialog + generic import logic
  utils/
    proteomics-types.ts   # EXISTING: SEMTYPE constants (no changes)
    column-detection.ts   # EXISTING: findColumn, findProteomicsColumns (no changes)
  package.ts              # MODIFIED: add importGenericMatrix menu entry, update importMaxQuant
  tests/
    parsers.ts            # EXISTING: MaxQuant parser tests (no changes needed)
    generic-parser.ts     # NEW: generic parser unit tests
```

### Pattern 1: Shared Utility Extraction
**What:** Extract reusable functions from maxquant-parser.ts into shared-utils.ts so both parsers share identical log2 transform and primary column logic.
**When to use:** Both MaxQuant and generic parsers need log2 transform and semicolon-delimited field extraction.

```typescript
// src/parsers/shared-utils.ts
import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';

/** Creates log2-transformed copies of the specified columns.
 * Original columns get SEMTYPE.INTENSITY. New log2 columns get SEMTYPE.INTENSITY
 * and names prefixed with 'log2(...)'. */
export function log2TransformColumns(df: DG.DataFrame, colNames: string[]): void {
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    col.semType = SEMTYPE.INTENSITY;
    const log2Name = `log2(${name})`;
    const log2Col = df.columns.addNewFloat(log2Name);
    for (let i = 0; i < df.rowCount; i++) {
      if (col.isNone(i))
        log2Col.set(i, DG.FLOAT_NULL);
      else {
        const val = col.get(i) as number;
        log2Col.set(i, val > 0 ? Math.log2(val) : DG.FLOAT_NULL);
      }
    }
    log2Col.semType = SEMTYPE.INTENSITY;
  }
}

/** Copies intensity columns with log2() prefix for already-transformed data.
 * No math transformation -- just copies values into log2-named columns
 * for downstream compatibility. */
export function copyAsLog2Columns(df: DG.DataFrame, colNames: string[]): void {
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    col.semType = SEMTYPE.INTENSITY;
    const log2Name = `log2(${name})`;
    const log2Col = df.columns.addNewFloat(log2Name);
    for (let i = 0; i < df.rowCount; i++) {
      if (col.isNone(i))
        log2Col.set(i, DG.FLOAT_NULL);
      else
        log2Col.set(i, col.get(i) as number);
    }
    log2Col.semType = SEMTYPE.INTENSITY;
  }
}

/** Creates a primary (first semicolon-delimited entry) column if semicolons
 * are detected in source data. Skips creation if no semicolons found. */
export function addPrimaryColumnIfNeeded(
  df: DG.DataFrame, sourceColName: string, newName: string, semType: string,
): void {
  const srcCol = df.col(sourceColName);
  if (!srcCol) return;
  let hasSemicolon = false;
  for (let i = 0; i < df.rowCount && !hasSemicolon; i++) {
    const val = srcCol.get(i);
    if (typeof val === 'string' && val.includes(';'))
      hasSemicolon = true;
  }
  if (!hasSemicolon) return;
  const newCol = df.columns.addNewString(newName);
  for (let i = 0; i < df.rowCount; i++) {
    const val = srcCol.get(i);
    if (typeof val === 'string' && val.length > 0) {
      const idx = val.indexOf(';');
      newCol.set(i, idx >= 0 ? val.substring(0, idx) : val);
    }
  }
  newCol.semType = semType;
}
```

### Pattern 2: Column Mapping Dialog
**What:** A single dialog with dataset input, column pickers, log2 toggle, hint text, and live preview grid.
**When to use:** Generic matrix import.

Key API patterns verified from codebase (js-api/ui.ts lines 1061-1067, 901-913, 1037-1039):

```typescript
// Single column picker with type filter
const proteinIdInput = ui.input.column('Protein ID', {
  table: df,
  filter: (col: DG.Column) => col.type === DG.COLUMN_TYPE.STRING,
  value: autoDetectedProteinCol,
  onValueChanged: (value) => updatePreview(),
});

// Multi-select column picker with available list and pre-selection
const intensityCols = ui.input.columns('Intensity Columns', {
  table: df,
  available: suggestedIntensityNames,
  value: autoSelectedColumns,
  onValueChanged: (value) => updatePreview(),
});

// Boolean toggle
const log2Toggle = ui.input.bool('Log2 Transform', {
  value: !detectedAsLog2,
  onValueChanged: (value) => updateHintMessage(),
});

// Dialog assembly (pattern from experiment-setup.ts)
ui.dialog('Import Generic Matrix')
  .add(datasetInput)
  .add(proteinIdInput)
  .add(geneNameInput)
  .add(intensityCols)
  .add(log2Toggle)
  .add(hintDiv)
  .add(previewContainer)
  .onOK(() => { /* finalize import */ })
  .show();
```

### Pattern 3: Log2 Auto-Detection Heuristic
**What:** Determine if intensity data is already log2-transformed by sampling values.
**When to use:** Pre-setting the log2 toggle and showing a hint message.

Recommended thresholds (Claude's discretion area):

```typescript
/** Detects whether intensity columns appear to be log2-transformed.
 * Samples up to 200 non-null values across intensity columns.
 * Heuristic: values in [0, 30] are log2-scale; values >= 1000 are raw scale. */
function detectLog2Status(
  df: DG.DataFrame, colNames: string[],
): {isLog2: boolean; message: string} {
  let inLog2Range = 0;  // values 0-30
  let inRawRange = 0;   // values >= 1000
  let total = 0;
  const maxSamples = 200;

  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    if (col.type !== DG.COLUMN_TYPE.FLOAT && col.type !== DG.COLUMN_TYPE.INT) continue;
    for (let i = 0; i < df.rowCount && total < maxSamples; i++) {
      if (col.isNone(i)) continue;
      const val = col.get(i) as number;
      if (val >= 0 && val <= 30) inLog2Range++;
      if (val >= 1000) inRawRange++;
      total++;
    }
    if (total >= maxSamples) break;
  }

  if (total === 0)
    return {isLog2: false, message: 'No intensity values found'};
  if (inRawRange / total > 0.5)
    return {isLog2: false, message: 'Data appears to be raw intensities'};
  if (inLog2Range / total > 0.8)
    return {isLog2: true, message: 'Data appears to be already log2-transformed'};
  return {isLog2: false, message: 'Could not determine data scale'};
}
```

### Pattern 4: Dataset Input Widget
**What:** Use Datagrok's built-in dataset input widget instead of bare OS file dialog.
**When to use:** Both MaxQuant and generic importers (backport for consistency).

The CONTEXT.md specifies the "built-in dataset input widget (table selector + file picker + folder browser + database icons)." Based on the API:

- `ui.input.table()` -- table selector dropdown for tables already open in the platform
- `ui.input.file()` -- file picker for platform file system
- `DG.Utils.openFile()` -- OS file dialog (current MaxQuant approach)

**Implementation approach:** The dialog should first let the user load a file. Two viable options:

1. **Two-step flow:** A "Select File" button using `DG.Utils.openFile()` that parses the file, then the column mapping dialog appears below with all inputs enabled.
2. **Integrated flow:** Use `ui.input.file()` at the top of the dialog, and when a file is selected, parse it and populate the column pickers below.

**Recommendation:** Use `DG.Utils.openFile()` for file selection since it is the proven pattern in this codebase. Parse the file immediately, then show the full column mapping dialog with all fields populated. For the "dataset input widget" aspect, this can be enhanced later if the exact composite widget API is found. The key constraint from CONTEXT.md is that both importers should use the same approach.

**Confidence:** MEDIUM for the exact widget component. HIGH for the overall two-step approach.

### Pattern 5: Delimiter Auto-Detection
**What:** Detect whether file is CSV or TSV before parsing.
**When to use:** When the file extension does not reliably indicate the delimiter.

```typescript
/** Detects delimiter by checking the first line for tab vs comma prevalence. */
function detectDelimiter(text: string): string {
  const firstLine = text.substring(0, text.indexOf('\n'));
  const tabs = (firstLine.match(/\t/g) || []).length;
  const commas = (firstLine.match(/,/g) || []).length;
  return tabs > commas ? '\t' : ',';
}
```

### Anti-Patterns to Avoid
- **Duplicating log2 logic:** Do NOT copy processIntensityColumns() into generic-parser.ts. Extract to shared-utils.ts and import from both parsers.
- **Hardcoding column names:** The generic parser must use user-selected column names, not hardcoded MaxQuant column names.
- **Modifying original columns in-place:** The downstream pipeline expects BOTH the original columns AND `log2()` copies. Never rename -- always create new columns.
- **Ignoring the `log2()` prefix contract:** Experiment-setup dialog filters for columns with `c.name.startsWith('log2(')` AND `c.semType === SEMTYPE.INTENSITY`. Both conditions must be met.
- **Making gene name required:** Many proteomics formats lack gene name columns. The pipeline handles this gracefully (volcano uses protein ID as fallback).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| CSV/TSV parsing | Custom delimiter parser | `DG.DataFrame.fromCsv(text, {delimiter})` | Handles quoting, encoding, type inference, null detection |
| Column picker UI | Custom dropdown with checkboxes | `ui.input.column()` / `ui.input.columns()` | Platform-native, handles filtering, consistent UX |
| Dialog layout | Raw HTML form | `ui.dialog()` with `.add()` chain | Platform-styled, handles OK/Cancel, keyboard shortcuts |
| Data grid preview | Custom HTML table | Grid viewer embedded via `.root` | Full Datagrok grid with scrolling, formatting |
| Delimiter detection | Complex regex heuristics | Simple tab-vs-comma count on first line | Sufficient for proteomics matrices (always one or the other) |

**Key insight:** The Datagrok platform already handles CSV parsing, column pickers, and dialog construction. The custom work is limited to (1) auto-suggestion heuristics, (2) log2 detection, and (3) wiring everything together.

## Common Pitfalls

### Pitfall 1: Column Names Not Starting with `log2(`
**What goes wrong:** Downstream pipeline silently ignores intensity columns because experiment-setup.ts, normalization.ts, and imputation.ts filter for `c.name.startsWith('log2(')`.
**Why it happens:** Generic parser might assign semType but forget the naming convention.
**How to avoid:** Always create `log2(originalName)` columns. If data is already log2, copy values as-is into `log2(name)` columns (use `copyAsLog2Columns`). If data is raw, apply Math.log2 transform (use `log2TransformColumns`).
**Warning signs:** Normalization/imputation/DE dialogs show empty column lists after generic import.

### Pitfall 2: Delimiter Auto-Detection Failure
**What goes wrong:** CSV files with comma delimiters parsed incorrectly, or TSV files with `.csv` extension.
**Why it happens:** `DG.DataFrame.fromCsv()` needs the correct delimiter hint.
**How to avoid:** Sniff the first line for tab vs comma count before parsing.
**Warning signs:** Single column with all data concatenated, or wrong column count.

### Pitfall 3: Intensity Column Type Coercion
**What goes wrong:** Large integer intensity values get parsed as wrong type by `fromCsv()`.
**Why it happens:** Same issue maxquant-parser solved with `buildIntensityColumnOptions()`.
**How to avoid:** After initial parse, verify user-selected intensity columns are numeric. If not, re-parse with `columnImportOptions` forcing those columns to double.
**Warning signs:** Non-numeric values in intensity columns, log2 transform producing all nulls.

### Pitfall 4: Preview Grid Not Updating
**What goes wrong:** Live preview does not reflect current column selections.
**Why it happens:** Grid viewer not refreshed when inputs change.
**How to avoid:** Use a container div. On each column selection change, create a new preview DataFrame (5 rows, selected columns only), create a new grid viewer, and replace the container's children. Keep it lightweight.
**Warning signs:** Preview shows stale columns.

### Pitfall 5: DataFrame Name Not Set
**What goes wrong:** Tab title shows "data" instead of meaningful file name.
**Why it happens:** `DG.DataFrame.fromCsv()` assigns a default name.
**How to avoid:** Set `df.name` explicitly from the source file name (strip extension).
**Warning signs:** Tab title shows generic name.

### Pitfall 6: Semicolon-Delimited IDs Not Split
**What goes wrong:** UniProt panel fails to look up proteins because IDs contain semicolons.
**Why it happens:** Many proteomics tools output semicolon-delimited protein group lists.
**How to avoid:** Use `addPrimaryColumnIfNeeded()` -- only creates primary column when semicolons are actually detected in the data.
**Warning signs:** UniProt panel shows "not found" for valid proteins.

## Code Examples

### Existing Dialog Pattern (from experiment-setup.ts)
```typescript
// Source: packages/Proteomics/src/analysis/experiment-setup.ts
const group1Cols = ui.input.columns('Group 1', {table: df, available: intensityColNames, value: []});
ui.dialog('Annotate Experiment')
  .add(group1Name)
  .add(group1Cols)
  .onOK(() => {
    const g1: DG.Column[] = group1Cols.value;
    // process...
  })
  .show();
```

### Column Input with Filter (from API)
```typescript
// Source: js-api/ui.ts lines 1061-1063
// IColumnInputInitOptions: { table, filter, value, onValueChanged, ... }
const proteinInput = ui.input.column('Protein ID', {
  table: df,
  filter: (col: DG.Column) => col.type === DG.COLUMN_TYPE.STRING,
  value: autoDetectedProteinCol,
});
```

### Columns Input with Available/Value (from API)
```typescript
// Source: js-api/ui.ts lines 1065-1067
// IColumnsInputInitOptions: { table, filter, available, checked, value, ... }
const intensityInput = ui.input.columns('Intensity Columns', {
  table: df,
  available: suggestedNames,
  value: autoSelectedCols,
});
```

### CSV Parsing with Delimiter (from maxquant-parser.ts)
```typescript
// Source: packages/Proteomics/src/parsers/maxquant-parser.ts line 148
const raw = DG.DataFrame.fromCsv(text, {delimiter: '\t', columnImportOptions});
```

### CsvImportOptions Type Definition
```typescript
// Source: js-api/src/const.ts lines 992-996
type CsvImportColumnOptions = {name: string, type?: string, semType?: string};
type CsvImportOptions = {
  delimiter?: string, decimalSeparator?: string, thousandSeparator?: string,
  headerRow?: boolean, columnFilterNames?: string[], columnFilterRegexp?: string,
  mergeDelimiters?: boolean, maxRows?: number, doublePrecision?: boolean,
  rowFilterTop?: number, rowFilterProb?: number, nullStrings?: string[],
  columnImportOptions?: CsvImportColumnOptions[]
};
```

### Downstream Column Discovery (critical contract)
```typescript
// Source: packages/Proteomics/src/analysis/normalization.ts lines 35-37
// Source: packages/Proteomics/src/analysis/experiment-setup.ts lines 31-33
// BOTH conditions must be met for downstream to find intensity columns:
const log2ColNames = df.columns.toList()
  .filter((c) => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('))
  .map((c) => c.name);
```

### Existing Log2 Transform (from maxquant-parser.ts)
```typescript
// Source: packages/Proteomics/src/parsers/maxquant-parser.ts lines 90-125
// This logic should be extracted to shared-utils.ts
const log2Name = `log2(${name})`;
const log2Col = df.columns.addNewFloat(log2Name);
for (let i = 0; i < df.rowCount; i++) {
  if (col.isNone(i))
    log2Col.set(i, DG.FLOAT_NULL);
  else {
    const val = col.get(i) as number;
    log2Col.set(i, val > 0 ? Math.log2(val) : DG.FLOAT_NULL);
  }
}
log2Col.semType = SEMTYPE.INTENSITY;
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `DG.Utils.openFile()` bare file dialog | Dataset input widget | Phase 6 decision | Both importers should use the richer widget |
| MaxQuant-only import | MaxQuant + Generic import paths | Phase 6 | Any proteomics CSV/TSV can enter the pipeline |
| No source tagging | `proteomics.source = 'maxquant'\|'generic'` | Phase 6 | Downstream can distinguish import origin |
| No DataFrame name from file | `df.name` set from file name | Phase 6 | Better UX in tab titles |

## Open Questions

1. **Dataset Input Widget API**
   - What we know: CONTEXT.md specifies the "built-in dataset input widget" with table selector, file picker, folder browser, database icons. `ui.input.table()` provides table selector. `ui.input.file()` provides file picker. `DG.Utils.openFile()` provides OS file dialog.
   - What's unclear: The exact Datagrok API call to create the full composite widget shown in the user's screenshot reference. It may be a single platform component or may require manual assembly of individual inputs.
   - Recommendation: Start with `DG.Utils.openFile()` as the proven pattern, and validate the composite widget API at implementation time. The core column mapping dialog is independent of how the file is loaded. Confidence: MEDIUM.

2. **Preview Grid in Dialog**
   - What we know: `ui.dialog().add()` accepts `HTMLElement | Widget | InputBase`. Grid viewers expose a `.root` HTMLElement.
   - What's unclear: Whether the grid viewer renders correctly inside a dialog container with fixed height.
   - Recommendation: Create a container div with fixed height (e.g., 150px), add grid viewer's `.root` inside it. Test at implementation time. Confidence: MEDIUM.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | @datagrok-libraries/test (Puppeteer-based) |
| Config file | packages/Proteomics/src/package-test.ts |
| Quick run command | `grok test --category "Generic Parser" --host localhost` |
| Full suite command | `grok test --host localhost` |

### Phase Requirements to Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| IMPORT-01 | Parse CSV/TSV text into DataFrame with correct column count | unit | `grok test --test "Generic: parses" --host localhost` | No -- Wave 0 |
| IMPORT-02 | Column mapping assigns correct semTypes to selected columns | unit | `grok test --test "Generic: assigns" --host localhost` | No -- Wave 0 |
| IMPORT-03 | Auto-suggestion finds protein ID and intensity columns by name patterns | unit | `grok test --test "Generic: auto" --host localhost` | No -- Wave 0 |
| IMPORT-04 | Data preview grid shows mapped columns | manual-only | Manual UI verification -- dialog interaction required | N/A |
| IMPORT-05 | Log2 transform produces correct values; skip-transform copies values | unit | `grok test --test "Generic: log2" --host localhost` | No -- Wave 0 |
| IMPORT-06 | Downstream pipeline (normalization, DE) works with generic import output | integration | `grok test --test "Generic: downstream" --host localhost` | No -- Wave 0 |

### Sampling Rate
- **Per task commit:** `grok test --category "Generic Parser" --host localhost`
- **Per wave merge:** `grok test --host localhost` (full suite including existing Parsers and Analysis tests)
- **Phase gate:** Full suite green before verification

### Wave 0 Gaps
- [ ] `src/tests/generic-parser.ts` -- unit tests for generic parser logic (IMPORT-01 through IMPORT-06)
- [ ] Import `./tests/generic-parser` in `package-test.ts`
- [ ] Test helper: `makeGenericCsv()` function for building test CSV strings with various column layouts

### Testable Logic (unit-testable without UI)
These functions can be tested in isolation without dialog interaction:
- **Auto-suggestion functions:** Given column names, identify protein ID and intensity candidates
- **Log2 detection heuristic:** Given a DataFrame with known value ranges, detect log2 vs raw
- **log2TransformColumns:** Given raw intensity values, verify log2 output + semType + naming
- **copyAsLog2Columns:** Given pre-transformed values, verify copy + semType + naming
- **addPrimaryColumnIfNeeded:** Given semicolon-delimited data, verify extraction; no semicolons = no new column
- **Delimiter detection:** Given CSV vs TSV text, verify correct delimiter returned
- **Source tagging:** Verify `proteomics.source` tag is set to `'generic'`
- **DataFrame naming:** Verify df.name reflects the file name
- **Downstream compatibility:** Feed generic-imported DataFrame into `medianNormalize()`, `runDifferentialExpression()` -- verify they find the columns

### Manual-Only Tests (require dialog interaction)
- IMPORT-04: Preview grid updates live as columns are selected/deselected
- Dialog layout: all fields visible at once, correct ordering
- Dataset input widget: file picker works correctly

## Sources

### Primary (HIGH confidence)
- `packages/Proteomics/src/parsers/maxquant-parser.ts` -- log2 transform, semType assignment, primary column patterns
- `packages/Proteomics/src/analysis/experiment-setup.ts` -- dialog construction with `ui.input.columns()`
- `packages/Proteomics/src/analysis/normalization.ts` -- downstream column filtering contract (`name.startsWith('log2(')` + `semType === SEMTYPE.INTENSITY`)
- `packages/Proteomics/src/utils/proteomics-types.ts` -- SEMTYPE constants
- `packages/Proteomics/src/utils/column-detection.ts` -- column discovery patterns
- `js-api/ui.ts` lines 780-1129 -- `ui.input.column()`, `ui.input.columns()`, `ui.input.bool()`, `ui.input.table()` APIs with full option interfaces
- `js-api/src/const.ts` lines 992-996 -- `CsvImportOptions` and `CsvImportColumnOptions` type definitions
- `js-api/src/utils.ts` line 204 -- `DG.Utils.openFile()` API
- `packages/Proteomics/src/tests/parsers.ts` -- existing test patterns for parser functions
- `packages/Proteomics/src/tests/analysis.ts` -- existing test patterns for analysis functions

### Secondary (MEDIUM confidence)
- Dataset input widget composite component -- referenced in CONTEXT.md with screenshot; exact API needs runtime validation

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all APIs verified in codebase, no external dependencies
- Architecture: HIGH - patterns directly extracted from existing maxquant-parser.ts and analysis dialogs
- Pitfalls: HIGH - identified from actual code contracts (log2 prefix requirement, semType filtering)
- Dataset input widget: MEDIUM - exact composite widget API needs runtime validation
- Log2 detection heuristic: MEDIUM - thresholds are reasonable but may need tuning with real-world files

**Research date:** 2026-03-05
**Valid until:** 2026-04-05 (stable -- Datagrok API, no fast-moving dependencies)
