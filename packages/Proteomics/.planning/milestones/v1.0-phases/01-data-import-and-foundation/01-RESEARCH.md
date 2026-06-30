# Phase 1: Data Import and Foundation - Research

**Researched:** 2026-02-28
**Domain:** MaxQuant proteinGroups.txt parsing, Datagrok DataFrame API, semantic types
**Confidence:** HIGH

## Summary

Phase 1 implements a MaxQuant proteinGroups.txt parser within the existing Datagrok package scaffold. The scaffold already provides the file structure (`src/parsers/`, `src/utils/`, `detectors.js`), semantic type constants (`SEMTYPE`), column detection utilities, menu item stubs, and the decorator polyfill pattern. The core work is implementing `parseMaxQuant()` in `maxquant-parser.ts`, wiring it to `importMaxQuant()` in `package.ts`, creating the demo dataset loader, and bundling a small public proteinGroups.txt file.

The technical domain is well-understood: MaxQuant's proteinGroups.txt is a tab-delimited file with stable column naming across versions 1.6.x through 2.0.x. Datagrok's `DG.DataFrame.fromCsv()` accepts a `delimiter` option for tab-separated files. Semantic types are assigned via `column.semType = 'value'`. No external libraries are needed beyond `datagrok-api`.

**Primary recommendation:** Implement the parser as a pure TypeScript function that reads TSV text, parses via `DG.DataFrame.fromCsv()` with `{delimiter: '\t'}`, filters rows, detects and log2-transforms intensity columns, and assigns semantic types directly. Wire to menu via existing stubs.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Use menu item pattern (`Proteomics | Import | MaxQuant...`) already scaffolded in package.ts
- Follow Datagrok file-handler pattern for auto-detection of proteinGroups.txt format
- Parser function signature already defined: `parseMaxQuant(file: DG.FileInfo): Promise<DG.DataFrame>`
- Filter contaminants (CON__ prefix or `+` in "Potential contaminant" column), reverse hits (REV__ prefix or `+` in "Reverse" column), and only-identified-by-site rows
- Handle both MaxQuant 1.x and 2.x column naming conventions
- Auto-detect intensity type by prefix pattern: `LFQ intensity`, `Intensity`, `iBAQ`, `Reporter intensity`
- Keep raw intensity columns and add log2-transformed columns alongside them
- Assign Proteomics-Intensity semantic type to intensity columns
- Use existing SEMTYPE constants from `src/utils/proteomics-types.ts`
- Existing detectors in `detectors.js` handle detection for subsequently opened tables
- Parser should explicitly assign semantic types on import (don't rely solely on detectors)
- Use a small, publicly available MaxQuant proteinGroups.txt (e.g., PXD000561 HeLa benchmark or similar)
- Bundle in `files/demo/` directory
- Keep small enough to load quickly (<5MB)
- Parse semicolon-delimited Protein IDs into primary ID (first entry)
- Parse semicolon-delimited Gene Names into primary gene name
- Scaffold already has the full file structure: parsers/, utils/, detectors.js

### Claude's Discretion
- Exact import dialog UI layout and options
- Column naming convention for log2-transformed columns
- How to display filtering summary (info balloon vs silent)
- Whether to auto-open a table view after import
- Demo dataset selection (any public MaxQuant output that's small and representative)
- Whether to show before/after protein count after filtering
- Error handling for malformed files

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| IMPORT-01 | Import MaxQuant proteinGroups.txt with automatic contaminant and reverse hit filtering | Parser implementation using `DG.DataFrame.fromCsv()` with tab delimiter, row filtering by "Potential contaminant", "Reverse", "Only identified by site" columns and ID prefix patterns |
| IMPORT-02 | Auto-detect intensity column type (LFQ, iBAQ, Reporter) and assign semantic types | Column name prefix matching (`LFQ intensity`, `iBAQ`, `Reporter intensity`, `Intensity`) with `column.semType` assignment |
| IMPORT-03 | Apply log2 transformation to intensity columns on import | Create new `log2(columnName)` columns using `Math.log2()`, handling zero/negative values by setting to NaN |
| IMPORT-04 | Assign semantic types to protein ID, gene name, log2FC, p-value, and intensity columns | Direct assignment via `column.semType = SEMTYPE.X` using constants from `proteomics-types.ts` |
| IMPORT-05 | Bundle a public MaxQuant demo dataset for testing and demonstrations | Place TSV file in `files/demo/`, load via `_package.files.readAsText('demo/...')` |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| datagrok-api | ^1.25.0 | DataFrame API, UI, shell, file I/O | Platform API -- the only way to interact with Datagrok |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| @datagrok-libraries/test | ^1.1.0 | Test framework (`category`, `test`, `expect`) | Already in devDependencies, used by `package-test.ts` |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `DG.DataFrame.fromCsv()` | Manual TSV line splitting | `fromCsv` handles type detection, quoting, edge cases; manual parsing is error-prone |
| Direct `column.semType =` | Relying solely on `detectors.js` | Detectors run on subsequently opened tables, not on import; explicit assignment ensures types are set immediately |

**Installation:**
No additional packages needed. All dependencies already in `package.json`.

## Architecture Patterns

### Recommended Project Structure
```
src/
├── package.ts              # Entry point -- wire importMaxQuant() and proteomicsDemo()
├── parsers/
│   └── maxquant-parser.ts  # Core parser: parseMaxQuant(file) -> DataFrame
├── utils/
│   ├── proteomics-types.ts # SEMTYPE constants (exists)
│   └── column-detection.ts # findColumn, findProteomicsColumns (exists)
├── tests/
│   └── parsers.ts          # Parser tests
└── package-test.ts         # Test entry (exists)
files/
└── demo/
    └── proteinGroups.txt   # Bundled demo dataset
```

### Pattern 1: File Import via Menu Item + DG.Utils.openFile
**What:** User clicks menu item, browser file picker opens, file content is read and parsed into a DataFrame
**When to use:** For the `importMaxQuant()` menu handler
**Example:**
```typescript
// In package.ts importMaxQuant()
static async importMaxQuant(): Promise<void> {
  DG.Utils.openFile({
    accept: '.txt,.tsv',
    open: async (file: File) => {
      const text = await file.text();
      const df = parseMaxQuantText(text);
      grok.shell.addTableView(df);
    },
  });
}
```

### Pattern 2: TSV Parsing with DG.DataFrame.fromCsv
**What:** Use Datagrok's built-in CSV parser with tab delimiter for the initial parse
**When to use:** First step of parseMaxQuant -- get raw data into a DataFrame
**Example:**
```typescript
// Parse tab-delimited proteinGroups.txt
const df = DG.DataFrame.fromCsv(tsvText, {delimiter: '\t'});
```

### Pattern 3: Row Filtering via BitSet
**What:** Use DataFrame filter/selection BitSet to mark rows for removal, then create a filtered copy
**When to use:** Removing contaminant, reverse, and only-identified-by-site rows
**Example:**
```typescript
// Filter rows where "Potential contaminant" column has "+"
const contCol = df.col('Potential contaminant');
if (contCol) {
  for (let i = 0; i < df.rowCount; i++) {
    if (contCol.get(i) === '+')
      df.filter.set(i, false);
  }
}
// Create filtered DataFrame
const filtered = df.clone(df.filter);
```

### Pattern 4: Semantic Type Assignment on Import
**What:** Explicitly assign semantic types to recognized columns after parsing
**When to use:** In the parser, after DataFrame is constructed
**Example:**
```typescript
import {SEMTYPE} from '../utils/proteomics-types';

const proteinIdCol = df.col('Protein IDs') ?? df.col('Majority protein IDs');
if (proteinIdCol)
  proteinIdCol.semType = SEMTYPE.PROTEIN_ID;

const geneCol = df.col('Gene names') ?? df.col('Gene name');
if (geneCol)
  geneCol.semType = SEMTYPE.GENE_SYMBOL;
```

### Pattern 5: Log2 Transformation as New Columns
**What:** Add new columns with log2-transformed values alongside originals
**When to use:** For intensity columns detected by prefix
**Example:**
```typescript
for (const col of df.columns.toList()) {
  if (isIntensityColumn(col.name)) {
    const log2Col = df.columns.addNewFloat(`log2(${col.name})`);
    for (let i = 0; i < df.rowCount; i++) {
      const val = col.get(i);
      log2Col.set(i, (val > 0) ? Math.log2(val) : DG.FLOAT_NULL);
    }
    log2Col.semType = SEMTYPE.INTENSITY;
  }
}
```

### Pattern 6: Demo Data Loading via _package.files
**What:** Load bundled demo data from the package's `files/` directory at runtime
**When to use:** For the proteomicsDemo() function
**Example:**
```typescript
// In package.ts
import {_package} from './package';

static async proteomicsDemo(): Promise<void> {
  const text = await _package.files.readAsText('demo/proteinGroups.txt');
  const df = parseMaxQuantText(text);
  grok.shell.addTableView(df);
}
```

### Pattern 7: Decorator + Metadata Comments
**What:** The scaffold uses the decorator polyfill pattern. Functions need `@grok.decorators.func()` or `@grok.decorators.demo()` decorators AND may need metadata comments for `grok api` generation.
**When to use:** All registered functions in package.ts
**Important:** The decorator polyfill (lines 13-36 of package.ts) must remain. Do NOT modify it.

### Anti-Patterns to Avoid
- **Modifying detectors.js for import:** Detectors run on generic table open, not import. The parser should assign semantic types directly -- detectors are a fallback for non-imported data.
- **Bundling the full MaxQuant output:** proteinGroups.txt can be huge. The demo must be trimmed to <5MB. Strip unnecessary columns and limit to ~500-1000 proteins.
- **Using `grok.data.loadTable(url)` for demo data:** This loads from a URL. Use `_package.files.readAsText()` which reads from the package's bundled `files/` directory.
- **Creating a custom TSV parser:** `DG.DataFrame.fromCsv()` with `{delimiter: '\t'}` handles all TSV parsing edge cases. Do not hand-roll line splitting.
- **Changing the parser function signature:** The CONTEXT specifies `parseMaxQuant(file: DG.FileInfo): Promise<DG.DataFrame>`. Additionally create a `parseMaxQuantText(text: string): DG.DataFrame` for when you already have the text content (demo data, tests).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| TSV parsing | Manual line/column splitting | `DG.DataFrame.fromCsv(text, {delimiter: '\t'})` | Handles quoting, type inference, missing values, large files |
| File picker dialog | Custom HTML input element | `DG.Utils.openFile({accept, open})` | Standard Datagrok pattern, handles browser differences |
| Table view creation | Manual DOM manipulation | `grok.shell.addTableView(df)` | Creates proper Datagrok table view with grid, toolbar, menus |
| Column type detection | Custom type sniffing | `DG.DataFrame.fromCsv()` auto-detection | Datagrok's parser infers int/float/string correctly |
| Info notification | Custom toast/modal | `grok.shell.info('message')` | Standard Datagrok balloon notification |

**Key insight:** The Datagrok API already provides all the building blocks. The parser is mostly orchestration -- parse TSV, filter rows, detect intensity columns, add log2 columns, assign semantic types.

## Common Pitfalls

### Pitfall 1: Zero/Negative Values in Log2 Transformation
**What goes wrong:** `Math.log2(0)` returns `-Infinity`, `Math.log2(-1)` returns `NaN`. These corrupt downstream analysis.
**Why it happens:** MaxQuant intensity values often contain 0 for undetected proteins.
**How to avoid:** Map zero and negative values to `DG.FLOAT_NULL` (Datagrok's null representation for float columns). Use `col.set(i, DG.FLOAT_NULL)` or check `isNaN(col.get(i))` before transformation.
**Warning signs:** Seeing `-Infinity` or extreme negative values in log2 columns.

### Pitfall 2: MaxQuant Column Name Variations
**What goes wrong:** Parser fails to find columns because of naming differences between MaxQuant versions.
**Why it happens:** MaxQuant 1.x uses "Protein IDs", 2.x may use "Protein IDs" or "Majority protein IDs" as primary. Case and whitespace can vary.
**How to avoid:** Use case-insensitive matching. Check multiple known column names for each concept. The existing `findColumn()` utility already supports name hints with `toLowerCase()`.
**Warning signs:** Semantic types not being assigned, null columns where data is expected.

### Pitfall 3: Semicolon-Delimited Multi-Value Fields
**What goes wrong:** Protein IDs like `P12345;Q67890;O11111` are treated as a single string instead of extracting the primary ID.
**Why it happens:** MaxQuant stores multiple protein accessions in semicolon-delimited format.
**How to avoid:** Parse semicolon-delimited fields and extract the first entry as the primary identifier. Optionally keep the full string in the original column and add a "Primary Protein ID" column.
**Warning signs:** UniProt lookup failing on semicolon-containing strings in later phases.

### Pitfall 4: Contaminant Prefix Mismatch
**What goes wrong:** Contaminants not filtered because the parser only checks one detection method.
**Why it happens:** MaxQuant marks contaminants two ways: (1) a `+` in the "Potential contaminant" column, AND (2) a `CON__` prefix in protein IDs. Some datasets only have one marker.
**How to avoid:** Check BOTH the dedicated column marker AND the protein ID prefix. Same for reverse hits (`REV__` prefix OR `+` in "Reverse" column).
**Warning signs:** CON__ or REV__ proteins remaining in filtered output.

### Pitfall 5: DG.FLOAT_NULL vs JavaScript NaN
**What goes wrong:** Using JavaScript `NaN` where Datagrok expects its own null representation, or vice versa.
**Why it happens:** Datagrok uses a special float null value internally.
**How to avoid:** When checking if a float column value is null, use `col.isNone(i)`. When setting null, use `col.set(i, DG.FLOAT_NULL)`. `DG.DataFrame.fromCsv()` handles this automatically for missing values in the source data.
**Warning signs:** Unexpected behavior in aggregations, filters, or visualizations with missing data.

### Pitfall 6: Large Demo File Bloating Package
**What goes wrong:** Demo proteinGroups.txt is too large, slowing package loading and publishing.
**Why it happens:** Full MaxQuant output can be 50-200MB with thousands of columns.
**How to avoid:** Curate the demo file: keep only essential columns (Protein IDs, Gene names, filtering columns, a few LFQ intensity columns), limit to ~500 proteins, target <1MB.
**Warning signs:** Package publish takes too long, demo loading is slow.

## Code Examples

### Complete Parser Flow (Recommended Implementation)
```typescript
// Source: Verified against Datagrok API (js-api/src/dataframe/, js-api/src/const.ts, js-api/src/utils.ts)

import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';

const INTENSITY_PREFIXES = ['lfq intensity', 'ibaq', 'reporter intensity', 'intensity'];
const FILTER_COLUMNS = {
  contaminant: ['Potential contaminant', 'Potential.contaminant'],
  reverse: ['Reverse'],
  onlyBySite: ['Only identified by site', 'Only.identified.by.site'],
};

/** Parses MaxQuant proteinGroups.txt from text content. */
export function parseMaxQuantText(text: string): DG.DataFrame {
  // Step 1: Parse TSV
  const raw = DG.DataFrame.fromCsv(text, {delimiter: '\t'});

  // Step 2: Filter rows (contaminants, reverse, only-by-site)
  raw.filter.init(true);  // all rows visible initially
  filterRows(raw);
  const df = raw.clone(raw.filter);
  df.name = 'proteinGroups';

  // Step 3: Assign semantic types
  assignSemTypes(df);

  // Step 4: Detect and log2-transform intensity columns
  addLog2Columns(df);

  return df;
}

function filterRows(df: DG.DataFrame): void {
  const beforeCount = df.rowCount;

  for (const [_key, names] of Object.entries(FILTER_COLUMNS)) {
    for (const name of names) {
      const col = df.col(name);
      if (col) {
        for (let i = 0; i < df.rowCount; i++) {
          if (col.get(i) === '+')
            df.filter.set(i, false);
        }
        break;
      }
    }
  }

  // Also check protein ID prefixes
  const idCol = df.col('Protein IDs') ?? df.col('Majority protein IDs');
  if (idCol) {
    for (let i = 0; i < df.rowCount; i++) {
      const id = idCol.get(i);
      if (id && (id.startsWith('CON__') || id.startsWith('REV__')))
        df.filter.set(i, false);
    }
  }
}
```

### Intensity Column Detection
```typescript
// Source: MaxQuant output documentation (cox-labs.github.io/coxdocs/output_tables.html)

function isIntensityColumn(name: string): boolean {
  const lower = name.toLowerCase();
  return INTENSITY_PREFIXES.some((prefix) => lower.startsWith(prefix));
}

function detectIntensityType(df: DG.DataFrame): string | null {
  const colNames = df.columns.names();
  if (colNames.some((n) => n.toLowerCase().startsWith('lfq intensity'))) return 'LFQ';
  if (colNames.some((n) => n.toLowerCase().startsWith('ibaq'))) return 'iBAQ';
  if (colNames.some((n) => n.toLowerCase().startsWith('reporter intensity'))) return 'Reporter';
  if (colNames.some((n) => n.toLowerCase().startsWith('intensity '))) return 'Intensity';
  return null;
}
```

### Log2 Transformation
```typescript
// Source: Verified against Datagrok API (js-api/src/dataframe/column.ts, js-api/src/const.ts)

function addLog2Columns(df: DG.DataFrame): void {
  for (const col of df.columns.toList()) {
    if (col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) {
      if (isIntensityColumn(col.name)) {
        const log2Col = df.columns.addNewFloat(`log2(${col.name})`);
        for (let i = 0; i < df.rowCount; i++) {
          if (col.isNone(i)) {
            log2Col.set(i, DG.FLOAT_NULL);
          } else {
            const val = col.get(i);
            log2Col.set(i, val > 0 ? Math.log2(val) : DG.FLOAT_NULL);
          }
        }
        log2Col.semType = SEMTYPE.INTENSITY;
      }
    }
  }
}
```

### Semantic Type Assignment
```typescript
// Source: Existing scaffold (src/utils/proteomics-types.ts, detectors.js)

function assignSemTypes(df: DG.DataFrame): void {
  // Protein IDs
  const proteinIdCol = df.col('Protein IDs') ?? df.col('Majority protein IDs');
  if (proteinIdCol) proteinIdCol.semType = SEMTYPE.PROTEIN_ID;

  // Gene names
  const geneCol = df.col('Gene names') ?? df.col('Gene name');
  if (geneCol) geneCol.semType = SEMTYPE.GENE_SYMBOL;

  // Intensity columns
  for (const col of df.columns.toList()) {
    if (isIntensityColumn(col.name) && (col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT))
      col.semType = SEMTYPE.INTENSITY;
  }
}
```

### File Import Handler
```typescript
// Source: Verified against DG.Utils.openFile (js-api/src/utils.ts lines 204-213)

static async importMaxQuant(): Promise<void> {
  DG.Utils.openFile({
    accept: '.txt,.tsv',
    open: async (file: File) => {
      const text = await file.text();
      const df = parseMaxQuantText(text);
      grok.shell.addTableView(df);
      grok.shell.info(`Imported ${df.rowCount} protein groups`);
    },
  });
}
```

### Demo Data Loading
```typescript
// Source: Common Datagrok pattern (Dendrogram, Bio packages use _package.files.readAsText)

static async proteomicsDemo(): Promise<void> {
  const text = await _package.files.readAsText('demo/proteinGroups.txt');
  const df = parseMaxQuantText(text);
  grok.shell.addTableView(df);
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `//name:` metadata comments only | `@grok.decorators.func()` decorator + metadata | Datagrok recent versions | Scaffold already uses decorator polyfill; follow this pattern |
| Manual npm link | `grok link` command | datagrok-tools update | Use `grok link` or `npm run link-all` for local dev |
| Custom file dialogs | `DG.Utils.openFile()` | Available in datagrok-api | Standard cross-browser file picker |

**Deprecated/outdated:**
- Direct TypeScript decorators without the polyfill: The scaffold includes a decorator polyfill (lines 13-36 of package.ts) because the platform may not provide `grok.decorators` yet. Do not remove or modify the polyfill.

## Open Questions

1. **Demo dataset source**
   - What we know: Need a small (<5MB) public MaxQuant proteinGroups.txt with LFQ intensity columns
   - What's unclear: Which specific dataset to use (PXD000561 HeLa or another)
   - Recommendation: Create a synthetic or heavily trimmed demo file with ~200-500 proteins, 6-8 LFQ intensity columns, representative gene names, and a few contaminant/reverse entries for testing filtering. This avoids licensing concerns and keeps the file small. Alternatively, search PRIDE for a CC0-licensed small dataset.

2. **parseMaxQuant signature vs parseMaxQuantText**
   - What we know: CONTEXT locks the signature as `parseMaxQuant(file: DG.FileInfo): Promise<DG.DataFrame>`. But for the menu handler and demo, we need to parse from text strings.
   - What's unclear: Whether `DG.FileInfo` can be used from `DG.Utils.openFile()` (it provides a browser `File`, not `DG.FileInfo`).
   - Recommendation: Create two functions: `parseMaxQuantText(text: string): DG.DataFrame` for the core logic, and `parseMaxQuant(file: DG.FileInfo): Promise<DG.DataFrame>` as a wrapper that reads file content then calls `parseMaxQuantText`. The menu handler uses `DG.Utils.openFile` + `parseMaxQuantText` directly. This also makes testing easier.

3. **Intensity column naming for log2 columns**
   - What we know: User decided to keep raw columns and add log2 alongside. Claude has discretion on naming convention.
   - Recommendation: Use `log2(Original Column Name)` format (e.g., `log2(LFQ intensity Sample1)`). This is clear, matches mathematical notation, and is distinguishable from raw columns.

## Sources

### Primary (HIGH confidence)
- Datagrok js-api source code: `js-api/src/const.ts` (CsvImportOptions with delimiter), `js-api/src/utils.ts` (DG.Utils.openFile), `js-api/src/dataframe/data-frame.ts` (DataFrame.fromCsv), `js-api/src/dataframe/column.ts` (column.semType setter)
- Existing package scaffold: `packages/Proteomics/src/` (all source files examined)
- Bio package file handler pattern: `packages/Bio/src/package.ts` (@grok.decorators.fileHandler usage)
- Dendrogram package demo data pattern: `packages/Dendrogram/src/demos/` (_package.files.readAsText usage)

### Secondary (MEDIUM confidence)
- [MaxQuant output tables documentation](https://cox-labs.github.io/coxdocs/output_tables.html) - Column names for proteinGroups.txt
- [wrProteo R package](https://rdrr.io/cran/wrProteo/man/readMaxQuantFile.html) - MaxQuant version compatibility (1.6.10.x to 2.0.x), filtering patterns
- [Emory Proteomics Core MaxQuant definitions](https://www.cores.emory.edu/eipc/_includes/documents/definitions-mmaxquant.pdf) - Column descriptions

### Tertiary (LOW confidence)
- MaxQuant version-specific column name differences: Could not find authoritative documentation listing exact differences between 1.6.x and 2.0.x column names. The wrProteo R package confirms format is "well conserved" between versions. Recommend testing with real files from both versions when available.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Datagrok API is well-documented in source code; no external libraries needed
- Architecture: HIGH - Parser pattern is straightforward; existing scaffold provides clear integration points
- Pitfalls: HIGH - Common proteomics data issues are well-known; Datagrok API behavior verified from source

**Research date:** 2026-02-28
**Valid until:** 2026-03-28 (stable domain, no fast-moving dependencies)
