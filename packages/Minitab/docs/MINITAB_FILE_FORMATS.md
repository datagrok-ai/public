# Minitab File Formats — Technical Specification

> **Purpose:** Developer reference for parsing, generating, and converting Minitab file formats within the Datagrok Minitab plugin. This document contains everything needed to implement file I/O without reference to Minitab itself.

---

## Overview of Minitab File Types

| Extension | Name | Generation | Structure | Priority |
|---|---|---|---|---|
| `.mpx` | Minitab Project | Minitab 19+ (2019–current) | ZIP / JSON | **P0** |
| `.mwx` | Minitab Worksheet | Minitab 19+ (2019–current) | ZIP / JSON | **P0** |
| `.mpj` | Minitab Project (legacy) | Minitab 15–18 | OLE2 / Binary | P2 |
| `.mtw` | Minitab Worksheet (legacy) | Up to Minitab 18 | OLE2 / Binary | P2 |
| `.mgf` | Minitab Graph | All versions | Binary (proprietary) | P3 |
| `.mtb` | Minitab Exec / Script | All versions | Plain text | P2 |
| `.mac` | Minitab Macro | All versions | Plain text | P2 |

---

## `.mpx` — Minitab Project File (Current Format)

### Container Structure

An `.mpx` file is a **standard ZIP archive**. Rename to `.zip` and unzip to inspect. Typical internal layout:

```
project.mpx (as ZIP)
├── project_metadata.json
├── sheets/
│   ├── 0/
│   │   └── sheet.json          ← Worksheet 0 data and column definitions
│   ├── 1/
│   │   └── sheet.json          ← Worksheet 1 (if multiple worksheets in project)
│   └── ...
├── commands/
│   ├── 1/
│   │   ├── command.json        ← Analysis command run by user
│   │   └── groups/
│   │       ├── 1/
│   │       │   └── group.json  ← Output group (e.g., session text block)
│   │       ├── 2/
│   │       │   └── group.json  ← Output group (e.g., graph)
│   │       └── ...
│   ├── 2/
│   │   └── command.json
│   └── ...
└── [optional: embedded graph/image data in subdirectories]
```

### `project_metadata.json` Schema

```json
{
  "version": "19.2",
  "minitabVersion": "20.3.0",
  "name": "Process Validation Study Q3",
  "createdDate": "2024-09-15T10:30:00Z",
  "modifiedDate": "2024-11-02T14:22:11Z",
  "activeSheetIndex": 0
}
```

Key fields:
- `version`: MPX format version (useful for forward-compat guards)
- `minitabVersion`: Minitab product version that created the file
- `name`: Project name (display in Datagrok)
- `activeSheetIndex`: Which worksheet was active when saved

### `sheet.json` Schema

This is the most important structure — it defines the worksheet data.

```json
{
  "name": "Batch Data",
  "description": "Q3 2024 tablet dissolution batches",
  "columns": [
    {
      "id": 1,
      "name": "Batch",
      "type": "Text",
      "description": "Batch identifier",
      "values": ["B001", "B001", "B002", "B002", "*", "B003"],
      "valueLabels": null
    },
    {
      "id": 2,
      "name": "Dissolution_%",
      "type": "Numeric",
      "description": "Dissolution at 45 min (%)",
      "values": [88.2, 91.4, 85.7, 89.3, "*", 92.1],
      "valueLabels": null,
      "formula": null,
      "dataDisplayType": "Numeric"
    },
    {
      "id": 3,
      "name": "Date",
      "type": "DateTime",
      "description": "",
      "values": ["2024-09-01T00:00:00", "2024-09-01T00:00:00", "*"],
      "format": "yyyy-MM-dd"
    }
  ],
  "constants": [],
  "matrices": [],
  "designObjects": []
}
```

**Column type mapping:**

| Minitab Type | JSON `type` value | Datagrok Type |
|---|---|---|
| Numeric | `"Numeric"` | `float64` |
| Text | `"Text"` | `string` |
| Date/Time | `"DateTime"` | `datetime` |

**Missing value convention:**
- Numeric: `"*"` (string in JSON array) → `null` / `NaN`
- Text: `""` (empty string) or `"*"` → `null`
- DateTime: `"*"` → `null`

**Constants:** Scalar values stored in the worksheet (Minitab `K1`, `K2`, etc.) — preserve as DataFrame metadata.

**Matrices:** 2D numeric matrices stored in the worksheet — import as separate DataFrames or metadata.

**Design objects:** Experimental designs stored in the worksheet (DOE factor levels, run order) — preserve as structured metadata for DOE analysis.

### `command.json` Schema (Session History)

```json
{
  "id": 5,
  "sessionCommand": "CAPABILITY C3;\n  USPEC 100;\n  LSPEC 80;\n  TARGET 90;\n  CONF 95.",
  "menuPath": "Stat > Quality Tools > Capability Analysis > Normal...",
  "timestamp": "2024-09-15T11:42:33Z",
  "outputGroupIds": [8, 9]
}
```

Key fields:
- `sessionCommand`: The raw Minitab command string — valuable for migration (can be used to re-run or transpile to R)
- `menuPath`: Human-readable navigation path that generated the command
- `timestamp`: When the analysis was run
- `outputGroupIds`: References to output groups containing results

### `group.json` Schema (Analysis Output)

```json
{
  "id": 8,
  "type": "SessionText",
  "content": "Process Capability Report for Dissolution_%\n\nLSL    Target    USL\n80     90        100\n\nCp: 1.47    CPL: 1.31    CPU: 1.63    Cpk: 1.31\nPp: 1.52    PPL: 1.35    PPU: 1.70    Ppk: 1.35\n"
}
```

```json
{
  "id": 9,
  "type": "Graph",
  "graphType": "CapabilitySixpack",
  "sheetId": 0,
  "columnRef": "Dissolution_%",
  "parameters": {
    "usl": 100,
    "lsl": 80,
    "target": 90,
    "distribution": "Normal"
  }
}
```

Output group types encountered in MPX:
- `"SessionText"` — raw text output from analysis
- `"Graph"` — graph with type identifier and parameters
- `"Table"` — structured output table (ANOVA, regression coefficients, etc.)

### Parsing Implementation (JavaScript / TypeScript)

```typescript
import JSZip from 'jszip';

interface MinitabColumn {
  id: number;
  name: string;
  type: 'Numeric' | 'Text' | 'DateTime';
  description: string;
  values: (number | string | null)[];
  formula?: string | null;
}

interface MinitabWorksheet {
  name: string;
  description: string;
  columns: MinitabColumn[];
}

async function parseMpx(file: File): Promise<{
  projectName: string;
  worksheets: MinitabWorksheet[];
  commands: any[];
}> {
  const zip = await JSZip.loadAsync(file);
  
  // 1. Read project metadata
  const metaJson = await zip.file('project_metadata.json')?.async('text');
  const metadata = metaJson ? JSON.parse(metaJson) : {};
  
  // 2. Find and parse all worksheets
  const worksheets: MinitabWorksheet[] = [];
  const sheetFiles = Object.keys(zip.files).filter(f => f.match(/^sheets\/\d+\/sheet\.json$/));
  
  for (const sheetPath of sheetFiles.sort()) {
    const sheetJson = await zip.file(sheetPath)?.async('text');
    if (!sheetJson) continue;
    const sheet = JSON.parse(sheetJson);
    
    // Normalize missing values
    const normalizedSheet: MinitabWorksheet = {
      name: sheet.name,
      description: sheet.description || '',
      columns: sheet.columns.map((col: any) => ({
        ...col,
        values: col.values.map((v: any) => (v === '*' || v === null) ? null : v)
      }))
    };
    worksheets.push(normalizedSheet);
  }
  
  // 3. Parse command history
  const commands: any[] = [];
  const cmdFiles = Object.keys(zip.files).filter(f => f.match(/^commands\/\d+\/command\.json$/));
  for (const cmdPath of cmdFiles.sort()) {
    const cmdJson = await zip.file(cmdPath)?.async('text');
    if (cmdJson) commands.push(JSON.parse(cmdJson));
  }
  
  return {
    projectName: metadata.name || file.name,
    worksheets,
    commands
  };
}

// Convert MinitabWorksheet to Datagrok DataFrame
function worksheetToDataFrame(ws: MinitabWorksheet): DG.DataFrame {
  const df = DG.DataFrame.create(ws.columns[0]?.values.length ?? 0);
  df.name = ws.name;
  
  for (const col of ws.columns) {
    let dgCol: DG.Column;
    switch (col.type) {
      case 'Numeric':
        dgCol = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, col.name, col.values as number[]);
        break;
      case 'Text':
        dgCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, col.name, col.values as string[]);
        break;
      case 'DateTime':
        dgCol = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, col.name, col.values as string[]);
        break;
      default:
        dgCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, col.name, col.values as string[]);
    }
    if (col.description) dgCol.setTag('description', col.description);
    df.columns.add(dgCol);
  }
  return df;
}
```

---

## `.mwx` — Minitab Worksheet File (Current Format)

`.mwx` files use the same ZIP/JSON structure as `.mpx` but contain only worksheet data — no project metadata, commands, or graphs.

### Internal Structure

```
worksheet.mwx (as ZIP)
└── sheets/
    └── 0/
        └── sheet.json          ← Single worksheet data
```

The `sheet.json` schema is identical to the one described above for MPX.

**Detection:** Check for presence of `sheets/0/sheet.json` without `project_metadata.json`.

---

## `.mpj` — Minitab Project File (Legacy, Minitab 15–18)

### Container Format

OLE2 Compound Document (Microsoft Compound File Binary / CFB format).

**Magic bytes (file signature):** `D0 CF 11 E0 A1 B1 1A E1`

**CLSID:** `de9268e3-30f8-11d0-9a73-00a024c92ff4`

### Internal Streams

| Stream Name | Contents |
|---|---|
| `Worksheets` | Binary-encoded worksheet data |
| `Compressed` | Compressed additional project data |
| `ProjectInfo` | Project metadata |
| `Charts` | Graph data |
| `UI` | User interface state |

### Parsing Approach

Use an OLE2/CFB library:
- JavaScript: `cfb` npm package (`npm install cfb`)
- Python: `olefile` package (`pip install olefile`)

The binary format within each stream is Minitab-proprietary and not publicly documented. The most practical parsing approach:

1. Extract the `Worksheets` stream using CFB library
2. Reverse-engineer the binary structure:
   - Look for column name strings (null-terminated or length-prefixed ASCII/UTF-16)
   - Locate arrays of float64 values (IEEE 754 little-endian)
   - Use `*` placeholder (special bit pattern or sentinel value) for missing numerics
3. Alternatively: invoke Minitab's own COM automation API to export to MWX, then parse MWX

**Recommended approach for V1:** Rather than full binary reverse-engineering, provide a conversion utility:
1. Install a lightweight Python script on the user's machine that uses `olefile` + Minitab COM automation
2. Script opens `.mpj`, exports each worksheet to CSV/MWX
3. Datagrok imports the CSV/MWX

---

## `.mtw` — Minitab Worksheet File (Legacy)

Same OLE2 structure as `.mpj`, contains only the worksheet streams. Same parsing challenges apply.

**Detection:** OLE2 magic bytes + file extension `.mtw`

---

## `.mtb` — Minitab Exec Script

### Format

Plain ASCII text file. No header, no delimiter.

```
# Monthly Process Capability Report
# Generated: 2024-03-01
# Author: QE Team

RETRIEVE 'C:\Data\BatchData_March.MPJ'
ZCHART 'Assay_%' 'Batch';
  MU 100;
  SIGMA 2.5;
  CONF 95;
  TEST.
ICHART 'Dissolution_%';
  CONF 95;
  TEST.
CAPABILITY 'Assay_%';
  USPEC 105;
  LSPEC 95;
  TARGET 100.
```

**Conventions:**
- Lines starting with `#` are comments (ignored by Minitab)
- Column names in single quotes: `'Column Name'`
- Constants referenced as `K1`, `K2`, etc.
- Subcommands indented with spaces and terminated with `.` or `;` (continuing)
- Commands terminated with `.`

### Command Language Reference (Key Commands)

| Command | Purpose | Example |
|---|---|---|
| `RETRIEVE` | Open a saved project/worksheet | `RETRIEVE 'file.mpj'` |
| `READ` | Import data into columns | `READ C1 C2 C3` |
| `LET` | Assign a value | `LET K1 = MEAN(C1)` |
| `COPY` | Copy data | `COPY C1 C2` |
| `DELETE` | Delete rows | `DELETE K1:K10 C1` |
| `ICHART` | Individuals chart (I-MR) | `ICHART C1` |
| `XCHART` | Xbar chart | `XCHART C1 C2` |
| `RCHART` | Range chart | `RCHART C1 C2` |
| `PCHART` | P chart | `PCHART C1 C2` |
| `CCHART` | C chart | `CCHART C1` |
| `UCHART` | U chart | `UCHART C1 C2` |
| `CUSUMCHART` | CUSUM chart | `CUSUMCHART C1` |
| `EWMACHART` | EWMA chart | `EWMACHART C1` |
| `CAPABILITY` | Process capability analysis | `CAPABILITY C1` + subcommands |
| `GRREPEAT` | Gage R&R (crossed) | `GRREPEAT C1 C2 C3` |
| `REGRESS` | Regression | `REGRESS C1 1 C2` |
| `ANOVA` | Analysis of variance | `ANOVA C1 = C2|C3` |
| `TTEST` | One-sample t-test | `TTEST 0 C1` |
| `TWOSAMPLE` | Two-sample t-test | `TWOSAMPLE C1 C2` |
| `HISTOGRAM` | Histogram | `HISTOGRAM C1` |
| `BOXPLOT` | Boxplot | `BOXPLOT C1` |
| `PLOT` | Scatterplot | `PLOT C1*C2` |
| `NOTE` | Write to session window | `NOTE Analysis complete` |
| `OUTFILE` | Redirect output to file | `OUTFILE 'output.txt'` |
| `GMACRO` | Begin global macro | `GMACRO\nMacroName` |
| `ENDMACRO` | End macro | `ENDMACRO` |

### Subcommands (Common Examples)

```
CAPABILITY C1;
  USPEC 105;        # Upper spec limit
  LSPEC 95;         # Lower spec limit
  TARGET 100;       # Target value
  CONF 95;          # Confidence level for CI
  WITHIN;           # Use within-subgroup variation (Cp/Cpk)
  OVERALL.          # Use overall variation (Pp/Ppk)

ICHART C1;
  TEST 1 2 3 4 5 6 7 8;   # Apply specified Western Electric rules
  CONF 95.

REGRESS C1 2 C2 C3;
  BRIEF 2;          # Amount of output: 1=minimal, 2=default, 3=verbose
  RESIDUALS C4;     # Store residuals in C4
  FITS C5.          # Store fitted values in C5
```

---

## `.mac` — Minitab Macro

### Global Macro (`.mac`)

```
GMACRO
CapReport
# Arguments: column to analyze, USL, LSL
MCOLUMN Col
MCONSTANT USL LSL

CAPABILITY 'Col';
  USPEC USL;
  LSPEC LSL;
  TARGET (USL+LSL)/2.

ENDMACRO
```

Invoked in Minitab: `%CapReport C1 105 95`

### Local Macro (`.mac`)

```
MACRO
Analyze x y z
MCOLUMN x y z
MCONSTANT N i

LET N = COUNT(x)
DO i = 1:N
  LET y(i) = x(i) * 2
  LET z(i) = SQRT(x(i))
ENDDO
ENDMACRO
```

### Control Flow Keywords

| Keyword | Purpose |
|---|---|
| `IF ... ELSEIF ... ELSE ... ENDIF` | Conditional |
| `DO ... ENDDO` | Count-controlled loop |
| `WHILE ... ENDWHILE` | Condition-controlled loop |
| `CALL macro args` | Invoke another macro |
| `RETURN` | Exit macro |
| `BREAK` | Exit loop |
| `NEXT` | Continue to next loop iteration |

---

## Column and Data Type Details

### Numeric Columns

- Stored as IEEE 754 double-precision float (64-bit) internally
- Missing value: `*` in display, special sentinel in binary formats
- Column subtype: `Continuous` (default) or `Discrete`

### Text Columns

- Variable-length strings
- Missing value: blank (`""`)

### Date/Time Columns

- Internally stored as float (days since an epoch, similar to Excel)
- Display format configurable (e.g., `MM/DD/YYYY`, `YYYY-MM-DD HH:MM:SS`)
- Common formats used in pharma: `MM/DD/YYYY` (US), `DD/MM/YYYY` (EU)

### Value Labels (Coded Text)

Numeric columns can have value labels — e.g., `1 = "Pass"`, `2 = "Fail"`. These are stored as a mapping in column metadata. Preserve as Datagrok column categories or tags.

### Column Formulas

Columns can have formulas (like Excel): `C3 = C1 + C2`. Stored as expression strings. Datagrok can re-evaluate these after import using its formula engine.

---

## Recommended NPM Packages

| Purpose | Package | Notes |
|---|---|---|
| ZIP parsing (MPX/MWX) | `jszip` | Browser-compatible; well-maintained |
| OLE2 parsing (MPJ/MTW) | `cfb` | Pure JS; handles compound files |
| JSON schema validation | `zod` | Type-safe parsing of sheet.json |
| Date parsing | `date-fns` | For DateTime column values |

### Detection Logic

```typescript
async function detectMinitabFormat(file: File): Promise<'mpx' | 'mwx' | 'mpj' | 'mtw' | 'mtb' | 'mac' | 'unknown'> {
  const ext = file.name.split('.').pop()?.toLowerCase();
  
  if (ext === 'mtb' || ext === 'mac') return ext;
  
  // Read first 8 bytes for magic number detection
  const header = await file.slice(0, 8).arrayBuffer();
  const bytes = new Uint8Array(header);
  
  // OLE2: D0 CF 11 E0 A1 B1 1A E1
  const isOle2 = bytes[0] === 0xD0 && bytes[1] === 0xCF && bytes[2] === 0x11 && bytes[3] === 0xE0;
  if (isOle2) return ext === 'mpj' ? 'mpj' : 'mtw';
  
  // ZIP: PK\x03\x04
  const isZip = bytes[0] === 0x50 && bytes[1] === 0x4B && bytes[2] === 0x03 && bytes[3] === 0x04;
  if (isZip) {
    // Distinguish MPX from MWX by presence of project_metadata.json
    const zip = await JSZip.loadAsync(file);
    const hasMeta = !!zip.file('project_metadata.json');
    return hasMeta ? 'mpx' : 'mwx';
  }
  
  return 'unknown';
}
```

---

## Writing `.mwx` (Export from Datagrok to Minitab)

To generate a valid `.mwx` file that Minitab 19+ can open:

```typescript
async function dataFrameToMwx(df: DG.DataFrame): Promise<Blob> {
  const zip = new JSZip();
  
  // Build sheet.json
  const columns = [];
  for (let i = 0; i < df.columns.length; i++) {
    const col = df.columns.byIndex(i);
    const values = [];
    for (let r = 0; r < df.rowCount; r++) {
      const val = col.get(r);
      // Represent missing values as Minitab convention
      if (val === null || val === undefined || (typeof val === 'number' && isNaN(val))) {
        values.push('*');
      } else {
        values.push(val);
      }
    }
    columns.push({
      id: i + 1,
      name: col.name,
      type: colTypeToMinitab(col.type),
      description: col.getTag('description') ?? '',
      values
    });
  }
  
  const sheetJson = JSON.stringify({
    name: df.name || 'Sheet',
    description: '',
    columns,
    constants: [],
    matrices: [],
    designObjects: []
  }, null, 2);
  
  zip.folder('sheets')?.folder('0')?.file('sheet.json', sheetJson);
  
  return await zip.generateAsync({ type: 'blob', compression: 'DEFLATE' });
}

function colTypeToMinitab(dgType: string): string {
  switch (dgType) {
    case DG.COLUMN_TYPE.FLOAT:
    case DG.COLUMN_TYPE.INT: return 'Numeric';
    case DG.COLUMN_TYPE.DATE_TIME: return 'DateTime';
    default: return 'Text';
  }
}
```

---

## Known Limitations and Edge Cases

1. **Multi-worksheet projects:** `.mpx` can contain multiple worksheets (multiple `sheets/N/` directories). Always import all of them and let the user choose or keep all.

2. **Design objects:** DOE factor-level structures are stored in `designObjects` within `sheet.json`. These require specialized handling to reconstruct the DOE metadata needed for analysis. For V1, import the data columns normally and flag that DOE metadata was detected.

3. **Column formulas:** Columns with active formulas should have the formula re-evaluated or at minimum the stored computed values imported. Flag formula columns in Datagrok column metadata.

4. **Encoded graphs in MPX:** Some MPX files embed graph image data (SVG, PNG) in subdirectories. For V1, skip these — offer to re-create analyses natively. For V2, render the embedded images in a read-only viewer.

5. **Legacy MPJ/MTW binary:** Minitab's internal binary format within OLE2 streams has changed across versions 13–18. A robust parser may need version-specific handling. The `.mpj`/`.mtw` parser should fail gracefully and prompt the user to save as `.mwx` in Minitab first if parsing fails.

6. **Character encoding:** Minitab stores text in UTF-8 for MPX/MWX. Legacy MPJ/MTW may use Windows ANSI code pages (CP1252 for Western Europe, CP932 for Japanese). Handle encoding detection when parsing legacy formats.

7. **Row count consistency:** Minitab worksheets allow columns of different lengths (shorter columns are padded with missing values). Normalize to the maximum column length when creating a Datagrok DataFrame.

---

*See also: `MINITAB_PLUGIN.md` for full feature and workflow documentation.*
