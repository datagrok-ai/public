---
name: define-semantic-type-detectors
version: 0.1.0
description: |
  Add a JS-only data-shape recognizer to a Datagrok package so the
  platform auto-tags matching columns on every table import — for
  custom column kinds identified by type, name pattern, or value
  statistics, not a single regex over string content. Produces a
  populated `detectors.js` at the package root and a column whose
  Semantic-type field is set the moment a user opens a table.
  Use when asked to "auto-tag a column as my custom kind when tables
  load", "register column recognition logic in my package", or "make
  the platform classify columns by data shape".
triggers:
  - auto-tag columns on table load
  - recognize column data shape
  - register column recognition logic
  - tag columns by name and stats
  - classify columns when import runs
  - assign a custom column kind
allowed-tools:
  - Read
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# define-semantic-type-detectors

## When to use

Your package owns a custom semantic type whose membership depends on
data shape — `col.type`, `col.name`, `col.min/max`, or a regex over
sampled values (`Molecule3D`, `Magnitude`, `nucleotides`, a 4-char
`PDB_ID`). Built-in detectors don't cover it, and a declarative regex
in `package.json` `meta.semanticTypes` is too weak. For regex-only
string-content detection (e.g. `CHEMBL\d+`), use `register-identifiers`
instead — no JS file needed.

## Prerequisites

- Familiarity with the column API you'll branch on (`col.type`,
  `col.name`, `col.min`, `col.max`, `col.categories`, `col.meta.units`).

## Steps

1. **Scaffold `detectors.js` — file name and class name are load-bearing.**
   ```bash
   grok add detector <SemType>          # e.g. grok add detector Nucleotides
   ```
   Creates `<PackageRoot>/detectors.js` (no `src/`) wrapping
   `class <PackageName>PackageDetectors extends DG.Package`. Both
   names are required — loader appends `PackageDetectors` to the
   package name (`DG-FACT-255`, `256`).
   ```javascript
   class SequencePackageDetectors extends DG.Package {
     //meta.role: semTypeDetector
     //input: column col
     //output: string semType
     detectNucleotides(col) { /* fill body in step 2 */ return null; }
   }
   ```

2. **Order cheap checks first — every detector runs on every column on every import** (`DG-FACT-257`).
   Lead with `col.type` and `col.name` regex; escalate to value
   sampling only after.
   ```javascript
   //meta.role: semTypeDetector
   //input: column col
   //output: string semType
   detectMagnitude(col) {
     if ((col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) &&
         col.name.toLowerCase() === 'magnitude' &&
         0 < col.min && col.max < 10) {
       col.semType = 'Magnitude';
       return col.semType;
     }
     return null;
   }
   ```

3. **Sample, don't iterate, for shape detection on values.**
   `DG.Detector.sampleCategories(col, check, min=5, max=10, ratio=1,
   minStringLength=1)` fail-closes if `col.type !== STRING` or
   categories < min (`DG-FACT-258`, `263`).
   ```javascript
   //meta.role: semTypeDetector
   //input: column col
   //output: string semType
   detectPdb(col) {
     if (DG.Detector.sampleCategories(col, (s) =>
         s.match(/^COMPND/m) && s.match(/^END/m) &&
         (s.match(/^ATOM/m) || s.match(/^HETATM/m)), 1)) {
       col.meta.units = 'pdb';
       return 'Molecule3D';
     }
     return null;
   }
   ```

4. **Disambiguate format variants with `col.meta.units`, not parallel semTypes**
   (`DG-FACT-264`). Set `col.meta.units = '<unit>'` alongside the
   returned semType; renderers key on the base type.

5. **Tag the auto-test** — three header tags configure the
   per-detector auto-test (`DG-FACT-261`):
   - `//meta.skipTest: <reason>` — skip the test.
   - `//meta.testData: <csv-path>` — fixture CSV.
   - `//meta.testDataColumnName: <name>` — assert which column matches.
   ```javascript
   //meta.role: semTypeDetector
   //input: column col
   //output: string semType
   //meta.testData: pdb_data.csv
   //meta.testDataColumnName: Molecule3D
   detectPdb(col) { /* body from step 3 */ }
   ```

6. **Build and publish.** `detectors.js` ships as a standalone
   artifact, NOT through webpack (`DG-FACT-255`) — syntax errors only
   surface at runtime. Smoke-test after publishing.
   ```bash
   npm install
   grok check                     # exits 0
   grok publish <host>            # add --release once stable
   ```

## Common failure modes

- **Detector never fires.** File or class misnamed — must be
  `<PackageRoot>/detectors.js` and `<PackageName>PackageDetectors`
  (`DG-FACT-255`, `256`).
- **No build error, recognition broken at runtime.** Webpack syntax
  warnings don't reach `detectors.js` — check the browser console
  after publish (`DG-FACT-255`).
- **`sampleCategories` always `false`.** Column is non-string; branch
  on `col.type === DG.TYPE.STRING` first (`DG-FACT-258`).
- **All-null column gets your semType stuck.** Direct `col.categories`
  loop with regex matching empty strings; use `sampleCategories` or
  guard `if (!s) continue;` (`DG-FACT-263`).
- **Two detectors match — wrong one wins.** FIRST non-null wins, no
  platform preference (`DG-FACT-262`). Tighten cheap checks.

## See also

- Source: `help/develop/how-to/functions/define-semantic-type-detectors.md`.
- Knowledge: `DG-FACT-255`–`264`.
- Related skills: `register-identifiers`, `context-actions`.
