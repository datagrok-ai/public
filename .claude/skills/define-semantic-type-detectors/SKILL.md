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

1. **Scaffold `detectors.js` with the CLI — file name and class name are load-bearing.**
   ```bash
   grok add detector <SemType>          # e.g. grok add detector Nucleotides
   ```
   Creates `<PackageRoot>/detectors.js` (no `src/`, no rename) and
   wraps it in `class <PackageName>PackageDetectors extends DG.Package`
   with a `detect<SemType>` stub. Both names are required: the loader
   resolves the class by appending `PackageDetectors` to the
   camel-cased package name (`DG-FACT-255`, `DG-FACT-256`). Expected
   skeleton:

   ```javascript
   class SequencePackageDetectors extends DG.Package {
     //meta.role: semTypeDetector
     //input: column col
     //output: string semType
     detectNucleotides(col) { /* fill body in step 2 */ return null; }
   }
   ```

2. **Order the cheap checks first; a detector runs on every column on every import.**
   The platform calls every registered detector on every column each
   time a table is opened (`DG-FACT-257`). Lead with `col.type ===
   DG.TYPE.STRING` (or `FLOAT` / `INT`) and a `col.name`-regex /
   startsWith — O(1). Escalate to value sampling only when those pass;
   `col.min`/`col.max` are precomputed. The header-comment annotations
   ARE the registration — `.js` files have no decorator surface.

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

3. **Sample, don't iterate, when shape detection needs to look at values.**
   `DG.Detector.sampleCategories(col, check, min=5, max=10, ratio=1, minStringLength=1)`
   samples up to `max` distinct categories and returns `true` only if
   `ratio * max` pass `check` (`DG-FACT-258`). Silent fail-closed: returns
   `false` if `col.type !== TYPE.STRING` or `categories.length < min`.
   The `minStringLength: 1` floor blocks empty categories
   (`DG-FACT-263`); if you iterate `col.categories` directly, guard
   `if (!s) continue;` or an all-null column gets the type stuck on
   it. Canonical pattern (`packages/BiostructureViewer/detectors.js:13-31`):

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

4. **Disambiguate format variants with `col.meta.units`, not parallel semTypes.**
   When one semType covers several physical formats (PDB vs PDBQT
   under `Molecule3D`, SMILES vs MolBlock under `Molecule`), set
   `col.meta.units = '<unit>'` alongside the returned semType
   (`DG-FACT-264`). Renderers and `context-actions` branch on the
   units tag. Do NOT invent `Molecule3D-pdb` / `Molecule3D-pdbqt` as
   separate semTypes — renderers key on the base type and won't fire.

5. **Tag the auto-test — every detector gets one and it runs in CI.**
   Three header tags configure the per-detector auto-test
   (`DG-FACT-261`); they affect the test only, not runtime:
   - `//meta.skipTest: <reason>` — skip the test (free-form reason).
   - `//meta.testData: <csv-path>` — fixture CSV; exactly one column
     in it must match.
   - `//meta.testDataColumnName: <name>` — assert the detector picks
     that column (catches false positives on siblings).

   Applied to `detectPdb` from step 3 (pattern matches the article at
   `define-semantic-type-detectors.md:122-127`):

   ```javascript
   //meta.role: semTypeDetector
   //input: column col
   //output: string semType
   //meta.testData: pdb_data.csv
   //meta.testDataColumnName: Molecule3D
   detectPdb(col) { /* body from step 3 */ }
   ```
   On `grok test`, the harness loads `pdb_data.csv`, runs `detectPdb`
   on every column, and asserts the column named `Molecule3D` (and
   only that one) returns non-null.

6. **Build and publish — `detectors.js` ships separately from the bundle.**
   ```bash
   npm install
   grok check                     # exits 0
   grok publish <host>            # add --release once stable
   ```
   `detectors.js` is uploaded as a standalone artifact, NOT bundled
   through webpack (`DG-FACT-255`) — syntax errors here produce no
   build warning, only a runtime class-load failure. Always smoke-test
   by opening a table after publishing.

## Common failure modes

- **Detector never fires.** File misnamed (must be
  `<PackageRoot>/detectors.js` — no `src/`, no rename) or class
  misnamed (must be `<PackageName>PackageDetectors`; `DG-FACT-255`,
  `DG-FACT-256`). The loader is a string match — a typo silently
  registers nothing.
- **No build error, recognition broken at runtime.** `detectors.js`
  is uploaded separately; webpack syntax-error warnings DO NOT
  surface here (`DG-FACT-255`). Check the browser console after
  `grok publish` — class-load failures appear there, not in `grok check`.
- **`sampleCategories` always returns `false` for obviously matching
  values.** The column is non-string — `sampleCategories` fail-closes
  silently for `col.type !== TYPE.STRING` (`DG-FACT-258`). Branch on
  `col.type === DG.TYPE.STRING` first.
- **All-null column gets your semType stuck on it.** A direct
  `col.categories` loop with a regex matching empty strings. Route
  through `DG.Detector.sampleCategories` (its `minStringLength: 1`
  floor blocks empties) or guard with `if (!s) continue;` (`DG-FACT-263`).
- **Two detectors match — the wrong one wins.** FIRST non-null wins
  (`DG-FACT-262`); platform detectors get NO preference. A too-loose
  detector that aliases an existing platform semType will shadow the
  stricter one. Tighten cheap checks, or return `null` for cases the
  platform already handles.

## See also

- Source: `help/develop/how-to/functions/define-semantic-type-detectors.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-255` … `DG-FACT-264` (file/class names, `sampleCategories`
  semantics, return forms, test tags, ordering, empty-value pitfall,
  `units` tag).
- Reference packages:
  - `packages/BiostructureViewer/detectors.js:1-44` — canonical
    `sampleCategories` + regex for `Molecule3D` with `pdb`/`pdbqt`
    units split, plus a type+name+length check for `PDB_ID`.
  - `packages/Bio/detectors.js` — large-package detector class for
    Macromolecule sequences with constants and helpers.
  - `tools/entity-template/sem-type-detector.js` — `grok add detector` template.
- Related skills:
  - `register-identifiers` — declarative regex-only path via
    `package.json` `meta.semanticTypes`; no JS file needed.
  - `context-actions` — dispatches on the `semType` your detector assigns.
