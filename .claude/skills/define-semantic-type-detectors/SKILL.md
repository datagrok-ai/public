---
name: define-semantic-type-detectors
description: Add a JS-only semantic-type detector to a Datagrok package so the platform tags matching columns on every table import
---

# define-semantic-type-detectors

## When to use

Your package owns a custom semantic type whose membership depends on
data shape — column type, name pattern, value statistics, regex over
sampled values (`Molecule3D`, `Magnitude`, `nucleotides`, a 4-char
`PDB_ID`). Built-in detectors don't cover it, and a declarative regex
in `package.json` `meta.semanticTypes` is too weak — you need JS to
inspect `col.type`, `col.min/max`, or sample categories. Triggers:
"tag this column as `<SemType>` automatically", "register a detector
for `<SemType>`". For declarative regex-only string-content detection
(e.g. `CHEMBL\d+`), use `register-identifiers` — no JS file needed.

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- Familiarity with the column API surface you'll branch on
  (`col.type`, `col.name`, `col.min`, `col.max`, `col.categories`,
  `col.meta.units`).

## Steps

1. **Scaffold `detectors.js` with the CLI — file name and class name are load-bearing.**
   ```bash
   grok add detector <SemType>          # e.g. grok add detector Nucleotides
   ```
   Creates `<PackageRoot>/detectors.js` (no `src/`, no rename) and
   wraps it in `class <PackageName>PackageDetectors extends DG.Package`
   with a `detect<SemType>` stub. Both names are required: the loader
   resolves the class by appending `PackageDetectors` to the
   camel-cased package name (`DG-FACT-255`, `DG-FACT-256`). Expected:

   ```javascript
   class SequencePackageDetectors extends DG.Package {
     //meta.role: semTypeDetector
     //input: column col
     //output: string semType
     detectNucleotides(col) {
       if (col.name.startsWith('nuc')) {
         col.semType = 'nucleotides';
         return col.semType;
       }
       return null;
     }
   }
   ```

2. **Order the cheap checks first; a detector runs on every column on every import.**
   The platform calls every registered detector on every column each
   time a table is opened (`DG-FACT-257`). Lead with
   `col.type === DG.TYPE.STRING` (or `FLOAT` / `INT`) and a
   `col.name`-regex / startsWith — these are O(1). Only escalate to
   value sampling when the cheap checks pass. `col.min` / `col.max`
   are precomputed; using them is fine. The header-comment annotations
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
   at least `ratio * max` of them pass `check` (`DG-FACT-258`). Two
   silent fail-closed gotchas: returns `false` immediately if
   `col.type !== TYPE.STRING`, and if `categories.length < min`. The
   `minStringLength: 1` floor blocks empty-string categories
   (`DG-FACT-263`); if you iterate `col.categories` directly, you
   must `if (!s) continue;` or an all-null column gets the type stuck
   on it for the session. Canonical pattern
   (`packages/BiostructureViewer/detectors.js:13-31`):

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
   units tag for parser selection. Do NOT invent `Molecule3D-pdb`
   and `Molecule3D-pdbqt` as separate semTypes — renderers won't
   fire because they're keyed on the base type.

5. **Tag the auto-test — every detector gets one and it runs in CI.**
   Three header tags configure the per-detector auto-test
   (`DG-FACT-261`); they affect the test only, not runtime:
   - `//meta.skipTest: <reason>` — skip the test (free-form reason).
   - `//meta.testData: <csv-path>` — fixture CSV; exactly one column
     in it must match.
   - `//meta.testDataColumnName: <name>` — assert the detector picks
     that column (catches false positives on siblings).

   ```javascript
   //meta.testData: pdb_data.csv
   //meta.testDataColumnName: Molecule3D
   ```

6. **Build and publish — `detectors.js` ships separately from the bundle.**
   ```bash
   npm install
   grok check                     # exits 0
   grok publish <host>            # add --release once stable
   ```
   `detectors.js` is uploaded as a standalone artifact, NOT bundled
   through webpack (`DG-FACT-255`) — a syntax error here produces no
   build warning and only manifests at runtime as a class-load
   failure. Always smoke-test by opening a table after publishing.

## Common failure modes

- **Detector never fires.** Either the file is misnamed (must be
  exactly `<PackageRoot>/detectors.js` — no `src/`, no rename) or the
  class is misnamed (must be exactly `<PackageName>PackageDetectors`;
  `DG-FACT-255`, `DG-FACT-256`). The loader is a string match;
  a typo silently registers nothing.
- **No build error, but detection broken at runtime.** `detectors.js`
  is uploaded separately, so webpack syntax-error warnings DO NOT
  surface here (`DG-FACT-255`). Open the browser console after
  `grok publish`; class-load failures appear there, not in
  `grok check`.
- **`sampleCategories` always returns `false` for a column with
  obviously matching values.** The column is non-string —
  `sampleCategories` fail-closes silently for
  `col.type !== TYPE.STRING` (`DG-FACT-258`). Branch on
  `col.type === DG.TYPE.STRING` first, or rely on the silent guard
  to skip non-string columns.
- **All-null column gets your semType stuck on it.** The detector
  iterates `col.categories` directly and the regex matches the empty
  string. Either route through `DG.Detector.sampleCategories`
  (`minStringLength: 1` blocks empties) or guard with
  `if (!s) continue;` (`DG-FACT-263`).
- **Two detectors match — the wrong one wins.** A column gets exactly
  one `semType`; FIRST non-null wins (`DG-FACT-262`). Standard
  detectors do NOT get preference. If your detector aliases an
  existing platform semType too loosely, it shadows the stricter one.
  Tighten cheap checks first, or return `null` for cases the platform
  handles.

## Verification

- `grok check` and `grok publish <host>` exit `0`.
- Open a table whose data should match: the column's "Semantic type"
  in the Properties panel shows the returned string.
- Open a table whose data should NOT match: the column's semantic
  type is empty or set by a different (stricter) detector.
- For variant disambiguation: the column's `meta.units` tag shows
  the expected per-variant unit.
- Browser console after import — no `<Pkg>PackageDetectors` errors.

## See also

- Source articles:
  - `help/develop/how-to/functions/define-semantic-type-detectors.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-255` … `DG-FACT-264` (file/class names, `sampleCategories`
  semantics, return forms, test tags, ordering, empty-value pitfall,
  `units` tag).
- Reference packages:
  - `packages/BiostructureViewer/detectors.js:1-44` — canonical
    `sampleCategories` + regex for `Molecule3D` with `pdb` / `pdbqt`
    units split, plus a type+name+length check for `PDB_ID`.
  - `packages/Bio/detectors.js` — large-package detector class for
    Macromolecule sequences with constants and helpers above the class.
  - `tools/entity-template/sem-type-detector.js` — the template
    `grok add detector` substitutes from.
- Related skills:
  - `register-identifiers` — declarative regex-only path via
    `package.json` `meta.semanticTypes`; use when no JS logic is needed.
  - `context-actions` — actions dispatch on the `semType` your
    detector assigns.
