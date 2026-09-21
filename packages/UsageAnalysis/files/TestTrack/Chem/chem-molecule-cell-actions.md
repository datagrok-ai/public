---
feature: chem
realizes_atlas:
  - chem.cp.molecule-cell-actions
realizes:
  - chem.cp.molecule-cell-actions
priority: p0
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
produced_from: atlas-driven
realized_as:
  - chem-molecule-cell-actions-spec.ts
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-08-20-chem-new-06
    timestamp: "2026-08-20T00:00:00Z"
    failure_keys: []
    review_round: 1
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-23-chem-molecule-cell-actions-helper
    timestamp: 2026-08-23T19:27:10Z
    spec_runs:
      - spec: chem-molecule-cell-actions-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 36
        failure_keys: []
        run_mode: headless-cold
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-molecule-cell-actions-r3
    timestamp: 2026-08-22T12:46:37Z
    failure_keys: []
expected_results:
  - anchor: Scenario 1 Step 1
    expectation: >-
      The molecule cell context menu lists Copy as SMILES, Copy as MOLFILE
      V2000, Copy as MOLFILE V3000, Copy as SMARTS, Copy as Image, Export as
      SVG, and Sort by similarity.
  - anchor: Scenario 1 Step 4-5
    expectation: >-
      The "Smiles copied to clipboard" balloon appears and the clipboard SMILES
      is non-empty, carries no newline, and round-trips to the same canonical
      molecule as the cell.
  - anchor: Scenario 1 Step 6
    expectation: >-
      The clipboard V2000 molfile begins with a blank header line, contains
      "V2000" and "M  END", and carries no "V3000" tag.
  - anchor: Scenario 1 Step 7
    expectation: >-
      The clipboard V3000 molfile contains "V3000" and "M  END" and differs from
      the V2000 clipboard text.
  - anchor: Scenario 1 Step 8
    expectation: >-
      The clipboard SMARTS contains a SMARTS token (bracket atom, "#", or bond
      primitive) and differs from the copied SMILES.
  - anchor: Scenario 1 Step 9
    expectation: >-
      The "Image copied to clipboard" balloon appears and the clipboard holds an
      image/png blob larger than 500 bytes.
  - anchor: Scenario 1 Step 10
    expectation: >-
      Export as SVG downloads a file named molecule.svg whose root element is
      <svg> — RDKit's renderer emits an XML prolog first, so the document does
      not literally begin with "<svg" — and which holds at least one drawing
      primitive (path, line, or circle), with no "Failed to export structure as
      SVG" warning.
  - anchor: Scenario 2 Step 1-2
    expectation: >-
      After Sort by similarity the query molecule (the right-clicked row) leads
      the grid at row 0 and is pinned as the sort anchor on the molecule column.
      The sort also leaves a hidden Morgan-fingerprint working column
      "~canonical_smiles.Morgan" on the table (excluded from CSV and binary
      export); no similarity score column is added to the table.
  - anchor: Scenario 2 Step 3 (scores)
    expectation: >-
      Tanimoto/Morgan similarity to the query, computed independently for the
      molecules at sorted rows 0-4, is non-increasing down the grid and equals
      1.0 at row 0.
  - anchor: Scenario 2 Step 3 (partial)
    expectation: >-
      The four clipboard texts — SMILES, MOLFILE V2000, MOLFILE V3000, SMARTS —
      are pairwise distinct.
---

# Chem — Molecule Cell Actions

## Setup

1. Open the demo dataset `System:DemoFiles/chem/smiles.csv` using "Open test data" from the TestTrack
   toolbar (or via File > Open and navigating to DemoFiles/chem/smiles.csv). Wait for the grid to
   render — molecule cells show structural drawings, not raw text.

2. Identify a molecule cell in the grid's canonical_smiles column whose structure is fully rendered
   (non-empty, non-null). Note the row index of that cell (e.g. row 0).

## Scenarios

### Scenario 1: Copy-as and Export actions on a molecule cell

Steps:
1. Right-click the molecule cell at row 0 in the canonical_smiles column. The context menu
   opens; confirm the menu contains the entries "Copy as SMILES", "Copy as MOLFILE V2000",
   "Copy as MOLFILE V3000", "Copy as SMARTS", "Copy as Image", "Export as SVG", and
   "Sort by similarity".

2. Close the context menu without selecting anything (press Escape).

3. Right-click the same molecule cell again to reopen the context menu.

4. Select "Copy as SMILES" from the context menu. Observe that the Datagrok info notification
   "Smiles copied to clipboard" appears in the status bar.

5. Read the clipboard text. Verify that it is a valid SMILES string: non-empty, contains no
   embedded newlines, and a round-trip through SMILES → MolBlock → canonical SMILES restores
   the same canonical SMILES.

6. Right-click the same molecule cell and select "Copy as MOLFILE V2000". Read the clipboard.
   Verify that the text starts with a blank-line header and contains "V2000" and "M  END".

7. Right-click the same molecule cell and select "Copy as MOLFILE V3000". Read the clipboard.
   Verify that the text contains "V3000" and "M  END", and that it differs from the V2000 text
   copied in the previous step.

8. Right-click the same molecule cell and select "Copy as SMARTS". Read the clipboard. Verify
   that the text contains at least one SMARTS token (square-bracket atom descriptor, "#", or
   a bond primitive such as "-", "=", "#") and is character-for-character different from the
   SMILES copied in Step 5.

9. Right-click the same molecule cell and select "Copy as Image". Observe that the Datagrok
   info notification "Image copied to clipboard" appears. Verify that the clipboard contains
   an image/png blob larger than 500 bytes.

10. Right-click the same molecule cell and select "Export as SVG". Observe that a file named
    "molecule.svg" is downloaded. Read the downloaded file content and verify its root element
    is `<svg>` — the renderer emits an XML prolog ahead of it, so the file does not literally
    begin with "<svg" — and that it contains at least one drawing primitive element (a `<path>`,
    `<line>`, or `<circle>` tag). Confirm the "Failed to export structure as SVG" warning does NOT appear.

Expected:
- Steps 5–8: each copy action places correctly formatted text in the clipboard; the four formats
  are pairwise distinct and all parse back to the same molecule identity.
- Step 9: the clipboard holds a non-trivial PNG image blob.
- Step 10: the downloaded SVG file is a valid SVG document with drawing content.

### Scenario 2: Sort by similarity reorders the grid

Steps:
1. Reload the dataset (close and reopen smiles.csv) so the grid is in its original order. Note
   the SMILES value of the molecule at row 2 in the canonical_smiles column — this is the query
   molecule.

2. Right-click the molecule cell at row 2 and select "Sort by similarity" from the context menu.
   Wait for the sort operation to complete (the "Sort by similarity" info notification appears
   or the grid order visibly changes).

3. Read the SMILES value now shown at row 0. Verify it equals the SMILES noted in Step 1 —
   confirming the query molecule was placed first. Then check the ordering itself: the sort
   attaches no score column to the table, so take the molecules now shown at rows 0–4, compute
   each one's Tanimoto/Morgan similarity to the query molecule (Chem | Similarities, or
   `grok.chem.getSimilarities`), and confirm the five values are non-increasing down the grid
   with row 0's value equal to 1.0.

4. Verify pairwise distinctness of the four format representations obtained in Scenario 1
   Steps 5–8 for the same molecule: SMILES != MOLFILE V2000 text, SMILES != MOLFILE V3000
   text, SMILES != SMARTS text, MOLFILE V2000 text != MOLFILE V3000 text.

Expected:
- The grid row order changes so the query molecule (row 2 before sort) leads at row 0 with
  similarity score 1.0.
- Similarity to the query, computed independently for the molecules at rows 0 through 4, is
  non-increasing.
- The table gains the hidden fingerprint-cache column `~canonical_smiles.Morgan` (not shown in the
  grid and excluded from CSV export), and gains no similarity score column.
- The four clipboard-format texts from Scenario 1 are pairwise distinct.

## Notes

- The "Copy as Image" assertion uses a byte-length floor (> 500 bytes) rather than pixel
  comparison — the rendered canvas may vary with DPI; a non-trivial blob length is sufficient
  to confirm the render completed without a null or empty result.
- The similarity monotonicity check requires reading at least 5 rows to distinguish a
  real sort from a no-op; confirming only that row 0 equals 1.0 is not sufficient on its own.
  It is checked on similarity values computed independently from the grid's row order, because
  the product attaches no score column (see the next note).
- What Sort by similarity leaves on the table (checked against source on 2026-08-22, after two
  earlier versions of this note got it wrong in opposite directions):
  - A hidden Morgan-fingerprint column IS added to the table being sorted. `sortBySimilarity`
    (`public/packages/Chem/src/package.ts#L2101-2131`) fingerprints the molecule column through
    `callChemSimilaritySearch` → `chemSimilaritySearch` → `chemSimilarityBitset`
    (`src/analysis/chem-similarity-viewer.ts#L286-307`) → `chemGetFingerprints`
    (`src/chem-searches.ts#L270-272`) → `getUint8ArrayFingerprints` → `invalidateAndSaveColumns`
    → `saveColumn` (`src/chem-searches.ts#L178-190`), which does
    `col.dataFrame.columns.getOrCreate('~' + col.name + '.' + tagName, BYTE_ARRAY)` and then sets
    `includeInCsvExport = false` and `includeInBinaryExport = false`. For this dataset the column
    is named `~canonical_smiles.Morgan` (the tag is the `Fingerprint.Morgan` enum value,
    `src/utils/chem-common.ts#L42-50`). It is a cache: a second sort on the same table reuses it.
    The sort path passes `returnSmiles = false`, so no `~canonical_smiles.canonicalSmiles` column
    appears alongside it.
  - The similarity SCORE is NOT added to the table. `similarityResultDf`
    (`src/analysis/chem-similarity-viewer.ts#L311-334`) builds `smiles` / `score` / `indexes` as
    three fresh columns of a SEPARATE DataFrame, and `sortBySimilarity` consumes only `indexes`,
    via `grid.setRowOrder` (`package.ts#L2125`). That is why Scenario 2 Step 3 recomputes the
    similarities instead of reading a column.
  - `col.temp` is not where the fingerprints live — it only carries the change-tracking flags
    (`src/chem-searches.ts#L110-134`).
- The "Convert mixture to smiles..." action is also available in the molecule cell context menu
  but targets columns with chemical mixture semantic type, not ordinary molecule columns.
  The smiles.csv dataset has a Molecule column, so that action does not appear and is not
  tested here.

## Automation notes

- All seven cell actions are triggered through the grid cell context menu (right-click opens a
  DOM menu whose items are dispatched as platform actions). Clipboard read, file download, and
  sort-result assertions all require a real browser context. Clipboard text is read via
  `page.evaluate(() => navigator.clipboard.readText())`; file downloads are intercepted
  natively by the test runner.
- Source: public/packages/Chem/src/package.ts#L2190 — copyAsSmiles and sibling
  copyAsMolfileV2000 / copyAsMolfileV3000 / copyAsSmarts / copyAsImage / exportAsSvg /
  sortBySimilarity static methods on PackageFunctions, all registered with meta.action keys
  so they appear in the Datagrok cell context menu for semType=Molecule columns.
