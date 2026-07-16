---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: []
realizes: []
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher.md
migration_date: 2026-05-11
source_text_fixes:
  - typo-smiles-csv-smiles-csv
  - keyboard-shortcut-ctrl-v-ctrl-v
  - step-numbering-glitch-1-2-3-4-5-6-7-8-10-11-9-12-13
  - 1608-paragraph-dropped
  - 2448-paragraph-split
  - step-13-backend-enumeration-scope-hard-coded
candidate_helpers:
  - helpers.playwright.chem.openSketcherFromCell
  - helpers.playwright.chem.sketcherHamburgerMenuClick
  - helpers.playwright.chem.sketcherTypeSmilesInInput
  - helpers.playwright.chem.sketcherPasteFromClipboard
  - helpers.playwright.chem.sketcherBackendSwitcherSelect
unresolved_ambiguities:
  - smiles-highlighted-csv-dataset-path
  - backend-switcher-hamburger-menu-label
  - sketcher-hamburger-menu-favorites-recent-list-selectors
  - alternate-make-these-structures-highlighted-entry-point-in-2448
  - 2448-coverage-type-smoke-vs-edge
  - pyramid-layer-bug-focused-chain-annotation-inert-after-migration
  - cross-cutting-bug-citations-not-propagated-as-related-bugs
scope_reductions:
  - id: SR-01
    check: A-STRUCT-02
    rationale: "A-STRUCT-02 carryforward (chain-level edge/perf coverage) — applies to sketcher.md only"
related_bugs: []
---

# Chem | Sketcher cell editor — Favorites / Recent / Copy / Paste round-trip

Canonical end-to-end exercise of the sketcher cell editor opened from a molecule cell on
`DemoFiles/chem/smiles.csv`. Walks the hamburger-menu **Favorites** add path, the **Recent**
list refresh after a SMILES edit, and the **Copy as SMILES** / **Copy as MOLBLOCK** clipboard
round-trip back into the molecular input field. The active sketcher backend for this scenario
is **OCL** — the platform fallback that's always available without needing an external
sketcher package.

Two bug-repro sub-blocks that used to live at the end of this scenario have moved: the #1608
block was dropped entirely, and the #2448 stereochemistry-preservation check now lives in the
companion file `sketcher-ui.md`. What remains here is the canonical sketcher modal walk —
Favorites / Recent / Copy / Paste / backend enumeration.

## Setup

1. **Provision linked dataset.** Bundled in the platform's System file shares (no external
   provisioning required): `System:DemoFiles/chem/smiles.csv` — SMILES-notation molecule
   column.
2. **Confirm Chem package is loaded** so that the sketcher cell editor opens on
   molecule-column double-click (`public/packages/Chem/src/utils/cell-editor.ts`). The
   **OCL** sketcher backend (`public/packages/Chem/src/open-chem/ocl-sketcher.ts`) is the
   platform fallback that must always be available even when no external sketcher package
   (Ketcher / Marvin / ChemDraw) is installed.
3. **Confirm clipboard access is granted in the test browser context.** The Copy as SMILES /
   Copy as MOLBLOCK steps round-trip through the system clipboard back into the sketcher's
   molecular input field via paste; playwright contexts must be granted clipboard read +
   write permissions for the host origin.
4. **No fixture consumed.** Per chain YAML `depends_on: []`. Dataset opened fresh inside the
   scenario.

## Scenarios

### Sketcher cell-editor walk on OCL backend

Canonical sketcher exercise on the active sketcher backend (OCL, the platform fallback). The
walk covers double-click open → hamburger-menu **Favorites** add → SMILES input from the
molecular field → **Recent** + **Favorites** list inspection on reopen → **Copy as SMILES** +
paste round-trip → **Copy as MOLBLOCK** + paste round-trip.

1. Open the dataset `System:DemoFiles/chem/smiles.csv`. Wait for the table view to render
   with the molecule column populated and the RDKit cell renderer applied.
2. Double-click any rendered molecule cell in the molecule column. The sketcher cell-editor
   modal opens with the picked molecule shown in the sketch canvas.
3. In the sketcher's hamburger menu, click **Favorites > Add to Favorites**. The current
   molecule is registered into the sketcher's Favorites list.
4. In the sketcher's molecular input field, type `C1CCCCC1` (cyclohexane) — type the SMILES
   directly, **not** via Ctrl+V paste. Press **Enter** to apply the SMILES to the sketch
   canvas. Click **OK**. The sketcher modal closes and the underlying cell is updated to
   the new molecule.
5. Reopen the sketcher by double-clicking the same (now updated) molecule cell. In the
   hamburger menu, inspect the **Recent** list and the **Favorites** list. Verify:
   - **Recent** contains the previously edited molecule (the cyclohexane SMILES applied at
     step 4).
   - **Favorites** contains the molecule that was added at step 3.
6. In the hamburger menu, click **Copy as SMILES**. The current sketcher molecule's SMILES
   representation is placed on the system clipboard.
7. Change the molecule in the sketcher (any edit — draw / modify / clear-then-redraw — so
   that the canvas no longer matches the SMILES on the clipboard).
8. With the molecular input field focused, press **Ctrl+V**. Press **Enter**. Verify: the
   sketch canvas shows the molecule whose SMILES was copied at step 6 (the clipboard
   paste round-trip restores the pre-edit state).
9. In the hamburger menu, click **Copy as MOLBLOCK**. The current sketcher molecule's
   MOLBLOCK (V2000) representation is placed on the system clipboard.
10. Change the molecule in the sketcher again (any edit so that the canvas no longer matches
    the MOLBLOCK on the clipboard).
11. With the molecular input field focused, press **Ctrl+V**. Verify: the sketch canvas
    shows the molecule whose MOLBLOCK was copied at step 9 (the clipboard paste round-trip
    restores that pre-edit state). Click **OK**. The sketcher modal closes and the
    underlying cell is updated.

### Backend enumeration (opportunistic)

Per chain `ui_coverage_responsibility` entry `chem-sketcher-backend-selection-hamburger-menu`:
the sketcher cell editor exposes a backend-switcher entry in its hamburger menu
("Sketcher type" / equivalent label) that lists every installed sketcher backend package
(OCL is always present; Ketcher / Marvin / ChemDraw are opportunistic — only available when
the corresponding plugin package is installed on the target test environment).

- **OCL**: required (platform fallback) — exercised in full by the previous scenario block.
- **Ketcher / Marvin / ChemDraw**: opportunistic. When the backend is installed on the
  test environment, the Automator MAY repeat steps 2-11 with that backend selected via the
  hamburger menu's backend switcher (assertion: the sketcher modal closes and reopens on
  the selected backend; the Favorites add + Recent list + Copy-as cycle + paste round-trip
  all complete on the alternate backend without console errors). When the backend is not
  installed, the assertion degrades to "the backend switcher entry is absent OR disabled —
  not silently broken".

The scenario passes if the OCL block (steps 1-11) completes cleanly; alternate-backend
re-runs are pass-extending, not pass-blocking. The Automator may surface alternate-backend
re-runs as `test.skip` / fixture-skipped branches when the host env lacks the
non-OCL backend package.

## Notes

- **Why only OCL is required.** OCL is the platform fallback and is always available, so it's
  the only backend this scenario requires. Ketcher / Marvin / ChemDraw each ship as a separate
  plugin package, so their availability is environment-dependent — the "Backend enumeration"
  block exercises them opportunistically rather than requiring them.
- **Sibling spec.** A Playwright spec already exists at `sketcher-spec.ts` ("Chem: Sketcher");
  it is extended / aligned to this scenario's Favorites + Recent + Copy + Paste +
  backend-enumeration walk.
