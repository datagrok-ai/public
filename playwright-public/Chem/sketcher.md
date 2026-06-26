---
feature: chem
sub_features_covered: [chem.sketcher, chem.sketcher.ocl, chem.sketcher.cell-editor, chem.actions.copy-smiles, chem.actions.copy-molfile-v2000, chem.actions.copy-as]
target_layer: playwright
coverage_type: regression
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
is **OCL** (`chem.sketcher.ocl`) — the platform fallback always available without an external
sketcher package per atlas.

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario (`depends_on: []`),
`classification: medium`, `pyramid_layer: bug-focused` (chain-level annotation) → migrated
`coverage_type: regression` per chain rev 2 directive (footer note (d): the #1608 bug-focused
sub-block is dropped, the #2448 sub-block is split into `sketcher-ui.md`; what remains in this
file is the canonical sketcher modal walk over Favorites / Recent / Copy / Paste / backend-
enumeration). `target_layer: playwright`, strategy `simple`. UI coverage owned
(`ui_coverage_delegated_to: null`) over `chem-sketcher-open-by-double-click`,
`chem-sketcher-favorites-add`, `chem-sketcher-recent-list`,
`chem-sketcher-molecular-input-enter`, `chem-sketcher-copy-as-smiles`,
`chem-sketcher-copy-as-molblock`, `chem-sketcher-paste-into-input`,
`chem-sketcher-backend-selection-hamburger-menu` per chain
`ui_coverage_responsibility`.

## Setup

1. **Provision linked dataset.** Bundled in the platform's System file shares (no external
   provisioning required): `System:DemoFiles/chem/smiles.csv` — SMILES-notation molecule
   column.
2. **Confirm Chem package is loaded** so that the sketcher cell editor opens on
   molecule-column double-click (atlas `chem.sketcher.cell-editor`,
   `public/packages/Chem/src/utils/cell-editor.ts`). The **OCL** sketcher backend
   (`chem.sketcher.ocl`, `public/packages/Chem/src/open-chem/ocl-sketcher.ts`) is the
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
   representation is placed on the system clipboard (atlas `chem.actions.copy-smiles`).
7. Change the molecule in the sketcher (any edit — draw / modify / clear-then-redraw — so
   that the canvas no longer matches the SMILES on the clipboard).
8. With the molecular input field focused, press **Ctrl+V**. Press **Enter**. Verify: the
   sketch canvas shows the molecule whose SMILES was copied at step 6 (the clipboard
   paste round-trip restores the pre-edit state).
9. In the hamburger menu, click **Copy as MOLBLOCK**. The current sketcher molecule's
   MOLBLOCK (V2000) representation is placed on the system clipboard (atlas
   `chem.actions.copy-molfile-v2000`).
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

- **`coverage_type: regression`** — canonical sketcher modal walk (Favorites add + Recent
  list + Copy-as cycle + paste round-trip) plus backend-enumeration loop. Not `smoke` (the
  section's smoke is `Advanced/scaffold-tree-functions.md` per chain
  `ui_coverage_plan.smoke_scenario`; sketcher's deep-dive isn't a viewer create+configure
  smoke). Not `edge` / `perf` (no specific failure-mode invariant or threshold being
  asserted after the #1608 drop + #2448 split). After the chain rev 2 source-text
  transformations (drop #1608 paragraph + split #2448 into `sketcher-ui.md`), what remains
  in this file is multi-action regression-of-the-set across the sketcher's Favorites /
  Recent / Copy / Paste surface — `regression` is the natural fit.
- **No JS API substitution.** Every entry in chain `ui_coverage_responsibility`
  (`chem-sketcher-open-by-double-click`, `chem-sketcher-favorites-add`,
  `chem-sketcher-recent-list`, `chem-sketcher-molecular-input-enter`,
  `chem-sketcher-copy-as-smiles`, `chem-sketcher-copy-as-molblock`,
  `chem-sketcher-paste-into-input`, `chem-sketcher-backend-selection-hamburger-menu`) is
  exercised via UI driving — double-click open + hamburger-menu walks + molecular input
  field typing + clipboard keyboard shortcuts. Per `pyramid_layer: bug-focused`. Invoking
  `chem.actions.copy-smiles` / `copy-molfile-v2000` (atlas role `action`) directly via
  `grok.functions.eval` would defeat the hamburger-menu + clipboard-roundtrip workflow
  this scenario asserts.
- **OCL is the only required backend per chain rev 2.** The original step 13 ("Repeat
  steps 2-9 for all available sketcher types for selection from hamburger menu") is
  preserved as a backend-enumeration loop (Scenarios > "Backend enumeration") but the
  enumeration scope is hard-coded to OCL-required + others-opportunistic per the chain
  rev 2 directive. Each non-OCL backend ships as a separate plugin package; their
  availability is environment-dependent.
- **#1608 dropped + #2448 split per chain rev 2.** The original body included two trailing
  bug-reproduction blocks under a "Check:" header — #1608 (large-substructure tooltip in
  Filter Panel sketch box) and #2448 (stereochemistry preservation on SMILES highlight).
  Per chain rev 2 footer note (d), #1608 is dropped entirely (no propagation, no bug-
  library entry); #2448 is split into a new ui-only scenario `sketcher-ui.md`
  (`target_layer: manual`, `coverage_type: smoke`) which carries the SMILES_highlighted.csv
  + Highlights + Sketch + stereo-preservation invariant. See migration report Decisions.
- **Existing sibling spec.** A test file already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher-spec.ts` (per
  `existing-test-index.yaml`: `test_name: "Chem: Sketcher"`, category `Chem`,
  `layer: playwright`). The Automator extends / aligns the existing spec to this
  migrated scenario's Favorites + Recent + Copy + Paste + backend-enumeration walk.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` + `closeAllViews` from
  `helpers-registry.yaml` cover the section harness. No registered helper currently
  abstracts the sketcher hamburger-menu walk or the clipboard-paste round-trip; candidate
  helpers surfaced in migration report Decisions.
- **Order in chain.** `order: 5` per source JSON. Independent scenario.
- **Source-text fixes silently applied during migration.** The original body had a step-
  numbering glitch (`1.` / `2.` / `3.` / `4.` / `5.` / `6.` / `7.` / `8.` / `10.` / `11.` /
  `9.` / `12.` / `13.` — gap at 9, doubled out-of-order 9/11). Renumbered cleanly to
  `1.` / `2.` / `3.` … in the migrated body. The dataset typo "smiles. csv" was
  normalized to `smiles.csv` and the keyboard-shortcut form "CTRL+V" was normalized to
  the canonical "Ctrl+V". The #1608 trailing paragraph was dropped (per chain rev 2
  directive); the #2448 trailing paragraph was lifted out into the new ui-only
  scenario `sketcher-ui.md`. See migration report Decisions § Source-text fixes.
