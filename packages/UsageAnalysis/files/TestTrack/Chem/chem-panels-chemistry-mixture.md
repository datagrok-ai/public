---
feature: chem
realizes_atlas:
  - chem.cp.panels-chemistry-mixture
realizes:
  - chem.cp.panels-chemistry-mixture
priority: p1
target_layer: playwright
coverage_type: regression
pyramid_layer: integration
produced_from: atlas-driven
realized_as:
  - chem-panels-chemistry-mixture-spec.ts
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-panels-chemistry-mixture-r2
    timestamp: 2026-08-22T13:42:07Z
    spec_runs:
      - spec: chem-panels-chemistry-mixture-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 73
        failure_keys: []
        run_mode: headless-cold
  a:
    verdict: PASS
    cycle_id: 2026-08-21-chem-new-02
    timestamp: "2026-08-21T12:00:00Z"
    failure_keys: []
    review_round: 1
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-panels-chemistry-mixture-r3
    timestamp: 2026-08-22T14:01:20Z
    failure_keys: []
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      The Chemistry | Mixture panel renders its component TABLE. Mixture draws a
      table and has no per-component panes of its own — those belong to
      MixtureTree and are asserted at Step 5 — so this step measures only the
      panel it names. All panel headers are read with the Chemistry GROUP pane
      expanded: a collapsed group renders none of its children, so a header list
      of ["Actions","Chemistry"] is the correct pre-expansion state and not
      evidence of a missing panel.
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      The Chemistry | MixtureTree context panel renders the mixfile version
      header and one accordion pane per component of the mixture; the pane count
      is greater than one for the multi- component row.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      For a molecular cell the expanded Chemistry group OFFERS the molecular
      panels and does NOT offer Mixture / MixtureTree. There is no
      "single-component Mixture" state to observe: those panels' parameters are
      typed ChemicalMixture, so for a Molecule cell they are not offered at all
      — the differentiation is in WHICH panels appear, not in how Mixture
      renders. The absence is asserted only alongside a positive control in the
      same expanded group (Gasteiger Partial Charges present), so that a group
      rendering nothing at all cannot satisfy the claim.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      Expanding the Chemistry | Gasteiger Partial Charges panel for a SMILES
      cell produces a GRAPHIC, not merely an instantiated panel: the graphics
      host is present with non-zero width and height on screen, its base64 PNG
      payload decodes to non-zero pixel dimensions, and the decoded image
      carries non-white paint (a non-empty floor — it proves the figure is not
      blank, not that a contour was drawn). The script also runs without a
      non-ambient console error, and the console channel is proved live by a
      deliberate error probe before any step relies on its silence.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      Repeating Step 3 on a molV2000 cell produces a graphic with the same
      properties — present host with non-zero on-screen box, decoded PNG of
      non-zero dimensions, non-white paint — and, crucially, a PNG payload whose
      hash DIFFERS from the SMILES render's. That difference is what makes the
      dual-format claim about the product's own output: identical bytes would
      mean the panel is still showing the previous binding. The panel also runs
      without a non-ambient console error.
---

# Chem — Chemistry Mixture and Gasteiger Partial Charges context panels

## Setup

1. Open the table `System:AppData/Chem/test_mixtures.csv` via the Files browser
   (Browse > Files > App Data > Chem > test_mixtures.csv). The file ships in the
   Chem package (`public/packages/Chem/files/test_mixtures.csv`) and has no copy
   under Demo Files. The table contains a `mixture` column detected as semType
   ChemicalMixture.
2. Open a second table `System:DemoFiles/chem/smiles.csv` via the Files browser.
   The table contains a `smiles` column detected as semType Molecule.

## Scenarios

### Scenario 1: Mixture and MixtureTree panels differentiate multi-component from single-component input

Steps:
1. In the `test_mixtures.csv` table view, click the first cell of the `mixture`
   column to select it. The Context Panel opens on the right side of the screen.
2. In the Context Panel, navigate to the Chemistry tab.
3. Verify that the Chemistry tab contains both the Mixture panel and the
   MixtureTree panel sub-sections.
   Note that Chemistry is a GROUP pane: while it is collapsed none of its children
   are rendered, so a header list of `["Actions","Chemistry"]` is the correct state
   before expansion and does not mean the panels are missing. Expand the group
   before reading which panels are offered.
4. Expand the Mixture sub-section and observe its content.
   The panel renders a component TABLE. It has no per-component panes of its own —
   those belong to MixtureTree and are checked in Step 5 — so this step verifies
   only that the Mixture panel is present and draws its table.
5. Expand the MixtureTree sub-section and observe its content.
   The panel must show the mixfileVersion header and render one collapsible pane
   per mixture component; the total pane count must be greater than one.
6. Switch to the `smiles.csv` table view and click the first cell of the `smiles`
   column. In the Context Panel, expand the Chemistry group and read which panels
   it offers. For a molecular cell the group must offer the molecular panels
   (Gasteiger Partial Charges among them) and must NOT offer Mixture or
   MixtureTree at all: those panels take a ChemicalMixture parameter, so there is
   no "single-component Mixture" rendering to observe. The differentiation is in
   WHICH panels are offered, not in how Mixture draws a lone molecule.

Expected:
- The Chemistry | Mixture panel for a multi-component mixture row is present and
  renders its component table.
- The Chemistry | MixtureTree panel shows the mixfileVersion line and one
  accordion pane per component with count > 1.
- For a plain SMILES cell the expanded Chemistry group offers the molecular
  panels and offers neither Mixture nor MixtureTree — confirming the panel set is
  chosen from the cell's semantic type. The absence is only meaningful next to the
  positive control that molecular panels ARE offered in the same expanded group.

### Scenario 2: Gasteiger Partial Charges panel renders without script error

Steps:
1. In the `smiles.csv` table view, click the first cell of the `smiles` column
   to select it. The Context Panel opens on the right side of the screen.
2. In the Context Panel, navigate to the Chemistry tab.
3. Locate the Gasteiger Partial Charges sub-section and expand it.
   The Python script `Chemistry | Gasteiger Partial Charges` runs for the selected
   SMILES and renders its charge contour graphic. Verify the graphic element is
   present and has non-zero dimensions, and confirm no script error appears in
   the browser console during or after rendering.
4. Click a molV2000-format cell in `System:AppData/Chem/mol1K.sdf` (open via
   Browse > Files > App Data > Chem > mol1K.sdf) and repeat the observation of
   Step 3 for the molblock input.
   The Gasteiger Partial Charges graphic must render without a console script
   error for the molV2000 format as well.

Expected:
- The Gasteiger Partial Charges panel renders its charge contour graphic for
  a SMILES cell; the graphic DOM element is present with non-zero dimensions
  and the browser console shows no script error from the panel script.
- The panel also renders for a molV2000 cell without a console error,
  confirming the script's dual-format branch (MolFromMolBlock vs MolFromSmiles)
  executes without fault.

## Notes

- target_layer rationale: both scenarios require clicking grid cells to trigger
  context-panel rendering and then inspecting DOM sub-sections (panel pane count,
  graphic element presence and dimensions, console error absence). These are
  DOM-driven assertions that playwright handles natively; apitest cannot observe
  context-panel widget rendering.
- Gasteiger Partial Charges is a Python script panel (`scripts/gasteiger_charges.py`,
  registered as `Chemistry | Gasteiger Partial Charges`). Step 4 of Scenario 2 needs
  `mol1K.sdf`, which lives under `System:AppData/Chem/` — `System:DemoFiles/chem/mol1K.sdf`
  is recorded as absent in `chem-import-export-formats-spec.ts#L181`.
- derived_from: public/packages/Chem/src/package.ts#L3087
  (mixtureWidget / mixtureTreeWidget panel registrations);
  public/packages/Chem/scripts/gasteiger_charges.py (Gasteiger Partial Charges script).
