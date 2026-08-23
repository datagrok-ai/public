---
feature: chem
realizes_atlas:
  - chem.cp.molecule-rendering-end-to-end
  - chem.int.render-feeds-search
realizes: []
priority: p0
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
related_bugs:
  - id: GROK-16870
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      The Molecule column region of the grid canvas carries coloured ink — RDKit
      draws heteroatoms in colour (O red, N blue) where a raw-SMILES text
      fallback is monochrome — so the cells are drawn as structures rather than
      as text, the region is painted rather than blank, and the monochrome
      molregno column of the same grid, read the same way in the same run,
      carries no coloured ink, so the reading is structure ink and not grid
      chrome. No console errors appear during rendering.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      Hovering over a Box Plot data point raises a tooltip whose molecule region
      is painted (a non-blank canvas inside the tooltip), the Box Plot survives
      the hover, and the browser console shows no errors from
      rdkit-cell-renderer.ts during hover.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      Hovering over a Scatter Plot data point raises a tooltip whose molecule
      region is painted, the viewer survives the hover, and no console errors
      appear.
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      The reaction column is detected as semType ChemicalReaction, the
      ChemicalReaction cell renderer is bound to it, and the region of the grid
      canvas occupied by the reaction column is painted rather than blank. No
      console errors appear. (Arrow layout — reactants left, products right — is
      not checked; nothing in the automated step can see where the arrow sits.)
  - anchor: "Scenario 3 Step 6"
    expectation: >-
      The mixture column is detected as semType ChemicalMixture, the
      ChemicalMixture cell renderer is bound to it, every cell parses into at
      least one component and at least one cell carries more than one component,
      and the region of the grid canvas occupied by the mixture column is
      painted rather than blank. No console errors appear.
realized_as:
  - molecule-rendering-spec.ts
scope_reductions:
  - id: SR-01
    check: F-STRUCT-COVERAGE-01
    rationale: |
      WITHDRAWN CLAIM, not a reduction of legitimate scope. The operator has
      ruled (2026-08-20) that the Histogram tooltip cannot and should not show a
      molecule, and that including the Histogram in this test was incorrect in
      the first place. Scenario 2 once hovered three viewers; the Histogram step
      is gone from both the spec and this body and must NOT be re-added by a
      future producer — it asserted something the product was never meant to do.

      The probe evidence is kept because it is what makes the ruling checkable
      rather than an assertion of taste. From a standalone node+chromium probe
      against dev on 2026-08-20 (not a spec run):
      - smiles.csv's first numeric column, molregno, puts all 972 rows into a
        single bar at the left edge of the Histogram. The viewer's ink is one
        narrow column at roughly 6% of the width, which is why the earlier
        bounding-box-fraction sweep never touched it and reported "no tooltip"
        deterministically — the stale Gate B diagnostic above describes that
        sweep, not a missing tooltip.
      - With the pointer actually on the bar the Histogram DOES raise a tooltip,
        but it is a bin aggregate ("molregno: 1480014..1481727.80 total: 972
        selected: 0 filtered: 972 972 rows") and it carries no canvas at all.
      - Re-probed with six numeric value columns (molregno,
        NumSaturatedHeterocycles, NumAliphaticHeterocycles,
        NumAliphaticCarbocycles, NOCount, NumSaturatedRings), each with several
        populated bars. Every hit had the same shape: bin aggregate, no canvas.

      So on the evidence the Histogram raises only bin/group tooltips and never a
      molecule tooltip. That is the product's design, per the ruling — not a
      defect, and not a coverage gap to close.

      Box Plot and Scatter Plot ARE covered and were verified to raise tooltips
      carrying real painted molecule canvases.

      Atlas discrepancy, recorded but NOT acted on here (reference artefacts are
      out of this scenario's scope; the operator has been asked separately). The
      atlas entry chem.cp.molecule-rendering-end-to-end reads: "Cross-viewer
      tooltip rendering (Box Plot, Scatter Plot, Grid, Histogram) must not crash
      the renderer on null DataFrame." That line asks the renderer MUST NOT
      CRASH; it does not ask for a molecule tooltip. Scenario 2 read it as the
      latter, and that reading was the error — so the atlas may well be correct
      as written while the scenario's interpretation of it was wrong, and the
      two are not necessarily in conflict. The same line also names Grid, which
      this scenario never covered as a tooltip hover either.
    verdict_status: PASS
gate_verdicts:
  a:
    verdict: SCOPE_REDUCTION
    cycle_id: gate-a-molecule-rendering-2026-08-22
    timestamp: 2026-08-22T00:00:00Z
    failure_keys: []
    review_round: 1
    scope_reduction_proposal: |
      Three prose-versus-spec mismatches, all of the same family (the body
      claims an actuation channel or an attribution the paired spec does not
      realize). Each is repaired by a stated reduction of claimed scope, not by
      re-authoring the test. If none is applied, the scenario should return as
      FAIL with failure_keys [A-CONT-01].

      (1) A-CONT-01, Scenario 2 Step 1. The step reads "add a Box Plot and a
      Scatter Plot viewer" and names no control anywhere. The spec adds both
      programmatically (molecule-rendering-spec.ts:147-151, tv.addViewer('Box
      plot') / tv.addViewer('Scatter plot')); no Add-viewer menu gesture is
      exercised. Reduction: state the channel the way Setup already does — e.g.
      "the viewers are added programmatically through the table-view API, so the
      Add > Viewers menu gesture is not covered here; the tooltip hover, which
      is the subject of this scenario, is a real pointer gesture
      (page.mouse.move)". The dependency is real and technical: the hover sweep
      needs deterministic viewer placement, and the spec's ink-targeted sweep
      reads the viewer canvas directly.

      (2) A-CONT-01, Setup step 1. The numbered step still prescribes "use the
      Browse sidebar > Demo Files > chem folder" while the disclosure two lines
      below correctly says the table is opened programmatically. The disclosure
      is honest and matches the spec (openChemTable uses
      grok.dapi.files.readCsv + grok.shell.addTableView, spec:11-38), but the
      step text contradicts it, and a future producer reading only the numbered
      step would author a Browse gesture. Reduction: delete the Browse
      parenthetical from the step, leaving the existing disclosure as the single
      statement of channel.

      (3) Expected-results correspondence, Scenario 3 Steps 3 and 6. Both
      expectations attribute paint to the typed cells ("the reaction cells paint
      rather than leaving the canvas blank", "the cells paint rather than
      leaving the canvas blank"), but the spec measures gridNonBlankPixels,
      which reads the largest grid canvas whole (spec:46-64, called at :190 and
      :214). A blank reaction or mixture column beside painted neighbouring
      columns satisfies it. Scenario 1 is not affected — it scopes its read to
      the column via moleculeColumnColouredPixels(page, 'canonical_smiles')
      (spec:74-98). Reduction: reword both expectations to what is measured —
      "the grid canvas carrying the column is painted" — matching the wording
      the file's own Automation notes already use ("that the canvas is
      painted"). The dependency is real: the column-scoped read that exists is
      colour-based, and whether ChemicalReaction / ChemicalMixture cells emit
      channel-separated ink is unverified, so demanding a column-scoped paint
      assertion here would introduce a new unprobed claim rather than tighten an
      existing one.

      Optional, same direction, not required for PASS: two stated expectations
      are asserted by the spec but carry no expected_results anchor — Scenario 1
      Step 1 (semType Molecule + Molecule renderer bound, asserted at
      spec:124-125) and Scenario 3 Step 2 (at least 5 genuine >> reaction rows,
      asserted at spec:188-189). Both anchor cleanly onto exactly one softStep
      ("Scenario 1 Step 1-2", "Scenario 3 Step 3"). This is the safe direction
      (asserted but unanchored, not claimed but unasserted), so it is a
      completeness note rather than a defect.
    claims:
      - check_id: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses as YAML. feature=chem, priority=p0,
          target_layer=playwright, coverage_type=smoke all present.
          realizes is the empty list, which is correct rather than missing:
          post-inversion realizes carries KG slugs, and the atlas gives
          validates_kg as the empty list for both ids this scenario names
          (chem.yaml:205 for the critical path, :481 for the interaction).
          The coverage channel is realizes_atlas, which carries
          chem.cp.molecule-rendering-end-to-end and chem.int.render-feeds-search
          and is the field compute-coverage.ts:236 reads; a p0 typed by
          realizes_atlas is not the untyped p2 case the checklist describes.
      - check_id: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body carries the level-two headings Setup, Scenarios and
          Automation notes, plus three level-three scenario headings.
      - check_id: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Scenario 1 has steps 1-3, Scenario 2 has steps 1-4, Scenario 3 has
          steps 1-6, and Setup has one numbered step.
      - check_id: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          No scenario is empty; the smallest, Scenario 1, carries three steps.
      - check_id: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer is playwright, which is in the permitted set, and a
          paired Playwright spec exists at molecule-rendering-spec.ts.
      - check_id: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type is smoke, from the test-kind enum.
      - check_id: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type smoke is declared once at frontmatter level and covers
          all three scenarios. It is a test-kind value, not a severity-axis
          value; the severity axis is carried separately by priority p0.
      - check_id: A-STRUCT-04
        status: PASS
        evidence: |
          The shared smiles.csv open is factored into Setup and Scenario 2
          explicitly continues from Scenario 1 rather than repeating it.
          Scenario 3 opens test-reactions.csv and test_mixtures.csv inline,
          which is correct — those are its own fixtures, not repeated setup.
      - check_id: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          pyramid_layer is ui-smoke and coverage_type is smoke, so the hard
          alignment rule is satisfied rather than vacuous.
      - check_id: A-CONT-01
        status: FAIL
        evidence: |
          Names are real throughout — canonical_smiles, reaction, mixture,
          System:DemoFiles/chem/smiles.csv, System:AppData/Chem/test-reactions.csv,
          System:AppData/Chem/test_mixtures.csv, rdkit-cell-renderer.ts — with no
          bracket placeholders or generic stand-ins, and all resolve against the
          spec. The actuation clause fails on two steps. Scenario 2 Step 1 says
          only "add a Box Plot and a Scatter Plot viewer", naming no control,
          while the spec adds them programmatically via tv.addViewer
          (spec:147-151). Setup step 1 names the opposite defect — it prescribes
          "use the Browse sidebar > Demo Files > chem folder" although the spec
          opens the file through grok.dapi.files.readCsv (spec:16-18), and the
          disclosure below the step correctly says so. Both are repaired by the
          stated reduction; see scope_reduction_proposal items 1 and 2.
      - check_id: A-BUG-01
        status: PASS
        evidence: |
          The atlas known_issues block uses the exists-false schema rather than
          the literal needed marker, and lists twelve chem bugs. The only entry
          in scope for this scenario, GROK-16870 ("Fix rdkit cell renderer error
          in scatter plot tooltip"), is addressed by both routes: it is carried
          in related_bugs frontmatter with status fixed, and Scenario 2 names it
          in its heading and asserts the exact crash signature, with the spec
          matching 'gS' and NullError inside its filter (spec:105). The other
          eleven belong to unrelated chem flows; a single-scenario gate cannot
          require this file to carry them, and whether the chain covers them is
          a Gate F question.
      - check_id: A-MERIT-01
        status: PASS
        evidence: |
          Two exclusions are recorded and neither is effort-based. The
          OpenChemLib renderer switch is dropped because the Renderer property
          was removed from the product (Automation notes, operator 2026-08-20),
          and the Histogram hover is withdrawn on an operator ruling backed by a
          probe across six numeric value columns, every hit a bin aggregate with
          no canvas (SR-01). Both cite a product-level dependency. The two
          not-covered notes inside the scenarios — the column-header icon and
          the reaction arrow layout — each state the technical reason (nothing
          in the step reads that channel) rather than pleading difficulty.
      - check_id: A-MERIT-02
        status: PASS
        evidence: |
          No TODO, no deferral, and no next-phase language anywhere in the body
          or frontmatter. Both removals are stated as permanent and as
          must-not-be-re-added, which is the opposite of a deferral. The one
          open item, the atlas discrepancy in SR-01, is explicitly recorded as
          out of this scenario's scope and routed to the operator rather than
          parked as future work on this file.
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-molecule-rendering
    timestamp: 2026-08-22T12:07:00Z
    spec_runs:
      - spec: molecule-rendering-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 26
        failure_keys: []
        run_mode: headless-cold
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-molecule-rendering
    timestamp: 2026-08-22T12:33:13Z
    failure_keys: []
---

# Chem — Molecule rendering end-to-end (RDKit / cross-viewer tooltip)

This scenario covers the backbone rendering path that every other Chem flow
depends on: SMILES/Molecule column cell rendering via the RDKit renderer
property, cross-viewer tooltip rendering of molecule cells in non-Chem viewers
(Box Plot, Scatter Plot), and the sibling ChemicalReaction and
ChemicalMixture cell renderers drawing their own cell types.

## Setup

1. Open `System:DemoFiles/chem/smiles.csv` (contains a SMILES column with
   known valid structures).

Actuation channel: the table is opened programmatically, so the file-open
gesture through Browse is not covered here. Rendering — the subject of this
scenario — is unaffected: the column is typed by the platform's own semantic
detection after the table is loaded, and the cell renderer binds from that.

## Scenarios

### Scenario 1: RDKit renderer on the Molecule column

Steps:
1. With `smiles.csv` open, confirm that the `canonical_smiles` column carries
   semType `Molecule` and that the `Molecule` cell renderer is bound to it.
   (The check reads the column model; the column-header icon itself is not
   inspected.)
2. Wait for the grid to finish painting its first screen of rows.
3. Verify that the Molecule column is drawn as structures rather than as text:
   the column region of the grid canvas is painted, and it carries coloured ink
   — RDKit draws heteroatoms in colour (O red, N blue) while a raw-SMILES text
   fallback is monochrome. Read the same colour measure over the monochrome
   `molregno` column of the same grid as a control: it must carry no coloured
   ink at all, which is what makes the reading over the Molecule column a
   statement about structure ink rather than about grid chrome. Open the browser
   developer console and confirm zero errors are reported by
   rdkit-cell-renderer.ts during rendering.
Expected:
- Scenario 1 Step 3: the Molecule column region is painted and carries coloured
  heteroatom ink, while the `molregno` control column carries none; console is
  clean (no rdkit-cell-renderer errors).

### Scenario 2: Cross-viewer tooltip rendering (GROK-16870 regression lock)

This scenario verifies that hovering over a data point in a non-Chem viewer
triggers the molecule-cell renderer in tooltip context without the null-DataFrame
crash (fixed in 1.22.0, locked here as regression guard).

Steps:
1. With `smiles.csv` open (continue from Scenario 1 with the RDKit renderer
   active), add a Box Plot and a Scatter Plot viewer. Each is added
   with its default axes — no axis is chosen by hand, so the scenario does not
   depend on any particular numeric column being present.

   Actuation channel: the viewers are added programmatically through the
   table-view API, so the Add > Viewers menu gesture is not covered here. The
   hover, which is the subject of this scenario, is a real pointer gesture —
   the sweep needs deterministic viewer placement and reads the viewer's own
   canvas to choose where to point.
2. Open the browser developer console, then hover the mouse over the Box Plot,
   sweeping across it until a tooltip appears.
3. Verify that a tooltip appears and that its molecule region is painted (the
   tooltip carries a canvas that is not blank). Confirm the Box Plot survives
   the hover, and that the console shows no errors from
   `rdkit-cell-renderer.ts` — specifically, no
   `NullError: method not found: 'gS' on null` message.
4. Sweep the mouse over the Scatter Plot the same way and verify a tooltip
   appears with a painted molecule region, the viewer survives, and no console
   errors appear.

Expected:
- Scenario 2 Step 3: a tooltip appears over the Box Plot with a painted
  molecule region; the viewer survives; no console errors from
  rdkit-cell-renderer.ts.
- Scenario 2 Step 4: same over the Scatter Plot.

### Scenario 3: ChemicalReaction and ChemicalMixture cell renderers

Steps:
1. Open `System:AppData/Chem/test-reactions.csv`, whose `reaction` column holds
   reaction SMILES and is detected as semType `ChemicalReaction`.
2. Wait for the grid to paint; the fixture must carry at least 5 genuine
   reaction rows (values containing the `>>` arrow) for the step to mean
   anything.
3. Verify that the `reaction` column carries semType `ChemicalReaction`, that
   the `ChemicalReaction` cell renderer is bound to it, and that the region of
   the grid canvas occupied by the `reaction` column is painted rather than
   blank. Check the browser console for no renderer errors.

   Not covered: where the arrow sits within a cell. Nothing in this step reads
   the layout, so "reactants on the left, products on the right" is not
   verified — only that the reaction renderer is bound and that the cells paint.
4. Open `System:AppData/Chem/test_mixtures.csv`, whose `mixture` column is
   detected as semType `ChemicalMixture`.
5. Wait for the grid to paint.
6. Verify that the `mixture` column carries semType `ChemicalMixture`, that the
   `ChemicalMixture` cell renderer is bound to it, that every cell parses into
   at least one component and at least one cell carries more than one
   component, and that the region of the grid canvas occupied by the `mixture`
   column is painted. Check the browser console for no renderer errors.

Expected:
- Scenario 3 Step 3: the reaction column is typed `ChemicalReaction` with its
  renderer bound and its own column region painted; console is clean.
- Scenario 3 Step 6: the mixture column is typed `ChemicalMixture` with its
  renderer bound, component counts parse (min ≥ 1, max > 1), and its own column
  region painted; console is clean.

## Automation notes

- The OpenChemLib renderer and the `Renderer` column property were REMOVED from
  the product (operator, 2026-08-20; no `Renderer` choice survives in
  `packages/Chem/src`). Scenario 1's former renderer-switch steps are deleted
  outright — this is a removed feature, not a deferred check, so it is NOT a
  scope reduction and must not be re-added by a future producer.
- Structure-versus-text discrimination is done by colour, not by pixel count. A
  raw-SMILES text fallback paints plenty of dark pixels, so a non-blank count
  says only "something was drawn". RDKit's palette draws heteroatoms in
  saturated colour (O ≈ rgb(255,0,0), N ≈ rgb(0,0,255)) while canvas text is
  monochrome, so a pixel whose channels differ from each other is ink only a
  structure renderer produces. The non-blank count is kept alongside it as a
  cheap fault guard against an unpainted or unreadable canvas — it does not
  distinguish structures from text and its message must not claim to.
- The colour measure reads the column strip from the top of the canvas, so it
  also covers the column-header band and the filter chrome the setup pins there.
  The same measure is therefore taken over `molregno`, a plain int column of the
  same grid, in the same run, and required to read zero. Without that control a
  single coloured chrome pixel would satisfy the Molecule assertion while the
  product rendered raw SMILES text — exactly the regression the step names. A
  standalone node+chromium probe against dev (2026-08-22) measured 856 coloured
  pixels for `canonical_smiles` and 0 for `molregno` and for every other numeric
  column on screen, so the chrome in the band is greyscale; the control asserts
  that in-run instead of relying on it.
- For Scenario 2: sweep the mouse over the viewer until `.d4-tooltip` appears,
  then read the largest canvas inside the tooltip and require it to be
  non-blank. Hover positions are taken from where the viewer's own canvas
  carries ink, densest area first, rather than from fractions of the viewer's
  bounding box: the box includes the header and the axis gutters, and a viewer
  whose ink sits in one narrow band is missed by any central fraction.
  "Non-blank" counts a pixel only when it is opaque as well as dark: an
  untouched canvas reads (0,0,0,0), so a colour-only rule scores a canvas that
  was never drawn on at its full area. Console errors during the hover are
  captured through the wrapped `console.error` channel installed when the table
  is opened; the assertion is that no message matching
  `rdkit-cell-renderer` appears. The tooltip check is
  the primary rung — without it the console check passes when no tooltip ever
  appeared.
- Scenario 2 covers Box Plot and Scatter Plot only. The Histogram tooltip cannot
  and should not show a molecule, and including the Histogram here was incorrect
  in the first place (operator, 2026-08-20; probed across six numeric value
  columns, every hit a bin aggregate carrying no canvas). Its hover steps are
  deleted outright — this is a withdrawn claim, not a deferred check, so it is
  NOT a scope reduction and must not be re-added by a future producer. SR-01
  keeps the ruling and the probe evidence behind it.
- For Scenario 3: the reaction and mixture steps assert semType, cell-renderer
  binding, component parsing, and that the typed column's own region of the
  canvas is painted — scoped to the grid column's x-range the way Scenario 1
  scopes its colour read, so a blank reaction or mixture column beside painted
  neighbours no longer passes. The floor stays a fault floor and is placed from
  measurement (probe against dev, 2026-08-22): the painted reaction strip reads
  6324 and the mixture strip 3963, against 400 for an equally sized blank strip
  of canvas below the last data row and at most 937 for a header band. Neither
  step reads glyph layout, so neither can carry an arrow-position or a raw-text
  claim.
