---
feature: chem
realizes_atlas:
  - chem.cp.substructure-search-with-filter
  - filters.cp.chem-and-bio-filters
realizes:
  - viewers.filters
priority: p0
target_layer: playwright
pyramid_layer: integration
coverage_type: regression
related_bugs:
  - id: GROK-18530
    status: fixed
  - id: GROK-18383
    status: fixed
  - id: github-3445
    status: fixed
  - id: GROK-20001
    status: fixed
  - id: GROK-14952
    status: fixed
realized_as:
  - chem-filters-spec.ts
scope_reductions:
  - id: SR-01
    check: E-TRACE-02
    rationale: |
      Step 9 asks for the drawn scaffold to be ROTATED inside the sketcher
      canvas. That is an in-canvas gesture on a third-party sketcher widget with
      no DOM handle and no Playwright primitive: there is no element to grab, no
      rotation API on the filter's saved state, and a raw mouse drag across the
      canvas draws a bond rather than rotating the fragment. The rotation is
      therefore NOT performed and nothing in the spec claims it.
      Two things stand in for it, and they are measured as separate claims.
      (1) The realignment path itself — the Align + Highlight checkboxes on the
      sketcher, which re-tag the column and repaint the rendered structure column
      (chem-substructure-filter.ts:373-382, :423-437) — measured as a real pixel
      delta on the grid canvas, under a sustained-zero negative control, and
      isolated from row-set change by asserting df.filter.trueCount is unmoved.
      (2) The in-sketcher edit that the rotation would have been an instance of —
      a further, DIFFERENT structure entered into an already-applied filter
      through the sketcher's SMILES field, which runs the same
      setValue -> setMolecule -> sketcher.onChanged path any drawing gesture runs
      (js-api/src/chem.ts:386-392, :653-656) and is therefore gated by
      "Filter as you draw" exactly as an in-canvas edit is. What that edit is
      NOT is a geometric transform of an unchanged molecule, so the alignment of
      the SAME structure under rotation remains unmeasured. The Step 9
      expectation text cites this entry.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Step 4's last bullet asks that switching through Exact, Similar, Included
      in and Not included in each CHANGE the displayed row count ("no mode
      silently stalls the filter"). That claim is dataset-dependent and
      unrealizable on spgi-100 with this probe: Exact, Similar and Included in
      all legitimately yield 0 rows, so three of the six modes produce the same
      count as one another and a "the count moved" assertion would fail on
      correct product behaviour — or, worse, would have to be weakened into
      something a stalled filter also satisfies.
      Asserted in its place are the dataset-INDEPENDENT invariants that a
      stalled mode cannot satisfy: Contains + Not contains == full row count,
      Included in + Not included in == full row count, Exact is a subset of both
      Contains and Included in, and Similar (at the default cutoff) contains
      every Exact match. Each of the six modes is still driven through the
      card's own search-type control and each is confirmed to be OFFERED by that
      control, so a build that drops or renames a mode fails by name. The prose
      bullet is kept because it states the user-facing intent; this entry records
      that its literal form is not what the spec measures.
    verdict_status: SCOPE_REDUCTION
  - id: SR-03
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Step 5's second bullet asks that selecting a scaffold drop the row count
      "below the Contains row count recorded in Step 4". With both cards armed
      the two filters AND together, so the combined count can legitimately EQUAL
      the Contains count whenever every row the substructure criterion keeps also
      carries the seeded scaffold — a strict inequality would then fail on
      correct behaviour, and which of the two holds is a property of the dataset,
      not of the product.
      Asserted in its place: with the Structure substructure card switched OFF
      (and the switch-off itself confirmed, so the measurement really is
      isolated) the scaffold card alone holds the row count strictly below the
      full row count and strictly above zero — which is what GROK-18383 is
      about, a card that is added but filters nothing — and, with the
      substructure card re-armed, the AND of both cards passes no more rows than
      the substructure card alone.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: "Scenario 1, Step 1:"
    expectation: >-
      spgi-100 opens with one Chem substructure card per molecular column
      already in the panel. Each of those cards goes down through its own X
      icon, and once they are gone the panel holds no substructure card and none
      comes back on its own — on the saved-state channel and in the DOM. Driving
      the Filter Panel hamburger menu Add Filter > Substructure Filter... then
      opens the "Select columns..." picker, which offers exactly the molecular
      columns; confirming it with every one checked leaves each of those columns
      carrying exactly one Chem substructure card, the Structure card among
      them, and grows the panel's card count by exactly the number of columns
      the picker reported checked.
  - anchor: "Scenario 1, Step 2"
    expectation: >-
      The card's own sketcher opens from the card body and exposes its own
      SMILES field; the probe typed into that field and committed arms the
      Structure card's filter.
  - anchor: "Scenario 1, Step 3"
    expectation: >-
      Fewer than 100 rows pass the filter once the substructure is committed,
      and at least one row still passes. The Filter Panel header's active-filter
      counter is present, carries text and reads exactly 1 — it reports how many
      cards are FILTERING, not how many rows pass, and the Structure card that
      just took the probe is the only filtering card. The Structure card's saved
      state carries the committed structure: a molblock whose counts line
      declares atoms, or the non-empty SMILES string when the probe was entered
      through the sketcher's SMILES field.
  - anchor: "Scenario 1, Step 4"
    expectation: >-
      The Structure card's own search-type control offers all six named search
      types, and every count read back lies between 0 and 100. Contains passes
      exactly the row count Step 2 measured for the same probe in the same mode,
      itself strictly between 0 and 100. Contains plus Not contains equals 100,
      and so does Included in plus Not included in — each complement pair
      partitions the dataset exactly. Exact is a subset of both Contains and
      Included in, and Similar at the default cutoff passes every Exact match.
      Returning the card to Contains reproduces the Contains count. See
      scope_reductions SR-02 for the per-mode "the count moved" clause, which is
      dataset-dependent and is not asserted.
  - anchor: "Scenario 1, Step 5"
    expectation: >-
      GROK-18383 — the Scaffold Tree Filter picker offers exactly the molecular
      columns, and confirming it leaves exactly one new scaffold-tree card per
      checked column relative to the card count Step 1 recorded, with no column
      carrying two. The scaffold reaches the card through the card's own
      controls: the "Sketch scaffolds manually" plus icon opens the "Add new
      scaffold..." sketcher, Add places the scaffold, and clicking the node it
      creates leaves that node CHECKED — an unchecked node is not a criterion.
      With the Structure substructure card switched off — the card reading
      disabled and its checkbox unchecked — the scaffold card alone holds the
      row count strictly below 100 and strictly above zero. Selecting the
      scaffold updates the card the menu added rather than creating a second
      one, and with the substructure card re-armed the two cards together pass
      no more rows than the substructure card alone. See scope_reductions SR-03
      for the strict "below the Contains count" form of the prose, which is
      dataset-dependent and is not asserted.
  - anchor: "Scenario 1, Step 6"
    expectation: >-
      GROK-20001 — with the panel's substructure cards taken down and the
      Structure column confirmed to hold no substructure state, right-clicking a
      molecular cell makes exactly that cell current: the grid's first visible
      row of the Structure column, pinned to the table row the product's own hit
      test reports for the aimed cell rather than to a fixed index, because the
      scaffold-tree criterion armed in Step 5 is still filtering here and the
      first visible row is not table row 0. Current Value > Use as filter then
      adds exactly ONE card, not two. The Structure column then carries exactly
      one substructure state, the new card carries the clicked cell's structure
      rather than an empty sketch, and the row count is below 100 and above zero
      — the clicked row itself still passes, the cell's own molecule matching at
      least itself.
  - anchor: "Scenario 1, Step 7"
    expectation: >-
      The panel rebuilt through Remove All, close and reopen comes back with
      exactly one substructure state and one substructure card body per
      molecular column, exactly one of them on Structure, and that card comes up
      with no structure, so this step's own probe is what arms it. GROK-18530 —
      switched through the card's own controls to Similar mode at a non-default
      0.6 cutoff, the original card reads Similar, carries the 0.6 cutoff and
      holds a molecule. The top menu View > Layout > Clone View then produces a
      second view holding its own filter group and its own substructure state
      object, and the clone's substructure card carries the same search type,
      the same similarity cutoff, the same fingerprint type and the same
      molecule as the original, with both cards displaying Similar. Row counts
      are not compared across the two views — see the Automation notes.
  - anchor: "Scenario 1, Step 8"
    expectation: >-
      Immediately before the reset the panel is filtering, a card holds a drawn
      structure, a card body is painting it, and a card is in a non-default
      search mode. After the Filter Panel header's reset icon is clicked, all
      100 rows pass again, the header's active-filter counter falls to 0, every
      card is still present — the card count and the number of saved
      substructure states are both unchanged, so the reset clears criteria
      rather than removing cards — no card's saved molblock declares any atoms
      any more, and every card is back in the default Contains mode. No card
      body still PAINTS a structure: read as non-background pixels on each
      body's own canvas, the painting is gone within about a second of the
      reset, and it is STILL gone when every body is re-read 10 s after the
      clearing run ends. Every body still offers a sketch area — a canvas or the
      "Sketch" placeholder — for that reading to be taken from, and every canvas
      that is still there is READABLE: a body whose canvas cannot be measured at
      all fails, rather than counting as one that paints nothing. The watch has
      a 60 s ceiling, which only a reset that never clears ever spends. Which of
      those two an emptied body ends up showing is recorded per card and
      deliberately not asserted: the armed card comes back blank rather than
      placeholdered, which is the open bug GROK-20739. The clone view opened in
      Step 7 is closed first and a single table view is left, so the reset acts
      on the original panel.
  - anchor: "Scenario 1, Step 9"
    expectation: >-
      github-3445 — two separate claims hold, and "Filter as you draw" gates one
      of them and not the other. CLAIM 1 — with a structure already applied and
      Align and Highlight OFF, turning them ON repaints the rendered structure
      column immediately, with the sketcher dialog still open: a positive pixel
      delta on the grid's own data canvas over a proven sustained-zero baseline,
      while df.filter.trueCount does not move. It holds with "Filter as you
      draw" cleared and with it checked. CLAIM 2 — with "Filter as you draw"
      CLEARED, a further, DIFFERENT structure entered into the already-filtering
      sketcher does not reach the grid while the dialog stays open: no repaint
      of the structure column and no movement of the row count, held across
      re-reads that span the window in which the checked half lands, and
      confirming with OK then applies it, moving the row set to a different,
      still-filtering count. With "Filter as you draw" CHECKED the same further
      edit reaches the grid LIVE with the dialog still open, repainting it as
      well as moving the row count, and lands on exactly the row set the cleared
      half reached only on OK. The cleared half's opening probe behaves the same
      way: entered with the checkbox already clear, it leaves the row set at the
      full 100 — held across the same re-reads on both the pixel and the
      row-count channel — until OK is pressed. In both halves the checkbox is
      set to the state under test before the probe molecule is entered and read
      back; both claims start from Align and Highlight OFF, read back as
      cleared; both halves measure on a filtered grid, strictly between zero and
      100 and on the same row set; and no new console errors appear in either.
      The in-canvas scaffold rotation itself is not performed — see
      scope_reductions SR-01.
  - anchor: "Scenario 1 Step 10"
    expectation: >-
      GROK-14952 — the same molecular column is filtered from two surfaces and
      the two must not diverge. With every substructure card taken down (a held
      zero on the saved-state channel and in the DOM) and the row count
      sustained at the full 100, the Structure grid column header's own
      column-options icon — sitting inside that column's header band — opens a
      popup titled Structure holding exactly ONE substructure filter body,
      OUTSIDE the Filter Panel, while the panel holds none and its filter group
      carries no Structure substructure state. The probe entered into the
      popup's own sketcher field narrows the table to a sustained count strictly
      between 0 and 100. The popup's own "Add filter" link then moves the
      criterion to the panel: the popup filter is gone from the document and the
      panel gains exactly one card carrying that structure, whose own sustained
      count is strictly between 0 and 100 and equals the popup count exactly.
      Dismissing the popup with a neutral click — the popup asserted on screen
      first, at a point computed from its own rectangle, and asserted gone
      afterwards — leaves the row count where it was. Switching the panel card
      off through its own enable checkbox returns the row count to exactly 100
      and HOLDS it there; a count still below 100 means a surviving popup-side
      filter is still contributing its mask. Re-arming the card restores exactly
      the count it held before.
  - anchor: "Scenario 1, Step 11"
    expectation: >-
      On a panel rebuilt from scratch (Remove All, close, reopen) with the
      auto-populated substructure cards taken down, so the Structure column
      carries no substructure state, and with at least one card on a
      non-molecular column already seated, dragging the Structure grid column
      header onto the panel as a real pointer drag — released inside the panel's
      own "Add filter" drop zone — leaves the panel holding exactly one more
      card. The card at index 0 is the Structure Chem substructure card, on the
      column name and on the filter source alike, and it arrives with NO
      structure in its saved state. The row count, read as a sustained value, is
      exactly what it was before the drag.
  - anchor: "Scenario 1, Step 12"
    expectation: >-
      With the substructure filter armed and narrowing strictly between 0 and
      100, a category on a non-molecular string column that keeps some
      substructure hits, drops others and also covers rows the substructure
      filter rejects is applied through the card's own category row. Both
      criteria armed, the visible row set is exactly the intersection of the two
      — a pinned count, not a bound — and strictly below each filter alone.
      Switching the substructure card off through its own checkbox leaves
      exactly the category's own row count, and switching it back on restores
      exactly the intersection, so neither criterion destroyed the other.
  - anchor: "Scenario 1, Step 13"
    expectation: >-
      With the Structure card filtering and carrying a drawn structure,
      switching it off through its own enable checkbox and then closing and
      reopening the panel brings the card back DISABLED — its own enable
      checkbox reads unchecked and the card carries its disabled class. The
      restored card still carries the structure drawn before the collapse, and
      the row count is still the switched-off one, held as a sustained reading.
      Switching the card back on reproduces exactly the row count it filtered to
      before the panel was collapsed.
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: direct-gate-a-2026-08-22
    timestamp: 2026-08-22T00:00:00Z
    review_round: 1
    failure_keys: []
    claims:
      - check_id: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses as YAML. feature=chem, priority=p0,
          target_layer=playwright, coverage_type=regression are all present.
          realizes carries one non-empty id (viewers.filters) and
          realizes_atlas carries the two mandatory atlas type-ids
          (chem.cp.substructure-search-with-filter,
          filters.cp.chem-and-bio-filters). The two-axis form matches the
          section convention verified against chem-descriptors-docker.md,
          chem-empty-input-analyses.md and chem-grok-12758.md, where
          realizes holds KG ids and realizes_atlas holds atlas ids. No
          produced_from is declared, so the non-empty realizes requirement
          is not even binding here; both lists are non-empty regardless.
      - check_id: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          The body carries the level-2 headings Setup, Scenarios and
          Automation notes, with the scenario itself under
          "### Scenario 1: Chem substructure filter ...".
      - check_id: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Scenario 1 carries thirteen numbered steps, 1. through 13., and
          the Setup section carries three of its own.
      - check_id: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          The single scenario is the largest in the section; no empty
          scenario heading exists.
      - check_id: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer reads playwright, and the paired
          chem-filters-spec.ts is a Playwright spec, so the declaration and
          the realization agree.
      - check_id: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type reads regression, which is in the test-kind enum.
          The scenario is bug-guard heavy (five fixed bugs named in the
          body), so regression is the right kind rather than smoke.
      - check_id: A-STRUCT-03
        status: PASS
        evidence: |
          The coverage_type label sits at file-frontmatter level and is a
          test-kind value, not a severity-axis one. priority=p0 is carried
          separately on its own field, so the two axes are not conflated.
      - check_id: A-STRUCT-04
        status: PASS
        evidence: |
          Table opening, the Chem-availability gate and the Filter Panel
          opening are factored into the "## Setup" section and are not
          restated in Scenario 1. The two in-scenario repeated
          preconditions (the Remove All / close / reopen rebuild used by
          Steps 7, 9 and 11) are stated once in the Automation notes and
          referenced from the steps, and the spec factors them into the
          shared rebuildFullSubstructurePanel and armStructureCard
          helpers.
      - check_id: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          PASS by vacuity. The frontmatter declares no pyramid_layer, so
          no ui-smoke alignment rule applies. Sibling Chem scenarios do
          carry pyramid_layer=ui-smoke; its absence here is consistent
          with a regression scenario and is not a Gate A defect.
      - check_id: A-CONT-01
        status: PASS
        evidence: |
          Every name in the body is real and checkable against the spec.
          The dataset is System AppData/Chem/tests/spgi-100.csv with its
          row count stated as 100 (CHEM_PATH, CHEM_FULL); the molecular
          column is Structure (CHEM_COL) and the body says explicitly that
          this dataset carries no SMILES column. Every actuation step
          names the concrete control and where it lives, not the effect —
          the Filter Panel hamburger Add Filter > Substructure Filter...
          and Add Filter > Scaffold Tree Filter..., the Select columns...
          picker with its All link and OK, each card's own X icon, the
          card's search-type control in the gear popup, the scaffold
          card's "Sketch scaffolds manually" plus icon with its "Add new
          scaffold..." sketcher and Add button, the grid cell's right-click
          Current Value > Use as filter, the top menu View > Layout >
          Clone View, the panel header's reset icon, the sketcher's Align,
          Highlight and "Filter as you draw" checkboxes, the card's own
          enable checkbox, the Structure grid column header's
          column-options icon and the popup's "Add filter" link, and a real
          pointer drag of the column header into the panel's "Add filter"
          drop zone. The one derived value is Step 12's second criterion,
          which is deliberately chosen from the data; the step states the
          selection rule in full (a category that keeps some substructure
          hits, drops others, and covers rows the substructure filter
          rejects) and pickIntersectingCategory implements exactly that
          rule and fails the step when no such category exists, so it is a
          specified derivation rather than a placeholder stand-in.
      - check_id: A-BUG-01
        status: PASS
        evidence: |
          The chem atlas known_issues list carries 12 entries, all with
          test_coverage.exists false. The one that falls inside this
          scenario's critical paths, GROK-14028 (Reset filter does not
          clear the sketcher input line), is addressed under clause (b) —
          Step 8 names the GROK-14028 family in the body and the spec
          asserts the clearing on the painted-pixel channel plus the saved
          state. The remaining 11 belong to other Chem critical paths
          (MMP, R-Group, Chemical Space, Descriptors, MPO, project save,
          scatterplot tooltip, scaffold-tree rendering) and are the
          business of other scenarios in the section; deciding whether the
          section as a whole covers them is a Gate F predicate and is not
          decided here. Separately, GROK-20739 is an OPEN bug held out of
          related_bugs by the known-open-bug mechanism and referenced in
          prose (Step 8 bullet plus two Automation notes), which is clause
          (b) and is also the section's documented convention (the same
          shape as GROK-20731 in add-remove-entry-points.md).
      - check_id: A-MERIT-01
        status: PASS
        evidence: |
          The three scope_reductions each cite a real technical or
          dataset-level dependency, never effort. SR-01 cites the absence
          of any Playwright primitive for an in-canvas sketcher rotation —
          no DOM handle, no rotation API on the saved state, and a raw drag
          draws a bond instead. SR-02 and SR-03 cite dataset dependence on
          spgi-100, where the literal prose form would fail on correct
          product behaviour, and each names the dataset-independent
          invariant asserted in its place. All three are already
          adjudicated with verdict_status SCOPE_REDUCTION on Gate E keys
          (E-TRACE-02, E-EXPECT-COVERAGE-01), not on any Gate A check.
      - check_id: A-MERIT-02
        status: PASS
        evidence: |
          No TODO, no "add later", no "next phase" anywhere in the body,
          the Automation notes or the Spec-must-keep list. The only
          deferrals are the three SR entries, each bound to a stated
          prerequisite rather than to a future intention.
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-filters-r2
    timestamp: 2026-08-22T13:19:30Z
    spec_runs:
      - spec: chem-filters-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 186
        run_mode: headless-cold
        failure_keys: []
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-filters-r2
    timestamp: 2026-08-22T13:57:02Z
    failure_keys: []
---

# Filters — Chem package filters

## Setup

The order below matters: the Filter Panel is the first thing that depends on Chem, so the
Chem-availability decision is taken **before** the panel is opened.

1. Open the file **System:AppData/Chem/tests/spgi-100.csv** (100 rows; the molecular column used
   throughout is `Structure` — there is no `SMILES` column on this dataset), WITHOUT opening the
   Filter Panel yet. Note the full row count (100 rows).
2. Confirm the Chem package is loaded and initialized. Four things are read, all off the table and
   none off the panel: the `Structure` column reports the `Molecule` semantic type; the Substructure
   Filter widget is registered; that filter factory is **applicable** (applying it is both what
   pulls the package in and the signal that it is initialized — see Automation notes); and the
   molecular-cell **Use as filter** action is registered. If any of the four is missing, the
   scenario cannot be observed at all — stop here and report it as **skipped**, naming
   **which of the four** was missing and what value was read for it. Do not treat it as a failure,
   and do not report a build-level conclusion ("Chem is not installed") that the four readings do
   not support.
3. Only then open the Filter Panel by clicking the funnel icon in the ribbon toolbar, and wait for
   at least one filter card to appear (about 2 s on a healthy build). Past the check in step 2, a
   panel that comes up with a complete header and zero cards is a **different** condition — Chem is
   registered and the panel is still empty — and must be reported as such, naming what was seen:
   whether the panel root exists, how many `.d4-filter` cards there are, and which columns report
   the `Molecule` semantic type.

## Scenarios

### Scenario 1: Chem substructure filter — search modes, scaffold tree, Use-as-filter entry point, clone-view similarity, reset and alignment, and panel/column-popup non-divergence

Steps:

1. The panel opens with a substructure filter card already on every molecular column, so remove
   those cards first — click the **X** icon on each substructure card in turn (not **Remove All**,
   which would clear the rest of the panel too) and confirm that no substructure card is left and
   that none comes back on its own. Only then open the viewer hamburger menu (the three-line or
   dots icon in the Filter Panel title bar), hover **Add Filter**, then click
   **Substructure Filter...**. A **Select columns...** dialog opens listing the molecular columns;
   check the columns to filter on (the **All** link checks every one) and press **OK**. Wait for
   the substructure filter card to appear in the panel on each checked column, and verify the card
   count grew by exactly the number of columns checked. Note the current filter card count.
   - Verify the number of cards removed equals the number of molecular columns — the panel really
     did open with one substructure card per molecular column, so their removal is a real gesture.
   - Verify the picker reports exactly the molecular columns checked, and that after OK every one
     of those columns carries exactly ONE substructure filter state — not zero, not two.

2. In the substructure filter card, place a probe molecule in the sketcher input area.
   The sketcher canvas cannot be driven by clicking individual atoms, but the sketcher dialog's
   SMILES field can: open the card's sketcher, type the probe SMILES into that field and commit it
   with Enter (confirming with **OK** if the dialog is still open). Then wait for the number of rows
   passing the filter to drop below 100. This confirms the filter is engaged before
   the search-type switch test begins. Seeding the filter state programmatically (see Automation
   notes) is the fallback if that field is unavailable.

3. Verify the filter is engaged on both channels at once. On the row-set channel: fewer than 100
   rows pass (the substructure filter is active and narrowing the dataset) and at least one row
   still passes — a probe matching nothing would leave the Step 4 partition arithmetic empty. On
   the header channel: the Filter Panel header's **active-filter counter** reads exactly **1**. That
   counter shows how many cards are currently filtering — it is not a row count — so read it as
   such: it must be present and carry non-empty text before its value is compared, and because it is
   debounced off the filter events it is polled until its text settles rather than read in the same
   turn as the sketch.
   - Verify the Structure card's saved state really carries the committed structure. "Carries a
     structure" is read as: a molblock whose counts line declares atoms, or — when the probe was
     typed into the sketcher's SMILES field, which is what the saved state then holds — a non-empty
     string. Never a raw atom count on its own, and never `is not null`: an empty sketch is a
     well-formed molblock with a `0 0 0` counts line and would pass a null check.

4. Switch across all six search-type options in the Substructure filter gear (or filter-menu)
   icon on the filter card (hover to reveal the icon if needed):
   - Verify each of the six options is actually OFFERED by the card's search-type control before
     selecting it, so a build that drops or renames a mode fails by name instead of silently
     leaving the previous mode in place.
   - Select **Contains** and note how many rows pass the filter.
   - Select **Not contains** and note how many rows pass the filter.
   - Verify that the two row counts add up to exactly 100 (the two modes partition the
     dataset exactly for the same molecule, with no row belonging to both or neither).
   - Also switch through **Exact**, **Similar**, **Included in**, and **Not included in**
     to confirm each option changes the displayed row count (no mode silently stalls the
     filter). See scope_reductions **SR-02**: the per-mode "the count moved" form of this clause is
     dataset-dependent — on spgi-100 Exact, Similar and Included in all legitimately return 0 — so
     what is measured instead are the invariants below, which a stalled mode cannot satisfy.
   - Verify **Included in** plus **Not included in** also add up to exactly 100 — the second
     complement pair partitions the dataset just as the first does.
   - Verify **Exact** is a subset of both **Contains** and **Included in**, and that **Similar** at
     the default cutoff includes every **Exact** match. These orderings are dataset-independent.
   - Verify every count read lies between 0 and 100 inclusive — a mode that returned a nonsense
     count would otherwise satisfy the sums by accident.
   - Return the card to **Contains** at the end and verify the row count reproduces the Contains
     count measured above, so the downstream steps run against a known, active filter.

5. With the Substructure filter still active, open the hamburger menu → **Add Filter** →
   **Scaffold Tree Filter...**, then check the columns in the **Select columns...** dialog and
   press **OK**. Wait for the new filter cards to appear in the panel.
   - Verify this picker too reports exactly the molecular columns checked, as the Substructure picker in
     Step 1 does — a picker that offered fewer columns would make the card-count delta below agree with a
     build that only half-added the filter.
   - Verify the panel gained exactly one new filter card per checked column relative to the
     count after Step 1, and that no column ended up with two (regression guard for
     GROK-18383: the scaffold tree filter card was sometimes not added at all).
   - Verify no card appeared between Step 1 and this step — the GROK-18383 delta is measured
     against the count recorded in Step 1, and a card that arrived in between would corrupt it.
   - Select a scaffold in the scaffold tree filter card, using the card's own controls: click its
     **Sketch scaffolds manually** (plus) icon, type the scaffold SMILES into the
     **Add new scaffold...** sketcher that opens, confirm with **Add**, then click the node that
     appears in the tree — the click is what selects it. Verify the node reads **checked**
     afterwards; an unchecked node is not a criterion and the card filters nothing. Verify the
     number of rows passing the filter drops further (below the Contains row count recorded in
     Step 4).
     See scope_reductions **SR-03**: with both cards armed the two criteria AND together and the
     combined count can legitimately equal the Contains count, so what is measured instead is the
     isolated effect described in the next two bullets.
   - Verify selecting the scaffold did not add a second card — the card count is unchanged, so the
     scaffold updated the card the menu created rather than creating another one.
   - Switch the **Structure** substructure card OFF through its own enable checkbox, and verify the
     switch-off really landed (the card reads disabled and its checkbox is unchecked) before
     measuring anything — otherwise the "isolated" reading is not isolated. Verify the scaffold card
     on its own then holds the row count strictly below 100 and strictly above zero: a card that was
     added but filters nothing is the same user-visible defect GROK-18383 describes.
   - Re-arm the substructure card and verify the AND of both cards passes no more rows than the
     substructure card alone.

6. Current Value > Use as filter (regression guard for GROK-20001):
   - Remove the substructure cards currently in the panel through each card's own **X** icon, and
     verify the molecular column holds no substructure state afterwards. The action behind the menu
     item is an update-or-add on the same column and filter type: with a card already present it
     would UPDATE that card and add nothing at all, so the duplicate-card regression would be
     unobservable and the step would pass on a broken build.
   - Note the current number of filter cards in the panel.
   - In the data grid, right-click a cell in the **Structure** column and choose
     **Current Value > Use as filter** from the context menu. Do not make that cell current through
     any other route first: after the right-click, and before choosing the menu item, verify the
     current cell really is the one that was clicked — the grid's first visible row of the
     **Structure** column. Read that cell's table row from the product's own hit test and compare
     the current row against it; a fixed index such as 0 would be wrong, because the scaffold-tree
     criterion armed in Step 5 is still filtering here, so the first visible row is some later table
     row. Otherwise a click that landed on the wrong cell still produces the expected structure and
     the check below tests nothing.
   - Wait for the Filter Panel to update.
   - Verify the panel now contains exactly one more filter card than before (GROK-20001:
     the action was creating two cards instead of one), and that the molecular column carries
     exactly ONE substructure state — the same claim read on the state channel as well as the DOM
     one.
   - Verify the new card carries the clicked cell's structure, not an empty sketch.
   - Verify the number of rows passing the filter is below 100 and above zero — the cell's
     own molecule must at least match itself.

7. First put the panel back into the configuration this scenario describes, so what Step 8 later
   resets is guaranteed rather than inherited from the steps above: use **Remove All** in the panel
   menu, close the Filter Panel, and reopen it from the ribbon funnel icon — every eligible column
   comes back represented by exactly one card. Verify that every molecular column carries exactly
   one substructure filter, that there is one substructure card per molecular column and exactly one
   of them on **Structure**, and that it comes up empty; then draw this step's own probe into it.
   In the original view, switch the substructure filter to **Similar** mode using the filter
   card's gear menu, and set a similarity cutoff that is NOT the default 0.8 (use 0.6, so a
   cutoff that gets dropped is visible rather than coincidentally right). Confirm the card
   really shows Similar at that cutoff. Only then clone the view through the product's own command,
   the top menu **View > Layout > Clone View**, and open the Filter Panel in the clone. That command
   saves the view's layout and rebuilds a view from it, and the filter-state save-and-restore it
   runs is the code the bug lived in — simply opening a second view over the same table skips it,
   and every comparison below would then be reading the original's own values back.
   - Before comparing any field, verify the clone really is a second view with its own filter group
     holding its own substructure state object, not the original's read twice.
   - Before cloning, verify the original card's saved state really reads Similar, really carries
     the non-default cutoff, and really holds a molecule. Without those three the clone comparison
     below would be comparing two unset values.
   - Verify the clone's substructure filter card came up in **Similar** mode carrying the
     same cutoff, the same fingerprint type and the same molecule as the original, and that
     both cards display Similar (regression guard for GROK-18530: the Similar cutoff was
     dropped on clone). Each compared value is guarded for presence first — the fingerprint is
     asserted to BE a string and the molecule asserted to be a real structure before either is
     compared — so "equal" can never mean "both missing".
   - Do NOT compare the two views' row counts. A clone is a second view over the same table
     and the same filter bitset, so the two counts are the same number by construction and
     that comparison cannot fail even while the bug is present.

8. Close the clone view opened in Step 7 and confirm a single table view is left, so the reset acts
   on the original panel. Then click the reset icon (the curved arrow icon) in
   the Filter Panel header bar. Wait for the panel to update.
   - Immediately BEFORE the reset, verify all three positive halves — otherwise every "nothing is
     left" clause below is satisfied by a panel that never held anything: a substructure card is
     present, at least one card holds a drawn structure, at least one card is in a non-default
     search mode, and the panel is actually filtering (row count below 100).
   - Verify all 100 rows pass the filter again.
   - Verify the Filter Panel header's **active-filter counter** falls to **0** — no card is filtering
     any more. Step 3 only ever read that counter at 1, which a stale or hardcoded indicator also
     satisfies; the reset is the contrast that tells them apart.
   - Verify the cards are all still present — the card count in the panel after the reset equals the
     count before it, and so does the number of saved substructure states. Reset clears criteria; it does
     not remove cards.
   - Verify the substructure filter card's sketcher input area is visually clear — the
     card body shows no drawn structure remaining after the reset. Read this on what the card body
     actually DRAWS, not only on the saved state and not on whether a drawing area is still there:
     before the reset at least one card body must be showing a drawn structure, and after it no card
     body may still be showing one, while every card still offers a sketch area to click into. Read
     it again about ten seconds later and it must still be clear — the clearing happens within a
     second or so, so a single reading taken straight after the reset says nothing about whether it
     lasts. A card body whose drawing area cannot be read at all is a FAILURE here, not a card
     showing nothing: an unreadable area and an empty one look the same to the eye and must not to
     the check. The regression being guarded (GROK-14028 family): the sketcher input line persisted after reset even
     when the row count returned to full. Do NOT require the emptied card to show its **Sketch**
     placeholder again — the card that carried the structure comes back blank instead, which is the
     open bug GROK-20739; record what each card body ends up showing, and do not fail on it.
   - Verify every card is back in the default **Contains** search mode — the non-default mode set
     in Step 7 must not survive the reset either.

9. github-3445 — scaffold rotation and realignment:
   - Rebuild the panel the same way Step 7 does — **Remove All**, close the Filter Panel, reopen it
     — and verify one substructure filter per molecular column came back, exactly one of them on
     **Structure**, so each half below enters its own probe from a clean card on a known panel.
   - Note any console errors present before the realignment step (to establish a baseline
     for comparison afterwards).
   - The rotation of the drawn scaffold inside the sketcher canvas is not performed — see
     scope_reductions **SR-01**. The two claims below stand in for it.
   - **Two separate claims are checked here, and they must not be run together.** The sketcher's
     third checkbox, **Filter as you draw**, gates the second one and not the first: it decides
     whether a NEW or CHANGED sketch reaches the grid before **OK**, and it has no say over the
     repaint of a structure that is already applied.
   - **Claim 1 — the realignment repaints on the spot.** With the substructure filter already
     applied and Align and Highlight OFF, turn Align and Highlight ON with the sketcher dialog
     still open. Verify the rendered structure column repaints immediately — a real, non-zero
     pixel change on the grid's own data canvas — and that the row count does not move. This is a
     positive assertion, not a "nothing must change": the repaint is the product doing its job, and
     it happens in either state of **Filter as you draw**. Run this claim in BOTH halves, once with
     the checkbox cleared and once with it checked, each over its own proven-still grid: "it is not
     gated by Filter as you draw" is only shown by measuring it in both states.
   - **Claim 2 — what the checkbox really gates.** Once the filter is applied and its structure is
     on screen, make a FURTHER edit in the sketcher: enter a DIFFERENT structure from the one
     already applied.
     - With **Filter as you draw** cleared, that edit must NOT reach the grid while the dialog
       stays open: no repaint of the structure column and no change to the row count, checked
       again after the window in which the checked half lands, so a slow update cannot pass as a
       withheld one. Confirming with **OK** must then apply it — the row count moves to a
       different, still-filtering value.
     - With **Filter as you draw** checked, the same further edit must reach the grid LIVE, with
       the dialog still open, and must land on exactly the row count the cleared half reached only
       after **OK**. Verify it repaints the grid as well as moving the row count, so "it reached the
       grid" is not a bare number. The two halves differ in timing, not in outcome.
   - The cleared half's OPENING probe is the same claim in miniature and is checked the same way:
     entered into a clean card with the checkbox already clear, it must leave the row count at 100
     until **OK** is pressed.
   - In each half the checkbox is set to the state under test **before the probe molecule is
     entered**, and the value is read back and confirmed — the checkbox is sticky across dialog
     openings, so an inherited state would make the observation meaningless. Setting it after the
     molecule has already been committed is not the case under test.
   - In each half, before measuring anything, verify the substructure filter is engaged (row count
     below 100) and still keeps rows on screen (above zero) — an empty grid can show neither an
     edit nor a repaint — and verify both halves start from the same row set.
   - Start each CLAIM — not just each half — from Align and Highlight OFF, and read both checkboxes
     back as cleared before measuring. Claim 1 leaves them ON, so a Claim 2 that ran straight after
     it would attribute Claim 1's alignment repaint to the withheld edit and report a product
     failure that is not one.
   - Before each measurement, take a negative control: let the grid go quiet, snapshot it, idle, and
     verify the pixel delta is exactly zero. Only after a SUSTAINED zero does a non-zero delta
     afterwards mean the gesture caused it. An unreadable canvas must fail this control, not pass it.
   - In both claims verify that no new console errors appeared compared to the baseline noted
     before the action (compare by error message text, excluding known ambient noise by name).

10. GROK-14952 — the Filter Panel card and the column-popup filter on the same molecular column
    must not diverge:
    - Take down every substructure card in the panel through each card's own **X** icon and confirm
      none is left, on the saved-state channel and in the DOM, and that the row count is back at 100
      and stays there. Both surfaces are then armed from nothing.
    - Hover the **Structure** grid column header until its **column options** icon (the small bars
      icon at the right of the header) appears, and click it. A popup opens on that column with a
      **Filter** section holding a substructure filter of its own.
      - Verify the icon that was clicked really belongs to the Structure column (it sits inside that
        column's header band) and that the popup is titled **Structure**.
      - Verify the popup holds exactly ONE substructure filter body and that it is OUTSIDE the
        Filter Panel, while the panel holds none and the panel's filter group carries no Structure
        substructure state. Without this the two readings below could be one filter read twice.
    - Type the probe SMILES into the popup filter's own sketcher field and commit it with Enter.
      Verify the row count drops strictly below 100 and stays above zero, and that it HOLDS on a
      re-read — a count equal to 100 on both surfaces would satisfy "the same" while proving nothing.
      Note this count.
    - Click the popup's own **Add filter** link, which moves the criterion into the Filter Panel.
      - Verify the popup's substructure filter is GONE from the page and the panel gained exactly
        one card, carrying that structure and exactly one Structure substructure state. This is what
        makes the next reading the panel's rather than the previous one taken again.
      - Verify the row count still holds strictly between 0 and 100 and equals the count the popup
        filter produced.
    - Dismiss the popup with a neutral click on empty space beside it, and verify it is really gone
      before concluding anything: verify the popup is on screen first, click a point worked out from
      where the popup actually sits rather than a fixed spot on the screen, verify the popup has
      disappeared, and only then verify the row count did not move. A click that dismissed nothing
      leaves the row count unmoved for the wrong reason.
    - Switch the panel card OFF through its own enable checkbox and verify the card reads disabled
      and the row count returns to 100 and STAYS there. A count still below 100 means a popup-side
      filter is still contributing its mask even though the panel card is off — the divergence
      GROK-14952 reported. Re-arm the card and verify the row count returns to the value noted above.

11. Dragging the molecular column header onto the panel adds an empty card at the top:
    - Rebuild the panel from scratch (**Remove All**, close the Filter Panel, reopen it) and then
      take down the substructure cards it auto-populated, so the **Structure** column carries no
      substructure state at all. A drop onto a column that already has a card would update that
      card and add nothing, and neither claim below could be observed.
    - Make sure the panel still holds at least one card on some other, non-molecular column — add
      one through the panel's column selector if it does not. On an empty panel the dropped card is
      first simply because it is the only one, and "it landed at the top" would be true by
      arithmetic rather than because the product placed it there.
    - Drag the **Structure** grid column header onto the Filter Panel with a real pointer drag:
      press on the header, move toward the panel until the panel puts up its **Add filter** drop
      zone, and release inside that zone.
    - Verify the panel gained exactly one card, that the new **Structure** substructure card sits
      at the **top** of the panel rather than appended at the bottom, and that it arrives
      **empty** — no structure in its saved state.
    - Verify the row count is unchanged from before the drag. A freshly dropped, empty filter must
      not narrow the table; a card that arrives already carrying a criterion has filtered the data
      before the user drew anything.

12. A substructure filter and a non-molecular filter intersect, and each toggles independently:
    - Commit the probe structure into the **Structure** card and note the row count it filters to.
      Confirm it is below 100 and above zero.
    - Choose the second criterion **from the data**: a category on some non-molecular column that
      keeps part of the substructure hits, drops another part, and also covers rows the
      substructure filter rejects. Picking a category blindly risks one whose intersection is the
      whole of either filter, which no product behaviour could violate.
    - Add a filter card for that column and click the chosen category row on the card.
    - Verify the visible rows are exactly the rows satisfying **both** criteria — the pre-computed
      intersection count — and that this count is strictly smaller than the substructure filter
      alone **and** strictly smaller than the category filter alone.
    - Switch the substructure card off through its own checkbox and verify the row count rises to
      exactly the category's own row count; switch it back on and verify the intersection returns.
      Each filter must act independently of the other.

13. A disabled substructure filter survives collapsing the panel, still disabled:
    - With the **Structure** card filtering and carrying a drawn structure, switch it off through
      its own enable checkbox and note the row count it returns to.
    - Close the Filter Panel and reopen it.
    - Verify the card is back, still reads **disabled** on both its own enable checkbox and the
      card itself, and still carries the structure that was drawn before. The card merely being
      present is not enough — a card restored switched **on** has silently re-applied a criterion
      the user had turned off.
    - Verify the row count still matches the switched-off count, so the reopen did not re-apply the
      criterion.
    - Switch the card back on and verify it filters to exactly the row count it held before the
      panel was collapsed, which is what shows the surviving criterion is the one that was saved.

Expected:
- Step 3: fewer than 100 rows pass and more than zero; the Filter Panel header's active-filter
  counter is present, carries text and reads exactly 1 — the count of cards that are filtering, not
  of rows that pass; and the Structure card's saved state carries the committed structure.
- Step 4: the card's search-type control offers all six modes, and every count read lies between 0
  and 100. Contains passes exactly the row count Step 2 measured for the same probe in the same
  mode, itself strictly between 0 and 100. The Contains row count plus the Not-contains row count
  equals exactly 100, and so does Included-in plus Not-included-in (partition correctness). Exact is
  a subset of Contains and of Included in, Similar includes every Exact match, and restoring
  Contains reproduces the Contains count.
- Step 5: the panel holds exactly one new filter card per column checked in the picker, relative to
  Step 1, and no column holds two (GROK-18383 guard); the scaffold reaches the card through the
  card's own plus icon, its Add new scaffold... sketcher and its Add button, and the tree node the
  click selects reads checked afterwards; and, with the Structure substructure card switched off,
  the scaffold card on its own holds the row count strictly below 100 while keeping at least one
  row.
- Step 6: with the panel's substructure cards taken down first, and the current cell after the
  right-click reading the cell that was clicked — the grid's first visible row of Structure, pinned
  to that cell's own table row rather than to a fixed index, since the Step 5 scaffold-tree
  criterion is still filtering and the top visible row is not table row 0 — Use as filter leaves the
  panel holding exactly one more filter
  card, not two (GROK-20001 guard); the column carries exactly one substructure state; the new card
  carries the clicked cell's structure; and the row count is below 100 and above zero.
- Step 7: the clone made through View > Layout > Clone View is a second view holding its own filter
  group and its own substructure state, and its substructure card carries the original's Similar
  mode, cutoff, fingerprint and molecule (GROK-18530 guard). The two views' row counts are not
  compared — see the Automation notes.
- Step 8: immediately before the reset the panel is filtering, a card holds a drawn structure, a
  card body is painting it, and a card is in a non-default search mode. After the reset all 100 rows
  pass again, the Filter Panel header's active-filter counter reads 0, the card count is unchanged,
  no card's saved structure declares any atoms any more, every card is back in the default Contains
  mode, and no card body still paints a structure — still true when the bodies are read again ten
  seconds later, with every body still offering a sketch area (a canvas or the Sketch placeholder)
  and every canvas that is still there readable.
- Step 9: both halves reach a filtered grid (row count below 100 and above zero, the same row set in
  each), with "Filter as you draw" set before the probe is entered, Align and Highlight off, and the
  grid pixel-stable beforehand. Then, as two separate claims: turning Align and Highlight on
  repaints the rendered structure column immediately, with the dialog still open, while the row
  count stays put — in both checkbox states (Claim 1); and a further, different structure entered
  into the already-applied sketcher is withheld from the grid — no repaint, no row-count movement,
  held across re-reads — until OK when "Filter as you draw" is cleared, but reaches the grid live
  with the dialog open when it is checked, landing on the same row count either way (Claim 2). The
  cleared half's opening probe is withheld in the same way and on the same two channels. No new
  console errors appear in either claim. The scaffold rotation itself is out of scope, see SR-01.
- Step 10: the popup filter and the panel card are different objects — the popup body sits outside
  the Filter Panel, which holds no substructure card and no substructure state of its own; each
  surface's row count is sustained and strictly between 0 and 100; the two counts are equal;
  switching the panel card off returns the row count to 100 and holds it there (GROK-14952 guard),
  and re-arming restores the filtered count.

## Automation notes

- **Keep-both pair with `bio-filters.md`.** Both realize `filters.cp.chem-and-bio-filters`, so a
  coverage pass will propose consolidating them; the disposition is KEEP BOTH. They were split because
  Gate B has no skip channel — a skipped Chem half beside a passing Bio half in one spec file reports
  the whole file as PASS, which is how this scenario ran zero steps for weeks while its file was green.
- **Setup — gate Chem BEFORE opening the panel.** On a build where Chem never registered no filter card
  ever appears, so opening the panel first burns the shared 15 s first-card barrier and hard-fails
  inside a helper before the skip can name which of the four readings was missing. Open the table
  without the panel, take the readings (all four come off the dataframe), then open the panel. The
  applicability reading is an awaited `apply()` on the `Chem:substructureFilter` factory — that call
  both pulls the package in and proves it initialized; `pkg.load({file: 'package.js'})` never resolved
  on dev. Report "Chem registered, panel empty" as its own named failure, never as "Chem is not
  installed".
- **Painting inside the sketcher canvas has no DOM row to target**, so every molecule and every
  scaffold reaches its filter through a real input field instead: the substructure sketcher's own
  SMILES field, and — for the scaffold tree — the card's plus icon, its **Add new scaffold...**
  sketcher and its **Add** button, followed by a click on the tree node. `fg.updateOrAdd` is NOT a
  stand-in for any of those gestures; it appears nowhere in the spec. The gestures under test are
  the search-type switch (Step 4), the scaffold selection (Step 5), the clone command (Step 7) and
  the sketcher edits (Step 9), all driven through the UI. `fg` is shorthand for
  `grok.shell.tv.getFiltersGroup()`.
- **The header's active-filter counter is read through a header-scoped selector, and polled.** The bare
  indicator class also matches every card's own indicator, so the selector is scoped to
  `.d4-filter-group-header`. The counter is debounced off the filter events, so it is polled until its
  text settles rather than read in the turn that changed the criteria, and it is asserted PRESENT with
  non-empty text before its value is compared — an absent or blank counter proves nothing.
- **Step 9's in-canvas scaffold rotation has no Playwright primitive.** There is no element to grab, no
  rotation API on the saved state, and a raw mouse drag across the sketcher canvas draws a bond instead
  of rotating the fragment, so the gesture is not driven and the two claims below it stand in for it.
  The reduction is recorded as scope_reductions SR-01 with reason "no Playwright primitive for
  in-canvas scaffold rotation"; do not re-add it as a step a manual executor is asked to attempt.
- **Every per-mode row count in Step 4 waits for that mode's own recompute.** Idle-plus-stable is not
  enough on its own: a mode whose recompute is dispatched on a later task leaves the card's loader
  hidden and the previous count standing, and the reading would then be the previous mode's count
  carried forward into the subset and partition claims. The barrier is compound — the card's saved
  state must report the mode that was just picked, the table must have run at least one NEW filtering
  pass since the switch (counted off the frame's own rows-filtered event), and the count must then
  HOLD across a re-read a settle later. Measured on dev: picking the mode a card already displays
  fires no change and therefore no pass, so the new-pass half is required only when the switch really
  changes the control; that case arises once, for Contains at the top of the loop, whose count is
  separately pinned to the count Step 2 measured for the same probe in the same mode.
- **Search-type and cutoff must be driven as real input events**, not by assigning `select.value` and
  dispatching a synthetic `change`: use Playwright's `selectOption` on the card's search-type select
  and `fill` on the cutoff editor. A synthetic assignment bypasses the very control the step names.
- **The card's enable checkbox is a real click too, and it needs the card HOVERED first.** The card's
  `.d4-filter-controls` block is `display: none` until the pointer is over the card, so the checkbox
  measures 0x0 and Playwright refuses to click it — measured on dev 2026-08-22: unhovered, the click
  times out and nothing moves; hovered, the same click toggles the card in both directions and the row
  count follows (14 -> 100 -> 14 on the Step 2 probe). The checkbox carries the central claim of Steps
  5, 10, 12 and 13, so it is driven the way the user drives it rather than through an in-page
  `cb.click()`. Read the resulting state off the input and the card's class; a READ needs no hover.
  **The card's X remove icon lives in that same hover-revealed block** (`viewers/filters.md:59`), so
  the card-removal helper hovers each card and clicks its X for real, one card at a time, asserting
  the panel's substructure-card count fell by exactly one per click. An in-page click on the icon
  reaches the handler without the icon ever being clickable, which proves the handler rather than
  the gesture Step 1 names — and that helper is also the entry state for Steps 6, 10 and 11.
- **The three sketcher checkboxes — Align, Highlight and "Filter as you draw" — are real clicks too,
  and they need NO hover.** They sit in `.chem-sketcher-filter-options` inside the sketcher dialog and
  are 13x28 and visible as soon as the dialog is up, so a plain Playwright click lands on all three.
  Measured on dev 2026-08-22, at the level of the consequence rather than the checkbox state: with
  "Filter as you draw" CLEARED by a real click the opening probe was withheld (100 with the dialog
  open, 14 after OK) and a further edit was withheld too (held at 14, 32 after OK); with it CHECKED by
  a real click the further edit propagated live with the dialog still open (32 -> 14); and turning
  Align + Highlight on by real click repainted the grid by 9070 px over a proven idle-zero control
  while the row count stayed at 14. The product registers the real click on all three, so none of them
  needs the in-page `cb.click()` the spec used to run.
- **How the Steps 7/9 preconditions rebuild the panel, and why not through the picker.** Driving
  Add Filter > Substructure Filter... on a panel that already carries substructure cards does not
  top the panel up. The Chem factory declares `allowMultipleFiltersForColumn: 'false'`
  (`public/packages/Chem/src/package.ts:278`), and the picker's OK handler
  (`core/client/d4/lib/src/viewers/filters/filters_core.dart:676-685`) returns out of its whole
  column loop at the FIRST checked column that already carries a filter of that type - the columns
  after it in column order get nothing, so a picker-driven rebuild leaves a partial,
  column-order-dependent panel. The preconditions therefore use the auto-population path instead:
  **Remove All, close the Filter Panel, reopen it**, which brings every eligible column back
  represented by exactly one card without the picker being involved at all (the same shape
  `add-remove-entry-points.md` rebuilds its panel with). The rebuild is UNCONDITIONAL and the
  configuration it produces is then asserted - one substructure state per molecular column, one
  substructure card body per molecular column, exactly one of them on `Structure` - so the panel
  Step 8 resets is the per-molecular-column panel the scenario describes rather than whatever the
  earlier steps left behind.
- **Measured on dev 2026-08-18: the reset DOES clear the drawn structure — read the card body on
  PAINTED PIXELS, never on the presence of a canvas element.** On the full per-molecular-column panel
  (seven substructure cards, a structure drawn on `Structure` only) the reset returns the row count to
  100, restores every card to Contains, empties every saved state, drops the header's active-filter
  counter to 0 — and the `Structure` card body stops showing its structure too. The per-card
  measurement: before the reset the `Structure` body carries a 204×102 canvas with 313 non-background
  pixels and a saved SMILES of 12 characters; after the reset the same canvas is 100×50 with ZERO
  non-background pixels and the saved state is a 75-character empty molblock declaring no atoms. The
  other six bodies hold their "Sketch" placeholder and no canvas throughout, before and after.
- **How fast the clearing is, and how the hold is measured — re-measured on dev 2026-08-22.** The
  clearing is effectively immediate: the first clean read came 50 ms after the reset click and three
  consecutive clean reads at 1074 ms, so the watch ends about a second in. The 60 s figure is a
  CEILING on that watch, not a duration anything is observed for, and an earlier version of this note
  read as though the cleared state had been watched for a minute. The persistence claim is therefore
  measured separately: the bodies are re-read once, 10 s after the watch ends, and must still be
  clear. Across 87 reads over 45 s, and again at +50 s, +60 s and +75 s, all seven bodies stayed
  clear, so 10 s is comfortably inside observed behaviour while still costing one wait.
- **An unreadable canvas must FAIL the clearing check, not pass it.** The per-body ink count is -1
  when the canvas cannot be measured (0×0, or no 2d context), and "no body is painting" is read as
  ink > 0 being nowhere — so -1 would satisfy it. The reset resizes this canvas (204×102 -> 100×50),
  so a body caught mid-teardown is a real state rather than a hypothetical. A body that still carries
  a canvas but cannot be read is asserted against directly, and such a sample is also excluded from
  the clearing watch's three-consecutive-clean run, so it can neither end the watch early nor pass
  the final reading. The sketch-area guard does not cover this: it counts a canvas ELEMENT, and a 0×0
  canvas is still an element. No unreadable sample was seen in the 2026-08-22 measurement; the guard
  is there because the failure it catches is silent.
  An earlier version of this note recorded the step as RED on a never-clears product defect, with 60 s
  persistence figures. That reading was wrong: it came from a channel that called a body "drawing a
  structure" whenever a `.chem-external-sketcher-canvas` ELEMENT was present, and that element survives
  the reset unpainted. The clearing claim is now read on pixels plus the saved state, and it passes.
- **What the reset really leaves behind on `Structure` is the missing "Sketch" placeholder — that is
  GROK-20739, an OPEN bug, and it is carried under the known-open-bug mechanism rather than asserted.**
  After the reset the `Structure` body holds a blank 100×50 canvas and NO "Sketch" link, so it reads
  visually empty with no visible affordance, while its six never-armed siblings show their placeholder
  again. Mechanically, the empty branch of `updateExtSketcherContent` (`js-api/src/chem.ts:432-451`) is
  what puts the link back and it does not win the last pass; the non-empty branch re-appends the (now
  unpainted) canvas instead. Full anatomy in `bug-library/filters.yaml`. Because the bug is open, the
  spec must NEVER assert that every body is back to its "Sketch" placeholder — that clause would be a
  known red, and it is a separate and much weaker claim than the clearing one in any case. What the spec
  does instead is LOG the per-card body shape on every run (placeholder vs blank canvas, with the ink
  count), so the residue stays visible and a fix shows up as a change in the log. GROK-20739 is
  deliberately kept OUT of `related_bugs`, as GROK-20731 is in `add-remove-entry-points.md`: open bugs
  are referenced in prose, not declared as regression guards. Do NOT make the Step 7/9 rebuild
  conditional again to avoid the configuration.
- **The grid's right-click context menu is driven by a local helper, not the shared panel-menu one.**
  Step 6 opens the cell's own context menu, which lands as a fresh popup on top of whatever is already
  on screen, so the helper scopes every lookup to the LAST popup; the shared Filter-Panel menu helper
  addresses the panel hamburger's menu instead and does not fit. The shared reset-icon selector IS
  imported rather than restated.
- **Step 6 must start from a column with no substructure card.** `Use as filter` is an `updateOrAdd` on
  the same column and type, so with a card present it updates rather than adds and the duplicate-card
  regression is unobservable. Remove the existing cards through their own X icons first.
- **Step 7 cutoff: drive 0.6, not the widget default 0.8**, so a dropped cutoff cannot read back as
  coincidentally right. Clone through the top menu leaf **View | Layout | Clone View**
  (`[name="div-View---Layout---Clone-View"]`, `layout.dart:97-114`) — it is the command the bug's own
  reproduction names. Do NOT substitute `grok.shell.addTableView(grok.shell.tv.dataFrame)`: that opens
  a second view without saving and restoring the layout, so the filter-state serialization GROK-18530
  lived in is never run and the clone's fields are the original's own values read back. The command
  attaches the new view to the SAME DataFrame by design, so the honest "these are two objects" guard
  is on the two views, their two filter groups and their two state objects, not on the DataFrame.
  For the same reason the two views' row counts are the same number by construction — they share one
  filter BitSet — so that comparison cannot fail and is not made.
- **Step 9 canvas reads must name the body canvas** — pass `canvasSelector: 'canvas[name="canvas"]'` to
  every quiet / snapshot / diff / wait-for-change call. The default first canvas under
  `[name="viewer-Grid"]` is a 10 px strip on which a realignment cannot appear, so a read taken there
  passes the withheld half vacuously and fails the repaint half at 0.
- **Step 9 console baseline:** capture error texts with `page.on('console', ...)` before each action and
  compare by message text afterwards, excluding ambient noise by name. Prove the channel is live once,
  at the top of the run, by emitting a marker error from the page and asserting the listener received
  it — four "no new console errors" clauses are all satisfied by a listener that receives nothing.
- **Step 10 popup surface.** The column-options icon is `[column_name]` inside the grid root
  (`grid_column_tools.dart:19-49`); it is display:none until the pointer enters that column's
  header cell, so hover the header first and only then click the visible one. The popup it opens is
  `.d4-popup-host`, an accordion whose **Filter** pane renders `FilterControl.createPopupFilter`
  (`column_meta.dart:66`) — a filter instance held in `GridCore.popupFilters`, NOT in the panel's
  filter group, which is why the two surfaces can be told apart in the DOM. The popup's chem filter
  carries its sketcher INLINE (`.grok-sketcher-input input`), with no dialog. Its **Add filter**
  link transfers the saved state into the filter group and kills the popup filter
  (`column_meta.dart:76-88`) — that teardown is the GROK-14952 fix, and the disable half of the step
  is what fails without it.
- Teardown: close probe tables in a `finally`, even on failure.
- Selectors, dialog anatomy, the search-type control, the busy-loader readiness barrier, saved-state
  field names, cell-menu geometry, and what "Filter as you draw" does and does not gate all live in
  `.claude/skills/grok-browser/references/viewers/filters.md`.

### Spec must keep

The Automator re-authors this spec from scratch on every cycle. Each item below was established by a
live run and regressing it re-opens a defect that was already paid for.

1. **The column-picker flow in every card-adding step** — leaf → picker → All → OK → picker
   confirmed gone. Not `fg.updateOrAdd` as a substitute for the menu gesture.
2. **Dialogs addressed unambiguously** — picker by `name=`, sketcher by its own content
   (`.d4-dialog:has(input[placeholder*="SMILES" i])`). No bare `.d4-dialog` anywhere.
3. **Every card read and write scoped to a named column** (Structure), never "the first
   `.chem-filter`".
4. **Step 1 must be able to fail**: remove the pre-existing substructure cards through each
   card's own X, prove a **held zero** (poll to zero on the state channel, wait, re-read, and
   check the DOM channel too), and only then drive the leaf and assert one card per molecular
   column — plus a delta assert that the card count rose by the number of columns checked.
   Asserting "one card per column" without clearing first is green by construction.
5. **"Carries a structure" tested as molblock-or-SMILES**, never as a raw atom count: an empty sketch
   is a well-formed molblock with a `0 0 0` counts line and passes a null check, while a probe entered
   through the SMILES field is stored as the raw SMILES and has no counts line to count atoms off.
6. **The test-level `chemReady` guard stays OUTSIDE every `softStep`, and it runs BEFORE the
   Filter Panel is opened.** A `test.skip()` inside a step wrapper is swallowed and the file
   reports a non-failure with zero assertions executed; a panel opened before the guard dies in
   the shared 15 s card barrier and never reaches the skip at all.
7. **The compound partition invariants in Step 4** — `Contains + Not contains == full`,
   `Included in + Not included in == full`, `Exact ⊆ Contains`, `Similar ⊇ Exact`. These are
   dataset-independent; do not replace them with absolute counts. (On spgi-100 the current probe
   yields Exact 0, Similar 0, Included in 0 — the invariants still hold and still discriminate.)
8. **Ambient console noise excluded BY NAME** in the Step 9 delta guard. Narrowing the filter to
   the bug's own signature makes the guard evergreen.
9. **Step 9 sets "Filter as you draw" BEFORE the probe is entered, in both halves, and measures on
   a filtered grid in both halves** (operator ruling 2026-08-17 — the checkbox is sticky across dialog
   openings, so an inherited state discriminates nothing). The engaged-filter reading is an assertion,
   not an assumption; an unfiltered grid must fail the step.
10. **Step 9 keeps its two claims apart.** Claim 1 (Align + Highlight repaint an already-applied
    structure immediately, dialog open, row count unmoved) is a POSITIVE pixel delta and must never
    be re-written as "must not repaint yet" — that form was measured wrong three times. Claim 2 (the
    checkbox gates a FURTHER edit, withheld until OK when cleared, live when checked) is what the
    checkbox actually controls, and its cleared half must be a HELD zero re-read after the window in
    which the checked half lands, on both the pixel and the row-count channel.

11. **Step 10 proves the two surfaces are distinct BEFORE comparing what they report**, and proves
    each count narrows. Comparing two counts that both equal the full row count, or comparing one
    number read twice, is satisfied by a build where only one filter ever exists. The disable half
    (panel card off must return the row count to the full one, sustained) is the only clause that
    fails when the popup filter survives the transfer, which is the bug itself.

12. **Step 8 reads the cleared sketch area on the CARD BODY, not only on the saved state**, and
    proves a drawn body was there first. The GROK-14028-family regression was a card that still
    SHOWED the structure while the row count had returned to full, so a saved-state-only check is
    blind to exactly the defect the step names. Read the body on PAINTED PIXELS — count the
    non-background pixels on the body's own canvas — never on whether a
    `.chem-external-sketcher-canvas` element is present: that element survives the reset unpainted, and
    the element-presence form of this check produced a deterministic false RED for three runs and a
    wrong product finding in these notes. Both ends are measured on the pixel channel: above zero
    before the reset, exactly zero after. The placeholder-vs-blank-canvas shape each body ends in is
    LOGGED, not asserted (see the Automation notes) — the armed card's missing placeholder is the OPEN
    bug GROK-20739, so asserting it green would be a known red, and asserting it as a failure here would
    re-import the same conflation between "no structure is drawn" and "the sketch area is gone".
    Two things about that pixel channel are load-bearing and must not be dropped by a re-author.
    **An UNREADABLE canvas must fail, not pass**: an unmeasurable canvas leaves the ink count at -1,
    which "no body has ink above zero" is satisfied by, so a body that still carries a canvas and
    cannot be read is asserted against directly and is excluded from the clearing watch's clean run.
    The element-presence sketch-area guard does not cover it — a 0×0 canvas is still an element.
    **And the persistence is measured, not asserted**: the watch ends about a second after the reset
    (measured 2026-08-22: first clean read at 50 ms, three consecutive at 1074 ms), so the bodies are
    re-read once 10 s later and must still be clear. Do not restate the 60 s ceiling on the watch as
    though it were an observation window the cleared state was watched across.

13. **Step 9 measures Claim 1 in both halves.** The claim is that the alignment repaint is NOT gated
    by "Filter as you draw"; measuring it only in the cleared half asserts the repaint happens, never
    that the checkbox has no say over it. The checked half's live edit is read on the pixel channel
    as well as the row count, and the cleared half's opening probe is itself a withheld-until-OK
    measurement.

14. **Step 7 clones through View | Layout | Clone View.** `grok.shell.addTableView` on the same table
    does not run the layout save-and-restore the bug lived in, so the clone's fields are the
    original's own values read back and the step cannot fail. The command shares the DataFrame by
    design; the two-objects guard belongs on the views, their filter groups and their state objects.

15. **Step 4 pins Contains before any partition arithmetic.** Every invariant in the step —
    both complement sums, both subset claims, the Similar-contains-Exact claim — is satisfied by a
    filter that matches every row, and the Contains count also bounds Step 5. Pin it to the
    narrowing count Step 2 measured for the same probe in the same mode.

16. **Every card-changing gesture is the card's own control.** The scaffold is added through the
    tree card's plus icon and Add new scaffold... sketcher and selected by clicking its node; the
    search type is switched with a real `selectOption`. `fg.updateOrAdd` and `select.value = …` are
    substitutions for the very surfaces GROK-18383 and Step 4 are about.

17. **The console listener gets a positive control.** A detached listener satisfies all four
    "no new console error" clauses; emit a marker error once and assert it arrives.

18. **The Steps 7 and 9 preconditions rebuild the panel UNCONDITIONALLY, through Remove All → close
    the panel → reopen it**, and then assert the configuration that came back: one substructure
    state on every molecular column, one substructure card body per molecular column, exactly one
    of them on Structure, and — for Step 7 — the Structure card armed by this step's own probe
    rather than by an inherited structure. Do not restore the old "rebuild only when the panel holds
    no substructure state" branch, and do not rebuild through
    Add Filter > Substructure Filter...: the picker's OK handler abandons the rest of its column
    selection at the first column that already carries a card (see the Automation note), which is
    why the auto-population path is used instead. The point of the unconditional rebuild is that
    Step 8's reset is measured on the per-molecular-column panel the scenario describes, guaranteed
    rather than inherited.

19. **The dead `.catch` fallback after `waitForCanvasChange` stays deleted.** The helper does not
    update its snapshot while polling, so the fallback recomputes the delta that just failed — it
    can only mask, never rescue.
