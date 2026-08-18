---
feature: filters
realizes_atlas:
  - filters.cp.chem-and-bio-filters
realizes:
  - viewers.filters
priority: p2
target_layer: playwright
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
  - anchor: "Scenario 1, Step 1"
    expectation: >-
      spgi-100 opens with one Chem substructure card per molecular column
      already in the panel, so those cards are first taken down one by one
      through each card's own X icon and the panel is confirmed — and confirmed
      to STAY — free of any substructure card, on both the saved-state and the
      DOM channel, before the menu is touched. Then driving the Filter Panel
      hamburger menu Add Filter > Substructure Filter... opens the "Select
      columns..." picker; confirming it with every offered column checked (the
      picker offers exactly the molecular columns) makes each of those columns
      carry exactly one Chem substructure filter card — the Structure card among
      them — the panel's card count grows by exactly the number of columns the
      picker reported checked, and that count is recorded as the baseline that
      Step 5 measures its new-card delta against.
  - anchor: "Scenario 1, Step 2"
    expectation: >-
      The probe reaches the filter through the card's OWN sketcher, not through a seeded state: the card's
      sketcher is opened from the card body, the dialog is confirmed to expose its own SMILES field, and the
      probe is typed into that field and committed. The programmatic `fg.updateOrAdd` seeding named in the
      Automation notes is the fallback for an unavailable field, never the path this step takes.
  - anchor: "Scenario 1, Step 3"
    expectation: >-
      BOTH channels are read in the same passage. On the state channel,
      df.filter.trueCount drops below the full spgi-100 row count (100) after
      the substructure is sketched and committed to the filter, and stays
      strictly above zero (a probe matching nothing would leave the Step 4
      partition arithmetic empty). On the header channel, the Filter Panel
      header's active-filter counter — read only through the header-scoped
      selector, since the bare class also matches every card's own indicator —
      reports how many cards are FILTERING, not how many rows pass, and it reads
      exactly 1: the Structure substructure card that just took the probe is the
      only filtering card. The counter is debounced off the filter events, so it
      is polled until its text settles rather than read in the turn that changed
      the criteria, and it is asserted PRESENT with non-empty text before its
      value is compared — an absent or blank counter proves nothing. The
      Structure card's saved state carries the committed structure — tested as
      molblock-with-atoms or, when the probe was entered through the sketcher's
      SMILES field, a non-empty SMILES string, never as a raw atom count.
  - anchor: "Scenario 1, Step 4"
    expectation: >-
      Each of the six search types named in the step is confirmed to be OFFERED
      by the Structure card's own search-type control before it is selected, so a
      build that drops or renames a mode fails by name rather than silently
      selecting nothing, and each resulting count is a real number between 0 and
      the full row count. Contains is PINNED first, to exactly the narrowing
      count Step 2 already proved for the same probe in the same mode (itself
      asserted to lie strictly between 0 and 100): every invariant below is
      satisfied by a filter that matches every row — Contains 100 plus Not
      contains 0 sums to 100, Exact 0 is a subset of everything — so without
      that pin the whole step passes on a filter that filters nothing, and the
      Contains count it exports would leave Step 5's comparison unbounded. The
      load-bearing claims are the dataset-independent
      compound invariants: trueCount(Contains) + trueCount(Not contains) equals
      the full row count (100) for the same sketched substructure, and so does
      trueCount(Included in) + trueCount(Not included in) — the two complement
      pairs partition the dataset exactly. Exact is a subset of both Contains and
      Included in, and Similar at the default cutoff contains every Exact match.
      Restoring Contains at the end reproduces the Contains count measured
      earlier in the step. See scope_reductions SR-02 for the per-mode
      "the count moved" clause, which is dataset-dependent and is not asserted.
  - anchor: "Scenario 1, Step 5"
    expectation: >-
      GROK-18383 — after adding a Scaffold Tree Filter from the hamburger Add
      Filter group and confirming its "Select columns..." picker (which is itself
      confirmed to offer exactly the molecular columns, as the Substructure
      picker in Step 1 is), exactly one
      new scaffold-tree card per confirmed column is present in the panel (the
      card count grew by exactly the number of columns the picker reported
      checked, relative to the count recorded in Step 1, and no column carries
      two). The scaffold is then placed through the card's OWN controls, never
      through a seeded state: the Structure scaffold-tree card is addressed by
      column name and card kind (exactly one such card must exist), its
      "Sketch scaffolds manually" plus icon opens the card's own
      "Add new scaffold..." sketcher, the scaffold SMILES is typed into that
      sketcher's field and confirmed with Add, and the node it creates is then
      clicked in the tree — the click is what selects it, and the node is
      asserted to read CHECKED afterwards, since an unchecked node is not a
      criterion at all. With the Structure substructure card switched off so the
      scaffold card's own effect is what is read, that selection drops
      df.filter.trueCount strictly below the full row count
      while keeping at least one row. The switch-off is itself confirmed (the
      card reads disabled and its checkbox unchecked), so the "isolated"
      measurement is proven isolated rather than assumed. Selecting the scaffold
      UPDATES the card the menu added — the card count is unchanged by it, so no
      second card was created — and with the substructure card
      re-armed the AND of both cards passes no more rows than the substructure
      card alone (see scope_reductions SR-03 for why the strict "below the
      Contains count" form of the prose is not what is measured).
  - anchor: "Scenario 1, Step 6"
    expectation: >-
      GROK-20001 — the pre-existing substructure cards are first taken down
      through their own X icons and the Structure column is confirmed to hold no
      substructure state, because the action is an updateOrAdd on the same column
      and type: with a card already present it would update it and add nothing,
      leaving the duplicate-card regression unobservable. The current cell is
      NOT set through the API before the right-click: after the right-click and
      before the menu leaf is driven, the current cell is asserted to be row 0 of
      the Structure column, so a click that landed elsewhere fails here instead
      of being covered up by a pre-set that would make the structure carried by
      the new card right whatever was clicked. After using Use as
      filter on a molecular cell the filter card count grows by exactly 1 (not
      2), the Structure column carries exactly ONE substructure state, that card
      carries the clicked cell's structure, and df.filter.trueCount drops below
      the full row count while staying above zero (the cell's own molecule must
      at least match itself).
  - anchor: "Scenario 1, Step 7"
    expectation: >-
      The step opens by rebuilding the panel UNCONDITIONALLY through the
      auto-population path — Remove All, close the Filter Panel, reopen it — and
      the configuration that comes back is asserted, not assumed: every
      molecular column carries exactly one substructure state, the panel renders
      one substructure card body per molecular column, exactly one of them sits
      on Structure, and that card comes up with NO structure, so the probe this
      step enters is what arms it. That is the per-molecular-column panel the
      scenario describes, and it is what Step 8's reset is then measured on.
      GROK-18530 — the original view's substructure card is switched to Similar
      mode with a NON-DEFAULT similarity cutoff (0.6, against the widget default
      of 0.8) through the card's own controls, and that mode and cutoff are
      confirmed on the original before the clone is made. The clone is produced
      by the product's OWN command — the top menu View > Layout > Clone View,
      which saves the view layout and rebuilds a view from it; that
      save-and-restore of the filter state IS the code the bug lived in, and
      opening a second view over the same table bypasses it entirely, leaving the
      comparison below reading the original's own values back. Before any field
      is compared, the clone is proven to be a SEPARATE view holding a SEPARATE
      filter group and a SEPARATE substructure state object from the original's.
      After cloning, the
      clone's substructure card carries the same search type (Similar), the same
      similarity cutoff, the same fingerprint type and the same molecule as the
      original, and both cards display Similar. Row counts are NOT compared
      across the two views: Clone View attaches the new view to the SAME
      DataFrame and therefore the same filter BitSet, so that comparison cannot
      fail.
  - anchor: "Scenario 1, Step 8"
    expectation: >-
      Immediately before the reset the panel is confirmed to be filtering, to
      hold a drawn structure and to be in a non-default search mode. After
      clicking the panel reset icon, df.filter.trueCount returns to the full row
      count (100), the cards are all still present — the surviving card count
      equals the pre-reset count exactly, on the DOM channel as well as on the
      saved-state one, so reset clears criteria rather than removing cards — no
      card's saved molblock declares any atoms any more and every card is back in
      the default Contains mode. The "the sketch area is cleared, not merely
      detached from filtering" half is read on the CARD BODY too, which is the
      channel the GROK-14028-family regression showed on, and it is read as
      PAINTED PIXELS on each body's own canvas rather than as the presence of a
      canvas element — the element survives the reset unpainted, so element
      presence reports a structure that is not on screen. Before the reset at
      least one card body PAINTS a structure (non-background pixels above zero);
      after it no card body paints one, and the saved-state channel agrees that
      no card still holds a structure. Every body is confirmed to still expose a
      sketch area — a canvas or the "Sketch" placeholder — so the pixel reading
      is not taken off nothing. WHICH of those two an emptied body comes back as
      is recorded per card and deliberately NOT asserted: the armed card comes
      back blank rather than placeholdered, which is the OPEN bug GROK-20739 and
      is carried under the known-open-bug mechanism, logged on every run rather
      than turned into a known red. The panel this is measured on is the
      full one Step 7 rebuilt — a substructure card on every molecular column —
      and the clearing is WATCHED rather than merely waited on: the card bodies
      are sampled every 500 ms across a deliberate 60 s observation window and
      the cleared reading must HOLD across three consecutive samples, so a
      failure reports how many bodies on which columns were still painting a
      structure and for how long, instead of a bare barrier timeout. The header's active-filter counter is
      read again once the panel settles and must have fallen to 0 — Step 3 reads
      that counter only at 1, a value a stale or hardcoded indicator also
      produces, so the reset supplies the contrast that discriminates. The clone view opened in Step 7 is
      closed first and the shell is confirmed down to a single table view, so the
      reset acts on the original panel.
  - anchor: "Scenario 1, Step 9"
    expectation: >-
      github-3445 (scope-reduced — the scope_reductions entry records that the
      in-canvas scaffold rotation has no Playwright primitive and is not
      performed). TWO SEPARATE CLAIMS hold, and they are not conflated: the
      "Filter as you draw" checkbox gates one of them and not the other.
      CLAIM 1 — realignment. With a structure already applied and Align and
      Highlight OFF, turning them ON repaints the rendered structure column
      IMMEDIATELY, with the sketcher dialog still open, while
      df.filter.trueCount does not move. This is asserted as a POSITIVE pixel
      delta on the grid's own data canvas over a proven sustained-zero baseline,
      not as a "must not change", and it holds in either checkbox state — which
      is MEASURED, not assumed: Claim 1 runs once with "Filter as you draw"
      CLEARED and once with it CHECKED, each over its own sustained-zero control,
      each with the row count asserted unmoved and no new console errors.
      CLAIM 2 — what the checkbox really gates: the propagation of a NEW or
      CHANGED sketch, never the repaint of an already-applied structure. With
      "Filter as you draw" CLEARED, entering a FURTHER, DIFFERENT structure into
      the already-filtering sketcher does not reach the grid while the dialog
      stays open — a held zero, re-read after the window in which the checked
      half lands, on both channels: no repaint of the structure column and no
      movement of df.filter.trueCount — and confirming with OK then applies it,
      moving the row set to a different, still-filtering count. With
      "Filter as you draw" CHECKED, the same further edit reaches the grid LIVE
      with the dialog still open — on the pixel channel as well as the row-count
      one, so "it reached the grid" is not a bare number — and lands on exactly
      the row set the cleared half reached only on OK, so the difference between
      the two halves is timing and nothing else. The cleared half's OPENING probe
      is measured the same way: entered with the checkbox already clear, it leaves
      the row set at the full 100 until OK is pressed.
      In both halves the checkbox is set to the state under test BEFORE the
      probe molecule is entered and read back, never inherited (it is sticky
      across dialog openings); both halves start from Align and Highlight OFF
      and BOTH checkboxes are read back as cleared at the start of EACH claim,
      including Claim 2, which follows a Claim 1 that left them on — a repaint
      Claim 2 inherited from Claim 1's alignment would otherwise be reported as a
      product failure that is not one; both halves
      measure on a FILTERED grid with df.filter.trueCount asserted strictly
      between zero and the full 100, and start from the same row set; and no new
      console errors are logged in either (baseline captured before each action;
      errors compared by message text with ambient noise excluded by name).
  - anchor: "Scenario 1 Step 10"
    expectation: >-
      GROK-14952 — the same molecular column is filtered from two surfaces and the two must not
      diverge. The step first proves the two surfaces are DIFFERENT objects: with every
      substructure card taken down (held zero on both the state and the DOM channel) and the row
      count sustained at the full 100, the Structure grid column header's own column-options icon
      is hovered into view and clicked — the icon is confirmed to sit inside that column's header
      band and the popup is confirmed to be titled Structure — and the popup is confirmed to hold
      exactly ONE substructure filter body, OUTSIDE the Filter Panel, while the panel holds none
      and the panel's filter group carries no Structure substructure state. The probe entered into
      the popup's own sketcher field must then really narrow: a SUSTAINED count strictly between 0
      and 100, so an agreement between two counts that both equal the full row count cannot satisfy
      the step. The popup's own "Add filter" link then moves the criterion to the panel, and the
      popup filter is asserted GONE from the document while the panel gains exactly one card
      carrying that structure — which is what makes the second reading the PANEL's rather than the
      first reading taken again. That panel count is sustained, is itself strictly between 0 and
      100, and equals the popup count exactly. Finally the bug's own mechanism: switching the panel
      card off through its own enable checkbox (confirmed disabled) must return the row count to
      exactly 100 and HOLD there — a count still below 100 means a surviving popup-side filter is
      still contributing its mask, which is precisely the divergence GROK-14952 reported — and
      re-arming the card restores exactly the count it held before. The popup dismissal in between
      is proven to have HAPPENED before it is called inert: the popup is asserted present, the
      neutral click lands at a point computed live from the popup's own rectangle rather than a
      fixed viewport coordinate, and the popup is asserted gone afterwards — only then is the
      row count re-read and required to hold, so a click that dismissed nothing on a shorter
      viewport can no longer pass the clause by inaction.
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
     current cell really is the one that was clicked (row 0 of **Structure**). Otherwise a click
     that landed on the wrong cell still produces the expected structure and the check below tests
     nothing.
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
     body may still be showing one, while every card still offers a sketch area to click into. The
     regression being guarded (GROK-14028 family): the sketcher input line persisted after reset even
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
   - Rotate the drawn scaffold within the sketcher widget (this is a canvas gesture inside
     the sketcher — if this gesture cannot be performed in a headless browser environment,
     record a scope reduction with reason "no Playwright primitive for in-canvas scaffold
     rotation" and drive the two claims below instead). The rotation itself is not performed —
     see scope_reductions **SR-01**.
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

Expected:
- Step 3: fewer than 100 rows pass but more than zero (the filter is active), the Filter Panel
  header's active-filter counter — present, non-empty and settled — reads exactly 1 filtering card,
  and the card's saved state carries the committed structure.
- Step 4: every one of the six modes is offered by the card's search-type control; Contains is
  first pinned to the narrowing row count Step 2 proved for the same probe in the same mode (itself
  strictly between 0 and 100), so a filter matching every row cannot satisfy the arithmetic below;
  the Contains row
  count plus the Not-contains row count equals exactly 100, and so does Included-in plus
  Not-included-in (partition correctness); Exact is a subset of Contains and of Included in;
  Similar includes every Exact match; restoring Contains reproduces the Contains count.
- Step 5: the panel gained exactly one new filter card per column checked in the picker, relative
  to Step 1, and no column gained two (GROK-18383 guard); the scaffold reaches the card through the
  card's own plus icon, Add new scaffold... sketcher and Add button, and is selected by clicking the
  resulting tree node, which is confirmed checked; and, with the Structure substructure card
  switched off, the scaffold card on its own still holds the row count strictly below 100 while
  keeping at least one row.
- Step 6: measured from a panel whose substructure cards were removed first, and with the current
  cell confirmed after the right-click to be the cell that was clicked rather than one set through
  any other route, the panel gained
  exactly one new filter card, not two (GROK-20001 guard), the column carries exactly one
  substructure state, the new card carries the clicked cell's structure, and the row count is below
  100 and above zero.
- Step 7: the clone is made through View > Layout > Clone View and is proven to be a second view
  with its own filter group and its own substructure state before anything is compared; the clone's
  substructure card carries the original's Similar mode, cutoff, fingerprint
  and molecule (GROK-18530 guard); row counts are not compared across the views.
- Step 8: all 100 rows pass the filter again AND the Filter Panel header's active-filter counter
  reads 0 AND the card count is unchanged AND no card's saved
  structure declares any atoms any more AND no card BODY still DRAWS a structure — measured on what is
  painted in the body, not on whether a drawing area is still present, with every card confirmed to
  still offer a sketch area — AND every card is back in the default Contains mode —
  asserted after confirming, right before the reset, that the panel was filtering, did hold a drawn
  structure, was painting that structure in a card body, and was in a non-default search mode.
- Step 10: the popup filter and the panel card are shown to be different objects (popup body
  outside the panel, zero substructure cards and zero substructure states inside it) before their
  row counts are compared; each surface's count is sustained and strictly between 0 and 100; the
  two counts are equal; switching the panel card off returns the row count to 100 and holds it
  there (GROK-14952 guard), and re-arming restores the filtered count.
- Step 9: each half sets "Filter as you draw" before its probe is entered, reaches a filtered grid
  (row count below 100 and above zero, the same row set in both halves), and starts from a proven
  pixel-stable grid (sustained-zero negative control) with Align and Highlight off. Then, as two
  separate claims: turning Align and Highlight on repaints the rendered structure column
  immediately, with the dialog still open, while the row count stays put — measured in BOTH checkbox
  states, which is what shows the alignment path is not gated by it (Claim 1); and a further,
  different structure entered into the already-applied sketcher is withheld from the grid — no
  repaint, no row-count movement, held across a re-read — until OK with "Filter as you draw"
  cleared, but reaches the grid live with the dialog open when it is checked, landing on the same
  row count either way (Claim 2). No new console errors appear in either claim. The scaffold
  rotation itself is out of scope, see SR-01.

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
  non-background pixels and the saved state is a 75-character empty molblock declaring no atoms, and
  it stays that way for the whole 60 s observation window. The other six bodies hold their "Sketch"
  placeholder and no canvas throughout, before and after.
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
