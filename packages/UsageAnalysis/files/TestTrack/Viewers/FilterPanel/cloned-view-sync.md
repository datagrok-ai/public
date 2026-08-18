---
feature: filters
realizes_atlas:
  - filters.cp.cloned-view-sync
  - filters.int.same-column-sync
realizes:
  - viewers.filters
priority: p2
target_layer: playwright
coverage_type: edge
related_bugs:
  - id: GROK-17569
    status: fixed
  - id: GROK-13582
    status: fixed
  - id: github-1984
    status: fixed
realized_as:
  - cloned-view-sync-spec.ts
scope_reductions:
  - id: SR-02
    check: E-TRACE-02
    rationale: |
      Setup Step 3 as written adds an AGE histogram card alongside RACE, and
      Scenario 1 Step 1 restates that ("RACE filtered, AGE and SEX cards
      present"). The spec's setup builds only the RACE and SEX categorical
      cards and asserts that card set exactly; the AGE histogram card is
      created later, inside Scenario 1 Step 10, where it is the different-type
      half of the type-scoping check and where its narrowing window is proven
      to bite. The cards named in Setup all exist by the time any step reads
      them, but the "AGE card present from Setup onwards" ordering of the
      source text is not covered by any test step.
    verdict_status: SCOPE_REDUCTION
  - id: SR-03
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Scenario 1 Steps 6-7 as originally written changed the SEX card's
      criterion by typing into its search input and then verified that "the row
      count agrees in both views". Two reductions apply. First, the criterion
      change is driven on
      the RACE categorical card instead: the SEX card the prose names is a
      column-free-text card, a filter type the product no longer builds (see
      the product-fact automation note), and RACE is the card the setup already
      armed with a real narrowing criterion. Second, the cross-view row-count
      agreement is deliberately NOT asserted — a clone is a second TableView
      over the same DataFrame and the same df.filter BitSet object, so
      comparing trueCount across the two views is A === A and no sync
      regression can falsify it. The step asserts instead the real row-effect
      delta of the criterion change, the original card's own criterion read
      back through getStates, and the clone's independently-computed
      isFiltering list. The Step 7 and Step 9 bodies no longer state the
      cross-view agreement as an expectation; this entry records why the
      original wording was dropped rather than realized.
    verdict_status: SCOPE_REDUCTION
  - id: SR-04
    check: E-TRACE-02
    rationale: |
      Scenario 4 as written runs the layout round-trip on the Scenario 1 clone
      view and takes the save/apply through the top-menu leaves View > Layout >
      Save to Gallery and View > Layout > Apply. The spec runs it on a
      purpose-built single view ("LayoutHost") and drives the round-trip
      through view.saveLayout() + grok.dapi.layouts.save/find + loadLayout. The
      round-trip semantics under test — panel restored, card set restored,
      per-card enabled state restored, trueCount back at the saved value after
      a proven perturbation — are all exercised, but neither the
      round-trip-from-a-clone context nor the two gallery menu leaves are
      covered by any test step.
    verdict_status: SCOPE_REDUCTION
  - id: SR-05
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Scenario 2 Step 3 as written asks for a NEW AGE filter to be added in the
      CLONE through the panel's column selector, producing "a second,
      independently-configured filter instance on the same column, left at its
      default missing-values setting". That premise is not constructible. Live
      measurement on dev: the picker opens, accepts AGE (the spec reads its
      search box back at the moment of the commit, so this is measured and not
      assumed), offers the row, and adds nothing — the clone's card count, its `fg.filters.length` and its
      AGE-card count all stay where they were. The same refusal occurs through
      a direct click on the picker's canvas row and through dragging the AGE
      grid header onto the panel (the drop was genuine: `body` gained
      `d4-drag d4-drag-column` and a `.d4-drop-zone` appeared). Source agrees:
      `core/client/d4/lib/src/viewers/filters/filters_core.dart:829-841` —
      `addDefaultFilter` dedupes on (column, default filter type) and merely
      calls `f.root.scrollIntoView()` when a match already exists, and all
      three entry points route through it. The step therefore asserts what the
      product does — the refusal itself, over a sustained window — and the
      github-1984 guard is reduced to: the missing-values option chosen from
      the AGE card's own menu SURVIVES the clone, read back from that same
      menu, with the row count in agreement. What is NOT covered is a second
      independently-configured instance on the same column being the thing that
      fails to overwrite the first.
    verdict_status: SCOPE_REDUCTION
  - id: SR-06
    check: E-TRACE-02
    rationale: |
      Scenario 1 Step 9 and Step 11, and Scenario 3 Step 4, set their
      categorical CRITERION through the platform API (the filter group's
      updateOrAdd) rather than by clicking category rows in a card. Scenario 3
      Step 4 names the gesture explicitly — "filtering column c from v1's Filter
      Panel". Step 7 is NOT in this reduction: it drives the criterion change
      through the RACE card's own category-name row, because the two actuation
      paths are observably different — a card click broadcasts the new criterion
      on the DataFrame event bus and every other card on the same (column,
      filter type) copies it, while updateOrAdd does not raise that broadcast at
      all. In Steps 9 and 11 the API call only ARMS a criterion whose mirroring
      is not what the step asserts (the enable/disable mirror, the locality of a
      card removal); in Scenario 3 Step 4 it arms the criterion whose
      NON-crossing to a second DataFrame is the assertion, so that step's
      isolation reading is the one that would change meaning if the click
      behaved differently. That reading has since been taken through the click
      as well, on dev 2026-08-18, against a fixture STRONGER than the step's own
      (the second frame was given its own categorical card on the very same
      shared column object, which the step does not build): the category-name
      click narrowed the host frame from 5850 to 3243 rows inside 82 ms, and the
      second frame stayed at all 5850 rows with its own card's criterion
      untouched at every 30 ms sample of a 6.1 s window. The two actuation paths
      agree, and the reason they must is that the criteria broadcast is fired on
      the SENDING frame's own eventBus (`data_frame.dart:31,50`), which no second
      DataFrame listens on however many column objects the two share. A
      categorical card's category list is a canvas grid
      (see the filters reference doc): a plain click on the checkbox column of
      an all-selected card trips the product's own inversion rule and ends up
      selecting only the clicked category, and Ctrl/Shift-click suppresses the
      inversion but also moves the row SELECTION, so neither is a clean
      criterion-only gesture that could set an ARBITRARY category set. Step 7
      needs only a single-category narrowing, which is exactly what the
      exclusive-select click produces, so it is driven that way; the steps that
      need a specific multi-category set are the ones substituted. What each
      step then asserts — the row-effect delta, the criterion read back through
      getStates on both the original and the clone, the isolation of the frame
      that merely shares a column object — is all read back from the product.
      What is NOT covered is the click that sets the criterion. The gestures
      that ARE driven from the DOM and not substituted: the Step 7 category-name
      row click, Clone View from the top menu, each card's own enable/disable
      checkbox, the card's X remove icon, the card's filter menu and its Missing
      values submenu, and the panel's column picker.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: "Scenario 1 - Step 5"
    expectation: >-
      The clone's own Filter Panel comes up with the same filter cards as the
      original and each card's product-computed (isActive, isFiltering) state
      matches the original's, with exactly the one narrowing card marked
      isFiltering. The Clone View menu leaf is proven actuated and adds exactly
      one distinctly-named TableView, whose name the rest of the run is
      addressed by rather than a hardcoded convention, and the clone's own
      [name="viewer-Filters"] root is present and laid out inside that view's
      own root element.
  - anchor: "Scenario 1 - Step 7"
    expectation: >-
      The criterion is changed in the original view by CLICKING a category-name
      row in that view's own RACE card — the exclusive-select gesture, driven on
      the card's canvas body in the original view's own root, not through the
      filter group's API. That click produces a real row-effect delta in
      df.filter.trueCount (a new, strictly positive value, not the one before
      the change) and narrows the original's own card to exactly one category,
      read back through getStates. The clone's OWN card is then interrogated for
      its OWN criterion, not merely for which card is filtering, and it must
      BECOME that same single category within a bounded wait: the criterion set
      by a card gesture is broadcast on the DataFrame event bus and copied into
      every other card on the same (column, filter type), so the clone's card
      follows. The wait is a poll that FAILS on timeout, not a settle followed by
      a snapshot, so a clone that never follows fails here and so does one that
      follows too late to be usable. Non-vacuity comes from the pre-state read
      before the click: the clone held the SAME, strictly WIDER category set as
      the original, so "the clone reads a single category" cannot be satisfied by
      a clone that never moved. The clone's independently-computed isFiltering
      list still names that same card as its single filtering one. Cross-view
      trueCount agreement is deliberately not asserted: see scope reduction
      SR-03.
  - anchor: "Scenario 1 - Step 9"
    expectation: >-
      Toggling one card's checkbox in the original view to disabled causes the
      same-named same-type card in the clone to follow (isActive == false in
      both, d4-filter-disabled present on both cards). The isActive mirror at
      filters_core.dart:280-291 fired across both FiltersCore instances on the
      same DataFrame. Both cards are read as enabled in both views BEFORE the
      toggle, so the mirror is a real change and not a pre-satisfied state, and
      the untouched RACE card is read as still enabled in both panels as the
      negative control — a panel that painted everything disabled would fail
      there instead of passing the positive half. The SEX card is narrowed to a
      single category before the toggle (it held every category until then, so
      switching it off could not have moved a row), and the disable is shown to
      lift that narrowing: df.filter.trueCount returns to exactly the value it
      held before the narrowing and is strictly above the narrowed one. That is
      a row-effect delta of the disable on the one shared BitSet, not a
      comparison of the two views' counts (SR-03).
  - anchor: "Scenario 1 - Step 10"
    expectation: >-
      A card on a DIFFERENT filter type for the same column in the original view
      does NOT follow the disabling of that column's card in the CLONE's filter
      group — the isActive mirror travels between the two groups over the shared
      DataFrame but is type-scoped and must not propagate across types. That
      card's isActive remains true and its contribution to df.filter.trueCount
      is unchanged. The different-type card is first shown to really narrow
      rows (its window drops trueCount strictly below the pre-card value and
      strictly above zero) and the cross-group toggle is shown to have really
      flipped isActive from true to false, so neither half of the "unchanged"
      read is the difference of two nils.
  - anchor: "Scenario 1 - Step 11"
    expectation: >-
      After the original's SEX card is switched back on and narrowed to one
      category (isFiltering true, row count noted), removing the SEX card in the
      CLONE through its own X icon leaves the clone's panel without that caption
      while the ORIGINAL still shows it, with isActive true, isFiltering true,
      the same single selected category, and df.filter.trueCount equal to the
      value noted before the removal. The re-arming is what makes the row-count
      half able to fail; the pairing with Step 9 is what makes the step distinct
      — enable/disable crosses between the groups, removal does not.
  - anchor: "Scenario 2 - Step 4"
    expectation: >-
      The Missing values option chosen inside the AGE card's OWN filter menu is
      still the selected one when that same menu is re-opened after the view has
      been cloned, and df.filter.trueCount still reads 5849. Both halves are read
      in the same step: the row count alone would be satisfied by a card whose
      menu state was silently reset but whose narrowing happened to coincide, and
      the glyph read alone would be satisfied by a menu that showed the right
      radio over a filter that had stopped applying. All THREE radios are read
      every time, not just the selected one — a one-sided clause would pass a
      regression that lit two of them at once. The read before anything is
      touched is the non-vacuity baseline: Keep missing values selected, the
      other two clear. The submenu is proven to open on HOVER (its nested
      container's computed display goes none to flex) and the option is proven to
      bite (5850 to exactly 5849, demog's AGE column holding exactly one missing
      value). The menu is never dismissed with Escape — Escape RESETS the card,
      unchecks it and pushes the count back to 5850, manufacturing a fake "the
      setting did not hold" reading — so each dismissal is a neutral click on
      inert chrome and is proven inert by re-reading the row count after it. The
      clone-side re-add of AGE asserts the product's REFUSAL: the picker is proven
      to open, proven to have ACCEPTED the name AGE (its search box is read back
      at the moment Enter commits the pick, and that text is what the picker
      commits), and the clone's card count, filter count and AGE-card count are
      unchanged across a sustained window, because addDefaultFilter dedupes on
      (column, default filter type) (scope reduction SR-05). The refusal is
      preceded by a POSITIVE CONTROL through the same picker: a column that has
      no card yet is picked first, the picker is read back as holding THAT
      column's name, the clone's card count is required to go up by exactly one,
      and the card that appeared is required to be CAPTIONED with the column
      that was asked for. The control column is deliberately the LAST card-free
      column rather than the first, so a picker that ignored the typed name and
      committed whichever row it had highlighted adds the wrong card and fails
      here. Without all of that, "the card set did not change" is also what a
      picker that commits nothing at all produces, and a picker that always
      commits its first row would pass the control while never offering AGE.
  - anchor: "Scenario 3 - Step 4"
    expectation: >-
      t1's df.filter.trueCount drops to the filtered value after filtering
      column c from v1's Filter Panel AND t2's df.filter.trueCount remains equal
      to t2.rowCount (the full row count) at EVERY sample of a sustained window
      taken after the host frame was filtered, so contamination arriving on a
      longer debounce than one settle cannot slip between two reads. The two frames are first shown to
      really share one column object (a frame built on a COPY of the column
      could not fail this check at all). Only the host frame of the
      FilterControl is narrowed; the other DataFrame that merely shares the
      column object is untouched. This is the fix verified by commit 1eb9170760
      — reading only "a table got filtered" is satisfied by the bug.
  - anchor: "Scenario 4 - Step 3"
    expectation: >-
      After re-applying the saved layout from the clone, the Filter Panel is
      present in the clone view, the card set matches what was saved (same card
      captions readable from .d4-filter-column-name), each card's checkbox state
      is restored, and df.filter.trueCount matches the value recorded at save
      time. Before the save, exactly one of the two cards is switched off and
      the product stamps that on both surfaces (getStates active false plus the
      d4-filter-disabled class on the card, with the other card true and
      undisabled), and switching it off is shown to have lifted its narrowing,
      so the saved trueCount is a real one-card value strictly above the
      two-card count and strictly below the full row count. Before the
      re-apply, the perturbation is proven: no [name="viewer-Filters"] is left
      in the document and df.filter.trueCount is demonstrably NOT at the saved
      value, so the restore assertion is able to fail. Panel presence after the
      re-apply is read BEFORE getFiltersGroup() is called, since that call
      would create the group and make "the panel is back" true by the test's
      own doing.
---

# Filters — Cloned View Synchronization

## Setup

1. Open the demog golden dataset (5850 rows). Name the DataFrame and its view "demog" explicitly —
   every by-name view lookup in this scenario resolves against that name, and the by-name resolver
   must throw when a name matches no view or more than one, so a drifted name reads as a failure
   rather than as a green run. Then clear any filter cards the view came up with and release the
   row filter.
   Verify before going on: the rename took, the table has 5850 rows, nothing is filtered
   (trueCount equals the full row count, so a filter leaked from an earlier run surfaces here), and
   the panel really is empty of cards — a teardown that threw on every card would leave cards behind
   and every "the clone inherited the same card set" check below would be about the wrong set.
2. Open the Filter Panel using the funnel icon in the ribbon toolbar.
3. Add a RACE categorical filter by clicking the column selector at the top of
   the Filter Panel and choosing RACE, then add an AGE histogram filter the same
   way. Both cards appear at the top of the panel.
4. Apply a filter on RACE by unchecking one category in the RACE card, so that
   the row count shown in the active-filter indicator drops below 5850 and the
   indicator reads 1. Note this filtered row count as the initial count.
5. Also add a column-free-text filter on SEX via the panel column selector.
   Leave it unconfigured so it is not yet filtering.
   Verify the setup state the whole of Scenario 1 compares the clone against is itself
   non-degenerate: both cards are present, the filtered row count is strictly below the full row
   count and strictly above zero, and exactly one card — the RACE one — is reported as filtering.
   The card set the spec asserts here is RACE and SEX only; the AGE card the prose above names is
   created later, at Scenario 1 Step 10 (scope reduction SR-02).

## Scenarios

### Scenario 1: Card state inheritance and same-column-sync on a cloned view

**Steps:**

1. With the panel configured (RACE filtered, AGE and SEX cards present but
   not yet filtering), click **View > Layout > Clone View** in the top menu.
   A new tab "demog copy" opens.
2. In the **clone view**, open the Filter Panel via the ribbon funnel icon.
3. Verify the clone's Filter Panel is open (the Filter Panel viewer is visible
   in the clone view).
4. Read the ordered list of column names shown on each card in the clone's
   Filter Panel.
5. Verify that the card column-name list and the enabled/disabled checkbox state
   for each card match the original view's Filter Panel exactly. Confirm the
   clone's active row count equals the initial count recorded in Setup Step 4.
6. In the **original view**, locate the SEX card (same column, same
   column-free-text type as in the clone). Change the SEX filter criterion —
   type a single letter into the filter's search input to narrow the matching
   rows. Verify the row count in the original view drops from the initial count.
7. In the **original view**, click a category name in its RACE card so the card narrows to that
   single category. Before the click, read both cards' criteria: the clone must be holding the same,
   wider set of categories the original had, otherwise what the clone shows afterwards would not be
   a move. Verify the click produced a real row-effect delta (a new, strictly positive filtered row
   count, not the value it held before) and that the original's own card now holds exactly the one
   category. Then read the **clone's** OWN card and wait, up to a bounded time, for its criterion to
   become that same single category — the clone follows a criterion set by a card gesture. The wait
   must fail if the clone never follows. Also read the clone's own (column, enabled, filtering) list
   and verify it still names that same card as its single filtering one. The two views' row counts
   are deliberately NOT compared with each other: a clone is a second view over the same table and
   the same row filter, so that comparison is A === A and no sync regression could fail it
   (scope reduction SR-03).
8. In the **original view**, first narrow the SEX card to a single category so it is really
   filtering — while it holds every category it constrains nothing, and switching it off could not
   move a single row, which would leave the mirror below asserted over a card that was doing
   nothing. Verify the narrowing bites: the row count drops strictly below the value before it and
   stays strictly above zero. Then disable the SEX card by unchecking its checkbox (the checkbox at
   the top of the SEX card). Verify the row count returns to exactly the value it held before the
   narrowing, strictly above the narrowed one — the disable lifted the card's criterion.
9. Switch to the **clone view**. Verify the clone's SEX card reads disabled too — the isActive
   mirror carried the change from the original's filter group to the clone's.
   Read both cards as enabled in both views BEFORE the toggle, so the mirror is a real change and
   not a state that was already satisfied. After the toggle, read the disabled state from each
   view's OWN panel root — a document-wide query can only ever see whichever panel happens to be
   laid out — and read the untouched RACE card in both panels as the negative control: it must
   still be enabled, so the positive half cannot be satisfied by a panel that paints everything
   disabled.
10. In the **original view**, add a SECOND card on a column of a DIFFERENT
    filter type — add an AGE histogram card via the panel column selector and
    set a range window on it so it actually narrows rows. Verify the card was really created
    (exactly one AGE histogram state, reported active) and that its window really bites — the row
    count drops strictly below the value before the card was added and stays strictly above zero.
    Note the row count.
    Then switch to the **clone view**, add an AGE categorical card there via the
    column selector, and uncheck that card's checkbox. Verify that card was created and read as
    enabled, and that the uncheck really flipped it from enabled to disabled — otherwise the mirror
    had nothing to carry and the check below is vacuous. Take the reference row count between
    adding the card and toggling it, so whatever the freshly-added card does to the row set is not
    charged to the toggle. Switch back to the
    **original view** and verify its AGE histogram card is still checked and the
    row count is unchanged — both against the reference taken just before the toggle and, as the
    second half, still strictly below the pre-histogram value, so a mirror that had lifted the
    window's narrowing would fail here too.

11. Removing a card is view-local, unlike switching it off. In the **original view**, switch the
    SEX card back on and narrow it to a single category so it is really filtering — removing a card
    that was not filtering could not move the row count either way, and the check below would hold
    even on a broken build. Note the row count. Then switch to the **clone view** and remove its SEX
    card with the X icon on the card, confirming the card really is gone from the clone's panel.

    Switch back to the **original view**. Verify: it still has its own SEX card; that card is still
    switched on and still filtering; its criterion is still the single category chosen above; and
    the row count is exactly the value noted before the removal. The failure this guards against is
    a removal that crosses to the other view and takes its criterion with it — the user removes a
    card in a clone and silently loses a filter in the view they were not looking at.

**Expected:**

- Step 11: after removing the SEX card in the clone, the original still shows a SEX card that is
  active, filtering, holding the same category, and the filtered row count is unchanged. Removal is
  scoped to the group it happens in, while the enable/disable mirror of Step 9 is not.
- Step 5: card set and per-card enabled/disabled state match between original
  and clone; clone's active row count equals the initial count.
- Step 7: clicking a category name in the original's RACE card produces a real row-effect delta,
  the original's own card reads back the single category, and the clone's own card BECOMES that same
  single category within a bounded wait, having held the wider set before — the same-column
  same-type criterion mirror fired without looping. The clone's independently-computed filter list
  still names that same card as its single filtering one. The two views' row counts are not compared
  with each other (SR-03).
- Step 8: narrowing SEX drops the row count strictly below the value before it and strictly above
  zero, and disabling the card returns it to exactly that pre-narrowing value.
- Step 9: clone's SEX card follows the disable — both cards are unchecked in
  both views.
- Step 10: the AGE categorical card disabled in the clone is NOT mirrored onto
  the original's AGE histogram card (different type, sync is type-scoped) — that
  card stays checked and its contribution to the row count is unchanged.

### Scenario 2: Missing-values option isolation between views (github-1984)

**Steps:**

1. Open the demog dataset with its Filter Panel and let the panel come up with the cards it
   creates by itself. Verify the starting state is the one the rest of the scenario reads
   against: 5850 rows, exactly one missing value in AGE (that single row is the whole row-count
   delta below), nothing filtered yet, and exactly one card captioned AGE in the laid-out panel.
2. Hover the AGE card so its filter-options indicator appears — it is hidden until the card is
   under the pointer — and click the indicator to open the card's own filter menu. Verify exactly
   one menu popup is open and that it carries a **Missing values** group.
3. Before touching anything, read all three **Missing values** options: **Keep missing values**
   must be the selected one and **Filter out missing values** and **Show only missing value** must
   both be clear. This is the baseline that makes the read-back at Step 4 mean something. All three
   are always read together — reading only the selected one would pass a build that lit two at once.
4. Hover the **Missing values** group (it opens on hover, not on click; verify it really was
   collapsed before the hover and expanded after) and click **Filter out missing values**. Verify
   the row count settles at exactly 5849 — the one missing AGE row and nothing else.
   Dismiss the menu by clicking inert chrome, never with **Escape**: Escape resets the card,
   unchecks it and puts the row count back to 5850, which reads exactly like the bug this scenario
   guards. Verify the dismissal was inert — the menu is gone and the row count is still 5849. Then
   re-open the menu and verify it now reads **Filter out missing values** as the selected option
   with the other two clear.

   Dismiss again and clone the view (**View > Layout > Clone View**). Verify the clone leaf was
   actuated, that exactly one new view appeared, and that the row count is still 5849.

   In the **clone view**, drive the panel's column selector with AGE. The product refuses to add a
   second AGE card: verify the picker really opened, and that the clone's card count, its filter
   count and its AGE-card count are all unchanged and stay unchanged over a sustained window.
   `addDefaultFilter` dedupes on (column, default filter type) and only scrolls the existing card
   into view, so a second independently-configured instance on the same column cannot be built
   here at all (scope reduction SR-05).

   Switch back to the **original view**. Verify the row count is still 5849, re-open the AGE card's
   menu and verify **Filter out missing values** is still the selected option with the other two
   clear, and that the row count is still 5849 while that menu is open. Dismiss neutrally.

**Expected:**

- Step 4: the original view's AGE card still reads **Filter out missing values** as its selected
  Missing values option in its own menu, and the row count is still 5849. Both halves are read in
  the same step — the row count alone would be satisfied by a card whose menu state was silently
  reset but whose narrowing happened to coincide, and the glyph alone would be satisfied by a menu
  that showed the right option over a filter that had stopped applying. The clone-side attempt to
  add a second AGE card is asserted as the refusal it is, not as a successful add.

### Scenario 3: Shared-column isolation across two independent DataFrames (GROK-13582)

**Steps:**

1. Open a second table that shares the SEX column with the already-open demog
   table. Open a grid view for this second table so it has its own table view.
   Prove the sharing before relying on it: write a probe tag through the first table's handle on
   the SEX column and read it back through the second table's handle. If the two tables did not
   really hold one column object — if the second table were built on a copy — they would be
   independent by construction and the "the second table was not filtered" check in Step 4 could
   never fail, whatever the filter code did. The probe is behavioural, so it does not rest on JS
   wrapper identity.
2. Open a Filter Panel in the **demog view** if not already open.
   Confirm both the demog view and the second table view each show their full
   row count (5850 rows each) in the active-filter indicator.
3. In the **demog Filter Panel**, apply a categorical criterion on SEX
   (uncheck one category). Verify the demog view's row count drops below 5850.
4. Verify the second table's row count is still equal to its full row count
   (5850) — the other table that merely shares the SEX column was NOT filtered.

**Expected:**

- Step 3: demog's active row count is below 5850.
- Step 4: the second table's row count is still 5850. Only the demog frame is
  narrowed; the second table is untouched. Both row counts must be read in the
  same step — reading only "a table got filtered" passes while the bug
  (GROK-13582) is present.

### Scenario 4: Layout round-trip from the clone view

**Steps:**

1. Return to the **clone view** from Scenario 1. Configure the clone's
   Filter Panel — apply an additional criterion on AGE by setting a numeric
   range window (min 30, max 60) via the AGE histogram card's range handles.
   Then uncheck the RACE card's checkbox, so exactly one card is switched off
   and one is left on. Note the clone's active row count after that.
   Verify before saving: the card that was switched off reads as off on both surfaces the product
   stamps it on (its own state reports inactive AND its card carries `d4-filter-disabled`), the card
   left on reads as on on both, and switching RACE off actually lifted its narrowing — the noted row
   count is strictly above the count with both cards active and still strictly below the full row
   count. Without that, "the checkbox state was restored" would be satisfied by a re-apply that
   restored nothing.
2. Save the layout from the clone view: top menu **View > Layout > Save to
   Gallery**, enter a probe name (e.g. "cloned-view-sync-probe"), confirm.
3. Perturb the clone: first widen the AGE range window so the active row count
   is no longer the value noted in Step 1, then close the Filter Panel using its
   close button in the panel titlebar, then add an unrelated viewer (e.g. Bar
   Chart from Add viewer) so both the live layout and the filtered row set
   genuinely differ from the saved ones.
   Verify the perturbation actually landed before re-applying anything: no Filter Panel root is
   left in the document, and the filtered row count is demonstrably NOT at the value noted in
   Step 1. Closing the panel and adding a viewer change the layout but leave the row set where it
   was, which is why the AGE window is widened first — otherwise the row-count comparison after the
   re-apply could never fail.
4. Re-apply the saved layout: top menu **View > Layout > Apply**,
   select "cloned-view-sync-probe".
5. Verify the Filter Panel is now open in the clone view, the card column
   names match what was saved, each card's checkbox state is restored, and
   the active row count equals the value noted in Step 1.
6. Delete the probe layout in teardown (finally block), even on failure.

**Expected:**

- Step 1: exactly one card reads as switched off on both surfaces the product stamps it on, the
  other reads as on, and switching it off lifted its narrowing — the noted row count is strictly
  above the both-cards-active count and strictly below the full row count.
- Step 3: the perturbation is proven before the re-apply — no Filter Panel root is left in the
  document, and the filtered row count is demonstrably not at the value noted in Step 1.
- Step 5: Filter Panel is present in the clone, the same card set is shown
  (column names match), the same checkbox states are restored, and the active
  row count matches the value noted at save time in Step 1. Panel presence is read before any
  call that would create a filter group, so "the panel is back" cannot be true by the check's own
  doing.

## Automation notes

Selectors, canvas geometry, the Missing-values menu anatomy, the card dedupe key and the
shared-DataFrame fact behind the clone are in
`.claude/skills/grok-browser/references/viewers/filters.md`.

- **`filters.int.same-column-sync` (Step 7)** needs both halves: the same-column same-type card
  mirrors the criterion, and a card of a DIFFERENT type on that column does NOT. Asserting only
  the positive half passes against a blanket broadcast.
- **The criterion mirror is raised by the card gesture, not by the filter group's API.** A click on
  a category-name row goes through `grid_filter_base.dart:111-125`, which fires
  `FILTER_CRITERIA_CHANGED` on the DATAFRAME event bus; `grid_filter_base.dart:71-81` listens on
  that bus in every other card and applies the sender's `saveState()` when the (column, filter type)
  match. `fg.updateOrAdd` reaches the same card through `add(updateIfExists: true)`
  (`filters_core.dart:369-385`), which calls `applyState` and requests a re-filter but raises no
  DataFrame-level criteria event — so a criterion set that way stays with the one card. Setting the
  criterion through the API and then asserting that the clone does not follow is therefore an
  artefact of the actuation, not a product fact. Measured on dev 2026-08-18: after a category-name
  click in the original, the clone's card read the new single category at the first sample taken,
  ~330 ms after the click, for every one of four categories tried; the same change made through
  `updateOrAdd` left the clone's card untouched over a 2.4 s window.
- **Scenario 3 Step 4 still arms its criterion through the API (SR-06), and BOTH actuation paths
  have now been measured and agree.** Its assertion is that a SECOND DataFrame merely sharing a
  column object is not filtered, and a non-crossing claim taken only through the weaker path would
  look the same whether the product isolates the frames or the API simply never broadcasts. Measured
  on dev 2026-08-18 through the category-name click, with the second frame given its own categorical
  card on the shared SEX column object so the criteria mirror had a card to land on: the host frame
  went 5850 → 3243 rows within 82 ms of the click, the second frame held all 5850 rows and its card
  held its full category set at every 30 ms sample of a 6.1 s window, and the two frames were first
  proven to hold one and the same column object. The isolation is a product fact, not an artefact of
  the actuation — the broadcast at `grid_filter_base.dart:111-125` goes to the SENDING frame's own
  `dataFrame.eventBus` (`data_frame.dart:31,50`), which the second frame's cards never subscribe to.
  Steps 9 and 11 are unaffected: they use the API only to arm a criterion, and assert the
  enable/disable mirror and the locality of a card removal, neither of which rides the criteria bus.
- **GROK-17569** (Steps 5 and 9): assert the per-card enabled/disabled state, not only the card
  set and the row count.
- **GROK-13582** (Scenario 3): read BOTH the demog row count and the second table's row count in
  the same step — "a table got filtered" passes while the bug is present.
- **github-1984** (Scenario 2): the option must be SET from the AGE card's own filter menu and
  READ BACK from that same menu. Setting it through `fg.updateOrAdd` is a different path and
  produces a divergence that looks like a reproduction but is not one. The bug does not
  reproduce, so nothing here is wrapped in `knownOpenBug()`.
- **Never dismiss the AGE card menu with Escape** — it disables the card and returns the row
  count to full while the setting stays set, which reads exactly like "the setting did not hold".
  Every dismissal in Scenario 2 is a neutral click on inert chrome, and the row count is re-read
  after each one so the dismissal is proven inert rather than assumed.
- **The AGE card is located by its caption, never by index** — a card added through the picker
  lands at index 0 (`appendToTop`) and shifts the rest. `[name="viewer-Filters"]` exists once per
  view after cloning, so every panel read is `:visible` or scoped to a named view's own root.
- **Do not re-attempt the "switch to categorical first, then re-add" route** for a second AGE
  card — the card it produces does not carry the premise the check needs.
- **`column-free-text` is a dead filter type.** Setup Step 5 and Scenario 1 Steps 6-9 still name a
  "SEX column-free-text" card in prose; the spec drives a SEX **categorical** card there, and
  Step 10 uses the AGE histogram + AGE categorical pair for its same-column two-type premise.
- **Active-filter indicator vs `isFiltering`** (Setup Step 4): the realized check reads
  `fg.filters[i].isFiltering`. Both are the product-computed `isActive && !allTrue` and the
  indicator is a debounced render of the same fact — a channel substitution, not a weakening.
- **Scenario 4 is realized as a single step** anchored `Scenario 4 - Step 3`: configure, save,
  perturb, re-apply and verify all happen inside it, because the perturbation proof and the
  restore assertion have to share the recorded save-time values.
- **GROK-14952** (molecular filter / panel sync) is realized by `chem-filters.md` Scenario 1
  Step 10, not here — it needs Chem init, a molecular dataset and the column-popup filter surface,
  none of which demog offers (operator decision, 2026-08-18).
  **GROK-18331** is deliberately absent from
  `related_bugs`: still In Progress, so it would have to be carried under the known-open-bug
  mechanism rather than asserted green.
- **Teardown:** close/delete the cloned view, the second table (Scenario 3) and the probe layout
  (Scenario 4) in `finally`, even on failure.
