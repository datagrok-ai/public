---
feature: filters
realizes_atlas:
  - filters.cp.hierarchical-and-combined-boolean
realizes:
  - viewers.filters
  - viewers.filters.bool-columns
priority: p2
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19968
    status: fixed
  - id: GROK-16528
    status: fixed
  - id: GROK-16488
    status: fixed
source_text_fixes: []
realized_as:
  - hierarchical-and-combined-boolean-spec.ts
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Step 7's expectation names "Caucasian" as the parent whose checkbox must
      read the indeterminate glyph. Under the SEX / RACE hierarchy configured in
      Step 5, RACE is the terminal level, so "Caucasian" is a LEAF and has no
      children to partially check — the three-state assert is unrealizable on
      that node. The spec asserts the same three-state contract one level up, on
      the SEX node "F": bring F to fully checked (3243 rows) and assert F reads
      the checked glyph, uncheck its Caucasian child, then assert the F
      substitute glyph reads the indeterminate codepoint U+F146 and the row count
      sits strictly between the single-descendant count and the full-F count. The
      tri-state semantics are preserved; only the node the glyph is read on
      differs from the wording.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-TRACE-02
    rationale: |
      Step 10 as written saves the layout AFTER the Step 9 column reorder. A
      reorder CLEARS the hierarchical criterion, so a save taken in that state
      captures nothing, while Step 11 expects the restored layout to carry
      caption "SEX / RACE" and the Caucasian-female count from Step 6 — the two
      steps as written are mutually inconsistent. The spec saves the layout in
      the coherent pre-reorder state instead (reset to SEX / RACE, check
      Caucasian only, 2823 rows, save), then performs the Step 9 reorder. Both
      numbered steps are exercised, but the save-in-post-reorder-state ordering
      of the source text is not covered by any test step.
    verdict_status: SCOPE_REDUCTION
  - id: SR-03
    check: E-TRACE-02
    rationale: |
      The hierarchical card's COLUMN CONFIGURATION is not driven through the
      card's column picker. Steps 5, 8b, 6-restore and 9 each set the hierarchy
      (SEX / RACE, then SEX / RACE / SEVERITY, then back, then RACE / SEX)
      through the platform API `window.grok_GridFilterBase_ApplyState`, followed
      by a row-filter request. The picker itself is a canvas-rendered column
      list whose row clicks do not register (see the filters reference doc), so
      the pick, the multi-select and the reorder cannot be actuated from the
      DOM. Everything the steps then assert — the card caption naming the
      levels in order, the tree rebuilding with the new root level, the row
      counts — is read back from the product. What is NOT covered is the
      picker gesture itself: choosing the columns and reordering them by hand.
      The tri-state checkbox toggles that carry Steps 6, 7 and 8b are NOT
      substituted; those are real clicks on the tree's own
      `input.d4-hierarchical-filter-checkbox` elements.
    verdict_status: SCOPE_REDUCTION
  - id: SR-04
    check: E-TRACE-02
    rationale: |
      The combined boolean card's per-column true/false checkboxes and its
      OR/AND mode switch are driven through the platform API
      `window.grok_GridFilterBase_ApplyState` in Steps 15, 16, 17 and 18,
      rather than by clicking the card. The card body is a canvas-rendered grid
      (see the filters reference doc), so neither the two checkboxes per
      boolean column nor the mode control is a DOM target. The effects the
      steps assert — the row count matching the union or the intersection
      computed from the two boolean columns' own data, the active-filter
      counter, the state surviving a project reopen and a layout round-trip —
      are all read back from the product. What is NOT covered is the toggling
      and mode-switching gestures themselves. The card's REMOVAL (its X icon),
      its creation through Add Filter > Combined Boolean, Remove All, and the
      panel close/reopen are all real DOM gestures and are not substituted.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: Step 5
    expectation: >-
      The hierarchical filter card caption reads "SEX / RACE" (two columns in
      order), confirming the column picker applied the selection.
      df.filter.trueCount equals grok.shell.tv.dataFrame.rowCount (no criterion
      checked yet).
  - anchor: Step 6
    expectation: >-
      The SEX="F" parent node is expanded and its RACE children are on screen.
      Under the SEX / RACE hierarchy RACE is the terminal level, so "Caucasian"
      is itself a LEAF and has no children — the same wording problem SR-01
      records for Step 7. The node is addressed by its caption path under F, not
      by [name="tree-expander-…"], which only ever names top-level nodes. After
      checking it, its own .d4-hierarchical-filter-checkbox-substitute renders
      the checked glyph U+F14A — read as a codepoint, because unchecked (U+F0C8)
      and indeterminate (U+F146) are BOTH font-weight 100 and a font-weight read
      cannot tell them apart. df.filter.trueCount is exactly 2823 — the
      Caucasian-female count, strictly positive and strictly below the full row
      count.
  - anchor: Step 7
    expectation: >-
      The SEX node "F" is first brought to fully checked (3243 rows) and asserted
      to render the checked glyph U+F14A, so the indeterminate read that follows
      is a transition and not a standing state. After unchecking its "Caucasian"
      child, the .d4-hierarchical-filter-checkbox-substitute on "F" renders the
      indeterminate glyph U+F146 and is explicitly not the unchecked U+F0C8 —
      the discrimination the earlier font-weight channel could not make, since
      unchecked and indeterminate are both font-weight 100 and an untouched node
      satisfied that read. df.filter.trueCount is 420, strictly above the
      descendant-only count and strictly below the all-F count — a three-state
      assert.
  - anchor: Step 8
    expectation: >-
      GROK-19968 — the tree search input is hidden until the card's search icon
      is clicked, and the typed term is confirmed by reading the input's own
      value back, so a keystroke that never reached the Dart onInput listener
      fails here instead of asserting against an untouched tree. After typing a
      search term, the visible tree nodes shrink to
      the matching nodes and their ancestors only: the matching caption
      ("Caucasian") is still visible while a non-matching sibling ("Asian") is
      gone, and the visible-node count is strictly lower than before typing.
      df.filter.trueCount is UNCHANGED from the pre-search value, which was
      itself a narrowed count (the tree search hides nodes, it does not filter
      rows). Clearing the search restores exactly the pre-search visible node
      set and the same row count.
  - anchor: Step 8b
    expectation: >-
      With the hierarchy set to SEX / RACE / SEVERITY, the card caption reads
      "SEX / RACE / SEVERITY" and checking the whole F branch gives the female row
      count (3243 on demog) and an F node with the checked glyph U+F14A. Every row
      count in this step is computed from the raw SEX / RACE / SEVERITY column
      values — never read back from df.filter, which would make the expectation
      circular — and each derived number is asserted to be a strict, non-empty
      subset before it is used. F and Caucasian are then expanded
      by caption path — a second-level expander cannot be addressed by
      [name="tree-expander-<value>"], whose name is built from the node's whole
      ancestor path text. The SEVERITY level under F / Caucasian is then asserted
      to list EXACTLY the SEVERITY values the data holds for that branch, derived
      from the raw columns like the counts and never read back from the tree
      (Critical, High, Low, Medium, None on demog today). Each generated level
      also renders one caption-less row — the bold column-header naming the level
      — which the tree reader reports separately instead of returning as an empty
      caption; its count and label are asserted (exactly one, reading "SEVERITY"),
      so a build that emits more such rows, or that strips the caption off a real
      value node, fails visibly rather than being filtered away. The load-bearing
      node is "Caucasian", asserted as a
      TRANSITION: it renders the checked glyph U+F14A BEFORE the grandchild
      "None" is unchecked and the indeterminate glyph U+F146 after. The uncheck
      is local — the level still lists the same values afterwards, every SEVERITY
      sibling in the derived set still reads U+F14A (the loop walks that proven
      set, so a sibling that vanished from the tree fails as a missing node rather
      than shortening the loop) and "None"
      reads the unchecked U+F0C8 — and df.filter.trueCount lands on the female row
      count minus the female / Caucasian / None rows (3243 - 1588 = 1655 on demog),
      strictly below the all-F count and strictly above 0. The baseline is the
      WHOLE F branch, so the subtrahend is the Caucasian / None slice of it, not a
      Caucasian-only baseline. The F node is asserted indeterminate too,
      as the propagation half, but F alone cannot fail: with SEX / RACE /
      SEVERITY, F reads indeterminate even while Caucasian is fully checked,
      because F's other RACE children (Other, Black, Asian) are unchecked.
  - anchor: Step 9
    expectation: >-
      GROK-16528 — after reordering columns in the picker (RACE first, then
      SEX), the card caption reads "RACE / SEX" and the top level of the tree
      now shows RACE values as root nodes.
  - anchor: Step 6-restore
    expectation: >-
      SR-02's coherent pre-save state — the card is reset to the SEX / RACE
      hierarchy and "Caucasian" alone is checked again, and the criterion is
      proven live before the save: df.filter.trueCount is exactly the
      Caucasian-female count 2823 recorded at Step 6, so the layout is captured
      while it is genuinely narrowing rather than in the criterion-free state a
      column reorder leaves behind.
  - anchor: Step 10
    expectation: >-
      Saving the layout through saveLayout + grok.dapi.layouts.save returns a
      layout carrying a non-empty id, so Step 11 has a real saved artefact to
      look up and re-apply rather than silently round-tripping nothing.
  - anchor: Step 11
    expectation: >-
      Layout round-trip (hierarchical half) — the panel is first brought back
      into a filtering state (a RACE root node checked in the reordered tree)
      and the narrowed df.filter.trueCount is recorded immediately before the
      close gesture: it is strictly positive and strictly below the full row
      count. Closing the panel really undoes it — df.filter.trueCount is back at
      the full row count, strictly above the recorded pre-close value, and no
      [name="viewer-Filters"] root remains in the DOM. After re-applying the
      saved layout, the hierarchical filter card is present, the caption reads
      "SEX / RACE", and df.filter.trueCount equals the recorded value from Step
      6.
  - anchor: Step 12
    expectation: >-
      The calculated column SEX_bool is created and its product-reported type
      is 'bool', and the table's own pre-existing boolean column CONTROL is
      still present — together the two boolean columns that
      showBoolCombinedFilter requires, asserted as a precondition so a missing
      combined boolean card in Step 14 cannot be blamed on the fixture.
  - anchor: Step 14
    expectation: >-
      Exactly one .d4-filter card with class d4-bool-combined-filter exists in
      the panel. The active-filter counter (.d4-filter-group-header
      .d4-filter-indicator) reads 0 (the combined card exists but nothing is
      toggled yet), and df.filter.trueCount is the full row count — the freshly
      auto-added card is not narrowing anything, which is what makes the Step 15
      drop attributable to the toggle.
  - anchor: Step 15
    expectation: >-
      After toggling the first boolean column's true-row in the combined card,
      df.filter.trueCount equals the number of rows whose first boolean column
      is true — the expected value is computed from that column's data, not
      hardcoded, and the table is checked to hold exactly two boolean columns
      that are neither all-true nor all-false. The active-filter counter reads 1
      once the panel's debounced counter refresh has run.
  - anchor: Step 16
    expectation: >-
      After toggling a second boolean column's row, df.filter.trueCount_OR
      equals the number of rows where either boolean column is true; after
      switching the card to AND mode with the same two criteria,
      df.filter.trueCount_AND equals the number of rows where both are true.
      Both expected values are computed from the two columns' data, and the
      intersection is checked to be strictly smaller than the union, so the AND
      mode is required to narrow — not merely "not widen".
  - anchor: Step 17
    expectation: >-
      GROK-16488 — after reopening the project, the combined boolean filter card
      is present in the Filter Panel. The reopened card is then put into a
      filtering state (df.filter.trueCount equals the first boolean column's
      true-row count) so that the removal has something to undo. After clicking
      the card's X icon (remove), the card count drops to zero,
      df.filter.trueCount returns to the full row count, and there is no console
      error (console-error delta across the removal is zero).
  - anchor: Step 18
    expectation: >-
      Any combined boolean card still on screen is first removed through the
      product — the card's own X icon — and the panel is asserted to hold zero
      such cards and zero bool-columns filters. The Add Filter > Combined
      Boolean leaf is then driven UNCONDITIONALLY (never guarded on "is one
      already there", which made the menu flow run on some executions only), and
      the count is asserted to go 0 -> 1: zero read immediately before the
      gesture, exactly one after it — exactly one, because a leaf that fired
      twice would leave two identical cards. Then the layout round-trip
      (combined boolean half) — the recorded pre-save count is a narrowed one,
      and closing the panel really undoes it (no combined boolean card remains
      and df.filter.trueCount is back at the full row count). After re-applying
      the saved layout, the combined boolean card is present again and
      df.filter.trueCount equals the value recorded before the layout was saved.
  - anchor: Step 19
    expectation: >-
      After Remove All the panel holds zero cards, zero
      .d4-bool-combined-filter cards, and df.filter.trueCount is the full row
      count — the perturbation is asserted before the recreation is. After
      closing the panel (its element leaves the DOM) and reopening it with the
      ribbon funnel icon, EXACTLY one .d4-bool-combined-filter card is present
      (not "at least one" — a rule firing twice would leave two), the table
      still carries more than one boolean column so the showBoolCombinedFilter
      precondition still holds, and df.filter.trueCount is the full row count
      because a freshly recreated card filters nothing.
---

# Filter Panel — Hierarchical and Combined Boolean Filters

## Setup

1. Open `System:DemoFiles/demog.csv` using the `openFile` helper. Wait for the table view to appear.
2. Prime the Filter Panel with `grok.shell.tv.getFiltersGroup()` (the `openTable`/`openFilterPanel`
   helper path). Wait for the Filter Panel to appear and for at least one filter card to render
   before proceeding, then hover the first card.
3. Record the total row count (expected: 5850 for demog).

## Scenarios

### Scenario 1: Hierarchical Filter — tree interaction and search

Steps:

1. Add a hierarchical filter through the panel's own hamburger menu: **Add Filter >
   Hierarchical**. Wait for exactly one hierarchical card to appear, showing the hierarchy tree.
   This is a real menu gesture, not a shortcut — adding the card through the filter group would
   reach the same card without the menu ever opening, and the menu is one of the flows this
   section is answerable for. Only the COLUMN configuration in the next steps stays programmatic,
   because its picker is a canvas-rendered column grid.
2. Record the current filtered row count (should equal the total row count — no criterion applied yet).
3. Apply a two-column hierarchy (SEX then RACE) by configuring the hierarchical filter card with
   columns SEX and RACE in that order, with all values enabled.
4. Wait for the card caption to update to reflect the configured columns.
5. Verify: the card caption reads "SEX / RACE". The filtered row count equals the total row count
   (no category checked yet).
6. Expand the "F" (female) node by clicking its expand arrow, and wait for its RACE children to
   appear. Then locate the "Caucasian" node **inside the "F" branch** — scoped by ancestor, since
   "Caucasian" also exists under "M" — and check its checkbox. Wait for the filtered row count to
   change.
   Verify: the filtered row count drops below the total row count and is strictly positive. Record
   this count as the Caucasian-female count. Verify also that the "Caucasian" node's own checkbox
   substitute renders the checked glyph U+F14A — under SEX / RACE, RACE is the terminal level, so
   this node is a leaf and a checked leaf must read checked. Read the glyph as a codepoint: the
   unchecked and indeterminate glyphs are both font-weight 100, so a font-weight read cannot make
   this distinction.
7. Bring "F" to fully checked and verify it renders the checked glyph U+F14A first — otherwise the
   indeterminate read below is a standing state rather than a transition. Then uncheck its
   "Caucasian" child and wait for the filtered row count to change.
   Verify: the "F" node's checkbox substitute renders the indeterminate glyph U+F146, and
   explicitly not the unchecked U+F0C8 (the two share font-weight 100, so an untouched node used to
   satisfy this check). The filtered row count is strictly above the single-descendant count and
   strictly below the all-F count.
8. GROK-19968 — Tree search. The search input is hidden until the search icon in the hierarchical
   filter card's header is clicked, so click that icon first, then type "Cau" in the search input.
   Confirm the typed term actually landed by reading the input's own value back — a keystroke that
   never reached the control would otherwise leave the assertions below reading an untouched tree.
   Wait for the visible node list to stabilize.
   Verify: the visible tree nodes shrink to the matches and their ancestors — "Caucasian" is still
   visible, the non-matching sibling "Asian" is not, and the visible node count is lower than
   before typing. The filtered row count is UNCHANGED from the value in Step 7 (which is a narrowed
   count, not the full table) — the tree search hides display nodes, it does not move the row-level
   filter. Clear the search input and verify the visible node set and the row count return to
   exactly what they were before typing.
8b. Three levels — the partial state has to travel further than one step. Reconfigure the card to
    the hierarchy SEX / RACE / SEVERITY and verify the caption names all three. Check the whole "F"
    branch and verify it reads as fully checked (U+F14A), then expand F and expand "Caucasian"
    under it, addressing both by caption path — the `[name="tree-expander-…"]` handle only names
    top-level nodes. Verify the SEVERITY children listed under F / Caucasian are EXACTLY the
    SEVERITY values the data holds for that branch — derive that set from the raw SEX / RACE /
    SEVERITY column values, the same way the counts are derived, and never read it back from the
    tree, which would accept whatever the tree happened to render. On demog today the set is
    Critical, High, Low, Medium, None. Each generated level also renders one caption-less row (the
    bold column header naming the level), which the tree reader reports separately rather than as an
    empty caption; verify there is exactly one and that it reads "SEVERITY", so a build that emits
    more of them — or that strips the caption off a real value node — is visible instead of being
    filtered away. "None" is among the derived values, so the grandchild this step unchecks exists.
    The load-bearing node is "Caucasian", and its state is asserted as a TRANSITION: read it
    BEFORE the grandchild uncheck and verify it is checked (U+F14A), then uncheck "None" and verify
    it is indeterminate (U+F146).
    Verify: the uncheck is local — the level still lists the same values, every SEVERITY sibling in
    the derived set still reads checked (walk the derived set, not the tree's own list, so a sibling
    that vanished fails as a missing node instead of shortening the loop) and "None"
    reads unchecked (U+F0C8); the "F" grandparent reads indeterminate too, which is the propagation
    half; and the filtered row count is the all-F count minus the rows that are female AND Caucasian
    AND SEVERITY "None" (3243 - 1588 = 1655 on demog), below the all-F count and above zero, so the
    tree state and the row set agree. The baseline here is the WHOLE "F" branch — not the
    Caucasian-only state of Step 6 — so the subtrahend is only the Caucasian / None slice of it.
    Both counts are computed from the raw SEX / RACE / SEVERITY column values rather than hardcoded,
    and never read back from the row filter, which would make the expectation circular; each derived
    number is checked to be strictly positive and strictly smaller than the all-F count before it is
    used, so a mis-derived zero cannot satisfy the assertion for free. Note that "F" alone would not
    fail here: under SEX / RACE / SEVERITY, F
    reads indeterminate even while Caucasian is fully checked, because F's other RACE children
    (Other, Black, Asian) are unchecked. Only the Caucasian transition discriminates.
9. GROK-16528 — Column reorder. Reconfigure the hierarchical filter card to place RACE first and
   SEX second (swap the column order). Wait for the card caption to update.
   Verify: the card caption reads "RACE / SEX". The top-level tree nodes are now RACE values
   (not SEX values).
6-restore. Before the layout can be saved the hierarchy has to be back in a state that is actually
    filtering. A column reorder CLEARS the criterion, so the state left by Step 9 would be captured
    empty, while Step 11 expects the restored layout to carry the caption "SEX / RACE" and the
    Caucasian-female count from Step 6. This is the ordering recorded as scope reduction SR-02 — the
    layout is saved in the coherent pre-reorder state and the Step 9 reorder is performed afterwards.
    Reconfigure the card back to the two-column SEX / RACE hierarchy, expand "F", and check
    "Caucasian" alone.
    Verify: the filtered row count is exactly the Caucasian-female count recorded in Step 6 (2823),
    so the criterion about to be saved is proven live rather than assumed.
10. Save the current layout using the **Save Layout** action. Record the current filtered row count
    before saving.
    Verify: the save returns a layout carrying a non-empty identifier — otherwise Step 11 would look
    up nothing and its round-trip would prove nothing.
11. The Step 9 reorder cleared the criterion, so first put the panel back into a filtering state:
    check the "Caucasian" root node of the reordered tree and record the resulting row count
    immediately before closing. Verify it is strictly positive and strictly below the total row
    count — otherwise the close has nothing to undo. Then perturb the live state: close the Filter
    Panel using its close button.
    Verify the perturbation actually happened: the filtered row count is back at the total row
    count (strictly above the pre-close value) and no Filter Panel root remains in the DOM.
    Re-apply the saved layout. Wait for the filter panel to rebuild and the filter cards to appear.
    Verify: the hierarchical filter card is present, the card caption reads "SEX / RACE" (the
    original column order from Step 5), and the filtered row count equals the Caucasian-female
    count recorded in Step 6.
    Delete the probe layout in cleanup, even on failure.

### Scenario 2: Combined Boolean Filter — AND/OR and project round-trip

Steps:

12. Open `System:DemoFiles/demog.csv` in a fresh view (or reset the current view's filters). To
    make the combined boolean filter appear, the table needs at least two boolean columns. Add a
    calculated boolean column named `SEX_bool` with the expression `${SEX} == "F"`.
    Wait for the new column to appear in the grid.
    Verify: the new column's product-reported type is `bool`, and the table's own pre-existing
    boolean column CONTROL is still present. Those two together are the precondition
    `showBoolCombinedFilter` needs, so a missing card in Step 14 cannot be blamed on the fixture.
13. Prime the Filter Panel with `grok.shell.tv.getFiltersGroup()` (the ribbon filter icon is the
    equivalent UI path; the panel is opened programmatically here). The panel should automatically
    add a combined boolean filter card when more than one boolean column exists.
14. Verify: exactly one combined boolean filter card is present in the panel. The active-filter
    counter shows "0" (or is hidden). Record the current filtered row count (should equal the full
    row count).
15. Toggle one boolean column's true-condition inside the combined card (true on, false off — a
    column with both values checked means "no constraint"). Wait for the filtered row count to
    change.
    Verify: the table has exactly two boolean columns; the filtered row count equals the number of
    rows whose first boolean column is true, computed from that column's data; the active-filter
    counter reads "1" (the counter is refreshed on a short debounce, so give it a moment).
    Record this count.
16. Toggle the second boolean column's row as well (add it to the filter). Record the count with
    both columns active in OR mode. Then switch the card to AND mode — the card body is
    canvas-rendered, so the mode is driven through the card's `mode` state rather than by clicking
    the AND/OR toggle. Wait for the filtered row count to change. Record the count in AND mode.
    Verify: the OR-mode count equals the number of rows where either boolean column is true and the
    AND-mode count equals the number of rows where both are true, both computed from the columns'
    data. Since the intersection is strictly smaller than the union on this table, the AND count
    must be strictly below the OR count — a mode switch that did nothing would fail here. Both
    counts are strictly positive and strictly below the full row count.
17. GROK-16488 — Save as project and reopen. First reset the combined boolean card to its no-op
    state (both the true and the false checkbox on, on both columns), so the reopened project starts
    at the full row count and the re-arming below is the only thing narrowing it. Record the current
    filtered row count before saving.
    Save the project using the ribbon Save button (use the `saveProjectViaUI` helper). Reopen the
    project. Wait for the Filter Panel to appear.
    Verify after reopen: the combined boolean filter card is present. Put it into a filtering state
    (toggle the first boolean column's true-condition) and verify the filtered row count drops to
    that column's true-row count, so the removal has something to undo. Click the card's remove icon
    (the X button on the card) to remove it. Wait for the card to disappear.
    Verify: the combined boolean card is gone. The filtered row count returns to the full row count.
    There is no new console error from the removal step.
    Delete the probe project in cleanup, even on failure.
18. Bring the combined boolean card back through the panel's own hamburger menu (**Add Filter >
    Combined Boolean**), the same real menu gesture Step 1 uses for the hierarchical card. The menu
    drive is unconditional: first remove any combined boolean card still on screen through the
    product — click the card's own X icon — and verify none remains (zero cards, zero bool-columns
    filters in the group). Then drive the leaf and verify the count went 0 -> 1: it was zero at the
    moment of the gesture and is exactly one afterwards, so the leaf demonstrably created the card
    rather than a pre-existing one satisfying the check. Exactly one, not "at least one" — a rule
    that fired twice would leave two. Guarding the drive on "is one already there" would leave the
    menu flow exercised on some executions only.
    Save a layout with the combined boolean filter active; the recorded row count must be a narrowed
    one. Then perturb the state (close the panel) and verify the perturbation actually happened —
    no combined boolean card remains and the row count is back at the full row count. Re-apply the
    saved layout and wait for the panel to rebuild.
    Verify: the combined boolean card is present and the filtered row count equals the value recorded
    before the layout was saved.
    Delete the probe layout in cleanup, even on failure.

19. Remove All, then reopen — the combined boolean card must come back on its own. Step 13/14 shows
    the card appearing the first time the panel is opened on a table with two boolean columns; this
    is the other half. Open the panel menu and click **Remove All**, and confirm the perturbation
    really happened: no cards at all remain, no combined boolean card remains, and every row is
    visible again. Then close the Filter Panel with its own close control, confirm it is gone, and
    reopen it with the ribbon funnel icon.
    Verify: exactly one combined boolean card is present — not "at least one", because a rule that
    fires twice would leave two identical cards; the table still has more than one boolean column,
    so a missing card would be the product's doing and not the fixture's; and the freshly recreated
    card is not filtering anything, so the full row count is back.

## Automation notes

Tree node addressing by caption path, the caption-less column-header row each generated level opens
with, the tri-state glyph codepoints, the per-card hover-only icons, the 50 ms group-counter debounce,
the console-error rule and the panel-close control — which is host-dependent, so BOTH candidate names
must be looked up inside the filter viewer's own `.panel-titlebar` — are in
`.claude/skills/grok-browser/references/viewers/filters.md`.

- **Canvas-driven preconditions.** The hierarchical filter's column picker
  (`.d4-column-selector-backdrop`) and the combined boolean card's row bodies are canvas-rendered,
  so the column configuration and the boolean row toggling (including the AND/OR mode) are driven
  through `grok_GridFilterBase_ApplyState`, followed by
  `grok.shell.tv.dataFrame.rows.requestFilter()` to trigger the BitSet recomputation. The gesture
  under test in both scenarios is the resulting card behaviour, not the picker or the toggle. Same
  reasoning for `SEX_bool`, created with `columns.addNewCalculated` rather than through the
  canvas-heavy formula dialog.
- **Setup:** the panel is primed by `openTable({withFilterPanel: true})` →
  `grok.shell.tv.getFiltersGroup()`, not by the ribbon `[name="icon-filter"]` click. Wait for
  `[name="viewer-Filters"]` and at least one `.d4-filter` card — the panel runs `look.auto()`
  asynchronously — then hover the first card.
- **Adding a card through the panel's own hamburger** (Scenario 1 step 1, and Step 18 for the
  combined boolean leaf) goes through the shared `drivePanelMenuLeaf` helper: it opens the panel
  hamburger, hovers the **Add Filter** group and clicks the leaf, resolving both labels inside the
  panel's own popup — the view's top menu carries identically-labelled items and precedes it in
  document order. Groups open on HOVER, and the rows carry no per-level `[name="div-…"]` handles, so
  `driveTopMenuLeaf` does not apply here. After the leaf, poll until the card is in the group rather
  than waiting a fixed time.
- **Step 8b derives its expectations from the raw SEX / RACE / SEVERITY columns**, never from the
  tree and never from `df.filter` — reading the sibling list out of the tree and asserting against
  the same tree is circular, and it was what let the caption-less header row through. Each derived
  number is checked to be strictly positive and strictly smaller than the all-F count before it is
  used. Expected on demog today: Critical, High, Low, Medium, None (10 / 319 / 552 / 354 / 1588),
  plus exactly one caption-less header row reading "SEVERITY".
- **Step 8 counts only tree nodes that carry `.d4-hierarchical-filter-caption-value` and are laid
  out** (`offsetParent !== null`) — the bold column-header items are never hidden by the search, so
  counting them dilutes the signal. The search round-trip is synchronous and display-based and needs
  no stabilization waits (measured 2223 / 2206 / 2219 ms across three attempts); do not relax the
  captions-equality assert on the restored node set to chase a hypothetical flake.
- **Step 11 close:** the close is a mandatory click (throw if the control is missing) and must drop
  `[name="viewer-Filters"]` from the DOM. Re-apply with `loadLayout` and wait on the layout-apply
  polling barrier.
- **Signal discipline.** `df.filter.trueCount` is the primary observable on both halves and is
  updated synchronously by `requestFilter()`. Canvas content is never read pixel-wise. The
  hierarchical card additionally uses DOM signals (visible node count, checkbox glyph, card
  caption); the combined boolean card uses DOM card presence only.
- **A boolean column counts as applied only when its true/false checkboxes differ** — both checked
  means "no constraint", both unchecked contributes nothing. The expected row counts are derived
  from the columns' own values rather than hardcoded.
- **Step 17 project round-trip:** `saveProjectViaUI`, then wait for the panel keyed to
  `[name="viewer-Filters"]` DOM presence (GROK-19152). The console-error claim window here is the
  X-click plus 800 ms — under ~1% of the ~85 s test — which is what makes the signature-only
  whitelist of the project-save emitters safe. Do not stretch that window to span a whole phase.
