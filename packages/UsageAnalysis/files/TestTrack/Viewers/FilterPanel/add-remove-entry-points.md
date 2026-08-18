---
feature: filters
realizes_atlas:
  - filters.cp.add-remove-entry-points
  - filters.cp.multi-value-filter
realizes:
  - viewers.filters
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19516
    status: fixed
  - id: GROK-18765
    status: fixed
  - id: GROK-12955
    status: fixed
realized_as:
  - add-remove-entry-points-spec.ts
scope_reductions:
  - id: SR-01
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      SETUP substitution, not substitution of the gesture under test.
      Narrowing a card to its largest category is done over the JS API
      (`fg.updateOrAdd({type: CATEGORICAL, column, selected: [top]})`)
      rather than by clicking the category row on the card canvas. The
      Automation notes already authorise this for Scenario 1 Steps 6-7;
      this entry extends the same declaration to the two later steps that
      use the identical technique for the identical purpose — Scenario 5
      Step 1 (RACE, to give the panel a criterion worth restoring across
      close/reopen) and Scenario 6 Step 3 (SEX, to give "Filter to
      Column..." an active filter to write out). Neither step asserts
      anything about how the criterion was applied: what they assert is
      the close/reopen restore and the boolean column the menu command
      writes, and both of those remain fully UI-driven. The category-row
      click itself is exercised end-to-end by the sibling
      panel-core-ladder scenario.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: "Setup — the zero-card baseline"
    expectation: >-
      After "Remove All" on the auto-populated panel, the panel carries exactly
      0 .d4-filter cards, df.filter.trueCount is the full row count (5850), and
      the group-header active-filter counter is hidden or absent. The card count
      of 0 is what makes every "the count is now exactly N" claim in Scenario 1
      a statement about cards this run added.
  - anchor: "Scenario 1 Step 1 — add path (a) header drag"
    expectation: >-
      After dragging the RACE grid column header onto the Filter Panel with the
      real mouse and releasing over the panel's drop zone, a RACE card exists
      and its caption is the FIRST entry of the ordered .d4-filter-column-name
      list, the card count is exactly 1 (it was 0 before the drag),
      df.filter.trueCount is unchanged at 5850, and the group-header
      active-filter counter is hidden or reads 0. The card must actually appear
      — an unchanged card count is the failure, not a tolerated outcome.
  - anchor: "Scenario 1 Step 2 — GROK-19516 pinned header drag"
    expectation: >-
      After pinning AGE and dragging the PINNED AGE column header onto the
      Filter Panel, an AGE card appears whose caption is the FIRST entry of the
      ordered caption list with RACE second, the card count is exactly 2,
      df.filter.trueCount is unchanged at 5850, and the counter is hidden or 0.
      GROK-19516 was the pinned source producing NO card at all, so the negative
      half is load-bearing: the count must have grown by exactly 1 and the AGE
      caption must be present. The pin itself is proven on both sides of the
      gesture, because a silently-unpinned column would make the whole guard
      vacuous: AGE's grid index is read together with the grid's frozen-column
      count (a column is pinned exactly when idx < frozenColumns), both readings
      must be numbers, AGE must NOT be pinned before the context-menu gesture,
      and it must be pinned after it.
  - anchor: "Scenario 1 Step 3 — add path (b) panel header picker"
    expectation: >-
      After using the panel header's column picker to add SEX, the SEX caption
      is FIRST in the ordered list (appendToTop), the card count is 3,
      df.filter.trueCount is unchanged at 5850, and the counter is hidden or 0.
  - anchor: "Scenario 1 Step 4 — add path (c) column properties Add filter"
    expectation: >-
      After clicking the "Add filter" action link in the Filter section of the
      DIS_POP column properties panel, a DIS_POP card appears FIRST in the
      ordered caption list, the card count is 4, df.filter.trueCount is
      unchanged at 5850, and the counter is hidden or 0. This path is a plain
      click, not a drag.
  - anchor: "Scenario 1 Step 5 — add path (d) context menu Expression"
    expectation: >-
      After adding an Expression filter from the panel context menu (Add Filter
      > Expression), a card carrying the Expression filter caption appears FIRST
      in the ordered list, that card really is an expression filter (its body
      carries a .d4-expression-filter editor, so a card that merely borrowed the
      caption fails), the card count is 5, df.filter.trueCount is unchanged at
      5850, and the counter is hidden or 0.
  - anchor: "Scenario 1 Step 6 — RACE card starts filtering"
    expectation: >-
      After narrowing the RACE card to a single category, df.filter.trueCount is
      strictly below 5850 AND strictly above 0, the group-header active-filter
      counter reads 1, and the card count is still 5. Record the reduced value
      as trueCount_race; a null/undefined reading must fail before it is
      compared.
  - anchor: "Scenario 1 Step 7 — GROK-18765 disable of a column-properties card"
    expectation: >-
      With the DIS_POP card (created in Step 4 through the column properties
      panel) narrowed so that df.filter.trueCount is strictly below
      trueCount_race and the counter reads 2, unchecking that card's own
      checkbox returns df.filter.trueCount to EXACTLY trueCount_race and the
      counter to 1, the card gains the d4-filter-disabled class, its filter
      reports isFiltering false, and the card count stays 5. Re-checking it
      restores both the lower trueCount and a counter of 2. GROK-18765 was the
      disable having no effect for a card created from this entry point, so the
      return of the row count is the half that carries the guard — an
      entry-point-agnostic disable check does not guard this bug.
  - anchor: "Scenario 2 Step 1 — GROK-12955 reorder drag"
    expectation: >-
      After dragging the bottom-most card by its caption label and releasing on
      the drop strip above the first card, the ordered caption list equals the
      pre-drag list with that one caption moved to position 0 and every other
      caption in its original relative order; the card count is unchanged;
      df.filter.trueCount is unchanged; the page still answers a round-trip
      evaluate after the drop; and the set of browser console error TEXTS
      contains no entry that was absent before the drag (compare texts, not
      counts; exclude the ambient "Permissions policy violation:
      compute-pressure" message by name). GROK-12955 was a freeze plus repeated
      "Cannot read properties of null (reading 'cy')" noise, so both the
      responsiveness check and the console-text delta are load-bearing — an
      order-only check passes while the freeze is present.
  - anchor: "Scenario 2 Step 2 — reorder revert"
    expectation: >-
      After dragging the same card back to the position it held before Step 1,
      the ordered caption list is IDENTICAL to the list recorded at the start of
      Step 1, the card count and df.filter.trueCount are unchanged, and again no
      console error text absent before the drags has appeared.
  - anchor: "Scenario 3 Step 1 — remove via card X icon"
    expectation: >-
      After clicking the X icon on the Expression filter card, the ordered
      caption list no longer contains the Expression filter caption, the card
      count dropped by exactly 1, and df.filter.trueCount is unchanged (that
      card was not filtering).
  - anchor: "Scenario 3 Step 2 — remove others"
    expectation: >-
      After hovering the group-header active-filter counter and clicking "Remove
      others" (removeNonActiveFilters), every non-filtering card's caption is
      gone from the ordered list while the RACE and DIS_POP captions — the two
      that are actively filtering — are both still present, and
      df.filter.trueCount is unchanged (removing non-filtering cards cannot move
      the filtered count).
  - anchor: "Scenario 3 Step 3 — Remove All"
    expectation: >-
      After clicking the panel context menu's Remove All, zero
      .d4-filter-column-name elements remain, df.filter.trueCount is the full
      row count (5850), and the active-filter counter is hidden or reads 0.
      Remove All fires immediately — no confirmation dialog appears and no OK
      click may follow it.
  - anchor: "Scenario 4 Step 1 — the cards exist before the columns are hidden"
    expectation: >-
      After emptying the panel with "Remove All", closing it and reopening it —
      the sequence that makes the panel build a fresh card set from the current
      column visibility rather than restore the set it saved — the ordered
      caption list contains BOTH RACE and AGE and df.filter.trueCount is the full
      row count (5850). "Remove All" must really leave zero cards, since a
      non-empty saved set turns the reopen into a restore. This is the
      perturbation proof for the whole scenario: without it, their absence in
      Step 3 would prove nothing, and it is read through the SAME path as
      Steps 3 and 4 so all three measure one behaviour.
  - anchor: "Scenario 4 Step 2 — the hide gesture really takes"
    expectation: >-
      For each of RACE and AGE in turn: the column reads visible in the grid
      BEFORE the gesture (the dialog's per-row check cell is a toggle, not a
      set, so a wrong pre-state would show the column instead of hiding it), the
      dialog's search box is narrowed to that column so it is row 0, its check
      cell is clicked in the dialog's canvas ColumnGrid, and
      grok.shell.tv.grid.columns.byName(c).visible then reads false while the
      dialog is still open — visibility applies live, so the click landing is
      what is asserted, naming the column, its visible state before and after,
      and the computed click point. The dialog is then closed and the column is
      still hidden. A gesture that did not land must fail here rather than
      downstream, and no JS visibility write may stand in for it.
  - anchor: "Scenario 4 Step 3 — hidden columns have no cards"
    expectation: >-
      After hiding RACE and AGE through the Order or Hide Columns dialog and
      making the panel rebuild its card set from column visibility ("Remove All",
      close, reopen — a bare close/reopen restores the saved card set instead and
      would keep the cards), neither caption is in the ordered caption list, the
      list is NOT empty and still contains SEX (a column that was not hidden),
      and df.filter.trueCount is still the full row count (5850) — hiding a
      column is not a filter. Both halves are read in the same step: the
      negative alone is also satisfied by a panel that built no cards at all.
  - anchor: "Scenario 4 Step 4 — showing the columns again restores their cards"
    expectation: >-
      For each of RACE and AGE in turn: the column reads hidden in the grid
      before the gesture, the dialog — reopened from SEX, since a hidden column
      has no header left to right-click — is narrowed to that column and its
      check cell clicked, and the column then reads visible while the dialog is
      still open, the failure message naming the column, its visible state
      before and after, and the computed click point. The panel is then rebuilt
      through the SAME path as Step 3 ("Remove All", close, reopen), and both
      captions are back in the ordered list with df.filter.trueCount still 5850,
      the failure message printing the caption list it actually got. Using the
      same path both times is what makes this step falsifiable: Step 3 measured
      that this path yields no card for a hidden column, and Step 3 left a panel
      with no RACE or AGE card for a restore to bring back, so a card appearing
      here can only come from the restored visibility.
  - anchor: "Scenario 5 Step 1 — three distinguishable card states"
    expectation: >-
      Narrowing RACE to a single category puts df.filter.trueCount strictly
      below 5850 and strictly above 0; that value is recorded as the RACE-only
      count. Unchecking the AGE card's own checkbox leaves it unchecked, and
      removing the SEX card with its X icon takes the SEX caption out of the
      ordered list. Neither of those two was filtering, so df.filter.trueCount
      is still EXACTLY the RACE-only count afterwards and the group-header
      active-filter counter reads 1. Three different states are needed because a
      single-state check cannot tell a panel that restored its state from one
      that merely rebuilt itself from scratch.
  - anchor: "Scenario 5 Step 2 — closing the panel releases its filtering"
    expectation: >-
      After the title-bar close, [name="viewer-Filters"] is gone from the DOM
      and df.filter.trueCount is the full row count (5850). This is the
      perturbation proof — without it the restore assertion in Step 3 could not
      fail.
  - anchor: "Scenario 5 Step 3 — reopening restores per-card state"
    expectation: >-
      After reopening through the ribbon funnel icon and waiting for the cards
      to rebuild, the RACE card is present and df.filter.trueCount equals the
      RACE-only value recorded in Step 1 exactly; the AGE card is present, its
      own enable/disable checkbox is FOUND (a missing checkbox fails the check
      rather than passing it) and is still unchecked; the SEX card, which was
      removed rather than disabled, is NOT recreated; and the header
      active-filter counter reads 1. The reopen here is deliberately a bare one,
      with no "Remove All" before it: the saved card set is non-empty, so the
      panel restores it instead of rebuilding, and that restore is the subject of
      this scenario.
  - anchor: "Scenario 6 Step 1 — Select Columns... drives the card set both ways"
    expectation: >-
      Starting from a panel emptied with "Remove All" — asserted at zero cards, because
      otherwise the card set produced below would not be attributable to the picker — and with
      the row count df.filter.trueCount recorded as this scenario's baseline (it must be a
      positive number; see the Automation note on why it is recorded rather than assumed to be
      5850), the hamburger's "Select Columns..." leaf opens the
      .d4-dialog[name="dialog-Select-columns..."] picker and its "N checked" label reads
      exactly 0, matching the empty panel (a reading of -1 means the label was not found at
      all and fails the check rather than passing it). Clicking the "All" link
      ([name="label-All"]) takes the label to exactly the table's column count, and OK closes
      the dialog and leaves the panel with exactly that many cards whose caption SET equals
      the table's column-name set, with df.filter.trueCount still at the recorded baseline —
      the committed cards are all unconfigured, so choosing which columns get cards is not
      filtering. Reopening the picker then reads back that same
      count, which is what proves the dialog reflects the live card set rather than a fixed
      default; clicking "None" takes the label to 0, OK empties the panel again, and the row
      count is still the baseline — none of those cards was filtering. Both
      directions are asserted in one step: an "All" that adds cards proves nothing on its own
      if a "None" that removes them is never shown to work.
  - anchor: "Scenario 6 Step 2 — Add Filter > Multi Value... really filters a multi-value column"
    expectation: >-
      From the empty panel Step 1 left, and with the filters group carrying ZERO filters whose
      filterType is 'multi-value' (asserted first, so the new one is attributable), the step opens
      a FRESH demog table view and closes the previous one. It has to: by this point in the run the
      original demog dataFrame still has "RACE: Caucasian" applied through a leaked filter
      subscriber that no panel gesture releases (see the Automation note), and this step's whole-table
      arithmetic cannot be measured on a table something else is already filtering. The fresh view is
      asserted clean before anything else happens — exactly one [name="viewer-Filters"] panel in the
      DOM, that panel emptied to zero cards, dataFrame.rows.filters EMPTY, and df.filter.trueCount
      the full 5850. A declared FIXTURE string column `MVF_TOKENS` is then added to it: row i takes the
      `;`-joined token list at index i mod 10 of a fixed 10-entry pattern, so the data is
      reproducible run to run. Before it is used, the distribution is asserted from the column's
      OWN values — 5850 rows covered, ZERO rows without a token, and for the pair
      `red` / `blue`: each token's row count strictly between 0 and 5850, their intersection
      strictly above 0 AND strictly below the smaller single count, their union strictly above
      the larger single count and strictly below 5850. Degenerate data that would make the
      filtering assertions unfalsifiable fails here, before any gesture.
      The hamburger's Add Filter > "Multi Value..." leaf then opens
      .d4-dialog[name="dialog-Select-column-and-separator-for-multi-value-filter"]. demog carries
      no column tagged with a MultiValueSeparator, so the categorical FALLBACK is what the dialog
      opens on: its column selector shows a non-empty column name that is NOT `MVF_TOKENS` (so the
      pick below is a real change), no "No categorical columns found in the table" balloon is on
      screen, and the dialog refuses to commit until a separator is given — its OK button carries
      the `disabled` class and its Separator editor carries `d4-invalid`. The column selector is
      then driven to `MVF_TOKENS` and the selector is asserted to READ `MVF_TOKENS` before the
      commit: a dialog that silently kept its fallback column would build a card over patient ids
      and every count below would be measuring the wrong thing. Typing `;` enables OK; the commit
      closes the dialog and leaves the panel with exactly 1 card captioned `MVF_TOKENS`, the
      filters group holding exactly one filter with filterType 'multi-value' on that column — a
      card that merely borrowed the caption fails — and df.filter.trueCount still at the full 5850,
      because a multi-value filter with nothing selected excludes no further row. The fresh card
      reports include [] and mode AND.
      The card's values are then ticked on its canvas body, and every expected row count is
      derived by splitting the column's own cells on `;` and testing token membership — never
      from df.filter and never from the card's own state. Each row index 0..3 is ticked and
      unticked once to build the row -> token map; that map must be exactly the set of distinct
      tokens the column carries, so a card rendering the wrong values, or fewer rows than tokens,
      fails there. Then: ticking `red`'s row alone puts df.filter.trueCount at exactly the derived
      `red` count and that count is strictly between 0 and 5850; ticking `blue`'s row ADDS to the
      selection rather than replacing it and, in the default AND mode, puts the count at exactly
      the derived intersection, strictly below the single-token count; clicking the card header's
      AND/OR control flips its label to OR and the filter's own mode to OR, leaves the ticked
      tokens alone, and puts the count at exactly the derived union, strictly ABOVE the AND count
      — so a mode switch that did nothing fails; clicking it again restores the label AND and the
      intersection count exactly. Unticking `blue` returns the count to `red`'s own count, and
      unticking `red` leaves include [] with all 5850 demog rows back.
      The fixture column and the multi-value filter are removed in a `finally`, and the table's
      column list is asserted back to what it was.
      This step realizes the atlas critical path `filters.cp.multi-value-filter`, promoted
      from rev 6's declared coverage gap by operator decision on 2026-08-18 (atlas rev 9),
      and it closes the refdoc flow `pcmdAddFilterMultiValue`.
  - anchor: "Scenario 6 Step 3 — Filter to Column... writes the filter into a boolean column"
    expectation: >-
      With the SEX card narrowed to its largest category so df.filter.trueCount is strictly
      between 0 and 5850 (both bounds asserted — an all-false or all-true column would satisfy
      the equality below otherwise), the panel hamburger's "Filter to Column..." leaf opens
      .d4-dialog[name="dialog-Filter-to-Column"]. Its Name editor arrives pre-filled from the
      active filter; that proposal must be non-empty, and the table must NOT already carry a
      column of that name, or the column appearing afterwards would prove nothing. OK closes
      the dialog and a column of that name appears on the table, its type is 'bool', and
      counting the TRUE values of that column itself — not re-reading the filter — gives
      exactly the row count the panel was showing when the command ran. A count of -1 (column
      missing or not boolean) fails before that comparison is made. df.filter.trueCount is
      unchanged, since saving the filter does not alter it. The column is removed in a
      `finally`, and the table's column list is asserted back to what it was.
---

# Filter Panel — Add, Reorder, and Remove Entry Points

## Setup

1. Open the demog golden dataset (`System:DemoFiles/demog.csv`).
2. Open the Filter Panel: click the funnel icon on the ribbon toolbar above the
   grid.
3. The panel does not open empty — it automatically creates a card for every column
   it considers suitable. Bring it back to an empty state before any checks start:
   open the Filter Panel menu and click "Remove All". Verify the panel now shows
   zero filter cards, no rows are excluded from the dataset (all 5850 rows remain
   visible in the grid), and the header active-filter counter is hidden or reads 0.
   This empty panel is the starting point for the scenarios below — the zero is what
   makes every "the panel now shows exactly N cards" claim below a statement about
   cards these steps added.

## Scenarios

### Scenario 1: All four ADD paths insert at the top (GROK-19516, GROK-18765)

Steps:

1. Drag the RACE column header from the grid onto the Filter Panel: press the
   mouse button on the header, move the pointer over the panel until the panel
   shows its "Add filter" drop hint, and release the button over that hint.
   Verify: a RACE card is now the first card in the panel, the panel shows exactly
   1 card, no rows are excluded from the dataset (a freshly added filter is not
   yet filtering), and the header active-filter counter is hidden or reads 0.
2. Confirm first that AGE is NOT already pinned — a column counts as pinned when its
   grid index falls inside the grid's frozen block — then right-click the AGE column
   header and choose "Pin Column" (a leaf of the menu's "Pin" group) so AGE moves into
   the pinned block at the left of the grid. Verify AGE really is pinned now: this drag
   is only a pinned-source drag if the pin took, and a silently unpinned column would
   make the whole GROK-19516 guard vacuous. Then drag the pinned AGE header onto the
   Filter Panel the same way. Verify: the AGE card is now the first card, RACE is
   second, the panel shows 2 cards, no rows are excluded, and the counter is still
   hidden or 0. This is the GROK-19516 guard — dragging from a pinned header must
   produce a card exactly as an unpinned one does, so a panel that still shows
   only the RACE card is the failure.
3. Use the panel header's column picker button to add the SEX column. Verify: the
   SEX card is now the first card, the panel shows 3 cards, no rows are excluded,
   and the counter is hidden or 0.
4. With the DIS_POP column selected so that the column properties panel is showing
   that column, open its Filter section and click the "Add filter" action link there
   — that click is the add path under check. Verify: a DIS_POP card is now the first
   card, the panel shows 4 cards, no rows are excluded, and the counter is hidden or
   0.
5. Open the Filter Panel menu from the hamburger icon in the panel's title bar and
   choose Add Filter > Expression. Verify: an Expression filter card is
   now the first card and it really is an expression filter — its body carries the
   expression-filter editor, so a card that merely borrowed the caption fails — the
   panel shows 5 cards, no rows are excluded, and the counter is hidden or 0.
6. Narrow the RACE card to a single category so it starts filtering. Verify: the
   number of rows visible in the grid is now below 5850 and above zero, the header
   active-filter counter reads 1, and the panel still shows 5 cards. Note the
   reduced row count.
7. Narrow the DIS_POP card — the one added from the column properties panel in
   Step 4 — so that it starts filtering as well. Verify: the visible row count
   drops below the count noted in Step 6 and the counter reads 2; note this second
   count. Then uncheck the DIS_POP card's own checkbox. Verify: the visible row
   count returns to exactly the Step 6 count, the counter reads 1, the DIS_POP
   card is shown as disabled and reports it is no longer filtering, and the panel
   still shows 5 cards. Re-check that checkbox and verify the second count and a
   counter of 2 come back. This is the GROK-18765 guard — the defect was the
   disable having no effect for a card created from the column properties panel,
   so the return of the row count is the half that carries it.

### Scenario 2: Reorder by dragging a card, and GROK-12955

Steps:

1. Record the ordered list of filter card captions from top to bottom, and the
   browser console errors present right now. Then drag the bottom-most card by its
   caption label and release it on the drop strip above the first card. Verify: the
   caption list is the recorded list with that one caption moved to the front and
   every other caption in its original relative order; the card count is unchanged;
   the number of rows visible in the grid is unchanged (a pure reorder changes no
   filter criteria); the page still responds to interaction; and no console error
   text has appeared that was not already present before the drag. GROK-12955 was a
   freeze accompanied by repeated null-dereference console noise, so the
   responsiveness check and the console-error comparison carry the guard — an order
   check alone passes while the freeze is present.
2. Drag the same card back to the position it held before Step 1, the same way.
   Verify: the caption list is identical to the list recorded at the start of Step
   1, the card count and the visible row count are unchanged, and again no console
   error text has appeared that was absent before the drags.

### Scenario 3: Three REMOVE paths

Steps:

1. Hover the Expression filter card to reveal its icons and click its X icon.
   Verify: the Expression caption is gone from the panel, the total card count
   dropped by exactly 1, and the number of rows visible in the grid is unchanged —
   that card was not filtering.
2. Hover the header active-filter counter until its context-action menu appears and
   click "Remove others". Verify: the cards that were not filtering are gone from
   the panel; the RACE and DIS_POP cards, the two that are actively filtering, both
   remain; and the number of rows visible in the grid is unchanged (removing
   non-filtering cards cannot change the filtered row count).
3. Open the Filter Panel context menu and click "Remove All". Verify: zero filter
   cards remain in the panel, all 5850 rows are visible in the grid again, and the
   active-filter counter is hidden or reads 0. Removing all filters happens
   immediately — there is no confirmation dialog to accept.

### Scenario 4: A hidden column has no filter card, and gets one back when shown again

Steps:

The panel treats "which cards do I show" in two different ways, and only one of them looks at
column visibility. When it has a saved set of cards — which it keeps whenever it is closed with
cards in it — reopening brings that exact set back, whatever the columns are doing now (that is
Scenario 5). Only when it has nothing saved does it work the set out from the table, and only then
does it skip the columns that are hidden in the grid. So every check below empties the panel with
"Remove All" first, and only then closes and reopens it; a bare close/reopen would restore the old
cards and measure nothing about visibility.

1. The panel fills itself with a card per suitable column. Bring it back to that state — click
   "Remove All", verify the panel really is empty, then close the panel and reopen it — and record
   the ordered caption list. Verify RACE and AGE both have a card right now, and that no rows are
   excluded from the dataset. Their cards have to exist before they are hidden, otherwise their
   absence afterwards would prove nothing; and reading them through the same close/reopen sequence
   the later steps use is what makes all three steps measure one behaviour.
2. Hide RACE, and then AGE, through the grid's own dialog: right-click a column header, choose
   "Order or Hide Columns...", type the column's name into the dialog's search box so the list
   narrows to that single row, and click the check box in that row of the dialog's column list.
   Verify before each click that the column is still visible — the row's check box toggles, so
   clicking it from the wrong state would show the column instead of hiding it — and verify right
   after the click, with the dialog still open, that the column really is no longer visible in the
   grid; visibility changes as soon as the row is unchecked. If it did not change, the click did
   not land and the rest of the scenario proves nothing. Then close the dialog. Do NOT use the
   single check box next to the search box: it acts on the whole filtered list and does not
   refresh when the search text changes, so a second use of it undoes the first. Once a column is
   hidden it has no header left to right-click, so the dialog is reopened from a column that is
   still visible.
3. Make the panel work out its cards again: click "Remove All", verify the panel is empty, then
   close and reopen it. Verify: neither RACE nor AGE has a card any more; the panel is NOT empty
   and SEX — a column that was not hidden — still has its card; and no rows are excluded from the
   dataset, because hiding a column is not a filter.
4. Show RACE and AGE again the same way — reopen the dialog from SEX, narrow the search to the
   column, and click that row's check box. Verify before each click that the column is still
   hidden, and after it that it is visible in the grid again. Then empty, close and reopen the
   panel exactly as in Step 3. Verify both cards are back and still no rows are excluded. Doing it
   the same way both times is what gives this step something to fail on: Step 3 showed this
   sequence produces no card for a hidden column, and it left a panel with no RACE or AGE card to
   bring back, so a card here can only be the result of the column being visible again.

### Scenario 5: Closing and reopening the panel keeps each card's state

Steps:

1. Put three cards into three distinguishable states: narrow the RACE card to a single category so
   it filters, and note the row count — verify it is below 5850 and above zero, so the card is
   really filtering; uncheck the AGE card's own checkbox so it is disabled but
   still present; and remove the SEX card with its X icon. Verify the row count is still the
   RACE-only value (neither disabling AGE nor removing SEX was filtering anything), that SEX is
   gone from the caption list, and that the header active-filter counter reads 1.
2. Close the Filter Panel with its title-bar close button. Verify all rows are visible again — the
   panel's filtering is released while it is closed. This has to be confirmed, not assumed: if
   nothing changed here, "the state came back" in the next step would be unfalsifiable.
3. Reopen the Filter Panel with the ribbon funnel icon and wait for its cards to come back. Do not
   empty the panel first, unlike Scenario 4 — the panel kept a saved set of cards here, and
   restoring that set is exactly what this scenario is about.

   Verify: the RACE card is back and the row count is exactly the value noted in Step 1; the AGE
   card is back, carries its own enable/disable checkbox, and that checkbox is still unchecked — a
   card that came back enabled means the disabled state was not kept, and a card with no checkbox
   at all is a failure of the check rather than a pass; the SEX card is NOT recreated, because it
   was removed rather than disabled; and the header counter reads 1.

### Scenario 6: Three more Filter Panel menu commands

Steps:

1. Empty the panel with "Remove All" and confirm it really is empty — the card set built below
   only means something against a zero. Note the number of rows the grid is showing right now
   and check it is above zero; that number, not 5850, is what the checks below compare against
   (see the Automation notes — Remove All does not always release the previous scenario's
   filtering). Then open the Filter Panel menu and click
   "Select Columns...". A picker dialog opens listing the table's columns with a check box each,
   a search box, and a running "N checked" label. Verify the label reads 0, matching the empty
   panel. Click the "All" link. Verify the label now counts every column in the table. Click OK.
   Verify: the panel now shows one card per column and nothing else — the set of card captions
   is exactly the set of column names — and the number of rows shown is unchanged, because the
   new cards are all unconfigured and choosing which columns get cards is not filtering. Now
   reopen the picker. Verify it opens
   showing that same count checked; that is what shows the dialog reads the live panel rather
   than a fixed default. Click the "None" link, verify the label goes to 0, and click OK.
   Verify the panel is empty again and the row count is still unchanged. Both directions are
   checked here on purpose: an "All" that
   fills the panel proves nothing on its own unless a "None" that empties it is also shown to
   work.
2. From that empty panel, first confirm no multi-value filter is present, so a new one can only
   come from this gesture. Then open demog again in a fresh table view and close the old one: the
   checks below count rows across the whole table, and by this point the original table still has
   Scenario 5's RACE criterion applied by something the panel can no longer reach (see the
   Automation notes), which would skew every count. Confirm the fresh view is genuinely clean —
   only one Filter Panel on screen, emptied to no cards, nothing at all filtering the table, and
   all 5850 rows visible — before going on.

   demog has no genuinely multi-value column, so give it one: add a
   string column `MVF_TOKENS` whose every row holds one to three `;`-separated tokens drawn from
   a four-word vocabulary, the value for each row worked out from its row number so the same run
   always produces the same data. Check the data before using it: all 5850 rows have at least one
   token, none is blank, and for the two tokens the checks below use (`red` and `blue`) each one
   picks out some but not all rows, the rows carrying both are fewer than the rows carrying
   either one alone, and the rows carrying at least one are more than either alone but still not
   all of them. Without that, the counts below could not tell a working filter from a broken one.

   Then open the Filter Panel menu and choose Add Filter > "Multi Value...". A dialog opens
   asking for a column and a value separator. The demog table has no column marked with a
   multi-value separator, so the dialog falls back to the first categorical column: verify it has
   resolved some column anyway (its column selector is not blank), that the column it landed on
   is NOT `MVF_TOKENS` — otherwise choosing it below would change nothing — that the
   "No categorical columns found in the table" message did NOT appear, and that the dialog
   refuses to commit while the separator is empty (OK disabled, separator field flagged invalid).
   Now open the dialog's column selector and choose `MVF_TOKENS`, and verify the selector really
   reads `MVF_TOKENS` before going any further: a dialog that quietly kept its own column would
   give a card full of patient ids and every count afterwards would be about the wrong data.
   Type `;` as the separator, verify OK becomes enabled, and click it. Verify: the panel now
   shows exactly one card, captioned `MVF_TOKENS`, that card really is a multi-value filter
   rather than a card that merely took the same column, the row count noted in Step 1 is
   unchanged because nothing is selected yet, and the card reports nothing included and AND mode.

   Now use the card. Its values are drawn on a small canvas grid, one row per token, and clicking
   a row's check box adds or removes that token. Every number checked below is worked out from the
   column's own cells — split each on `;` and see which tokens it holds — never from the filter
   itself and never from the card. First tick and untick each of the four rows in turn to learn
   which token each row stands for, and verify those four are exactly the four tokens the column
   contains. Then: tick `red`'s row alone and verify the grid shows exactly the rows whose own
   value contains `red`, and that this is more than none and fewer than all 5850. Tick `blue`'s
   row as well and verify it was added to the selection rather than replacing `red`, and that the
   grid now shows exactly the rows carrying BOTH tokens — fewer than with `red` alone, because the
   card is in AND mode. Click the AND/OR control in the card's header: verify the label now reads
   OR, the filter reports OR mode, the same two tokens are still ticked, and the grid shows
   exactly the rows carrying EITHER token — more than the AND count, so a mode control that did
   nothing fails here. Click it again and verify AND and the smaller count come back exactly.
   Finally untick `blue` — the count returns to `red`'s own — and untick `red`, and verify nothing
   is included any more and all 5850 rows are visible again. Then remove the `MVF_TOKENS` column
   and its filter and verify the table's column list is back to what it was.
3. Narrow the SEX card to its largest category so the panel is filtering, and note the row
   count — verify it is above zero and below 5850, since an all-false or an all-true column
   would satisfy the equality below if it were not. Open the Filter Panel menu and click
   "Filter to Column...". A small dialog opens with a Name field already filled in from the
   filter that is currently active. Verify the proposed name is not blank and that the table
   does not already have a column with that name — otherwise the column appearing afterwards
   would prove nothing. Click OK. Verify: a column of that name now exists on the table, it is
   a boolean column, and counting its own true values — reading the column, not the filter —
   gives exactly the row count noted above; and the panel is still filtering the same rows.
   Finally remove that column again and verify the table's column list is back to what it was.

## Automation notes

- Signals: the ordered card-caption list, `df.filter.trueCount`, and the active-filter counter.
  Selectors, the caption and counter selector traps, dialog anatomy, menu-leaf names, canvas
  geometry, the Order or Hide Columns and Select Columns picker mechanics, the console-error rule
  and readiness barriers for everything below live in
  `.claude/skills/grok-browser/references/viewers/filters.md`.
- Scenario 1 Steps 6-7 narrowing:
  `fg.updateOrAdd({type: 'categorical', column: 'RACE', selected: [...]})` and the same shape for
  `DIS_POP`. The keys are `column` and `selected`; `DG.FILTER_TYPE` resolves only inside a
  `page.evaluate` body, so spell the value out as a string in the Node process. Record the resulting
  counts as `trueCount_race` and `trueCount_race_dispop`, guarding both against null/undefined at the
  point of capture, before any comparison.
- Scenario 2 Step 1: the GROK-12955 guard is the responsiveness round-trip
  (`await page.evaluate('1+1')` after the drop) plus the console-error text delta — an order-only check
  passes while the freeze is present.
- Scenario 4 empties the panel with "Remove All" and asserts it really is empty before every
  close/reopen, because a panel with a non-empty saved card set RESTORES that set instead of rebuilding
  it from column visibility, and only a rebuild consults the grid. Steps 1, 3 and 4 all read through
  that same sequence, so the three measure one behaviour.
- Scenario 4 hide/show gesture: the "Order or Hide Columns..." dialog's per-row check cell, driven as
  a real click — never a JS write to `visible` as a stand-in for the gesture.
- Scenario 5 close/reopen: close with the title-bar control and reopen with the ribbon funnel icon. Do
  not substitute `getFiltersGroup().close()` / `getFiltersGroup()` — the second call CREATES a panel
  when none exists, which manufactures its own green.
- Scenario 6 Step 1 records the row count after its opening "Remove All" as this scenario's baseline.
  Do NOT restore a hardcoded 5850, and do NOT paper over the difference with a `rows.requestFilter()`
  call — that would hide the behaviour rather than test it. The residue is GROK-20731, an OPEN bug
  (full anatomy in `bug-library/filters.yaml`): "Remove All" over a card RESTORED from a saved look
  leaves a live `onRowsFiltering` subscriber that no panel gesture can reach. Asserting that the picker
  and the Multi Value add leave THAT number alone is the claim these steps are actually making.
- Scenario 6 Step 1 drives the picker's All / None links rather than per-column picks, and compares
  committed captions as SORTED SETS, not order.
- Scenario 6 Step 2 opens a FRESH demog view and closes the previous one, because its arithmetic is
  whole-table and by this point the run's dataFrame carries the GROK-20731 residue. Deriving the
  expectations against the leftover would mean reading `df.filter` to build an expectation, which is
  circular, and hardcoding "RACE: Caucasian" would couple the step to a defect. The fresh view is
  asserted clean before anything else happens — exactly one `[name="viewer-Filters"]`, zero cards after
  "Remove All", `dataFrame.rows.filters` EMPTY (the leak detector), and `df.filter.trueCount === 5850`.
  Step 3 runs against that fresh view and is unaffected: its assertions are all self-derived.
- Scenario 6 Step 2 FIXTURE COLUMN. demog has no multi-value column, so the step declares one through
  the JS API (`df.columns.add(DG.Column.fromStrings('MVF_TOKENS', values))`) — the gesture under test is
  the filter, not column creation. The generation rule is deterministic and index-derived, with no
  `Math.random`, so a failure message's numbers are reproducible: row `i` takes
  `MVF_PATTERN[i % 10].join(';')` where `MVF_PATTERN` is
  `[[red,green],[red],[green,blue],[red,green,blue],[blue],[red,amber],[green],[red,green],[amber],[blue,amber]]`.
  demog's 5850 rows are an exact multiple of 10, so the distribution is exactly `red 2925, green 2925,
  blue 2340, amber 1755`, with `red AND blue` = 585, `red OR blue` = 4680, and no empty row. Those five
  numbers are distinct and strictly ordered (585 < 2340 < 2925 < 4680 < 5850), which is what makes the
  single-token, AND, OR and clear assertions each falsifiable. The spec does NOT trust that arithmetic:
  it recomputes every count by splitting the column's own cells on `;` and asserts the ordering
  properties BEFORE driving anything, so a fixture that silently degenerated fails at the fixture, not
  at the filter. Cleanup is a `finally` that removes both the multi-value filter and the column, then
  asserts the table's column list equals the pre-step list — the filter is removed too, so Step 3's
  `Filter to Column...` name proposal is derived from the SEX card alone.
- Scenario 6 Step 2 ticks the CHECK-BOX column (x ≈ left+10) throughout: a click in the name column
  REPLACES the selection instead of adding to it, which would make the two-token assertions impossible.
  Row order is not relied on — each row is ticked and unticked once to build the row -> token map, which
  is then asserted to equal the derived token set exactly.
- Scenario 6 Step 3 counts the TRUE values of the boolean column itself; deriving the expected number
  from the filter would be circular, and the "the column did not exist before" guard is what makes that
  count falsifiable. The column is removed in a `finally` and the table's column list asserted back to
  what it was — every later scenario would otherwise see a 12th column, in particular Step 1's "one card
  per column" assertion if the steps are ever reordered.

Spec must keep:

- Every drag is driven with the real Playwright mouse and never with synthetic HTML5 drag-and-drop:
  press on the source, make SEVERAL small moves until the accumulated distance passes the start
  threshold, wait for the drop indication, move the cursor into it, and release the button OVER it —
  the completing `onMouseUp` is bound to the drop zone, not to the source. The two start thresholds and
  the two drop-indication shapes are in the refdoc drag-and-drop section; the load-bearing consequence
  here is that the ADD drag waits on the `.d4-drop-zone` element appended to `document.body`, while the
  REORDER drag never gets that class and must wait on `document.body.d4-drag` plus the card strip's own
  colour change.
- The waiver class `gesture-uncontrollable-headless` is FORBIDDEN for all of these drag steps. The
  mechanism is plain mouse events on `document` with no `DataTransfer` anywhere, and the grid header's
  lack of a DOM handle is a coordinate problem — the press point is computed from the column's screen
  rect — not an impossible gesture. If a step still cannot be driven, diagnose and report it rather than
  writing it off.
