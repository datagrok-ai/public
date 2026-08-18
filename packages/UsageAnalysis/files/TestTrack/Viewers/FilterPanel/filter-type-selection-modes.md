---
feature: filters
realizes_atlas:
  - filters.cp.filter-type-and-selection-modes
realizes:
  - viewers.filters
  - viewers.filters.categorical
  - viewers.filters.histogram
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19915
    status: fixed
  - id: GROK-20270
    status: fixed
  - id: GROK-19878
    status: fixed
  - id: GROK-20242
    status: fixed
  - id: GROK-19537
    status: fixed
  - id: GROK-19897
    status: fixed
realized_as:
  - filter-type-selection-modes-spec.ts
scope_reductions:
  - id: SR-01
    check: E-TRACE-02
    rationale: |
      The categorical card body is a GridCore created inside the filter
      (grid_filter_base.dart:340) rather than a DG viewer, so its row pitch,
      header-band height and column boundaries are not readable from the JS
      API and not derivable from the DOM — the body is a single canvas with
      no per-row elements. The spec therefore computes candidate row centres
      from pinned constants (header band 10px, pitch 27px, centre offset
      13px, checkbox x 10px, name x 60px) instead of measuring them.
      Nothing is asserted on those numbers. Every click is accepted only by
      the OUTCOME it produces — the card must report exactly one checked
      category, with a non-empty label whose own row count is above zero —
      so a wrong pitch makes the step FAIL on the probe table it prints,
      never pass on a click that landed somewhere else. Rows are likewise
      addressed by the accepted index, not by an assumed category.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-TRACE-02
    rationale: |
      Step 11 asks for the zero-variance probe column to be added "in the
      grid" and Step 12 for it to be removed through the column header's
      context menu. Both are driven instead through the dataframe API
      (`columns.addNewCalculated` / `columns.remove`), because the column
      under test is a fixture for the filter card, not the gesture under
      test: the ticket the step guards (GROK-19897) is about the filter
      card's layout computation on degenerate data. The filter-side
      gestures the step is really about — creating the card and clicking
      its body — ARE driven through the UI: the card is added with the
      panel header's column selector, and the click is dispatched on the
      card's own canvas only after that canvas is asserted present with a
      non-zero box.
    verdict_status: SCOPE_REDUCTION
  - id: SR-03
    check: E-TRACE-02
    rationale: |
      Step 10 asks for the Asian rows to be selected in the grid and deleted
      through Edit > Delete selected rows. The spec instead selects them
      with `df.selection.init` and invokes the platform's own
      `CmdRemoveSelectedRows` function — the same command the menu item
      calls. The gesture the GROK-19537 guard is actually about, the Ctrl+Z
      undo, IS driven for real: the main grid is clicked to take focus and
      Control+z is pressed through the keyboard. The card-refresh claim is
      addressed to the CARD through its own geometry: every click lands on
      one of the card's rendered rows, and what comes back is the filter
      state that click produced (getStates(col,'categorical')[0].selected)
      rather than the column's category list. A card still rendering its
      stale three-row list therefore has no row that can put the restored
      category into that state, and the step fails on the row-to-category
      map it prints.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: Step 1
    expectation: >-
      The AGE filter card is now a categorical filter (its body shows category
      rows rather than a histogram with a slider). The per-card filter menu icon
      is visible on the card header immediately after the switch — GROK-19915
      guard: the icon must be present without any further click into the card
      body. df.filter.trueCount is unchanged at the full demog row count (5850);
      the active-filter counter is hidden or reads 0 (the card has no category
      checked yet).
  - anchor: Step 2
    expectation: >-
      Both categories are checked by CLICKING their rows in the card body — the
      first by clicking its category name (an exclusive select), the second by
      clicking its checkbox (a toggle that adds to the selection). Neither row is
      chosen by position: the card lists the BLANK category first on this data
      (demog has exactly one row with an empty AGE), so a row taken by index can
      select the empty label whose own row count is 0, which makes every count
      comparison below unfalsifiable. Both rows are therefore chosen by OUTCOME —
      candidate rows are clicked in order and a row is accepted only when the
      card reports exactly one checked category, that label is non-empty and its
      own row count is above zero; the second row must additionally yield a
      DIFFERENT category with a DIFFERENT own count so Step 4 can discriminate.
      If no row in the scanned range is usable, the step FAILS and prints the
      probe table of every row tried and the state it produced. Every click
      lands on a real category row: the click y is a measured row centre, below
      the card grid's header band, because a click in that band reaches no row
      and changes nothing. Each click is proved to have landed: the card's
      checked-category list must actually name the clicked category, and a click
      that leaves the selection unmoved fails the step rather than falling back
      to a JS-API assignment. The checkbox click is also proved not to have
      taken the invert() branch — a card left unfiltered at 5850 with no
      selection fails rather than reading as "nothing happened".
      df.filter.trueCount drops below 5850 and equals the sum of the two
      categories' own row counts, computed independently off the AGE column. The
      card's checked-category list (fg.getStates('AGE', 'categorical')[0]
      .selected) holds exactly those two categories. The active-filter counter
      (.d4-filter-group-header .d4-filter-indicator) reads 1. Record the reduced
      value as trueCount_two_cats (it is the combined row count for the two
      checked categories). Failure messages carry the read values — the clicked
      categories, their own counts, and the pre/post card state — so a miss is
      diagnosable from the run report alone.
  - anchor: Step 3
    expectation: >-
      The card is in Radio mode. The category selection collapses to exactly ONE
      category (the first checked one survives; the second is dropped —
      cat_filter.dart:68-72). df.filter.trueCount equals that single category's
      own row count, which is strictly less than trueCount_two_cats, AND the
      card's own checked-category list (fg.getStates('AGE', 'categorical')[0]
      .selected) holds exactly that one category — the count alone does not say
      WHICH category survived. The Select
      all, Deselect all, and Invert all items are ABSENT from the card's
      indicator menu — GROK-20270 and GROK-19878 guard: the three items are
      omitted entirely in Radio mode (grid_filter_base.dart:137-141); assert
      their absence by opening the indicator menu and verifying none of the
      three labels appears.
  - anchor: Step 4
    expectation: >-
      The second category is selected by a REAL CLICK on its row in the card
      body, not by assigning the selection through the JS API: only the click
      enters the Radio replacement branch (grid_filter_base.dart:405-414), so an
      assigned selection could not fail on the behaviour this step names. The
      card's pre-click state is captured in the SAME reading as the click and
      asserted to be a known NARROWING state — the first category alone, at its
      own row count, not the full 5850 — and the click is then proved to have
      landed by a pre/post DELTA on that reading: both the checked-category list
      and df.filter.trueCount must move. Nothing is compared against a value
      carried over from an earlier step, so a no-op click fails here instead of
      matching the state Radio mode was already holding. Clicking the second
      category replaces the first rather than adding to it. df.filter.trueCount
      equals the second category's own row count (not the union of the two, and
      not 5850), and the card's checked-category list holds exactly the second
      category — one entry, not two. If the selection moves but the row count
      does not follow it, the step FAILS and reports both readings — that
      combination is never accepted here. The active-filter counter still
      reads 1.
  - anchor: Step 5
    expectation: >-
      The card is back in Multi-Select mode, and all three batch items — Select
      all, Deselect all and Invert all — are present again in the card's
      indicator menu. That is the other half of the Step 3 absence check: their
      absence in Radio mode only means something if they are demonstrably there
      in Multi-Select. The ordered list of card captions read from
      .d4-filter-column-name still shows the AGE-derived categorical card as
      expected. df.filter.trueCount is still at the second category's row count;
      the active-filter counter reads 1.
  - anchor: Step 6
    expectation: >-
      The matched-category count inside the card changes to reflect only the
      categories matching the typed fragment (the match is a substring regexp —
      GROK-20242 anchor: a mid-word fragment must match). df.filter.trueCount
      drops or stays dropped from Step 5, consistent with enableFoundCategories
      narrowing the main dataframe (grid_filter_base.dart:181). The value must
      differ from the pre-search trueCount unless the search matches all
      categories. The surviving rows are also the RIGHT rows: the set of distinct
      AGE labels still passing equals the set of labels containing the fragment,
      counted in both directions — no passing row's label lacks the fragment and
      no excluded row's label carries it. The fragment is checked to match at
      least one and fewer than all categories BEFORE the comparison, so a search
      that did nothing could not satisfy it.
  - anchor: Step 7
    expectation: >-
      The pasted block is accepted as a LIST of category terms, not as one
      literal search string. df.filter.trueCount equals the UNION of the rows of
      the distinct pasted category values — for categories of a single column
      that union is the arithmetic sum of their individual row counts. Negative
      half, both directions: the value is NOT the full demog row count (5850)
      and NOT the row count of any one pasted value taken alone — it is strictly
      greater than the largest individual pasted value's count and strictly less
      than 5850. The trailing separator (line break or comma) contributes no
      condition: the count obtained with the trailing separator equals the count
      obtained without it, and neither collapses to 0 or to 5850. The
      active-filter counter (.d4-filter-group-header .d4-filter-indicator) reads
      1. GROK-20242 anchor: multi-line paste of a compound value list.
  - anchor: Step 8
    expectation: >-
      trueCount_select_all equals the full demog row count (5850).
      trueCount_deselect_all equals 0. trueCount_invert equals the complement of
      trueCount_deselect_all (5850 - trueCount_deselect_all), because inverting
      the empty selection restores the full set — so trueCount_invert equals
      trueCount_select_all and only trueCount_deselect_all differs from the full
      row count. The active-filter counter follows the same three readings: it is
      hidden after Select all (the card excludes nothing), reads 1 after Deselect
      all (the card excludes everything), and is hidden again after Invert all.
  - anchor: Step 9
    expectation: >-
      The AGE filter card is a histogram again (a numerical card with a slider
      or histogram body). df.filter.trueCount returns to the full demog row
      count (5850). The active-filter counter is hidden or reads 0. This is the
      round-trip revert rung — the type switch back to histogram must undo the
      categorical selection.
  - anchor: Step 10
    expectation: >-
      After Ctrl+Z restores the deleted rows, the categorical filter card for
      RACE re-lists its categories including any category that was absent while
      those rows were deleted. The category count returned must equal the count
      before the deletion. df.filter.trueCount is consistent with the restored
      row set. GROK-19537 guard: the undo direction must refresh the card, not
      only the delete direction.
  - anchor: Step 11
    expectation: >-
      A zero-variance constant column is added, a filter is created on it, and
      clicking it produces no new console errors. The console-error delta across
      this step (measured by capturing error texts before the step and comparing
      by identity) is zero. df.filter.trueCount stays plausible. GROK-19897
      guard: the degenerate-data layout computation must not throw.
  - anchor: Step 12
    expectation: >-
      The probe column is removed from the table and is absent from the Filter
      Panel. No error is produced during removal. df.filter.trueCount is
      unchanged from before the probe column was added.
---

# Filters — Filter type switching and category selection modes

Verifies that the Filter Panel correctly handles switching a filter card between
numeric histogram and categorical modes, Radio vs Multi-Select category selection,
in-card category search, batch category operations (Select all / Deselect all /
Invert all), round-trip type revert, and two regression guards: GROK-19537
(categorical filter not refreshed after Ctrl+Z row restore) and GROK-19897
(no console errors on a zero-variance column filter).

## Setup

1. Open `System:DemoFiles/demog.csv` and ensure the table is the active table
   view. Record its full row count (expected: 5850).
2. Open the Filter Panel using the funnel icon on the toolbar ribbon. Wait for
   the panel to fully appear before proceeding.
3. Using the column selector at the top of the Filter Panel (the plus button or
   column-combo box), add the `AGE` column. A numerical (histogram) filter card
   for AGE appears in the panel. Verify the card count is 1 and the total row
   count shown in the table is still 5850 (adding an untouched card does not filter).

## Scenarios

### Scenario 1: Type switch — numeric histogram to categorical, then back

#### Step 1 — Switch the AGE card to categorical and verify the menu icon

On the AGE filter card, locate the **filter menu** icon in the card header (the
three-dot or gear icon). Click it to open the card's filter menu and select
**Switch to categorical filter**. The card is destroyed and rebuilt as a
categorical card.

Immediately after the switch — without clicking anywhere inside the new card body —
verify that the filter menu icon is present on the rebuilt card header.

Verify the row count shown in the table is unchanged at 5850 and the active-filter
counter on the panel header is hidden or reads 0 (no category is checked yet).

Expected (Step 1): the AGE categorical card is visible, the filter menu icon is
present on the header without any extra click, and the row count is unchanged.

#### Step 2 — Check two discrete categories and verify the AND-combination

Inside the categorical AGE card, check two category values by **clicking their rows
in the card body** — the category list is painted on canvas, so the rows are reached
by clicking the card's overlay at the row position (the measured header band, row
pitch and column boundaries are in the Filter Panel UI reference). A click on a row's name
cell checks that category alone. A click on a row's checkbox is a toggle in
Multi-Select, so it adds a category and both remain checked simultaneously.

**Do not pick the rows by position.** The card lists the blank category first on this
data — demog has exactly one row whose AGE is empty — so the top row selects the empty
label, whose own row count is 0, and every count comparison built on it proves nothing.
Choose both rows by their OUTCOME instead: click candidate rows in order on the name
cell and accept the first one whose resulting selection is exactly one category, with a
non-empty label and an own row count above zero. Then find the second the same way,
below the first, requiring a different category whose own row count also differs — Step
4 discriminates on that count. If no row in the scanned range is usable, the step fails
and reports every row tried together with the state it produced.

Every click must land on a real category row. The top of the card body is the internal
grid's **header band**, not the first category — a click there does nothing at all and
looks exactly like a failed gesture. Address rows by their measured centres, never by a
small y that happens to be near the top of the card.

Confirm each click actually landed before drawing any conclusion from it — the card's
checked-category list must name the category that was clicked. A click that leaves the
selection unmoved is an actuation failure and fails the step; do not substitute an
API-assigned selection for it. Guard the second look-alike no-op too: clicking the
checkbox of the only currently checked category flips the card back to "everything
selected", where the row count returns to 5850 and the selection disappears. Reading
only the selection, that is indistinguishable from nothing having happened, so assert
the row count as well.

Note each clicked category's own row count, worked out independently from the AGE
column, and record the filtered row count as **trueCount_two_cats**.

Verify the filtered row count equals the sum of the two categories' own row counts,
that the card's checked-category list holds exactly those two categories, and that
the active-filter counter on the panel header reads 1.

Expected (Step 2): the two categories are checked by clicking their rows, the filtered
row count is the sum of their own counts and strictly below 5850, and the counter
reads 1.

#### Step 3 — Switch to Radio mode and assert single-selection invariant

Open the card's indicator or filter menu and navigate to **Mode > Radio** (or
the equivalent Radio / combo-box option) to switch the card to Radio mode.

Open the card's indicator menu and verify that **Select all**, **Deselect all**,
and **Invert all** are absent from the menu. Do not invoke them — assert their
absence by the fact that none of those labels appears in the opened menu.

Verify that the multi-category selection collapsed to exactly one category
(the first checked one): the filtered row count equals that single category's row
count, which is strictly less than **trueCount_two_cats**. Read the criterion back
as well as the count — the card's own checked-category list must hold exactly that
one category. A row count on its own does not say which category survived.

Expected (Step 3): Radio mode is active, the batch operation items are absent from
the indicator menu, and the filtered row count equals one category's row count.

#### Step 4 — Click a second category in Radio mode and assert replacement

In Radio mode, click a different category row in the card — the same overlay click on
the row's checkbox position used in Step 2. This click is the gesture under test: only
a click goes through the Radio replacement path, so the selection must NOT be assigned
through the API here. In Radio mode clicking a category replaces the previous selection
rather than adding to it.

Before the click, read the card's checked-category list and the filtered row count, and
confirm the card is in a **known narrowing state**: the first category alone, at that
category's own row count, not the full 5850. Then click, and confirm the click landed by
comparing against those two just-read values — both the checked category and the row
count must move. Do not compare against a number recorded in an earlier step: in Radio
mode the state on entry can already equal the state you expect afterwards, and such a
check passes without the click doing anything.

Verify the filtered row count equals the clicked category's own row count, NOT the
union of the two categories and NOT 5850. Verify the card's checked-category list holds
exactly one entry and that entry is the category just clicked. Verify the active-filter
counter reads 1.

If the checked category moves but the filtered row count does not follow it, the step
fails and reports both readings. A 2026-08-17 recon run saw exactly that; the operator drove
both entry states by hand on 2026-08-18 and the count followed the selection in both, so no
product defect is claimed. The guard stays as a regression guard and must not be relaxed.

Expected (Step 4): the second category replaces the first; the filtered row count
equals the second category's count, not their union.

#### Step 5 — Return to Multi-Select mode

Open the card's indicator or filter menu and navigate back to **Mode > Multi-Select**
(the default mode) to return the card to multi-category selection.

Verify the card is in Multi-Select mode. Open the card's indicator menu again and
verify that **Select all**, **Deselect all** and **Invert all** are all present —
this is the other half of the Step 3 check, since their absence in Radio mode only
means something if they are demonstrably offered in Multi-Select. Verify the
filtered row count is still the second category's row count from Step 4. Verify the
active-filter counter reads 1.

Expected (Step 5): Multi-Select mode is active, the three batch items are back in
the indicator menu, and the row count is unchanged from Step 4.

#### Step 6 — In-card category search moves the main dataframe

Clear any checked categories in the AGE categorical card (click Select all to
check all, then clear the search box if needed). Record the filtered row count
after the clear.

In the card's search box, type a partial text fragment that matches a subset of the
age-bucket labels (for example, a digit fragment). Observe the filtered row count
after typing.

Verify the matched-category count inside the card changes. Verify the filtered row
count changes from the pre-search value, because the in-card search under the
"Filter by search results" checkbox default narrows the main dataframe.

Then check WHICH rows survived, not only how many: the set of distinct category
labels still passing must equal the set of labels containing the typed fragment,
counted in both directions — no passing row's label may lack the fragment, and no
row whose label carries the fragment may have been excluded. Before that comparison,
confirm the fragment matches at least one and fewer than all categories; otherwise
a search that did nothing at all would satisfy it.

Expected (Step 6): the filtered row count moves when the in-card search box is used
(the default "Filter by search results" checkbox is on); the match is a substring
and a mid-word fragment must match at least one category; and the surviving rows are
exactly the rows of the matching categories.

#### Step 7 — Paste a multi-line list of values into the card search

Clear the in-card search box and uncheck any checked categories, so the card lists
all of its categories and the whole table is unfiltered again. Record the row count
after the clear (expected: 5850).

Pick two or three category values that this card actually lists and that are not
substrings of one another. Note the row count of each of them separately.

Now **paste** into the card's search box — do not type it character by character —
a single block of text holding those values, one per line, ending with a trailing
line break. This is what a user gets when copying a few cells of a column out of
the table. How the paste gesture is reproduced is in the Filter Panel UI reference.

Verify the filtered row count equals the union of the rows of the pasted values —
their individual row counts added up — and is neither the count of just one of them
nor the full 5850. Verify the active-filter counter on the panel header reads 1.

Verify the trailing line break did not create an empty condition: repeat the paste
without the trailing line break and confirm the filtered row count is the same as
with it (and is still neither 0 nor 5850).

Expected (Step 7): pasting several category values separated by line breaks filters
to the union of those values, a trailing separator changes nothing, and the
active-filter counter reads 1.

#### Step 8 — Batch operations in Multi-Select mode

Clear the in-card search box (so all categories are visible). Open the card's
indicator menu. Select **Select all**. Record the filtered row count as
**trueCount_select_all**.

Open the indicator menu. Select **Deselect all**. Record the filtered row count as
**trueCount_deselect_all**.

Open the indicator menu. Select **Invert all**. Record the filtered row count as
**trueCount_invert**.

Verify **trueCount_select_all** equals 5850 (all rows pass when every category is
checked). Verify **trueCount_deselect_all** equals 0 (no category checked, so no row
passes). Verify **trueCount_invert** is the complement of **trueCount_deselect_all** —
5850 − 0 = 5850, because inverting an empty selection restores the full set. So
**trueCount_select_all** and **trueCount_invert** are both 5850 and only
**trueCount_deselect_all** differs; the three readings are not three distinct numbers.

Verify the active-filter counter on the panel header follows the same three states:
hidden after Select all (the card excludes nothing, so nothing is filtering), reading 1
after Deselect all (the card excludes everything), and hidden again after Invert all.

Expected (Step 8): Select all at 5850, Deselect all at 0, Invert all at the complement
of Deselect all (5850); the counter is hidden, then 1, then hidden again.

#### Step 9 — Switch back to histogram and verify round-trip revert

Open the card's filter menu and select **Switch to histogram filter**.
The card is rebuilt as a numerical histogram card.

Verify the filtered row count returns to the full demog row count (5850). Verify
the active-filter counter is hidden or reads 0.

Expected (Step 9): the AGE card is a histogram, and the filtered row count is 5850.

### Scenario 2: GROK-19537 — categorical filter refreshes after Ctrl+Z row restore

#### Setup for Scenario 2

Using the column selector, add the `RACE` column to the Filter Panel. A categorical
RACE card appears. Verify the card count grows by 1 and the filtered row count is
unchanged.

Open the card's indicator menu and select **Select all** to ensure all RACE
categories are checked. Record the current RACE category count (the number of
distinct categories shown in the card).

#### Step 10 — Delete rows, verify card updates, undo, verify card restores

In the grid, select all rows belonging to a RACE category that has a small count
(for example, `Asian` with 72 rows). Delete those rows using the grid's row-delete
action or menu (Edit > Delete selected rows). Verify the RACE card's category list
no longer shows the deleted category.

Press Ctrl+Z to undo the deletion. Wait for the undo to complete.

Verify the RACE card's category list is restored and the deleted category reappears.
Verify the category count returned to the pre-deletion count. Verify the filtered
row count is consistent with the restored row set.

Expected (Step 10): after Ctrl+Z, the RACE card re-lists all original categories
including the one that was absent while the rows were deleted; the category count
matches the pre-deletion count.

### Scenario 3: GROK-19897 — zero-variance column produces no console errors

#### Step 11 — Add a constant calculated column and filter it

Capture the current console error texts as a baseline (record each error's message
string before this step). The baseline captures any ambient errors already present.

In the grid, add a calculated column whose value is the same constant in every row
(for example, the expression `1` or a literal string). Name it `probe_constant`.

Using the column selector, add a filter for `probe_constant`. Click a value in the
`probe_constant` filter card.

Compare the console errors after this step against the baseline captured before.
Assert the delta is zero: no new error messages whose text was absent from the
baseline appear after the click. Compare by error message identity, not by count,
so ambient noise is excluded.

Expected (Step 11): no new console errors appear across this step; the zero-variance
column filter renders and responds to interaction without throwing.

#### Step 12 — Remove the probe column (teardown for Scenario 3)

Remove the `probe_constant` column from the table (right-click the column header >
Remove column, or use the equivalent menu action). Verify the `probe_constant`
filter card is removed from the Filter Panel. Verify the filtered row count is
unchanged from its value before the probe column was added.

Expected (Step 12): the probe column is absent from the table and the panel, and
no error is produced during removal.

## Teardown

In the finally block (runs even on failure):
- Remove the `probe_constant` column if it was added and not yet removed.
- Close the demog table view.

## Automation notes

Panel opening, the column selector, the card canvas row geometry, the two look-alike no-ops,
addressing category rows by the outcome of the click, the demog RACE category sizes, the in-card
search defaults, the multi-line paste branch and the console-error rule are in
`.claude/skills/grok-browser/references/viewers/filters.md`.

- **Panel opening:** the funnel icon on the toolbar ribbon. The Toolbox entry may not work (see
  chain `unresolved_ambiguities`).
- **Reaching a checked-category state:** by real overlay clicks only.
  `fg.updateOrAdd({..., selected})` never enters the Radio replacement branch, so an assigned
  selection would leave Step 4 unable to fail on the behaviour it names. Each step therefore
  asserts a pre/post DELTA rather than a value an earlier step left behind, and the wait after a
  click is a barrier on that reading (leave once it has moved and held), not a fixed sleep.
- **Picking the two rows Steps 3-5 use:** rows are accepted by the outcome of the click, and the
  search for the SECOND row continues BELOW the first and additionally requires a different
  category with a different own count, because Step 4 discriminates on that count. Each accepted
  row's index and computed y are recorded so later steps address the same row. There is no API
  fallback and no shifting to a fixed index: if no row in the scanned range is usable, the step
  fails with the probed rows and their observed states in the message.
- **Step 4 establishes a known narrowing pre-state before the click**, and if the count then still
  does not follow the selection the step fails and names both readings. Do not weaken it. The
  operator drove both entry states by hand on 2026-08-18 and found no defect.
- **Data (demog, 5850 rows):** the categorical view of AGE shows individual integer values or age
  buckets depending on the auto-range resolution — pick two small-count values to maximize the
  shift. `Asian` is the deletion category in Step 10 (smallest RACE count, sharpest shift).
- **Type switch:** driven by the per-card hover-only icons, never by setting the filter type through
  the API. The menu items are "Switch to categorical filter" and "Switch to histogram filter"
  (`filters_core.dart:827-843`).
- **Multi-line paste (Step 7)** is the GROK-20242 gesture and must NOT be waived as "clipboard is
  not automatable" — a single-line comma-separated paste does not exercise the ticket at all. The
  values pasted must not be substrings of one another, and their per-value counts are computed off
  the AGE column at run time so the union is compared against an independently derived number.
- **Row deletion for GROK-19537:** the grid's row-delete action (select rows, then Edit > Delete
  selected rows or equivalent); Ctrl+Z triggers `cat_filter.dart:110-138`.
- **Derived from:** `cat_filter.dart:63-78` (Radio mode), `filters_core.dart:827-843` (type-switch
  icons), `grid_filter_base.dart:355-363` (Radio replacement), `:137-141` (batch-ops guard).
