---
feature: filters
realizes_atlas:
  - filters.cp.panel-core-ladder
  - filters.int.and-combination
  - filters.int.master-active-toggle
  - filters.int.active-counter-counts-filtering-only
  - filters.int.header-search-hides-cards
  - filters.int.esc-toggles-not-resets
realizes:
  - viewers.filters
  - viewers.filters.categorical
  - viewers.filters.histogram
priority: p0
target_layer: playwright
coverage_type: smoke
pyramid_layer: ui-smoke
related_bugs:
  - id: GROK-16677
    status: fixed
  - id: GROK-19152
    status: fixed
  - id: github-1103
    status: fixed
  - id: github-2307
    status: fixed
realized_as:
  - panel-core-ladder-spec.ts
scope_reductions:
  - id: SR-01
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      SETUP substitution, not substitution of the gesture under test. In
      Step 10, Step 11-negative and Step 13 the two-filter state that the
      persistence rungs round-trip is established over the JS API
      (`fg.updateOrAdd({type: CATEGORICAL, column: 'RACE', selected: [...]})`
      and `fg.updateOrAdd({type: 'histogram', column: 'AGE', min, max})`),
      not by re-driving the Step 3 canvas category click and the Step 5
      min/max inputs. The anchors for those steps say only that "the
      two-filter state is re-established"; what they assert is the SAVE
      (View > Layout > Save to Gallery, driven through the real menu leaf)
      and the ROUND TRIP, both of which remain fully UI-driven. The card
      gestures themselves are asserted end-to-end in Steps 3 and 5 of the
      same run, so nothing is left unexercised — only unrepeated.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      Layout RE-APPLY is over the JS API (`grok.dapi.layouts.find(id)` →
      `tv.loadLayout(saved)`) in Step 11 and Step 11-negative. The Layout
      menu offers no "apply this named layout" leaf — its commands are Save
      to Gallery, Download, Open Gallery, Clone View and Clear — and the
      gallery panel's thumbnails carry nothing that tells one saved layout
      from another, so the layout this run saved cannot be selected through
      the UI without guessing. The SAVE half of the round trip is driven
      through the real menu leaf and is what the anchors name; the re-apply
      is the unreachable half and is declared here rather than left implied
      by the "after re-applying the saved layout" phrasing.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: Setup
    expectation: >-
      Before the ladder starts, df.filter.trueCount is the full demog row count
      (5850). After the panel is emptied through its own "Remove All" menu leaf,
      the panel carries ZERO cards (both the caption list and the .d4-filter
      count read empty) and df.filter.trueCount is still 5850. A single
      BYSTANDER card (SEX) is then added through the panel header column
      selector and proven: the caption list is exactly ['SEX'], the card count
      is 1, df.filter.trueCount is still 5850, and the header active-filter
      counter block is PRESENT in the DOM but not on screen — it is hidden by
      `visibility` while nothing is filtering, so its computed style is what is
      read, never its rect. SEX is never touched again, which is what makes
      Step 6's "the summary does not list SEX" a claim about a card that is
      really in the panel.
  - anchor: Step 1
    expectation: >-
      After a real-mouse drag of the RACE grid column header onto the panel's
      "Add filter" drop zone, a card whose .d4-filter-column-name caption reads
      exactly "RACE" is in the panel. The drag is trusted page.mouse input from
      a press point computed off the grid overlay rect; an unchanged panel is
      the failure, not a tolerated outcome.
  - anchor: Step 4
    expectation: >-
      Adding AGE through the panel header column selector genuinely ADDS the
      card: the caption list does not contain AGE before the gesture, the header
      exposes a column-selector element with a non-zero box, the picker's search
      input holds exactly "AGE" before Enter commits it (a dropped keystroke
      must fail here rather than silently commit another column), and the AGE
      caption appears in the list afterwards.
  - anchor: Step 10
    expectation: >-
      The two-filter state is re-established (RACE narrowed to the clicked
      category, AGE windowed) and the header active-filter counter reads exactly
      2 before the layout is saved. The View > Layout > Save to Gallery leaf is
      really driven (driveTopMenuLeaf returns true) and the silent, auto-named
      save is located by diffing this user's applicable layouts: exactly ONE new
      layout id must appear, and the step throws rather than guessing — and
      deletes nothing — if more than one does.
  - anchor: Step 3
    expectation: >-
      After clicking the Black category row on the RACE card, the card's
      selection is exactly that one category — a plain click narrows the filter
      to the clicked category instead of dropping it — df.filter.trueCount drops
      below the full demog row count (5850) and the header active-filter counter
      (.d4-filter-group-header .d4-filter-indicator) reads exactly 1. This value
      is kept as the categorical-only row count that Step 5 is measured against.
  - anchor: Step 5
    expectation: >-
      After adding the AGE card through the panel header column selector and
      typing the upper end of the window into the card's max input field,
      committing each entry with Enter (an in-range end, then an out-of-range
      one, then an in-range one again), df.filter.trueCount is STRICTLY BELOW
      the categorical-only value recorded in Step 3 and the header active-filter
      counter reads exactly 2. The range control still answers the max field
      after the out-of-range entry — the row count moves again (github-2307
      regression guard).
  - anchor: Step 8c
    expectation: >-
      Esc on the focused panel flips the group off rather than resetting it:
      df.filter.trueCount goes to the full demog count (5850) and the panel root
      carries d4-filters-disabled, while BOTH criteria survive — the RACE card
      still reads exactly the clicked category and the AGE card still reports a
      non-null upper bound. A second Esc returns df.filter.trueCount to the
      exact pre-Esc value, drops d4-filters-disabled and brings the header
      counter back to 1. The surviving criteria and the exact restored value are
      what distinguish TOGGLE_FILTER_REQUEST from the reset icon's
      RESET_FILTER_REQUEST asserted in Step 9.
  - anchor: Step 6
    expectation: >-
      The counter tooltip is raised the way a user raises it — the pointer is
      parked away, the counter's box is re-read at that moment, and the pointer
      is walked onto it — and the tooltip that comes up must be a FRESH one (not
      hidden by its computed style, carrying text, and carrying text different
      from what the reused tooltip div held before the hover). Its table cell
      list is non-empty, and it lists exactly the two filtering columns (RACE
      and AGE) with their active criteria. The RACE criterion cell is non-empty,
      names the clicked category, and names none of the column's other
      categories (taken from the live column, not hard-coded). The AGE criterion
      cell is non-empty. No non-filtering card appears in the summary table —
      the SEX bystander added in Setup is named explicitly, so the clause can
      fail.
  - anchor: Step 7
    expectation: >-
      After unchecking the RACE filter card's own checkbox, df.filter.trueCount
      rises to the value the AGE-only filter produces (strictly above the
      two-filter value and strictly below the full count). The header
      active-filter counter drops to exactly 1. The RACE card is still present
      in the panel, carries the product's d4-filter-disabled class, and its
      criterion is intact — read back through the filter group, its selection is
      still exactly the clicked category. Re-checking the card returns
      df.filter.trueCount to EXACTLY the two-filter value it was disabled from
      and the counter to 2 (a criterion that came back diminished or widened
      could not reproduce that number); unchecking it again returns the exact
      AGE-only value and a counter of 1, leaving the panel with one disabled
      card for Steps 8a-8b.
  - anchor: Step 8a
    expectation: >-
      After clicking the master "Turn filters on/off" checkbox off,
      df.filter.trueCount returns to the full row count (5850), the panel root
      carries the d4-filters-disabled CSS class
      ([name="viewer-Filters"].d4-filters-disabled), and every card is still
      present BY CAPTION — RACE, AGE and the untouched SEX bystander (the master
      toggle disables filters, it does not drop cards) — with the criteria
      intact: read back through the filter group, the categorical card still
      selects exactly the clicked category, and the numerical card still reports
      a non-null min and max with max > min and max strictly below the AGE
      column's own maximum (a window reset to the full range reads back as null
      here, not as numbers that happen to bracket the data).
  - anchor: Step 8b
    expectation: >-
      After clicking the master checkbox back on, df.filter.trueCount returns to
      the exact pre-toggle AGE-only value recorded in Step 7. The
      d4-filters-disabled class is absent. The RACE card checkbox remains
      unchecked (its stashed isActive state is restored, not reset to on).
  - anchor: Step 9
    expectation: >-
      Typing a card-caption fragment into the panel header search (with real
      keys, into the input the header search icon reveals) hides the
      non-matching cards — the AGE card, the one that is actually filtering at
      this point, among them — while df.filter.trueCount and the header
      active-filter counter both stay exactly where they were (the header search
      is presentational; a hidden card keeps filtering). The counter is read
      before the search and must be non-empty, else the comparison across the
      search would be vacuous. Closing the search shows every card again and
      leaves df.filter.trueCount unchanged. After clicking the header reset
      icon, df.filter.trueCount is the full row count (5850), no confirmation
      dialog appears, the header active-filter counter is hidden or reads 0, the
      RACE, AGE and SEX cards are all present by caption, every card
      carries an enable/disable checkbox and all of them are checked, no card is
      left hidden behind the search, and the header search box is either
      collapsed or open and empty — the input must be present in the DOM for the
      check to count.
  - anchor: Step 11
    expectation: >-
      The perturbation is proven first — after the reset plus the unrelated
      viewer, df.filter.trueCount differs from the saved value, so "it came
      back" is falsifiable. After re-applying the saved layout (and waiting on
      the panel's STATE, never a clock), the Filter Panel is present in the
      view, the RACE and AGE filter cards are present, and the criteria saved in
      Step 10 are read back through the filter group: RACE selects exactly the
      clicked category, and AGE reports a non-null max strictly below the AGE
      column's own maximum. Every card's checkbox is in the state it held at
      save time — the checkboxes are COUNTED first (non-zero, and one per card)
      so that "all of nothing" cannot pass — and df.filter.trueCount equals the
      value recorded before saving. (github-1103 regression guard: the
      individual card checkboxes come back to their saved state, not
      blanket-on.)
  - anchor: Step 11-negative
    expectation: >-
      After the GROK-16677 rung — reset filters, type a search into the panel
      header, then save and re-apply the layout — df.filter.trueCount is the
      full row count (5850), the re-applied panel's header search input is
      present and empty, and no card is left hidden (no residue serialized from
      the prior search).
  - anchor: Step 13
    expectation: >-
      After reopening the saved project, the Filter Panel opens automatically
      (GROK-19152: the readiness barrier is keyed to the panel appearing in the
      DOM with a non-zero box, not a fixed delay), the RACE and AGE filter cards
      are present, their criteria are read back through the filter group of the
      reopened view (RACE selects exactly the clicked category; AGE reports a
      non-null max strictly below the AGE column's own maximum), and
      df.filter.trueCount equals the value recorded before saving the project.
---

# Filter Panel — Panel Core Ladder (p0)

Validates the full Filter Panel core as an accumulating ladder on the demog golden dataset:
categorical and numerical filtering composing via AND, the master toggle with stash/restore,
the header reset, and the mandatory p0 persistence peak (layout round-trip and project
round-trip). Every assertion is a row-count change visible in the status bar, the header
active-filter counter, or a panel state change visible in the panel itself. Covers the core
ladder critical path, the AND-combination interaction, the master toggle interaction, and the
active-counter-counts-filtering-only interaction.

## Setup

1. Open the demog dataset (System:DemoFiles/demog.csv) — 5850 rows, 11 columns.
2. Open the Filter Panel by clicking the funnel icon in the ribbon toolbar. Wait until the
   Filter Panel is visible with a non-zero size.
3. Confirm no filter is active: the status bar should show all 5850 rows selected before
   any filter is applied.
4. Empty the panel completely — open the Filter Panel menu from the hamburger icon in the
   panel's title bar and click "Remove All" — so that Step 1 and Step 4 genuinely add their
   cards rather than reuse cards that were already there. Dropping just the RACE and AGE
   cards by name is not durable: the panel re-populates its column list on later mutations
   and the two cards come back. "Remove All" fires immediately, with no confirmation dialog.
   Verify the panel really is empty before going on — zero card captions, zero filter cards,
   and all 5850 rows still in the dataset. A clear that silently did not happen would
   otherwise surface three steps later as a cascade of unrelated failures.
5. Add one BYSTANDER card — SEX — through the panel header column selector (the same gesture
   Step 4 uses for AGE, never the drag, which is Step 1's own subject). Verify the panel now
   holds exactly that one card, that no rows are excluded (the bystander does not filter),
   and that the header active-filter counter block, which is only shown while something is
   filtering, is present in the DOM but not on screen. SEX is never touched again by the
   ladder, which is what makes Step 6's "the summary does not list SEX" a claim about a card
   that is really in the panel rather than one that is true merely because the panel was
   emptied.

## Scenarios

### Scenario 1: Accumulating ladder — categorical then numerical, then toggle and reset

**Step 1.** Drag the RACE column header from the demog grid onto the Filter Panel's
"Add filter" drop target. Wait for a new filter card labeled "RACE" to appear in the panel.

**Step 2.** Click the row of the category labeled "Black" in the RACE filter card list. A plain
click on a category row narrows the filter down to that one category: the clicked category stays
selected and all the others are dropped from the selection.

**Step 3.** Verify:
- The RACE card now selects exactly one category — the one that was clicked. A falling row count
  on its own cannot tell this apart from the opposite reading (the clicked category dropped and
  the rest kept), so the selection itself has to be checked.
- The row count shown in the status bar drops below 5850.
- The header active-filter counter reads exactly "1".

Record this value as the categorical-only row count.

**Step 4.** Add an AGE numerical filter by clicking the column selector in the Filter Panel
header, waiting for the column picker to appear, typing `AGE` into the picker's search box
and committing the choice with Enter. Wait for the AGE filter card to appear.

**Step 5.** Set a min/max window on the AGE filter:
- Open the AGE card's filter menu by clicking the filter indicator on the card and choose
  **Min / max** to reveal the min and max input fields.
- Enter a maximum age value of `60` in the max field.
- Then deliberately enter an out-of-range maximum value (`999`) in the max field —
  this is the github-2307 regression step.
- After the out-of-range entry, enter a mid-range maximum (`55`) in the same field to confirm
  the range control still answers.

Verify:
- The row count is strictly below the categorical-only value recorded in Step 3.
- The header active-filter counter reads exactly "2".
- The row count changed again after the out-of-range entry (regression guard for github-2307:
  an out-of-range end value used to lock the range control).

Record this value as the two-filter row count.

**Step 6.** Hover over the header active-filter counter to reveal its tooltip summary table.
Park the pointer away first and note whatever the tooltip element already holds — the product
reuses one tooltip div and leaves the previous hover's markup in it — then walk the pointer onto
the counter. Verify:
- A tooltip really came up: it is not hidden by its computed style, it carries text, and that
  text differs from what it held before this hover. Residual content from an earlier hover must
  not be mistaken for a fresh summary.
- The tooltip's table cell list is non-empty. A tooltip that came up empty would satisfy the
  "does not contain" clause below without proving anything.
- The summary lists exactly RACE and AGE with their criteria. Each of the two criterion cells is
  non-empty, and the RACE criterion names the clicked category and none of the other RACE
  categories (taken from the live column, not from a hard-coded list).
- No other column appears in it — in particular the SEX bystander, which is genuinely in the
  panel and is not filtering, is absent.

**Step 7.** Disable the RACE filter card by unchecking its own enable/disable checkbox
(the checkbox at the top of the RACE filter card, always visible without hovering).

Verify:
- The row count rises from the two-filter value to a value strictly above it and strictly
  below 5850. Record as the AGE-only row count.
- The header active-filter counter reads exactly "1".
- The RACE card is still present in the panel with its criterion intact (the card caption
  still reads "RACE", and the card still selects exactly the clicked category).
- The RACE card is shown in its disabled state (it carries the product's disabled marker).

Then round-trip the card checkbox to prove the criterion was kept rather than merely
re-derivable: re-check it and verify the row count lands on exactly the two-filter value it was
disabled from and the counter reads "2" again — a criterion that came back diminished or widened
could not reproduce that number. Uncheck it once more and verify the exact AGE-only value and a
counter of "1" return, so Steps 8a-8b act on a panel that holds one disabled card.

**Step 8a — master toggle off.** Click the master "Turn filters on/off" checkbox in the
Filter Panel header to switch it off.

Verify:
- The row count returns to 5850.
- The panel visually shows the disabled state (all filter cards are grayed out or marked
  inactive).
- Every card is still visible in the panel, identified by caption — RACE, AGE, and the
  untouched SEX bystander: the master toggle disables filtering, it does not drop cards.
- The RACE and AGE criteria are intact: RACE still selects just the clicked category, and AGE
  still reports both bounds, with its upper bound above its lower one and below the AGE
  column's own maximum — a window reset to the full range would report no bounds at all.

**Step 8b — master toggle back on.** Click the master checkbox again to re-enable it.

Verify:
- The row count returns to exactly the AGE-only value recorded in Step 7 (the per-card
  active state is restored, not reset to all-on).
- The panel returns to its normal enabled appearance.
- The RACE card's own checkbox remains unchecked (it was disabled when the master was
  toggled, so its saved state is unchecked).

**Step 8c — Esc toggles, it does not reset.** Click an empty spot on the Filter Panel so the panel
has focus, then press Esc.

Verify:
- No rows are excluded from the dataset any more, and the panel shows its disabled state.
- Every criterion is still on its card: RACE still selects just the clicked category, and AGE still
  holds a window narrower than the column. This is what separates Esc from the header reset icon in
  Step 9 — Esc flips the group off and leaves the criteria alone, the reset destroys them. A build
  that wired Esc to reset passes a row-count-only check and fails right here.

Press Esc again.

Verify: the row count returns to exactly the value it had before the first Esc, the panel looks
enabled again, and the header active-filter counter reads "1". A reset could never produce that
exact value back.

**Step 9 — header search, then reset.** First read the header active-filter counter and confirm
it is not empty — without a value to compare, the "counter unchanged" check below would be
vacuous. Then open the panel header search (click the search icon in the Filter Panel header)
and type `RACE` into the box that appears, with real keys.

Verify:
- Fewer cards are shown than before, and every card still shown matches what was typed.
- The AGE card — the one that is actually filtering at this point — is among the hidden ones.
- The row count is unchanged — the header search hides cards, it does not filter data.
- The header active-filter counter is unchanged too. Hiding a card that is actively filtering does
  not stop it filtering; without this half the search could be a filter in disguise on a build where
  hiding a card silently drops its criterion.

Close the search by clicking the search icon again and verify every card is shown once more and
the row count is still what it was — collapsing the box undoes the hiding and nothing else.

Then click the header reset icon (the circular arrow icon in the Filter Panel header). The
reset fires immediately with no confirmation dialog.

Verify:
- The row count is 5850.
- No confirmation dialog is on screen.
- The header active-filter counter is absent or reads "0".
- The RACE, AGE and SEX cards are all still present in the panel, identified by their captions —
  the reset clears criteria, it does not remove cards.
- Every card carries its own enable/disable checkbox and all of them are checked again
  (regression guard for github-1103: count the checkboxes and assert their states, not only
  the row count, because the bug passes the row-count check while the per-card checkboxes are
  wrong — and a check that finds no checkbox at all is not a pass).
- No card is left hidden behind the search, and the header search box is either closed or
  open and empty. The search box must be present to be judged: a missing box is a failure of
  the check, not an empty value.

### Scenario 2: Persistence peak — layout round-trip (mandatory close of p0)

**Step 10 — set up state and save layout.** Bring the panel back to the two-filter state:
RACE narrowed to Black alone, and AGE limited by a maximum of 60. What is under test in this
scenario is the layout round-trip, so how that state is reached does not matter here — the
gestures themselves are already covered by Steps 2–5. Record the resulting row count as the
saved row count. Verify the header active-filter counter reads "2".

Save the current layout via **View > Layout > Save to Gallery**. The command saves silently
and names the layout after the table, so there is no name to type and no dialog to confirm —
identify the layout that was just saved by comparing the layouts available for this table
before and after the save. Exactly one new layout of your own must appear; if more than one
does, stop rather than guess which is yours.

**Step 10b — GROK-16677 negative rung.** Before the re-apply, also run the negative
variant in this same step:
1. Reset all filters (click the header reset icon) — verify the row count is 5850 and that
   the cards are still there (the reset clears criteria, it does not remove cards).
2. Open the panel header search and type `RACE` — verify the non-matching cards are hidden
   and the row count does not move. This is the residue the layout must not carry.
3. Save the layout again via **View > Layout > Save to Gallery**, identifying it the same
   way, and re-apply it after perturbing the filters.
4. Verify: the row count is 5850 (the reset state was saved, not the prior criterion), the
   panel header search box is present and empty, and no card is hidden.

Order note: this negative rung and the main round-trip of Step 11 are independent of each
other, and the automation runs the main rung FIRST, then this one — so there is nothing to
restore between them and no re-apply of the two-filter layout happens here. Read the two
rungs as a pair that may run in either order; each establishes its own state before saving,
which is what makes that possible. Both probe layouts are deleted in cleanup regardless of
which ran first.

**Step 11 — perturb and re-apply.** To make the re-apply meaningful, reset the filters and
add an unrelated viewer (for example, a Bar Chart). The live layout now differs from the one
that was saved — verify that before going on: the row count must no longer equal the saved
row count, otherwise "it came back" is unfalsifiable.

Re-apply the saved layout. The Layout menu offers no "apply this named layout" command — its
commands are Save to Gallery, Download, Open Gallery, Clone View and Clear, and the gallery
panel's thumbnails give no way to tell one saved layout from another — so the re-apply is
performed on the saved layout directly rather than through a menu.

Wait for the panel to finish rebuilding before checking anything: wait until the panel has
its cards back and the row count has stopped changing, not for a fixed number of seconds.

Verify:
- The Filter Panel is present in the view.
- The RACE and AGE filter cards are present, and their criteria are read back off the cards:
  RACE selects exactly the clicked category and nothing else, and AGE still reports an upper
  bound below the AGE column's own maximum. Card captions alone would also be produced by a
  panel that rebuilt itself empty.
- Every card carries its own enable/disable checkbox and all of them are checked (regression
  guard for github-1103: count the checkboxes and assert their states, not only the row
  count, because the bug passes the row-count check while the per-card checkboxes are wrong).
- The row count equals the saved row count recorded in Step 10.

Clean up: delete both probe layouts even if the scenario fails.

### Scenario 3: Persistence peak — project round-trip (mandatory close of p0)

**Step 12 — configure and save project.** From the filtered state (RACE narrowed to Black,
AGE max ~60, two-filter row count from Step 10), save the current view as a project (see
Automation notes for the mechanism this scenario actually scripts) under the name
"probe-filters-ladder-project".

**Step 13 — close and reopen.** Close all views. Then reopen the project from the
**Projects** browser. Wait until the Filter Panel is visible — regression guard for
GROK-19152: the readiness check must key on the panel appearing, not on a fixed delay.

Verify:
- The Filter Panel is open.
- The RACE and AGE filter cards are present, and their criteria are read back off the cards in
  the reopened view: RACE selects exactly the clicked category, and AGE still reports an upper
  bound below the AGE column's own maximum.
- The row count equals the saved row count recorded in Step 10.

Clean up: delete the saved project even if the scenario fails.

## Automation notes

Selectors, the column-picker typed-name commit, the min/max inputs, the header-search mechanics,
the `View > Layout` leaf inventory, the layout-apply readiness barrier and the counter tooltip
anatomy are in `.claude/skills/grok-browser/references/viewers/filters.md`.

- **Panel opening is setup-only and deliberately not the ribbon click.** The shared opener
  (`helpers/viewers.ts`, `openTable(..., {withFilterPanel: true})` → `openFilterPanel`) primes the
  panel over the JS API, waits for the first `.d4-filter` and hovers it. Nothing in the ladder
  asserts anything about the ribbon, so this is not a substitution that needs a fallback record.
- **Category row click (Step 6):** the column MUST be driven through the UI gesture or the exact
  JS-API equivalent of that click; an `applyState()` call is programmatic actuation, not a UI
  gesture, and is not an acceptable substitution here.
- **Primary signal discipline:** assert `df.filter.trueCount` on BOTH halves of every filtering
  step. The row count alone is not enough for the categorical rung — it falls both when the filter
  narrows to the clicked category and when that category is dropped — so the step also reads
  `getStates(column, 'categorical').selected` back and requires it to be exactly the clicked
  category. Never read the canvas option lists pixel-wise; assert their effect plus the filter's
  own computed `isFiltering`.
- **Card presence (Steps 8a, 9, 11) is asserted by CAPTION** (`.d4-filter-column-name` text
  contains RACE and AGE), never by `.d4-filter` count — demog auto-creates a card per suitable
  column, so a non-zero card count survives the destruction of exactly the two cards this
  scenario is about.
- **Column selector (Step 4):** Setup removes the auto-created AGE card first, so this step
  genuinely adds one.
- **Min/max (Step 5):** the out-of-range entry (999) is the github-2307 regression step; the guard
  is that a further in-range entry in the same field still moves the row count, i.e. the control
  was not locked. The histogram slider itself is canvas and is NOT dragged.
- **Reset (Step 9), github-1103 guard:** count the per-card checkboxes before judging them —
  assert the number found is non-zero AND equals the card count, and only then that all are
  checked.
- **Panel header search (Steps 9, 11-negative):** assert BOTH halves — fewer cards visible and
  `df.filter.trueCount` unchanged. After collapsing the box the honest check is "the input exists
  AND (it is collapsed OR its value is empty)", with the two states distinguished explicitly; a
  `querySelector(...)?.value ?? ''` reading cannot fail and must not be used.
- **Layout (Steps 10, 11):** SAVE goes through the real menu leaf via `driveTopMenuLeaf`
  (helpers/viewers.ts:767); APPLY is over the JS API. Perturb by resetting the filters and adding
  an unrelated viewer.
- **Readiness after a layout apply (Steps 11, 11-negative):** do not sleep — both rungs share the
  one polling barrier.
- **Project round-trip (Steps 12, 13):** save with `saveProjectViaApi` (`helpers/projects.ts`) —
  `DG.Project.create()` + `df.getTableInfo()` + `tv.getInfo()` as project children,
  `dapi.tables.uploadDataFrame`/`save`, `dapi.views.save(viewInfo)` before
  `dapi.projects.save(project)` — no ribbon click, no dialogs, no fixed sleeps or polling; the
  round-trip assertions (RACE/AGE filter card presence and criteria, row count) are ordinary
  Filters-viewer state, which `tv.getInfo()`'s `ViewInfo` restores the same as the ribbon Save
  does (verified live: a TableView carrying Grid + Pie chart + Filters round-trips all three
  viewers correctly via this path). `saveProjectViaUI` (real ribbon Save) is still required when
  a round-trip asserts state outside `@Prop`/layout serialization or the Save/Share dialog UI
  itself is under test. Reopen via `grok.dapi.projects.find(projectId)` → `project.open()`.
  GROK-19152 guard: the readiness barrier keys on `[name="viewer-Filters"]` becoming visible with
  a non-zero bounding box, never on a fixed delay — a fixed sleep satisfies readiness on small
  datasets and hides the bug on large ones.
- **Counter tooltip (Step 6):** assert the tooltip's cell list is non-empty before asserting what
  it does not contain. `fg.getFilterSummary()` returns the same table without proving it ever
  reached the screen.
- **Interaction coverage:** `and-combination` ← Step 5 (two active filters, row count strictly
  below either alone); `master-active-toggle` ← Steps 8a-8b (stash/restore with the per-card
  checkbox surviving); `active-counter-counts-filtering-only` ← Step 7 (counter 2 → 1 on disable
  while card count stays 2); `esc-toggles-not-resets` ← Step 8c (both criteria survive, exact
  pre-Esc count returns on the second Esc); `header-search-hides-cards` ← Step 9 first half.
- **Cleanup:** in `finally`, delete both probe layouts by the ids captured from the two gallery
  saves and call `deleteProjectWithCleanup` (helpers/projects.ts:1243) with the saved projectId,
  even on failure.
---
{
  "order": 1,
  "datasets": ["System:DemoFiles/demog.csv"]
}
