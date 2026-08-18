---
feature: filters
realizes_atlas:
  - filters.cp.save-and-reapply-state
realizes:
  - viewers.filters
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-20386
    status: fixed
realized_as:
  - save-and-reapply-state-spec.ts
scope_reductions:
  - id: SR-01
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      CRITERION setup over the JS API, not substitution of the gesture under
      test. Both filter CARDS are created through the real panel header column
      picker (mouse click on the header combobox, typed column name into
      `input.d4-column-selector-search-input`, Enter), and the card is asserted
      to have appeared FROM that commit — a picker that stopped adding cards
      fails the spec. What is not driven through the UI is the AGE card's
      numeric RANGE: it is set with `v.applyNumericFilter` (min 30, max 60)
      instead of through the card's Min / max inputs, and the Step 5
      perturbation likewise uses `v.applyNumericFilter` /
      `v.applyCategoricalFilter`. The RACE criterion, the one this scenario
      round-trips visually, IS set by a real canvas click on the category row
      (Step 2). The min/max input gestures are asserted end to end in
      panel-core-ladder Step 4/5 on the same build.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      The beer table in Step 8 is opened over the JS API
      (`grok.dapi.files.readCsv` + `grok.shell.addTableView`) rather than
      through the file browser. The step's subject is the type-matching rule of
      the Save or Apply submenu on a second, differently-typed table; how that
      table got on screen is setup. The submenu itself is opened and read
      through the real panel hamburger, scoped to the ACTIVE view's root, and
      the absence claim is now guarded by a positive control (the always-offered
      `Save...` leaf must be enumerated) so a submenu that failed to populate
      cannot read as "the state is correctly not offered".
    verdict_status: SCOPE_REDUCTION
  - id: SR-03
    check: E-ASSERTION-STRENGTH-01
    rationale: |
      The Teardown step asserts the spec's OWN cleanup, not product behaviour:
      the `finally` deletes the probe key from the `filter-states` localStorage
      entry and the step re-reads that entry to confirm the deletion. It is kept
      as a leak guard — a named state surviving the run would poison later runs
      of this same spec — and it is no longer able to pass on a corrupt store:
      the JSON parse is unguarded, so an unparseable `filter-states` fails the
      step instead of reading as "not leaked". No product claim is made by it.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: Setup
    expectation: >-
      The panel starts empty AND present: the Filter Panel's group header is on
      screen exactly once, and inside it no .d4-filter card is left after the
      reset — so "no cards" cannot be satisfied by the panel itself being gone,
      and the counter read at Steps 3 and 7 can only be produced by the two
      cards this scenario configures. The full demog row count is read from the
      runtime, is a number greater than zero, and every row is unfiltered.
  - anchor: Step 1
    expectation: >-
      Activating the panel header's column combobox opens the column picker
      popup (.d4-column-selector-backdrop); RACE is typed into the picker's
      search box and COMMITTED with Enter, and the card captioned RACE (read
      from .d4-filter-column-name) appears as a result of that commit — a picker
      that opens but adds nothing fails here. The popup is gone again,
      df.filter.trueCount still equals the full demog row count, and the header
      active-filter counter is present in the DOM and either hidden or reading 0
      — an absent counter element is a failure, not an empty reading.
  - anchor: Step 2
    expectation: >-
      The real click on the Asian category row MOVES df.filter.trueCount — the
      value after the click differs from the value before it, which is what
      proves the gesture landed on the row rather than missing it. The resulting
      count is strictly between 0 and the full demog row count, and the header
      active-filter counter reads exactly 1.
  - anchor: Step 3
    expectation: >-
      df.filter.trueCount is strictly below the full demog row count (5850) and
      above zero, and the header active-filter counter (.d4-filter-group-header
      .d4-filter-indicator) reads exactly 2. The counter tooltip summary table
      is recorded as the criteria baseline for Step 7 and is asserted non-empty
      HERE, at the point the baseline is taken: it lists exactly RACE and AGE,
      the RACE summary names the Asian category, and the AGE summary is a
      rendered numeric range of the form [min,max].
  - anchor: Step 4
    expectation: >-
      A named state entry is stored in window.localStorage under the
      'filter-states' key for the probe name — the entry itself must exist and
      parse. The Save dialog closes, and no NEW console error text appears
      across the step: the error texts accumulated before the save are compared
      by identity against those after it, with the standing ambient message
      ("Permissions policy violation: compute-pressure") excluded by name.
  - anchor: Step 5
    expectation: >-
      df.filter.trueCount moves to a value DIFFERENT from the value recorded in
      Step 3, confirming the modification is in effect before the re-apply.
  - anchor: Step 6
    expectation: >-
      The Save or Apply submenu, opened through the panel hamburger of the
      ACTIVE view and read scoped to it, enumerates its always-offered Save...
      leaf, and among its leaves is the saved probe name — confirming the state
      is offered (getApplicableStates matched the column types).
  - anchor: Step 7
    expectation: >-
      df.filter.trueCount RETURNS to the exact value recorded in Step 3. The
      active-filter counter (.d4-filter-group-header .d4-filter-indicator) reads
      2. The counter tooltip summary table lists exactly the same two columns
      (RACE and AGE) with the same criteria as at Step 3 — compared against the
      summary recorded at Step 3, entry for entry — and no other demog column
      appears in it. GROK-20386 guard: no
      "Cannot fire new event" error appears in the console across this step (the
      re-apply is driven through the real menu, not through applyState over JS).
  - anchor: Step 8
    expectation: >-
      The Save or Apply submenu on the beer table is confirmed to have
      ENUMERATED — its always-offered Save... leaf is among the leaves read —
      and the probe name is absent from it. No error is thrown; the menu is
      simply limited to states whose column types match.
  - anchor: Teardown
    expectation: >-
      The probe key is removed from window.localStorage so no named state leaks
      into later runs. This is a leak guard on the spec's own cleanup, not a
      product claim (see scope_reductions SR-03); a 'filter-states' entry that
      does not parse fails it rather than reading as "not leaked".
---

# Filters — Save and re-apply named filter state

Realizes `filters.cp.save-and-reapply-state`: configure at least two filters,
save the configuration under a probe name via the panel context menu, perturb the
state so the re-apply has observable work to do, then re-apply the saved name
through the real menu and assert the row count, the active-filter counter, and the
counter tooltip summary all return to the recorded values.

GROK-20386 (applying filters from the hamburger menu got broken) lives on the
re-apply step; driving through the real menu, rather than through the JS API, is
the regression guard.

## Setup

1. Open `System:DemoFiles/demog.csv` and ensure the table is the active table view.
   Record its full row count (expected: 5850).
2. Open the Filter Panel. The recommended UI path is the funnel icon on the toolbar
   ribbon. Wait for the panel to finish loading before proceeding.
3. Remove every filter card the panel opened with, and set the row filter back to
   all-true. The active-filter counter is asserted to read exactly 2 later on, so
   the two cards this scenario configures must be the only ones present.

   Expected (Setup): no filter card is left in the panel, and the full demog row
   count read from the runtime is a number greater than zero.

4. Choose a probe name that is unique per run (for example, append the current
   timestamp to a fixed prefix). This name is used in the Save dialog in Step 4 and
   in the re-apply step (Step 7) and must be cleaned up in the Teardown block.

## Scenarios

### Scenario 1: Save a named state, perturb, re-apply, assert round-trip

#### Step 1 — Add a categorical filter (RACE)

In the Filter Panel header, click the column selector (the plus / column picker at
the top of the panel). Choose the `RACE` column. A new categorical filter card for
RACE appears in the panel.

Verify the column picker popup actually opened on that activation, then dismiss it.
The column pick itself is made programmatically — the picker draws its rows on a
canvas with no per-row handle and the row a column occupies shifts with the cards
already present, so a blind geometric click could add the wrong column. This step
therefore exercises the entry point, not the column pick; see the Automation notes.

Verify the RACE card caption reads `RACE`. Verify the filtered row count is still
the full row count (no categories are checked yet) and the active-filter counter is
hidden or reads 0 (an added-but-untouched card does not count as active filtering).

#### Step 2 — Select one category in the RACE card

In the RACE card, click the row of the `Asian` category in the card's category list.
The card body paints its category rows rather than laying them out as separate elements,
so the click has to be aimed at the row's position inside the card body — the geometry
for that aim is in the Automation notes. This step must be a real click landing on the
category row; it is the interface gesture this scenario is here to exercise.

Only if the click cannot be delivered at all may the selection be set through the filter
state instead. That fallback must be reported as such: the run notes then say the category
was set programmatically, and the step is NOT to be described as a selection made through
the card.

Verify the click MOVED the filtered row count: the value read right after the click
must differ from the value read right before it. That difference is what proves the
gesture landed on the row — a click that misses leaves the count where it was, and
every count assertion below would then pass on an unfiltered table.

Verify the filtered row count drops below the full row count (to the `Asian` category's
row count, strictly less than 5850, and above zero) and the active-filter counter reads 1.

Expected (Step 2): the count before the click differs from the count after it; the
resulting count is strictly between 0 and the full row count; the counter reads 1.

#### Step 3 — Add a numerical filter (AGE) and narrow the range

Using the column selector in the panel header, add the `AGE` column. A numerical
(histogram) filter card for AGE appears in the panel.

Set a min/max window that excludes some rows — for example, restrict to ages 30–60
(driven programmatically via the range slider; see Automation notes).

Record the resulting filtered row count as **trueCount_saved** and note that the
active-filter counter in the panel header must read **2**.

Hover the counter and record its tooltip summary table as **summary_saved** — the
criteria baseline the re-apply is checked against in Step 7. Verify it here, where it
is taken: it must list exactly `RACE` and `AGE`, the RACE entry must name `Asian`, and
the AGE entry must be a rendered numeric range of the form `[min,max]`. A baseline
that is silently empty would make the Step 7 comparison satisfiable by an equally
empty tooltip.

Expected (Step 3): the filtered row count is strictly below 5850 and above 0, the
counter reads 2, and summary_saved is non-empty and well-formed as described.

#### Step 4 — Save the state under the probe name via the panel context menu

Right-click the Filter Panel header (or click the hamburger / gear icon) to open the
panel context menu. Navigate to **Save or Apply** > **Save…**.

In the save dialog, type the probe name from Setup step 4 and confirm. The dialog
closes.

Verify the probe name is now present among the stored named states (inspect the
browser's local storage for the filter-states entry; the probe name must appear
there). Verify the dialog is gone from the DOM, that the save added no platform
warning, and that no browser console error was emitted across this step.

Expected (Step 4): a named state entry exists in local storage for the probe name,
the dialog has closed, and neither a platform warning nor a console error is
recorded.

#### Step 5 — Perturb the state so the re-apply has visible work to do

Change the active configuration: widen the AGE range to include all rows and select
a different RACE category (for example, `Black`). Both changes are driven
programmatically via the filter group (see Automation notes).

Verify the filtered row count has changed to a value **different from trueCount_saved**
recorded in Step 3.

Expected (Step 5): the filtered row count is not equal to trueCount_saved.

#### Step 6 — Confirm the saved state is offered on the demog table

Open the panel context menu again. Navigate to **Save or Apply**. Verify the probe
name appears as a menu item in that group. The platform only offers a saved state on
a table whose columns match the saved column types, so its presence here confirms
both that the save succeeded and that the state is compatible with the current table.

Expected (Step 6): the probe name is visible as a menu item under **Save or Apply**.

#### Step 7 — Re-apply the saved state through the real menu and assert the round-trip

Click the probe name item in the **Save or Apply** submenu. This triggers the saved
state to be re-applied through the real panel context menu (not through the JS API —
driving through the API would bypass the GROK-20386 regression surface).

Hover the active-filter counter in the panel header and read its tooltip summary
table. The tooltip lists each active filter column and its criteria.

Expected (Step 7):
- The filtered row count equals trueCount_saved (exact value, not an approximation).
- The active-filter counter reads **2**.
- The tooltip summary table lists exactly **RACE** and **AGE** with the criteria from
  Step 3 (no extra columns, no missing columns). It is compared entry for entry against
  **summary_saved**, and no other demog column name may appear among its captions — a
  restore that hits the right row count through some other combination of criteria
  fails here.
- No `"Cannot fire new event"` text appears in the browser console across this step
  (the GROK-20386 regression guard: the re-apply fires a filter event and must not
  re-enter an already-firing event controller).

### Scenario 2: Saved state is not offered on a mismatched table (negative path)

#### Setup for Scenario 2

After Scenario 1 completes (probe name is in local storage), open a second table that
does NOT contain a `RACE` or `AGE` column — for example, `System:DemoFiles/beer.csv`.
Make that table the active view.

#### Step 8 — Confirm the saved state is NOT offered on the mismatched table

Open the Filter Panel for the beer table. Open the panel context menu and navigate to
**Save or Apply**. The probe name saved in Scenario 1 must NOT appear in the submenu
(or the **Save or Apply** group must be absent / empty).

This verifies that the platform's column-type guard works: a state saved for `RACE`
(categorical) and `AGE` (numerical) is silently filtered out for a table that lacks
those column types.

Expected (Step 8): the probe name is absent from the **Save or Apply** submenu on the
beer table. No error is thrown; the menu is simply limited to states whose column types
match.

## Teardown

In the finally block (runs even on failure):
- Remove the probe name entry from the named filter states in local storage.
- Close the beer table view if Scenario 2 was reached.
- Close the demog table view.

Expected (Teardown): the probe key is removed from local storage so it does not bleed
into later runs.

## Automation notes

Column-combobox mechanics, categorical canvas row geometry, the counter tooltip anatomy and its
50 ms debounce, the caption and counter selector traps, the demog RACE category sizes, the
no-confirmation-dialog fact for Reset and Remove All, the `filter-states` localStorage key and its
teardown, and the Save-or-Apply source anchors including GROK-20386 are in
`.claude/skills/grok-browser/references/viewers/filters.md`.

- **The column PICK is programmatic and must be reported as such.** Steps 1 and 3 actuate the
  combobox for real and assert the picker popup opened, then create the card through
  `fg.updateOrAdd({...})`: the picker's rows are canvas and a row's index depends on which columns
  already have cards, so a blind geometric click can add the wrong column and silently break the
  "counter reads 2" assertion. Only the entry point is exercised through the interface.
- **Step 2 is different — a real click on the card canvas, with no default fallback.** If the
  click genuinely cannot be delivered, the fallback
  `fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']})` must
  be labelled as programmatic in both code and run report; describing it as a selection made
  through the card is a false claim about what ran.
- **`Asian` is the deliberate choice** — demog's smallest RACE category (72 of 5850), so the
  row-count shift is the sharpest and cannot be confused with any other selection. There is no
  `White` category in demog.
- **Perturbation constants (Step 5):** AGE `min: 0, max: 200` (all rows) and a different category,
  `selected: ['Black']` — the perturbation has to be visibly different from the saved state or the
  re-apply has no work to prove.
- **The counter is polled, never read once**, and must hold `2` after Steps 3 and 7.
- **Probe hygiene:** the state name is `filter-state-${Date.now()}`, unique per run, so parallel
  runners cannot collide on it.

Spec must keep:
- The full demog size is read from the runtime (`grok.shell.tv.dataFrame.rowCount`), never hardcoded; in prose it is 5850 — never 500.
