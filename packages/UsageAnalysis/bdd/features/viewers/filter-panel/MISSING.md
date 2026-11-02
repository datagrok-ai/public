# Filter panel — what the features cannot say yet

## Closed issues the regression round of 26 Oct 2026 could not claim

Written during the closed-issue regression round of 26 Oct 2026 (bucket b5, filtering), against the
panel's "Automation surface" (`core/client/d4/lib/src/viewers/filters/CLAUDE.md`). What could be said
went into `filter-panel-regressions.feature` (public#790, #791, #3059, #1464, #319, #2501, #2943) and
the panel's own journeys (public#1986 ladder; #3027, core#479, #1889 categorical; #1343 range-handles;
#377 persistence; #1837 cloned-view). public#2364 is an ApiTests `FilterGroup` test (the min and max
fields clamp a typed end); public#1103 is filter-panel-cloned-view's per-card layout claim and the
ladder's master toggle. Chem's substructure-card candidates of the same bucket are in
`Chem/bdd/features/filters/`.

### Still open

- **A click keeps the clicked category after the column's custom order changed**
  ([public#576](https://github.com/datagrok-ai/public/issues/576)). The custom order is set in the
  grid header's `Sort > Custom...` dialog (`editColumnCategoryOrder`,
  `core/client/d4/lib/src/common/column_category_ordering.dart`), which reorders the categories by
  dragging; no step moves an item of that dialog, and setting the column's order tag through the API
  would skip the dialog the issue goes through. Needs a gesture for the dialog's list (`user drags the
  "Black" item above the "Asian" item in the category order dialog`, or the dialog's own up/down
  icons named) — a library step, as other ordered lists want it. Then:

  ```gherkin
  Scenario: A click keeps the clicked category after the column's custom order changed (public#576)
    When user adds a card for "RACE" to the filter panel
    And user picks "Sort > Custom..." from the context menu of the "header RACE" area of grid
    And user drags the "Other" item above the "Asian" item in the category order dialog
    And user clicks on OK button in category order dialog
    Then the "categories of RACE" reading of filter panel should be "Other, Asian, Black, Caucasian"
    When user clicks on the "category Black of RACE" area of filter panel
    Then the filter should pass exactly the rows where "RACE" is "Black"
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And no errors should have been logged
  ```

  (The `categories of` reading lists the card's `nameCol` in its own order — check that it follows
  the custom order before relying on it; otherwise the claim is `category Other of RACE` lying above
  `category Asian of RACE`.)

- **Hovering a category highlights its rows in a scatter plot**
  ([public#2344](https://github.com/datagrok-ai/public/issues/2344)). A hover on a category sets the
  table's `mouseOverRowFunc` (`cat_filter.dart`), and nothing reports that group: the library reads
  only the single `mouse-over row`, and `scatter plot viewer should have repainted` is one pixel, not
  the 27 Black rows. Needs a core signal — a `mouse-over rows` count on the table's status (or the
  scatter plot's `highlighted rows` reading, as Peptides' status provider has one) computed from
  `rows.mouseOverRowFunc` when asked. Then:

  ```gherkin
  Scenario: Hovering a category highlights its rows in a scatter plot (public#2344)
    Given user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    When user adds a card for "RACE" to the filter panel
    And user hovers over the "category Black of RACE" area of filter panel
    Then the "highlighted rows" reading of scatter plot viewer should be 27
    When user moves the pointer away from grid
    Then the "highlighted rows" reading of scatter plot viewer should be 0
    And no errors should have been logged
  ```

- **Typing in the panel's add-filter column picker narrows its list**
  ([public#605](https://github.com/datagrok-ai/public/issues/605)). The picker is the header's
  `ColumnComboBox` named `div-column-combobox-add-filter` (`filters_core.dart`). `user types … into the
  … column selector of …` types and picks (`typeInColumnGrid`), and `… should offer the columns …` opens
  the picker and reads all of it with nothing typed; no step types into an open picker and reads what
  is left. Needs a library step in `widgets.ts`, e.g. `the column picker of {element} should offer only
  {string} after typing {string}` (open, type into its search, read the rows of its grid, dismiss).
  Then:

  ```gherkin
  Scenario: The add-filter column picker narrows its list as the user types (public#605)
    When user hovers over filter panel
    Then the column picker of add filter of filter panel should offer only "HEIGHT" after typing "HEI"
    And the filter panel should have 0 filters
    And no errors should have been logged
  ```

- **A column made multi-value from Column Properties counts each value in its card**
  ([public#3874](https://github.com/datagrok-ai/public/issues/3874)). The tag goes in through the
  Column Properties dialog's Tags pane, and no step adds a tag row there; the card's counts are
  `count <label> of <caption>` areas only, not readings, so "shows a count" has nothing to compare.
  Needs a step for a tag row of the dialog (`user adds the tag "<name>" = "<value>" in the Column
  Properties dialog`) and a `counts of <caption>` reading next to `categories of` in
  `filters_status.dart` (`_addGridFilter` already walks `filteredCountCol`). Then (spgi-100, Series):

  ```gherkin
  Scenario: A column made multi-value from Column Properties counts each value in its card (public#3874)
    When user picks "Column Properties..." from the context menu of the "header Series" area of grid
    And user adds the tag ".multi-value-separator" = "|" in the Column Properties dialog
    And user clicks on OK button in Column Properties dialog
    And user adds a card for "Series" to the filter panel
    Then the "type of Series" reading of filter panel should be "multi-value"
    And the "counts of Series" reading of filter panel should contain "Pyrrolidines: 12"
    And no errors should have been logged
  ```

  (the count is a placeholder — measure it in the run that first passes.)

- **The second card's search icon itself moves the focus to its search**
  ([public#2359](https://github.com/datagrok-ai/public/issues/2359)). Every `Icon` puts the focus back
  where it was on mouse down one `Timer.run` later (`icons.dart`, GROK-8211); Playwright's click
  releases within that tick, so the restore moves the caret back to the first card's search (6 of 6
  headless clicks; 6 of 6 clicks held 120 ms keep it). Needs a hand-length press on the library's
  icon click before `search of "DIS_POP" filter card should be focused` and `user types "RA" at the
  caret` can say it. Typing into the second search by name does not reach the bug, so no scenario
  claims the issue now.

- **The single-value row keeps its close and menu icons** (core#479, second half).
  `filter-panel-categorical.feature` claims the row (`card <caption>` area) after the cards and the
  card itself hidden; the row's icons are not parts the panel reports, so "still offers its close" is
  unclaimed. Needs `close of "<caption>" single-value row` — a kind for the row
  (`_singleValueRow`, `filters_core.dart`) with its close and menu icons as parts.

## Edge cases (now `@edge-cases` scenarios of the filter-panel features)

Added with the filter-panel edge cases (edge-case round of 2026-10-26, localhost), now `@edge-cases` scenarios of their own, after the master scenarios of the filter-panel features whose topic they are; each known failure is a scenario of its own: the ones of a journey feature are in `filter-panel-known-failures.feature` while their tickets are open (core#1174 from `combined-boolean-filter.feature`, core#1223 from `filter-panel-persistence.feature` and `hierarchical-filter.feature`, core#1224 from `filter-panel-saved-states.feature`), the others stay in their feature: `combined-boolean-filter.feature`, `expression-filter-all-columns.feature`, `filter-panel-categorical.feature`, `filter-panel-entry-points.feature`, `filter-panel-expression-and-text.feature`, `filter-panel-linked-tables.feature`, `filter-panel-persistence.feature`, `filter-panel-regressions.feature`, `filter-panel-saved-states.feature` and `hierarchical-filter.feature`. Only the scenarios were written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

#### Steps requested

1. **The validation state of a range card's min / max fields.** Min equal to max, or above it, makes
   the card's input validator say "Max should be greater than min" (`histogram_core.dart`), and `abc`
   makes it say nothing. No step reads a card field's invalid state or its message. Wanted: the card's
   min/max inputs reporting `aria-invalid` and the message, then `"max" field of the "AGE" filter card
   should be invalid`.
2. **The free-text input's invalid state.** A free-text rule that does not parse is refused on Enter
   and the input gets `d4-invalid` with an error tooltip (`text_filter.dart`). The scenario claims only
   that no rule is added and the rows stay. Wanted: the state of the Expression card's search input.

#### Waiting for core#1174

3. **Expression and Combined Boolean cards, a column they use removed** while the panel has no card of
   that column (with the column's own card present, its removal refilters the table and the rows come
   back): the table stays filtered. Four known failures: a Combined Boolean flag in
   `filter-panel-known-failures.feature`, and in `expression-filter-all-columns.feature` an Expression rule
   on one column (HEIGHT is null, 128 rows), a Numeric Columns rule (884 rows, 131 once HEIGHT is gone)
   and an All Columns comparison.

#### Decided (2026-10-28)

5. **A range card's min typed above its max**: the present handling is right (the fields are flagged
   with a message); claimed once a step reads the flag (item 1).
7. **An expression rule naming a removed column**: covered by the rule "a filter goes with its
   column".
10. **The text card's `!malt` and `^hop`**: right as it is. At Fuzzyness 0 the term is literal text
    (a plain "contains"); above 0 the card searches with Fuse.js extended search, where `!` is "does not
    contain" and `^` is "starts with" (`text_filter.dart`).

#### API or Dart tests to write

12. A card on a numeric column whose values are all empty (Dart regTest with log capture).
13. A number typed with a decimal comma under a comma-decimal locale (`NumericMatcher`, Dart unit test).

### Won't do

- **A layout naming a column the table no longer has** on a card of one column (`createFilter(null)`):
  not reproduced on localhost (the other cards come back, no error), and a removed calculated column
  is recreated by the layout itself. A hierarchical card naming a removed level column is another
  matter: the card should come back without that level (core#1223), two known failures in
  `filter-panel-known-failures.feature` (from `filter-panel-persistence.feature` and `hierarchical-filter.feature`).
- **Card summary at 19 / 20 categories**: a cosmetic threshold of the summary text.
- **Already covered**: all rows filtered out (`filter-panel-categorical.feature`), Filter out missing
  values (`filter-panel-cloned-view.feature`), a column with one value (`filter-panel-categorical.feature`).
