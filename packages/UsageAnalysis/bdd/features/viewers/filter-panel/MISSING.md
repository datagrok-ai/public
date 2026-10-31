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
