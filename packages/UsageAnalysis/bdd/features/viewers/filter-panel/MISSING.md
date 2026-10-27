# Filter panel — what the features cannot say yet

## Closed issues the regression round of 26 Oct 2026 could not claim

Written during the closed-issue regression round of 26 Oct 2026 (bucket b5, filtering), against the
panel's "Automation surface" (`core/client/d4/lib/src/viewers/filters/CLAUDE.md`). What could be said
went into `filter-panel-regressions.feature` (public#790, #2364, #1343, #2359, #1986, #3027, core#479,
public#2501, #2943) and `filter-panel-layout-regressions.feature` (public#3059, #1837, #1103, #791,
#377, #1464, #1889, #319). Chem's substructure-card candidates of the same bucket are in
`Chem/bdd/features/filters/`.

### Regressions found 26 Oct 2026

- **An invalid viewer Filter formula logs a RangeError** ([public#344](https://github.com/datagrok-ai/public/issues/344)).
  The page keeps working and the formula is ignored (the scatter plot keeps showing all 1000 rows), but
  every time the formula is set the console gets an unhandled
  `RangeError (index): Index out of range: no indices are valid: 0`, thrown at
  `packages/d4/src/viewer_base/data_frame_viewer.dart 474:26 DataFrameViewer.createCombinedRowFilter`:
  `getNewColumn(dataFrame, filterFormula, …)` answers an empty list for a formula that does not parse,
  and `columns[0]` is read before the `column is BoolColumn` check that already handles a formula
  that is not boolean. Measured in three runs of the scenario below and once in a scratch Playwright
  probe on the local master stand (26 Oct 2026). Steps by hand: open demog, add a scatter plot, open
  its properties, set Filter to `${RACE} = "Black"` (one `=`) and press Enter; the console shows the
  RangeError. Expected (the issue): the filter is not applied or shows no data, and nothing breaks —
  the first half holds, the error is the regression. Then:

  ```gherkin
  Scenario: An invalid viewer Filter formula leaves the page working (public#344)
    Given user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    Then scatter plot viewer should show 1000 rows
    When user sets "Filter" property of scatter plot viewer to '${RACE} = "Black"'
    Then all rows should pass the filter
    And scatter plot viewer should show 1000 rows
    And no errors should have been logged
    When user sets "Filter" property of scatter plot viewer to '${RACE} == "Black"'
    Then scatter plot viewer should show 27 rows
    And all rows should pass the filter
    When user sets "Filter" property of scatter plot viewer to ""
    Then scatter plot viewer should show 1000 rows
    And no errors should have been logged
  ```

  It fails on the first `no errors should have been logged`. It goes back into
  `filter-panel-regressions.feature` (demog-1000, X and Y on AGE and WEIGHT, which have no blanks —
  the default axes include HEIGHT and draw 872) once the error is gone.

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
  ([public#2359](https://github.com/datagrok-ai/public/issues/2359), first half).
  `filter-panel-regressions.feature` opens both searches and types into the second one by name; it
  does not claim that the icon click alone hands it the caret. Every `Icon` puts the focus back where
  it was on mouse down, one `Timer.run` later (`core/client/d4/lib/src/widgets/icons.dart`, GROK-8211).
  A hand click holds the button longer than that, so the restore runs first and the search icon's
  click then focuses its search; Playwright's click releases within the same tick, the icon focuses
  the DIS_POP search and the restore moves the caret back to RACE 13 ms later (traced 27 Oct 2026:
  6 of 6 headless clicks lose it, 6 of 6 clicks held 120 ms keep it). Needs a press that lasts as long
  as a hand click — a `delay` on the library's icon click — before `search of "DIS_POP" filter card
  should be focused` and `user types "RA" at the caret` can say it.

- **The single-value row keeps its close and menu icons** (core#479, second half).
  `filter-panel-regressions.feature` claims the row (`card <caption>` area) after the cards and the
  card itself hidden; the row's icons are not parts the panel reports, so "still offers its close" is
  unclaimed. Needs `close of "<caption>" single-value row` — a kind for the row
  (`_singleValueRow`, `filters_core.dart`) with its close and menu icons as parts.
