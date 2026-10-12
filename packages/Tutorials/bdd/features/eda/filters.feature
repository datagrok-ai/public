@tutorials @serial @realizes:tutorials.filters
Feature: The Filters tutorial
  Walks Exploratory Data Analysis > Filters from its card to the end: categorical filtering by a
  click and by the cursor keys, the category set inverted, rows selected from a count, two filters
  combined, the numeric card's bins hovered and selected, a range typed into Min / max, the filter
  saved as a column and as a named state, and the state restored. Each step is claimed as ticked
  and as done — the rows that pass, the categories each card holds, the rows selected, the new
  column, the restored range.
  Translated from playwright-tests/e2e/tutorials/filters.test.ts, which found the categories by a
  row pitch measured in pixels at one screen size. The filter panel reports its categories, and the
  numeric card's histogram is a viewer that reports its bins, so neither needs a signal of its own.

  Fixed in the tutorial for this translation: the row-count and bin-selection steps completed on
  any selection event; the indicator and reset hints were captured when their step began.
  The named filter state lives in the browser's storage and is removed before and after the walk.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "recentViewerSettings" user settings are put back at feature end
    And the "Filters" tutorial is not completed yet
    And no saved filter state "AGE: [40,60]" is kept, now or when the feature ends
    And the Tutorials app is open

  Scenario: A learner completes the Filters tutorial
    When user starts the "Filters" tutorial
    Then the tutorial progress should be 1 of 17
    When user clicks on scatter-plot icon in toolbox
    Then the tutorial step "Open scatter plot" should be done
    When user clicks on histogram icon in toolbox
    Then the tutorial step "Open histogram" should be done
    Then filter icon in toolbox should be hinted
    When user clicks on filter icon in toolbox
    Then the tutorial step "Open filters" should be done
    And filter panel should be visible

    When user clicks on the "category AS of DIS_POP" area of filter panel
    Then the tutorial step "Click on the \"AS\" label within the \"DIS_POP\" filter" should be done
    And the filter should pass exactly the rows where "DIS_POP" is "AS"
    When user presses ArrowDown in "DIS_POP" filter card
    Then the "selected categories of DIS_POP" reading of filter panel should be "Indigestion"
    When user presses ArrowDown in "DIS_POP" filter card
    Then the "selected categories of DIS_POP" reading of filter panel should be "PsA"
    When user presses ArrowDown in "DIS_POP" filter card
    Then the "selected categories of DIS_POP" reading of filter panel should be "Psoriasis"
    When user presses ArrowDown in "DIS_POP" filter card
    Then the "selected categories of DIS_POP" reading of filter panel should be "RA"
    And the tutorial step "Filter the dataset by the most common disease" should be done
    And 2550 rows should pass the filter

    When user picks "Invert all" from the indicator menu of the "DIS_POP" filter card
    Then the tutorial step "Invert the category set for the \"DIS_POP\" filter" should be done
    And the filter should pass exactly the rows where "DIS_POP" is one of "AS, Indigestion, PsA, Psoriasis, UC"

    When user clicks on the "count AS of DIS_POP" area of filter panel
    Then the tutorial step "Click on a non-empty row count" should be done
    And only rows where "DIS_POP" is "AS" should be selected

    When user clicks on the "category F of SEX" area of filter panel
    And user clicks on the "category Asian of RACE" area of filter panel
    # a name click applies "only this one" on a debounce that reads the card's row when it fires: claimed
    # before the next click on the same card moves that row
    Then the "selected categories of RACE" reading of filter panel should be "Asian"
    When user clicks on the "checkbox Black of RACE" area of filter panel
    Then the tutorial step "Filter the dataset to only females of Asian or Black origin" should be done
    And the "selected categories of SEX" reading of filter panel should be "F"
    And the "selected categories of RACE" reading of filter panel should be "Asian, Black"

    # the panel's header icons show while the panel is hovered
    When user hovers over filter panel
    Then "Reset filter" icon in filter panel should be hinted
    When user clicks on "Reset filter" icon in filter panel
    Then the tutorial step "Reset the filter" should be done
    And all rows should pass the filter

    When user takes a snapshot of scatter plot viewer
    And user hovers over the "bin 3" area of histogram viewer in filter panel
    Then the tutorial step "Hover over the histogram bins" should be done
    # the rows under the hovered bin are highlighted on the other viewers
    And scatter plot viewer should have repainted
    # the rows of the count clicked before are still selected: the bin's own selection is the change
    And the tutorial step "Select one of the histogram bins" should not be done yet
    When user remembers the "rows selected" reading of scatter plot viewer
    And user clicks on the "bin 3" area of histogram viewer in filter panel
    Then the tutorial step "Select one of the histogram bins" should be done
    And the "rows selected" reading of scatter plot viewer should not be as remembered
    And the "selected bin 3" area of histogram viewer in filter panel should be painted

    When user clicks on the "cell 4 of AGE" area of grid
    Then the tutorial step "Change the current row in the spreadsheet" should be done
    And row 4 should be current

    When user picks "Min / max" from the indicator menu of the "AGE" filter card
    And user enters "40" into the min field of the "AGE" filter card
    And user enters "60" into the max field of the "AGE" filter card
    Then the tutorial step "Find records for people aged 40 to 60" should be done
    # 2989 patients are 40 to 60, and a range filter keeps the one row with no AGE
    And 2990 rows should pass the filter
    And no rows where "AGE" is "39" should pass the filter
    And no rows where "AGE" is "61" should pass the filter

    When user picks "Filter to Column..." from the filter panel menu
    Then the tutorial step "Save the current filter as a column" should be done
    When user clicks on OK button in dialog
    Then the table should have a column "AGE: [40,60]"

    When user picks "Save or Apply > Save..." from the filter panel menu
    Then the tutorial step "Save the filter configuration as \"AGE: [40,60]\"" should be done
    When user clicks on OK button in dialog

    When user hovers over filter panel
    And user clicks on "Reset filter" icon in filter panel
    Then the tutorial step "Reset the filter" should be done 2 times
    And all rows should pass the filter
    When user picks "Save or Apply > AGE: [40,60]" from the filter panel menu
    Then the tutorial step "Restore the filter state" should be done
    And 2990 rows should pass the filter
    And no rows where "AGE" is "61" should pass the filter

    And the "Filters" tutorial should be completed
    And the tutorial should have listed 17 steps
    And the tutorial progress should be 17 of 17
    And no hint should be shown
    And no errors should have been logged
