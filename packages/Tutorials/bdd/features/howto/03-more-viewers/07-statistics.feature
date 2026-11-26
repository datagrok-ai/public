@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The Statistics viewer from the toolbox on demog (the toolbox
  then hidden with its sidebar icon and the viewer dragged by its title bar to the bottom, on camera, so every
  statistic fits): one row per column with its counts and statistics, the nulls of
  HEIGHT hovered; the row menu adds a statistic (sum) and a histogram split by SEX; then the settings icon: a
  Filter formula (Data) that summarizes only the 878 rows with AGE over 60. Open as table is left out: it writes
  the STARTED dates as raw numbers (see the report). demog is opened off camera; nothing is made on the server.

  Scenario: Summarize every column with the Statistics viewer
    The Statistics viewer gives one row per column: how many values, how many missing, and the usual statistics.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "Statistics: a summary of every column at once"
    And user clicks on statistics icon on toolbox
    Then statistics viewer should be added to the open tableview
    # caption: Hide the toolbox for more room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: Drag its title bar to the bottom, so every statistic fits
    When user docks statistics viewer to the bottom edge of the view
    Then statistics viewer should be docked along the bottom edge of the view
    When the video says "HEIGHT is missing in 751 rows: see nulls" at the top
    And user hovers over the "cell nulls of HEIGHT" area of statistics viewer
    Then the "nulls of HEIGHT" reading of statistics viewer should be "751"
    # caption: Right-click a row to add a statistic
    When user picks "Statistics > sum" from the context menu of the "row AGE" area of statistics viewer
    Then the "stats" reading of statistics viewer should contain "sum"
    # caption: or a histogram of each column, split by a category
    When user picks "Histograms > SEX" from the context menu of the "row AGE" area of statistics viewer
    Then the "histogram columns" reading of statistics viewer should be 1
    When user closes the context menu
    # caption: The settings hold the rest
    And user hovers over statistics viewer
    And user clicks on settings icon of statistics viewer
    Then the context panel should show "Statistics"
    # caption: Data: a filter formula summarizes only some rows
    When user types "${AGE} > 60" into "Filter" property in context panel
    And user presses Enter
    Then the "rows shown" reading of statistics viewer should be 878
    And all rows should pass the filter
    And no errors should have been logged
