@viewers @realizes:charts.viewer.sunburst
Feature: Sunburst segments, clicks, filtering, empty values and the project round trip
  The Sunburst draws a categorical hierarchy as nested rings, one segment per distinct value under
  its parent. Every claim is the viewer's own report: a segment is the `segment <path>` area (its
  path from the centre joined by " | ", `segment F | Caucasian`), and its row count the
  `rows of segment <path>` reading; `segments`, `hierarchy columns`, `on click` and `include nulls`
  are readings too. The viewer is added from the ribbon's Add viewer gallery and its hierarchy is
  picked in the Select columns dialog of the Hierarchy property, as a user does it. The dialog keeps
  a column that is already in the hierarchy in its place and appends the ones checked after it, so
  RACE with SEX checked again is RACE, SEX (12 segments). The md's last step, SEX dragged back above
  RACE, has no gesture: the dialog's rows cannot be dragged (a drag starts a column drag-out), and the
  checked columns come back in the order of its rows, the current hierarchy's on top.
  The table is demog-1000 under a run's own name, which the layout the Layouts pane saves takes.
  Counts on demog-1000: SEX F 553 / M 447; F | Caucasian 480, F | Other 48, F | Black 18,
  F | Asian 7, M | Caucasian 416, M | Other 14, M | Black 9, M | Asian 8.
  Translated from the TestTrack case Charts/sunburst. In ae.csv AESTDTC is a string column (its
  values are not parsed as dates), so the date column claimed missing is AEENDTC, as in the md. An
  empty value's sector is `segment S_PART | (empty)` (github-2992).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset keeping the first 1000 rows as "Sunburst-{time}"
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Sunburst" card in "Add Viewer" dialog
    Then sunburst viewer should be visible
    And "Add Viewer" dialog should be absent
    When user picks "Properties..." from the context menu of the "view" area of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user toggles the "SEX" column in the column list of "Select columns..." dialog
    And user toggles the "RACE" column in the column list of "Select columns..." dialog
    Then "2 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Select columns..." dialog should be absent
    And the "hierarchy columns" reading of sunburst viewer should be "SEX, RACE"
    And the "segments" reading of sunburst viewer should be 10

  Scenario: The hierarchy draws one segment per value with its row count
    Then the "segment names" reading of sunburst viewer should contain "F"
    And the "segment names" reading of sunburst viewer should contain "M"
    And the "segment names" reading of sunburst viewer should contain "F | Caucasian"
    And the "segment names" reading of sunburst viewer should contain "M | Asian"
    And the "rows of segment F" reading of sunburst viewer should be 553
    And the "rows of segment M" reading of sunburst viewer should be 447
    And the "rows of segment F | Caucasian" reading of sunburst viewer should be 480
    And the "rows of segment M | Asian" reading of sunburst viewer should be 8
    When user hovers over the "segment F | Black" area of sunburst viewer
    Then tooltip should contain text "18"
    And tooltip should contain text "Black"
    When user moves the pointer away from sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "SEX" column in the column list of "Select columns..." dialog
    Then "1 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the "hierarchy columns" reading of sunburst viewer should be "RACE"
    And the "segments" reading of sunburst viewer should be 4
    And the "rows of segment Caucasian" reading of sunburst viewer should be 896
    When user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "SEX" column in the column list of "Select columns..." dialog
    Then "2 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the "hierarchy columns" reading of sunburst viewer should be "RACE, SEX"
    And the "segments" reading of sunburst viewer should be 12
    And the "rows of segment Caucasian | F" reading of sunburst viewer should be 480
    And sunburst viewer should report no error
    And no errors should have been logged

  Scenario: Only categorical columns can build the hierarchy (github-2954, GROK-18010)
    ae.csv, opened from Browse, has 35 columns; the dialog lists the string and boolean ones only:
    the date column AEENDTC and the numeric AESEQ and AESTDY are not in its list, while a search for
    AESEV finds it. Then a Sunburst on spgi-100 draws Core and R101 (GROK-18010).
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Charts tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---Charts---ae.csv tree node inside browse tree
    Then the "ae" view should be current
    When user clicks on "Add viewer" icon
    And user clicks on first "Sunburst" card in "Add Viewer" dialog
    Then sunburst viewer should be bound to table "ae"
    When user picks "Properties..." from the context menu of the "view" area of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    And the "rows" reading of grid in "Select columns..." dialog should be 30
    And the column list of "Select columns..." dialog should not list "AEENDTC"
    And the column list of "Select columns..." dialog should not list "AESEQ"
    And the column list of "Select columns..." dialog should not list "AESTDY"
    When user types "AESEV" into "Search" input in "Select columns..." dialog
    Then the column list of "Select columns..." dialog should be exactly "AESEV"
    When user clicks on "CANCEL" button in "Select columns..." dialog
    Then "Select columns..." dialog should be absent
    Given user opens spgi dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Sunburst" card in "Add Viewer" dialog
    Then sunburst viewer should be bound to table "spgi-100"
    When user picks "Properties..." from the context menu of the "view" area of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user types "Core" into "Search" input in "Select columns..." dialog
    And user toggles the "Core" column in the column list of "Select columns..." dialog
    And user types "R101" into "Search" input in "Select columns..." dialog
    And user toggles the "R101" column in the column list of "Select columns..." dialog
    Then "2 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the "hierarchy columns" reading of sunburst viewer should be "Core, R101"
    And the "segments" reading of sunburst viewer should be 14
    And the "segment names" reading of sunburst viewer should contain "(empty)"
    And sunburst viewer should report no error
    And no errors should have been logged

  Scenario: Click, Ctrl+Click, Shift+Click and Ctrl+Shift+Click select segments
    When user presses Escape in grid
    Then no rows should be selected
    When user clicks on the "segment F" area of sunburst viewer
    Then 553 rows should be selected
    And only rows where "SEX" is "F" should be selected
    When user clicks on the "segment M | Asian" area of sunburst viewer holding Control
    Then 561 rows should be selected
    When user clicks on the "segment M | Asian" area of sunburst viewer holding Control+Shift
    Then 553 rows should be selected
    And only rows where "SEX" is "F" should be selected
    When user clicks on the "segment M | Other" area of sunburst viewer holding Shift
    Then 567 rows should be selected
    When user clicks on the "segment F | Black" area of sunburst viewer
    Then 18 rows should be selected
    And no rows where "SEX" is "M" should be selected
    And no rows where "RACE" is "Caucasian" should be selected
    And no rows where "RACE" is "Other" should be selected
    And no rows where "RACE" is "Asian" should be selected
    And no errors should have been logged

  Scenario: On Click = Filter; a double-click on empty space and Reset View clear it (github-2994, github-3097)
    Given "Misc" category in context panel is expanded
    When user presses Escape in grid
    And user selects "Filter" in "On Click" property in context panel
    Then the "on click" reading of sunburst viewer should be "Filter"
    And "Row Source" property in context panel should contain text "All"
    When user clicks on the "segment F | Black" area of sunburst viewer
    Then 18 rows should pass the filter
    And the "segments" reading of sunburst viewer should be 10
    When user double-clicks on empty plot space of sunburst viewer
    Then 1000 rows should pass the filter
    When user clicks on the "segment M" area of sunburst viewer
    Then 447 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "M"
    When user picks "Reset View" from the context menu of the "view" area of sunburst viewer
    Then 1000 rows should pass the filter
    When user selects "Select" in "On Click" property in context panel
    Then the "on click" reading of sunburst viewer should be "Select"
    And "Row Source" property in context panel should contain text "Filtered"
    And no errors should have been logged

  Scenario: Empty values — Include Nulls, and a click on the empty segment (github-2992)
    spgi-100: Stereo Category R_ONE 36, S_ACHIR 34, S_UNKN 18, S_PART 10, S_ABS 2; Series is empty
    in 3 S_PART rows, 3 R_ONE rows and 1 S_UNKN row.
    Given user opens spgi dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Sunburst" card in "Add Viewer" dialog
    Then sunburst viewer should be bound to table "spgi-100"
    When user picks "Properties..." from the context menu of the "view" area of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user types "Stereo Category" into "Search" input in "Select columns..." dialog
    And user toggles the "Stereo Category" column in the column list of "Select columns..." dialog
    And user types "Series" into "Search" input in "Select columns..." dialog
    And user toggles the "Series" column in the column list of "Select columns..." dialog
    Then "2 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the "hierarchy columns" reading of sunburst viewer should be "Stereo Category, Series"
    And the "include nulls" reading of sunburst viewer should be "true"
    And the "segments" reading of sunburst viewer should be 17
    And the "rows of segment S_PART" reading of sunburst viewer should be 10
    And the "rows of segment R_ONE" reading of sunburst viewer should be 36
    And the "rows of segment S_PART | (empty)" reading of sunburst viewer should be 3
    When user presses Escape in grid
    And user clicks on the "segment S_PART | (empty)" area of sunburst viewer
    Then 3 rows should be selected
    And no rows where "Stereo Category" is "R_ONE" should be selected
    And no rows where "Stereo Category" is "S_ACHIR" should be selected
    And no rows where "Stereo Category" is "S_UNKN" should be selected
    And no rows where "Stereo Category" is "S_ABS" should be selected
    When user presses Escape in grid
    Given "Value" category in context panel is expanded
    When user unchecks "Include Nulls" property in context panel
    Then the "include nulls" reading of sunburst viewer should be "false"
    And the "segments" reading of sunburst viewer should be 14
    And the "rows of segment S_PART" reading of sunburst viewer should be 7
    And the "rows of segment R_ONE" reading of sunburst viewer should be 33
    And the "rows of segment S_UNKN" reading of sunburst viewer should be 17
    And sunburst viewer should not have a "segment S_PART | (empty)" area
    When user checks "Include Nulls" property in context panel
    Then the "segments" reading of sunburst viewer should be 17
    And the "rows of segment S_PART" reading of sunburst viewer should be 10
    And no errors should have been logged

  Scenario: The viewer's filter and the Filter Panel combine (GROK-15543)
    Given "Misc" category in context panel is expanded
    When user selects "Filter" in "On Click" property in context panel
    And user clicks on the "segment F" area of sunburst viewer
    Then 553 rows should pass the filter
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user clicks on the "category Caucasian of RACE" area of filter panel
    Then 480 rows should pass the filter
    And no rows where "SEX" is "M" should pass the filter
    When user hovers over "RACE" filter card
    And user clicks on close of "RACE" filter card
    Then 553 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "F"
    When user double-clicks on empty plot space of sunburst viewer
    Then 1000 rows should pass the filter
    And no errors should have been logged

  Scenario: A scatter plot's legend color picker does not strip the Sunburst's colors (github-3412)
    Given user opens spgi dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Sunburst" card in "Add Viewer" dialog
    Then sunburst viewer should be bound to table "spgi-100"
    When user picks "Properties..." from the context menu of the "view" area of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user types "Stereo Category" into "Search" input in "Select columns..." dialog
    And user toggles the "Stereo Category" column in the column list of "Select columns..." dialog
    Then "1 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the "hierarchy columns" reading of sunburst viewer should be "Stereo Category"
    And the "segment R_ONE" area of sunburst viewer should contain the color "#1F77B4"
    And the "segment S_UNKN" area of sunburst viewer should contain the color "#9467BD"
    And the "segment S_UNKN" area of sunburst viewer should not contain the color "#1F77B4"
    When user clicks on "Add viewer" icon
    And user clicks on first "Scatter plot" card in "Add Viewer" dialog
    Then scatter plot viewer should be visible
    When user clicks on grid
    And user clicks on settings icon of scatter plot viewer
    Given "Color" category in context panel is expanded
    When user selects "Stereo Category" in "Color" property in context panel
    Then "Color" property of scatter plot viewer should be "Stereo Category"
    Then the legend of scatter plot viewer should list 5 items
    When user hovers over "R_ONE" legend item in legend of scatter plot viewer
    And user clicks on last color picker icon
    Then "R_ONE" dialog should be visible
    When user clicks on CANCEL button in "R_ONE" dialog
    Then "R_ONE" dialog should be absent
    When user hovers over the "segment S_ACHIR" area of sunburst viewer
    And user moves the pointer away from sunburst viewer
    Then the "segment R_ONE" area of sunburst viewer should contain the color "#1F77B4"
    And the "segment S_UNKN" area of sunburst viewer should contain the color "#9467BD"
    And the "segment S_UNKN" area of sunburst viewer should not contain the color "#1F77B4"
    When user clicks on close icon of sunburst viewer
    Then the open tableview should have 0 sunburst viewers
    When user clicks on "Add viewer" icon
    And user clicks on first "Sunburst" card in "Add Viewer" dialog
    Then sunburst viewer should be bound to table "spgi-100"
    When user picks "Properties..." from the context menu of the "view" area of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user types "Stereo Category" into "Search" input in "Select columns..." dialog
    And user toggles the "Stereo Category" column in the column list of "Select columns..." dialog
    Then "1 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the "hierarchy columns" reading of sunburst viewer should be "Stereo Category"
    And the "segment R_ONE" area of sunburst viewer should contain the color "#1F77B4"
    And the "segment S_UNKN" area of sunburst viewer should contain the color "#9467BD"
    And the "segment S_UNKN" area of sunburst viewer should not contain the color "#1F77B4"
    And no errors should have been logged

  Scenario: Hierarchy and On Click survive a project save and reopen
    Given no project named "SunburstRoundTrip{time}" is on the server
    And "Misc" category in context panel is expanded
    When user selects "Filter" in "On Click" property in context panel
    Then the "on click" reading of sunburst viewer should be "Filter"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "SunburstRoundTrip{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share SunburstRoundTrip{time}" dialog should be visible
    When user clicks on CANCEL button in "Share SunburstRoundTrip{time}" dialog
    Then the "Share SunburstRoundTrip{time}" dialog should close
    And 1 project named "SunburstRoundTrip{time}" should be on the server
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "SunburstRoundTrip{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on SunburstRoundTrip{time} gallery card
    Then the task bar should have finished "Opening project"
    And the "Sunburst-{time}" view should be current
    And sunburst viewer should be bound to table "Sunburst-{time}"
    And the "hierarchy columns" reading of sunburst viewer should be "SEX, RACE"
    And the "on click" reading of sunburst viewer should be "Filter"
    And the "segments" reading of sunburst viewer should be 10
    And no errors should have been logged

  Scenario: A layout saved from the Layouts pane restores the hierarchy
    Given the layouts named "Sunburst-{time}" are deleted when the feature ends
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    When user clicks on Save button in layouts pane
    Then "Sunburst-{time}" layout card should be visible
    When user picks "Properties..." from the context menu of the "view" area of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "SEX" column in the column list of "Select columns..." dialog
    Then "1 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the "hierarchy columns" reading of sunburst viewer should be "RACE"
    And the "segments" reading of sunburst viewer should be 4
    When user clicks on "Sunburst-{time}" layout card
    Then the "hierarchy columns" reading of sunburst viewer should be "SEX, RACE"
    And the "segments" reading of sunburst viewer should be 10
    And no errors should have been logged
