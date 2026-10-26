# Trellis plot — closed issues the regression round of 26 Oct 2026 could not claim

Kept apart from `MISSING.md` (the regression round of closed issues; the edge-case round keeps its own notes there).

From the closed-issues round of 26 Oct 2026
(`trellis-plot-regressions.feature` holds public#2144, public#3617 and public#1547):

- **A double-click on the trellis's whitespace no longer resets its selection or filter (core#932).**
  The removed handler (`root.onDoubleClick`, `_isWhitespace` in `trellis_plot_core.dart` before
  ec9f2efe25) counted a cell's own padding as whitespace, and that padding is also the cell's click
  target (the `cell <x> | <y>` corner band is the cell div). A plain double-click there is two cell
  clicks, which select the cell and then clear it (a click on the only selected cell clears), so the
  old reset and the new behaviour end in the same state. No point of the trellis is whitespace and
  not a cell (probed on demog-1000, SEX by RACE: every point of the charts grid is a cell, a label,
  a selector or a slider), so `user double-clicks on empty plot space` has nothing to aim at. Needs a
  double-click holding a modifier (Control toggles twice and keeps the selection; the old handler
  reset it), `user double-clicks on the {string} area of {widget} holding {key}`:

  ```gherkin
  Scenario: A double-click on a selected cell's whitespace keeps the selection (core#932)
    Given user sets "On Click" property of trellis plot viewer to "Select"
    When user clicks on the "cell F | Caucasian" area of trellis plot viewer
    And user clicks on the "cell M | Asian" area of trellis plot viewer holding Control
    Then 488 rows should be selected
    When user double-clicks on the "cell F | Black" area of trellis plot viewer holding Control
    Then 488 rows should be selected
    And the "selected cells" reading of trellis plot viewer should be "F | Caucasian, M | Asian"
    And no errors should have been logged
  ```

- **Inner box plot settings survive re-running the view's parameterised query (public#3692).** No step
  runs a query into a table view with its toolbox inputs, edits an input and clicks Refresh. Needs
  `user runs the query {string} into a table view` and `user enters {string} into the {string} input of
  the view's toolbox and refreshes` (the query fixture already exists: `a query … is on the server:`):

  ```gherkin
  Scenario: The inner box plot's settings survive a new query input (public#3692)
    Given a query "TrellisParam{time}" on "PostgresDemo" is on the server:
      """
      --input: int maxId = 10
      select * from employees where employee_id <= @maxId
      """
    When user runs the query "TrellisParam{time}" into a table view
    And user adds a trellis plot viewer
    And user picks "Box plot" in the viewer selector of trellis plot viewer
    And user sets "showStatistics" inner property of trellis plot viewer to "true"
    And user enters "20" into the "maxId" input of the view's toolbox and refreshes
    Then "showStatistics" inner property of trellis plot viewer should be "true"
    And no errors should have been logged
  ```

- **Show Empty Bars off gives every inner bar chart its own categories (public#3616).** The trellis
  reports a cell's picture (`cell signature`), not how many bars its inner bar chart drew. Needs a
  per-cell `bars of <xcat> | <ycat>` reading from the inner bar chart's own status (its `bars`
  reading), so the claim is the count of the categories present in that cell:

  ```gherkin
  Scenario: Show Empty Bars off draws only each cell's own categories (public#3616)
    Given user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX |
    When user picks "Bar chart" in the viewer selector of trellis plot viewer
    And user sets "showEmptyBars" inner property of trellis plot viewer to "false"
    Then the "bars of F" reading of trellis plot viewer should be <F categories, measured>
  ```
