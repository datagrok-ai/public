@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). A walk down the table view's Data menu with one line per
  item, every one of them running in the browser: Aggregate Rows, Unpivot, the four two-table commands (Join, Link,
  Append, Compare: greyed out while only one table is open), Batch Edit, Split Column, Split Column by RegExp,
  Anonymize, Categorize and Compare Columns (PowerPack). Curves at the bottom is left out: its Calculate MSR runs
  a Python script. Nothing is picked; demog is opened off camera and nothing is made on the server.

  Scenario: What is in the Data menu
    The Data menu reshapes, combines and cleans the table you are looking at.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "Data: commands that reshape and combine tables" at the top
    # caption: Aggregate Rows: summarize groups of rows
    And user points at "Data > Aggregate Rows..." in the top menu
    # caption: Unpivot: turn columns into rows
    And user points at "Data > Unpivot..." in the top menu
    # caption: Join, Link, Append, Compare: these need a second open table
    And user points at "Data > Join Tables..." in the top menu
    Then "Join Tables..." menu item should be disabled
    # caption: Batch Edit: one value into many rows at once
    When user points at "Data > Batch Edit..." in the top menu
    # caption: Split Column: text into several columns, by a separator
    And user points at "Data > Split Column..." in the top menu
    # caption: Split Column by RegExp: pull out parts that match a pattern
    And user points at "Data > Split Column by RegExp..." in the top menu
    # caption: Anonymize: replace values with made-up ones
    And user points at "Data > Anonymize..." in the top menu
    # caption: Categorize: bin numbers into named ranges
    And user points at "Data > Categorize..." in the top menu
    # caption: Compare Columns: row-by-row differences of two columns
    And user points at "Data > Compare Columns..." in the top menu
    Then "Compare Columns..." menu item should be enabled
    And no errors should have been logged
