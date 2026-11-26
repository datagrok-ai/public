@guide @tutorial
Feature: How-to · Basics
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A column header's menu formats and colors the column, on demog:
  Format shows every number format with the first value written in it (HEIGHT in scientific notation, WEIGHT with
  two digits after the point), the same item of a date column lists date formats (STARTED as dd.MM.yyyy), Color
  Coding > Linear shades AGE from low to high, and Color Coding > Categorical gives each RACE its own color. These
  items are radio choices, which leave the menu open: Escape closes it after each pick so the column shows. demog is
  opened off camera; the formats and colors stay in the table only, nothing is made on the server.

  Scenario: Format and color-code columns
    Number and date formats from a column's menu, then linear and categorical color coding.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    # caption: Right-click a header, then Format: each choice shows a value in it
    When user picks "Format > scientific (1.60E2)" from the context menu of the "header HEIGHT" area of grid
    And user presses Escape
    Then context menu should be absent
    And "HEIGHT" column should have tag "format" equal to "scientific"
    # caption: At most two digits after the point
    When user picks "Format > max two digits after comma (73.2)" from the context menu of the "header WEIGHT" area of grid
    And user presses Escape
    Then context menu should be absent
    And "WEIGHT" column should have tag "format" equal to "max two digits after comma"
    # caption: A date column offers date formats
    When user picks "Format > dd.MM.yyyy (11.09.1990)" from the context menu of the "header STARTED" area of grid
    And user presses Escape
    Then context menu should be absent
    And "STARTED" column should have tag "format" equal to "dd.MM.yyyy"
    # caption: Color Coding > Linear: shades from low to high
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user presses Escape
    Then context menu should be absent
    And "AGE" column should be color-coded linearly
    # caption: Color Coding > Categorical: a color for each category
    When user picks "Color Coding > Categorical" from the context menu of the "header RACE" area of grid
    And user presses Escape
    Then context menu should be absent
    And "RACE" column should be color-coded categorically
    And no errors should have been logged
