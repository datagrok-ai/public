@guide @help:visualize/viewers
Feature: Reorder columns by dragging and hide the ones you do not need
  A guide: the answer to "how do I reorder columns by drag-and-drop, and how do I hide some
  columns?". A column moves when its header is dragged onto another header: it lands right after
  the column it was dropped on. A column is hidden from the header's context menu (Hide), or from
  the Order or Hide Columns... dialog of the same menu, where every column has a checkbox and the
  names can be dragged into a new order. Demo: demog, whose HEIGHT column is dragged onto AGE, then
  WEIGHT is hidden from the header menu and STARTED from the dialog.

  Scenario: Drag a column header to a new place, then hide two columns
    Given user is logged in
    And user opens demog dataset
    When user drags the "header HEIGHT" area of grid to the "header AGE" area
    Then the "column order" reading of grid should include the text "AGE, HEIGHT, SEX"
    When user picks "Hide" from the context menu of the "header WEIGHT" area of grid
    Then grid should not have a "header WEIGHT" area
    When user picks "Order or Hide Columns..." from the context menu of the "header AGE" area of grid
    Then Order or Hide Columns dialog should be visible
    When user toggles the "STARTED" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then grid should not have a "header STARTED" area
    And the "column order" reading of grid should not include the text "WEIGHT"
