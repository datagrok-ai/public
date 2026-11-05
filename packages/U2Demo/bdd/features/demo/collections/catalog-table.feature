@demo @realizes:u2.catalog-table
Feature: Catalog table
  A catalog table is found by its name, its rows by their key or first cell, its headers as column headers by
  their text. A double click on a header sorts up, down and back to the default order; the group icon's menu
  groups the rows into sections that count them and fold; the eye's menu adds a column at the end; a row's actions
  are also its context menu. The table reports its sort, grouping, columns and groups.

  Background:
    Given the catalog table choices kept under "u2demo-campaigns" are cleared now and put back at feature end
    And user opens the "Catalog table" demo page

  Scenario: A double click on a header sorts up, down, then back to the default order
    Then campaigns catalog table should have 50 items in total
    And campaigns catalog table should report "sort" as ""
    When user double-clicks "Molecules" column header in campaigns catalog table
    Then campaigns catalog table should report "sort" as "Molecules asc"
    And campaigns catalog table should list the items "KHL-4, KHL-14, SKP-9" in this order
    When user double-clicks "Molecules" column header in campaigns catalog table
    Then campaigns catalog table should report "sort" as "Molecules desc"
    And campaigns catalog table should list the items "SKP-9, KHL-14, KHL-4" in this order
    And value of sort readout should have text "Molecules descending"
    When user double-clicks "Molecules" column header in campaigns catalog table
    Then campaigns catalog table should report "sort" as ""
    And value of sort readout should have text "(last modified)"

  Scenario: The group menu groups the rows into counted sections that fold
    When user clicks on campaigns-group icon button
    And user clicks on "Template" menu item
    Then campaigns catalog table should report "groupBy" as "Template"
    And campaigns catalog table should report "groups" as "Fragment expansion (10), Quick sketchpad (17), Kinase hit-to-lead (14), GPCR series (9)"
    And item "Quick sketchpad" of campaigns catalog table should be expanded
    When user collapses "Quick sketchpad" table row in campaigns catalog table
    Then item "Quick sketchpad" of campaigns catalog table should be collapsed
    And campaigns catalog table should report "collapsed" as "SKP"
    And campaigns catalog table should report "rows" as "33"
    When user expands "Quick sketchpad" table row in campaigns catalog table
    Then campaigns catalog table should report "rows" as "50"

  Scenario: The row under the pointer is marked as the platform grid marks it, lighter than a section header
    When user clicks on campaigns-group icon button
    And user clicks on "Template" menu item
    And user hovers over "SKP-9" table row in campaigns catalog table
    Then the row under the pointer in campaigns catalog table should be marked as the platform grid marks it

  Scenario: The eye's menu adds a column at the end
    When user clicks on campaigns-columns icon button
    And user clicks on "Last Modified" menu item
    Then campaigns catalog table should report "columns" as "Name,Code,Created,Author,Last Modified by,Molecules,Status,Last Modified"
    And "Last Modified" column header in campaigns catalog table should be visible

  Scenario: A row's actions show once its access is known, and are its context menu; a view-only row has none
    Then "Delete campaign" icon button in "SKP-9" table row in campaigns catalog table should be visible
    And "Delete campaign" icon button in "FRAG-9" table row in campaigns catalog table should be absent
    When user right-clicks on "SKP-9" table row in campaigns catalog table
    Then "Delete campaign" menu item should be visible
    And "Manage campaign permissions" menu item should be visible
