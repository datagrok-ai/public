@journey
Feature: Apps > Plates: the search form follows its template, a template opens in its editor, Create warns with no plate
  Plates fills its Apps node with Create, Search plates, Search wells, Search analyses and Templates
  (`src/plates/plates-app.ts`); the node itself opens the list of plates, with their properties as columns. A search
  view builds its property inputs from the template picked
  (`plates-search-view.ts`); a numeric one marks a criterion that is not one as invalid, and an emptied one as valid
  again. EDIT opens a template in its editor, and CREATE with no plate imported only
  warns. Each scenario opens its view from the tree again: the app's views hide the Browse panel. Search
  plates and Search analyses opening with their inputs: browse/browse-apps-matrix.feature. Left out
  (MISSING.md): the plate list's rows and the search results, CREATE and template writes, readers and the search
  syntax (package tests).

  Background:
    Given user is logged in
    And the "Plates" package is installed
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Plates tree node inside browse tree is expanded

  Scenario: Apps > Plates opens the list of plates with their properties as columns
    When user clicks on Apps---Plates tree node inside browse tree
    Then the "Plates" view should be current
    And the table should have a column "Imaging device"
    And no errors should have been logged

  Scenario: Search Wells checks a numeric criterion and swaps its well inputs for the template picked
    Given the browse panel is open
    And Apps---Plates tree node inside browse tree is expanded
    When user clicks on Apps---Plates---Search-wells tree node inside browse tree
    Then the "Search Wells" view should be current
    And "Well cell count" input should be visible
    When user enters "abc" into "Well cell count" input
    Then "Well cell count" input should be invalid
    When user clears "Well cell count" input
    Then "Well cell count" input should be valid
    When user selects "Dose-response" in "Template" input
    Then "Concentration" input should be visible
    And "Well cell count" input should be hidden
    And no errors should have been logged

  Scenario: EDIT opens the template in its editor
    Given the browse panel is open
    And Apps---Plates tree node inside browse tree is expanded
    When user clicks on Apps---Plates---Templates tree node inside browse tree
    Then the "Templates" view should be current
    When user clicks on first "EDIT" button
    Then the "Templates / Cell counting" view should be current
    And "Name" input should have value "Cell counting"
    And no errors should have been logged

  Scenario: CREATE with no plate imported only warns
    Given the browse panel is open
    And Apps---Plates tree node inside browse tree is expanded
    When user clicks on Apps---Plates---Create tree node inside browse tree
    Then the "Create Plate" view should be current
    When user clicks on "CREATE" button
    Then a warning balloon containing "No active plate is selected." should have been shown
    And no errors should have been logged
