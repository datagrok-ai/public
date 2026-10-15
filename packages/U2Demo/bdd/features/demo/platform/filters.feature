@demo @realizes:u2.filter-builder @realizes:u2.filter-query-input @realizes:u2.tabs
Feature: Filter builder
  The builder and the query input edit one filter: a condition built in the builder shows in the
  query at once and filters the rows on Apply, and the query outlives a visit to another tab. The
  domain tab counts the rows of the platform's users table, which every stand has.

  Background:
    Given user opens the "Filters" demo page

  Scenario: A condition built in the builder filters the dataframe on Apply
    When user clicks on "Add condition" button
    And user selects "city" in property of filter row
    And user selects "equals" in operator of filter row
    And user enters "Basel" into value of filter row
    Then Query input should have value 'city = "Basel"'
    When user clicks on Apply button
    Then value of rows readout should have text "2 of 6"
    When user clicks on "Entity type" tab
    Then tabs should show the "Entity type" tab
    When user clicks on DataFrame tab
    Then Query input should have value 'city = "Basel"'
    And value of rows readout should have text "2 of 6"

  Scenario: A condition on the users table counts its rows
    When user clicks on "Domain table" tab
    And user clicks on Open button
    And user clicks on "Add condition" button
    And user selects "Login" in property of filter row
    And user selects "equals" in operator of filter row
    And user enters "admin" into value of filter row
    Then Query input should have value 'login = "admin"'
    When user clicks on Count button
    Then value of rows readout should have text "1"
