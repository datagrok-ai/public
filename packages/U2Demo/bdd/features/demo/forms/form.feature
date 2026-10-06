@demo @realizes:u2.form
Feature: Forms
  Inputs resolve by label inside the form; a data table fills several of them at once.

  Background:
    Given user opens the "Form" demo page

  Scenario: Validity aggregates across the form and Submit reports the values
    Then value of "form validity" readout should have text "Required"
    And error of "Last name" input should have text "Required"
    When user types "Lovelace" into "Last name" input
    Then value of "form validity" readout should have text "valid"
    When user clicks on Submit button
    Then notification should contain text "\"first\":\"Ada\""

  Scenario: Filling and resetting
    When user fills in:
      | "First name" input | Grace  |
      | "Last name" input  | Hopper |
      | Age input          | 85     |
      | Role input         | Admin  |
      | Subscribe checkbox | no     |
    Then Age input should have value "85"
    And Role input should have value "Admin"
    And Subscribe checkbox should be unchecked
    When user clicks on Reset button
    Then "First name" input should have value "Ada"
    And "Last name" input should be empty
    And Role input should have value "Editor"
    And Subscribe checkbox should be checked
