@browse @realizes:views.browse
Feature: The quick filters of a Browse gallery
  The filter panel a gallery opens from its toolbar, with the quick filters on top. Translated from
  the manual cases Browse-Filter-02 and -03 (browse_manual_tests2.md section 11,
  playwright-public/browse/filter.test.ts), whose old specs clicked a tag and claimed nothing.

  A quick filter writes its query into the gallery's search ("created < 1w") and the gallery reloads
  on it; All clears the query and brings the counter back to what it was. Whether the filter
  narrows the list depends on the stand — a local stand whose demo files were copied this week keeps
  all of them, and the query is recursive, so it lists the files of the subfolders too — so the
  claim is the query and the way back, not a smaller count. Each scenario ends on All, the state the next one
  expects. Browse-Filter-01 (no nameless filter property, GROK-19691) and -04 (Users by group) are in
  users-view.feature as far as the Users filter panel goes; -06 (a search that finds nothing) there
  and in spaces-search.feature. Browse-Filter-05 (a filter and a search together, GROK-19690) is not
  written: the gallery search is fuzzy, so the number a search keeps depends on every login on the
  stand, and no count it leaves is one a claim could state.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: Created recently puts its query in the search and All takes it back
    Given Files tree node inside browse tree is expanded
    When user clicks on Files---Demo tree node inside browse tree
    Then the "Demo" view should be current
    When user clicks on "Toggle filters" icon
    And user remembers the gallery counter
    And user clicks on "Created recently" tag
    Then gallery search should have value "created < 1w"
    When user clicks on "All" tag
    Then gallery search should have value ""
    And the gallery counter should not be lower than remembered
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Dockers gallery offers Used by me and Created recently
    Given Platform tree node inside browse tree is expanded
    When user clicks on Platform---Dockers tree node inside browse tree
    Then the "Dockers" view should be current
    When user clicks on "Toggle filters" icon
    Then "Used by me" tag should be visible
    And "Created recently" tag should be visible
    When user remembers the gallery counter
    And user clicks on "Used by me" tag
    Then gallery search should have value "usedBy = @current"
    When user clicks on "All" tag
    Then gallery search should have value ""
    And the gallery counter should not be lower than remembered
    And no errors should have been logged
    And no error or warning balloon should have been shown
