@realizes:powerpack.search.power-pack @realizes:powerpack.view.welcome
Feature: What the Home search finds
  The "Search everywhere" box of the Home page searches as the user types, half a second after the
  typing stops, and lists what it finds under a header per kind: an entity by its name (a script
  under Scripts), the projects that carry a tag by "#tag", and nothing at all for a text nothing
  matches. The results are read once the search is over — the results host is aria-busy from the
  first keystroke until the last search path has answered — so an empty answer is claimed empty,
  not caught before it arrives. The script and the project are the feature's own, made through the
  API and deleted at its end, so no answer depends on what the stand holds: the old spec searched
  "tutorials" and "#demo", which a minimal stand has none of. Translated from playwright-public
  browse/search.test.ts (Browse-Search-01..03) and Browse-Nav-03 of browse/nav.test.ts: the Home icon
  of the Browse toolbar brings back the Home page with its search box.

  Background:
    Given user is logged in

  Scenario: The Home icon of the Browse toolbar brings back the Home page with its search box
    Given the browse panel is open
    And user clicks on "Open text" icon inside browse toolbar
    And the "Import text" view should be current
    And home search should be hidden
    When user clicks on "Home" icon inside browse toolbar
    Then the "Home" view should be current
    And home search should be visible
    And no errors should have been logged

  Scenario: A script is found by its name
    Given a script "BddPowerSearch{time}" is on the server:
      """
      //language: javascript
      //output: int x
      x = 1;
      """
    When user types "BddPowerSearch{time}" into home search
    Then the search should have finished
    And the "Scripts" search results should list "BddPowerSearch{time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A tag search finds the projects that carry the tag
    Given a project named "bdd-power-search-{time}" tagged "bddsearch{time}" is on the server
    When user types "#bddsearch{time}" into home search
    Then the search should have finished
    And the "Projects" search results should list "bdd-power-search-{time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A text nothing matches finds nothing
    When user types "zzqxnomatchqq" into home search
    Then the search should have finished
    And the search results should be empty
    And home search should have value "zzqxnomatchqq"
    And no errors should have been logged
    And no error or warning balloon should have been shown
