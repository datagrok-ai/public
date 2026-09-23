@realizes:powerpack.search.power-pack @realizes:powerpack.view.welcome
Feature: Enter in the Home search box
  Pressing Enter in the "Search everywhere" box of the Home page (GROK-18656: Enter on "QA" threw).
  The search starts on the text itself, half a second after the typing stops; Enter acts on the
  suggestion menu under the box. The case GROK-18656 was about is Enter with that menu shown and no
  suggestion highlighted, so every query that brings up suggestions waits for them, claims that none
  is highlighted, and only then presses Enter; a query with no suggestions presses Enter once its
  search is over and no menu is shown. Every query then finishes its search with nothing logged and
  no error balloon, and shows a known answer. The answer is read once the search is over — the
  results say so, `aria-busy` cleared when the last search path has answered — so a query that finds
  nothing is claimed to find nothing, not caught before its answer arrives. Translated from TestTrack
  PowerPack/power-search-enter.md.

  What each query shows on a stand is what the stand holds: the categories claimed are the ones every
  stand has (functions and help pages). "Project[0-9]+" is not read as a pattern, and "1+1" is not
  evaluated (a formula is evaluated only with parentheses): both finish with nothing shown at all. Of
  the suggestions, "PDB ID, e.g. 4AKZ" puts 4AKZ into the box, and the search then shows the PDB entry
  of 4AKZ (a link to it, found through data.rcsb.org, so the stand needs to reach it); the New Users
  suggestions carry no value of their own and leave the box as it was, so the suggestion walked
  through is the PDB one. The md's "dem" brings up no suggestion on a stand; "DG" brings up two, which
  is what the arrow keys need.

  Not translated: "the search input retains focus or a sensible follow-up view is loaded" — the
  box's focus is not claimed; the view that stays current is the Home page with its results.

  Background:
    Given user is logged in

  Scenario Outline: Enter on "<query>" with its suggestions shown and none highlighted finds functions and help pages
    When user types "<query>" into home search
    Then the search suggestions should be "<suggestions>"
    And the highlighted search suggestion should be "none"
    When user presses Enter in home search
    Then the search should have finished
    And home search should have value "<query>"
    And home widgets panel should be hidden
    And the page address should contain "search?q="
    And the search results should list the categories "Functions, Help"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | query | suggestions                                                                                                                      |
      | QA    | PDB ID, e.g. 4AKZ                                                                                                                |
      | new   | New Users Today \| New users This Month \| New users This Year \| New users last 3 months \| New user last 7 days \| New users yesterday |
      | user  | DGUSER-{User Login}                                                                                                              |

  Scenario: Enter on "a", with no suggestion shown, finds functions and help pages
    When user types "a" into home search
    Then the search should have finished
    And no search suggestion should be shown
    When user presses Enter in home search
    Then home search should have value "a"
    And the search results should list the categories "Functions, Help"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Enter on "1+1" with its suggestion shown and none highlighted finds nothing
    When user types "1+1" into home search
    Then the search suggestions should be "PDB ID, e.g. 4AKZ"
    And the highlighted search suggestion should be "none"
    When user presses Enter in home search
    Then the search should have finished
    And home widgets panel should be hidden
    And the search results should show nothing
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Enter on "Project[0-9]+", with no suggestion shown, finds nothing
    When user types "Project[0-9]+" into home search
    Then the search should have finished
    And no search suggestion should be shown
    When user presses Enter in home search
    Then home widgets panel should be hidden
    And the search results should show nothing
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: "new" lists the Add New Column function
    When user types "new" into home search
    Then the search should have finished
    And the "Functions" category of the search results should list "Add New Column"
    And no errors should have been logged

  Scenario: The arrow keys walk the suggestions, and Enter takes the highlighted one
    When user types "DG" into home search
    Then the search suggestions should be "DGUSER-{User Login} | PDB ID, e.g. 4AKZ"
    And the highlighted search suggestion should be "none"
    When user presses ArrowDown in home search
    Then the highlighted search suggestion should be "DGUSER-{User Login}"
    When user presses ArrowDown in home search
    Then the highlighted search suggestion should be "PDB ID, e.g. 4AKZ"
    When user presses ArrowUp in home search
    Then the highlighted search suggestion should be "DGUSER-{User Login}"
    When user presses ArrowDown in home search
    And user presses Enter in home search
    Then home search should have value "4AKZ"
    And the page address should contain "search?q=4AKZ"
    And the search should have finished
    And "PDB: 4AKZ" link in home search results should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
