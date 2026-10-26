@browse @edge-cases @realizes:views.browse
Feature: Gallery search edge cases: characters of the search syntax typed as text
  Guards for the search box of a Browse gallery when what is typed is part of the search syntax
  rather than a name: an opening parenthesis, a double quote, and a condition with no value
  ("name = "). The search matches names (help/datagrok/concepts/project/space.md, "Searching
  spaces"), and no space of this feature has such a character in its name, so each search empties
  the list and logs no error, and clearing the box brings the space back. spaces-search claims the
  search by a name and by a name nothing carries.

  The space is this feature's own, made through the API, named with the run's time, and deleted when
  the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And a space named "BDDEdgeSearch{time}" is on the server
    And Spaces tree node inside browse tree is expanded
    When user clicks on Spaces tree node inside browse tree
    And user enters "BDDEdgeSearch{time}" into gallery search
    Then BDDEdgeSearch{time} link in gallery should be visible

  Scenario Outline: A search of <what> empties the list with no error, and clearing it brings the space back
    When user enters "<text>" into gallery search
    Then BDDEdgeSearch{time} link in gallery should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clears gallery search
    Then BDDEdgeSearch{time} link in gallery should be visible

    Examples:
      | what                      | text                       |
      | an opening parenthesis    | (                          |
      | a double quote            | \"                         |
      | a condition with no value | name =                     |
