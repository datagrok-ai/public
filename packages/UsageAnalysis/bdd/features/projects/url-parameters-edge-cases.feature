@realizes:views.projects
Feature: URL parameters edge cases: a value of the wrong type
  Guards for the URL Parameters of a dashboard saved from a query with two parameters, an integer
  "lim" and a string "typeName" (made through the API on the Datagrok database: writing a query is
  the database features' subject). The Save dialog lists both, each with its alias.

  A value of the wrong type in the project's address (?lim=abc) is refused with a balloon that
  names the parameter.

  The project's address is /p/<namespace>.<name>, the namespace being the running account's login
  (`projects-url-parameters`). Names are letters and digits only and carry the run's time; the query
  and the project are removed when the feature starts and ends.

  Background:
    Given user is logged in
    And the context panel is open
    And simple mode is off
    And the browse panel is open
    And no project named "BDDUrlEdgeP{time}" is on the server
    And a query "BDDUrlEdgeQ{time}" on "System:Datagrok" is on the server:
      """
      --input: int lim = 2
      --input: string typeName = "Project"
      select name from entity_types where name <> @typeName order by name limit @lim
      """
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDUrlEdgeQ{time} tree node inside browse tree and the query it runs completes
    Then the "BDDUrlEdgeQ{time}" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDUrlEdgeP{time}" into Name text input in "Save project" dialog
    And user clicks on "URL Parameters" button in "BDDUrlEdgeQ{time}" project table in "Save project" dialog
    Then first "URL alias" text input in "Save project" dialog should have value "lim"

  @edge-cases
  Scenario: A value of the wrong type in the address is refused with a balloon naming the parameter
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDUrlEdgeP{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDUrlEdgeP{time}" dialog
    Then the "Share BDDUrlEdgeP{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user opens the address "/p/{login}.BDDUrlEdgeP{time}?lim=abc"
    Then an error or warning balloon matching "lim" should have been shown
