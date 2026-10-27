@realizes:views.projects
Feature: URL parameters edge cases: a value of the wrong type, an emptied alias, one alias for two parameters
  Guards for the URL Parameters of a dashboard saved from a query with two parameters, an integer
  "lim" and a string "typeName" (made through the API on the Datagrok database: writing a query is
  the database features' subject). The Save dialog lists both, each with its alias.

  - A value of the wrong type in the project's address (?lim=abc) is refused with a balloon that
    names the parameter.
  - An emptied alias falls back to the parameter's own name (`syncMapping` in
    shell_project_url_params.dart), so the Share link keeps ?lim=. Today clearing it logs
    "NullError: method not found" and the Share link line goes: a known failure, no ticket yet.
  - Two parameters given one alias are expected to be refused when the dashboard is saved, with a
    message. Today the mapping is keyed by the alias and the second parameter silently replaces the
    first (`syncMapping`); that scenario is a known failure with no ticket yet (`@ticket-pending`).

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
    When user opens the address "/p/Admin.BDDUrlEdgeP{time}?lim=abc"
    Then an error or warning balloon matching "lim" should have been shown

  @edge-cases @known-failure:59 @ticket-pending
  Scenario: An emptied alias falls back to the parameter's own name
    When user clears first "URL alias" text input in "Save project" dialog
    Then no errors should have been logged
    And "Save project" dialog should contain text "?lim=2"
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And 0 projects named "BDDUrlEdgeP{time}" should be on the server

  @edge-cases @known-failure:71 @ticket-pending
  Scenario: Two parameters given one alias are refused when the dashboard is saved
    Today the second parameter silently replaces the first under the shared alias and the dashboard
    is saved. No ticket yet. Remove the tag when the save is refused with a message.
    When user enters "typeName" into first "URL alias" text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then an error or warning balloon matching "alias" should have been shown
    And "Save project" dialog should be visible
    And 0 projects named "BDDUrlEdgeP{time}" should be on the server
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
