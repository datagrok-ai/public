@skip @journey @serial @realizes:powerpack.db-explorer
Feature: Identifiers of a database connection
  TODO: development only. The feature is set aside and does not run in CI: it saves and removes
  identifiers configurations in the one file every configuration of the stand is kept in, and after a
  save the page it runs on keeps the handlers it registered. Take @skip off while working on Configure
  Identifiers or the db-explorer library, run it on a local stand, and put the tag back before
  committing.

  Configure Identifiers... in the context menu of a database connection asks for the schema the
  identifiers live in and opens the editor on it. An identifier is a semantic type bound to a column
  of a table; the editor previews a row of the table, and Save keeps the configuration on the server,
  where it shows under the connection as "Identifiers - <schema>", opens again from there, and is
  removed from its context menu.

  The connection is the platform's own Postgres database (System:Datagrok), which every stand has.
  Postgres has no catalogs to choose from, so the dialog asks for the schema alone; a database with
  catalogs is identifiers-catalogs.feature. What an identifier does once the page has loaded it - the
  detection of the semantic type, the card and the context panel - is claimed there, on its own
  database; here the claims are the dialog, the editor and the saved configuration. The configuration
  the feature saves is removed when it ends.

  Background:
    Given user is logged in
    And no identifiers configuration of the "System:Datagrok" connection is on the server
    And the browse panel is open

  Scenario: The dialog asks for the schema alone and opens on public
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Configure Identifiers..." from the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then "Select primary schema for Identifiers Configuration" dialog should be visible
    And Schema input in "Select primary schema for Identifiers Configuration" dialog should have value "public"
    And Catalog input in "Select primary schema for Identifiers Configuration" dialog should be absent
    When user clicks on OK button in "Select primary schema for Identifiers Configuration" dialog
    Then the "Datagrok.public Identifiers" view should be current
    And identifiers editor should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An identifier is previewed from its table and saved
    When user clicks on add identifier icon
    Then "Add Identifier" dialog should be visible
    When user enters "BDD_EVENT_TYPE" into Semantic Type input in "Add Identifier" dialog
    And user selects "event_types" in Table input in "Add Identifier" dialog
    And user selects "friendly_name" in Column input in "Add Identifier" dialog
    Then identifier preview should contain text "friendly_name"
    When user clicks on ADD button in "Add Identifier" dialog
    Then identifiers editor should contain text "BDD_EVENT_TYPE"
    When user clicks on Save button in identifiers editor
    Then an info balloon containing "Identifier configuration saved" should have been shown
    And the identifiers configurations of the "System:Datagrok" connection on the server should cover "public"
    And the saved identifiers of "public" in the "System:Datagrok" connection should be "BDD_EVENT_TYPE"
    And no errors should have been logged

  Scenario: The saved configuration is listed under the connection and opens from there
    Given the browse panel is open
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree
    Then identifiers editor should contain text "BDD_EVENT_TYPE"
    And no errors should have been logged

  Scenario: Configure Identifiers on the same schema opens the saved configuration
    When user picks "Configure Identifiers..." from the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then "Select primary schema for Identifiers Configuration" dialog should contain text "public already has identifiers. OK opens them for editing."
    When user clicks on OK button in "Select primary schema for Identifiers Configuration" dialog
    Then an info balloon containing "already exists" should have been shown
    And the "Datagrok.public Identifiers" view should be current
    And identifiers editor should contain text "BDD_EVENT_TYPE"

  Scenario: Remove takes the configuration off the connection and the server
    Given the browse panel is open
    When user picks "Remove..." from the context menu of Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree
    Then "Remove Identifier Configuration" dialog should be visible
    When user clicks on OK button in "Remove Identifier Configuration" dialog
    Then an info balloon containing "Identifier configuration removed" should have been shown
    And Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree should be absent
    And the identifiers configurations of the "System:Datagrok" connection on the server should cover ""
    And no errors should have been logged
