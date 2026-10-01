@browse @realizes:views.browse
Feature: The Platform and Databases sections of the Browse tree
  The administrative section and the database browser: what each lists, and that a node opens the
  view it names. Translated from the manual cases Browse-Platform-01, -02 and Browse-DB-01, -02
  (playwright-public/browse/platform.test.ts, db.test.ts, browse_manual_tests2.md sections 9 and 10).

  The two "lists" scenarios are tagged @full-stand: they name what a full stand carries, and a
  minimal stack has fewer providers and fewer Platform sections. Everything else here holds on any
  stand — the old spec claimed only Postgres, and Plugins/Credentials/Functions/Users/Groups/Roles,
  for the same reason — or is gated on what it needs: a package, or the database of a connection
  outside the stand.

  Browse-DB-04 (a saved query opens from the tree) is claimed on the Orders query of the Dbtests
  package's NorthwindTest connection, which the queries features run on too; Browse-DB-05 on a
  table of its public schema. Browse-DB-03 (the Schema Browser) is claimed in
  connections-schema.feature. Browse-DB-06 (CHEMBL > Browse > Summary, GROK-16857) opens the
  CHEMBL package's own browse views, which need the CHEMBL database — a stand without it would
  report the database, not the tree — so it is not written. The Platform matrix (Node-Platform-01)
  opens every administrative node and claims the view it names, not only the absence of errors.

  The package manager (Browse/package-manager.md 1 and 3, GROK-16545) and local-deploy-ui.md 4 are
  claimed read-only: the Plugins view lists the packages the stand installed, and a package's menu
  offers no Delete. The version switch and the install and uninstall cases (2, 4-8) change a package
  every other feature uses and are not written; local-deploy-ui.md 1-3, 5 and 6 are what a fresh
  deployment does once, not something a shared stand can repeat.

  Browse-Tree-07 (a node the user may not read) is claimed through the second account: the first
  account's own connection is not in that user's tree at all, beside the shared Datagrok one.
  Browse-Platform-03 (Platform is hidden from a user without the privilege) is not: the second
  account holds whatever the stand gives users — the local stand gives them Platform — so the
  claim would describe the stand's defaults, not the rule.

  A node below the top level is named by its full tree path ("Files---Demo"), which is what the
  platform writes into its own `name` attribute. Several sections carry a node called Demo, Files
  or App Data, and a bare name matches whichever of them another feature happened to leave
  open: the tree remembers its expanded set per user, across features and across runs.

  Background:
    Given user is logged in
    And the browse panel is open

  @full-stand
  Scenario: The Platform section lists what an administrator manages
    Given Platform tree node inside browse tree is expanded
    Then the following elements should be visible:
      | Platform---Plugins tree node inside browse tree           |
      | Platform---Credentials tree node inside browse tree       |
      | Platform---Functions tree node inside browse tree         |
      | Platform---Users tree node inside browse tree             |
      | Platform---Groups tree node inside browse tree            |
      | Platform---Roles tree node inside browse tree             |
      | Platform---Predictive-models tree node inside browse tree |
      | Platform---Dockers tree node inside browse tree           |
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A Platform node opens the view it is named after
    Given Platform tree node inside browse tree is expanded
    When user clicks on Platform---Users tree node inside browse tree
    Then the "Users" view should be current
    When user clicks on Platform---Groups tree node inside browse tree
    Then the "Groups" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @full-stand
  Scenario: The Databases section lists the connected providers
    Given Databases tree node inside browse tree is expanded
    Then the following elements should be visible:
      | Databases---Postgres tree node inside browse tree |
      | Databases---MySQL tree node inside browse tree    |
      | Databases---Oracle tree node inside browse tree   |
      | Databases---MariaDB tree node inside browse tree  |
      | Databases---MS-SQL tree node inside browse tree   |
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A provider opens down to its connections
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    Then Databases---Postgres---Datagrok tree node inside browse tree should be visible
    And Databases---Postgres---CHEMBL tree node inside browse tree should be visible
    When user collapses Databases---Postgres tree node inside browse tree
    Then Databases---Postgres---Datagrok tree node inside browse tree should be hidden
    And Databases---Postgres tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # expanding the connection reads its schema from db.datagrok.ai, outside the stand (PostgresTest is its grok name:
  # two more connections are called NorthwindTest)
  Scenario: A saved query opens from the tree with its details
    Given the context panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And the stand has a reachable "PostgresTest" connection
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Orders tree node inside browse tree
    Then the context panel should show "Orders"
    And the page address should contain "/func/Dbtests.PostgresOrders"
    And "Script" accordion header in context panel should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A table of a schema shows its details
    Given the context panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And the stand has a reachable "PostgresTest" connection
    And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---orders tree node inside browse tree
    Then the context panel should show "orders"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: Platform > <node> opens the <view> view
    Given the context panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on Platform---<node> tree node inside browse tree
    Then the "<view>" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | node              | view         |
      | Plugins           | Plugins      |
      | Credentials       | Credentials  |
      | Functions         | Functions    |
      | Roles             | Roles        |
      | MCP-Servers       | MCP Servers  |
      | Predictive-models | Models       |
      | Sync              | Sync         |
      | Layouts           | View layouts |
      | URL-Aliases       | URL Aliases  |
      | Settings          | Settings     |
      | Sticky-Meta       | Schemas      |

  # the tree lists Notebooks only where a package registers the notebooks functions
  Scenario: Platform > Notebooks opens the Notebooks view
    Given the "Notebooks" package is installed
    And the context panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on Platform---Notebooks tree node inside browse tree
    Then the "Notebooks" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # Browse-Tree-07, through the second account: a plain user who was given nothing of what the
  # first account owns does not get the first account's connection in the tree at all
  Scenario: Another user does not see a connection nobody shared with them
    Given a "Postgres" connection named "BDD-Browse-Private-{run}" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    Then Databases---Postgres---BDD-Browse-Private-{run} tree node inside browse tree should be visible
    When user signs in as the sharing user
    And the browse panel is open
    Then the sharing user should be signed in
    When Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    Then Databases---Postgres---Datagrok tree node inside browse tree should be visible
    And Databases---Postgres---BDD-Browse-Private-{run} tree node inside browse tree should be absent
    And no errors should have been logged
    When user signs in as themselves again
    Then the running account should be signed in

  Scenario: The Plugins view lists the installed packages and offers no Delete for one
    Given the "PowerPack" package is installed
    And the "Tutorials" package is installed
    And Platform tree node inside browse tree is expanded
    When user clicks on Platform---Plugins tree node inside browse tree
    Then the "Plugins" view should be current
    And the following elements should be visible:
      | "PowerPack" gallery card |
      | "Tutorials" gallery card |
    When user opens the context menu of "PowerPack" gallery card
    Then the open menu should list "Uninstall"
    And the open menu should not list "Delete"
    And the open menu should not list "Delete..."
    When user closes the context menu
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  # the Dockers view warns in a balloon on a stand whose Grok Spawner is off
  Scenario: Platform > Dockers opens the Dockers view
    Given the stand runs the "Grok Spawner" service
    And the context panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on Platform---Dockers tree node inside browse tree
    Then the "Dockers" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown
