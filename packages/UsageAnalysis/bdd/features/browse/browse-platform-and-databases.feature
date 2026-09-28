@browse @realizes:views.browse
Feature: The Platform and Databases sections of the Browse tree
  The administrative section and the database browser: what each lists, and that a node opens the
  view it names. Translated from the manual cases Browse-Platform-01, -02 and Browse-DB-01, -02
  (playwright-public/browse/platform.test.ts, db.test.ts, browse_manual_tests2.md sections 9 and 10).

  The two "lists" scenarios are tagged @full-stand: they name what a full stand carries, and a
  minimal stack has fewer providers and fewer Platform sections. Everything else here holds on any
  stand — the old spec claimed only Postgres, and Plugins/Credentials/Functions/Users/Groups/Roles,
  for the same reason.

  Browse-Platform-03 (Platform is hidden from a user without the privilege) needs the second
  account. Browse-DB-03 (Schema Browser) is writable — this file already assumes the same CHEMBL
  connection — and is simply not written yet. Browse-DB-04 (a saved query) names a query only some
  stands carry. Browse-DB-05 asserted only that the pane held some text, which holds on a stale
  pane; Browse-DB-06 never clicked Browse > Summary at all, so it did not test GROK-16857.

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
