@browse @realizes:views.browse
Feature: The My stuff section of the Browse tree
  The section that gathers what belongs to the current user. Translated from the manual case
  Browse-MyStuff-01 (playwright-public/browse/mystuff.test.ts, browse_manual_tests2.md section 4).

  Recent and Favorites reload when they are opened again: the group goes loading → loaded, which the
  expand step waits for, so Browse-MyStuff-02 (a project just opened is in Recent), -03 and Fav-01
  (an entity added from its context menu is in Favorites, and out of it when the same item is
  picked again — "Add to favorites > Only for me" is a check since GROK-21108) are claimed on
  a reopened group. The favorite is a fixture connection the feature makes and deletes, taken out of
  the account's favorites before and after. -05 (Add to favorites from My Files, GROK-19848) needs
  the user's home share, which a stand names itself ("My files" on dev) and a stand without home
  storage does not have, so it is @full-stand. -06 (a new script lands in My stuff) is claimed in
  browse-navigation.feature, where Refresh brings it in. -04 (Shared with me grouped by who shared)
  needs the second account to own and share something first, which no step can make it do.

  Recent, Favorites and Shared with me are the section's own nodes. The rest are the buckets of
  the user's personal project (project_meta.dart bucketOf: Connections, Files, Dashboards,
  Scripts, Tables, Flows, Spaces), and a bucket is there only while the account owns something
  of that kind, so only the three fixed nodes are claimed. A bucket carries the same name as the
  section elsewhere in the tree; the full path tells them apart.

  A node below the top level is named by its full tree path ("Files---Demo"), which is what the
  platform writes into its own `name` attribute. Several sections carry a node called Demo, Files
  or App Data, and a bare name matches whichever of them another feature happened to leave
  open: the tree remembers its expanded set per user, across features and across runs.

  Background:
    Given user is logged in
    And the browse panel is open
    And My stuff tree node inside browse tree is expanded

  Scenario: My stuff gathers the user's own things
    Then the following elements should be visible:
      | My-stuff---Recent tree node inside browse tree          |
      | My-stuff---Favorites tree node inside browse tree       |
      | My-stuff---Shared-with-me tree node inside browse tree  |
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The section closes on its own twistie
    Given My-stuff---Recent tree node inside browse tree should be visible
    When user collapses My stuff tree node inside browse tree
    Then My-stuff---Recent tree node inside browse tree should be hidden
    And My-stuff---Favorites tree node inside browse tree should be hidden
    And My stuff tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # lowercase: a project saved through the JS API is listed by camelCaseToWords of its name
  # ("BDD-Browse-Recent" shows as "BD D- Browse- Recent")
  Scenario: A project just opened is in Recent
    Given no project named "bdd-browse-recent-{run}" is on the server
    And user opens demog-1000 dataset
    When user saves the current view as project "bdd-browse-recent-{run}"
    And user closes all views
    And user opens the "bdd-browse-recent-{run}" project
    Then "bdd-browse-recent-{run}" should be among the recently used entities on the server
    # the project brings its Toolbox, docked as a tab over Browse
    And user closes all views
    And the toolbox pane is hidden
    And the browse panel is open
    And My stuff tree node inside browse tree is expanded
    When user collapses My-stuff---Recent tree node inside browse tree
    And user expands My-stuff---Recent tree node inside browse tree
    Then My-stuff---Recent---bdd-browse-recent-{run} tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An entity added to favorites from its menu is in Favorites, and out when picked again
    Given a "Postgres" connection named "BDD-Browse-Fav-{run}" is on the server
    And "BDD-Browse-Fav-{run}" is not in favorites
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user picks "Add to favorites > Only for me" from the context menu of Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree
    Then "BDD-Browse-Fav-{run}" should be in favorites on the server
    When user collapses My-stuff---Favorites tree node inside browse tree
    And user expands My-stuff---Favorites tree node inside browse tree
    Then My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree should be visible
    When user picks "Add to favorites > Only for me" from the context menu of Databases---Postgres---BDD-Browse-Fav-{run} tree node inside browse tree
    Then "BDD-Browse-Fav-{run}" should not be in favorites on the server
    When user collapses My-stuff---Favorites tree node inside browse tree
    And user expands My-stuff---Favorites tree node inside browse tree
    Then My-stuff---Favorites---BDD-Browse-Fav-{run} tree node inside browse tree should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @full-stand
  Scenario: The home share is added to favorites from its menu
    Given "My files" is not in favorites
    And Files tree node inside browse tree is expanded
    When user picks "Add to favorites > Only for me" from the context menu of Files---My-files tree node inside browse tree
    Then "My files" should be in favorites on the server
    When user collapses My-stuff---Favorites tree node inside browse tree
    And user expands My-stuff---Favorites tree node inside browse tree
    Then My-stuff---Favorites---My-files tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
