@browse @realizes:views.browse
Feature: The Apps and Dashboards sections of the Browse tree
  Opening an application from the tree, what the Dashboards node opens, and a dashboard that
  survives a save and a reopen. Translated from the manual cases Browse-Apps-01, -02, -03 and
  Browse-Dash-01, -02 (playwright-public/browse/apps.test.ts, dash.test.ts,
  browse_manual_tests2.md sections 6 and 8).

  Browse-Apps-03 (the tooltip and the details of an application, GROK-19638) is claimed on an app of
  the Chem package: the tooltip names the app, what it does and the package it comes from. The
  five-second half of Browse-Apps-01 (GROK-20032) is claimed on the list reopened from a collapsed
  Apps node.

  Browse-Dash-02 and -03 (a dashboard opened from the list, and the context panel following from one
  dashboard to the next, GROK-19934) are claimed on the dashboards the Chem package ships
  (chemical_space_demo, demo_activity_cliffs) — fixtures wherever Chem is, which the suite already
  needs for spgi — besides the dashboard this feature saves itself and removes again.

  The Model Hub (Browse-ModelHub-01..04, GROK-17896, GROK-19740, GROK-19965, GROK-19628) is claimed on
  a JavaScript model the feature saves and deletes again: the catalog lists it, Uncategorized opens
  to it, a click previews it, a hover explains it, a double click keeps its view through the next
  click, and its menu offers Run.

  The Misc group is what a full stand carries: a stand with no ungrouped application has none.

  The Tutorials app itself (its RUN and the tracks it docks) is the Tutorials package's own feature.

  Background:
    Given user is logged in
    And the browse panel is open

  @full-stand
  Scenario: The Apps section lists the installed applications
    Given Apps tree node inside browse tree is expanded
    Then Tutorials tree node inside browse tree should be visible
    And Misc tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An application opens from the tree and becomes the current object
    Given Apps tree node inside browse tree is expanded
    When user clicks on Tutorials tree node inside browse tree
    Then the "Tutorials" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # Browse-ModelHub-01: the Model Catalog opening from Apps is where GROK-17896 / GROK-17664 were.
  Scenario: The Model Hub opens from the Compute group
    Given the "Compute2" package is installed
    And Apps tree node inside browse tree is expanded
    And Apps---Compute tree node inside browse tree is expanded
    When user clicks on Apps---Compute---Model-Hub tree node inside browse tree
    Then the "Model Hub" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Dashboards node opens the list of projects
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    And the page address should contain "/projects"
    # the gallery's container is built with the view, so its presence says nothing: the claim is
    # on a card that has to have come from the server
    And card of gallery should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A dashboard saved from a table view opens again with its viewers
    Given user opens demog-1000 dataset
    And user adds a scatter plot viewer
    When user saves the current view as project "BDD-Browse-Dash-{run}"
    And user closes all views
    And user opens the "BDD-Browse-Dash-{run}" project
    Then the current view should hold at least 2 viewers
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The application list comes back within five seconds of opening Apps
    Given the "Compute2" package is installed
    And the "Chem" package is installed
    And user collapses Apps tree node inside browse tree
    And Apps---Compute tree node inside browse tree should be hidden
    When user expands Apps tree node inside browse tree
    Then Apps---Compute tree node inside browse tree should become visible within 5 seconds
    And Apps---Chem tree node inside browse tree should be visible
    And no errors should have been logged

  Scenario: An application's tooltip says what it does and which package it comes from
    Given the "Chem" package is installed
    And Apps tree node inside browse tree is expanded
    And Apps---Chem tree node inside browse tree is expanded
    And Apps---Chem---Reactions tree node inside browse tree is expanded
    When user hovers over Apps---Chem---Reactions---Reaction-Enumerator tree node inside browse tree
    Then tooltip should contain text "Reaction Enumerator"
    And tooltip should contain text "Forward-reaction library enumeration"
    And tooltip should contain text "Package"
    And tooltip should contain text "Chem"
    And no errors should have been logged

  Scenario: A dashboard the Chem package ships opens from the Dashboards gallery with its viewers
    Given the "Chem" package is installed
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    When user types "chemical_space_demo" into gallery search
    And user double-clicks on "ChemicalSpaceDemo" gallery card
    Then the current view should hold at least 2 viewers
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # each dashboard is searched by its own name: a shorter search brings up every dashboard of the stand that matches
  Scenario: The context panel follows from one dashboard to the next
    Given the "Chem" package is installed
    And the context panel is open
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    When user types "chemical_space_demo" into gallery search
    And user clicks on "ChemicalSpaceDemo" gallery card
    Then the context panel should show "chemical_space_demo"
    When user types "demo_activity_cliffs" into gallery search
    And user clicks on "DemoActivityCliffs" gallery card
    Then the context panel should show "demo_activity_cliffs"
    And context panel should not contain text "chemical_space_demo"
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Rule: A model the feature saves is in the Model Hub
    Background:
      Given the "Compute2" package is installed
      And a script "BddBrowseModel" is on the server:
        """
        //language: javascript
        //meta.role: model
        //description: A model a BDD feature saved
        //input: int x = 1
        //output: int result
        result = x + 1;
        """
      And Apps tree node inside browse tree is expanded
      And Apps---Compute tree node inside browse tree is expanded

    Scenario: The Model Hub catalog lists the model
      When user clicks on Apps---Compute---Model-Hub tree node inside browse tree
      Then the "Model Hub" view should be current
      And "BddBrowseModel" link should be visible
      And no errors should have been logged

    Scenario: Uncategorized opens to the model, a hover explains it and a click previews it
      Given Apps---Compute---Model-Hub tree node inside browse tree is expanded
      When user expands Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree
      Then Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree should be visible
      When user hovers over Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree
      Then tooltip should contain text "A model a BDD feature saved"
      When user clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree
      Then the "BddBrowseModel preview" view should be current
      And no errors should have been logged
      And no error or warning balloon should have been shown

    Scenario: A double click keeps the model's view and its menu offers Run
      Given Apps---Compute---Model-Hub tree node inside browse tree is expanded
      And Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree is expanded
      When user double-clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree
      Then the "BddBrowseModel preview" view should be current
      When user clicks on Dashboards tree node inside browse tree
      Then Projects view should be visible
      And "BddBrowseModel preview" view should be present
      # the click on Dashboards moved the tree: the model's group is opened again before its menu
      Given Apps---Compute---Model-Hub---Uncategorized tree node inside browse tree is expanded
      And Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree should be visible
      When user right-clicks on Apps---Compute---Model-Hub---Uncategorized---BddBrowseModel tree node inside browse tree
      Then the open menu should list "Run..."
      When user closes the context menu
      Then no errors should have been logged
      And no error or warning balloon should have been shown
