@browse @realizes:views.browse
Feature: The Apps and Dashboards sections of the Browse tree
  Opening an application from the tree, what the Dashboards node opens, and a dashboard that
  survives a save and a reopen. Translated from the manual cases Browse-Apps-01, -02, -03 and
  Browse-Dash-01, -02 (playwright-public/browse/apps.test.ts, dash.test.ts,
  browse_manual_tests2.md sections 6 and 8).

  Browse-Apps-01 also asks that the list appear within five seconds (GROK-20032). That half is not
  translated: a wall-clock threshold on a shared stand reports the stand's load, not the product's,
  and the suite has no vocabulary for it. What is claimed is the part that would actually have
  caught the defect — the list arrives and holds the applications it should.

  Browse-Dash-02 asks for a dashboard opened from the list. The stand's own dashboards are not a
  fixture — dev carries over two hundred of them, created by whoever — so the dashboard this
  feature opens is one it saves itself and removes again at the end. It is opened by name rather
  than from the gallery, so the list-to-view step of the manual case is not claimed.
  Browse-Dash-03 (the context panel must not repeat the previous dashboard's content, GROK-19934)
  needs two fixtures with distinguishable panes and is left for the round that can tell them apart.

  The Model Hub scenario claims that the view opens, and no more: its catalog is the Compute
  package's own gallery, not the platform card gallery the shared vocabulary names, so a catalog
  that opens empty (the GROK-17896 shape) is not caught here.

  Browse-ModelHub-02, -03 and -04 are NOT translated, and that is a real gap: they carried
  GROK-19740, GROK-19965 and GROK-19628, the reproductions playwright-public/browse/KNOWN_BUGS.md
  is built around. Clicking and double-clicking a model in the tree is claimed by nothing here.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: The Apps section lists the installed applications
    Given Apps tree node inside browse tree is expanded
    Then Tutorials tree node inside browse tree should be visible
    And Chem tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An application opens from the tree and becomes the current object
    Given Apps tree node inside browse tree is expanded
    When user clicks on Tutorials tree node inside browse tree
    Then the "Tutorials" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Hovering an application does not break its details
    Given the context panel is open
    And Apps tree node inside browse tree is expanded
    When user hovers over Tutorials tree node inside browse tree
    And user clicks on Tutorials tree node inside browse tree
    Then the context panel should show "Tutorials"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # Browse-ModelHub-01: the Model Catalog opening from Apps is where GROK-17896 / GROK-17664 were.
  # Tagged so a stand without the Compute package can exclude it rather than go red.
  @compute
  Scenario: The Model Hub opens from the Compute group
    Given Apps tree node inside browse tree is expanded
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
