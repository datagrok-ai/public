@journey @serial @realizes:views.projects @realizes:views.functions
Feature: A dashboard on a query with a parameter: Toolbox > Source, URL parameters and saved values
  A query with a string parameter (typeName = "Project") is made on the Datagrok database's
  entity_types with New SQL Query... and run from the Browse tree. Before the save, the copy icon of
  Toolbox > Source gives a /func/ link that runs the query; in the Save dialog the table's URL
  Parameters expose typeName, its alias is changed to "type", and after the save the same icon gives
  the dashboard's /p/ link. Reopened, the dashboard takes a new value typed into Source (the link
  follows it before REFRESH; the copy icon puts it on the clipboard and turns into a check), saves
  it into the project and a third value into a copy independently, and unticking the parameter under
  the sliders icon of Source asks for a save ("Save dashboard to apply changes") that clears the hint. Translated
  from the TestTrack case Projects/project-url-parameters.

  Parked (see the request document): every step that opens an address — the copied links in a new
  tab, the project's URL from Links... with ?type=, ?typeName= and no parameter, and after the
  parameter is taken out (steps 5, 7-10, 12 and 16 of the md in part); the typeName check box of the
  Save dialog's URL Parameters (a bare check box no element reaches) and the check mark of typeName
  in the sliders icon's menu (`checked` does not read a Dart menu item's aria-checked).

  Known failures: no ticket (the md's note asks to record it) — right after the first save the Source
  link carries the parameter's own name (?typeName=Project), not the alias "type" set in the Save
  dialog; after a reopen it carries the alias. No ticket either: once typeName is unticked under
  the sliders icon (the "Save dashboard to apply changes" hint shows), the Source link still ends
  with ?type=Script. GROK-20930 — the Save dialog's "Share link:" line does not follow the name typed
  into the dialog; GROK-20929 — the sliders icon of Source appears only after the project is
  reopened, not right after the first save. Each is a scenario of its own claiming what the md
  expects.

  entity_types has no "Package" row on localhost (the md allows another value of its name column
  then), so "Script" stands for the md's "Package"; "Project" and "Script" are core entity types.
  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time. The query and the two projects are removed when the feature starts and ends. It is serial:
  the Dashboards search is shared with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And no project named "BDDUrlParamProj{time}" is on the server
    And no project named "BDDUrlParamCopy{time}" is on the server
    And no query named "BDDUrlParamQuery{time}" is on the server

  Scenario: A query with a typeName parameter is made on entity_types
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "New SQL Query..." from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the current view should be a DataQueryView view
    When user enters "BDDUrlParamQuery{time}" into Name input
    And user replaces the code of code editor with '--input: string typeName = "Project"'
    And user appends "select name from entity_types where name = @typeName" to code editor
    And user clicks on Save button
    Then 1 query named "BDDUrlParamQuery{time}" should be on the server
    When user closes the current view
    Then no errors should have been logged

  Scenario: The query runs from the tree with its default value
    Given the browse panel is open
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user double-clicks on Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree
    Then the "BDDUrlParamQuery{time}" view should be current
    And the "rows" reading of grid should be 1
    And the "text of cell 1 of name" reading of grid should be "Project"

  Scenario: Before the save, the Source link runs the query
    Given the toolbox pane is shown
    Then Source pane in toolbox should be visible
    When user hovers over copy icon in Source pane in toolbox
    Then tooltip should contain text "/func/"
    And tooltip should contain text "BDDUrlParamQuery{time}"
    And tooltip should contain text "typeName=%22Project%22"
    And tooltip should contain text "run=true"

  Scenario: The Save dialog exposes typeName under the table's URL Parameters
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "BDDUrlParamQuery{time}" project table in "Save project" dialog should be visible
    And Data sync switch in "BDDUrlParamQuery{time}" project table in "Save project" dialog should be checked
    When user enters "BDDUrlParamProj{time}" into Name text input in "Save project" dialog
    And user clicks on "URL Parameters" button in "BDDUrlParamQuery{time}" project table in "Save project" dialog
    Then "URL alias" text input in "Save project" dialog should have value "typeName"

  @known-failure @realizes:GROK-20930
  Scenario: The Share link line follows the name typed into the dialog
    Then "Save project" dialog should contain text ".BDDUrlParamProj{time}?typeName="

  Scenario: The alias is changed to "type" and the dashboard is saved
    When user enters "type" into "URL alias" text input in "Save project" dialog
    Then "Save project" dialog should contain text "?type=Project"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDUrlParamProj{time}" uploaded' should have been shown
    And "Share BDDUrlParamProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDUrlParamProj{time}" dialog
    Then the "Share BDDUrlParamProj{time}" dialog should close

  Scenario: After the save, the Source link opens the dashboard without run=true
    When user moves the pointer away from Source pane in toolbox
    And user hovers over copy icon in Source pane in toolbox
    Then tooltip should contain text "/p/"
    And tooltip should contain text ".BDDUrlParamProj{time}?"
    And tooltip should not contain text "run=true"

  @known-failure
  Scenario: Right after the save, the Source link already carries the alias
    Then tooltip should contain text ".BDDUrlParamProj{time}?type=Project"

  @known-failure @realizes:GROK-20929
  Scenario: Right after the save, Source offers the sliders icon
    Then sliders-h icon in Source pane in toolbox should be visible
    When user hovers over sliders-h icon in Source pane in toolbox
    Then tooltip should contain text "Choose which parameters the dashboard link carries"

  Scenario: Reopened, Source offers the sliders icon, and the link follows a new value before REFRESH
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDUrlParamProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDUrlParamProj{time} gallery card
    Then the "BDDUrlParamQuery{time}" view should be current
    And the "rows" reading of grid should be 1
    Given the toolbox pane is shown
    Then sliders-h icon in Source pane in toolbox should be visible
    When user enters "Script" into "Type Name" input in Source pane in toolbox
    And user hovers over copy icon in Source pane in toolbox
    Then tooltip should contain text "Copy link:"
    And tooltip should contain text ".BDDUrlParamProj{time}?type=Script"
    When user clicks on REFRESH button in Source pane in toolbox
    Then the "text of cell 1 of name" reading of grid should be "Script"
    And the "rows" reading of grid should be 1
    When user clicks on copy icon in Source pane in toolbox
    Then check icon in Source pane in toolbox should be visible
    And the clipboard should contain text ".BDDUrlParamProj{time}?type=Script"

  Scenario: The new value is saved into the project (GROK-20006)
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on "Creation script" button in "BDDUrlParamQuery{time}" project table in "Save project" dialog
    Then "BDDUrlParamQuery{time}" project table in "Save project" dialog should contain text '"Script"'
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDUrlParamProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDUrlParamProj{time} gallery card
    Then the "BDDUrlParamQuery{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Script"
    And the "rows" reading of grid should be 1

  Scenario: A copy saved with a third value leaves the original's value alone
    Given the toolbox pane is shown
    When user enters "Project" into "Type Name" input in Source pane in toolbox
    And user clicks on REFRESH button in Source pane in toolbox
    Then the "text of cell 1 of name" reading of grid should be "Project"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    And user enters "BDDUrlParamCopy{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDUrlParamCopy{time}" uploaded' should have been shown
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDUrlParamCopy{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDUrlParamCopy{time} gallery card
    Then the "BDDUrlParamQuery{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Project"
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDUrlParamProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDUrlParamProj{time} gallery card
    Then the "BDDUrlParamQuery{time}" view should be current
    And the "text of cell 1 of name" reading of grid should be "Script"

  Scenario: The sliders icon takes the parameter out of the link, and the save applies it
    Given the toolbox pane is shown
    When user clicks on sliders-h icon in Source pane in toolbox
    Then the open menu should list "typeName"
    When user picks "typeName" from the open menu
    Then "Save dashboard to apply changes" text in toolbox should be visible

  @known-failure
  Scenario: With the parameter taken out, the Source link no longer carries it
    When user moves the pointer away from Source pane in toolbox
    And user hovers over copy icon in Source pane in toolbox
    Then tooltip should contain text "/p/"
    And tooltip should not contain text "?type="

  Scenario: The save applies the change and the hint goes
    When user moves the pointer away from Source pane in toolbox
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And sliders-h icon in Source pane in toolbox should be visible
    And "Save dashboard to apply changes" text in toolbox should be hidden

  Scenario: The query and the projects are deleted
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDUrlParam" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Delete Project" from the context menu of BDDUrlParamProj{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user picks "Delete Project" from the context menu of BDDUrlParamCopy{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDUrlParamProj{time}" should be on the server
    And 0 projects named "BDDUrlParamCopy{time}" should be on the server
    When user clicks on "Refresh" icon inside browse toolbar
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user collapses Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree
    And user picks "Delete" from the context menu of Databases---Postgres---Datagrok---BDDUrlParamQuery{time} tree node inside browse tree
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 queries named "BDDUrlParamQuery{time}" should be on the server
