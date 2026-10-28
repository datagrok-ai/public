@journey
Feature: Hit Design template editor: name and key checks, fields, stages, submit function, layout, CREATE and CANCEL
  The "+" icon of the Template picker opens the New Template view, its crumb leading back to Hit Design: Details
  (Name, Key, Layout), Additional fields (rows of Name, Type, Required, added with "Add field"), Stages (rows added
  with "Add stage"), Compute, Submit, and CANCEL and CREATE in a footer that stays in view. A name or a key that is
  empty or taken shows its error under the input, and CREATE is disabled while it does; CREATE writes <name>.json
  with the fields and stages, a row still being typed in included, the submit function and the layout, and comes
  back to the landing with the new template selected; CANCEL writes nothing, and the crumb back closes the editor,
  asking first once something was typed. Every scenario starts on a freshly opened app.

  The descriptors cannot be listed in this feature (some stands list them through a service of Chem's): the Compute
  section says why, and the rest of the editor works without them. Picking descriptors is a package test.

  A template "BDD Editor Existing" is made for the feature through the JS API (bindings/landing.ts) to collide with;
  the templates the scenarios name go at feature end. The layout file is the bdd project's fixtures/bar-chart.layout
  (a grid and a bar chart of Stage), opened from disk.

  Background:
    Given user is logged in
    And the descriptors cannot be listed, put back at feature end
    And the "Hit Design" templates are on the server:
      | name                       | key          |
      | BDD Editor Existing {time} | BDDEDX{time} |
    And the "Hit Design" templates named "BDD Editor New {time}, BDD Editor Cancelled {time}" are removed now and at feature end
    And user opens the "Hit Design" app

  Scenario: The editor has Details, Additional fields, Stages, Compute and Submit, with CANCEL and CREATE
    When user clicks on "New Template" icon
    Then the "New Template" view should be current
    And "Hit Design" link should be visible
    And template name input should be visible
    And template key input should be visible
    And Layout input should be visible
    And "Submit function" input should be visible
    And the current view should show the text "Additional fields"
    And the current view should show the text "Stages"
    And the current view should show the text "Compute"
    And the current view should show the text "Could not list the descriptors: the descriptors are not available"
    And Cancel button should be visible
    And Create button should be visible
    When user clicks on Cancel button
    Then the "Hit Design" view should be current

  Scenario: A name or a key that is empty or taken shows its error, a free one clears it
    When user closes all views
    And user opens the "Hit Design" app
    And user clicks on "New Template" icon
    And user types "BDD Editor Existing {time}" into template name input
    Then the template editor should show the error "A template with this name already exists"
    When user types "BDD Editor Free {time}" into template name input
    Then the template editor should not show the error "A template with this name already exists"
    When user types "BDDEDX{time}" into template key input
    Then the template editor should show the error "BDDEDX{time} is used by BDD Editor Existing {time}"
    When user types "BDDEDF{time}" into template key input
    Then the template editor should not show the error "BDDEDX{time} is used by BDD Editor Existing {time}"
    When user clears template name input
    Then the template editor should show the error "Value can't be empty"
    When user clicks on Cancel button
    Then the "Hit Design" view should be current

  Scenario: While the name is taken, CREATE is disabled
    When user closes all views
    And user opens the "Hit Design" app
    And user clicks on "New Template" icon
    Then Create button should be disabled
    When user types "BDD Editor Existing {time}" into template name input
    And user types "BDDEDF{time}" into template key input
    Then the template editor should show the error "A template with this name already exists"
    And Create button should be disabled
    When user types "BDD Editor Free {time}" into template name input
    Then Create button should be enabled
    When user clicks on Cancel button
    Then the "Hit Design" view should be current

  Scenario: CREATE writes the fields, stages, submit function and layout, the rows being typed included, and comes back with the template selected
    When user closes all views
    And user opens the "Hit Design" app
    And user clicks on "New Template" icon
    And user types "BDD Editor New {time}" into template name input
    And user types "BDDEDN{time}" into template key input
    And user uploads "fixtures/bar-chart.layout" through "Open file" icon in Layout input
    And user adds the required field "Owner" of type "String" in the template editor
    Then the template editor should list the fields "Owner (String, required)"
    When user types the field "Due" of type "Date" in the template editor without adding it
    And user adds the stage "Design" in the template editor
    Then the template editor should list the stages "Design"
    When user types the stage "Make" in the template editor without adding it
    And user selects "Demo File Submit" in "Submit function" input
    And user clicks on Create button
    Then the "Hit Design" view should be current
    And Template input should have the value "BDD Editor New {time}"
    And "Owner" input should be visible
    And the "Hit Design" template "BDD Editor New {time}" on the server should have the key "BDDEDN{time}"
    And the "Hit Design" template "BDD Editor New {time}" on the server should have the fields "Owner (String, required), Due (Date)"
    And the "Hit Design" template "BDD Editor New {time}" on the server should have the stages "Design, Make"
    And the "Hit Design" template "BDD Editor New {time}" on the server should have the submit "HitTriage:demoFileSubmit"
    And the "Hit Design" template "BDD Editor New {time}" on the server should have the layout "a layout"
    And no error or warning balloon should have been shown

  Scenario: CANCEL closes the editor and writes no template
    When user closes all views
    And user opens the "Hit Design" app
    And user clicks on "New Template" icon
    And user types "BDD Editor Cancelled {time}" into template name input
    And user types "BDDEDC{time}" into template key input
    And user clicks on Cancel button
    Then the "Hit Design" view should be current
    And Template input should offer the choice "BDD Editor Existing {time}"
    And Template input should not offer the choice "BDD Editor Cancelled {time}"
    And the "Hit Design" template "BDD Editor Cancelled {time}" should not be on the server

  Scenario: The crumb back to Hit Design closes the editor, and asks first once something was typed
    When user closes all views
    And user opens the "Hit Design" app
    And user clicks on "New Template" icon
    And user clicks on "Hit Design" link
    Then the "Hit Design" view should be current
    And the "New Template" view should not be open
    When user clicks on "New Template" icon
    And user types "BDD Editor Cancelled {time}" into template name input
    And user clicks on "Hit Design" link
    Then "Discard the new template?" dialog should be visible
    When user clicks on Cancel button in "Discard the new template?" dialog
    Then the "New Template" view should be current
    And template name input should have the value "BDD Editor Cancelled {time}"
    When user clicks on "Hit Design" link
    And user clicks on Discard button in "Discard the new template?" dialog
    Then the "Hit Design" view should be current
    And the "New Template" view should not be open
    And the "Hit Design" template "BDD Editor Cancelled {time}" should not be on the server
