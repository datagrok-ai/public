@journey @realizes:chem.cp.mpo-profile-crud @realizes:chem.int.mpo-profile-sync
Feature: The MPO profile editor, from the Browse tree to a deleted profile
  Browse > Apps > Chem > MPO profiles opens the Manage Profiles list. A profile created there keeps
  the property name typed into it key by key (GROK-20918), is saved under the name typed on its tab,
  and appears in the list. It clones into a copy that can be saved (GROK-20832), renames to a free
  name without a "Replace profile" prompt (GROK-20822), and both are deleted from the list's own row
  menu, which is where the walk ends with the list holding neither.

  A table opened after the profile tab is offered in that tab's Dataset list (GROK-20821), and a
  data-driven profile is saved under the name and the description typed on the tab, with the pMPO
  model file beside it (GROK-20830). GROK-20830 was reported against a "Save model" dialog prefilled
  with the dataset name; since the profiles moved to the mpo domain table there is no such dialog,
  and what the defect lost — the typed name and description on the stored profile — is what is
  claimed here. The data-driven scenario needs the EDA package on the stand.

  The property name is typed once, at a person's pace, and read back: the library's own typing
  retypes until a field holds the text, which a field that drops keys on its first rename would
  survive. The profile title swallows the select-all key, so a title is retyped after its text is
  selected with the pointer.

  Every profile the user saves during the feature is deleted when it ends, whatever name it carries,
  together with the pMPO model file a data-driven save writes to System:AppData/EDA/pmpo — the
  scenarios delete their own profiles through the interface, and this is what is left if one fails.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And no MPO profile named "BDD MPO {run}, BDD MPO {run} (Copy), BDD MPO renamed {run}, BDD MPO data-driven {run}" is on the server
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Chem tree node inside browse tree is expanded

  Scenario: The app opens from the Browse tree
    When user clicks on Apps---Chem---MPO-profiles tree node inside browse tree
    Then the "MPO Profiles" view should be current
    And "Create profile" button should be visible
    And 0 MPO profiles named "BDD MPO {run}" should be on the server
    And no errors should have been logged

  Scenario: A property name typed key by key stays in its field
    When user clicks on "Create profile" button
    Then the "Untitled Profile" view should be current
    When user clicks on "+ Add Property" link
    Then first MPO property should have value "NewProperty1"
    When user types "HeavyAtomCount" key by key into first MPO property
    Then first MPO property should have value "HeavyAtomCount"
    And no errors should have been logged

  Scenario: Save stores the profile under the name typed on the tab
    When user selects the text of MPO profile title
    And user types "BDD MPO {run}" into MPO profile title
    And user clicks on "Save" button
    Then an info balloon containing "BDD MPO {run}" should have been shown
    And 1 MPO profile named "BDD MPO {run}" should be on the server
    And the MPO profile "BDD MPO {run}" should have the properties "HeavyAtomCount"
    And the "BDD MPO {run}" view should be current
    And no errors should have been logged

  Scenario: The list holds the saved profile
    When user switches to the "MPO Profiles" view
    Then "BDD MPO {run}" MPO profile should be visible
    And no errors should have been logged

  Scenario: A clone can be saved and then exists
    When user clicks on actions of "BDD MPO {run}" MPO profile
    And user picks "Clone" from the open menu
    Then the "BDD MPO {run} (Copy)" view should be current
    And 0 MPO profiles named "BDD MPO {run} (Copy)" should be on the server
    And "Save" button should be enabled
    When user clicks on "Save" button
    Then an info balloon containing "BDD MPO {run} (Copy)" should have been shown
    And 1 MPO profile named "BDD MPO {run} (Copy)" should be on the server
    And the MPO profile "BDD MPO {run} (Copy)" should have the properties "HeavyAtomCount"
    When user switches to the "MPO Profiles" view
    Then "BDD MPO {run} (Copy)" MPO profile should be visible
    And no errors should have been logged

  Scenario: A profile renamed to a free name is saved without a replace prompt
    When user clicks on actions of "BDD MPO {run}" MPO profile
    And user picks "Edit" from the open menu
    Then the "BDD MPO {run}" view should be current
    When user selects the text of MPO profile title
    And user types "BDD MPO renamed {run}" into MPO profile title
    And user clicks on "Save" button
    Then an info balloon containing "BDD MPO renamed {run}" should have been shown
    And 1 MPO profile named "BDD MPO renamed {run}" should be on the server
    And 0 MPO profiles named "BDD MPO {run}" should be on the server
    And the "BDD MPO renamed {run}" view should be current
    And no errors should have been logged

  Scenario: A table opened after the profile tab is offered in its Dataset list
    When user switches to the "MPO Profiles" view
    And user clicks on "Create profile" button
    Then Dataset input should have value ""
    Given user opens drugs-props-train dataset
    And user switches to the "Untitled Profile" view
    Then Dataset input should offer "drugs-props-train"
    And no errors should have been logged

  Scenario: A data-driven profile is stored under the name and description typed on the tab
    When user selects "Data-driven" in Method input
    And user selects "drugs-props-train" in Dataset input
    Then "Training data" heading should be visible
    When user selects the text of MPO profile title
    And user types "BDD MPO data-driven {run}" into MPO profile title
    And user types "made from the CNS column" into MPO profile description
    And user clicks on "Save" button
    Then an info balloon containing "BDD MPO data-driven {run}" should have been shown
    And 1 MPO profile named "BDD MPO data-driven {run}" should be on the server
    And the MPO profile "BDD MPO data-driven {run}" should have the description "made from the CNS column"
    And the pMPO model file of "BDD MPO data-driven {run}" should hold its name and the description "made from the CNS column"
    And no errors should have been logged

  Scenario: Delete in the row menu removes the profile from the list and the server
    When user switches to the "MPO Profiles" view
    And user clicks on actions of "BDD MPO renamed {run}" MPO profile
    And user picks "Delete" from the open menu
    Then "Delete profile" dialog should be visible
    When user clicks on OK button in "Delete profile" dialog
    Then 0 MPO profiles named "BDD MPO renamed {run}" should be on the server
    And "BDD MPO renamed {run}" MPO profile should be absent
    When user clicks on actions of "BDD MPO {run} (Copy)" MPO profile
    And user picks "Delete" from the open menu
    And user clicks on OK button in "Delete profile" dialog
    Then 0 MPO profiles named "BDD MPO {run} (Copy)" should be on the server
    And "BDD MPO {run} (Copy)" MPO profile should be absent
    And no errors should have been logged
