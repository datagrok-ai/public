@journey
Feature: PeptiHit landing and new campaign: the campaigns table, a campaign opened by its code, and START with Helm
  PeptiHit is Hit Design for peptides: the same campaigns table under its own hero, a campaign opened with its code
  in the address and the Helm column first (the name of its view is not claimed), and START making a campaign with
  Helm, Molecule, Stage, V-iD and the HELM editor open on its first row. All but the first two start on a fresh app.

  Made for the feature through the JS API (bindings/landing.ts): the campaigns "BDD PH Alpha" (status "BDD
  Synthesis"), "BDD PH Beta" and "BDD PH Locked" (the second account's, which the running account may only view), and
  one template per START, each with an optional field of its own ("Start note" …) so the form shows which template it
  was built for before START is pressed (the form is rebuilt after the template file is read), none with compute,
  "BDD PH Layout" with a scatter plot of MolWt in its layout. The campaigns the templates' keys start go with them at
  feature end, with their rows in the hitdesign database. Needs Bio and Helm on the stand.

  The Campaign Name names the new campaign, which takes nothing over from the campaign opened before it, not even from
  one the user may only view. Opening a campaign from the landing reached without the crumb closes the previous
  campaign's view. A new campaign from a template whose layout has a MolWt scatter plot starts that plot empty until
  the first calculation adds the column, with no 'x: Column "MolWt" does not exist.' in it; the scatter plot does not
  report that message as its error (getWidgetStatus().error is null), so the claim reads the plot's text.

  Background:
    Given user is logged in
    And the "Bio" package is installed
    And the "Helm" package is installed
    And the campaigns tables start from their defaults, put back at feature end
    And the "PeptiHit" templates are on the server:
      | name                 | key          | fields              | stages       | layout                |
      | BDD PH Start {time}  | BDDPHS{time} | Start note: String  | Design, Test |                       |
      | BDD PH Named {time}  | BDDPHN{time} | Named note: String  | Design       |                       |
      | BDD PH Fresh {time}  | BDDPHF{time} | Fresh note: String  | Design       |                       |
      | BDD PH After {time}  | BDDPHA{time} | After note: String  | Design       |                       |
      | BDD PH Layout {time} | BDDPHY{time} | Layout note: String | Design       | scatter plot of MolWt |
    And the "PeptiHit" campaigns are on the server:
      | code           | name                 | status               | created   | molecules | author       | can edit |
      | BDDPHL{time}-1 | BDD PH Alpha {time}  | BDD Synthesis {time} | 2099/1/14 | 3         | me           | everyone |
      | BDDPHL{time}-2 | BDD PH Beta {time}   | In Progress          | 2099/1/13 | 2         | me           | everyone |
      | BDDPHL{time}-3 | BDD PH Locked {time} | In Progress          | 2099/1/12 | 2         | sharing user | author   |
    And the campaigns table of "PeptiHit" is sorted by "Created" descending
    And user opens the "PeptiHit" app

  Scenario: The PeptiHit landing lists the campaigns as Hit Design does, under its own hero
    Then the "PeptiHit" view should be current
    And the current view should show the text "Sketch Helm sequences in the spreadsheet"
    And the campaigns table should have the columns "Name, Code, Created, Author, Last Modified by, Molecules, Status"
    And the "Status" column of the "BDD PH Alpha {time}" campaign row should be "BDD Synthesis {time}"
    And permissions icon of "BDD PH Alpha {time}" campaign row should be visible
    And permissions icon of "BDD PH Locked {time}" campaign row should be hidden

  Scenario: A campaign's name opens it with the Helm column first, and the crumb comes back
    When user clicks on "BDD PH Alpha {time}" link
    Then the page address should contain "campaignId=BDDPHL{time}-1"
    And 1 campaign view should be open
    And the table should have the columns "Helm, Molecule, Stage, V-iD"
    And "Helm" column should have semantic type "Macromolecule"
    When user clicks on "PeptiHit" link
    Then the "PeptiHit" view should be current

  Scenario: START makes a campaign with Helm, Molecule, Stage and V-iD and opens the HELM editor on it
    When user closes all views
    And user opens the "PeptiHit" app
    And user selects "BDD PH Start {time}" in Template input
    Then "Start note" input should be visible
    When user clicks on Start button
    Then HELM editor dialog should be visible
    And the "BDDPHS{time}-1" view should be current
    And the page address should contain "campaignId=BDDPHS{time}-1"
    And the table should have the columns "Helm, Molecule, Stage, V-iD"
    And the value of "Stage" column in row 1 should be "Design"
    And "Helm" column should have semantic type "Macromolecule"
    And the "PeptiHit" campaign "BDDPHS{time}-1" should be on the server
    When user clicks on CANCEL button in HELM editor dialog
    Then HELM editor dialog should be absent

  Scenario: The Campaign Name typed before START names the new campaign
    When user closes all views
    And user opens the "PeptiHit" app
    And user selects "BDD PH Named {time}" in Template input
    Then "Named note" input should be visible
    When user types "BDD PH Named Run {time}" into "Campaign Name" input
    And user clicks on Start button
    Then HELM editor dialog should be visible
    And the "PeptiHit" campaign "BDDPHN{time}-1" should be on the server
    And the "PeptiHit" campaign "BDDPHN{time}-1" on the server should have the name "BDD PH Named Run {time}"
    When user clicks on CANCEL button in HELM editor dialog

  Scenario: A new campaign does not take over the campaign opened before it
    When user closes all views
    And user opens the "PeptiHit" app
    And user clicks on "BDD PH Alpha {time}" link
    Then the page address should contain "campaignId=BDDPHL{time}-1"
    When user clicks on "PeptiHit" link
    And user selects "BDD PH Fresh {time}" in Template input
    Then "Fresh note" input should be visible
    When user clicks on Start button
    Then HELM editor dialog should be visible
    And the "PeptiHit" campaign "BDDPHF{time}-1" should be on the server
    And the "PeptiHit" campaign "BDDPHF{time}-1" on the server should have the status "In Progress"
    When user clicks on CANCEL button in HELM editor dialog

  Scenario: After a campaign the user may only view, START still makes and opens the new campaign
    When user closes all views
    And user opens the "PeptiHit" app
    And user clicks on "BDD PH Locked {time}" link
    Then the page address should contain "campaignId=BDDPHL{time}-3"
    When user clicks on "PeptiHit" link
    And user selects "BDD PH After {time}" in Template input
    Then "After note" input should be visible
    When user clicks on Start button
    Then the "PeptiHit" campaign "BDDPHA{time}-1" should be on the server
    And the "BDDPHA{time}-1" view should be current
    And no error or warning balloon should have been shown

  Scenario: Opening a campaign from the landing closes the campaign opened before it
    When user closes all views
    And user opens the "PeptiHit" app
    And user clicks on "BDD PH Beta {time}" link
    Then the page address should contain "campaignId=BDDPHL{time}-2"
    When user switches to the "PeptiHit" view
    And user clicks on "BDD PH Locked {time}" link
    Then the page address should contain "campaignId=BDDPHL{time}-3"
    And 1 campaign view should be open

  Scenario: A new campaign from a template with a layout opens with viewers that work
    When user closes all views
    And user opens the "PeptiHit" app
    And user selects "BDD PH Layout {time}" in Template input
    Then "Layout note" input should be visible
    When user clicks on Start button
    Then HELM editor dialog should be visible
    When user clicks on CANCEL button in HELM editor dialog
    Then the "BDDPHY{time}-1" view should be current
    And scatter plot viewer should be visible
    And scatter plot viewer should not contain text "does not exist"
