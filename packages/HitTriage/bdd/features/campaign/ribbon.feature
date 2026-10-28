Feature: The Hit Design campaign view's ribbon and crumbs
  A campaign opens as a platform table view whose ribbon ends with the package's own group: Refresh, Add new row,
  Merge, Calculate additional properties, Edit campaign permissions (with edit permission), Apply template layout
  (when the template has one), then the Progress tracker (when it has stages) and Submit... buttons. The platform's
  Save is taken out. The view's name slot shows the crumbs "Hit Design › <campaign>"; going back by the first one is
  crumbs.feature.

  The group is pinned to the ribbon's end, so at 1366 px with the toolbox and the context panel open the primary
  buttons stay in view; below a width the small icons fold into a "More actions" menu whose items are the icons'
  tooltips ("..." after the ones that open a dialog).

  The campaign is made for the feature through the JS API and opened by its link (bindings/campaign.ts); its template
  carries the package's shipped layout file, so the Apply template layout icon shows.

  Background:
    Given user is logged in
    And the Hit Design campaign "BDD Ribbon" with the stages "Idea, Design, Make" and a table of:
      | Molecule | Stage  | V-iD | MW    |
      | CCO      | Idea   |      | 46.07 |
      | c1ccccc1 | Design |      | 78.11 |
    When user opens the campaign "BDD Ribbon" by its link

  Scenario: The ribbon ends with the package's actions, and the crumbs name the app and the campaign
    Then the last group of the ribbon should list "Refresh, Add new row, Merge another table into the campaign, Calculate additional properties, Edit campaign permissions, Apply template layout (Loaded from mounted file storage), Progress tracker, Submit..."
    And "Save" button in toolbar should be hidden
    And "Hit Design" breadcrumb should be visible
    And "BDD Ribbon" breadcrumb should be visible
    And the page address should contain "campaignId=BDDBDDRIBBON-1"
    When user hovers over "Refresh" icon in toolbar
    Then tooltip should have text "Refresh Campaign"
    And no error or warning balloon should have been shown

  Scenario: Apply template layout loads the template's layout file into the view
    Then the open tableview should have 0 Bar chart viewers
    When user clicks on "Apply template layout (Loaded from mounted file storage)" icon in toolbar
    Then the open tableview should have 4 Bar chart viewers
    And the table should have 2 rows

  Scenario: At 1366 px the primary buttons stay in view and the package's icons are one menu away
    Given the toolbox pane is shown
    And the context panel is open
    When the window is 1366 by 768 px
    Then "Submit..." button in toolbar should not be clipped
    And "Progress tracker" button in toolbar should not be clipped
    And "Merge another table into the campaign" icon in toolbar should be hidden
    When user clicks on "More actions" icon in toolbar
    Then the open menu should list "Refresh Campaign"
    And the open menu should list "Merge another table into the campaign..."
    And the open menu should list "Edit campaign permissions..."
