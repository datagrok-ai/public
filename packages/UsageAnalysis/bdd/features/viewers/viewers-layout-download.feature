@viewers @realizes:viewers.chrome
Feature: The layout of a docked view downloads as a .layout file
  `viewers-docking-ui.md` ends by downloading the docked arrangement: View > Layout > Download saves
  the current view's layout as a `.layout` file, the one the md then drops back onto the view.
  `viewers-docking.feature` docks and restores through a saved layout but never downloads one. The
  file is claimed once, by its name and by holding the viewers the view had docked; what else a
  layout stores is the layout API's test. Dropping the file back from the desktop is not
  translated: a browser page cannot be handed a file dragged in from outside it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer
    And user adds a bar chart viewer
    Then the open tableview should have 1 histogram viewer
    And the open tableview should have 1 bar chart viewer

  Scenario: View > Layout > Download saves the view's layout with its viewers
    When user docks histogram viewer to the right edge of the view
    Then histogram viewer should be docked along the right edge of the view
    Given user watches downloads
    When user picks "View > Layout > Download" from the top menu
    Then a file matching "\.layout$" should have been downloaded
    And the downloaded file should contain "Histogram"
    And the downloaded file should contain "Bar chart"
    And no errors should have been logged
