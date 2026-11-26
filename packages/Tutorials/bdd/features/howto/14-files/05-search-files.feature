@guide @tutorial
Feature: How-to · Files
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Files found by name in a folder: Browse > Files > Demo clicked
  (its folder view), "aspirin" typed into the folder's search box finds the three aspirin files of its subfolders
  (chem/mol/aspirin.mol, chem/sdf/aspirin-2d.sdf, chem/sdf/aspirin-3d.sdf), the Search in subfolders icon switched
  off finds none of them (Demo itself holds no aspirin file) and on again, and a double-click on a result opens it
  (aspirin-2d.sdf, one molecule). Searching from Home is 01-getting-started/08's. Nothing is made on the server;
  the subfolder switch is the view's own, and the search box is cleared by the view closing.

  Scenario: Search files
    Find a file by its name in a folder and all the folders inside it.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When user clicks on browse tab
    And user expands Files tree node inside browse tree
    And user clicks on Files---Demo tree node inside browse tree
    Then the "Demo" view should be current
    When the video says "A folder has a search box: type part of a name"
    And user types "aspirin" into gallery search
    Then aspirin-2d.sdf link in gallery should be visible
    And aspirin-3d.sdf link in gallery should be visible
    And aspirin.mol link in gallery should be visible
    When the video says "It finds files in the folders inside it too"
    And user hovers over "Search in subfolders" icon inside gallery toolbar
    And user clicks on "Search in subfolders" icon inside gallery toolbar
    Then aspirin-2d.sdf link in gallery should be absent
    When the video says "Switched off, it looks in this folder only"
    And user clicks on "Search in subfolders" icon inside gallery toolbar
    Then aspirin-2d.sdf link in gallery should be visible
    When the video says "Double-click a file to open it"
    And user double-clicks on aspirin-2d.sdf link in gallery
    Then the table should have 1 row
    And no errors should have been logged
