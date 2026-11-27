@guide @tutorial
Feature: How-to · Molecules
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Activity cliffs: Chem > Analyze > Activity Cliffs... on a
  table with an Activity column plots the molecules by similarity (UMAP, in the browser) with a line between each
  pair of similar molecules whose activity differs a lot; Show only cliffs keeps the pairs, hovering a line shows
  both molecules, and the cliffs button lists every pair in a table. activity_cliffs (Chem's demo of 269 kinase
  inhibitors) is opened off camera with the grid showing chembl_tid, smiles and Activity; nothing is made on the
  server.

  Scenario: Find activity cliffs
    Find pairs of similar molecules with very different activity.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens activity_cliffs dataset
    And the grid shows only the columns "chembl_tid, smiles, Activity"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Chem > Analyze > Activity Cliffs
    When user picks "Chem > Analyze > Activity Cliffs..." from the top menu
    Then "Activity Cliffs" dialog should be visible
    # caption: Check the molecule and activity columns, then OK
    When user clicks on OK button in "Activity Cliffs" dialog
    Then scatter plot viewer should be visible
    And the "cliffs" reading of scatter plot viewer should be 20
    When the video says "Lines join similar molecules whose activity differs a lot"
    # caption: Show only cliffs hides the rest
    And user toggles "Show only cliffs" input in scatter plot viewer
    Then the "only cliffs" reading of scatter plot viewer should be "true"
    # caption: Point at a line to see the pair
    When user hovers over the first free "line" area of scatter plot viewer
    Then tooltip should be visible
    When the video says "The tooltip shows both molecules and their activities" at the top
    Given the tooltip is hidden
    # caption: The cliffs button lists every pair
    When user clicks on "20 cliffs" button in scatter plot viewer
    Then the open tableview should have 2 grid viewers
    # caption: Drag its top border up to see more pairs
    When user drags the top border of second grid viewer by 250 pixels up
    And the video says "Each row is one pair: both molecules and how far their activities differ"
    Then the table should have 270 rows
    And no errors should have been logged
