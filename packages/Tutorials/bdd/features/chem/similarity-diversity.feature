@tutorials @serial @realizes:tutorials.similarity-diversity-search
Feature: The Similarity and Diversity Search tutorial
  Walks Cheminformatics > Similarity and Diversity Search from its card to the end: both search
  viewers from the Chem menu, a card of each made current, the similarity viewer's reference locked,
  a reference molecule pasted into the sketcher, the reference explored, a property column added to
  the diversity cards and colour-coded. Each step is claimed as ticked and as done — the viewers,
  the row a clicked card showed made current, Follow Current Row off, the search run again around
  the pasted molecule, the column on the cards, the colour coding.
  Translated from playwright-tests/e2e/tutorials/similarity-diversity.test.ts, which clicked the cards
  at guessed offsets and the column picker's canvas by pixels. RDKit runs in the browser.

  Fixed in the tutorial for this translation: the gear step completed on a click on any gear of the page; the Follow Current Row
  and Molecule Properties steps waited for texts in the context panel ("1 / 31") and now read the
  viewers' own properties; the Edit hint looked for a class the viewer does not carry
  (`.similarity-search-edit`, the icon is `chem-similarity-search-edit`), so it showed nothing; the
  menu and icon hints were captured when their step began.
  The context panel drops a new current object within 2 s of a property edit and within 1 s of a
  scripted `grok.shell.o =` (by design, GROK-21024), which a feature is faster than: the cards are
  clicked beside the drawing (a click on the drawing makes the molecule current that way), a grid
  cell click releases the property-edit guard before Explore, and the 1 s after Explore — which
  nothing releases or reports — is waited out before the diversity link.
  Not claimed: the panes of the explored molecule's context panel — several are built by server-side
  scripts or outside lookups, which a feature does not run.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Similarity and Diversity Search" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Similarity and Diversity Search tutorial
    When user starts the "Similarity and Diversity Search" tutorial
    Then the tutorial progress should be 1 of 12
    When user picks "Chem > Search > Similarity Search..." from the top menu
    Then the tutorial step "On the Top Menu, click Chem > Search > Similarity Search..." should be done
    And Chem Similarity Search viewer should be visible
    When user picks "Chem > Search > Diversity Search..." from the top menu
    Then the tutorial step "Next, click Chem > Search > Diversity Search..." should be done
    And Chem Diversity Search viewer should be visible

    When user clicks on card 2 of Chem Similarity Search viewer beside its drawing
    Then the tutorial step "On the Most similar structures viewer, click the molecule next to the reference molecule" should be done
    And the current row should be the row of the clicked card
    When user clicks on card 3 of Chem Diversity Search viewer beside its drawing
    Then the tutorial step "Now, click any molecule in the diversity viewer" should be done
    And the current row should be the row of the clicked card
    # the table view puts the moved current cell into the panel 750 ms later, over a viewer made current meanwhile;
    # after another molecule nothing changes that a claim could read (the shell keeps an object of the same type)
    When user waits 1 second

    When user hovers over Chem Similarity Search viewer
    And user clicks on settings icon of Chem Similarity Search viewer
    Then the tutorial step "Hover over similarity viewer and click gear icon in the right top corner of the viewer to open settings" should be done
    And the context panel should show "Chem Similarity Search"
    When user expands Misc category
    And user unchecks "Follow Current Row" property
    Then the tutorial step "Under Misc, clear the Follow Current Row checkbox" should be done
    And "followCurrentRow" property of Chem Similarity Search viewer should be "false"

    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user clicks on "Edit" icon in Chem Similarity Search viewer
    Then the tutorial step "On the reference molecule, click the Edit icon" should be done
    And sketcher dialog should be visible
    When user types "CNc1nc(Nc2ccc(Br)cc2)nc(N)c1[N+](=O)[O-]" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then the tutorial step "Set new reference molecule" should be done
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered

    When user clicks on a "smiles" cell of grid other than the current one
    Then the context panel should show the current cell
    When user hovers over card 1 of Chem Similarity Search viewer
    And user clicks on "More" icon in Chem Similarity Search viewer
    And user picks "Explore" from the open menu
    Then the tutorial step "Hover over the reference molecule, click the More icon, and then Explore" should be done

    # Explore made the molecule current through grok.shell.o, which drops the next change within a second
    # (freezeCurrentObjectUntil) and announces nothing when the second is over
    When user waits 1 second
    And user clicks on "Tanimoto, Morgan" link in Chem Diversity Search viewer
    Then the tutorial step "In the top right corner of the diversity viewer, click Tanimoto, Morgan" should be done
    And the context panel should show "Chem Diversity Search"
    When user clicks on "..." button in "Molecule Properties" property
    Then "Select columns..." dialog should be visible
    When user types "NumValenceElectrons" into "Search" input in "Select columns..." dialog
    And user toggles the "NumValenceElectrons" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the tutorial step "Select NumValenceElectrons column" should be done
    And the "card properties" reading of Chem Diversity Search viewer should be "NumValenceElectrons"

    When user clicks on the "header smiles" area of grid
    And user moves the current cell of grid to the "NumValenceElectrons" column
    And user picks "Color Coding > Linear" from the context menu of the "header NumValenceElectrons" area of grid
    Then the tutorial step "Add Color Coding for NumValenceElectrons column" should be done
    And "NumValenceElectrons" column should be color-coded linearly

    And the "Similarity and Diversity Search" tutorial should be completed
    And the tutorial should have listed 12 steps
    And the tutorial progress should be 12 of 12
    And no hint should be shown
    And no errors should have been logged
