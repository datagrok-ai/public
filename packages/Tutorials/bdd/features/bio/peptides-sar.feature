@tutorials @serial @realizes:tutorials.peptides-sar
Feature: The Peptides SAR tutorial
  Walks Bioinformatics > Peptides SAR from its card to the end: the SAR analysis launched from the Bio
  menu with a similarity threshold of 90, the monomer grid and its WebLogo header, the Sequence
  Variability Map in both modes, the Most Potent Residues and Logo Summary Table settings, a pie chart
  aggregation added to the clusters, and the analysis re-run with a dendrogram. Each step is claimed as
  ticked and as done — the analysis ready (peptides-sar-ready), the tooltip, the selection and its
  clearing, the viewers made current, the mutation cliff pairs, the map's mode, the aggregated column.
  Translated from playwright-tests/e2e/tutorials/peptides-sar.test.ts, which swept the Sequence
  Variability Map's canvas for a non-empty cell; the map and the grid's WebLogo report their cells as
  areas. Needs Bio, Peptides and EDA (MCL); the clustering runs in the browser.

  Fixed in the tutorial for this translation: the counter was 18 for 24 steps (one more where WebGPU is
  available, which headless browsers have not); seven steps completed on their own after 30 minutes
  with nothing done.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Peptides SAR" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Peptides SAR tutorial
    When user starts the "Peptides SAR" tutorial
    Then the tutorial progress should be 1 of 24
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then the tutorial step "On the Top Menu, click Bio > Analyze > SAR..." should be done
    And "Analyze Peptides" dialog should be visible
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    Then the tutorial step "Click the Gear Icon (⚙)" should be done
    When user enters "90" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Then the tutorial step "Set Similarity Threshold to 90" should be done
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the tutorial step "Click OK to start analysis" should be done
    And the "peptides-sar-ready" custom event should have fired
    And the tutorial step "Wait for analysis to complete" should be done
    When user clicks on NEXT button in hint popup
    Then the tutorial step "Click NEXT to proceed" should be done

    When user hovers over the "cell 1 of 2" area of grid
    Then the tutorial step "Hover over a monomer cell in the main table grid" should be done
    And tooltip should be visible
    When user clicks on the "Aca at 3" area of grid
    Then the tutorial step "Click any WebLogo letter" should be done
    And some rows should be selected
    When user clicks on NEXT button in hint popup
    Then the tutorial step "Click NEXT to proceed" should be done 2 times
    When user presses Escape
    Then the tutorial step "Press Esc to clear selection" should be done
    And no rows should be selected

    When user hovers over Sequence Variability Map viewer
    And user clicks on settings icon of Sequence Variability Map viewer
    Then the tutorial step "Open SVM settings (gear)" should be done
    And the context panel should show "Sequence Variability Map"
    When user clicks on NEXT button in hint popup
    Then the tutorial step "Click NEXT to proceed" should be done 3 times

    When user clicks on the "cell 1Nal at 1" area of Sequence Variability Map viewer
    Then the tutorial step "Click a Mutation Cliffs cell" should be done
    And "Mutation Cliffs pairs" pane in context panel should be visible
    When user clicks on NEXT button in hint popup
    Then the tutorial step "Click NEXT to proceed" should be done 4 times
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    Then the tutorial step "Switch SVM mode to Invariant Map" should be done
    And the "mode" reading of Sequence Variability Map viewer should be "Invariant Map"
    When user clicks on NEXT button in hint popup
    Then the tutorial step "Click NEXT to proceed" should be done 5 times

    When user hovers over Most Potent Residues viewer
    And user clicks on settings icon of Most Potent Residues viewer
    Then the tutorial step "Open Most Potent Residues settings (gear)" should be done
    And the context panel should show "Most Potent Residues"
    When user clicks on NEXT button in hint popup
    Then the tutorial step "Click NEXT to proceed" should be done 6 times
    When user hovers over Logo Summary Table viewer
    And user clicks on settings icon of Logo Summary Table viewer
    Then the tutorial step "Open Logo Summary Table settings (gear)" should be done
    And the context panel should show "Logo Summary Table"
    When user clicks on "..." button in "Columns" property
    Then "Select columns..." dialog should be visible
    When user types "14" into "Search" input in "Select columns..." dialog
    And user toggles the "14" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the tutorial step "Add pie chart aggregation for position 14" should be done
    And "columns" property of Logo Summary Table viewer should be "14"
    When user clicks on NEXT button in hint popup
    Then the tutorial step "Scroll horizontally in Logo Summary Table. Click NEXT to proceed to next step" should be done

    When user clicks on "Peptides analysis settings" icon
    Then the tutorial step "Click the Wrench to update analysis configuration" should be done
    And "Peptides settings" dialog should be visible
    When user checks "Dendrogram" input in "Peptides settings" dialog
    Then the tutorial step "Check Dendrogram" should be done
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Peptides settings" dialog
    Then the tutorial step "Click OK to re-run analysis" should be done
    And the "peptides-sar-ready" custom event should have fired
    When user clicks on OK button in hint popup

    And the "Peptides SAR" tutorial should be completed
    And the tutorial should have listed 24 steps
    And the tutorial progress should be 24 of 24
    And no hint should be shown
    And no errors should have been logged
