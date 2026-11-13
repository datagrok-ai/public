@journey @realizes:chem.cp.sar-matrix
Feature: SAR Matrix Summary on component columns
  The SAR Matrix PROTAC demo project: 2,792 degraders whose warhead, linker and E3 ligand are columns of
  their own, with the linker as the core and the warhead across the matrix columns. Every scenario starts
  by putting its own tab on screen, and puts back what it changed. A component's best value opens its
  ranking, effects against the average compound: a row's compound count selects exactly the compounds it
  counts, and its grid icon opens that series on the SAR Matrix tab. A component's best swap opens the
  measured swaps, each with its weakest pair, whose grid icon selects the compounds carrying the new
  value. Lower is better re-reads every ranking and the key series without a rebuild. The dialog leaves
  Scaling at none, unavailable, for the negative solubility and offers -lg again for a positive activity.
  Another component across the columns rebuilds the matrices; a scaling the activity cannot take then leaves
  no matrix and every tab says why; Scaling none and the warhead back across the columns, set together,
  build the matrices again.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens the "SarMatrixProtacDemo" project
    And SAR Matrix Viewer viewer should have built its matrices

  Scenario: A component's ranking selects and opens what each row counts
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    And user clicks on value of Warhead component row
    Then Summary evidence should contain text "Warhead — effect vs the average compound"
    And value of first evidence row should have text "+1.65"
    And count badge of first evidence row should have text "8 cpd"
    When user clicks on count badge of first evidence row
    Then 8 rows should be selected
    When user clicks on open icon of first evidence row
    Then the "tab" reading of SAR Matrix Viewer viewer should be "SAR Matrix"
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    Then the "series" reading of SAR Matrix Viewer viewer should be the text of name of first evidence row
    And no errors should have been logged

  Scenario: A component's best swap selects the compounds with the new value
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    And user clicks on swap of Warhead component row
    Then Summary evidence should contain text "Swapping Warhead — what it bought"
    And first evidence row should contain text "· 4 pairs"
    And value of first evidence row should have text "+4.14"
    When user clicks on open icon of first evidence row
    Then 5 rows should be selected
    And no errors should have been logged

  Scenario: Lower is better re-reads every ranking without a rebuild
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    And user remembers the "key series" reading of SAR Matrix Viewer viewer
    And user sets "Activity direction" property of SAR Matrix Viewer viewer to "Lower is better"
    Then value of Warhead component row should contain text "−2.45"
    And the "key series" reading of SAR Matrix Viewer viewer should not be as remembered
    And the "restored" reading of SAR Matrix Viewer viewer should be "true"
    When user sets "Activity direction" property of SAR Matrix Viewer viewer to "Higher is better"
    Then value of Warhead component row should contain text "+1.65"
    And the "key series" reading of SAR Matrix Viewer viewer should be as remembered
    And the "restored" reading of SAR Matrix Viewer viewer should be "true"
    And no errors should have been logged

  Scenario: The dialog offers -lg again after an activity it cannot scale
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    When user selects "Solubility logS (pred)" in Activity input in "SAR Matrix" dialog
    Then Scaling input in "SAR Matrix" dialog should be disabled
    And Scaling input in "SAR Matrix" dialog should have value "none"
    When user selects "Average Mass" in Activity input in "SAR Matrix" dialog
    Then Scaling input in "SAR Matrix" dialog should be enabled
    And Scaling input in "SAR Matrix" dialog should have value "-lg"
    When user clicks on CANCEL button in "SAR Matrix" dialog
    Then "SAR Matrix" dialog should be absent
    And no errors should have been logged

  Scenario: Another component across the columns rebuilds, a scaling it cannot take leaves none, and both go back at once
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    And user selects "E3 Ligand" in Columns input in SAR Matrix Viewer viewer
    Then SAR Matrix Viewer viewer should have built its matrices
    And the "axis" reading of SAR Matrix Viewer viewer should be "E3 Ligand"
    And the "matrices" reading of SAR Matrix Viewer viewer should be 71
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    Then setup line should contain text "columns: E3 Ligand"
    And Warhead component row should be visible
    When user sets "Scaling" property of SAR Matrix Viewer viewer to "lg"
    Then SAR Matrix Viewer viewer should have finished building
    And a warning balloon containing "cannot be scaled by lg" should have been shown
    And the "matrices" reading of SAR Matrix Viewer viewer should be 0
    And SAR Matrix note should have text "No SAR matrices found: no activity value can be scaled by lg. Set Scaling to none."
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    Then SAR Matrix note should have text "No SAR matrices found: no activity value can be scaled by lg. Set Scaling to none."
    When user clicks on SAR Transfer tab in SAR Matrix Viewer viewer
    Then SAR Matrix note should contain text "No SAR matrices to compare"
    When user sets properties of SAR Matrix Viewer viewer:
      | Scaling        | none    |
      | Matrix columns | Warhead |
    Then SAR Matrix Viewer viewer should have built its matrices
    And the "matrices" reading of SAR Matrix Viewer viewer should be 70
    And the "axis" reading of SAR Matrix Viewer viewer should be "Warhead"
    And no errors should have been logged
