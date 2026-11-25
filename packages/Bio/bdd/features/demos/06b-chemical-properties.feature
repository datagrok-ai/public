@guide @demo
Feature: 2 · Every sequence is a molecule · Chemical properties
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). The first 63 peptides of the BILN_complex sample — disulfide-rich hormones and
  toxins, head-to-tail and side-chain rings, multi-chain insulins, lipidated branches, scaffolds — made molecules
  off camera (Bio | Transform | To Atomic Level, its defaults). Chem > Calculate > Chemical Properties then
  measures them as it measures small molecules (Crux on this stand, RDKit's definitions, in the browser): the
  molecular weight, ticked by default, the polar surface area and LogP, ticked here. A scatter plot of polar
  surface area against weight, colored by kind and on log axes, is set up off camera under a note, with the
  rule-of-five limits (MW 500, TPSA 140) drawn as lines: every peptide lies beyond them. A click on cyclosporin's
  marker, the classic oral drug beyond the rule of five, makes its row current, and the pointer leaves the plot.
  Set up off camera: the conversion, the columns, their formats, the plot and its place. Left out: the dialog's
  other calculators and properties.

  Scenario: Peptide drugs measured like small molecules
    Chem measures every peptide drug like a small molecule: weight, polar surface, LogP.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the Bio package is initialized
    And user opens BILN_complex dataset keeping the first 63 rows
    And the "biln" column is converted to atomic level
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "name, category, biln, molfile(biln)"
    And the "biln" column of the grid is 520 px wide
    And the "molfile(biln)" column of the grid is 360 px wide
    And the rows of the grid are 90 px tall
    When user picks "Chem > Calculate > Chemical Properties..." from the top menu
    Then "Chemical Properties" dialog should be visible
    And Molecules input in "Chemical Properties" dialog should contain text "molfile(biln)"
    When user checks "Chemical Properties" calculator in "Chemical Properties" dialog
    Then MW input in "Chemical Properties" dialog should be checked
    When user checks PSA input in "Chemical Properties" dialog
    And user checks "Log P" input in "Chemical Properties" dialog
    And user clicks on OK button in "Chemical Properties" dialog
    Then "Chemical Properties" dialog should be hidden
    And a new column "MW" should have been added
    And a new column "PSA" should have been added
    And a new column "LogP" should have been added
    And "MW" column should have no missing values
    And every value of "MW" column should lie between 500 and 10000
    Given user sets the format of "MW" column in the grid to "0"
    And user sets the format of "PSA" column in the grid to "0"
    And user sets the format of "LogP" column in the grid to "0.0"
    And the grid shows only the columns "name, category, MW, PSA, LogP, biln"
    And the monomer names of the "biln" column are cut at 8 characters
    When the video says "Weight, polar surface and LogP for every peptide" at the top
    Given the "biln" column of the grid is 300 px wide
    And the rows of the grid are 40 px tall
    And user adds a scatter plot viewer with:
      | xColumnName     | MW          |
      | yColumnName     | PSA         |
      | colorColumnName | category    |
      | xAxisType       | logarithmic |
      | yAxisType       | logarithmic |
      | xMin            | 300         |
      | yMin            | 100         |
      | legendPosition  | Top         |
    And user sets "formulaLines" property of scatter plot viewer to '[{"type":"line","formula":"${MW} = 500","title":"MW 500","color":"#e15759"},{"type":"line","formula":"${PSA} = 140","title":"TPSA 140","color":"#e15759"}]'
    And no row of the table is current
    When user clicks on the marker of row 35 of scatter plot viewer
    Then "name" of the current row should be "cyclosporin_like"
    When user moves the pointer away from scatter plot viewer
    And the video says "Cyclosporin and its peers: beyond the rule of five" at the top
    Then no error or warning balloon should have been shown
