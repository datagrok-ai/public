@journey @realizes:bio.render.msa-header
Feature: The WebLogo header of long sequence columns
  A long sequence column of a known alphabet gets a header above its cells: the conservation and
  the WebLogo of every position over a ruler that names the positions. The ruler and the WebLogo
  tooltip take the names from the column's position names, so a region extracted from an antibody
  numbering's aligned column is ruled in the scheme's positions — IMGT CDR3 from 105, not from 1 —
  whether Extract Region or the aligned cell's Extract CDR3 cut it out. Every long column of the
  table gets its header on its own.

  Not translated: the header's navigation along the positions (drag, wheel, arrow keys) and the
  Sequence Position Statistics viewer a click on a position docks; the ruler of the aligned column
  itself starts in the flank before the scheme's first position, and how long that flank is depends
  on the fixture's rows, so its names are claimed through the regions extracted from it.

  Background:
    Given user is logged in
    And user opens antibodies dataset keeping the first 40 rows as "antibodies"
    And the Bio package is initialized
    Then "AntibodyHC" column should have semantic type "Macromolecule"
    And "AntibodyLC" column should have semantic type "Macromolecule"

  Scenario: Both chains get the header, ruled from 1
    Then the "header tracks of AntibodyHC" reading of grid should be "Conservation, WebLogo"
    And the "header tracks of AntibodyLC" reading of grid should be "Conservation, WebLogo"
    And the "header positions of AntibodyHC" reading of grid should be "1, 10, 20"
    And the "header positions of AntibodyLC" reading of grid should be "1, 10, 20"

  Scenario: Extract Region cuts CDR3 out of the aligned column with the scheme's positions
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    And the "header tracks of AntibodyHC (aligned)" reading of grid should be "Conservation, WebLogo, Annotations"
    When user picks "Bio > Calculate > Extract Region..." from the top menu
    And user selects "AntibodyHC (aligned)" in Sequence input in "Get Sequence Region" dialog
    And user selects "CDR3: 105-117" in Region input in "Get Sequence Region" dialog
    And user enters "CDR3" into "Column name" input in "Get Sequence Region" dialog
    And user clicks on OK button in "Get Sequence Region" dialog
    Then a new column "CDR3" should have been added
    And "CDR3" column should have tag "aligned" equal to "SEQ.MSA"
    When user scrolls the mouse wheel down 5 times over the "row header 1" area of grid holding Shift
    Then the "header positions of CDR3" reading of grid should be "105, 110"
    When user hovers over the "position 110 of CDR3 header" area of grid
    Then tooltip should contain text "Position: 110"
    And no error or warning balloon should have been shown

  Scenario: The aligned cell's Extract CDR3 keeps the scheme's positions too
    When user scrolls the mouse wheel up 5 times over the "row header 1" area of grid holding Shift
    And user picks "Annotations > Extract CDR3 as Column" from the context menu of the "cell 1 of AntibodyHC (aligned)" area of grid
    Then a new column "AntibodyHC (aligned)(CDR3)" should have been added
    When user scrolls the mouse wheel down 5 times over the "row header 1" area of grid holding Shift
    Then the "header positions of AntibodyHC (aligned)(CDR3)" reading of grid should be "105, 110"
    And no error or warning balloon should have been shown
    And no errors should have been logged
