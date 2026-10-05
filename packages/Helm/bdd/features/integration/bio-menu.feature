@journey @realizes:helm.cell.helm
Feature: The Bio menu on the HELM showcase leaves the Helm renderer intact
  The showcase mixes linear, cyclic and modified peptides, nucleic acids, CHEM and BLOB polymers and
  their conjugates in one HELM column. Each Bio command that takes a standalone HELM column runs on it
  and produces its result, and when the sweep is over the column is still a Macromolecule painted by
  the Helm renderer. Translated from TestTrack General/helm-bio-menu-integration (.md, -spec.ts).

  The commands' own behaviour on small HELM files is claimed in Bio's features (analyze, transform,
  search) and Dendrogram's; this is the showcase case: every notation in one column, one command after
  another on the same table.

  Not translated, and why: MSA — on a HELM column its only engine is PepSeA, a Docker container (the
  lead's rule; the dialog itself is claimed in Bio's analyze/msa-helm-dialog); Activity Cliffs — the
  showcase has no numeric column to take as the activity; Apply Numbering Scheme, Manage Annotations,
  PolyTool Convert and Enumerate HELM — the old spec claimed nothing for them, and Bio's annotate
  features and SequenceTranslator own them. The md's per-command "no Helm-related balloon" is claimed
  stricter, as no error or warning balloon of any kind; the md's tolerance for a non-Helm warning was
  not needed on the showcase. The cells of the search viewers' own result grids are not claimed
  (MISSING.md). The table has 53 rows, not the 55 of the old md and spec.

  Background:
    Given user is logged in
    And the Helm package is initialized
    And user opens helm-showcase dataset
    Then the table should have 53 rows
    And "HELM" column should have semantic type "Macromolecule"

  Scenario: Composition docks a WebLogo over the showcase column
    When user picks "Bio > Analyze > Composition" from the top menu
    Then the top menu command should have completed
    And WebLogo viewer should be visible
    And "Sequence Column Name" property of WebLogo viewer should be "HELM"
    And WebLogo viewer should be painted
    And WebLogo viewer should have a "position 1" area
    And the "rows shown" reading of WebLogo viewer should be 53
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on close icon of WebLogo viewer
    Then WebLogo viewer should be absent

  Scenario: Sequence Space embeds the showcase column
    When user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    And editor of Column input in "Sequence Space" dialog should have text "HELM"
    When user clicks on OK button in "Sequence Space" dialog
    Then the "Sequence Space" dialog should close
    And the top menu command should have completed
    And a new column "Embed_X_1" should have been added
    And a new column "Embed_Y_1" should have been added
    And "Embed_X_1" column should have no missing values
    And "Embed_X_1" column should have at least 10 distinct values
    And "Embed_Y_1" column should have no missing values
    And scatter plot viewer should be visible
    And "X" property of scatter plot viewer should be "Embed_X_1"
    And scatter plot viewer should show 53 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on close icon of scatter plot viewer
    And user removes "Embed_X_1" column
    And user removes "Embed_Y_1" column
    Then scatter plot viewer should be absent

  Scenario: Hierarchical Clustering attaches a tree with a leaf for every row
    Given user watches the task bar
    When user picks "Bio > Analyze > Hierarchical Clustering..." from the top menu
    Then "Hierarchical Clustering" dialog should be visible
    When user clicks on OK button in "Hierarchical Clustering" dialog
    Then "Hierarchical Clustering" dialog should be hidden
    And the task bar should have finished "Creating dendrogram"
    And the "tree leaves" reading of grid should be 53
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Convert Sequence Notation writes the column in the notation asked for
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    Then "Convert Sequence Notation" dialog should be visible
    And "Convert Sequence Notation" dialog should contain text "Current notation: helm"
    When user selects "separator" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then 1 new column should have been added
    And a new column "separator(HELM)" should have been added
    And "separator(HELM)" column should have units "separator"
    And the value of "separator(HELM)" column in row 2 should be "A-C-D-E-F-G-H-I-K-L"
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user removes "separator(HELM)" column

  Scenario: Extract Region cuts a HELM region
    When user picks "Bio > Calculate > Extract Region..." from the top menu
    Then "Get Sequence Region" dialog should be visible
    When user selects "1" in Start input in "Get Sequence Region" dialog
    And user selects "2" in End input in "Get Sequence Region" dialog
    And user enters "region 1-2" into "Column name" input in "Get Sequence Region" dialog
    And user clicks on OK button in "Get Sequence Region" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And "region 1-2" column should have units "helm"
    And the value of "region 1-2" column in row 5 should be "PEPTIDE1{A.A}$$$$"
    And the value of "region 1-2" column in row 2 should be "PEPTIDE1{A.C}$$$$"
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user removes "region 1-2" column

  Scenario: Scan Liabilities annotates the showcase column
    When user picks "Bio > Annotate > Scan Liabilities..." from the top menu
    Then "Scan Sequence Liabilities" dialog should be visible
    When user clicks on OK button in "Scan Sequence Liabilities" dialog
    Then a new column "~HELM_annotations" should have been added
    And some value of "~HELM_annotations" column should contain "oxid-m"
    And some value of "~HELM_annotations" column should contain "oxid-w"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Similarity Search lists the neighbours of a HELM row
    When user picks "Bio > Search > Similarity Search" from the top menu
    Then the top menu command should have completed
    And "Sequence Similarity Search" viewer should be visible
    And the "source column" reading of "Sequence Similarity Search" viewer should be "HELM"
    And the "neighbours" reading of "Sequence Similarity Search" viewer should be 11
    And the "target row" reading of "Sequence Similarity Search" viewer should be 0
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on close icon of "Sequence Similarity Search" viewer
    Then "Sequence Similarity Search" viewer should be absent

  Scenario: Diversity Search picks a varied HELM subset
    When user picks "Bio > Search > Diversity Search" from the top menu
    Then the top menu command should have completed
    And "Sequence Diversity Search" viewer should be visible
    And the "source column" reading of "Sequence Diversity Search" viewer should be "HELM"
    And the "subset size" reading of "Sequence Diversity Search" viewer should be 10
    And the "distinct sequences" reading of "Sequence Diversity Search" viewer should be 10
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user clicks on close icon of "Sequence Diversity Search" viewer
    Then "Sequence Diversity Search" viewer should be absent

  Scenario: Subsequence Search adds a filter on the HELM column
    When user picks "Bio > Search > Subsequence Search ..." from the top menu
    Then filters viewer should be visible
    And the filter panel should have 1 filter
    And the filter panel should have a filter on "HELM" column
    And the "type of HELM" reading of filter panel should be "Bio:bioSubstructureFilter"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Split to Monomers gives a Monomer column per position
    When user picks "Bio > Transform > Split to Monomers..." from the top menu
    Then "Split to Monomers" dialog should be visible
    And editor of Sequence input in "Split to Monomers" dialog should have text "HELM"
    When user clicks on OK button in "Split to Monomers" dialog
    Then the top menu command should have completed
    And "1" column should have semantic type "Monomer"
    And the value of "1" column in row 2 should be "A"
    And the value of "10" column in row 2 should be "L"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: After the sweep the column is still a HELM column painted by Helm
    When user moves the current cell of grid to the "HELM" column
    # the hues are read off row 2, so no current-row tint can supply them
    Then the "current row" reading of grid should be 1
    And the table should have 53 rows
    And "HELM" column should have semantic type "Macromolecule"
    And "HELM" column should have units "helm"
    And "HELM" column should have tag "cell.renderer" equal to "helm"
    And the "cell type of HELM" reading of grid should be "helm"
    And the "cell 2 of HELM" area of grid should be painted in at least 3 colors
    And no errors should have been logged
    And no error or warning balloon should have been shown
