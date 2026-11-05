@realizes:bio.calculate.extract-region
Feature: The Get Region pane cuts a region out of a sequence column from the context panel
  A Macromolecule column's context panel carries Bioinformatics > Get Region (`Bio:getRegionPanel`): the
  Start, End and Column name inputs of the Extract Region dialog, proposed for the whole sequence, and a
  Get Region button that adds the region as a new sequence column (help `datagrok/solutions/domains/bio/bio.md`,
  "Get region": "expand the Get Region section. Fill in start and end positions of the region of interest, and
  name for the output column. A new column containing sequences of the region of interest is added"). The
  column is made current by a click on its header, as the Actions pane features do.

  Runs on the first nine sequences of filter_FASTA (38 positions; row 1 starts MDYKET), as
  transform/convert.feature does. Left out on purpose:
  - the Extract Region dialog itself: transform/convert.feature and transform/other-notations.feature;
  - the Region choice: it shows only on a column with annotated regions, and its regions are
    `src/tests/seq-handler-get-region-tests.ts`;
  - an MSA pane: none exists (MISSING.md, "Bio panes").
  Nothing is put on the server.

  Scenario: The pane proposes the whole range, and Get Region adds the chosen region as a sequence column
    Given user is logged in
    And user opens filter_FASTA dataset keeping the first 9 rows
    And the Bio package is initialized
    And the context panel is open
    When user clicks on the "header fasta" area of grid
    Then the context panel should show "fasta"
    When user expands Bioinformatics accordion header in context panel
    And user expands "Get Region" accordion header in context panel
    Then Start input in "Get Region" pane in context panel should have value "1"
    And End input in "Get Region" pane in context panel should have value "38"
    And "Column name" input in "Get Region" pane in context panel should have value "fasta: (1-38)"
    When user selects "3" in Start input in "Get Region" pane in context panel
    And user selects "6" in End input in "Get Region" pane in context panel
    Then "Column name" input in "Get Region" pane in context panel should have value "fasta: (3-6)"
    When user clicks on last "Get Region" button in "Get Region" pane in context panel
    Then the table should have a column "fasta: (3-6)"
    And "fasta: (3-6)" column should have semantic type "Macromolecule"
    And the value of "fasta: (3-6)" column in row 1 should be "YKET"
    And every value of "fasta: (3-6)" column should match "^[A-Z]{4}$"
    And no errors should have been logged
