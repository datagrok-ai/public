@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Bio's sample FASTA file, opened from Browse (Files > App Data >
  Bio > samples > FASTA.fasta), becomes a table whose sequence column Datagrok recognizes as a Macromolecule in FASTA
  notation and draws residue by residue; the pointer on a residue names it; the context panel shows the column
  (notation, alphabet), then one sequence's composition. With the panel closed again, Bio > Transform > Convert
  Sequence Notation writes the column as separator text in a new column, separator(sequence), the same residues under
  the pointer; the dialog's list offers separator, HELM and BILN (HELM, drawn as a chain of small circles in a column
  this narrow, did not read on screen). Off camera, before the first step: the Browse tree opened down to
  Bio, the context panel closed. Nothing is made on the server. Merged from the former clips "Open sequences in a
  table" and "Convert between sequence notations".

  Scenario: Open, inspect and convert sequences
    Open a FASTA file as a table, see what Datagrok knows about the sequences, and write them in another notation.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Bio tree node inside browse tree is expanded
    And the context panel is closed
    # caption: Sample sequences: Files > App Data > Bio > samples
    When user expands Files---App-Data---Bio---samples tree node inside browse tree
    # caption: A FASTA file opens as a table
    And user double-clicks on Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree
    Then the "FASTA" view should be current
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "fasta"
    When the video says "Recognized as sequences: every residue in its color"
    # caption: Drag the column's edge to see more of each sequence
    And user drags the "column resizer sequence" area of grid by 200 pixels to the right
    And user hovers over monomer 6 of the "cell 2 of sequence" area of grid
    # caption: The context panel: what Datagrok knows about the column
    And user clicks on "Context Panel" status bar toggle
    And user clicks on the "header sequence" area of grid
    Then the context panel should show "sequence"
    # caption: Details: FASTA notation, an alphabet of peptides
    When user clicks on "Details" accordion header in context panel
    # caption: Click a sequence for panes about that one
    And user clicks on the "cell 3 of sequence" area of grid
    Then the context panel should show the current cell
    # caption: Composition: how often each residue occurs
    When user clicks on "Composition analysis" accordion header in context panel
    # caption: Close the panel to make room
    And user clicks on "Context Panel" status bar toggle
    # caption: Bio > Transform > Convert Sequence Notation
    And user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    Then "Convert Sequence Notation" dialog should be visible
    # caption: Pick a notation: here, separator text
    When user selects "separator" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then a new column "separator(sequence)" should have been added
    And "separator(sequence)" column should have units "separator"
    When the video says "A new column: the same residues, with dashes between them"
    And user hovers over monomer 6 of the "cell 2 of separator(sequence)" area of grid
    Then no errors should have been logged
