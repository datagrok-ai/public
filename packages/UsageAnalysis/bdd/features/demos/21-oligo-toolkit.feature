@guide @demo
Feature: 3 · Oligonucleotides · Oligo Toolkit
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). SequenceTranslator's Oligo Toolkit, opened off camera: the Translator reads a
  21-mer siRNA sense strand typed in Axolabs notation (2'-OMe lowercase, 2'-F "f", phosphorothioate "s") and lists
  it in BioSpring, HELM, Mermade12 and plain nucleotides; the Pattern designer turns the example pattern to
  2'-O-methyl and marks two sense positions 2'-fluoro (its translation example follows), its phosphorothioate
  linkages drawn as stars. Left out for now: the structures — the Translator's structure box, hidden off camera,
  and the Structure tool's duplex — which this stand draws as two corner-to-corner lines (a Chem rendering
  regression since Nov 12; they drew right on Nov 9, and go back in once it is fixed). Nothing is saved: the
  pattern is edited, never saved, so the account's patterns stay as they were.

  Scenario: One siRNA from notation to pattern
    One strand in every vendor notation, then its modification pattern designed.

    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens the "Oligo Toolkit" app on its "Oligo Translator" view
    And the browse panel is closed
    And the structure box of the "Oligo Translator" view is hidden
    When user types "csAfsgCfaUfuGfaAfgGfcUfaCfaAfusUfsu" into text area
    Then "Single sequence format" choice input should have the value "Axolabs"
    And "Nucleotides" table row should contain the text "CAGCAUUGAAGGCUACAAUUU"
    When the video says "One strand, every vendor notation"
    And user clicks on PATTERN tab
    And user selects "2'-O-Methyl" in "Sequence basis" choice input
    And user clicks on "Edit strands" button
    And user selects "2'-Fluoro" in "Sense modification 2" choice input in "Edit strands" dialog
    And user selects "2'-Fluoro" in "Sense modification 4" choice input in "Edit strands" dialog
    And user clicks on OK button in "Edit strands" dialog
    Then the "Edit strands" dialog should close
    When the video says "2'-fluoro at positions 2 and 4"
    Then no error or warning balloon should have been shown
