@guide @demo @viewport:1440x810
Feature: Bioinformatics · Oligo Toolkit
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). SequenceTranslator's Oligo Toolkit, opened off camera: the Translator reads
  a 21-mer siRNA sense strand typed in Axolabs notation (2'-OMe lowercase, 2'-F "f", phosphorothioate "s") and
  lists it in BioSpring, HELM, Mermade12 and plain nucleotides, with its structure; the Pattern designer turns
  the example pattern to 2'-O-methyl and marks three sense positions 2'-fluoro (its translation example follows),
  its phosphorothioate linkages drawn as stars; the Structure tool links the sense strand to its antisense
  (typed 5' to 3', drawn 3' to 5' under the sense, base paired) and opens the duplex full size. Nothing is saved:
  the pattern is edited, never saved, so the account's patterns stay as they were.

  Scenario: One siRNA from notation to duplex
    One siRNA strand in every vendor notation, its modification pattern designed position by position, and the duplex drawn atom by atom.

    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens the "Oligo Toolkit" app on its "Oligo Translator" view
    When user types "csAfsgCfaUfuGfaAfgGfcUfaCfaAfusUfsu" into text area
    Then "Single sequence format" choice input should have the value "Axolabs"
    And "Nucleotides" table row should contain the text "CAGCAUUGAAGGCUACAAUUU"
    When the video says "One strand, every vendor notation"
    And user clicks on PATTERN tab
    And user selects "2'-O-Methyl" in "Sequence basis" choice input
    And user clicks on "Edit strands" button
    And user selects "2'-Fluoro" in "Sense modification 2" choice input in "Edit strands" dialog
    And user selects "2'-Fluoro" in "Sense modification 4" choice input in "Edit strands" dialog
    And user selects "2'-Fluoro" in "Sense modification 6" choice input in "Edit strands" dialog
    And user clicks on OK button in "Edit strands" dialog
    Then the "Edit strands" dialog should close
    When user clicks on STRUCTURE tab
    And user types "csAfsgCfaUfuGfaAfgGfcUfaCfaAfusUfsu" into SS text area
    And user types "asAfsaUfuGfuAfgCfcUfuCfaAfuGfcsUfsg" into AS text area
    And user clicks on the structure picture of the "Oligo Structure" view
    Then no error or warning balloon should have been shown
