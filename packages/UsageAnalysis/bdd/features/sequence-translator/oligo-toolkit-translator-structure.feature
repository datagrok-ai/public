@realizes:sequencetranslator.app.oligo-toolkit
Feature: Oligo Toolkit: the Translator and Structure tools
  The Oligo Toolkit app (Browse > Apps > Peptides > Oligo Toolkit) holds the TRANSLATOR, PATTERN and
  STRUCTURE tools as tabs, and the browser address follows the tab. The Translator detects the
  format of a typed oligonucleotide and lists its translations (a click on one copies it), marks
  in red the codes it cannot read and lists no translation for them, opens its monomer library as
  a read-only table, and converts a whole table column in bulk. The Structure tool refuses to save an SDF without a sense strand, and
  saves one with it. Translated from the TestTrack case
  SequenceTranslator/oligo-toolkit-translator-structure.

  The typed HELM of Blocks B and C is plain RNA, `RNA1{r(A)p.r(C)p.r(G)p.r(U)}$$$$`; the Nucleotides
  translation is read from the HELM branches themselves, so it does not depend on the monomer
  library. The GROK-20958 input of Block C is `RNA1{r(A)p.r(C)p.[meI]}$$$$` — detected as HELM, with
  the peptide monomer meI that the oligo library does not hold. The md's "no error at any
  intermediate state" is claimed on one partial HELM, `RNA1{r(A)p.r(`, typed on its own: the
  Translator reacts after a 300 ms pause, so a string typed at once is only processed as a whole.
  The prefilled Axolabs sample `Afcgacsu` (Nucleotides `ACGACU`) carries the detection and copy
  checks that pass on 1.11.5.

  GROK-20959 is fixed: the refusal names the sequence. It is claimed in one step, since a balloon claim takes
  the balloon it matched and a second claim on it finds none. GROK-20806 (a second conversion of the same
  column) is fixed: both columns stay, the second named "(2)"; Convert writes HELM.

  The translations themselves (the Axolabs sample's, the typed HELMs', none for an unreadable code), the
  values of the bulk conversion and the GROK-20958 refresh are SequenceTranslator's package tests
  'Formats support', 'Translator: bulk' and 'Translator: app'.

  The single-sequence format selector has no caption; it is named "Single sequence format". Not
  claimed: whether the structure picture drew a molecule — the picture is Chem's renderer, tested by
  the packages' own tests.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Apps tree node inside browse tree is expanded
    And Apps---Peptides tree node inside browse tree is expanded
    And Apps---Peptides---Oligo-Toolkit tree node inside browse tree is expanded
    When user double-clicks Apps---Peptides---Oligo-Toolkit---Oligo-Toolkit tree node inside browse tree
    Then the "Oligo Translator" view should be current
    And TRANSLATOR tab should be selected
    And text area should have the value "Afcgacsu"
    When user hovers over "Single sequence" heading
    Then tooltip should be hidden

  Scenario: The toolkit opens on the Translator, each tab builds its tool, and the address follows the tab
    Then "Single sequence" heading should be visible
    And "Bulk" heading should be visible
    And the page address should contain "/OligoToolkit/Translator"
    When user clicks on PATTERN tab
    Then "Load" heading should be visible
    And "Edit" heading should be visible
    And the page address should contain "/OligoToolkit/Pattern"
    When user clicks on STRUCTURE tab
    Then SS text area should be visible
    And AS text area should be visible
    And AS2 text area should be visible
    And "Save SDF" button should be visible
    And the page address should contain "/OligoToolkit/Structure"
    When user clicks on TRANSLATOR tab
    Then "Single sequence" heading should be visible
    And the page address should contain "/OligoToolkit/Translator"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The prefilled Axolabs sequence is detected as Axolabs
    Then "Single sequence format" choice input should have the value "Axolabs"
    And no errors should have been logged

  Scenario: A click on a translation copies it
    When user clicks on "ACGACU" link
    Then an info balloon containing "Copied" should have been shown
    And the clipboard should have the text "ACGACU"

  Scenario: A typed HELM switches the format selector to HELM and is translated to nucleotides without an error (GROK-20958)
    When user types "RNA1{r(A)p.r(C)p.r(G)p.r(U)}$$$$" into text area
    Then "Single sequence format" choice input should have the value "HELM"
    When user clicks on "ACGU" link
    Then an info balloon containing "Copied" should have been shown
    And the clipboard should have the text "ACGU"
    And no errors should have been logged

  Scenario: A HELM with a monomer the oligo library lacks refreshes the translations without an error, also while typed (GROK-20958, GROK-19926)
    When user types "RNA1{r(A)p.r(" into text area
    Then "Single sequence format" choice input should have the value "HELM"
    And "Nucleotides" table row should not contain the text "ACGACU"
    And no errors should have been logged
    When user types "RNA1{r(A)p.r(C)p.r(G)p.r(U)}$$$$" into text area
    Then "Nucleotides" table row should contain the text "ACGU"
    When user types "RNA1{r(A)p.r(C)p.[meI]}$$$$" into text area
    Then "Single sequence format" choice input should have the value "HELM"
    And "Nucleotides" table row should not contain the text "ACGU"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Codes the Translator cannot read are marked red and leave no translations
    When user types "AfcgacsuQQ" into text area
    Then "Single sequence format" choice input should have the value "Axolabs"
    And the "-webkit-text-fill-color" style of "QQ" text should contain "rgb(255, 0, 0)"
    And "Nucleotides" table row should be absent
    When user types "Afcgacsu" into text area
    Then "Nucleotides" table row should contain the text "ACGACU"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The Oligo Translator's monomer library opens as a read-only table
    When user double-clicks Apps---Peptides---Oligo-Toolkit---Oligo-Translator tree node inside browse tree
    Then the "Oligo Translator" view should be current
    When user clicks on "View monomer library" icon
    Then the "Monomer Library" view should be current
    And "symbol" column should have at least 10 distinct values
    And "molfile" column should have semantic type "Molecule"
    When user double-clicks on the "cell 1 of symbol" area of grid
    Then a warning balloon containing "read-only" should have been shown
    And cell editor should be hidden
    And no errors should have been logged

  Scenario: Bulk converts a table column to Nucleotides and to HELM
    Given Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded
    When user double-clicks Files---App-Data---SequenceTranslator---samples---bulk-translation-axolabs.csv tree node inside browse tree
    Then the "bulk-translation-axolabs" view should be current
    When user clicks on the tab of the "Oligo Translator" view
    Then Table input should have the value "bulk-translation-axolabs"
    And Sequence input should have the value "AxolabsSequences"
    And "Input format" input should have the value "Axolabs"
    And "Output format" input should have the value "Nucleotides"
    When user clicks on Convert button
    Then the "bulk-translation-axolabs" view should be current
    And the table should have a column "AxolabsSequences (Nucleotides)"
    And "AxolabsSequences (Nucleotides)" column should have semantic type "Macromolecule"
    When user clicks on the tab of the "Oligo Translator" view
    When user selects "HELM" in "Output format" input
    And user clicks on Convert button
    Then the "bulk-translation-axolabs" view should be current
    And the table should have a column "AxolabsSequences (HELM)"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Converting the same column twice keeps both columns under different names (GROK-20806)
    Given Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded
    When user double-clicks Files---App-Data---SequenceTranslator---samples---bulk-translation-axolabs.csv tree node inside browse tree
    Then the "bulk-translation-axolabs" view should be current
    When user clicks on the tab of the "Oligo Translator" view
    When user clicks on Convert button
    Then the table should have a column "AxolabsSequences (HELM)"
    When user clicks on the tab of the "Oligo Translator" view
    When user clicks on Convert button
    And user clicks on the tab of the "bulk-translation-axolabs" view
    Then the table should have the columns "AxolabsSequences, AxolabsSequences (HELM), AxolabsSequences (HELM) (2)"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Save SDF without a sense strand is refused, and with one saves an SDF record
    Given user watches downloads
    When user clicks on STRUCTURE tab
    And user types "Afcgacsu" into AS text area
    And user clicks on "Save SDF" button
    Then a warning balloon containing "Enter SENSE_STRAND and optionally ANTISENSE_STRAND/AS2 to save SDF" should have been shown
    And no file should have been downloaded
    When user clears AS text area
    And user types "Afcgacsu" into SS text area
    And user downloads a file through "Save SDF" button
    Then a file matching "^SequenceTranslator-\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}\.sdf$" should have been downloaded
    And the downloaded file should contain 1 occurrences of "$$$$"
    And the downloaded file should contain "M  END"
    And no errors should have been logged

  Scenario: Save SDF with an antisense strand that cannot be converted is refused by name (GROK-20959)
    Given user watches downloads
    When user clicks on STRUCTURE tab
    And user types "Afcgacsu" into SS text area
    And user types "NOTASEQUENCE" into AS text area
    And user clicks on "Save SDF" button
    Then a warning balloon containing "Unable to save SDF: Unable to detect the format of the sequence 'NOTASEQUENCE'" should have been shown
    And no file should have been downloaded
    And no errors should have been logged
