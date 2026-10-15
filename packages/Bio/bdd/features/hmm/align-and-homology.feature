@journey @realizes:hmm.align.msa-engine @realizes:hmm.align.bio-dialog @realizes:hmm.homology.top-menu @realizes:hmm.logo.weblogo @realizes:hmm.build.top-menu
Feature: HMM: profile alignment, homology search, profile logos and built profiles from the Bio menu
  The HMM package (HMMER reimplemented in the browser) adds its commands to Bio's menu and its engine
  to Bio's MSA dialog. On a protein column the dialog offers Kalign first and the profile HMM
  aligner (hmmalign) beside it, which aligns the chains to the best antibody germline model; Homology
  Search (phmmer) takes the current row as its query; Profile HMM Logo docks Bio's WebLogo drawing a
  family of the built-in Pfam library; Build Profile HMM builds a profile from the selected rows of
  the alignment, searches the table with it and shows its logo, whose letters select the rows of the
  alignment that have them: the valine of position 2 is in 35 of the 40 chains. Translated from the
  HMM package's playwright/hmm-align-homology.test.ts.

  What the programs compute (the alignment, scores and E-values, logo heights, the built model) is
  HMM's own package test against the C programs (src/tests/tier2-tests.ts, search-tests.ts); the
  feature claims the dialogs and what each OK adds. Every model comes with the package: the built-in
  Pfam library and the germline models; a Pfam accession outside the library would be fetched from
  InterPro, so none is used. The table is the first 40 antibodies.

  Background:
    Given user is logged in
    And the "HMM" package is installed
    And user opens antibodies dataset keeping the first 40 rows as "antibodies"
    And the Bio package is initialized
    Then "AntibodyHC" column should have semantic type "Macromolecule"

  Scenario: The MSA dialog offers hmmalign beside Kalign and aligns the chains to the germline profile
    When user picks "Bio > Analyze > MSA..." from the top menu
    Then MSA dialog should be visible
    And Engine input in MSA dialog should have value "Kalign"
    And Engine input in MSA dialog should offer the choice "Profile HMM (hmmalign)"
    When user selects "Profile HMM (hmmalign)" in Engine input in MSA dialog
    And user selects "Antibody/TCR germline (auto)" in Profile input in MSA dialog
    And user clicks on OK button in MSA dialog
    Then MSA dialog should be hidden
    And a new column "msa(AntibodyHC)" should have been added
    And "msa(AntibodyHC)" column should have units "fasta"
    And "msa(AntibodyHC)" column should have tag "aligned" equal to "SEQ.MSA"
    And every value of "msa(AntibodyHC)" column should have the same length
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Homology Search takes the current row as its query
    When user makes row 3 current
    And user picks "Bio > Search > Homology Search (phmmer)..." from the top menu
    Then "Homology Search (phmmer)" dialog should be visible
    And Row input in "Homology Search (phmmer)" dialog should have value "3"
    When user clicks on OK button in "Homology Search (phmmer)" dialog
    Then a new column "phmmer row 3 E-value" should have been added
    And a new column "phmmer row 3 score" should have been added
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Profile HMM Logo docks a WebLogo of a family of the built-in library
    When user picks "Bio > Analyze > Profile HMM Logo..." from the top menu
    Then "Profile HMM Logo" dialog should be visible
    When user selects "Pfam family (built-in)" in "Profile from" input in "Profile HMM Logo" dialog
    And user selects "C1-set" in Family input in "Profile HMM Logo" dialog
    And user clicks on OK button in "Profile HMM Logo" dialog
    Then WebLogo viewer should be visible
    And the "positions shown" reading of WebLogo viewer should be at least 20
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Build Profile HMM from the selected rows searches the table and shows a logo whose letters select rows
    When user selects the first 15 rows
    And user picks "Bio > Search > Build Profile HMM..." from the top menu
    Then "Build Profile HMM" dialog should be visible
    And Rows input in "Build Profile HMM" dialog should have value "Selected rows"
    When user enters "heavy15" into Name input in "Build Profile HMM" dialog
    And user clicks on OK button in "Build Profile HMM" dialog
    Then a new column "heavy15 E-value" should have been added
    # the family logo of the scenario before is the first
    And second WebLogo viewer should be visible
    When user selects no rows
    # caption: Click the valine at position 2 of the profile
    And user clicks on the "monomer V at position 2" area of second WebLogo viewer
    Then 35 rows should be selected
    And no error or warning balloon should have been shown
    And no errors should have been logged
