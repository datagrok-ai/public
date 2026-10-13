@journey @realizes:hmm.domains.top-menu @realizes:hmm.search.top-menu @realizes:hmm.panels.anarci @realizes:hmm.panels.pfam
Feature: HMM: protein domains, profile search and the sequence panels
  Find Domains (HMMER) annotates the Pfam domains of a protein column from the package's built-in
  library and adds a column naming each row's domains: the antibody chains carry the V-set and C1-set
  immunoglobulin domains. Profile HMM Search searches the column with one family, V-set by default,
  and adds its score, E-value and significance columns, named after it. A cell of an antibody chain
  shows, on the context panel's Bioinformatics group, the ANARCI pane (chain, germlines, IMGT CDRs)
  and the Pfam Domains pane; it is read first, while the chain column is in view (the searches'
  columns scroll the grid to them). Translated from the HMM package's playwright/hmm-domains-search.test.ts, whose panel
  step called the panel functions directly; here the panes are opened from the cell.

  What the searches compute (which rows hold which domains, scores, E-values, clan competition) is
  HMM's own package test against the C programs (src/tests/search-tests.ts). Every model is the
  package's own built-in library; an accession outside it would be fetched from InterPro, so none is
  used. The table is the first 40 antibodies.

  Background:
    Given user is logged in
    And the "HMM" package is installed
    And user opens antibodies dataset keeping the first 40 rows as "antibodies"
    And the Bio package is initialized
    Then "AntibodyHC" column should have semantic type "Macromolecule"

  Scenario: An antibody chain shows its ANARCI and Pfam Domains panes
    Given the context panel is open
    When user clicks on the "cell 2 of AntibodyHC" area of grid
    Then the "current row" reading of grid should be 2
    When user expands Bioinformatics accordion header in context panel
    And user expands ANARCI accordion header in context panel
    Then "ANARCI" pane in context panel should contain the text "Heavy"
    And "ANARCI" pane in context panel should contain the text "IGHV"
    And "ANARCI" pane in context panel should contain the text "CDR3 (IMGT)"
    When user expands "Pfam Domains" accordion header in context panel
    Then "Pfam Domains" pane in context panel should contain the text "V-set"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Find Domains names the immunoglobulin domains of the chains
    When user picks "Bio > Annotate > Find Domains (HMMER)..." from the top menu
    Then "Find Domains (HMMER)" dialog should be visible
    And editor of Sequence input in "Find Domains (HMMER)" dialog should have text "AntibodyHC"
    When user clicks on OK button in "Find Domains (HMMER)" dialog
    Then a new column "AntibodyHC domains" should have been added
    And some value of "AntibodyHC domains" column should contain "V-set"
    And some value of "AntibodyHC domains" column should contain "C1-set"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Profile HMM Search with the built-in V-set family adds its columns
    When user picks "Bio > Search > Profile HMM Search..." from the top menu
    Then "Profile HMM Search" dialog should be visible
    When user clicks on OK button in "Profile HMM Search" dialog
    Then a new column "V-set significant" should have been added
    And a new column "V-set E-value" should have been added
    And a new column "V-set score" should have been added
    And no error or warning balloon should have been shown
    And no errors should have been logged
