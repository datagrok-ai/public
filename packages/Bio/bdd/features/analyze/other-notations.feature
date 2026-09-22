@realizes:bio.analyze.sequence-space @realizes:bio.analyze.activity-cliffs @realizes:bio.analyze.composition
Feature: The Analyze commands on HELM and MSA columns
  Sequence Space, Activity Cliffs and Composition take a HELM or an aligned (MSA) column as well
  as a fasta one; a HELM column is split into monomers by the HELM parser and a separator column
  by its separator before anything is computed, so each notation is its own path. The fasta runs,
  with their dialogs and their edited reruns, are in sequence-space, activity-cliffs and
  composition; these are the other two notations. A HELM run takes longer than the command's
  own completion: the dialog stays open, its OK greyed, until the embedding is written, so the
  claims start when the dialog closes.

  Not translated, and why: the manual cases name the tests/filter_* files for every command, but
  filter_HELM has no activity column (so Activity Cliffs cannot start on it) and four rows (too
  few for an embedding), so the embeddings and the cliffs run on the first 100 rows of the
  package's samples/HELM.csv and samples/MSA.csv, which carry Activity. Composition, which needs
  no activity, runs on the filter files the cases name.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario Outline: Sequence Space embeds a <notation> column
    Given user opens <dataset> dataset keeping the first 100 rows
    Then "<column>" column should have units "<units>"
    When user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    And editor of Column input in "Sequence Space" dialog should have text "<column>"
    When user clicks on OK button in "Sequence Space" dialog
    Then the "Sequence Space" dialog should close
    And the top menu command should have completed
    And a new column "Embed_X_1" should have been added
    And a new column "Embed_Y_1" should have been added
    And "Embed_X_1" column should have no missing values
    And "Embed_X_1" column should have at least 10 distinct values
    And scatter plot viewer should be visible
    And "X" property of scatter plot viewer should be "Embed_X_1"
    And scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    And no errors should have been logged

    Examples:
      | notation | dataset     | column | units     |
      | HELM     | HELM_sample | HELM   | helm      |
      | MSA      | MSA_sample  | MSA    | separator |

  Scenario Outline: Activity Cliffs finds cliffs in a <notation> column
    Given user opens <dataset> dataset keeping the first 100 rows
    When user picks "Bio > Analyze > Activity Cliffs..." from the top menu
    Then "Sequence Activity Cliffs" dialog should be visible
    And editor of Column input in "Sequence Activity Cliffs" dialog should have text "<column>"
    When user selects "Activity" in Activities input in "Sequence Activity Cliffs" dialog
    And user clicks on OK button in "Sequence Activity Cliffs" dialog
    Then the "Sequence Activity Cliffs" dialog should close
    And the top menu command should have completed
    And a new column "Embed_X_1" should have been added
    And a new column matching "sali|SALI" should have been added
    And scatter plot viewer should be visible
    And title of scatter plot viewer should have text "Activity cliffs"
    And the activity cliffs plot should report at least 1 cliff
    And scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    And no errors should have been logged

    Examples:
      | notation | dataset     | column |
      | HELM     | HELM_sample | HELM   |
      | MSA      | MSA_sample  | MSA    |

  Scenario: Composition docks a WebLogo over a HELM column
    Given user opens filter_HELM dataset
    When user picks "Bio > Analyze > Composition" from the top menu
    Then the top menu command should have completed
    And WebLogo viewer should be visible
    And WebLogo viewer should be bound to table "filter_HELM"
    And "Sequence Column Name" property of WebLogo viewer should be "HELM string"
    And WebLogo viewer should be painted
    And WebLogo viewer should have a "position 1" area
    And the "rows shown" reading of WebLogo viewer should be 4
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Composition on an aligned column selects the rows of a multi-letter monomer
    Given user opens filter_MSA dataset
    Then "MSA" column should have units "separator"
    When user picks "Bio > Analyze > Composition" from the top menu
    Then the top menu command should have completed
    And WebLogo viewer should be visible
    And "Sequence Column Name" property of WebLogo viewer should be "MSA"
    And WebLogo viewer should be painted
    And WebLogo viewer should have a "monomer meI at position 1" area
    And the "positions shown" reading of WebLogo viewer should be at least 15
    And no rows should be selected
    When user clicks on the "monomer meI at position 1" area of WebLogo viewer
    Then some rows should be selected
    And only rows where "MSA" starts with "meI/" should be selected
    And no error or warning balloon should have been shown
    And no errors should have been logged
