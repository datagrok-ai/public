@realizes:bio.analyze.msa
Feature: The MSA dialog on a HELM column, before the PepSeA engine runs
  On a HELM column Bio | Analyze | MSA... offers the non-canonical engine, PepSeA, with its own gap
  penalties behind the Alignment parameters toggle, and a cluster column to align each cluster on its
  own; the same dialog is offered on a HELM table that came back from a saved project. Translated from
  TestTrack General/pepsea-spec.ts (and its Bio copy) and General/bio-lifecycle-pepsea-container-spec.ts,
  the dialog parts analyze/msa does not claim (the cluster column and Gap Open on a HELM column, the
  dialog after a project round trip).

  Not translated, and why: pressing OK — on a HELM column it runs PepSeA in its Docker container (the
  lead's rule), so the aligned column, its renderer, the per-cluster widths and the container's status
  and eviction checks of both specs are out; the old random cluster column (RandBetween over four rows,
  which can give one cluster) is a deterministic length parity here. The lifecycle spec's claim of a
  second, kalign engine is stale: PepSeA is the only engine registered for non-canonical sequences, and
  kalign is the dialog's other mode, for canonical columns. The JS API checks of the lifecycle spec
  (initBio, getSeqHelper, saving the project through the API) are not UI.

  The project the second scenario saves is removed now and at feature end.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: The cluster column and both PepSeA gap penalties are offered on a HELM column
    Given user opens filter_HELM dataset
    When user adds a calculated column "Clusters" with formula "Length(${HELM string}) % 2"
    Then "Clusters" column should have type "int"
    And "Clusters" column should have at least 2 distinct values
    When user picks "Bio > Analyze > MSA..." from the top menu
    Then MSA dialog should be visible
    And Engine input in MSA dialog should have value "PepSeA"
    And Clusters input in MSA dialog should be visible
    When user selects "Clusters" in Clusters input in MSA dialog
    Then editor of Clusters input in MSA dialog should have text "Clusters"
    And "Gap Open" input in MSA dialog should have value "1.53"
    And "Gap Extend" input in MSA dialog should have value "0"
    When user clicks on "Alignment parameters" button in MSA dialog
    Then "Gap Open" input in MSA dialog should be hidden
    And "Gap Extend" input in MSA dialog should be hidden
    When user clicks on "Alignment parameters" button in MSA dialog
    Then "Gap Open" input in MSA dialog should be visible
    And "Gap Extend" input in MSA dialog should be visible
    When user clicks on CANCEL button in MSA dialog
    Then MSA dialog should be hidden
    And no new column should have been added
    And no errors should have been logged

  Scenario: A HELM table reopened from a project offers the same engine
    Given the user's own project "bdd-bio-msa-helm-{run}" is removed now and at feature end
    And user opens filter_HELM dataset
    When user saves the current view as project "bdd-bio-msa-helm-{run}"
    Then 1 project named "bdd-bio-msa-helm-{run}" should be on the server
    When user closes all views
    And user opens the "bdd-bio-msa-helm-{run}" project and waits for its table
    Then "HELM string" column should have units "helm"
    When user picks "Bio > Analyze > MSA..." from the top menu
    Then MSA dialog should be visible
    And editor of Sequence input in MSA dialog should have text "HELM string"
    And Engine input in MSA dialog should have value "PepSeA"
    And Method input in MSA dialog should have value "mafft --auto"
    When user clicks on CANCEL button in MSA dialog
    Then MSA dialog should be hidden
    And no new column should have been added
    And no errors should have been logged
