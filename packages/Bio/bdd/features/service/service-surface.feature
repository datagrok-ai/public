@journey @realizes:bio.service-surface
Feature: The service surface other packages call
  Helm, Peptides, Dendrogram and others start by calling Bio's service getters; the platform
  holds those calls until Bio has initialized, so they resolve to the initialized singletons —
  never before, never to a stub. A sequence handler is per column and knows its notation.

  Background:
    Given user is logged in

  Scenario: The singletons resolve with the methods their consumers use
    When user calls "Bio:getSeqHelper" function
    Then the result should have methods "getSeqHandler, getSeqMonomers, helmToAtomicLevel, setUnitsToFastaColumn"
    When user calls "Bio:getMonomerLibHelper" function
    Then the result should have methods "getMonomerLib, awaitLoaded"
    When user calls "Bio:getBioLib" function
    Then the result should have methods "getMonomer, getMonomerSymbolsByType, getPolymerTypes"
    And no error or warning balloon should have been shown

  Scenario: A sequence handler is per column and reports the column's notation
    Given user opens filter_HELM dataset
    Then "HELM string" column should have units "helm"
    When user calls "Bio:getSeqHandler" function with:
      | sequence | column:HELM string |
    Then the result should have methods "getSplitter, getRegion, convert"
    And the result should have a "notation" of "helm"
    Given user opens filter_FASTA dataset
    When user calls "Bio:getSeqHandler" function with:
      | sequence | column:fasta |
    Then the result should have a "notation" of "fasta"
    And no errors should have been logged

  Scenario: The HELM monomer list is exactly the monomers of the column
    When user switches to the "filter_HELM" table view
    And user calls "Bio:getHelmMonomers" function with:
      | sequence | column:HELM string |
    Then the result should be a list of 1 or more items
    And the result should be exactly the monomers of "HELM string" column
    And no errors should have been logged
