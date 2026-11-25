@footguns @batch-5
Feature: Foot-guns in the reaction enumerator's filters
  The enumerator drops a blocking SMARTS (a template row's blocking_fg) or an exclusion SMARTS that the engine
  does not read, without a word (enumerate.ts:213-218, filters.ts:182-194): a typo turns a filter off and the
  products it should have removed come out. So in Crux mode a SMARTS no engine reads is a block (Enumerate
  disabled, the reason beside it, as a dialog's footer would say it); a filter that matches nothing where it
  filters is a warn. Chem's demo exclusion list now writes its nitro `[N+](=O)[O-]`: the `[NX3](=O)=O` it had, a
  nitro written with a five-valent N, matches no nitro group as either engine holds one (0 of 648 nitro compounds
  in Chem's smi10K.csv, both engines), so a run with it must say so. Expected products are RDKit's on the given
  building blocks.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"
    And user opens the Reaction Enumerator from "Chem | Reactions"

  @catalogue-TRF-05
  Scenario: A blocking SMARTS no engine reads blocks Enumerate and says which
    Given the enumerator's templates are:
      | reaction_smarts                                                         | blocking_fg | reaction_name |
      | [C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]                          | [N;H2       | Amide         |
    And the enumerator's building blocks are "CC(=O)O, NCCN, CCN"
    Then the Enumerate button is disabled
    And the words beside the Enumerate button name "[N;H2" and say it cannot be read
    # the fix is written "[N&H2]", since the enumerator cuts "[N;H2]" at ";" (next scenario)
    When user changes the blocking SMARTS of template 1 to "[N&H2]"
    Then the Enumerate button is enabled

  @catalogue-TRF-06
  Scenario: An exclusion SMARTS no engine reads blocks Enumerate and says which
    Given the enumerator's templates are:
      | reaction_smarts                                                         | blocking_fg | reaction_name |
      | [C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]                          |             | Amide         |
    And the enumerator's building blocks are "CC(=O)O, CCN"
    And the enumerator's exclusion SMARTS are "C(=O)Cl], S(=O)(=O)[#6]"
    Then the Enumerate button is disabled
    And the words beside the Enumerate button name "C(=O)Cl]" and say it cannot be read
    And the words beside the Enumerate button do not name "S(=O)(=O)[#6]"

  @catalogue-TRF-07
  Scenario: A blocking SMARTS that matches none of the building blocks is a warn before the run
    Given the enumerator's templates are:
      | reaction_smarts                                                         | blocking_fg   | reaction_name |
      | [C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]                          | [NX3](=O)=O   | Amide         |
    And the enumerator's building blocks are "O=C(O)c1ccc([N+](=O)[O-])cc1, CCN"
    Then the Enumerate button is enabled
    And the words beside the Enumerate button warn that "[NX3](=O)=O" matches none of the building blocks

  @catalogue-UI-03 @catalogue-TRF-05 @catalogue-TRF-07
  Scenario: RDKit mode: the enumerator behaves as before
    Given the Chemical Engine of this page is "RDKit"
    And the enumerator's templates are:
      | reaction_smarts                                                         | blocking_fg | reaction_name |
      | [C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]                          | [N;H2       | Amide         |
    And the enumerator's building blocks are "CC(=O)O, NCCN, CCN"
    Then the Enumerate button is enabled
    And no foot-gun words show beside the Enumerate button
