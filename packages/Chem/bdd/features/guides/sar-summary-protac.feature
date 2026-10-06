@guide
Feature: SAR Matrix Summary on a PROTAC set
  The Summary tab reads a decomposition the table already holds. A degrader library is the case it
  was designed for: every compound is a warhead, a linker and an E3 ligand, and those columns mean
  the same thing in every row, so a component can be ranked across the whole set rather than only
  inside one series. Demo: 2,792 degraders from PROTAC-PatentDB over 71 linkers, 632 warheads and
  12 E3 ligands. The endpoint is predicted Caco-2 permeability, higher being more permeable, which
  is a property rather than a potency.

  One run ranks every component. A single additive fit over the core column and the R-group columns at
  once puts every value of every component on one scale, cross-validated and refitted on each half of
  the table before anything is ranked. Only the column axis (here E3 Ligand) also gets per-series
  evidence — outcome strips, matched pairs and differences against each series' own reference — and the
  Method segment offers the other columns as the axis of the next run.

  Scenario: E3 ligands over linkers
    Given user is logged in
    And simple mode is off
    And user opens protac-degraders dataset
    Given user watches the task bar
    When user calls "Chem:SarMatrixAnalysis" function with:
      | table               | table                             |
      | molecules           | column:Compound                   |
      | activity            | column:Caco-2 Permeability (pred) |
      | scaling             | none                              |
      | activityDirection   | Higher is better                  |
      | fragmentCutoff      | 0.4                               |
      | fragmentationLevels | 3                                 |
      | predictVirtual      | true                              |
      | useMcsAnchors       | false                             |
      | seriesColumn        |                                   |
      | coreColumn          | column:Linker                     |
      | fragmentColumns     | columns:Warhead,E3 Ligand         |
      | columnAxis          | E3 Ligand                         |
    Then the task bar should have finished "Building SAR matrices"
    When user clicks on "Summary" tab
    Then tab panel should contain text "What to change"
    And tab panel should contain text "Best measured swap"
    When user clicks on "Effects" summary segment
    And user clicks on "Linker" effects tab
    Then tab panel should contain text "Linker — offsets from the additive fit"
    When user clicks on "Measured in series" effects tab
    Then tab panel should contain text "Swapping E3 Ligand"
    When user clicks on "Method" summary segment
    Then tab panel should contain text "Trust gate"
    And no error or warning balloon should have been shown

  Scenario: every component in one run
    Given user is logged in
    And simple mode is off
    And user opens protac-degraders dataset
    Given user watches the task bar
    When user calls "Chem:SarMatrixAnalysis" function with:
      | table               | table                             |
      | molecules           | column:Compound                   |
      | activity            | column:Caco-2 Permeability (pred) |
      | scaling             | none                              |
      | activityDirection   | Higher is better                  |
      | fragmentCutoff      | 0.4                               |
      | fragmentationLevels | 3                                 |
      | predictVirtual      | true                              |
      | useMcsAnchors       | false                             |
      | seriesColumn        |                                   |
      | coreColumn          | column:Linker                     |
      | fragmentColumns     | columns:Warhead,E3 Ligand         |
      | columnAxis          | E3 Ligand                         |
    Then the task bar should have finished "Building SAR matrices"
    When user clicks on "Summary" tab
    And user clicks on "Effects" summary segment
    Then tab panel should contain text "moves Solubility logS (pred) most: E3 Ligand"
    And tab panel should contain text "E3 Ligand — offsets from the additive fit"
    And tab panel should contain text "Bimodal:"
    And tab panel should contain text "Cross-validated R²"
    When user clicks on "Warhead" effects tab
    Then tab panel should contain text "Warhead — offsets from the additive fit"
    When user clicks on "Linker" effects tab
    Then tab panel should contain text "Linker — offsets from the additive fit"
    When user clicks on "Method" summary segment
    Then tab panel should contain text "Every component is ranked on Effects"
    And no error or warning balloon should have been shown

  Scenario: warheads over linkers
    Given user is logged in
    And simple mode is off
    And user opens protac-degraders dataset
    Given user watches the task bar
    When user calls "Chem:SarMatrixAnalysis" function with:
      | table               | table                             |
      | molecules           | column:Compound                   |
      | activity            | column:Caco-2 Permeability (pred) |
      | scaling             | none                              |
      | activityDirection   | Higher is better                  |
      | fragmentCutoff      | 0.4                               |
      | fragmentationLevels | 3                                 |
      | predictVirtual      | true                              |
      | useMcsAnchors       | false                             |
      | seriesColumn        |                                   |
      | coreColumn          | column:Linker                     |
      | fragmentColumns     | columns:Warhead,E3 Ligand         |
      | columnAxis          | Warhead                           |
    Then the task bar should have finished "Building SAR matrices"
    When user clicks on "Summary" tab
    Then tab panel should contain text "What to change"
    When user clicks on "Effects" summary segment
    And user clicks on "Measured in series" effects tab
    Then tab panel should contain text "Swapping Warhead"
    And no error or warning balloon should have been shown
