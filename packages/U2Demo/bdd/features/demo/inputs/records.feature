@demo @realizes:u2.records-input @realizes:u2.checklist-input
Feature: Records and an ordered checklist
  A records input's rows are its items, keyed by their first field: Add field adds an empty row, whose empty
  name the input's verdict names by row until it is filled; a row's remove button takes it out. Stages move with
  Alt+↑/↓. A checklist's value is its checked items in the order shown.

  Background:
    Given user opens the "Records" demo page

  Scenario: A new row is named in the verdict until its name is filled, and goes with its remove button
    Then Fields records input should have 2 items in total
    And Fields records input should list the items "Target, Dose" in this order
    And value of problem readout should have text "none"
    When user clicks on "Add field" link
    Then Fields records input should have 3 items in total
    And value of problem readout should have text "Field 3: Name: Value can't be empty"
    When user types "Chemist" into last text input in Fields records input
    Then value of problem readout should have text "none"
    And Fields records input should list the items "Target, Dose, Chemist" in this order
    When user clicks on last Remove icon button in Fields records input
    Then Fields records input should have 2 items in total

  Scenario: Stages move with Alt+↑/↓
    When user presses Alt+ArrowDown in "Design" item in Stages records input
    Then value of stages readout should have text "Synthesis, Design, Assay"
    And Stages records input should list the items "Synthesis, Design, Assay" in this order

  Scenario: The checklist's value is its checked items in the order shown
    Then value of run readout should have text "Toxicity risks, Descriptors"
    And Run checklist should have 4 items in total
    And item "Toxicity risks" of Run checklist should be checked
    And item "Solubility" of Run checklist should be unchecked
    When user clicks on "Solubility" item in Run checklist
    Then value of run readout should have text "Toxicity risks, Descriptors, Solubility"
    When user presses Alt+ArrowUp in "Solubility" item in Run checklist
    Then value of run readout should have text "Toxicity risks, Solubility, Descriptors"
    And value of order readout should have text "Toxicity risks, Solubility, Descriptors, Structural alerts"
