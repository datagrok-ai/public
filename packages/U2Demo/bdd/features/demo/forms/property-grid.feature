@demo @realizes:u2.property-grid
Feature: The property grid
  A row is a "property" by its name; a "category" header folds its rows.

  Background:
    Given user opens the "Property grid" demo page

  Scenario: Editing a property replaces the value record
    Then editor of number input in opacity property should be as wide as editor of text input in title property
    And editor of number input in width property should be as wide as editor of text input in title property
    And value of "last change" readout should have text "(none)"
    When user types "Heatmap" into title property
    Then value of "last change" readout should have text "title → Heatmap"
    And value of values readout should contain text "\"title\":\"Heatmap\""
    When user unchecks showLegend property
    Then value of "last change" readout should have text "showLegend → false"
    When user selects "top" in position property
    Then value of "last change" readout should have text "position → top"
    When user enters "0.5" into opacity property
    Then value of "last change" readout should have text "opacity → 0.5"

  Scenario: Categories collapse
    Then Layout category should be expanded
    And width property should be visible
    When user collapses Layout category
    Then Layout category should be collapsed
    And width property should be hidden
    When user expands Layout category
    Then width property should be visible
