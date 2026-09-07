@demo
Feature: Scoping and overrides
  "X in Y" resolves X inside Y; a registered whole phrase wins over composition; ordinals pick
  among several matches. The workbench vocabulary is a context that "user opens the MSA workbench"
  enters: its names and the generic kinds resolve inside the workbench first, while platform
  names such as "toolbox" keep their platform meaning.

  Background:
    Given user opens the MSA workbench

  Scenario: The same button name in two scopes
    When user clicks on save button in alignment panel
    Then notification should contain text "Saved from the form"
    When user clicks on save button in toolbar
    Then notification should contain text "Saved from the toolbar"

  Scenario: A registered whole phrase overrides composition
    When user clicks on save button inside toolbar
    Then notification should contain text "Saved from the toolbar"

  Scenario: Generic kinds inside a registered scope
    When user clicks on run button in toolbar
    Then MSA dialog should be visible
    When user closes MSA dialog
    Then MSA dialog should be hidden

  Scenario: Ordinals among the matches
    When user clicks on run button in toolbar
    And user clicks on OK button in MSA dialog
    Then second item in aligned sequences list should contain text "2"
    And last item in aligned sequences list should contain text "5"

  Scenario: A toggle reveals a panel
    Then settings panel should be hidden
    When user clicks on settings button in toolbar
    Then settings panel should be visible
    And theme input in settings panel should have value "light"

  Scenario: Platform names keep their platform meaning inside the workbench
    Then toolbar should be visible
    And MSA workbench should be visible
    And toolbox should be visible
    And browse tab should be visible
