@viewers @realizes:viewers.chrome
Feature: Viewer title and description
  The title and the description every viewer inherits from the base viewer: shown when set, placed
  where Description Position says, gone under Never, and cleared again. One outline over the
  viewers whose own features used to carry the same scenario; a viewer with a chrome of its own
  (the tile viewer's above/below content, the trellis plot's description slot) keeps its own.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario Outline: <viewer> shows and clears its title and description
    Given user adds a <viewer> viewer
    When user sets properties of <viewer> viewer:
      | Show Title | true         |
      | Title      | Demographics |
    Then title of <viewer> viewer should have text "Demographics"
    When user sets properties of <viewer> viewer:
      | Description                 | By race |
      | Description Visibility Mode | Always  |
    Then description of <viewer> viewer should have text "By race"
    When user sets "Description Position" property of <viewer> viewer to "Bottom"
    Then description of <viewer> viewer should be visible
    When user sets "Description Position" property of <viewer> viewer to "Left"
    Then description of <viewer> viewer should be visible
    When user sets "Description Position" property of <viewer> viewer to "Right"
    Then description of <viewer> viewer should be visible
    When user sets "Description Visibility Mode" property of <viewer> viewer to "Never"
    Then description of <viewer> viewer should be absent
    When user sets properties of <viewer> viewer:
      | Show Title                  | false |
      | Title                       |       |
      | Description                 |       |
      | Description Visibility Mode | Auto  |
      | Description Position        | Top   |
    Then description of <viewer> viewer should be absent
    And no errors should have been logged

    Examples:
      | viewer    |
      | bar chart |
      | box plot  |
      | histogram |
      | pc plot   |
