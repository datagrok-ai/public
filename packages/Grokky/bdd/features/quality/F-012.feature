@serial
Feature: Ordering a new viewer of a column leaves the other viewers of that column in their order
  pilot-1120 F-012: with a bar chart drawing a column's categories heavy, medium, light (custom sort,
  following the column's category order), the user asked for a box plot of the same column light,
  medium, heavy. The assistant reordered the shared column, so the bar chart flipped unasked, or switched
  the bar chart to sort by value. The box plot has its own categoryOrder ({"manual": [...]}), which leaves
  the bar chart alone. Both checks read the categories as drawn, from the viewers' hit areas; the counts
  (light 4, heavy 3, medium 2) keep a sort by value from drawing heavy, medium, light.

  These scenarios reach a Docker container and a paid model: they run only with GROKKY_BDD_LIVE=1.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens a table "patients" with:
      | weight class | age |
      | light        | 25  |
      | light        | 31  |
      | light        | 44  |
      | light        | 52  |
      | heavy        | 38  |
      | heavy        | 61  |
      | heavy        | 47  |
      | medium       | 29  |
      | medium       | 55  |
    And the AI panel is open on a new chat

  Scenario: A box plot ordered light first keeps the bar chart's heavy-first order
    Given the categories of "weight class" column are ordered "heavy, medium, light"
    And user adds a bar chart viewer with:
      | splitColumnName | weight class |
      | valueColumnName | age          |
      | valueAggrType   | count        |
      | barSortType     | custom       |
      | barSortOrder    | asc          |
    And the "Bar chart" viewer should draw its categories in the order "heavy, medium, light"
    When user asks the AI panel "Now also a box plot of age by weight class, light first, then medium, then heavy."
    Then the "Box plot" viewer should draw its categories in the order "light, medium, heavy"
    And the "Bar chart" viewer should draw its categories in the order "heavy, medium, light"
    And a screenshot "F-012-box-and-bar-order" of the AI panel is saved
