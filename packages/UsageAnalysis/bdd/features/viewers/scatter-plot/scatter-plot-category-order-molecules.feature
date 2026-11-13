@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot molecule categories ordered by a computed property
  A categorical axis of molecules is ordered by a property a package computes for them: Aggregated
  value on the molecule column itself offers Chem's Molecular Property choices (MW first), and the
  markers move to their core's place by molecular weight once the values arrive. Ascending MW is the
  claim that the computed order is not the by-name one (descending MW gives the by-name places; the
  direction is the Dart test `d4 | Axes | Semantic Category Order`).
  The plot reports a render pending while Chem computes the values (about a second for the seven
  cores), so the claims after the pick wait for them.
  spgi-100 with Core (six scaffolds with R-group attachment points, seven rows without one) on X and
  Series on Y, so the dialog starts on the first column it can order by, Structure. Row 8 holds the
  triazole core (the lightest), row 1 CC1(N…)CCN(…)C1 and row 2 the heavier CC1(CN…)CCN(…)C1; by
  name they come row 2, row 1, row 8. The numbers the property gives are Chem's package test
  (`Molecular Property`) and the Dart test `d4 | Axes | Semantic Category Order`. One journey; the
  scenario puts the order back to Default.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And user opens spgi dataset
    And user adds a scatter plot viewer with:
      | X | Core   |
      | Y | Series |
    Then the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 8" area

  Scenario: Molecular weight orders the cores
    When user picks "Order..." from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "Order categories: Core" dialog should be visible
    When user selects "Aggregated value" in "Order by" input in "Order categories: Core" dialog
    Then "Column" input in "Order categories: Core" dialog should contain text "Structure"
    When user picks "Core" in the column selector of the "Order categories: Core" dialog
    Then "Property" input in "Order categories: Core" dialog should be visible
    And "Property" input in "Order categories: Core" dialog should have value "MW"
    And "Aggregation" input in "Order categories: Core" dialog should be hidden
    When user selects "Ascending" in "Direction" input in "Order categories: Core" dialog
    Then the "marker of row 8" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And the "marker of row 1" area of scatter plot viewer should lie to the left of the "marker of row 2" area
    When user selects "Default" in "Order by" input in "Order categories: Core" dialog
    And user clicks on "CLOSE" button in "Order categories: Core" dialog
    Then the "marker of row 2" area of scatter plot viewer should lie to the left of the "marker of row 1" area
    And no error or warning balloon should have been shown
    And no errors should have been logged
