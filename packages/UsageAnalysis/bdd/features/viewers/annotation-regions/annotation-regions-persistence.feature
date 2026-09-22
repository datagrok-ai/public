@viewers @realizes:viewers.scatter-plot
Feature: Annotation regions persist with the view
  The regions are part of the viewer's look, so a layout saved from the view and loaded back
  brings them with it, drawn and titled. On a fresh view, since a journey that has swapped axes
  and renamed columns first cannot say which state its layout captured.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: The regions and their titles survive a layout round trip
    Given user adds a scatter plot viewer with:
      | xColumnName       | AGE    |
      | yColumnName       | WEIGHT |
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Outer","area":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{"type":"area","x":"AGE","y":"WEIGHT","header":"Inner","area":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}] |
    Then the "regions shown" reading of scatter plot viewer should be 2
    And scatter plot viewer should have a "region Outer title" area
    When user saves the layout of the current table view
    And user loads the saved layout
    Then the "viewer regions" reading of scatter plot viewer should be 2
    And the "regions shown" reading of scatter plot viewer should be 2
    And scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should have a "region Outer title" area
    And scatter plot viewer should have a "region Inner title" area
    And no errors should have been logged

  Scenario: The regions survive a layout saved through the server
    Given user adds a scatter plot viewer with:
      | xColumnName       | AGE    |
      | yColumnName       | WEIGHT |
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Outer","area":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]}] |
    Then the "regions shown" reading of scatter plot viewer should be 1
    When user saves the layout of the current table view to the server
    And user loads the saved layout
    Then the "viewer regions" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should have a "region Outer title" area
    And no errors should have been logged
