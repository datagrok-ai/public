@guide @help:visualize/viewers
Feature: Draw a custom annotation region on a scatter plot
  A guide: the answer to "how do I create custom annotation regions on a scatter plot?". A region
  is drawn by hand: the plot's context menu has Tools > Draw Annotation Region, which puts the plot
  into drawing mode, and dragging across the plot draws the rectangle. The Formula Lines dialog
  then opens on the new region (its title, color and opacity); OK keeps it. The region lives in the
  viewer's Annotation Regions property and is drawn on every repaint. Demo: demog-1000, AGE vs
  WEIGHT.

  Scenario: Draw an annotation region on a scatter plot
    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | WEIGHT |
    When user picks "Tools > Draw Annotation Region" from the context menu of scatter plot viewer
    Then the "region drawing mode" reading of scatter plot viewer should be "true"
    When user drags across the "view" area of scatter plot viewer
    Then scatter plot viewer should have a "region 1" area
    And the "region 1" area of scatter plot viewer should be painted
    When user clicks OK button in "Formula Lines" dialog
    Then the "regions shown" reading of scatter plot viewer should be 1
    And "annotationRegions" property of scatter plot viewer should contain "area"
