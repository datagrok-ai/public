@journey @viewers @realizes:viewers.correlation-plot
Feature: Correlation plot — cell colours, the cell tooltip switch, the column labels and a selection
  What `correlation-plot-ui.md` and the remaining anchors of `correlation-plot.md` ask that the other
  three features do not: a negative coefficient is drawn blue and a positive one red, a near-zero one
  nearly white; the Tooltip > Visible menu item that mirrors Show Tooltip; a column with no variance
  hovers as "N/A" (GROK-12586); Col Header Font repaints the column labels (GROK-16818); and a
  selection drawn into the diagonal histograms (core#1087, fixed). The title and the
  description are every viewer's, in viewer-chrome.feature.
  Each colour is claimed twice: the cell's resolved background (`color of cell …`, measured) and the
  pixels inside the cell's region on the canvas, where the old spec sampled hard-coded offsets.
  Not translated: Grid > Grid Color Coding > Apply to text (GROK-19052) — the coloured text is drawn
  on the inner grid's canvas and nothing reports it, so the scenario would be a lone "no errors".

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a correlation plot viewer
    Then the "cells" reading of correlation plot viewer should be 16
    And the "rows shown" reading of correlation plot viewer should be 1000
    And correlation plot viewer should be painted

  Scenario: A negative coefficient is blue, a positive one red, and one near zero almost white
    Then the "correlation of HEIGHT and AGE" reading of correlation plot viewer should be between -0.2349 and -0.2348
    And the "correlation of WEIGHT and HEIGHT" reading of correlation plot viewer should be between 0.4124 and 0.4125
    And the "color of cell HEIGHT x AGE" reading of correlation plot viewer should be "#c3c3ff"
    And the "color of cell WEIGHT x HEIGHT" reading of correlation plot viewer should be "#ff9595"
    And the "color of cell STARTED x AGE" reading of correlation plot viewer should be "#fcfcff"
    And the "cell HEIGHT x AGE" area of correlation plot viewer should contain the color "#C3C3FF"
    And the "cell HEIGHT x AGE" area of correlation plot viewer should not contain the color "#FF9595"
    And the "cell WEIGHT x HEIGHT" area of correlation plot viewer should contain the color "#FF9595"
    And the "cell WEIGHT x HEIGHT" area of correlation plot viewer should not contain the color "#C3C3FF"
    And no errors should have been logged

  Scenario: Tooltip > Visible in the cell menu turns Show Tooltip off and back on, and the tooltip returns
    # After any context menu the grid shows no tooltip until the next press (the table's own grid
    # too), so the last menu pick is followed by a click on a cell before the hover. That the
    # tooltip stays away while Visible is off is not claimed: see MISSING.md.
    When user hovers over the "cell HEIGHT x AGE" area of correlation plot viewer
    Then tooltip should contain the text "Pearson R: -0.235"
    When user moves the pointer away from correlation plot viewer
    And user picks "Tooltip > Visible" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then "showTooltip" property of correlation plot viewer should be "false"
    When user picks "Tooltip > Visible" from the context menu of the "cell HEIGHT x AGE" area of correlation plot viewer
    Then "showTooltip" property of correlation plot viewer should be "true"
    When user clicks on the "cell STARTED x AGE" area of correlation plot viewer
    And user hovers over the "cell WEIGHT x HEIGHT" area of correlation plot viewer
    Then tooltip should contain the text "Pearson R: 0.412"
    When user moves the pointer away from correlation plot viewer
    Then no errors should have been logged

  Scenario: A column with no variance hovers as N/A instead of a number (GROK-12586)
    Given user adds a calculated column "FLAT" with formula "0"
    When user sets properties of correlation plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED, FLAT |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED, FLAT |
    Then the "cells" reading of correlation plot viewer should be 25
    When user hovers over the "cell FLAT x AGE" area of correlation plot viewer
    Then exactly one tooltip should be shown
    And tooltip should contain the text "N/A"
    And tooltip should not contain the text "NaN"
    When user moves the pointer away from correlation plot viewer
    And user sets properties of correlation plot viewer:
      | xColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
      | yColumnNames | AGE, HEIGHT, WEIGHT, STARTED |
    And user removes "FLAT" column
    Then the "cells" reading of correlation plot viewer should be 16
    And no errors should have been logged

  Scenario: Col Header Font redraws the column labels and moves no coefficient (GROK-16818)
    When user sets "colHeaderFont" property of correlation plot viewer to 'normal normal 18px "Roboto"'
    Then the "header HEIGHT" area of correlation plot viewer should have repainted
    And the "text of cell HEIGHT x AGE" reading of correlation plot viewer should be "-0.23"
    And the "cells" reading of correlation plot viewer should be 16
    When user sets "colHeaderFont" property of correlation plot viewer to 'normal normal 13px "Roboto"'
    Then the "header HEIGHT" area of correlation plot viewer should have repainted
    And no errors should have been logged

  Scenario: A selection is drawn into the diagonal histograms
    Then the "cell AGE x AGE" area of correlation plot viewer should not contain the color "#FF8C00"
    When user selects the first 20 rows
    Then 20 rows should be selected
    And the "cell AGE x AGE" area of correlation plot viewer should have repainted
    And the "cell AGE x AGE" area of correlation plot viewer should contain the color "#FF8C00"
    And no errors should have been logged
