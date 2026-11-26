@guide @tutorial
Feature: How-to · Annotations
  For the maintainer: a Markup viewer as a note beside a chart, on demog-1000 with a scatter plot of WEIGHT against
  HEIGHT (added off camera). The Markup viewer comes from the Add viewer gallery (searched "markup"); its context
  menu's Edit content... opens the Edit dialog, where the note is typed in Markdown with column references, as a
  heading so it reads at video size: ## Patient ${USUBJID}: age ${AGE}, ${WEIGHT} kg, ${HEIGHT} cm. A reference shows the current row's value, so a
  click on the marker of row 11 (X0273T21000900001: 46, 118.90 kg, 168.974 cm, the heavy point on the right) makes
  that row current and the note shows that patient. demog-1000 is opened and the plot added off camera; the note
  lives in the view only, nothing is made on the server.

  Scenario: Add a note that explains a chart
    A Markup viewer puts text next to your charts, and the text can show values of the row you click.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    When the video says "A Markup viewer is a note next to your charts"
    # caption: Add viewer, then search for Markup
    And user clicks on "Add viewer" icon in toolbar
    And user types "markup" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Markup" card in "Add Viewer" dialog
    Then markup viewer should be visible
    # caption: Right-click it: Edit content...
    When user picks "Edit content..." from the context menu of markup viewer
    Then Edit dialog should be visible
    When user enters "## Patient ${USUBJID}: age ${AGE}, ${WEIGHT} kg, ${HEIGHT} cm" into input in Edit dialog
    # caption: Markdown text; ${COLUMN} shows the current row's value
    And user hovers over input in Edit dialog
    And user clicks on "OK" button in Edit dialog
    Then Edit dialog should be absent
    And "content" property of markup viewer should be "## Patient ${USUBJID}: age ${AGE}, ${WEIGHT} kg, ${HEIGHT} cm"
    # caption: Click a point: the note shows that patient
    When user clicks on the "marker of row 11" area of scatter plot viewer
    # caption: Age 46, 118.90 kg, 168.974 cm: the point you clicked
    Then row 11 should be current
    And markup viewer should contain text "Patient X0273T21000900001: age 46, 118.90 kg, 168.974 cm"
    And no errors should have been logged
