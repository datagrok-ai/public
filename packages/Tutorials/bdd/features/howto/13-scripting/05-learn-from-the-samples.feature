@guide @tutorial @viewport:1440x810
Feature: How-to · Scripting
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The JavaScript samples that come with the platform (ApiSamples):
  Browse > Platform > Functions > JS API Samples, grid > color-coding > color-coding opens as an unsaved script;
  its height line is changed (green up to 160 cm instead of 170) and the script is run: it opens a 10,000-row demo
  table ("demog 10000", random each run) with height colored by those ranges, age on a red-yellow-green scale,
  weight on its own scale and two races in their own colors. Starts with Browse open and the toolbox and the
  context panel hidden. Nothing is made on the server: the sample is not saved, its table lives in the page. Filmed
  in a 1440x810 window, so the code reads larger.

  Scenario: Learn from the sample scripts
    Datagrok comes with short JavaScript samples for common tasks. Open one, change it, run it.

    Given user is logged in
    And simple mode is off
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    # caption: JS API Samples: short working scripts, by topic
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Functions" tree node inside browse tree
    And user expands "Platform > Functions > JS API Samples" tree node inside browse tree
    And user expands "Platform > Functions > JS API Samples > grid" tree node inside browse tree
    And user expands "Platform > Functions > JS API Samples > grid > color-coding" tree node inside browse tree
    # caption: A sample opens as a script you can change and run
    And user clicks on "Platform > Functions > JS API Samples > grid > color-coding > color-coding" tree node inside browse tree
    Then code editor should contain the text "t.col('height').meta.colors.setConditional"
    # caption: It colors the columns of a demo table, one line per column
    When user hovers over the text "// Conditional color-coding for numerical columns" in code editor
    # caption: Change it: green up to 160 cm, dark from 160 to 190
    And user replaces the line starting with "t.col('height')" in code editor with "t.col('height').meta.colors.setConditional({'20-160': '#00FF00', '160-190': '#220505'});"
    # caption: Run it: this play icon, or F5
    And user hovers over "Run script (F5)" icon
    # caption: Its last line opens the table in a view
    And user hovers over the text "addTableView" in code editor
    And user presses F5
    Then the "demog 10000" view should be current
    And "height" column should be color-coded conditionally
    And "age" column should be color-coded linearly
    # caption: The script opened a demo table with its colors
    When user hovers over the "header height" area of grid
    Then no errors should have been logged
