Feature: A GraphPad Prism file opens from Browse
  Translated from TestTrack Curves/curves-data-to-curves.md, scenario 3. Curves 1.12.0 added a viewer
  and a handler for .pzfx files: a double-click in Browse, or Open in the context menu of the file,
  previews the file in a view of its name.
  Not claimed: the curve column of the preview grid draws charts — the renderer reports no readings.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Curves tree node inside browse tree is expanded

  Scenario: A double-click previews the file in a view of its name
    When user double-clicks on Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree
    Then the "dose-response-4pl.pzfx" view should be current
    And no errors should have been logged

  Scenario: Open from the context menu opens the file
    When user picks "Open" from the context menu of Files---App-Data---Curves---dose-response-4pl.pzfx tree node inside browse tree
    Then the "dose-response-4pl.pzfx" view should be current
    And no errors should have been logged
