@journey @GROK-20507
Feature: The file drop overlay always goes away
  While files are dragged over the window the platform lays its "Incoming! Drop your CSV files to open
  them locally" layer over the whole workspace. A drag whose end the page never hears (the dragleave or
  the drop lost while the UI was busy) used to leave that layer over everything for good (GROK-20507):
  Escape pressed while the grid had the focus never reached it, and nothing else took it away. Escape
  now goes to it first, and a pointer move takes it away too, since no pointer events arrive while a
  real drag is on: a move means the drag is over. A file dropped after those dismissals still opens.
  Translated from playwright-public shell/file-drop-overlay.test.ts.

  The lost end is made in the page: a dragenter carrying a file, and nothing after it. A click on the
  layer, the old spec's third case, is a pointer move first wherever a person makes it, so the move
  claims it. The drop is the browser's own drag of a file from disk (the xlsx-open feature's step),
  which claims the layer showing and then going on the drop.

  Background:
    Given user is logged in
    And user opens cars dataset

  Scenario: Escape takes the layer away while the grid has the focus
    When user clicks on the "cell 2 of model" area of grid
    And a file dragged from disk enters grid and the drag's end never reaches the page
    Then file drop overlay should be visible
    When user presses Escape
    Then file drop overlay should be absent
    And no errors should have been logged

  Scenario: A pointer move takes the layer away
    When a file dragged from disk enters grid and the drag's end never reaches the page
    Then file drop overlay should be visible
    When user moves the pointer over file drop overlay
    Then file drop overlay should be absent
    And no errors should have been logged

  Scenario: A file dropped after the layer was dismissed still opens
    When user drops the "fixtures/cars-small.csv" file of the project onto grid
    Then the "cars-small" table view should open with 3 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
