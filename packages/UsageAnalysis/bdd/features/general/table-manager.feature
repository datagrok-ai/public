Feature: The Table Manager is docked and closed with Alt+T
  Alt+T docks the Table Manager, a grid of the open tables in a panel titled Tables, and closes it
  again. Translated from TestTrack General/table-manager.md, table-manager-spec.ts and the
  manual-only table-manager-ui.md (tried at the user's request).

  Not translated, and why: the md and the old spec open the manager from View | Tables, which the
  View menu no longer lists (GROK-21116; Alt+T is the command's shortcut).
  Everything about the manager's rows — the tables it lists, a click on a row making that table's view
  current and the table the current object, Open as table, several rows selected for the "N tables"
  submenu, Show > All adding and removing the attribute columns — needs the rows and cells of the
  manager's grid, and the widget on that element is the TableManager, which reports no readings and no
  areas (MISSING.md). The old spec's check that each name was "rendered" searched the whole page's
  text, and its toggle check read back the setter it had just called; neither is kept. Nothing is put
  on the server; the pane's shown state is the browser's (grok-settings), and the scenario closes it.

  Background:
    Given user is logged in
    And user opens demog dataset
    And user opens smiles dataset
    And user opens spgi dataset
    Then the open table views should be exactly "demog, smiles, spgi-100"

  Scenario: Alt+T docks the Table Manager and closes it again
    Then "Tables" dock panel should be absent
    When user presses Alt+T
    Then "Tables" dock panel should be visible
    And Grid viewer in "Tables" dock panel should be visible
    When user presses Alt+T
    Then "Tables" dock panel should be absent
    And no errors should have been logged
