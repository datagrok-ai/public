@viewers @realizes:viewers.tile-viewer
Feature: Tile viewer — a card designed on one table, then rebound to another
  The manual probe of `tile-viewer-ui.md` ("Edit Form with a different source table"): a card is
  designed on spgi-100, and the viewer's Table is then switched to demog-1000. Switching the table
  always rebuilds the card for the new table's columns (`tile-viewer.md`, "Multiple table
  switching"), so a person should see demog's fields with demog's values and nothing left of the
  spgi card — no field bound to a column demog does not have, no empty binding. The lanes column is
  set because `Edit Form...` is reached only through the viewer menu, which needs a region no card
  covers (see `tile-viewer-lanes.feature`).

  Fixtures: spgi-100 (Stereo Category R_ONE 36, S_ACHIR 34, S_UNKN 18, S_PART 10, S_ABS 2; the
  auto-generated card holds Idea Created, Primary Series Name, Scaffold Names), demog-1000 (the
  auto-generated card holds ten of its eleven columns, all but SEVERITY, in an order that changes
  from run to run — so the 10 is that set, and the list is claimed by what it holds).

  Scenario: The designed spgi card gives way to demog's own card when the table changes
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens spgi dataset
    And user switches to the "spgi-100" table view
    And user adds a tile viewer with:
      | Lanes Column Name | Stereo Category |
    Then tile viewer should be visible
    And the "lane names" reading of tile viewer should be "R_ONE, S_ABS, S_ACHIR, S_PART, S_UNKN"
    When user picks "Edit Form..." from the viewer menu of tile viewer
    Then form designer should be visible
    And the form designer should show value fields "Idea-Created, Primary-Series-Name, Scaffold-Names"
    When user deletes the "Scaffold-Names" value field in the form designer
    And user clicks on "CLOSE AND APPLY" button
    Then form designer should be absent
    And the "form designed" reading of tile viewer should be "true"
    And the "fields shown" reading of tile viewer should be 2
    And the "fields" reading of tile viewer should contain "Idea Created"
    And the "fields" reading of tile viewer should not contain "Scaffold Names"
    When user sets "Table" property of tile viewer to "demog-1000"
    Then tile viewer should be bound to table "demog-1000"
    And the "table" reading of tile viewer should be "demog-1000"
    And tile viewer should show 1000 rows
    And the "fields shown" reading of tile viewer should be 10
    And the "fields" reading of tile viewer should contain "SEX"
    And the "fields" reading of tile viewer should not contain "Primary Series Name"
    And the "fields" reading of tile viewer should not contain "Idea Created"
    And the "SEX of row 1" reading of tile viewer should be "F"
    And the "AGE of row 1" reading of tile viewer should be "26"
    And the "RACE of row 1" reading of tile viewer should be "Caucasian"
    And no errors should have been logged
