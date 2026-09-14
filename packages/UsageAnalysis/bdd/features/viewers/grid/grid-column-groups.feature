@journey @viewers @realizes:viewers.grid
Feature: Grid column groups
  A column group is made from the Context Panel: select the columns by Control+clicking their
  headers, and Actions > Group columns... asks for the group's name, description and colour. The
  grid then draws a band above the member headers and reports it as the area `group <name>`, which
  spans the visible columns of the group and exists only while the group does; the band's name is
  drawn in the group's colour, and every member column carries the group's name in its `group`
  tag. A mouse-down on the band selects the group's columns. The guards are the ones
  `grid-dialogs-groups.md` scenario 4 names: a Shift+click across grouped headers (GROK-17505), a
  click on the group's band (which selects its columns), a second group, the first one's band
  selected and Escape (GROK-17442, GROK-18213), and the groups after a round-trip (GROK-17441):
  both groups with their names, their columns and their colours after a layout saved to the server
  and loaded over a view whose state was thrown away, and both groups with their names and their
  columns after a project, saved and opened with the library's API steps (operator decision D1).
  The colours are the band's own lettering, `#FF0000` and `#0000FF` against the `#4a4a49` an
  uncoloured band is drawn in.
  One journey on demog-1000, with two groups of neighbouring columns, Person (AGE, SEX) and Body
  (HEIGHT, WEIGHT), where the TestTrack spec groups AGE with HEIGHT and WEIGHT with SEX: for a
  group whose columns are not side by side the grid draws one band per run of members but reports
  one `group <name>` area from the first member to the last, over the columns between them too, so
  a click in the middle of that area lands on a column outside the group and selects every
  ungrouped column instead (measured on dev, reported to the operator).

  The first group goes through one layout round-trip before the second group is made. Expanding a
  pane of the Context Panel stamps `AppEvents.propertyEdited`, and for the next two seconds the
  platform drops every change of the current object for good (`setCurrentObject`, events.dart): a
  second selection made in that window leaves the panel on the first one, and its Group columns...
  groups the columns that are no longer selected (measured on dev: the whole first scenario takes
  1.3 s; reported to the operator). After the round-trip the panes are already expanded, nothing
  stamps that time again, and the panel follows the selection. The round-trip is a layout and not
  the project, because a project loses the colours (below) and both colours have to be on the
  grid when the layout of the second round-trip is saved.

  The colours are claimed on the layout round-trip, where they survive: a group's colour lives in
  the table tag `.columnGroups`, which the table view keeps in its layout's user data
  (`TableView.SYNC_TAGS`) and re-applies from it. The project scenarios claim the groups, their
  names and their columns and say nothing about colour: the library's `user saves the current view
  as project` puts only the layout's view state in the project and loses the user data the tag
  lives in, so the reopened bands are drawn in the default header colour - a request to the
  developers, filed apart from this feature. The `.columnGroups` table tag itself is not read (no
  step reads a table tag; the band, its colour and the column tag carry the claim). Of the
  `grid-ui.md` "Column Groups"
  checklist, collapsing and expanding a group, renaming a group in place and Ungroup from the header
  menu are not in the product: the core draws the band, selects on it and shows a tooltip
  (`grid_column_groups.dart`), the only command is Group columns... (`column_commands.dart`), and the
  header menu has no group item. Running Group columns... again on the same columns is how a group
  gets a new name and colour, and running it with an empty name is how it is ungrouped
  (`column_commands.dart` removes the `group` tag); the last two scenarios claim both. GROK-17443
  (a click on the band over columns that belong to no group) is not claimed: the grid reports an
  area only for a group, so the empty part of the band has no name to click. The Context Panel is
  checked by the column names it lists, not by the current object: a list of columns has no name
  for `the context panel should show`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And grid should not have a "group Person" area

  Scenario: Group columns... from the Context Panel draws a band in the group's colour
    When user clicks on the "header AGE" area of grid holding Control
    And user clicks on the "header SEX" area of grid holding Control
    Then columns "AGE, SEX" should be selected
    Given the context panel is open
    And Items accordion header in context panel is expanded
    Then context panel should contain text "AGE"
    And context panel should contain text "SEX"
    Given Actions accordion header in context panel is expanded
    When user clicks on "Group columns..." text in context panel
    Then Group input in dialog should be visible
    When user types "Person" into Group input in dialog
    And user types "#FF0000" into Color input in dialog
    And user clicks on OK button in dialog
    Then dialog should be hidden
    And grid should have a "group Person" area
    And "AGE" column should have tag "group" equal to "Person"
    And "SEX" column should have tag "group" equal to "Person"
    And the "group Person" area of grid should contain the color "#FF0000"
    And no errors should have been logged

  Scenario: Shift+click across the group, a click on its band and Escape raise no error
    When user clicks on the "header AGE" area of grid
    And user clicks on the "header SEX" area of grid holding Shift
    Then no errors should have been logged
    When user presses Escape
    Then no columns should be selected
    When user clicks on the "group Person" area of grid
    Then columns "AGE, SEX" should be selected
    And no errors should have been logged
    When user presses Escape
    Then no columns should be selected
    And no errors should have been logged

  Scenario: A layout saved to the server brings the group back over a fresh view of the table
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And grid should not have a "group Person" area
    And "AGE" column should have tag "group" equal to ""
    When user loads the saved layout
    Then grid should have a "group Person" area
    And "AGE" column should have tag "group" equal to "Person"
    And "SEX" column should have tag "group" equal to "Person"
    And no errors should have been logged

  Scenario: A second group gets its own band and colour; the first one's band selected and Escape raise no error
    When user clicks on the "header HEIGHT" area of grid holding Control
    And user clicks on the "header WEIGHT" area of grid holding Control
    Then columns "HEIGHT, WEIGHT" should be selected
    Given the context panel is open
    Then context panel should contain text "HEIGHT"
    And context panel should contain text "WEIGHT"
    Given Actions accordion header in context panel is expanded
    When user clicks on "Group columns..." text in context panel
    And user types "Body" into Group input in dialog
    And user types "#0000FF" into Color input in dialog
    And user clicks on OK button in dialog
    Then dialog should be hidden
    And grid should have a "group Body" area
    And grid should have a "group Person" area
    And "HEIGHT" column should have tag "group" equal to "Body"
    And "AGE" column should have tag "group" equal to "Person"
    And the "group Body" area of grid should contain the color "#0000FF"
    And the "group Body" area of grid should not contain the color "#FF0000"
    When user clicks on the "group Person" area of grid
    Then columns "AGE, SEX" should be selected
    When user presses Escape
    Then no columns should be selected
    And no errors should have been logged

  Scenario: A layout saved to the server restores the groups, their columns and their colours
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    Then grid should show 1000 rows
    And grid should not have a "group Person" area
    And grid should not have a "group Body" area
    And "AGE" column should have tag "group" equal to ""
    When user loads the saved layout
    Then grid should have a "group Person" area
    And grid should have a "group Body" area
    And "AGE" column should have tag "group" equal to "Person"
    And "SEX" column should have tag "group" equal to "Person"
    And "HEIGHT" column should have tag "group" equal to "Body"
    And "WEIGHT" column should have tag "group" equal to "Body"
    And the "group Person" area of grid should contain the color "#FF0000"
    And the "group Body" area of grid should contain the color "#0000FF"
    And the "group Person" area of grid should not contain the color "#0000FF"
    And no errors should have been logged

  Scenario: Both groups and their columns come back from a project
    When user saves the current view as project "zz-grid-column-groups-2"
    And user closes all views
    And user opens the "zz-grid-column-groups-2" project
    Then grid should show 1000 rows
    And grid should have a "group Person" area
    And grid should have a "group Body" area
    And "AGE" column should have tag "group" equal to "Person"
    And "WEIGHT" column should have tag "group" equal to "Body"
    And no errors should have been logged

  Scenario: Grouping the same columns again gives the band a new name and colour
    When user clicks on the "group Person" area of grid
    Then columns "AGE, SEX" should be selected
    Given the context panel is open
    Then context panel should contain text "AGE"
    And context panel should contain text "SEX"
    Given Actions accordion header in context panel is expanded
    When user clicks on "Group columns..." text in context panel
    And user types "Measures" into Group input in dialog
    And user types "#d62728" into Color input in dialog
    And user clicks on OK button in dialog
    Then dialog should be hidden
    And grid should have a "group Measures" area
    And grid should not have a "group Person" area
    And the "group Measures" area of grid should contain the color "#d62728"
    And "SEX" column should have tag "group" equal to "Measures"
    And "AGE" column should have tag "group" equal to "Measures"
    When user presses Escape
    Then no errors should have been logged

  Scenario: Grouping with an empty name ungroups the columns
    When user clicks on the "group Measures" area of grid
    Then columns "AGE, SEX" should be selected
    Given the context panel is open
    Then context panel should contain text "AGE"
    Given Actions accordion header in context panel is expanded
    When user clicks on "Group columns..." text in context panel
    And user clears Group input in dialog
    And user clicks on OK button in dialog
    Then dialog should be hidden
    And grid should not have a "group Measures" area
    And grid should have a "group Body" area
    When user presses Escape
    Then no errors should have been logged
