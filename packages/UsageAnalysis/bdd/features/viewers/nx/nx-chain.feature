@journey @viewers
Feature: The NX chain: linked tables, calculated columns, viewers, legends and filters over saved projects
  The five NX scenarios of TestTrack (linking, calc-columns, viewers-for-linked-tables,
  formula-lines-and-legend, filtering) as one journey in their order: each part opens the project
  the one before saved — NxProject, NxProjectCalcColumns, NxProjectViewers, NxProjectFormulaLegend,
  NxProjectFiltering, each named with the run's suffix so parallel runs do not meet — and the last
  scenario deletes all five. The Background lists the server's projects once and deletes whatever
  an earlier run of the same names left; at feature end each of the five is looked up by name and
  deleted again, and the last scenario proves them gone from a listing of every project.
  SPGI (3624 rows), SPGI-linked1 (3624) and SPGI-linked2 (224) are opened from the Files tree, so
  the tables carry their creation script and every save is made through the ribbon's Save dialog
  with Data sync on; the dialog's "Save a copy" makes each next project, and the first save closes
  the Share dialog it opens. Link Tables links SPGI to SPGI-linked1 by Id / Concept Id (selection
  to selection) and SPGI-linked1 to SPGI-linked2 by Sample Name and link column 1-3 (selection to
  filter); what the links carry is checked row by row against the keys of the rows selected in the
  table the link starts from, not only as counts (5 rows selected in SPGI select 9 in SPGI-linked1
  and leave 12 SPGI-linked2 rows; all rows 191 and 183; no selection all 224).
  Readings, not pictures: a viewer's "rows shown" against the rows of its table that pass the
  filter (for a scatter plot, those it can place on its axes), its axis range against its filtered
  rows (a check that fails when the filtered rows span most of both axes, since they cannot tell a
  zoom from none), the properties set read back, the columns the grid pins, the Formula Lines list,
  the "formula lines" reading and — for a line chart — the "formula line <title>" areas it reports
  for a line it drew, legend placement, item counts and colors, the Scaffold Tree card's own
  readings (read through the filter the panel holds: the viewer inside a filter card has no viewer
  name of its own), and what each view's filter panel filters by. The chain waits for the Chem
  filters (no substructure card searching, every checked scaffold node counted) before it counts.
  The stand needs the ApiTests datasets (System:AppData/ApiTests/datasets) and a Chem package
  whose Scaffold Tree and substructure card report their readings (published from this branch);
  dev runs an older Chem, so the filtering part runs on localhost only. The journey takes about
  five minutes on localhost.

  Done through the API, and why: the projects are reopened by their friendly name through the
  project API (the Dashboards gallery search finds no name that holds "-", and a run-suffixed name
  does); Close All is the shell's closeAll; a layout is saved and applied through the view and kept
  by name in the page (a layout saved to the gallery would be offered to every parallel run's SPGI,
  and the md does not name the gesture); a grid is scrolled to a column (scrollToCell) before its
  header menu is opened; viewer settings the md lists without a place are set as properties and
  read back (logarithmic axes, Zoom and Filter, tables, color columns, legend visibility and
  position, Row Source, label columns, the line charts' X and Y columns — a line chart on dates
  covers its X selector with the time-unit one); the conditional bins and the linear scheme of
  Chem Space X are written to the column after Color Coding > Conditional and > Linear are picked
  from its header menu (the color-coding editor is not driven).
  How the md's words are read: "not freezing, cannot be broken, check the legend" is every viewer
  of the view reporting a status with no error once its render has settled, plus the legend
  readings. The first bar chart is stacked, so a click on the Pyrrolidines bar filters to one of its
  segments: fewer than its 874 rows, none of another series. The rows selected for the scatter
  plots are the R_ONE | Chemist 2 segment of the second bar chart: their SPGI-linked2 rows cover a
  narrow part of both axes, which is what lets pack and zoom be seen. After Escape the md expects the
  third scatter plot to be empty: with the link as Link Tables makes it ("Filter All On No Rows
  Selected" off, link.dart) the SPGI-linked2 plot shows all 224 rows and the SPGI-linked1 plot (Row
  Source Selected) is the empty one; with that option switched on in the link through Link Tables,
  Escape leaves SPGI-linked2 no row, as the md expects, and the option is switched off again.
  "Stack them one over another" is the histogram, bar chart, pie chart, PC plot, box plot and
  trellis plot dragged by their title bars onto each other's dock compass as the tabs of one panel;
  each legend is read on its tab after switching to it, and again after the saved layout is applied
  over an earlier one and a move. "Other tabs" in the filtering part is the grid of another view
  showing the rows the shared table lets through. The first view's panel is emptied before the
  Scaffold Tree filter is added, so its card is on screen; "Use as filter" is taken on the drawn
  molecule with the longest SMILES; Filter type Categorical is checked on a Core card added to the
  clone's own panel; turning the master switch off lets more rows through, not all, since the other
  views' panels keep filtering the shared table.
  Not translated: resizing the grid when the layout is changed (closing the filter panel is the
  change); a new linear color scheme set on the first scatter plot itself (the scheme is changed on
  the column later, as the md's next steps do); the regression line drawn — with Row Source
  FilteredSelected and nothing selected the plot draws no rows, so the property is claimed; a
  scatter plot's formula lines drawn — the scatter plot reports no area for a formula line (the
  line chart does), so for scatter plots the "formula lines" reading (items shown, not drawn) and
  the look are claimed; the trellis plot's legend and the box plot's (neither draws one here); the
  Context Panel as the place a Legend Position is read.
  Known failures: a layout saved before Chemical Space X was renamed, applied after the rename,
  puts the old column name back into the "${Spec} result" formula, and a project saved afterwards
  warns on reopen that the column is missing (localhost and dev); the calc-columns layout is
  therefore applied before the layout this part saved, which puts the renamed formula back. A
  molecule used as a filter (Current Value > Use as filter) in a view while the first view's
  Scaffold Tree filters the same Structure column puts a card that reports the molecule and
  "filtering", yet the rows that pass do not contain the molecule (localhost; dev runs an older
  Chem). Not claimed, because of it: how many rows pass after the master switch goes off and on
  and after NxProjectFiltering is reopened — every view's panel reports the same cards and
  summaries, but the count differs from run to run (9 before; 9 or 1 after).

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the molecule sketcher is "OpenChemLib"
    And the projects "NxProject-{run}, NxProjectCalcColumns-{run}, NxProjectViewers-{run}, NxProjectFormulaLegend-{run}, NxProjectFiltering-{run}" are deleted now and when the feature ends
    And the browse panel is open

  Scenario: Linking - SPGI, SPGI-linked1 and SPGI-linked2 open from the Files tree with data sync
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---chem tree node inside browse tree is expanded
    When user double-clicks on Files---Demo---chem---SPGI.csv tree node inside browse tree
    Then the "SPGI" table view should open with 3624 rows
    When user clicks on browse tab
    Given Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---ApiTests tree node inside browse tree is expanded
    And Files---App-Data---ApiTests---datasets tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---ApiTests---datasets---SPGI-linked1.csv tree node inside browse tree
    Then the "SPGI-linked1" table view should open with 3624 rows
    When user clicks on browse tab
    And user double-clicks on Files---App-Data---ApiTests---datasets---SPGI-linked2.csv tree node inside browse tree
    Then the "SPGI-linked2" table view should open with 224 rows
    And no errors should have been logged

  Scenario: Linking - Link Tables links SPGI to SPGI-linked1 and SPGI-linked1 to SPGI-linked2
    When user clicks on the tab of the "SPGI" view
    And user picks "Data > Link Tables..." from the top menu
    And user sets the tables of the Link Tables dialog to "SPGI" and "SPGI-linked1"
    And user sets key columns 1 of the Link Tables dialog to "Id" and "Concept Id"
    And user selects "selection to selection" in Link Type input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "SPGI -> SPGI-linked1" text in "Link Tables" dialog should be visible
    When user clicks on "New Link" text in "Link Tables" dialog
    And user sets the tables of the Link Tables dialog to "SPGI-linked1" and "SPGI-linked2"
    And user sets key columns 1 of the Link Tables dialog to "Sample Name" and "Sample Name"
    And user clicks on Add button in "Link Tables" dialog
    And user sets key columns 2 of the Link Tables dialog to "link column 1" and "link column 1"
    And user clicks on Add button in "Link Tables" dialog
    And user sets key columns 3 of the Link Tables dialog to "link column 2" and "link column 2"
    And user clicks on Add button in "Link Tables" dialog
    And user sets key columns 4 of the Link Tables dialog to "link column 3" and "link column 3"
    And user selects "selection to filter" in Link Type input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "SPGI-linked1 -> SPGI-linked2" text in "Link Tables" dialog should be visible
    When user clicks on CLOSE button in "Link Tables" dialog
    Then 224 rows of table "SPGI-linked2" should pass the filter
    And no errors should have been logged

  Scenario: Linking - a line chart on SPGI-linked2 follows the selection made in SPGI
    When user clicks on the tab of the "SPGI" view
    And user clicks on line chart icon on toolbox
    Then the open tableview should have 1 line chart viewer
    When user clicks on settings icon of line chart viewer
    Given "Data" category in context panel is expanded
    When user selects "SPGI-linked2" in "Table" property in context panel
    Then line chart viewer should be bound to table "SPGI-linked2"
    When user enters "${link column 3}==\"v ii\" && ${link column 1} <30" in "Filter" property in context panel
    And user adds "link column 2" to the splits of line chart viewer
    And user selects "link column 1" in "Overview" property in context panel
    And user picks "Value1" in the "x" column selector of line chart viewer
    Given "Y Axis" category in context panel is expanded
    When user selects "logarithmic" in "Y Axis Type" property in context panel
    Then properties of line chart viewer should be:
      | Table        | SPGI-linked2                                          |
      | Filter       | ${link column 3}=="v ii" && ${link column 1} <30      |
      | Split Column Names | link column 2                                   |
      | Overview     | link column 1                                         |
      | X            | Value1                                                |
      | Y Axis Type  | logarithmic                                           |
    And line chart viewer should show 97 rows
    And line chart viewer should show the rows of table "SPGI-linked2" that pass the filter where "link column 3" is "v ii" and "link column 1" is below 30
    And line chart viewer should report no error
    And no errors should have been logged

  Scenario: Linking - selecting rows in SPGI narrows the line chart and clearing the selection restores it
    When user drags the "row header 1" area of grid to the "row header 5" area
    Then 5 rows should be selected
    And 9 rows of table "SPGI-linked1" should be selected
    And the rows selected in table "SPGI-linked1" should be exactly those matching the rows selected in table "SPGI" on "Concept Id" = "Id"
    And 12 rows of table "SPGI-linked2" should pass the filter
    And the rows of table "SPGI-linked2" that pass the filter should be exactly those matching the rows selected in table "SPGI-linked1" on "Sample Name, link column 1, link column 2, link column 3" = "Sample Name, link column 1, link column 2, link column 3"
    And line chart viewer should show 3 rows
    And line chart viewer should show the rows of table "SPGI-linked2" that pass the filter where "link column 3" is "v ii" and "link column 1" is below 30
    When user clicks on the "cell 1 of Id" area of grid
    And user presses Control+A
    Then all rows should be selected
    And 191 rows of table "SPGI-linked1" should be selected
    And 183 rows of table "SPGI-linked2" should pass the filter
    And line chart viewer should show 76 rows
    And line chart viewer should show the rows of table "SPGI-linked2" that pass the filter where "link column 3" is "v ii" and "link column 1" is below 30
    When user presses Escape
    Then no rows should be selected
    And 0 rows of table "SPGI-linked1" should be selected
    And 224 rows of table "SPGI-linked2" should pass the filter
    And line chart viewer should show 97 rows
    And line chart viewer should report no error
    And no errors should have been logged

  Scenario: Linking - the tables, their links and the line chart are saved with data sync as NxProject
    When user opens the Save project dialog from the ribbon
    Then the Save project dialog should save the tables "SPGI, SPGI-linked1, SPGI-linked2" with data sync
    When user types "NxProject-{run}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    And "Share NxProject-{run}" dialog should be visible
    When user presses Escape
    Then "Share NxProject-{run}" dialog should be absent
    And 1 project named "NxProject-{run}" should be on the server
    And the project saved as "NxProject-{run}" should link "SPGI -> SPGI-linked1 by Id = Concept Id as selection to selection; SPGI-linked1 -> SPGI-linked2 by Sample Name, link column 1, link column 2, link column 3 = Sample Name, link column 1, link column 2, link column 3 as selection to filter"
    And no errors should have been logged
    When user closes all views

  Scenario: Calc columns - the NxProject project opens with its three tables, the link and the line chart
    When user opens the project saved as "NxProject-{run}"
    Then table "SPGI" should be open
    And no error or warning balloon should have been shown
    And table "SPGI-linked1" should be open
    And table "SPGI-linked2" should be open
    When user clicks on the tab of the "SPGI" view
    Then line chart viewer should be bound to table "SPGI-linked2"
    And "Filter" property of line chart viewer should be "${link column 3}==\"v ii\" && ${link column 1} <30"
    And line chart viewer should show 97 rows
    When user drags the "row header 1" area of grid to the "row header 5" area
    Then 9 rows of table "SPGI-linked1" should be selected
    And 12 rows of table "SPGI-linked2" should pass the filter
    And line chart viewer should show 3 rows
    When user presses Escape
    Then no rows should be selected
    And 224 rows of table "SPGI-linked2" should pass the filter
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Calc columns - Add New Column warns about if(num, qnum) and adds the column the second formula makes
    When user clicks on the tab of the "SPGI" view
    And user picks "View > Layout > Clone View" from the top menu
    Then the "SPGI copy" view should be current
    When user picks "Edit > Add New Column..." from the top menu
    Then "Add New Column" dialog should be visible
    When user pastes "if(${Whole blood assay 1} != null, ${Whole blood assay 1}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}==\"Dog\", 30.9, if(${Species}==\"Monkey\", 43.6, if(${Species}==\"Minipig\", 39, null)))))*100" into code editor in "Add New Column" dialog
    Then "Add New Column" dialog should contain text "If function params types (qnum, num) do not match and cannot be casted to each other"
    When user pastes "if(${NIBR logP} != null, ${NIBR logP}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}==\"Dog\", 30.9, if(${Species}==\"Monkey\", 43.6, if(${Species}==\"Minipig\", 39, null)))))*100" into code editor in "Add New Column" dialog
    Then "Add New Column" dialog should not contain text "If function params types"
    When user types "${Species} result" into new column name input
    And user clicks on OK button in "Add New Column" dialog
    Then the "Add New Column" dialog should close
    And the table should have a column "${Species} result"
    And "${Species} result" column should have type "double"
    And "${Species} result" column should have tag "formula" equal to "if(${NIBR logP} != null, ${NIBR logP}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}==\"Dog\", 30.9, if(${Species}==\"Monkey\", 43.6, if(${Species}==\"Minipig\", 39, null)))))*100"
    And no errors should have been logged

  Scenario: Calc columns - Order or Hide Columns leaves six columns and the column menu adds the result filter
    When user picks "Order or Hide Columns..." from the context menu of the "header ${Species} result" area of grid
    Then Order or Hide Columns dialog should be visible
    When user unchecks all columns checkbox
    And user types "NIBR logP" into "Search" input in Order or Hide Columns dialog
    And user toggles the "NIBR logP" column in the column list of Order or Hide Columns dialog
    And user types "Whole blood assay 1" into "Search" input in Order or Hide Columns dialog
    And user toggles the "Whole blood assay 1" column in the column list of Order or Hide Columns dialog
    And user types "Chemical Space X" into "Search" input in Order or Hide Columns dialog
    And user toggles the "Chemical Space X" column in the column list of Order or Hide Columns dialog
    And user types "Average Mass" into "Search" input in Order or Hide Columns dialog
    And user toggles the "Average Mass" column in the column list of Order or Hide Columns dialog
    And user types "Species" into "Search" input in Order or Hide Columns dialog
    And user toggles the "Species" column in the column list of Order or Hide Columns dialog
    And user types "${Species} result" into "Search" input in Order or Hide Columns dialog
    And user toggles the "${Species} result" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then the "column order" reading of grid should be "Chemical Space X, Average Mass, NIBR logP, Whole blood assay 1, Species, ${Species} result"
    When user hovers over the "header ${Species} result" area of grid
    And user clicks on "Column options" icon in grid
    Then column popup should be visible
    When user clicks on "Add filter" action in column popup
    Then "${Species} result" filter card should be visible
    And no errors should have been logged

  Scenario: Calc columns - filtering out missing results keeps the 30 rows that have one
    When user picks "Missing values > Filter out missing values" from the indicator menu of the "${Species} result" filter card
    Then 30 rows should pass the filter
    And "${Species} result" column should have missing values
    And grid should show 30 rows
    And no errors should have been logged

  Scenario: Calc columns - renaming Species renames the result column, and an edited value is recalculated
    When user picks "Column Properties..." from the context menu of the "header Species" area of grid
    Then "Species" dialog should be visible
    When user types "Spec" into "New name:" input in "Species" dialog
    And user clicks on OK button in "Species" dialog
    Then the table should have a column "Spec"
    And the table should have a column "${Spec} result"
    And the table should not have a column "${Species} result"
    And "${Spec} result" column should have tag "formula" equal to "if(${NIBR logP} != null, ${NIBR logP}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Spec}, 'Rat') || Contains(${Spec}, 'Rat Legacy'), 80, if(Contains(${Spec}, 'Mouse'), 125, if(${Spec}==\"Dog\", 30.9, if(${Spec}==\"Monkey\", 43.6, if(${Spec}==\"Minipig\", 39, null)))))*100"
    And the "column order" reading of grid should be "Chemical Space X, Average Mass, NIBR logP, Whole blood assay 1, Spec, ${Spec} result"
    And the value of "${Spec} result" column in row 962 should be "4.353702545166016"
    When user double-clicks on the "cell 962 of NIBR logP" area of grid
    Then cell editor should be visible
    When user presses Control+A in cell editor
    And user types "8" into cell editor
    And user presses Enter in cell editor
    Then the value of "NIBR logP" column in row 962 should be "8"
    And the value of "${Spec} result" column in row 962 should be "10"
    And no errors should have been logged

  Scenario: Calc columns - a clone of the view keeps the six columns and the filter, and gets a pie chart summary column
    When user picks "View > Layout > Clone View" from the top menu
    Then the "SPGI copy copy" view should be current
    And the "column order" reading of grid should be "Chemical Space X, Average Mass, NIBR logP, Whole blood assay 1, Spec, ${Spec} result"
    And grid should show 30 rows
    And "${Spec} result" filter card should be visible
    When user picks "Add > Summary Columns > Pie Chart" from the context menu of the "cell 962 of Spec" area of grid
    Then the "column order" reading of grid should include the text "Pie"
    And no errors should have been logged

  Scenario: Calc columns - a saved layout brings the filter panel back after it was closed
    When user saves the layout of the current table view as "calc-columns"
    And user clicks on filter icon in toolbar
    Then filter panel should be absent
    And 30 rows should pass the filter
    When user applies the layout "calc-columns" to the current table view
    Then "${Spec} result" filter card should be visible
    And 30 rows should pass the filter
    And the "column order" reading of grid should include the text "Spec, ${Spec} result"
    And the "column order" reading of grid should include the text "Pie"
    And no errors should have been logged

  Scenario: Calc columns - Table > Add view opens every SPGI column without a filter panel, filtered like the other views
    When user picks "Table > Add View" from the context menu of the current view tab
    Then filter panel should be absent
    And grid should show 30 rows
    And the "column order" reading of grid should include the text "Id, Structure, CAST Idea ID"
    And the "column order" reading of grid should include the text "TK6 MNT, Spec, Strain"
    And the "column order" reading of grid should include the text "Compound Status, ${Spec} result"
    And no errors should have been logged

  Scenario: Calc columns - the new view's panel takes Structure and categorical cards, Structure moves to the top, and Reset filter clears every view
    When user clicks on filter icon in toolbar
    And user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    When user picks "Select Columns..." from the filter panel menu
    And user types "Series" into "Search" input in "Select columns..." dialog
    And user toggles the "Series" column in the column list of "Select columns..." dialog
    And user types "Stereo Category" into "Search" input in "Select columns..." dialog
    And user toggles the "Stereo Category" column in the column list of "Select columns..." dialog
    And user types "Structure" into "Search" input in "Select columns..." dialog
    And user toggles the "Structure" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "cards" reading of filter panel should be "Series, Stereo Category, Structure"
    When user drags the "Structure" filter card above the "Series" filter card
    Then the "cards" reading of filter panel should be "Structure, Series, Stereo Category"
    And 30 rows should pass the filter
    When user picks "View > Reset Filter" from the top menu
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: Calc columns - a copy of the project is saved as NxProjectCalcColumns
    When user opens the Save project dialog from the ribbon
    When user clicks on "Save a copy" text in "Save project" dialog
    And user types "NxProjectCalcColumns-{run}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    And 1 project named "NxProjectCalcColumns-{run}" should be on the server
    And 1 project named "NxProject-{run}" should be on the server
    When user closes all views

  Scenario: Viewers - the NxProjectCalcColumns project opens with its views and a clone of the last one pins two columns and colors Chemist 521 by category
    When user opens the project saved as "NxProjectCalcColumns-{run}"
    Then table "SPGI" should be open
    And no error or warning balloon should have been shown
    And table "SPGI-linked2" should be open
    And no errors should have been logged
    When user switches to the last table view of "SPGI"
    And user picks "View > Layout > Clone View" from the top menu
    And user picks "Pin > Pin Column" from the context menu of the "header Structure" area of grid
    And user picks "Pin > Pin Column" from the context menu of the "header Id" area of grid
    Then "Frozen Columns" property of grid should be "3"
    And the grid should pin the columns "Structure, Id"
    When user scrolls the grid to the "Chemist 521" column
    And user picks "Color Coding > Categorical" from the context menu of the "header Chemist 521" area of grid
    Then "Chemist 521" column should be color-coded categorically
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Viewers - three scatter plots on the three tables, two bar charts, and Chemist 27 recolored from the legend
    When user clicks on scatter plot icon on toolbox
    And user clicks on scatter plot icon on toolbox
    And user clicks on scatter plot icon on toolbox
    Then the open tableview should have 3 scatter plot viewers
    When user sets properties of first scatter plot viewer:
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | Series                  |
    And user sets properties of second scatter plot viewer:
      | Table           | SPGI-linked1            |
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | Primary Series Name     |
      | Row Source      | Selected                |
    And user sets properties of third scatter plot viewer:
      | Table           | SPGI-linked2            |
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | link column 2           |
    Then properties of first scatter plot viewer should be:
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | Series                  |
    And properties of second scatter plot viewer should be:
      | Table           | SPGI-linked1            |
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | Primary Series Name     |
      | Row Source      | Selected                |
    And properties of third scatter plot viewer should be:
      | Table           | SPGI-linked2            |
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | link column 2           |
    And first scatter plot viewer should be bound to table "SPGI"
    And second scatter plot viewer should be bound to table "SPGI-linked1"
    And third scatter plot viewer should be bound to table "SPGI-linked2"
    And the legend of first scatter plot viewer should be docked
    And the legend of third scatter plot viewer should be docked
    When user clicks on bar chart icon on toolbox
    And user sets properties of first bar chart viewer:
      | Split Column | Primary Series Name |
      | Stack Column | Scaffold Names      |
      | On Click     | Filter              |
    Then properties of first bar chart viewer should be:
      | Split Column | Primary Series Name |
      | Stack Column | Scaffold Names      |
      | On Click     | Filter              |
    When user clicks on bar chart icon on toolbox
    And user sets properties of second bar chart viewer:
      | Value Column    | Average Mass |
      | Value Aggr Type | sum          |
      | Stack Column    | Chemist 521  |
    Then the open tableview should have 2 bar chart viewers
    And properties of second bar chart viewer should be:
      | Value Column    | Average Mass |
      | Value Aggr Type | sum          |
      | Stack Column    | Chemist 521  |
    When user hovers over "Chemist 27" legend item in legend of second bar chart viewer
    And user clicks on color picker icon
    Then "Chemist 27" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    And user clicks on OK button in "Chemist 27" dialog
    Then "Chemist 27" dialog should be absent
    And the categorical color of "Chemist 27" in "Chemist 521" column should be "#9467BD"
    And the "Chemist 27" item in the legend of second bar chart viewer should be colored "#9467BD"
    And "Chemist 521" column should be color-coded categorically
    And no viewer of the current view should report an error
    And no errors should have been logged

  Scenario: Viewers - a bar filters SPGI and another selects in it, and the scatter plots follow through the links
    When user clicks on the "bar Pyrrolidines" area of first bar chart viewer
    Then fewer than 874 rows should pass the filter
    And no rows where "Primary Series Name" is "Aminopiperidines" should pass the filter
    And no rows where "Primary Series Name" is "Triazoles" should pass the filter
    And first scatter plot viewer should show every row that passes the filter of its table
    And the "rows shown" reading of first scatter plot viewer should be at least 100
    And first scatter plot viewer should be zoomed to the rows that pass the filter
    When user clicks on the "bar R_ONE | Chemist 2" area of second bar chart viewer
    Then some rows should be selected
    And every selected row should pass the filter
    And no rows where "Stereo Category" is "S_PART" should be selected
    And no rows where "Chemist 521" is "Chemist 27" should be selected
    And the rows selected in table "SPGI-linked1" should be exactly those matching the rows selected in table "SPGI" on "Concept Id" = "Id"
    And the rows of table "SPGI-linked2" that pass the filter should be exactly those matching the rows selected in table "SPGI-linked1" on "Sample Name, link column 1, link column 2, link column 3" = "Sample Name, link column 1, link column 2, link column 3"
    And second scatter plot viewer should show the selected rows of its table
    And third scatter plot viewer should show every row that passes the filter of its table
    And third scatter plot viewer should be zoomed to the rows that pass the filter
    And no viewer of the current view should report an error
    And no errors should have been logged

  Scenario: Viewers - a clone of the view keeps the viewers, the pinned columns and the filtering
    When user picks "View > Layout > Clone View" from the top menu
    Then the "SPGI (2) copy copy" view should be current
    And the open tableview should have 3 scatter plot viewers
    And the open tableview should have 2 bar chart viewers
    And "Frozen Columns" property of grid should be "3"
    And the grid should pin the columns "Structure, Id"
    And first scatter plot viewer should be bound to table "SPGI"
    And second scatter plot viewer should be bound to table "SPGI-linked1"
    And third scatter plot viewer should be bound to table "SPGI-linked2"
    And properties of first scatter plot viewer should be:
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | Series                  |
    And properties of second scatter plot viewer should be:
      | Table           | SPGI-linked1            |
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | Primary Series Name     |
      | Row Source      | Selected                |
    And properties of third scatter plot viewer should be:
      | Table           | SPGI-linked2            |
      | Zoom And Filter | pack and zoom by filter |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Color           | link column 2           |
    And fewer than 874 rows should pass the filter
    And no rows where "Primary Series Name" is "Triazoles" should pass the filter
    And first scatter plot viewer should show every row that passes the filter of its table
    And the rows selected in table "SPGI-linked1" should be exactly those matching the rows selected in table "SPGI" on "Concept Id" = "Id"
    And the rows of table "SPGI-linked2" that pass the filter should be exactly those matching the rows selected in table "SPGI-linked1" on "Sample Name, link column 1, link column 2, link column 3" = "Sample Name, link column 1, link column 2, link column 3"
    And second scatter plot viewer should show the selected rows of its table
    And third scatter plot viewer should show every row that passes the filter of its table
    And no viewer of the current view should report an error
    And no errors should have been logged

  Scenario: Viewers - a double-click on the bar chart resets the filtering and Escape the selection, on every view
    When user clicks on the tab of the "SPGI (2) copy" view
    And user double-clicks on empty plot space of first bar chart viewer
    Then all rows should pass the filter
    And first scatter plot viewer should show every row that passes the filter of its table
    When user presses Escape
    Then no rows should be selected
    And 0 rows of table "SPGI-linked1" should be selected
    And the "rows shown" reading of second scatter plot viewer should be 0
    And 224 rows of table "SPGI-linked2" should pass the filter
    And third scatter plot viewer should show every row that passes the filter of its table
    When user clicks on the tab of the "SPGI (2) copy copy" view
    Then all rows should pass the filter
    And no rows should be selected
    And the "rows shown" reading of second scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: Viewers - with Filter All On No Rows Selected switched on in the link, Escape empties the SPGI-linked2 plot
    When user picks "Data > Link Tables..." from the top menu
    And user clicks on "SPGI-linked1 -> SPGI-linked2" text in "Link Tables" dialog
    And user checks "Filter All On No Rows Selected" input in "Link Tables" dialog
    Then "Filter All On No Rows Selected" input in "Link Tables" dialog should be checked
    When user clicks on CLOSE button in "Link Tables" dialog
    And user clicks on the "bar R_ONE | Chemist 2" area of second bar chart viewer
    Then some rows should be selected
    And the rows of table "SPGI-linked2" that pass the filter should be exactly those matching the rows selected in table "SPGI-linked1" on "Sample Name, link column 1, link column 2, link column 3" = "Sample Name, link column 1, link column 2, link column 3"
    When user presses Escape
    Then no rows should be selected
    And 0 rows of table "SPGI-linked2" should pass the filter
    And the "rows shown" reading of third scatter plot viewer should be 0
    When user picks "Data > Link Tables..." from the top menu
    And user clicks on "SPGI-linked1 -> SPGI-linked2" text in "Link Tables" dialog
    And user unchecks "Filter All On No Rows Selected" input in "Link Tables" dialog
    And user clicks on CLOSE button in "Link Tables" dialog
    And user clicks on the "bar R_ONE | Chemist 2" area of second bar chart viewer
    And user presses Escape
    Then 224 rows of table "SPGI-linked2" should pass the filter
    And no errors should have been logged

  Scenario: Viewers - a copy of the project is saved as NxProjectViewers
    When user opens the Save project dialog from the ribbon
    When user clicks on "Save a copy" text in "Save project" dialog
    And user types "NxProjectViewers-{run}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    And 1 project named "NxProjectViewers-{run}" should be on the server
    And 1 project named "NxProjectCalcColumns-{run}" should be on the server
    When user closes all views

  Scenario: Formula lines - the NxProjectViewers project opens, and a new view gets a qualified-number column
    When user opens the project saved as "NxProjectViewers-{run}"
    Then table "SPGI" should be open
    And no error or warning balloon should have been shown
    And the open tableview should have 3 scatter plot viewers
    And no viewer of the current view should report an error
    And no errors should have been logged
    When user picks "Table > Add View" from the context menu of the current view tab
    And user picks "Edit > Add New Column..." from the top menu
    And user pastes "Qnum((100-${Chemical Space Y})/100, if(qualifier(${Chemical Space X})==\">\", \"<\", \"=\"))" into code editor in "Add New Column" dialog
    And user types "${Chemical Space Y} ${Chemical Space X}" into new column name input
    And user clicks on OK button in "Add New Column" dialog
    Then the "Add New Column" dialog should close
    And the table should have a column "${Chemical Space Y} ${Chemical Space X}"
    And "${Chemical Space Y} ${Chemical Space X}" column should have type "qnum"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Formula lines - a scatter plot gets three formula lines, a horizontal and a vertical one, labels and a color column
    When user clicks on scatter plot icon on toolbox
    Then the open tableview should have 1 scatter plot viewer
    When user picks "Chemical Space X" in the "x" column selector of scatter plot viewer
    And user sets properties of scatter plot viewer:
      | Y Column        | Chemical Space Y        |
      | X Axis Type     | logarithmic             |
      | Y Axis Type     | logarithmic             |
      | Zoom And Filter | pack and zoom by filter |
    And user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user adds the formula line "${Chemical Space Y} = 0.5* ${Chemical Space X}* ${Chemical Space X} - 1.5 * ${Chemical Space X} -1" in the Formula Lines dialog
    And user adds the formula line "${Chemical Space Y} = 0.1* ${Chemical Space X}* ${Chemical Space X} +  ${Chemical Space X}" in the Formula Lines dialog
    And user adds the formula line "${Chemical Space Y} = 0.1* ${Chemical Space X}" in the Formula Lines dialog
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Horizontal" from the open menu
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Line - Vertical" from the open menu
    And user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And the "formula lines" reading of scatter plot viewer should be 5
    And "formulaLines" property of scatter plot viewer should contain "${Chemical Space Y} = 0.1* ${Chemical Space X}"
    When user sets properties of scatter plot viewer:
      | Label Columns   | Structure, Id |
      | showLabelsFor   | Selected      |
    And user picks "Chemical Space X" in the "color" column selector of scatter plot viewer
    Then "colorColumnName" property of scatter plot viewer should be "Chemical Space X"
    And no viewer of the current view should report an error
    And no errors should have been logged

  Scenario: Formula lines - Pick Up / Apply copies the look to a second scatter plot, whose lines are then edited, doubled, hidden and shown again
    When user clicks on scatter plot icon on toolbox
    Then the open tableview should have 2 scatter plot viewers
    When user picks "Pick Up / Apply > Pick Up" from the context menu of first scatter plot viewer
    And user picks "Pick Up / Apply > Apply" from the context menu of second scatter plot viewer
    Then "xColumnName" property of second scatter plot viewer should be "Chemical Space X"
    And "colorColumnName" property of second scatter plot viewer should be "Chemical Space X"
    And the "formula lines" reading of second scatter plot viewer should be 5
    When user picks "Tools > Formula Lines..." from the context menu of second scatter plot viewer
    Then the "Formula Lines" dialog should be visible
    When user changes formula line 3 in the Formula Lines dialog to "${Chemical Space Y} = 0.3* ${Chemical Space X}"
    And user adds the formula line "${Chemical Space Y} = 0.1* ${Chemical Space X}" in the Formula Lines dialog
    And user sets the range of the selected formula line to "1" .. "3"
    And user adds the formula line "${Chemical Space Y} = 0.1* ${Chemical Space X}" in the Formula Lines dialog
    And user sets the range of the selected formula line to "5" .. "20"
    And user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And the "formula lines" reading of second scatter plot viewer should be 7
    And "formulaLines" property of second scatter plot viewer should contain "${Chemical Space Y} = 0.3* ${Chemical Space X}"
    And the formula lines of second scatter plot viewer should hold "${Chemical Space Y} = 0.1* ${Chemical Space X}" over the ranges "1..3; 5..20"
    And the "formula lines" reading of first scatter plot viewer should be 5
    When user picks "Tools > Formula Lines..." from the context menu of second scatter plot viewer
    And user clicks on the "cell 1 of visible" area of grid in "Formula Lines" dialog
    And user clicks on the "cell 2 of visible" area of grid in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of second scatter plot viewer should be 5
    When user picks "Tools > Formula Lines..." from the context menu of second scatter plot viewer
    And user clicks on the "cell 1 of visible" area of grid in "Formula Lines" dialog
    And user clicks on the "cell 2 of visible" area of grid in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "formula lines" reading of second scatter plot viewer should be 7
    When user sets "Row Source" property of second scatter plot viewer to "FilteredSelected"
    And user picks "Tools > Show Regression Line" from the context menu of second scatter plot viewer
    Then "showRegressionLine" property of second scatter plot viewer should be "true"
    And the "rows shown" reading of second scatter plot viewer should be 0
    When user picks "Series" in the "color" column selector of second scatter plot viewer
    Then "colorColumnName" property of second scatter plot viewer should be "Series"
    And no viewer of the current view should report an error
    And no errors should have been logged

  Scenario: Formula lines - a line chart split six ways stays whole with and without Multiaxes, and gets two dataframe lines
    When user clicks on line chart icon on toolbox
    Then the open tableview should have 1 line chart viewer
    When user adds "Series" to the splits of line chart viewer
    And user adds "Stereo Category" to the splits of line chart viewer
    And user adds "Core" to the splits of line chart viewer
    And user adds "R1" to the splits of line chart viewer
    And user adds "R2" to the splits of line chart viewer
    And user adds "R3" to the splits of line chart viewer
    Then "Split Column Names" property of line chart viewer should be "Series, Stereo Category, Core, R1, R2, R3"
    And the "split columns" reading of line chart viewer should be 6
    And line chart viewer should report no error
    When user sets "Multi Axis" property of line chart viewer to "true"
    Then the "multi axis" reading of line chart viewer should be "true"
    And no viewer of the current view should report an error
    When user sets "Multi Axis" property of line chart viewer to "false"
    Then the "multi axis" reading of line chart viewer should be "false"
    When user sets "xColumnName" property of line chart viewer to "Chemical Space X"
    And user sets "Y Column Names" property of line chart viewer to "Average Mass, Chemical Space Y"
    Then the "y columns" reading of line chart viewer should be "Average Mass, Chemical Space Y"
    And line chart viewer should report no error
    When user picks "Tools > Formula Lines..." from the context menu of line chart viewer
    Then the "Formula Lines" dialog should be visible
    When user clicks on "DataFrame" text in "Formula Lines" dialog
    And user adds the formula line "${Average Mass} = 0.75* ${Chemical Space X}* ${Chemical Space X} - 4 * ${Chemical Space X} -1+300" in the Formula Lines dialog
    And user types "#ff0000" into Color input in "Formula Lines" dialog
    And user selects "dashed" in Style input in "Formula Lines" dialog
    And user adds the formula line "${Chemical Space Y} = 0.75* ${Chemical Space X}* ${Chemical Space X} - 4 * ${Chemical Space X} " in the Formula Lines dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And the ".formula-lines" tag of the table should contain "#ff0000"
    And the ".formula-lines" tag of the table should contain "dashed"
    And the "formula lines" reading of line chart viewer should be 2
    And line chart viewer should draw 2 formula lines
    And no viewer of the current view should report an error
    And no errors should have been logged

  Scenario: Formula lines - a second line chart draws the dataframe lines of its Y columns, with Multiaxes and on a logarithmic axis
    When user clicks on line chart icon on toolbox
    Then the open tableview should have 2 line chart viewers
    When user sets properties of second line chart viewer:
      | xColumnName  | Chemical Space X               |
      | yColumnNames | Chemical Space Y, Average Mass |
    Then the "formula lines" reading of second line chart viewer should be 2
    And second line chart viewer should draw 2 formula lines
    When user adds "Series" to the splits of second line chart viewer
    And user sets "Multi Axis" property of second line chart viewer to "true"
    Then the "multi axis" reading of second line chart viewer should be "true"
    And the legend of second line chart viewer should list 14 items
    And the legend of second line chart viewer should be docked
    And second line chart viewer should report no error
    When user sets "yColumnNames" property of second line chart viewer to "Average Mass"
    Then the "formula lines" reading of second line chart viewer should be 1
    And second line chart viewer should draw 1 formula line
    When user sets "Y Axis Type" property of second line chart viewer to "logarithmic"
    Then "Y Axis Type" property of second line chart viewer should be "logarithmic"
    And second line chart viewer should draw 1 formula line
    And second line chart viewer should report no error
    And no errors should have been logged

  Scenario: Formula lines - renaming Chemical Space X and Y renames the calculated column and the formula lines
    When user scrolls the grid to the "Chemical Space X" column
    And user picks "Column Properties..." from the context menu of the "header Chemical Space X" area of grid
    And user types "Chem Space X" into "New name:" input in "Chemical Space X" dialog
    And user clicks on OK button in "Chemical Space X" dialog
    And user scrolls the grid to the "Chemical Space Y" column
    And user picks "Column Properties..." from the context menu of the "header Chemical Space Y" area of grid
    And user types "Chem Space Y" into "New name:" input in "Chemical Space Y" dialog
    And user clicks on OK button in "Chemical Space Y" dialog
    Then the table should have a column "${Chem Space Y} ${Chem Space X}"
    And the table should not have a column "${Chemical Space Y} ${Chemical Space X}"
    And "${Chem Space Y} ${Chem Space X}" column should have tag "formula" equal to "Qnum((100-${Chem Space Y})/100, if(qualifier(${Chem Space X})==\">\", \"<\", \"=\"))"
    And "${Spec} result" column should have tag "formula" equal to "if(${NIBR logP} != null, ${NIBR logP}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chem Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Spec}, 'Rat') || Contains(${Spec}, 'Rat Legacy'), 80, if(Contains(${Spec}, 'Mouse'), 125, if(${Spec}==\"Dog\", 30.9, if(${Spec}==\"Monkey\", 43.6, if(${Spec}==\"Minipig\", 39, null)))))*100"
    And "formulaLines" property of first scatter plot viewer should contain "${Chem Space Y} = 0.1* ${Chem Space X}"
    And no formula line of first scatter plot viewer should name "Chemical Space"
    And no formula line of second scatter plot viewer should name "Chemical Space"
    And the ".formula-lines" tag of the table should contain "${Chem Space Y} = 0.75* ${Chem Space X}"
    And no formula line of the table should name "Chemical Space"
    And the "formula lines" reading of second scatter plot viewer should be 8
    And the "formula lines" reading of first scatter plot viewer should be 6
    And the "formula lines" reading of second line chart viewer should be 1
    And second line chart viewer should draw 1 formula line
    And no viewer of the current view should report an error
    And no errors should have been logged

  Scenario: Formula lines - six more viewers get legends in different places and are stacked as the tabs of one panel
    When user clicks on histogram icon on toolbox
    And user clicks on bar chart icon on toolbox
    And user clicks on pie chart icon on toolbox
    And user clicks on trellis plot icon on toolbox
    And user clicks on pc plot icon on toolbox
    And user clicks on box plot icon on toolbox
    Then the open tableview should have 1 histogram viewer
    And the open tableview should have 1 box plot viewer
    When user sets properties of histogram viewer:
      | splitColumnName  | Stereo Category |
      | legendVisibility | Always          |
      | legendPosition   | Left            |
    And user sets properties of bar chart viewer:
      | splitColumnName  | Series          |
      | stackColumnName  | Stereo Category |
      | legendVisibility | Always          |
      | legendPosition   | Top             |
    And user sets properties of pie chart viewer:
      | categoryColumnName | Stereo Category |
      | legendVisibility   | Always          |
      | legendPosition     | Right           |
    And user sets properties of pc plot viewer:
      | colorColumnName  | Stereo Category |
      | legendVisibility | Always          |
      | legendPosition   | Bottom          |
    And user stacks bar chart viewer onto histogram viewer as a tab
    And user stacks pie chart viewer onto bar chart viewer as a tab
    And user stacks pc plot viewer onto pie chart viewer as a tab
    And user stacks box plot viewer onto pc plot viewer as a tab
    And user stacks trellis plot viewer onto box plot viewer as a tab
    Then the viewers "Histogram, Bar chart, Pie chart, PC Plot, Box plot, Trellis plot" should be the tabs of one panel
    When user switches the tabbed panel to the "Histogram" tab
    Then the legend of histogram viewer should be on the left
    And the legend of histogram viewer should list 5 items
    When user switches the tabbed panel to the "Bar chart" tab
    Then the legend of bar chart viewer should be on the top
    And the legend of bar chart viewer should list 5 items
    When user switches the tabbed panel to the "Pie chart" tab
    Then the legend of pie chart viewer should be on the right
    And the legend of pie chart viewer should list 5 items
    When user switches the tabbed panel to the "PC Plot" tab
    Then the legend of pc plot viewer should be on the bottom
    And the legend of pc plot viewer should list 5 items
    When user sets "legendVisibility" property of pc plot viewer to "Never"
    Then the legend of pc plot viewer should be placed nowhere
    When user switches the tabbed panel to the "Histogram" tab
    Then the legend of histogram viewer should be on the left
    And no viewer of the current view should report an error
    And no errors should have been logged

  Scenario: Formula lines - Chem Space X colored by conditions and then linearly colors the scatter plot the same way
    When user scrolls the grid to the "Chem Space X" column
    And user picks "Color Coding > Conditional" from the context menu of the "header Chem Space X" area of grid
    And user closes the context menu
    Then "Chem Space X" column should be color-coded conditionally
    When user sets "legendVisibility" property of first scatter plot viewer to "Always"
    And user colors "Chem Space X" column conditionally:
      | <0.1    | #1F77B4 |
      | 0.1-0.5 | #2CA02C |
      | 0.5-1   | #FFBB78 |
      | 1-5     | #D62728 |
      | >5      | #9467BD |
    Then the legend of first scatter plot viewer should list 5 items
    And the "view" area of first scatter plot viewer should contain the color "#9467BD"
    When user picks "Color Coding > Linear" from the context menu of the "header Chem Space X" area of grid
    And user closes the context menu
    Then "Chem Space X" column should be color-coded linearly
    When user colors "Chem Space X" column linearly from "#00FF00" to "#FF00FF"
    Then the "view" area of first scatter plot viewer should contain the color "#FF00FF"
    And the "view" area of first scatter plot viewer should not contain the color "#9467BD"
    And no errors should have been logged

  Scenario: Formula lines - a color picked in the pie chart legend reaches the other legends of Stereo Category
    When user switches the tabbed panel to the "Pie chart" tab
    And user hovers over "R_ONE" legend item in legend of pie chart viewer
    And user clicks on color picker icon
    Then "R_ONE" dialog should be visible
    When user picks the color "#17BECF" in the color picker dialog
    And user clicks on OK button in "R_ONE" dialog
    Then "R_ONE" dialog should be absent
    And the categorical color of "R_ONE" in "Stereo Category" column should be "#17BECF"
    And the "R_ONE" item in the legend of pie chart viewer should be colored "#17BECF"
    When user switches the tabbed panel to the "Histogram" tab
    Then the "R_ONE" item in the legend of histogram viewer should be colored "#17BECF"
    When user switches the tabbed panel to the "Bar chart" tab
    Then the "R_ONE" item in the legend of bar chart viewer should be colored "#17BECF"
    And no errors should have been logged

  Scenario: Formula lines - a saved layout brings the tabbed viewers and their legends back after an earlier layout and a move
    When user saves the layout of the current table view as "formula-legend"
    And user applies the layout "calc-columns" to the current table view
    Then the open tableview should have 0 histogram viewers
    And the open tableview should have 0 scatter plot viewers
    And "${Spec} result" filter card should be visible
    And no panel of the current view should hold tabs
    And no viewer of the current view should report an error
    And no errors should have been logged
    When user applies the layout "formula-legend" to the current table view
    Then the viewers "Histogram, Bar chart, Pie chart, PC Plot, Box plot, Trellis plot" should be the tabs of one panel
    When user switches the tabbed panel to the "Histogram" tab
    And user docks histogram viewer to the right edge of the view
    Then histogram viewer should be docked along the right edge of the view
    And the viewers "Bar chart, Pie chart, PC Plot, Box plot, Trellis plot" should be the tabs of one panel
    When user applies the layout "formula-legend" to the current table view
    Then the viewers "Histogram, Bar chart, Pie chart, PC Plot, Box plot, Trellis plot" should be the tabs of one panel
    When user switches the tabbed panel to the "Histogram" tab
    Then the legend of histogram viewer should be on the left
    And the legend of histogram viewer should list 5 items
    When user switches the tabbed panel to the "Bar chart" tab
    Then the legend of bar chart viewer should be on the top
    And the legend of bar chart viewer should list 5 items
    When user switches the tabbed panel to the "Pie chart" tab
    Then the legend of pie chart viewer should be on the right
    And the legend of pie chart viewer should list 5 items
    When user switches the tabbed panel to the "PC Plot" tab
    Then "legendVisibility" property of pc plot viewer should be "Never"
    When user switches the tabbed panel to the "Box plot" tab
    Then no viewer of the current view should report an error
    And the "formula lines" reading of first scatter plot viewer should be 6
    And no errors should have been logged

  Scenario: Formula lines - the scatter plots zoom, pan, select and filter by their legend
    When user remembers the value range of first scatter plot viewer
    And user scrolls the mouse wheel up 3 times over the "view" area of first scatter plot viewer
    Then first scatter plot viewer should show a narrower value range than before
    When user remembers the "x axis min" reading of first scatter plot viewer
    And user drags the "view" area of first scatter plot viewer by 60 pixels to the right
    Then the "x axis min" reading of first scatter plot viewer should not be as remembered
    When user double-clicks on empty plot space of first scatter plot viewer
    Then first scatter plot viewer should show the remembered value range
    When user drags a selection box over the "view" area of first scatter plot viewer
    Then some rows should be selected
    And the "rows selected" reading of first scatter plot viewer should be at least 1
    And the "labels shown" reading of first scatter plot viewer should be at least 1
    When user sets "legendVisibility" property of second scatter plot viewer to "Always"
    When user remembers the "rows shown" reading of second scatter plot viewer
    And user clicks on "Triazoles" item in the legend of second scatter plot viewer
    Then the "rows shown" reading of second scatter plot viewer should be lower than remembered
    And no errors should have been logged

  Scenario: Formula lines - the layout is saved again, changed and applied, and a copy of the project is saved as NxProjectFormulaLegend
    When user switches the tabbed panel to the "Histogram" tab
    And user saves the layout of the current table view as "formula-legend-final"
    And user remembers the place of the "view" area of histogram viewer
    And user docks histogram viewer to the right edge of the view
    Then histogram viewer should be docked along the right edge of the view
    When user applies the layout "formula-legend-final" to the current table view
    Then the viewers "Histogram, Bar chart, Pie chart, PC Plot, Box plot, Trellis plot" should be the tabs of one panel
    When user switches the tabbed panel to the "Histogram" tab
    Then histogram viewer should not be docked along the right edge of the view
    And the "view" area of histogram viewer should be placed as remembered
    And the legend of histogram viewer should be on the left
    And no viewer of the current view should report an error
    When user opens the Save project dialog from the ribbon
    When user clicks on "Save a copy" text in "Save project" dialog
    And user types "NxProjectFormulaLegend-{run}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    And 1 project named "NxProjectFormulaLegend-{run}" should be on the server
    When user closes all views

  Scenario: Filtering - the NxProjectFormulaLegend project opens, and a scaffold tree filter takes an old tree and an NX tree
    When user opens the project saved as "NxProjectFormulaLegend-{run}"
    Then table "SPGI" should be open
    And no error or warning balloon should have been shown
    And no viewer of the current view should report an error
    And no errors should have been logged
    When user clicks on the tab of the "SPGI" view
    And user clicks on filter icon in toolbar
    And user picks "Remove All" from the viewer menu of filter panel
    And user picks "Add Filter | Scaffold Tree Filter..." from the viewer menu of filter panel
    Then "Select columns..." dialog should be visible
    When user types "Structure" into "Search" input in "Select columns..." dialog
    And user toggles the "Structure" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "nodes" reading of the scaffold tree filter should be 0
    When user uploads "fixtures/nx/scaffold-tree-for-testing.tree" through "Upload saved tree file" icon inside filter panel
    Then the "nodes" reading of the scaffold tree filter should be at least 1
    And no node of the scaffold tree filter should be colored
    When user uploads "fixtures/nx/scaffold-tree-for-nx-testing.tree" through "Upload saved tree file" icon inside filter panel
    Then the "colored nodes" reading of the scaffold tree filter should be at least 1
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Filtering - checked scaffolds filter every view, the layout keeps them, and closing the panel releases them
    Given the "hits of node 1" reading of the scaffold tree filter should be at least 1
    And the "hits of node 2" reading of the scaffold tree filter should be at least 1
    When user remembers how many rows pass the filter
    And user clicks on the "checkbox of node 1" area of the scaffold tree filter
    Then the "checked nodes" reading of the scaffold tree filter should be 1
    And fewer rows than remembered should pass the filter
    And the "rows kept" reading of the scaffold tree filter should be at least 1
    When user clicks on the "checkbox of node 2" area of the scaffold tree filter
    Then the "checked nodes" reading of the scaffold tree filter should be 2
    And the "bit operation" reading of the scaffold tree filter should be "OR"
    When user remembers the "rows kept" reading of the scaffold tree filter
    And user sets the scaffold tree filter to combine the checked scaffolds with "AND"
    Then the "bit operation" reading of the scaffold tree filter should be "AND"
    And the "rows kept" reading of the scaffold tree filter should be lower than remembered
    When user remembers how many rows pass the filter
    And user clicks on the tab of the "SPGI copy" view
    Then grid should show every row that passes the filter of its table
    When user clicks on the tab of the "SPGI" view
    And user saves the layout of the current table view as "filtering"
    And user clicks on filter icon in toolbar
    Then filter panel should be absent
    And more rows than remembered should pass the filter
    When user applies the layout "filtering" to the current table view
    Then as many rows as remembered should pass the filter
    And the "checked nodes" reading of the scaffold tree filter should be 2
    And no errors should have been logged

  Scenario: Filtering - a molecule used as a filter puts a Structure card on the panel, and every view shows the rows left
    When user clicks on the tab of the "SPGI (2)" view
    And user remembers how many rows pass the filter
    And user picks "Current Value > Use as filter" from the context menu of the drawn cell of "Structure" column with the longest value
    Then filter panel should be visible
    And the "type of Structure" reading of filter panel should be "Chem:substructureFilter"
    And the "filtering of Structure" reading of filter panel should be "true"
    And the Chem filters of every view should have finished computing
    And no more rows than remembered should pass the filter
    When user clicks on the tab of the "SPGI" view
    Then grid should show every row that passes the filter of its table
    And no errors should have been logged

  @known-failure
  Scenario: Filtering - the rows a molecule used as a filter lets through contain that molecule
    Then every row that passes the filter should contain the molecule of the cell picked

  Scenario: Filtering - a structure changed in a clone's card reaches the first view's card
    When user clicks on the tab of the "SPGI (2)" view
    And user picks "View > Layout > Clone View" from the top menu
    And user clicks on the structure drawn in the "Structure" filter card
    Then sketcher dialog should be visible
    When user clears molecule input of sketcher dialog
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then the "structure of Structure" reading of filter panel should be "c1ccncc1"
    When user clicks on the tab of the "SPGI (2)" view
    Then the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And no errors should have been logged

  Scenario: Filtering - Filter type Categorical in the Chemistry Rendering pane makes the Core card categorical
    When user scrolls the grid to the "Core" column
    And user clicks on the "header Core" area of grid
    Then the context panel should show "Core"
    When user sets Filter type to "Categorical" in the Rendering pane of the context panel
    And user adds a card for "Core" to the filter panel
    Then the "cards" reading of filter panel should include the text "Core"
    And the "type of Core" reading of filter panel should be "categorical"
    And no errors should have been logged

  Scenario: Filtering - in a clone the master switch and a card's own switch turn their filtering off and on
    When user picks "View > Layout > Clone View" from the top menu
    And user remembers how many rows pass the filter
    And user clicks on the "category S_PART of Stereo Category" area of filter panel
    Then fewer rows than remembered should pass the filter
    And no rows where "Stereo Category" is "R_ONE" should pass the filter
    When user remembers how many rows pass the filter
    And user hovers over filter panel
    And user unchecks master of filter panel
    Then the "active" reading of filter panel should be "false"
    And more rows than remembered should pass the filter
    When user remembers how many rows pass the filter
    And user hovers over filter panel
    And user checks master of filter panel
    Then the "active" reading of filter panel should be "true"
    And fewer rows than remembered should pass the filter
    And no rows where "Stereo Category" is "R_ONE" should pass the filter
    When user remembers how many rows pass the filter
    And user hovers over "Stereo Category" filter card
    And user unchecks checkbox of "Stereo Category" filter card
    Then "Stereo Category" filter card should be disabled
    And more rows than remembered should pass the filter
    And no errors should have been logged

  Scenario: Filtering - the layout saved by the calc-columns scenario applies its filtering
    When user applies the layout "calc-columns" to the current table view
    Then "${Spec} result" filter card should be visible
    And no viewer of the current view should report an error
    And no errors should have been logged

  @known-failure
  Scenario: Filtering - after a layout saved before Chemical Space X was renamed, the calculated column's formula names only existing columns
    Then every column the formula of "${Spec} result" column refers to should exist
    And "${Spec} result" column should have tag "formula" equal to "if(${NIBR logP} != null, ${NIBR logP}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chem Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Spec}, 'Rat') || Contains(${Spec}, 'Rat Legacy'), 80, if(Contains(${Spec}, 'Mouse'), 125, if(${Spec}==\"Dog\", 30.9, if(${Spec}==\"Monkey\", 43.6, if(${Spec}==\"Minipig\", 39, null)))))*100"

  Scenario: Filtering - the layout saved in this scenario applies its filtering, and a copy is saved as NxProjectFiltering
    When user applies the layout "filtering" to the current table view
    Then the "checked nodes" reading of the scaffold tree filter should be 2
    And "${Spec} result" column should have tag "formula" equal to "if(${NIBR logP} != null, ${NIBR logP}, if(${Route Admin}==\"PO\", ${Whole blood assay 1} / ${Chem Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Spec}, 'Rat') || Contains(${Spec}, 'Rat Legacy'), 80, if(Contains(${Spec}, 'Mouse'), 125, if(${Spec}==\"Dog\", 30.9, if(${Spec}==\"Monkey\", 43.6, if(${Spec}==\"Minipig\", 39, null)))))*100"
    And no viewer of the current view should report an error
    And no errors should have been logged
    And the Chem filters of every view should have finished computing
    When user remembers what the filter panel of every view filters by
    And user opens the Save project dialog from the ribbon
    When user clicks on "Save a copy" text in "Save project" dialog
    And user types "NxProjectFiltering-{run}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    And 1 project named "NxProjectFiltering-{run}" should be on the server
    When user closes all views

  Scenario: Filtering - the NxProjectFiltering project opens with its filtering
    When user opens the project saved as "NxProjectFiltering-{run}"
    Then table "SPGI" should be open
    And no error or warning balloon should have been shown
    And the filter panel of every view should filter by what was remembered
    And the Chem filters of every view should have finished computing
    And no viewer of the current view should report an error
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Cleanup - the five projects of the chain are deleted
    When user closes all views
    And user deletes the projects "NxProject-{run}, NxProjectCalcColumns-{run}, NxProjectViewers-{run}, NxProjectFormulaLegend-{run}, NxProjectFiltering-{run}"
    Then none of the projects "NxProject-{run}, NxProjectCalcColumns-{run}, NxProjectViewers-{run}, NxProjectFormulaLegend-{run}, NxProjectFiltering-{run}" should be on the server
