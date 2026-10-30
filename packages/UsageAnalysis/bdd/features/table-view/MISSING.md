# Table view features — what is still missing

One section per group of the table-view round (30 Oct 2026). Each lists the scenarios that need a step
the library lacks (as Gherkin, with the signal it would read), the questions for the help or the product,
and the findings.

## Data and Select menus

Features: `data-join-tables.feature`, `data-unpivot.feature`, `data-text-to-columns.feature`,
`select-menu.feature`, `edit-add-rows.feature` (30 Oct 2026, local stand).

### Needs a step (Data, Select)

- **Select > Selection to Top** sorts the grid so that the selected rows come first (`selection_to_top.dart`).
  No step reads the grid's row order (every row step reads the table, which the sort leaves as it was).
  The step must read `grid.tableRowToGrid` / the grid's sorted order, e.g. `the grid should show rows where
  "<col>" is "<v>" first`:

  ```gherkin
  Scenario: Selection to Top brings the selected rows to the top of the grid
    Given user selects rows where "name" is "c"
    When user picks "Select > Selection to Top" from the top menu
    Then the first grid row should show "c" in "name"
    And the value of "name" column in row 5 should be "c"
    And no errors should have been logged
  ```

- **Join Tables: picking a key pair other than the default.** The library's key-pair step is scoped to the
  Link Tables dialog (`user sets key columns {int} of the Link Tables dialog to …`, `[name="dialog-Link-Tables"]`);
  the Join Tables dialog has the same `div-selectKeyCol{side}Row{n}` pickers. The step should take the
  dialog title (`… of the {string} dialog`), so the Join features can add a second pair (composite keys) and
  pick keys of different types (the test page's "different types of key columns" → OK greyed out with a
  reason):

  ```gherkin
  Scenario: Key columns of different types disable OK and say why
    When user picks "Data > Join Tables..." from the top menu
    And user sets key columns 1 of the "Join Tables" dialog to "name" and "id"
    Then OK button in Join Tables dialog should be disabled
    And Join Tables dialog should contain text "<the reason>"
    And no errors should have been logged
  ```

### Questions (Data, Select)

- **Join Tables, "OK disabled while keys are missing"**: the first key pair has no remove button (only the
  `+`), so the keys can never be all missing from the UI; with no pair (`keysGrid` empty) the code joins by
  row index ("Joining by row indexes…") with OK enabled. Is the row-index join reachable at all, and should it
  be claimed?
- **Join Tables, In-place**: a left join done in place renames the left table to "left+right" (measured). The
  help says the result "replaces the content of the left table" and says nothing about the name. Intended?
  The feature claims the measured name.
- **Unpivot: the role inputs carry no label.** The per-column Copy/Merge selects are named
  `input-<column>` but have no label, aria-label or title, so a phrase can reach them only by position
  (`3rd input in Unpivot dialog`). A screen reader hears an unnamed combo box. Should each select take its
  column's name as its accessible name (the text shown next to it)? Then the feature reads
  `"math" input in Unpivot dialog`.
- **Split Column and Unpivot: undo.** Neither help page (`text-to-columns.md`, `unpivot.md`) nor
  `datagrok/navigation/undo.md` says whether Ctrl+Z removes the columns Split Column adds. Not claimed.

### Findings 30 Oct 2026 (Data, Select)

- **Split Column: the default Split by is ";," and splits nothing** ([core#1236](https://github.com/datagrok-ai/core/issues/1236)). Steps: open a
  table with a text column holding "a;b"; Data > Split Column...; leave Split by at its default
  ";,", no By chars; OK. Measured: the default is taken as one two-character delimiter, so "a;b" is not split
  (col1 = "a;b"). Expected: the default splits on either character. The help (`text-to-columns.md`) does not say.

## Column operations, undo, properties

Features: `undo-redo.feature`, `change-type.feature`, `column-properties.feature`, `table-properties.feature`.

### Needs a step (column operations)

- **Change column type: the Convertible / Not convertible / Empty counts.** The dialog shows each count as a
  bare `label.ui-input-label` beside the category's label, with no name or role tying the two together, so no
  phrase reads "Not convertible: 1". The step needs a readout (`[role="status"]` with `<label>` + `<output>`)
  per row in `column_commands.dart` `category()`, then:
  ```gherkin
  When user picks "Change Type..." from the context menu of the "header code" area of grid
  And user selects "int" in "New type" input in "Change column type" dialog
  Then the value of "Convertible" readout in "Change column type" dialog should be "2"
  And the value of "Not convertible" readout in "Change column type" dialog should be "1"
  And the value of "Empty" readout in "Change column type" dialog should be "1"
  ```

### Questions (column operations)

- **Table tags and description.** `datagrok/concepts/table.md` says "each table can be annotated with
  attributes", but no dialog or context-panel pane edits a table's tags or description (the context panel shows
  General and Actions; Actions has Rename... and Clone). Is there meant to be a Table Properties editor? Not
  claimed.

- **Ctrl+Z after closing another view reopens that view, ahead of this view's own change: current behaviour,
  not a bug** ([core#1234](https://github.com/datagrok-ai/core/issues/1234)). In table view "undo" remove a column,
  open and close a second table's view, then Ctrl+Z in "undo": the closed view comes back first, and a second
  Ctrl+Z restores the column. Davit Rizhinashvili, 30 Oct: undo is global, and table views and view commands are
  not scoped yet; this is the expected behaviour today, scoping to the table view is planned. The help's "never one
  you performed somewhere else" (`datagrok/navigation/undo.md`) overstates it. Not claimed either way until the
  scoping lands; `undo-redo.feature` switches views without closing one, so it does not depend on it.

### Findings 30 Oct 2026 (column operations)

- **Edit > Add New Column: Ctrl+Z does not remove the added column** ([public#4461](https://github.com/datagrok-ai/public/issues/4461)). Steps: open a table with columns
  name, city, weight; Edit > Add New Column...; formula `${weight} * 2`, name `double`; OK; click a grid cell;
  press Ctrl+Z (or Edit > Undo). Measured: the column stays (`name, city, weight, double`),
  `grok.shell.canUndo` is false right after OK. Expected: the column is removed, and Ctrl+Shift+Z adds it back.
  Help: `datagrok/navigation/undo.md`, "What can be undone", Columns: "adding or editing a formula column".
  Cause seen in the source: the core Dart dialog (`edit_add_new_column.dart`) pushes an `UndoRecord` "Add
  column ...", but the menu opens PowerPack's `AddNewColumnDialog` (`dialogs/add-new-column.ts`), which runs
  `AddNewColumn` / `EditColumnFormula` and pushes no record. Failed in runs r2, w1a and a direct probe. The
  scenario, as it was in `undo-redo.feature`:
  ```gherkin
  Scenario: A formula column added through Add New Column is removed by Ctrl+Z and added again by redo
    When user picks "Edit > Add New Column..." from the top menu
    And user pastes "${weight} * 2" into code editor in "Add New Column" dialog
    And user types "double" into column name input
    And user clicks on OK button in "Add New Column" dialog
    Then the "Add New Column" dialog should close
    And the value of "double" column in row 2 should be "4.5"
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Z
    Then the table should not have a column "double"
    And the table should have 3 columns
    When user presses Control+Shift+Z
    Then the value of "double" column in row 3 should be "6"
    When user presses Control+Z
    Then the table should not have a column "double"
    And no errors should have been logged
  ```
- **Change column type: the "select" icon beside each count has no glyph** ([core#1220](https://github.com/datagrok-ai/core/issues/1220)). Steps: open a table with a string
  column code = 10, 20, n/a, (empty); right-click its header > Change Type...; New type int. Measured: each of
  the three rows (Convertible / Not convertible / Empty) shows Previous, Next and Filter icons, and a fourth
  button whose `<i aria-label="Select All">` has a 0×0 box — nothing is drawn, so there is nothing to see or
  click on its glyph. Cause: `core/client/d4/lib/src/widgets/icons.dart:30`
  `Icons.selectAll` uses `faSolid('th-lagre', …)` — a typo for `th-large`. Expected: a visible select icon
  that selects the rows of its category. Help: `transform/change-column-type.md`: "Next to each count, use the
  icons to go to the previous or next such row, filter the rows, or select them." Failed in runs r2 and w1a
  (Playwright: "element is not visible") and a direct probe. The claim, as it was in `change-type.feature`:
  ```gherkin
    When user clicks on third "Select All" icon in "Change column type" dialog
    Then 1 row should be selected
    And all rows where "name" is "d" should be selected
  ```

## Shortcuts, Column Manager, Files, Plugins

Features: `shortcuts.feature`, `column-manager.feature` (this folder), `../browse/files-lifecycle.feature`
(30 Oct 2026, local stand). No plugins feature: see the questions below.

### Needs a step (shortcuts, Column Manager, Files)

- **Ctrl+Shift+R shows "Ctrl+Z to undo".** `cmdViewLayoutClear` (`commands/view/layout.dart`) calls `info()`,
  which writes the text into the top-right `label.grok-info` for 5 s (`features/info.dart`), not a balloon, so
  `an info balloon containing …` never sees it. The same label carries "N rows removed", "Block copied" and the
  viewer-close hint. The step reads that label's text while it is shown:
  ```gherkin
    When user presses Control+Shift+R
    Then the info label should show "Ctrl+Z to undo"
  ```
- **F11 and Alt+F (full screen).** F11 requests full screen for the document (`cmdViewFullScreen`), Alt+F for the
  current viewer ("Full screen" viewer command). Nothing reads `document.fullscreenElement`, and the page is
  shared with the next features, so a step must also leave full screen at scenario end:
  ```gherkin
    When user presses F11
    Then the page should be in full screen
    When user presses Escape
    Then the page should not be in full screen
  ```
- **F1: the help follows the current object.** The help panel shows the page of what is current (a viewer, a
  column); on a dev stack no help pages are served, so the claim needs `the stand serves the help pages` and a
  step that reads the help panel's page (`the help panel should show the "<page>" page`).
- **Column Manager: a command applied to "all selected tables".** `column-manager.md`: "If a context menu is open
  when multiple columns are selected, you can choose to apply commands to either the current table or all
  selected tables." That needs several tables selected in the Tables pane (Alt+T), which no step selects
  (`user selects the {string} table in the Tables pane` holding Control).
- **Ctrl+L: closing the suggestion panels again.** Ctrl+L does not toggle: `cmdViewLayoutLoad` docks two new
  panels on every press, into the shell's dock (not the view's), so they outlive the table view. Their title-bar
  close button is `div.panel-titlebar-button-close[name="Close"]` (`panel_dock_container.dart`), which no kind
  reaches: `user closes "Layout suggestions" dock panel` falls back to Escape and the panel stays (run rf1).
  `shortcuts.feature` therefore leaves them open. The dock panel kind needs a `close icon` part, then:
  ```gherkin
    When user clicks on the close icon of "Layout suggestions" dock panel
    And user clicks on the close icon of "Viewer suggestions" dock panel
    Then "Layout suggestions" dock panel should be absent
    And "Viewer suggestions" dock panel should be absent
  ```

### Questions (shortcuts, Column Manager, Files)

- **Ctrl+L is "Open Gallery" in the help** (`shortcuts.md`: "Table View | Layouts | Open Gallery"), and the
  command's top menu path is View > Layout > Open Gallery, but what it does is dock two panels, "Viewer
  suggestions" and "Layout suggestions" (`cmdViewLayoutLoad`); there is no gallery view. Should the help say
  "Show layout and viewer suggestions"? `shortcuts.feature` claims the two panels.
- **`column-manager.md` says Column Manager is "available via View | Columns or Add | Column Viewer".** View >
  Columns (Alt+C) docks it; the viewer is registered (`ViewerDescriptor.columnManager`), but no Add menu
  item named Column Viewer was checked. Should the feature claim adding it as a viewer? Not claimed. (`shortcuts.md` also
  spells it "Colum Manager".)
- **Plugins: switching a package's version and Uninstall are not written** (`browse/plugins-versions.feature`
  was planned). Measured on the local stand: a package card's menu (Browse > Platform > Plugins) offers
  Uninstall; in its Versions pane the current version has no menu, an installed older version offers "Set
  current" and "Uninstall", a version not installed offers "Install", the account's debug version offers
  "Delete". Why not written:
  - Uninstall removes the package from the stand; putting it back means installing it from the npm registry (an
    outside service, never in a feature) or republishing it from source (no step, and minutes of build). It
    cannot be restored within the feature.
  - "Set current" can be undone by "Set current" on the version that was current, but it changes the version
    every account and every parallel feature loads (a package other features use), the admin account the suite
    runs as is served its debug version anyway (so nothing visible changes for it but the icon), and it needs a
    stand where a second version is installed (the CI stand installs one). The outcome (which version is
    current) is a server state — an API/datlas test of the packages endpoint.
  Should a version switch be a feature at all? If yes, on what fixture package (one no feature and no user
  needs, with two versions published by the feature's own setup), and with which step to publish it?

## Batch Edit, Anonymize, Categorize, Compare Tables

Features: `batch-edit.feature`, `anonymize.feature`, `categorize.feature`, `compare-tables.feature` (30 Oct 2026,
local stand). All four dialogs are core (`xamgle`): `views/batch_edit_view.dart`, `commands/tools/tools_anonymize_data.dart`,
`commands/data/data_commands.dart` + `tools/data_categorizer.dart`, `commands/data/compare_tables.dart`.

### Needs a step or a name in core (Batch Edit, Anonymize, Categorize, Compare)

- **Categorize: typed patterns.** The Pattern and Category editors of each row are bare `ui.textInput`s with no
  name, label or aria-label (only a tooltip), so no phrase reaches them; the feature drives the dialog through
  Bin... only. With an `aria-label="Pattern"` / `"Category"` per row (`data_categorizer.dart`, the `catTable`
  renderer), the `text input` kind would find them as `first "Pattern" text input`:

  ```gherkin
  Scenario: A regular expression pattern on a text column marks the matching rows with its category
    When user picks "Data > Categorize..." from the top menu
    And user selects "city" in "Categorize" input in Categorize Data dialog
    And user enters "Os" into first "Pattern" text input in Categorize Data dialog
    And user enters "North" into first "Category" text input in Categorize Data dialog
    Then Categorize Data dialog should contain text "Matching: 2"
    When user clicks on "SAVE COLUMN" button in Categorize Data dialog
    Then the value of "category" column in row 1 should be "North"
    And the value of "category" column in row 2 should be ""
    And no errors should have been logged
  ```

- **Anonymize: the "category" column.** A text column is mapped onto the categories of another open table with
  exactly one string column by picking that table in the column's "category" cell of the dialog's grid (a choice
  cell editor). No step picks a choice in a grid cell inside a dialog; it needs one, e.g.
  `user picks {string} in the {string} area of grid in {element}` (double-click the hit area, choose in the
  cell editor's list):

  ```gherkin
  Scenario: A category table maps a text column onto its categories
    Given user opens a table "towns" with:
      | town |
      | X    |
      | Y    |
    And user switches to the "people" table view
    When user picks "Data > Anonymize..." from the top menu
    And user toggles the "city" column in the column list of Anonymize Data dialog
    And user picks "towns" in the "cell 2 of category" area of grid in Anonymize Data dialog
    And user clicks on OK button in Anonymize Data dialog
    Then the value of "city" column in row 1 should be "X"
    And the value of "city" column in row 2 should be "Y"
    And no errors should have been logged
  ```

### Questions (Batch Edit, Anonymize, Categorize, Compare)

- **Batch Edit: undo.** Neither `transform/batch-edit.md` nor `datagrok/navigation/undo.md` says whether Ctrl+Z
  takes back a value written by the batch editor (`Column.setAllValues`). Not claimed.
- **Categorize: undo.** The help and `undo.md` say nothing about Ctrl+Z after SAVE COLUMN / SAVE MATRIX /
  APPLY COLOR CODING. Not claimed.

### Findings 30 Oct 2026 (Batch Edit, Anonymize, Categorize, Compare)

- **Anonymize: an invalid Number randomization factor throws instead of being refused** ([core#1225](https://github.com/datagrok-ai/core/issues/1225)). Steps: open any table
  with a number column (e.g. demog); Data > Anonymize...; check AGE in the column list; type `abc` into Number
  randomization factor; click OK. Measured: OK stays enabled and the input is not marked invalid; the click logs
  `FormatException: Invalid double abc` (console error, `double.parse` in `anonymize()`), the dialog stays open
  with no message, and the data is unchanged. Expected: the factor validated like other numeric dialog inputs
  (the Categorize "Number of bins" dialog marks its input `d4-invalid` and disables OK), so OK cannot be pressed,
  and no error is logged. Help: `transform/anonymize-data.md` says only that "the scale of the change depends on
  the specified number randomization factor r"; the test page (`anonymize-data-test.md`) asks to try "different
  values … (negatives, integers, floats)". Seen in a direct probe and a BDD run (the dialog never closed). The
  scenario, as it was drafted:
  ```gherkin
  Scenario: An invalid randomization factor disables OK
    When user picks "Data > Anonymize..." from the top menu
    And user toggles the "w" column in the column list of Anonymize Data dialog
    And user enters "abc" into "Number randomization factor" input in Anonymize Data dialog
    Then OK button in Anonymize Data dialog should be disabled
    And no errors should have been logged
  ```
- **Compare Tables: new and deleted rows are painted in each other's color** ([core#1231](https://github.com/datagrok-ai/core/issues/1231), confirmed on dev 30 Oct).
  Steps: open t1 (name, city: Ann/Oslo, Bob/Rome, Cid/Kyiv) and t2 (Ann/Oslo, Bob/Paris, Dan/Lviv); Data >
  Compare Tables, key name, OK. Measured: in the diff grid Cid (only in t1, deleted) is green (#CFFFC5) and Dan
  (only in t2, new) is pink (#FFB3B0). Expected: the reverse — `compare_tables.dart:210-222` defines
  `newRowColor` green and `deletedRowColor` pink but paints the `// new row` branch with `deletedRowColor` and the
  `// deleted row` branch with `newRowColor`; table 1 is the old one (the changed cell's tooltip says "Old value:"
  with t1's value). The help says only "color-coded … added, removed, or changed". `compare-tables.feature`
  claims the two colors differ, not which is which.
- **Batch Edit: the editor has no title** ([core#1235](https://github.com/datagrok-ai/core/issues/1235)). Steps: open any table; Data > Batch Edit....
  Measured: the docked panel shows no caption, and `grok.shell.views` lists it with the name `""`, so a phrase
  can only reach its rows unscoped (`"city" property`). Expected: a caption, "Batch Edit".
