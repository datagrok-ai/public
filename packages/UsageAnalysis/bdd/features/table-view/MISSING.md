# Table view features — what is still missing

What the table-view features (`undo-redo`, `data-menu`, `select-menu`, `shortcuts`, `../browse/files-lifecycle`)
cannot claim until a step, a reading or a name exists. The findings of the 30 Oct 2026 round are filed:
[public#4461](https://github.com/datagrok-ai/public/issues/4461) (Add New Column not undone),
[core#1220](https://github.com/datagrok-ai/core/issues/1220) (Change Type select icon has no glyph),
[core#1225](https://github.com/datagrok-ai/core/issues/1225) (Anonymize: a non-numeric factor throws),
[core#1231](https://github.com/datagrok-ai/core/issues/1231) (Compare Tables: new and deleted rows swap colors),
[core#1235](https://github.com/datagrok-ai/core/issues/1235) (Batch Edit has no title),
[core#1236](https://github.com/datagrok-ai/core/issues/1236) (Split Column's default `;,` splits nothing),
[public#4465](https://github.com/datagrok-ai/public/issues/4465) (help pages: wrong menu paths, table metadata).

## Needs a step

- **Join Tables: picking a key pair.** `user sets key columns {int} of the Link Tables dialog to …` is scoped to
  `[name="dialog-Link-Tables"]`; the Join Tables dialog has the same `div-selectKeyCol{side}Row{n}` pickers. With
  the dialog title as a parameter (`… of the {string} dialog`), keys of different types can be claimed:

  ```gherkin
    When user picks "Data > Join Tables..." from the top menu
    And user sets key columns 1 of the "Join Tables" dialog to "name" and "id"
    Then OK button in Join Tables dialog should be disabled
    And Join Tables dialog should contain text "<the reason>"
  ```

- **Ctrl+Shift+R shows "Ctrl+Z to undo".** `cmdViewLayoutClear` (`commands/view/layout.dart`) writes it into the
  top-right `label.grok-info` for 5 s (`features/info.dart`), not a balloon; the same label carries "N rows removed"
  and "Block copied". A step reads that label while it is shown: `Then the info label should show "Ctrl+Z to undo"`.

- **F1: the help follows the current object.** Needs `the stand serves the help pages` and a step that reads the
  help panel's page (`the help panel should show the "<page>" page`).

- **Column Manager: a command applied to "all selected tables"** (`column-manager.md`). Needs several tables
  selected in the Tables pane (Alt+T), which no step does (`user selects the {string} table in the Tables pane`
  holding Control).

- **Ctrl+L: closing the suggestion panels.** Ctrl+L docks "Layout suggestions" and "Viewer suggestions" into the
  shell's dock (not the view's) on every press, so they outlive the table view and the feature. Their close button
  is `div.panel-titlebar-button-close[name="Close"]` (`panel_dock_container.dart`), which no kind reaches. With a
  `close icon` part on the dock panel kind, the Ctrl+L scenario goes back into `shortcuts.feature`:

  ```gherkin
    When user presses Control+L
    Then "Layout suggestions" dock panel should be visible
    And "Viewer suggestions" dock panel should be visible
    When user clicks on the close icon of "Layout suggestions" dock panel
    And user clicks on the close icon of "Viewer suggestions" dock panel
    Then "Layout suggestions" dock panel should be absent
  ```

## Needs a name in core

- **Unpivot: the role inputs carry no label.** The per-column Copy/Merge selects are `input-<column>` with no
  label, aria-label or title; `data-menu.feature` reaches them by position (`3rd input in Unpivot dialog`). With
  the column's name as the accessible name: `"math" input in Unpivot dialog`.

- **Change column type: the Convertible / Not convertible / Empty counts.** Each count is a bare
  `label.ui-input-label` with no name tying it to its category. A readout per row (`[role="status"]` with
  `<label>` + `<output>`, `column_commands.dart` `category()`) lets `undo-redo.feature` claim
  `the value of "Not convertible" readout in "Change column type" dialog should be "1"`.

- **Categorize: typed patterns.** The Pattern and Category editors of each row are bare `ui.textInput`s with
  only a tooltip (`data_categorizer.dart`, the `catTable` renderer). With `aria-label="Pattern"` / `"Category"`:
  `When user enters "Os" into first "Pattern" text input in Categorize Data dialog`.
