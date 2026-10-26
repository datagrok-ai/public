# Forms viewer — what the features cannot say yet

Written during the viewer-gaps round of 23 Oct 2026, against `files/TestTrack/Viewers/FormsViewer/`
and the viewer's "Automation surface" (`public/libraries/utils/src/viewers/CLAUDE.md`). New this round:
`forms-cards-drawn.feature` (`forms.md` Step 4c — every selected row gets a card).

## Still open

- **The field takes the grid column's alignment** (`forms.md` Scenario 1, Step 1b, GROK-19339). The
  reading exists — `align of <COL> of <label>` — and `forms-presentation.feature` claims the default
  (`start`). What is missing is a step that sets a grid column's horizontal alignment the way a person
  does (the column's Style / format pane, or `user sets the alignment of {string} column of grid to
  {word}` writing `GridColumn.cellStyle.horzAlign`). Then:

  ```gherkin
  Scenario: A centred grid column centres its field on every card
    When user sets the alignment of "AGE" column of grid to "center"
    Then the "align of AGE of card 1" reading of forms viewer should be "center"
    And the "align of SEX of card 1" reading of forms viewer should be "start"
    And no errors should have been logged
  ```

- **"Molecules occupy the left part of the card"** (`forms-ui.md`, exploratory) is not a gap but a
  question: the card stacks its fields in one column (Structure above Primary Series Name, same x), so
  the md's "text fields to its right" does not describe the viewer. Nothing to write until the owner
  says which layout is intended.

- **Each card draws that row's molecule** (`forms.md` Scenario 4, Step 4c) is not a gap: what a cell
  renderer draws, on a card or in the grid, is the renderer's package test, not a feature (the
  library's "What never becomes a feature"). `forms-cards-drawn.feature` claims which row each card
  shows (`record of card <n>`).

## Edge cases (`forms-edge-cases.feature`)

Added with `forms-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Questions for owners (candidate defect, no ticket)**: observed on localhost: fields AGE, HEIGHT,
   SEX; layout saved; HEIGHT removed (the field goes, `fields` = "AGE, SEX"); the saved layout loaded:
   `fields` is "AGE, HEIGHT, SEX" again, naming a column the table does not have. A column that leaves
   the table otherwise takes its field with it (`forms-fields.feature`). Is the field kept on purpose?
   (A calculated column cannot show it: the layout recreates the column.)
2. **Questions for owners**: a pinned row filtered out (stays pinned or hides); a row pinned by an
   empty cell (stored as the text "null").
3. **Package tests to write (PowerGrid)**: the Data Manipulation set the Dart viewers get (remove all
   columns, all values empty, every row filtered out) for the Forms viewer, which is a JS viewer the
   Dart regTest does not run; `fixedRowNumbers` set before Show Current Row.
