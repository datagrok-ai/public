# Forms viewer — what the features cannot say yet

Written during the viewer-gaps round of 23 Oct 2026, against `files/TestTrack/Viewers/FormsViewer/`
and the viewer's "Automation surface" (`public/libraries/utils/src/viewers/CLAUDE.md`). New this round:
`forms-cards-drawn.feature` (`forms.md` Step 4c — each selected row's card draws its molecule and curve).

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
