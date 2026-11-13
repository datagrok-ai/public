# Forms viewer — what the features cannot say yet

Written during the viewer-gaps round of 23 Oct 2026, against `files/TestTrack/Viewers/FormsViewer/`
and the viewer's "Automation surface" (`public/libraries/utils/src/viewers/CLAUDE.md`); trimmed in the BDD PR
round of 8 Nov 2026.

## Still open

- **public#4452** — a forms layout naming a column the table no longer has throws "TypeError: Cannot read properties of
  null (reading 'semType')" on every render (`buildForm`); claim when fixed: the layout loads with fields "AGE, SEX".
- **The field takes the grid column's alignment** (`forms.md` Step 1b, GROK-19339): the `align of <COL> of <label>`
  reading exists; needs a step that sets a grid column's horizontal alignment the way a person does.

## Won't do

- **Each card draws that row's molecule** (`forms.md` Step 4c): what a cell renderer draws is the renderer's package
  test; `forms-cards-drawn.feature` claims which row each card shows.
- A pinned row under a filter that hides it, a pinned row deleted: PowerGrid `Forms viewer` "A pinned row stays pinned
  under a filter that hides it and leaves the pinned pane once deleted".
- No rows, rows deleted, all columns removed or empty, `fixedRowNumbers` before Show Current Row, "molecules on the
  left of the card": the shared viewer rules or no gap.
