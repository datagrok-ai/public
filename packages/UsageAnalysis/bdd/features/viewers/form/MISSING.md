# Form — what the features cannot say yet

Written 2026-10-22 when `form.feature` took the manual sections "Column selector", "Column changes
reaction" and "Table switching" of TestTrack `Viewers/Form/form.md`. No library, binding or code
change was made.

## Still open

1. **The text a field's label shows.** The rename of AGE to AGE_NEW is claimed through the
   `label AGE_NEW` area, which is keyed by the column the label belongs to, not by what it shows.
   Wanted: a `label text of <column>` reading in the form's status, then the feature claims the
   label reads "AGE_NEW".

## Edge cases (`form-edge-cases.feature`)

Added with `form-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Decided and claimed in `form-edge-cases`**:
   - A file that is not a form (a CSV file, a saved table layout) is refused with the error balloon
     "<file> is not a valid form file" and the form keeps its fields (core#1138, fixed).
   - Text typed into a numeric field of an editable form leaves the cell as it was: core#1192,
     fixed (the cell was emptied on the first keystroke).
   - A column removed from the table takes its field off the form. A default form does; a form built
     from chosen columns (`columnNames`) keeps an empty field: a known failure with no ticket yet.
     The scenario that claimed the empty field is gone from `form.feature`.
2. **Fixture**: the reported case of core#1138 is a file holding `{}`; a dedicated
   `fixtures/form-not-a-form.json` would say it more plainly.
3. **Steps requested**: the arrows at the first and the last row flash the arrow (`ui.noMoreMatches`)
   and report nothing; the scenario claims only that the row stays. Every row filtered out leaves
   empty fields (the no-rows rule of `viewer-rules-edge-cases.feature`); claiming it needs a reading
   of a field's value (`value of <field>`), which the form does not report.
4. **Dart test to write**: metadata changed ten times, then one keystroke in edit mode, writes the cell
   once (the field listeners are added again on every metadata change).
