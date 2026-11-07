# Form — what the features cannot say yet

Written 2026-10-22 when `form.feature` took the manual sections "Column selector", "Column changes
reaction" and "Table switching" of TestTrack `Viewers/Form/form.md`. No library, binding or code
change was made.

## Still open

1. **The text a field's label shows.** The rename of AGE to AGE_NEW is claimed through the
   `label AGE_NEW` area, which is keyed by the column the label belongs to, not by what it shows.
   Wanted: a `label text of <column>` reading in the form's status, then the feature claims the
   label reads "AGE_NEW".

## Edge cases (now in `form.feature`)

Added in the edge-case round of 2026-10-26 (localhost); its scenarios are now `@edge-cases` scenarios of their own, after the master scenarios of `form.feature`. Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Decided and claimed in `form.feature`**:
   - A file that is not a form (a CSV file, a saved table layout) is refused with the error balloon
     "<file> is not a valid form file" and the form keeps its fields (core#1138, fixed).
   - Text typed into a numeric field of an editable form leaves the cell as it was: core#1192,
     fixed (the cell was emptied on the first keystroke).
   - On a form built from chosen columns a column removed from the table keeps its field, empty
     (core#1214, fixed): claimed in `form.feature`.
   - A column removed from the table under a default form takes its field off (DEMOG removed: 10
     fields are shown, none for DEMOG, and no error is reported or logged).
2. **Fixture**: the reported case of core#1138 is a file holding `{}`; a dedicated
   `fixtures/form-not-a-form.json` would say it more plainly.
3. **Steps requested**: the arrows at the first and the last row flash the arrow (`ui.noMoreMatches`)
   and report nothing; the scenario claims only that the row stays. Every row filtered out leaves
   empty fields; claiming it needs a reading of a field's value (`value of <field>`), which the form
   does not report.
4. **API or Dart tests to write**: metadata changed ten times, then one keystroke in edit mode, writes
   the cell once (the field listeners are added again on every metadata change).
