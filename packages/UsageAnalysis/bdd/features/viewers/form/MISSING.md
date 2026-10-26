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

1. **Fixture.** The two `@ticket-pending` Load form from file scenarios use files the project already
   keeps (`fixtures/cars-small.csv`, `fixtures/nx/spgi-legend-position-old.layout`). The reported case
   is a file holding `{}`; a dedicated `fixtures/form-not-a-form.json` would say it more plainly (adding
   a fixture was outside this round).
2. **Questions for owners**:
   - The message for a file that is not a form, and whether the form keeps its fields.
   - Editable form: `abc` or `-` typed into the numeric AGE field writes an empty value into the cell
     on every keystroke. Refuse, or clear the cell?
   - Editable form: typing into the field of a removed column (the field stays on a designed form).
   - Every row filtered out: empty fields, or the last row's values.
3. **Steps requested**: the arrows at the first and the last row flash the arrow (`ui.noMoreMatches`)
   and report nothing; the scenario claims only that the row stays.
4. **Dart test to write**: metadata changed ten times, then one keystroke in edit mode, writes the cell
   once (the field listeners are added again on every metadata change).
