# Form — what the features cannot say yet

Written 2026-10-22 when `form.feature` took the manual sections "Column selector", "Column changes
reaction" and "Table switching" of TestTrack `Viewers/Form/form.md`; trimmed in the BDD PR round of 8 Nov 2026.

## Still open

- **The text a field's label shows**: the rename of AGE to AGE_NEW is claimed through the `label AGE_NEW` area, keyed by
  the column; needs a `label text of <column>` reading in the form's status.

## Won't do

- The arrows at the first and the last row (a no-op), a CSV loaded as a form (the layout-file case is claimed), every
  row filtered out, ten metadata changes followed by one keystroke: low value or not UI checks.
- A column removed under a default form: the shared viewer rule (Dart `Column removal`, `Data Manipulation`).
