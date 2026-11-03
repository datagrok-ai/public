# Form — what the features cannot say yet

Written 2026-10-22 when `form.feature` took the manual sections "Column selector", "Column changes
reaction" and "Table switching" of TestTrack `Viewers/Form/form.md`. No library, binding or code
change was made.

## Still open

1. **The text a field's label shows.** The rename of AGE to AGE_NEW is claimed through the
   `label AGE_NEW` area, which is keyed by the column the label belongs to, not by what it shows.
   Wanted: a `label text of <column>` reading in the form's status, then the feature claims the
   label reads "AGE_NEW".
