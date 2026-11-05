# Peptides SAR features — what is still missing

## Peptide Generation

Feature: `peptide-generation.feature` (30 Oct 2026, local stand).

### Questions (Peptide Generation)

- **The default Activity column.** Added to the peptides dataset (ID, AlignedSequence, IC50), the viewer takes
  **ID** as its Activity (`onTableAttached`: the first numerical column whose name contains "activity", else
  the first numerical column), so it generates peptides "for highest activity" ranked by row ID. The SAR dialog
  proposes IC50 on the same table. The help (`visualize/viewers/peptide-generation.md`) does not say how the
  column is picked. Should the viewer skip integer ID-like columns, or take the column the SAR dialog proposes?
  The feature sets IC50 itself and claims nothing about the default.

### Not written (Peptide Generation)

- **Filtering**: the help says the statistics, predictions and peptides "are computed using only the filtered
  rows". No visible number in the viewer is bound to change with a filter (the count stays at Peptide count),
  and the predicted values are the package test's job; a claim needs a reading such as `rows used` in a
  `getWidgetStatus` the viewer does not publish yet (the other SAR viewers do).
- **The basis cell tooltip** (the activity distribution of the monomer at that position): the viewer reports no
  hit areas, so no step can hover a basis cell by name.
