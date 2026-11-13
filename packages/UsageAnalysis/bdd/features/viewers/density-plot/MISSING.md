# Density plot — what the features cannot say yet

Written during the viewer gap review of 23 Oct 2026 against `density_plot/CLAUDE.md` "Automation
surface"; trimmed in the BDD PR round of 8 Nov 2026.

## Resolved (core#1048)

- Bins 0, −5 or null and a 0 × 0 root: Dart `Density plot | Invalid bins keep the view…` and `A property keeps a
  value inside its declared bounds`.
- A logarithmic X axis bins only the positive values, and a click selects the rows a bin holds in either shape: Dart
  `Density plot | A click selects the rows a bin holds for every shape and axis…`.
- Listeners across table swaps: Dart `Table swap | Density plot` and `Lifecycle`.

## Still open

- **The bin selector on screen, and gone with Show Bin Selector off** (`density-plot-ui.md`): the viewer publishes it as
  the `bin selector` part; needs a library step `the {string} part of {widget} should be visible` / `… hidden`.

## Won't do

- A click on a rectangular bin (core#1048): Dart `Density plot | A click selects the rows a bin holds…`.
- Rectangular bins after a zoom (core#1048): Dart `Density plot | Rectangular bins over a narrowed view bin every row
  unless Bin To Range is on`.
- A column name with markup characters in the bin tooltip (core#1048): Dart `Data shows as text, not markup…`.
- The bin tooltip's range on a logarithmic axis: not worth a reading.
