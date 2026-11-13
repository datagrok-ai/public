# Histogram — what the features cannot say yet

Written during the bar chart / box plot / histogram / pie / line chart / density plot gap round of
23 Oct 2026, against the histogram's "Automation surface" (`core/client/d4/lib/src/viewers/histogram/CLAUDE.md`)
and the TestTrack checklists (`files/TestTrack/Viewers/Histogram/`); trimmed in the BDD PR round of 8 Nov 2026.

## Resolved

- The distribution lines follow a filter made elsewhere (core#998): Dart `Histogram | Distribution lines follow the
  rows shown…`.

## Still open

- **X Axis Height (core#1086)**: with Auto Layout off the bars get shorter (claimed), but the drawn axis strip is not
  reported; needs an `x axis strip height` reading from the axis render in `histogram_core.dart`.
- **The Show Values counts follow the value column (GROK-19759)**: the `bin labels` strip is reported, not the text drawn
  in it; needs `bin label <n>` readings from the Show Values loop in `histogram_core.dart`.

## Won't do

- A maximum below the minimum on a date column: the rule and its message are claimed on AGE ("A maximum below the
  minimum collapses the filter and reports it"); a date typed into the range input is claimed on STARTED in the first
  scenario.
- Bins 201 typed in the panel (core#1200): Dart `A property keeps a value inside its declared bounds`.
- Every row filtered out under Normalize To Filter (core#998): Dart `Data Manipulation | Histogram`, ApiTests
  `AI: Viewers: Histogram extras` "renders on a table with no rows (core#998)".
- Min equal to max in the range inputs, the selectors at 200 × 200, a bound past the dates of a date column, a layout
  without its Value column: not worth a scenario, or the shared viewer rules.
