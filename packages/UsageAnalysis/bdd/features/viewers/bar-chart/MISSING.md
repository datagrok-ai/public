# Bar chart — what the features cannot say yet

Written during the bar chart / box plot / histogram / pie / line chart / density plot gap round of
23 Oct 2026, against `files/TestTrack/Viewers/BarChart/bar-chart-ui.md`; trimmed in the BDD PR round of 8 Nov 2026.

## Still open

- **Scrolling with the range slider** (`bar-chart-ui.md` "Scrolling (range sliders)"): the vertical slider of zoomed
  categories (`vertScroll` in `bar_chart_core.dart`) is not in the status; needs `scroll slider`, `scroll min handle`
  and `scroll max handle` areas, as the histogram reports `range slider`.
- **The mouse wheel over zoomed categories** scrolls a 600-pixel step by about 11 pixels (`_initPanning` passes
  `deltaY / barHeight` to `scrollByPixels`): needs a ruling before a scenario.
- **public#463, a column picked in the on-chart value selector under count**: the step clicks the selector as a whole
  and opens the aggregation list; it needs to open the column half (`.d4-column-selector-column`) when there is one.

## Won't do

- Out-of-bounds values typed into a property (core#1200): the property rule — prop_gen `Property.clamp…`, Dart `A
  property keeps a value inside its declared bounds`. Text that does not parse: the generated setter puts the
  property's default back (`v ?? <default>`); not a UI claim.
- A Split column removed under a Filter click (core#1174, fixed): Dart `Bar chart click filter goes when the column it
  was made on is removed (core#1174)`.
- A date-stacked Filter click released by a click on empty space after another setting (core#974): Dart `Bar chart
  filter click on a bar stacked by a date column outlives another setting…`.
- A category with markup characters in the bar tooltip (core#974): Dart `Data shows as text, not markup…`.
- Stack Maps the editor cannot make, 200-character names, frame budgets, corner radius −5, every row filtered with Show
  Filtered Rows, zoom under a filter, empty-bar sorting, the Stack column removed: not UI checks (Dart `Data
  Manipulation`, `Non-finite values`, `Column removal`).
