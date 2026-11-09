# Row Source — what the features cannot say yet

## Won't do

- Closing the viewer whose group is hovered empties MouseOverGroup (core#1227): xamgle `d4 | Viewers | Fixed issues | Row
  Source MouseOverGroup empties when the group is released` and ddt `highlight clears with the mouse-over group`.
- Row Source and Filter through a saved layout for the seven viewers: what a layout stores is an API test, and the
  viewers' own persistence features restore their settings.
- A Filter formula that does not parse, FilteredSelected with a selection the filter hides, two viewers sharing one
  formula column: `row-source.feature` and the Dart viewer tests; a renamed Filter formula is ApiTests `AI: Viewers:
  Shared rules`.
