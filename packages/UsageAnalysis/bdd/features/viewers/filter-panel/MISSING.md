# Filter panel edge cases — what the features cannot say yet

Added with `filter-panel-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

## Still open

### Steps requested

1. **The validation state of a range card's min / max fields.** Min equal to max, or above it, makes
   the card's input validator say "Max should be greater than min" (`histogram_core.dart`), and `abc`
   makes it say nothing. No step reads a card field's invalid state or its message. Wanted: the card's
   min/max inputs reporting `aria-invalid` and the message, then `"max" field of the "AGE" filter card
   should be invalid`.
2. **The free-text input's invalid state.** A free-text rule that does not parse is refused on Enter
   and the input gets `d4-invalid` with an error tooltip (`text_filter.dart`). The scenario claims only
   that no rule is added and the rows stay. Wanted: the state of the Expression card's search input.

### Waiting for core#1174

3. **Expression and Combined Boolean cards, a column they use removed** while the panel has no card of
   that column (with the column's own card present, its removal refilters the table and the rows come
   back): the table stays filtered. Three known failures in `filter-panel-edge-cases`.

### Decided, scenarios to write (2026-10-28)

4. **The categorical search with `(`, `.`, `+`, `[`** finds the literal text (core#1134, fixed, no error any more).
5. **A range card's min typed above its max**: the present handling is right (the fields are flagged
   with a message); claimed once a step reads the flag (item 1).
6. **Hierarchical card when a level column is removed**, and a layout naming a removed level column:
   the level is dropped.
7. **An expression rule naming a removed column**: covered by the rule "a filter goes with its
   column".
8. **A saved filter state saved twice under one name** asks before it overwrites; a state that cannot
   be applied to the table is not offered.
9. **Combined Boolean on a table with no boolean column**: the menu item is disabled.

10. **The text card's `!malt` and `^hop`**: right as it is. At Fuzzyness 0 the term is literal text
    (a plain "contains"); above 0 the card searches with Fuse.js extended search, where `!` is "does not
    contain" and `^` is "starts with" (`text_filter.dart`).
11. **Linked tables, filter to filter: the source view closed** while its filter narrows the target:
    the target's filter is released.

### API or Dart tests to write

12. A card on a numeric column whose values are all empty (Dart regTest with log capture).
13. A number typed with a decimal comma under a comma-decimal locale (`NumericMatcher`, Dart unit test).

## Won't do

- **A layout naming a column the table no longer has** (`createFilter(null)`): not reproduced on
  localhost (the other cards come back, no error), and a removed calculated column is recreated by the
  layout itself.
- **Card summary at 19 / 20 categories**: a cosmetic threshold of the summary text.
- **Already covered**: all rows filtered out (`filter-panel-categorical.feature`), Filter out missing
  values (`filter-panel-cloned-view.feature`), a column with one value (`filter-panel-categorical.feature`).
