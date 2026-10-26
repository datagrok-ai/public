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

### Questions for owners

3. **Combined boolean card, the kept flag's column removed.** Observed on localhost: Combined Boolean
   card, click CONTROL (6 rows pass), remove the CONTROL column. The card lists no flag, but the table
   stays filtered (not every row passes). A card on a removed column otherwise goes and releases its
   rows (`filter-panel-edge-cases.feature`). Is this a defect? It has no ticket; the scenario waits
   for the answer.
4. **What the categorical search does with `(`, `.`, `+`, `[`.** core#1134 is the error; whether the
   search should match the literal text or treat the word as a regular expression is not documented
   (only `,` and `*` are search syntax).
5. **A range card's min typed above its max**: which rows pass while the validator flags the fields.
6. **Hierarchical card when a level column is removed**, and a layout naming a removed level column:
   drop the level, drop the card, or fall back to auto-detection.
7. **An expression rule naming a removed column** (`${HEIGHT} < 150`, then HEIGHT removed): drop the
   rule or flag it.
8. **A saved filter state saved twice under one name** (overwritten silently today), and a state
   applied to a table of the same shape that lacks the saved category.
9. **The text card's extended-search operators** (`!malt`, `^hop`) at fuzziness above 0: operators or
   literal text (help filters.md says only "fuzzy search").
10. **Combined Boolean added to a table with no boolean column**: a message, an empty card, or a
    disabled menu item.
11. **Linked tables: the source view closed while its filter narrows the target**: release or keep.

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
