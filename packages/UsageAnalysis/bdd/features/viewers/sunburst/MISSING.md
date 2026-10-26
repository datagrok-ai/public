# Sunburst edge cases — what the features cannot say yet

Written for `sunburst-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data shapes,
columns removed or renamed under the open viewer, layouts naming a missing column, property bounds,
closing the viewer. Nothing here asks for a library, binding or code change made by that feature; it
lists what the feature could not claim and why.

## Still open

### 1. Steps requested

None.

### 2. Questions for owners

1. **Filter mode, the clicked segment's column removed.** After a click on F (553 rows pass) and
   SEX removed, the table stays filtered to F while the sunburst shows RACE only, and nothing on it
   says a filter is active; a double-click on empty space clears it (claimed). Should the removal
   itself release the sunburst's filter? (`sunburst-viewer.ts:238-247` keeps `viewerFilter`.)
2. **A hierarchy column renamed**: the viewer has no rename handler and the level silently
   disappears; follow the rename or drop the level?
3. **One level.** The help says at least two levels are needed (sunburst.md:30); the viewer draws
   one. Which is right?
4. **The 500-category limit.** The message says "fewer than 500 unique categories"
   (`sunburst-viewer.ts:423`), the check keeps columns with at most 500 (`:412`), and it runs only
   while at least 500 rows pass the filter (`:407`), so the hierarchy changes as the user filters.
   Where is the boundary, and should it depend on the filter?
5. **More than 20 hierarchy columns** are cut to 20 silently (`:429`): a message, or a limit in the
   dialog?
6. **Closing the sunburst in Filter mode** leaves the table filtered to the clicked segment (pinned
   as a known failure, `@ticket-pending`): the expectation that closing gives every row back is
   taken from core#233, where the calendar keeping its filter after close is a defect. Confirm it is
   a defect for the sunburst too, and file it.
7. **A layout naming a missing hierarchy column** keeps the remaining levels (claimed: RACE); is
   that the intended fallback, rather than the default columns (public#2979)?

### 3. API or Dart tests to write

1. Charts package test: columns with 500 and 501 categories on a table of at least 500 rows, then a
   filter below 500 rows.
2. Charts package test: 21 hierarchy columns.

## Won't do

Nothing.
