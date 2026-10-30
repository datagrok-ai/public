# Sunburst edge cases — what the features cannot say yet

Written for `sunburst-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data shapes,
columns removed or renamed under the open viewer, layouts naming a missing column, property bounds,
closing the viewer. Nothing here asks for a library, binding or code change made by that feature; it
lists what the feature could not claim and why.

## Still open

### 1. Steps requested

None.

### 2. Decided (the rules of `viewer-rules-edge-cases.feature`)

1. **Filter mode, the clicked segment's column removed**: the removal releases the sunburst's
   filter. Today the table stays filtered to F after SEX is removed until the next refilter, as for
   the bar chart and the tree: core#1174, a known failure.
2. **A hierarchy column renamed** is followed: claimed (SEX, ORIGIN).
3. **Closing the sunburst in Filter mode** gives every row back: public#4425, fixed, claimed.
4. **A layout naming a missing hierarchy column** keeps the remaining levels (RACE), claimed; no rows
   is a row of the no-rows outline.

5. **The 500-category limit** does not depend on the filter: a column with more than 500
   categories is refused whatever the filter. Today the check runs only while at least 500 rows pass
   (`sunburst-viewer.ts:407`): a known failure with no ticket yet, the category-limit outline of
   `viewer-rules-edge-cases.feature`.
6. **More than 20 hierarchy columns** are cut to 20, and the viewer says that part was dropped.
   Today the cut is silent (`:429`). Steps requested: the note has no wording yet, and the
   `hierarchy columns` reading reports the property, not the levels drawn; wanted a `levels`
   reading, then: 21 categorical columns, claim 20 levels and the note.

### 3. Questions for owners

1. **One level.** The help says at least two levels are needed (sunburst.md:30); the viewer draws
   one. Which is right?

### 4. API or Dart tests to write

1. Charts package test: columns with 500 and 501 categories on a table of at least 500 rows, then a
   filter below 500 rows.
2. Charts package test: 21 hierarchy columns.

## Won't do

Nothing.
