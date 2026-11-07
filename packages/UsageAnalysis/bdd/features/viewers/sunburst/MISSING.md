# Sunburst — what the features cannot say yet

## Edge cases

The edge cases of the Sunburst (26 Oct 2026) that no scenario claims, and why. The claimed
ones are `@edge-cases` scenarios of their own, after the master scenarios of `sunburst.feature`. Nothing here asks for a library, binding or code change.

### Still open

#### 1. Steps requested

None.

#### 2. Decided (the rules every viewer shares, listed in `viewers/MISSING.md`)

1. **Filter mode, the clicked segment's column removed**: the removal releases the sunburst's
   filter. Today the table stays filtered to F after SEX is removed until the next refilter, as for
   the bar chart and the tree: core#1174, a known failure: "Removing the column of a clicked segment
   gives back the rows it filtered out (core#1174)", the last scenario of `sunburst.feature`.
2. **A hierarchy column renamed** is followed (SEX, ORIGIN): API test 5.
3. **Closing the sunburst in Filter mode** gives every row back: public#4425, fixed, claimed in an `@edge-cases` scenario of its own.
4. **A layout naming a missing hierarchy column** is applied and the sunburst draws without the
   column: API test 7. No rows: API test 6.

5. **The 500-category limit** counts the rows that pass the filter: the check runs only while at
   least 500 rows pass (`CATEGORIES_NUMBER` in `sunburst-viewer.ts`): right as it is (2026-10-30), API tests 1 and 8.
6. **More than 20 hierarchy columns** are cut to 20, and the viewer says that part was dropped.
   Today the cut is silent (`MAXIMUM_COLUMN_NUMBER`). Steps requested: the note has no wording yet, and the
   `hierarchy columns` reading reports the property, not the levels drawn; wanted a `levels`
   reading, then: 21 categorical columns, claim 20 levels and the note.

#### 3. Questions for owners

1. **One level.** The help says at least two levels are needed (sunburst.md:30); the viewer draws
   one. Which is right?

#### 4. API or Dart tests to write

1. Charts package test: columns with 500 and 501 categories on a table of at least 500 rows, then a
   filter below 500 rows.
2. Charts package test: 21 hierarchy columns.
3. An empty hierarchy shows "The Sunburst viewer requires at least one categorical column with fewer
   than 500 unique categories"; SEX, RACE set back gives 10 segments without the message; no error.
4. SEX removed from the table: the hierarchy is RACE, 4 segments, 896 rows in Caucasian, no error.
5. RACE renamed to ORIGIN: the hierarchy is SEX, ORIGIN, no viewer error, no error logged.
6. Every row filtered out (AGE 500–600): 0 rows shown, no viewer error; the filter reset gives 1000.
7. A layout saved with SEX, RACE and loaded after SEX is removed: the sunburst is shown, the hierarchy
   setting stays SEX, RACE, no error.
8. Hierarchy USUBJID shows the 500-category message; filtered to SEX M (447 rows) the message goes;
   no error.

Moved out of BDD (no user gesture); requested from ApiTests.

### Won't do

Nothing.
