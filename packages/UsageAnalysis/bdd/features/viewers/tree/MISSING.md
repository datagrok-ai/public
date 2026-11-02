# Tree — what the features cannot say yet

## Edge cases

The edge cases of the Tree (26 Oct 2026) that no scenario claims, and why. The claimed ones are
`@edge-cases` scenarios of their own, after the master scenarios of `tree-on-click-filter.feature` and `tree.feature`. Nothing here asks for a library,
binding or code change.

### Still open

#### 1. Steps requested

None.

#### 2. Decided (the rules every viewer shares, listed in `viewers/MISSING.md`)

1. **Every row filtered out**: an empty tree, no message, no error; the branches come back: API test 4.
2. **The Size and Color columns removed** are dropped, with no other column picked. The removal
   handler cleans the hierarchy only (`:629-633`), and a removed Size or Color column throws on every
   render: public#4423, fixed, API test 5.
3. **Values above the declared maximum typed in the panel** (Initial Tree Depth 6, Font Size 31) are
   not taken: core#1200, known failures in a scenario outline of `tree.feature`, one Examples
   row per property. A negative Label
   Rotate has no declared minimum (`:130`), so the rule does not reach it.
4. **Closing the tree in Filter mode** gives every row back (public#4425); the removed column of the
   clicked branch should too (core#1174, a known failure). Both are claimed in
   `tree-on-click-filter.feature`. A renamed hierarchy column and a layout naming a
   missing hierarchy column: API tests 1 and 6.

5. **A column with more than 500 categories** is refused only while at least 500 rows pass
   (`tree-viewer.ts:921-927`); under a filter of fewer rows the tree is drawn: right as it is
   (2026-10-30): API test 8.

#### 3. Questions for owners

None.

#### 4. API or Dart tests to write

1. RACE renamed to ORIGIN, with every row passing and under a filter of SEX M (447 rows): the
   hierarchy is CONTROL, SEX, ORIGIN, the tree keeps drawing, no error.
2. RACE removed from the table: the hierarchy is CONTROL, SEX, the RACE level is gone, no error.
3. An empty value is the `(empty)` branch only while Include Nulls is on (the sunburst claims its
   twin in `sunburst.feature`, Empty values).
4. Every row filtered out (AGE 500–600): 0 rows shown, no viewer error; the filter reset gives 1000.
5. The Size (HEIGHT) or Color (AGE) column removed: the property no longer names it, no viewer error,
   no error logged (public#4423, fixed).
6. A layout saved with CONTROL, SEX, RACE and loaded after RACE is removed: the tree is shown, the
   hierarchy setting stays CONTROL, SEX, RACE, no error.
7. Initial Tree Depth 6 and Font Size 31 set through the API are not taken (3 and 12 stay): known
   failure core#1200, still open.
8. Hierarchy USUBJID shows the 500-category message; filtered to SEX M (447 rows) the message goes;
   no error.

Moved out of BDD (no user gesture); requested from ApiTests.

### Won't do

Nothing.
