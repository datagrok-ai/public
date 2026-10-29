# Tree edge cases — what the features cannot say yet

Written for `tree-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data shapes,
columns removed or renamed under the open viewer, layouts naming a missing column, property bounds,
closing the viewer. Nothing here asks for a library, binding or code change made by that feature; it
lists what the feature could not claim and why.

## Still open

### 1. Steps requested

None.

### 2. Decided (the rules of `viewer-rules-edge-cases.feature`)

1. **Every row filtered out**: an empty tree, no message, no error; the branches come back. A row of
   the no-rows outline.
2. **The Size and Color columns removed** are dropped, with no other column picked. The removal
   handler cleans the hierarchy only (`:629-633`), and a removed Size or Color column throws on every
   render: public#4423, known failures in the removal outline.
3. **Values above the declared maximum typed in the panel** (Initial Tree Depth 6, Font Size 31) are
   not taken: known failures with no ticket yet in the out-of-bounds outline. A negative Label
   Rotate has no declared minimum (`:130`), so the rule does not reach it.
4. **Closing the tree in Filter mode** gives every row back (public#4425); a renamed hierarchy column
   with every row passing, a layout naming a missing hierarchy column and the removed column of the
   clicked branch are rows of the outlines too.

5. **A column with more than 500 categories** is refused whatever the filter. Today the check runs
   only while at least 500 rows pass (`tree-viewer.ts:921-927`), so under a filter of fewer rows the
   message turns into a tree with one branch per category: a known failure with no ticket yet, the
   category-limit outline of `viewer-rules-edge-cases.feature`.

### 3. Questions for owners

None.

### 4. API or Dart tests to write

None.

## Won't do

Nothing.
