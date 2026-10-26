# Word cloud edge cases — what the features cannot say yet

Written for `word-cloud-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data
shapes, columns removed or renamed under the open viewer, layouts naming a missing column, property
bounds, closing the viewer. Nothing here asks for a library, binding or code change made by that
feature; it lists what the feature could not claim and why.

## Still open

### 1. Steps requested

None.

### 2. Questions for owners

1. **The cloud over no rows**: an empty cloud or a message? (No error and the words coming back are
   claimed.)
2. **After the word column is removed** (core#1136, pinned): pick another string column, or show a
   message?
3. **The word column renamed**: follow the rename?
4. **Empty strings** are counted as a word with no visible text (`word-cloud-viewer.ts:217-221`):
   skip them, or draw "(empty)" as the Tree does?
5. **Inverted ranges**: Min Text Size above Max Text Size, Min Rotation above Max Rotation (no
   bounds, `:76-80`): swap, reject, or accept?

### 3. API or Dart tests to write

1. Charts package test: word columns with exactly 500 and 501 distinct values (the gate is `> 500`,
   `:175`; the feature side has only the 1000-value case).
2. Charts package test: values over 50 characters (non-Latin and emoji too) are drawn cut with "…"
   and a click selects the rows of the full value.

## Won't do

Nothing.
