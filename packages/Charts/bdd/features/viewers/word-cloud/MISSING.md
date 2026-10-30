# Word cloud edge cases — what the features cannot say yet

Written for `word-cloud-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data
shapes, columns removed or renamed under the open viewer, layouts naming a missing column, property
bounds, closing the viewer. Nothing here asks for a library, binding or code change made by that
feature; it lists what the feature could not claim and why.

## Still open

### 1. Steps requested

None.

### 2. Decided (the rules of `viewer-rules-edge-cases.feature`)

1. **The cloud over no rows**: an empty cloud, no message, no error. A row of the no-rows outline.
2. **The word column renamed** is followed: a row of the rename outline.
3. **Empty strings** are counted as a word with no visible text (`word-cloud-viewer.ts:217-221`):
   right as it is (2026-10-30), claimed in `viewer-rules-edge-cases.feature`.
4. **Rotation Step 0 typed in the context panel** is not taken: public#4424, a known failure in the
   out-of-bounds outline.
5. **Inverted ranges** (Min Text Size above Max Text Size, Min Rotation above Max Rotation): swapped.
   Steps requested: the cloud reports the font of its first laid-out word only (`font`), not the
   size of each word, so no row claims which word is drawn largest.
6. **After the word column is removed** no other column is picked. The scenario stays in
   `word-cloud-edge-cases.feature` (core#1136, claims no error); what the cloud shows with no word
   column, the requires-message or an empty cloud, is not said by the rules.

### 3. API or Dart tests to write

1. Charts package test: word columns with exactly 500 and 501 distinct values (the gate is `> 500`,
   `:175`; the feature side has only the 1000-value case).
2. Charts package test: values over 50 characters (non-Latin and emoji too) are drawn cut with "…"
   and a click selects the rows of the full value.

## Won't do

Nothing.
