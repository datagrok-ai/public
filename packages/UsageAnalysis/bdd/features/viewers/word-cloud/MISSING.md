# Word cloud — what the features cannot say yet

## Edge cases

The edge cases of the word cloud (26 Oct 2026) that no feature claims, and why. Nothing here asks
for a library, binding or code change.

### Still open

#### 1. Steps requested

None.

#### 2. Decided (the rules every viewer shares, listed in `viewers/MISSING.md`)

1. **The cloud over no rows**: an empty cloud, no message, no error: API test 3.
2. **The word column renamed** is followed: API test 4.
3. **Empty strings** are counted as a word with no visible text (`word-cloud-viewer.ts:217-221`):
   right as it is (2026-10-30), API test 5.
4. **Rotation Step 0 typed in the context panel** is taken, and the cloud still draws every word (public#4424, fixed; the render clamps the step to 1):
   claimed in "Rotation Step 0 typed in the context panel still draws every word (public#4424)", the last
   scenario of `word-cloud-appearance.feature`, under a rotation range of 0 to 90 so that the step is used.
5. **Inverted ranges** (Min Text Size above Max Text Size, Min Rotation above Max Rotation): swapped.
   Steps requested: the cloud reports the font of its first laid-out word only (`font`), not the
   size of each word, so no row claims which word is drawn largest.
6. **After the word column is removed** no other column is picked. API test 6 claims no error
   (core#1136); what the cloud shows with no word column, the requires-message or an empty cloud, is
   not said by the rules.

#### 3. API or Dart tests to write

1. Charts package test: word columns with exactly 500 and 501 distinct values (the gate is `> 500`,
   `:175`; the feature side has only the 1000-value case).
2. Charts package test: values over 50 characters (non-Latin and emoji too) are drawn cut with "…"
   and a click selects the rows of the full value.
3. Every row filtered out (AGE 500–600): 0 rows shown, no viewer error; the filter reset gives 1000.
4. RACE renamed to ORIGIN: the word column is ORIGIN, no viewer error, no error logged.
5. A word column with an empty value among three rows: 3 rows shown, 3 words, no error.
6. The word column (RACE) removed, then a filter (AGE 30–50, 494 rows): no error is logged
   (core#1136, fixed).

Moved out of BDD (no user gesture); requested from ApiTests.

### Won't do

Nothing.
