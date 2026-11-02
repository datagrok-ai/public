# Sankey and Chord — what the features cannot say yet

## Edge cases

The edge cases of the Sankey and Chord viewers that no feature states, with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

### Still open

**Steps requested**

1. The Charts package on the local stand reports no Sankey or Chord readings. With `nodes`, `links` and `categories` on the stand a check can claim counts, not only the messages and a category label.

**Decided** (the rules every viewer shares, listed in `viewers/MISSING.md`)

1. A blank Source value is a node named `(empty)`. No outline row: the Sankey needs a numeric Value column the outline's three-column table lacks, and no reading lists its node names yet.
2. The Source / From column removed: no other column is picked and the viewer shows its requires-message. The removal check (API test 3) should claim the message (`The Sankey viewer requires…`) as well as no error.

3. Sankey Source set through the property to a column with more than 50 categories (demog USUBJID) is cut as on the table, and the viewer says so (the silent-cut rule). Today the 50-category rule is checked only when the viewer picks its columns. Steps requested: the note has no wording yet, and the Sankey reports no readings on the stand (Steps requested, 1).

4. A Sankey Value column with negative values, or with only zeros, is refused with a message as the pie chart refuses them (`sum(NEG) contains negative values`, `sum(ZERO): all values are 0`). public#4462, still open: today the Sankey draws its node labels piled at the top with no links and no message. The Sankey's own wording is not fixed, so a check can read only the words the pie's messages use (`contains negative values`, `all values are 0`); with the wording fixed it can read the whole message.

**API or Dart tests to write**

1. A Sankey with Source and Target both SEX shows "The graph contains cycles. Please remove circular dependencies." and logs no error; Target RACE draws it again without the message.
2. A Chord with From SEX and To USUBJID (more than 1500 pairs) shows "Too many categories to render." and logs no error; To RACE draws it again without the message.
3. The Sankey's Source or the Chord's From column (SEX) removed from the table, then a filter (AGE 30–50, 494 rows): no error is logged (public#4420, public#4421, fixed).
4. A Sankey Value column of negative values (`-${AGE}`) or of zeros: the viewer refuses to draw and says "contains negative values" / "all values are 0"; Value AGE draws it again. Known failure public#4462, still open.

Moved out of BDD (no user gesture); requested from ApiTests.

### Won't do

- Chord Include Nulls: help (`chord.md`) describes it as showing categories without connections; no edge beyond the existing feature.
