# Sankey and Chord edge cases — what the features cannot say yet

Written for `sankey-chord-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Steps requested**

1. The Charts package on the local stand reports no Sankey or Chord readings, so the feature claims the messages and a category label drawn after the change. With `nodes`, `links` and `categories` on the stand the claims can be counts.

**Decided** (the rules of `viewer-rules-edge-cases.feature`)

1. A blank Source value is a node named `(empty)`. No outline row: the Sankey needs a numeric Value column the outline's three-column table lacks, and no reading lists its node names yet.
2. The Source / From column removed: no other column is picked and the viewer shows its requires-message. The two scenarios (public#4420, public#4421, fixed) claim only "no error"; they should also claim the message (`"The Sankey viewer requires…" text in sankey viewer` visible).

3. Sankey Source set through the property to a column with more than 50 categories (demog USUBJID) is cut as on the table, and the viewer says so (the silent-cut rule). Today the 50-category rule is checked only when the viewer picks its columns. Steps requested: the note has no wording yet, and the Sankey reports no readings on the stand (Steps requested, 1).

4. A Sankey Value column with negative values, or with only zeros, is refused with a message as the pie chart refuses them (`sum(NEG) contains negative values`, `sum(ZERO): all values are 0`). In the feature as two known failures: today the Sankey draws its node labels piled at the top with no links and no message. The Sankey's own wording is not fixed, so the claims read only the words the pie's messages use (`contains negative values`, `all values are 0`); with the wording fixed they can read the whole message.

## Won't do

- Chord Include Nulls: help (`chord.md`) describes it as showing categories without connections; no edge beyond the existing feature.
