# Sankey and Chord edge cases — what the features cannot say yet

Written for `sankey-chord-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Steps requested**

1. The Charts package on the local stand reports no Sankey or Chord readings, so the feature claims the messages and a category label drawn after the change. With `nodes`, `links` and `categories` on the stand the claims can be counts.

**Questions for owners**

1. Sankey Source set through the property to a column with more than 50 categories (demog SUBJ): the 50-category rule is checked on the table only. Refuse the column?
2. Sankey Value column with negative values, or all zeros: refuse with a message as the pie does, or draw?
3. Sankey with a blank Source value: a node with an empty name, or skipped?
4. The Source / From column removed from the table: fall back to another column or show the requires-message? (The feature claims only "no error", tagged `@ticket-pending`: both viewers log "No such column".)

## Won't do

- Chord Include Nulls: help (`chord.md`) describes it as showing categories without connections; no edge beyond the existing feature.
