# Sankey and Chord — what the features cannot say yet

## Won't do

- The cycles message (Source and Target both SEX), the Chord's "Too many categories to render.", a Value column of
  negatives or zeros (public#4462, fixed), a removed Source or From column then a filter (public#4420, public#4421,
  fixed), a blank Source value as `(empty)`: Charts `Viewers` and `Viewers: no rows`, the rules every viewer shares
  (`../MISSING.md`); the messages are cosmetic.
- Chord Include Nulls: help (`chord.md`) describes it as showing categories without connections; no edge beyond the
  existing feature.
