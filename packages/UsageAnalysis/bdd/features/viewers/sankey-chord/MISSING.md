# Sankey and Chord — what the features cannot say yet

## Click selection (viewer tails round, 2026-11-05)

Claimed in `sankey-chord.feature`: a click on a Sankey link selects the rows of that pair, a
Ctrl+click on a node adds the rows of that category, a plain click on a node replaces the
selection; a click on a Chord category selects its rows, Ctrl+click adds another category, a plain
click replaces the selection.

### Found on the run

- **A Sankey node click selects rows the filter hides** ([public#4514](https://github.com/datagrok-ai/public/issues/4514)). Under the Filter Panel's SEX = F (553 rows),
  a Ctrl+click on `node Asian` selects all 15 Asian rows of demog-1000, the 8 male ones included,
  not the 7 that pass the filter (measured on localhost, 2026-11-05). The Sankey's click predicate
  (`rowMatchesColumnValues`, `Charts/src/viewers/sankey/sankey.ts`) does not read the filter; the
  Chord's (`chord-viewer.ts`) starts with `this.filter.get(i)`. The scenario claims only that Asian
  rows are selected and Other rows are not; once fixed, claim `487 rows should be selected` and
  `every selected row should pass the filter` after the Ctrl+click, and `7 rows should be selected`
  after the plain click.
