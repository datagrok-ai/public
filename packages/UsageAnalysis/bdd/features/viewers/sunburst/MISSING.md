# Sunburst — what the features cannot say yet

## Closed issues the regression round of 26 Oct 2026 could not claim

### Regressions found 26 Oct 2026

- **Rows selected in the grid are not drawn on the sunburst (public#2262, closed).** public#2262 asked
  that a selection made in the grid or in a viewer light up only the segments it belongs to. On
  localhost (core master; Charts 1.8.1, with a debug Charts build published there for the admin account) the sunburst draws nothing of a
  grid selection.

  Steps: open demog-1000; add a Sunburst (its default hierarchy is SEX, CONTROL, RACE); in the grid,
  select the rows where SEX is F (Select > … or a Shift-click on the F cells) — 553 rows selected.

  Measured (two `grok-bdd` runs and one in-page probe): the sunburst's canvas does not change at all —
  166 px in the selection colour before and 166 after (the probe's pixel histogram of the canvas is
  identical before and after the selection). Expected: the F segment (and the segments under it)
  drawn in the selection colour.

  Likely cause: `TreeUtils.toForestAsync` (`packages/Charts/src/utils/tree-utils.ts`) builds the
  tree with `linkSelection = false` since 3b60ca785f (#2954, Aug 2024), so the aggregated table's
  selection is never linked to the source table and `markSelectedNodes` finds no selected path. The
  Tree viewer, which links it, draws the same grid selection (`tree-regressions.feature`).

  The scenario, to go back into `sunburst-regressions.feature` once fixed:

  ```gherkin
  Scenario: Rows selected in the grid light up their segments (public#2262)
    Given user opens demog-1000 dataset
    And user adds a sunburst viewer
    Then sunburst viewer should be painted
    When user takes a snapshot of sunburst viewer
    And user selects rows where "SEX" is "F"
    Then 553 rows should be selected
    And sunburst viewer should show more selection highlight than before
    When user takes a snapshot of sunburst viewer
    And user clears the row selection
    Then no rows should be selected
    And sunburst viewer should show less selection highlight than before
    And no errors should have been logged
  ```
