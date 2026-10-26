# Scatter plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `scatterplot/CLAUDE.md` "Automation
surface", the TestTrack `scatterplot-*.md` and `scatter-plot-ui.md`. The "(no value)" legend click is
in `scatter-plot-legend-no-value.feature`; the rejected line formula (GROK-17282) in
PowerPack's `formula-lines-dialog.feature`.

## Needs a reading

- **Structures on the X axis** (`scatter-plot-ui.md`, "Structures as markers and as an axis", items
  7-8, spgi-100). With a molecule column on X the axis draws structures instead of tick labels;
  the status reports the axis column and range, not what the ticks are drawn as. Needs an
  `x axis labels drawn as` reading (`text` / `structure`, from the axis painter's renderer), claimed
  after the pick and again after a saved layout is applied.
