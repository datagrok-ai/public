---
feature: linechart
sub_features_covered: [viewers.line-chart]
target_layer: manual-only
pyramid_layer: regression
ui_companion_of: legend-color-and-persistence.md
companion_steps: [9, 10, 11, 12, 13]
manual_execution_notes: |
  Project save/reopen round-trip for per-category color persistence (GROK-17278)
  and markers-legend presence on project reopen (GROK-19825). Split to manual
  because grok.dapi.projects.save(grok.shell.project) currently fails on dev with
  "Unable to add entity to the project" for a file-backed TableView (a dev-server
  regression, verified 2026-07-18 via chrome-devtools MCP), so the project
  round-trip cannot be driven headless; only the SAVE ribbon dialog completes it.
  The layout round-trip for the same color signal is covered automatically in
  legend-color-and-persistence-spec.ts.
---

# Line Chart — project save/reopen persistence (manual-only companion)

Companion to `legend-color-and-persistence.md`. Covers Scenario 2 Steps 9-13
(the project round-trip), which the automated spec cannot drive because the
JS-API project save is a current dev-server regression.

## Manual steps

1. With the split Line Chart from the parent Setup active (X = `Chemical Space X`,
   Y = `Chemical Space Y`, split = `Stereo Category`), change one category's color
   (right-click the legend entry, Change Color, or the property-panel color picker).
   Note the chosen color.
2. Click the **SAVE** ribbon button, choose **Save Project**, give it a distinct
   name (e.g. `linechart-color-persist-test`), and save.
3. Close all open tables and views.
4. Reopen the project from the **Projects** panel (Browse tree).
5. Locate the Line Chart viewer in the reopened project layout.
6. Verify: the per-category color reads back the same value chosen in Step 1
   (GROK-17278 — project serialization preserves per-category color assignments).
7. Verify: the markers legend is present in the reopened viewer (GROK-19825 —
   the legend must not disappear on project reopen).
8. Verify: no console or page errors occurred during Steps 2-7.
9. Delete the saved project (cleanup).

## Notes

- The color property is stored on the split column's metadata
  (`col('Stereo Category').meta.colors`, categorical). Read back a category's
  color in the console via
  `grok.shell.tv.dataFrame.col('Stereo Category').meta.colors.getColor(rowIdx, col)`.
- The markers legend is the DOM node `[name="viewer-Line-chart"] [name="legend"]`.
