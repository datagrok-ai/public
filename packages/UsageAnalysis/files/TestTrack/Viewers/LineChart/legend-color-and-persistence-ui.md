---
feature: linechart
sub_features_covered: [viewers.line-chart]
target_layer: manual-only
pyramid_layer: regression
ui_companion_of: legend-color-and-persistence.md
companion_steps: [13]
manual_execution_notes: |
  Markers-legend PRESENCE on project reopen (GROK-19825) only. Manual because a
  df-only project restores the TABLE but not the viewer LAYOUT headless, so the
  legend DOM node is not reconstructed — the visual "is the legend there after
  reopen" check needs a human eye on a real reopened layout.

  CORRECTION 2026-07-21: the earlier note ("grok.dapi.projects.save fails on dev
  with a dev-server regression") was a MIS-DIAGNOSIS. The JS-API project save
  fails with project_relations_entity_id_fkey ONLY when the table has no
  server-side entity; calling grok.dapi.tables.uploadDataFrame(df) first makes
  save/find/reopen/delete work headless (verified 2026-07-21). Consequently the
  per-category color persistence (GROK-17278) is now covered AUTOMATICALLY via a
  full project save/close/reopen in legend-color-and-persistence-spec.ts
  (S2 Steps 9-12) — no longer manual. Only the viewer-layout-dependent legend
  presence (Step 13) remains here.
---

# Line Chart — markers legend on project reopen (manual-only companion)

Companion to `legend-color-and-persistence.md`. Covers Scenario 2 Step 13 — the
visual "markers legend is present after a project reopen" check (GROK-19825).
The per-category color persistence (GROK-17278, Steps 9-12) is now AUTOMATED in
`legend-color-and-persistence-spec.ts` via a full project save/close/reopen; only
this legend-presence check stays manual, because a headless df-only project
reopen restores the table but not the viewer layout (no legend DOM node to read).

## Manual steps

1. With the split Line Chart from the parent Setup active (X = `Chemical Space X`,
   Y = `Chemical Space Y`, split = `Stereo Category`) and a per-category color set,
   click the **SAVE** ribbon button, choose **Save Project**, give it a distinct
   name (e.g. `linechart-color-persist-test`), and save.
2. Close all open tables and views.
3. Reopen the project from the **Projects** panel (Browse tree) and locate the
   Line Chart viewer in the restored layout.
4. Verify: the markers legend is present in the reopened viewer (GROK-19825 —
   the legend must not disappear on project reopen).
5. Verify: no console or page errors occurred during Steps 1-4.
6. Delete the saved project (cleanup).

## Notes

- The color property is stored on the split column's metadata
  (`col('Stereo Category').meta.colors`, categorical). Read back a category's
  color in the console via
  `grok.shell.tv.dataFrame.col('Stereo Category').meta.colors.getColor(rowIdx, col)`.
- The markers legend is the DOM node `[name="viewer-Line-chart"] [name="legend"]`.
