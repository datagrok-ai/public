---
feature: powerpack
realizes_atlas: [powerpack.cp.formula-lines-dialog]
priority: p0
target_layer: manual-only
coverage_type: regression
related_bugs:
  - github-2487
  - github-671
  - github-2747
manual_only_reason: |
  All checks are visual judgments on canvas-rendered formula lines and the
  live preview inside the Formula Lines editor — no DOM state to assert on.
---

# Formula Lines — regression checks

Each section guards a fixed GitHub issue. All scenarios start with:
1. Close all
2. Open SPGI

## Datetime columns on scatter plot ([#2487](https://github.com/datagrok-ai/public/issues/2487))

1. Add a Scatter plot; set X to Synthesis Date (datetime), Y to Average Mass
2. Right-click the canvas and select Tools > Formula Lines
3. Click ADD NEW and add a vertical line bound to the datetime X axis
4. Verify the editor accepts the datetime column (no error, preview shows the line)
5. Click OK — the line renders on the scatter plot at the expected date position
6. Cleanup: reopen the editor and delete the created line

## Editor preview follows viewer configuration ([#671](https://github.com/datagrok-ai/public/issues/671))

1. Add a Scatter plot with numeric X and Y
2. Open Tools > Formula Lines and add a horizontal line — the preview shows it
3. Close the editor; change the scatter plot's X and Y columns
4. Reopen the Formula Lines editor — the preview reflects the new axes, and the
   existing line is drawn against them (not against the stale configuration)
5. In the editor, select each line in the list — the preview updates to
   highlight the selected line
6. Cleanup: delete the created line

## Consistency between scatter plot and line chart ([#2747](https://github.com/datagrok-ai/public/issues/2747))

1. Add a Scatter plot and a Line chart on the same columns (same X and Y)
2. Open Tools > Formula Lines on the scatter plot; on the Data frame tab add a
   horizontal line at a fixed Y value
3. Verify the line renders on both viewers at the same position
4. On each viewer, open the Formula Lines editor and verify the axes shown in
   the preview match that viewer's actual axes and can be changed from the dialog
5. Cleanup: delete the created line

---
{
  "order": 1,
  "datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
