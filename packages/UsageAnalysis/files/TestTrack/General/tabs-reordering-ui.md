---
feature: general
target_layer: manual-only
coverage_type: regression
produced_from: split
original_path: public/packages/UsageAnalysis/files/TestTrack/General/tabs-reordering.md
split_date: 2026-06-16
related_bugs: []
manual_only_reason: |
  The core of this scenario is drag-and-drop reordering of table tabs, which is
  implemented by the custom dock_spawn docking layer (core/client/libs/dock_spawn).
  Playwright dragTo / mouse-drag against dock_spawn tab headers is unreliable
  (no stable DOM drop target; the reorder is driven by pixel-level pointer-move
  events on a canvas-backed tab strip). Step 3 (order persistence after
  project save/reopen) cannot be automated independently either: verifying that
  a *modified* order survives a reload requires first establishing that order
  via the manual drag in steps 1-2. Distinct from Viewers/grid-ui.md, which
  covers column reordering on the grid canvas — this is tab reordering.
---

### Tabs Ordering and Persistence

1. Open multiple datasets.
- Drag and drop two different datasets into the workspace.
- Open two additional random datasets from the demo files.
- Expected result: All four datasets should open without errors.
2. Rearrange datasets by drag-and-drop.
- Hover over a dataset name, press and hold the left mouse button, then drag datasets left and right to reorder them randomly.
- Verify the following:
  - All datasets should be movable without errors.
  - Dropping a dataset in a new position should correctly update the order.
  - A dataset can be returned to its original position without issues.
  - Adding a new dataset should not affect the modified order of existing datasets.
3. Verify order persistence after project reload
- Memorize the modified dataset order.
- Save the project.
- Files > Close all.
- Reopen the project.
- Expected result: The dataset order should remain exactly as it was before saving.
4. Additional Notes:
- Ensure that dataset movement is smooth and does not cause UI glitches.
- Verify that datasets snap into place correctly when repositioned.
- Check that closing and reopening datasets does not reset the order.

---
{
  "order": 8
}