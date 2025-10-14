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