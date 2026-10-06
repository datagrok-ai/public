---
feature: projects
realizes: [views.projects]
target_layer: manual-only
coverage_type: regression
companion_to: projects-copy-clone.md
manual_only_reason: |
  Thumbnail render quality and tile-layout consistency, and the preserved
  workspace layout on reopen, are visual judgments without an explicit
  per-element criterion.
---

# Projects — tile rendering and view-state preservation (manual)

Two visual checks around project tiles and Save Copy: how project tiles render
in Dashboards, and whether a copy saved with Personal View Customizations
brings the customized view back on reopen.

Self-contained: the scenario creates its own project and deletes it in Cleanup.

## Create the test project

1. Close all, open `System:DemoFiles/demog.csv`
2. Add a Scatter plot and a Histogram to the view
3. Select **File > Save Project**, name it `test_copy_clone`, click OK — the
   project is saved; cancel the Share dialog if it opens
4. Go to **Browse > Dashboards** — the `test_copy_clone` tile is listed

## Tile rendering in Dashboards

1. Locate the `test_copy_clone` tile among the other project tiles
2. Its thumbnail renders as a visible preview image — not a blank box, not a
   broken-image icon, and the same size as neighbouring tiles
3. The project name below the thumbnail is fully readable at the default panel
   width — no clipping or truncation artifacts
4. The tile's metadata (creation date, owner) renders in line with the
   neighbouring tiles, nothing shifted or overlapping
5. Hover the tiles and scroll the gallery — all tiles keep the same dimensions,
   nothing overflows, is cropped, or jitters

## Personal View Customizations preserved on reopen

1. Open `test_copy_clone` from Dashboards
2. Customize the view:
   - filter `AGE > 50` in the Filter Panel
   - sort the grid by HEIGHT descending
   - hide the DIS_POP column
   - drag the Histogram to a different dock position
3. Select **File > Save Copy**, choose the **Personal View Customizations**
   mode, name the copy `test_copy_clone_pvc`, click OK — the copy is saved and
   appears in Dashboards
4. Close all
5. Reopen `test_copy_clone_pvc` from Dashboards
6. The customized view comes back: the `AGE > 50` filter is applied (the
   filtered row count matches what step 2 produced), the grid is still sorted
   by HEIGHT descending, DIS_POP is still hidden, and the Histogram sits in the
   dock position it was moved to
7. The rest of the workspace matches the pre-save state — panel positions and
   viewer sizes, with no drift

## Cleanup

Delete both projects created by this run — `test_copy_clone` and
`test_copy_clone_pvc` (Browse > Dashboards, right-click the tile > **Delete**).
After a partial run, delete whichever of them exists.

---
{
  "order": 3,
  "datasets": ["System:DemoFiles/demog.csv"]
}
