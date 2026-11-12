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

Two visual checks: how project tiles render in **Dashboards**, and
whether personal view customizations bring the customized view back on
reopen, down to the viewer layout.

Self-contained: the scenario creates its own project and deletes it in
Cleanup.

## Create the test project

1. Right-click the left sidebar and select **Close All**.
2. Go to **Browse > Files > Demo** and double-click `demog.csv`.
3. In **Toolbox > Viewers**, click **Scatter plot**, then **Histogram**.
4. Click **SAVE** on the ribbon, enter `testCopyClone` as the name,
   click **OK**. In the **Share** dialog, click **CANCEL**.
5. Go to **Browse > Dashboards** and type `testCopyClone` into the
   search box — the `testCopyClone` tile is listed.

## Tile rendering in Dashboards

1. Locate the `testCopyClone` tile among the other project tiles.
2. Its thumbnail renders as a visible preview image — not a blank box,
   not a broken-image icon, and the same size as neighbouring tiles.
3. The project name below the thumbnail is fully readable at the default
   panel width — no clipping or truncation artifacts.
4. The date under the name renders in line with the neighbouring tiles,
   nothing shifted or overlapping.
5. Hover the tiles and scroll the gallery — all tiles keep the same
   dimensions, nothing overflows, is cropped, or jitters.

## Personal view customizations preserved on reopen

1. Right-click the left sidebar and select **Close All**. In **Browse >
   Dashboards**, double-click the `testCopyClone` tile.
2. Customize the view:
   - click the filter icon on the ribbon and, in the `SEX` filter, clear
     the check box next to `F`; note the **Filtered** count in the
     status bar;
   - right-click the `HEIGHT` column header and choose **Sort >
     Descending**;
   - right-click any column header, choose **Order or Hide Columns...**,
     uncheck `DIS_POP` and click **CLOSE**;
   - drag the histogram by its header to a different dock position.
3. Click **SAVE** on the ribbon and select **Save personal view
   customizations**. The **Name** field is greyed out and reads
   `testCopyClone`. Click **OK**.
4. Right-click the left sidebar and select **Close All**.
5. In **Browse > Dashboards**, double-click the `testCopyClone` tile —
   no new tile appeared, and the balloon *Project has personal view
   customizations* is shown.
6. The customized view comes back: the `SEX` filter keeps only `M` with
   the **Filtered** count noted in step 2, the grid is sorted by
   `HEIGHT` descending, `DIS_POP` is hidden, and the histogram sits in
   the dock position it was moved to.
7. The rest of the workspace matches the pre-save state — panel
   positions and viewer sizes, with no drift.

## Cleanup

In **Browse > Dashboards**, right-click the `testCopyClone` tile, choose
**Delete Project** and click **DELETE**.

---
{
  "order": 3,
  "datasets": ["System:DemoFiles/demog.csv"]
}
