---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.save-copy-with-link-mode]
realizes: [views.projects]
realized_as:
  - complex-save-copy-spec.ts
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
related_bugs: []
---

# Complex — Save a copy without Data sync, then turn Data sync back on

Checks the Save-a-copy round trip:

1. Save a project.
2. Save a copy of it with **Data sync OFF** (a static snapshot).
3. Reopen the copy.
4. Save it again with **Data sync ON**, so it reads the file again on
   every open.

The original project must not change.

## Setup

1. Log in as the test user.
2. Projects in this test: `saveCopyOriginal` and `saveCopyNoSync`.
3. The **Save project** dialog shows a **CREATION SCRIPT** block under
   a table only when its Data sync is ON.

## Scenario

1. **Save the original.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.
   - Click **SAVE** on the ribbon.
   - Enter `saveCopyOriginal` as the name.
   - Leave **Data sync** ON.
   - Click **OK**.
   - **Verify:** the balloon *Project "saveCopyOriginal" uploaded.*
     appears.
   - In the **Share** dialog, click **CANCEL**.

2. **Save a copy without Data sync.**
   - Click **SAVE** on the ribbon.
   - In the **Save project** dialog, select **Save a copy**.
   - **Verify:** the name changes to *Copy of saveCopyOriginal*.
   - **Verify:** a **Link / Clone** choice with **Clone** selected
     appears next to `demog`.
   - Enter `saveCopyNoSync` as the name.
   - Switch **Data sync** OFF.
   - Click **OK**.
   - **Verify:** the balloon *Project "saveCopyNoSync" uploaded.*
     appears.
   - On the left sidebar, click the **Dashboards** icon.
   - **Verify:** the open project node is `saveCopyNoSync`.

3. **Reopen the copy.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `saveCopy` into the search box.
   - Click the refresh icon.
   - **Verify:** the tiles `saveCopyOriginal` and `saveCopyNoSync` are
     shown.
   - Double-click the `saveCopyNoSync` tile.
   - **Verify:** the `demog` view opens.
   - Click **SAVE** on the ribbon.
   - **Verify:** the dialog shows no **CREATION SCRIPT**.
   - Click **CANCEL**.

4. **Turn Data sync on for the copy.**
   - Click **SAVE** on the ribbon.
   - Leave **Save original project** selected.
   - Switch **Data sync** ON.
   - Click **OK**.

5. **Check the copy.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `saveCopyNoSync` into the search box.
   - Click the refresh icon.
   - Double-click the `saveCopyNoSync` tile.
   - Click **SAVE** on the ribbon.
   - Expand **CREATION SCRIPT**.
   - **Verify:** it shows `OpenFile("System:DemoFiles/demog.csv")`.
   - Click **CANCEL**.

6. **Check the original.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `saveCopyOriginal` into the search box.
   - Click the refresh icon.
   - Double-click the `saveCopyOriginal` tile.
   - **Verify:** the `demog` view opens.
   - Click **SAVE** on the ribbon.
   - **Verify:** the dialog shows **CREATION SCRIPT**.
   - Click **CANCEL**.

7. **Cleanup.**
   - Right-click the left sidebar and select **Close All**.
   - In **Browse > Dashboards**, right-click the `saveCopyOriginal` tile
     and choose **Delete Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.
   - Right-click the `saveCopyNoSync` tile and choose **Delete
     Project**.
   - Click **DELETE**.

## Expected results

- **Save a copy** creates a separate project and leaves the original
  unchanged.
- A copy saved with Data sync OFF opens as a static snapshot.
- Re-saving the copy with Data sync ON makes it a live, file-backed
  project.
