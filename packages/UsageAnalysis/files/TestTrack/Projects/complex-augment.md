---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.upload-save-reopen-golden]
realizes: [views.projects]
realized_as:
  - complex-augment-spec.ts
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
related_bugs: []
---

# Complex — Augment project (drag-drop add a table)

Checks that a table can be added to a project that is already saved,
and that the added table is still there after the project is reopened.

## Setup

1. Log in as the test user.
2. The project in this test is `augmentTest`.

## Scenario

1. **Save a one-table project.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.
   - Click **SAVE** on the ribbon.
   - In the **Save project** dialog, enter `augmentTest` as the name.
   - Leave **Data sync** ON for `demog`.
   - Click **OK**.
   - **Verify:** the balloon *Project "augmentTest" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.

2. **Open a second table.**
   - Go to **Browse > Files > Demo**.
   - Double-click `iris.csv`.

3. **Drag the table onto the project.**
   - On the left sidebar, click the **Dashboards** icon.
   - **Verify:** the panel lists **New Dashboard** with `iris`, and the
     `augmentTest` node with `demog`.
   - Collapse the `augmentTest` node.
   - Drag `iris` from **New Dashboard** onto the `augmentTest` node.
   - **Verify:** the **Move entity** dialog opens with the heading
     *to …:AugmentTest project* and the row `iris`.
   - Click **YES**.
   - **Verify:** `iris` is listed under the `augmentTest` node.

4. **Save the project.**
   - Click **SAVE** next to the `augmentTest` node.
   - **Verify:** the **Save project** dialog lists `iris` and `demog`.
   - Leave **Save original project** selected.
   - Click **OK**.

5. **Reopen the project.**
   - Right-click the left sidebar and select **Close All**.
   - **Verify:** the `demog` and `iris` views are closed.
   - Go to **Browse > Dashboards**.
   - Type `augmentTest` into the search box.
   - Click the `augmentTest` tile.
   - In the **Context Panel**, expand **Content**.
   - **Verify:** **Content** lists `demog` and `iris`.
   - Double-click the `augmentTest` tile.
   - **Verify:** the views `demog` and `iris` open.
   - **Verify:** the status bar of `iris` shows **Rows: 150** and
     **Columns: 6**.
   - **Verify:** no error balloon appears.

6. **Cleanup.**
   - Right-click the left sidebar and select **Close All**.
   - Right-click the `augmentTest` tile and choose **Delete Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.
   - Click the refresh icon next to the search box.
   - **Verify:** the `augmentTest` tile is gone.

## Expected results

- A table dropped on an open project and then saved becomes part of
  that project.
- After close and reopen, the project opens both tables with their
  data.
