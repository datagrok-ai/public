---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.upload-save-reopen-golden]
realizes: [views.projects, file.menu.save.tables-as-project]
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

# Complex — Augment project (leave a table out, then drag-drop add it)

Every table row in the **Save project** dialog has a cross icon that
leaves the table out of the project, and a plus icon that brings it
back. The project is first saved with one of two open tables left out,
and reopens without it. Then that table is added to the saved project
by drag and drop, and it is still there after the project is reopened.

## Setup

1. Log in as the test user.
2. The project in this test is `augmentTest`.

## Scenario

1. **Open two tables.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.
   - Double-click `iris.csv`.

2. **Save the project with `iris` left out.**
   - Click **SAVE** on the ribbon.
   - **Verify:** the **Save project** dialog lists `demog` and `iris`,
     each with a cross icon at the end of its row.
   - Hover the cross icon in the `iris` row.
   - **Verify:** the tooltip reads *Exclude table from the project.*
   - Click it.
   - **Verify:** the `iris` row is greyed out, and a plus icon (tooltip
     *Include table to the project.*) replaces the cross.
   - Click the plus icon.
   - **Verify:** the `iris` row is active again, with the cross icon.
   - Click the cross icon in the `iris` row again.
   - Enter `augmentTest` as the name.
   - Leave **Data sync** ON for `demog`.
   - Click **OK**.
   - **Verify:** the balloon *Project "augmentTest" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.

3. **The project holds only `demog`.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `augmentTest` into the search box.
   - Click the refresh icon.
   - Click the `augmentTest` tile.
   - In the **Context Panel**, expand **Content**.
   - **Verify:** **Content** lists `demog` and not `iris`.
   - Double-click the `augmentTest` tile.
   - **Verify:** only the `demog` view opens, with 5,850 rows.

4. **Open the second table again.**
   - Go to **Browse > Files > Demo**.
   - Double-click `iris.csv`.

5. **Drag the table onto the project.**
   - On the left sidebar, click the **Dashboards** icon.
   - **Verify:** the panel lists **New Dashboard** with `iris`, and the
     `augmentTest` node with `demog`.
   - Collapse the `augmentTest` node.
   - Drag `iris` from **New Dashboard** onto the `augmentTest` node.
   - **Verify:** the **Move entity** dialog opens with the heading
     *to …:AugmentTest project* and the row `iris`.
   - Click **YES**.
   - **Verify:** `iris` is listed under the `augmentTest` node.

6. **Save the project.**
   - Click **SAVE** next to the `augmentTest` node.
   - **Verify:** the **Save project** dialog lists `iris` and `demog`.
   - **Verify:** **Data sync** is OFF for `iris` and ON for `demog`. A
     table dropped into an open project arrives with Data sync OFF.
   - Switch **Data sync** ON for `iris`.
   - Expand **CREATION SCRIPT** under `iris`.
   - **Verify:** it shows `OpenFile("System:DemoFiles/iris.csv")`.
   - Leave **Save original project** selected.
   - Click **OK**.

7. **Reopen the project.**
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
   - Click **SAVE** on the ribbon.
   - **Verify:** **Data sync** is ON for both `demog` and `iris`.
   - Click **CANCEL**.

8. **Cleanup.**
   - Right-click the left sidebar and select **Close All**.
   - Right-click the `augmentTest` tile and choose **Delete Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.
   - Click the refresh icon next to the search box.
   - **Verify:** the `augmentTest` tile is gone.

## Expected results

- The cross icon leaves a table out of the project and the plus icon
  brings it back before saving.
- A table left out is not saved; its view is not saved either.
- A table dropped on an open project and then saved becomes part of
  that project.
- A dropped table arrives with Data sync OFF; once it is switched on,
  the table is saved with its creation script.
- After close and reopen, the project opens both tables with their
  data, and both keep Data sync ON.

## Automation notes

- The icons, their tooltips and the greyed row are read from
  `project_entity_move.dart`; that the view of a left-out table is not
  saved follows from the same file (a view takes its table's action).
