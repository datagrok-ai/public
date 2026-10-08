---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.int.derive-then-save-inside-project]
realizes: [views.projects, data.menu.join-tables]
realized_as:
  - projects-lifecycle-derived-spec.ts
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-derived.md
migration_date: 2026-05-04
related_bugs:
  - GROK-19103
---

# Projects — lifecycle of a project with derived tables

A project holds a source table plus three tables built from it: a
pivot, an aggregate and a join. It is saved and reopened, and every
derived table must come back with its rows. This also covers
GROK-19103: the join result must stay in the current workspace and
must not be saved as a separate, broken project.

## Setup

1. Log in as the test user.
2. The project in this test is `lifecycleDerived`.

## Scenario

1. **Open the source.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.

2. **Pivot.**
   - In **Toolbox > Viewers**, click **Pivot table**.
   - Set **Group by** to `RACE`.
   - Set **Pivot** to `SEX`.
   - Set **Aggregate** to `avg(AGE)`.
   - Click **ADD** in the top-right corner of the viewer.
   - **Verify:** the view *demog aggregation* opens.
   - Write down its row count from the status bar.

3. **Aggregate.**
   - Click the `demog` view tab.
   - Open **Data > Aggregate Rows...**.
   - **Verify:** the **Pivot table** panel opens.
   - Set **Group by** to `DIS_POP`.
   - Set **Aggregate** to `count(USUBJID)`.
   - Clear **Pivot**.
   - Click **ADD**.
   - **Verify:** a new view with the aggregated table opens with 6 rows,
     one per `DIS_POP` value.

4. **Join (GROK-19103).**
   - Open **Data > Join Tables...**.
   - In **Tables**, set the first table to `demog`.
   - Set the second table to *demog aggregation*.
   - Set both **Key Columns** to `RACE`.
   - Set **Join Type** to `inner`.
   - Click **OK**.
   - **Verify:** the join result opens as a new view.
   - Write down its row count.
   - On the left sidebar, click the **Dashboards** icon.
   - **Verify:** all four tables are listed under **New Dashboard**.
   - **Verify:** no other project node is listed.

5. **Save.**
   - Click **SAVE** on the ribbon.
   - **Verify:** each of the four tables has its own **Data sync**
     toggle and **CREATION SCRIPT**.
   - Enter `lifecycleDerived` as the name.
   - Leave **Data sync** ON for every table.
   - Click **OK**.
   - In the **Share** dialog, click **CANCEL**.

6. **Reopen.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDerived` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDerived` tile.
   - **Verify:** four views open: `demog`, the pivot, the aggregate and
     the join.
   - **Verify:** `demog` has 5,850 rows, the aggregate has 6 rows, and
     the pivot and the join have the row counts written down in steps
     2 and 4.
   - **Verify:** no error balloon appears.
   - Right-click the left sidebar and select **Close All**.

7. **Cleanup.**
   - Right-click the `lifecycleDerived` tile and choose **Delete
     Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.

## Expected results

- Pivot, aggregate and join results are saved in the project and
  reopen with the same rows.
- The join result is part of the project, not a separate project.
