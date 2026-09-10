---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.cp.rename-dependent-entity-reopen, projects.op.rename_external_dep, projects.op.share_with_recipient_open, projects.op.rename_project]
realizes: [views.projects]
realized_as:
  - projects-lifecycle-query-spec.ts
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-query.md
migration_date: 2026-05-04
related_bugs:
  - github-3550
---

# Projects — lifecycle of a query-based project

A project built from your own saved query is saved, shared, and then
the **query** is renamed. The project must still open with its data
(github-3550).

## Setup

1. Two accounts: the **owner** (test user) and a **second user** who
   can access the **NorthwindTest** Postgres connection.
2. Names in this test: query `lifecycleQuery`, project
   `lifecycleQueryProj`.

## Scenario

1. **Create the query.**
   - Go to **Browse > Databases > Postgres > NorthwindTest > Schemas >
     public**.
   - Right-click `orders` and choose **New SQL Query...**.
   - **Verify:** the query editor opens with `select * from
     public.orders`.
   - Replace the **Name** `orders` with `lifecycleQuery`.
   - Click **Save** on the ribbon.
   - Close the query editor.

2. **Open its result as a table.**
   - If `lifecycleQuery` is not shown under **NorthwindTest**, click
     **Load more** at the end of the list.
   - Double-click `lifecycleQuery`.
   - **Verify:** the view `lifecycleQuery` opens with 830 rows.

3. **Save the project.**
   - Click **SAVE** on the ribbon.
   - **Verify:** the **Save project** dialog lists only
     `lifecycleQuery`.
   - Enter `lifecycleQueryProj` as the name.
   - Leave **Data sync** ON.
   - Click **OK**.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

4. **Share the project only.**
   - Go to **Browse > Dashboards**.
   - Right-click the `lifecycleQueryProj` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Switch **Send notifications** off.
   - Click **OK**.

5. **The recipient opens it.**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleQueryProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleQueryProj` tile.
   - **Verify:** the table opens with 830 rows.
   - **Verify:** no error dialog appears.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

6. **Rename the query.**
   - Go to **Browse > Databases > Postgres > NorthwindTest**.
   - Right-click `lifecycleQuery` and choose **Rename...**.
   - In the **Rename dataquery** dialog, change the name to
     `lifecycleQueryRenamed`.
   - Click **OK**.

7. **The owner reopens the project (github-3550).**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleQueryProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleQueryProj` tile.
   - **Verify:** the table opens with 830 rows.
   - **Verify:** no error dialog appears.
   - Right-click the left sidebar and select **Close All**.

8. **The recipient reopens the project (github-3550).**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleQueryProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleQueryProj` tile.
   - **Verify:** the table opens with 830 rows.
   - **Verify:** no error dialog appears.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

9. **Rename the project.**
   - Go to **Browse > Dashboards**.
   - Right-click the `lifecycleQueryProj` tile and choose **Rename...**.
   - Change the name to `lifecycleQueryProjRenamed`.
   - Click **OK**.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleQueryProjRenamed` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleQueryProjRenamed` tile.
   - **Verify:** the table opens with 830 rows.
   - Right-click the left sidebar and select **Close All**.

10. **Cleanup.**
    - Right-click the `lifecycleQueryProjRenamed` tile and choose
      **Delete Project**.
    - Click **DELETE**.
    - Wait until the dialog closes.
    - Under **NorthwindTest**, right-click `lifecycleQueryRenamed` and
      choose **Delete**.
    - **Verify:** the dialog says *Delete query
      "lifecycleQueryRenamed"?*.
    - Click **DELETE**.

## Expected results

- A query-based project reopens and reloads the query result.
- The recipient can open the shared project.
- After the query is renamed, the project still opens with its data.
