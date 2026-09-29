---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.cp.rename-dependent-entity-reopen, projects.op.rename_external_dep, projects.op.share_with_recipient_open]
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
(github-3550). Then the query's SQL is changed: the next open must show
the new result, which proves that Data sync runs the query again
instead of showing the saved rows.

## Setup

1. Two accounts: the **owner** (test user) and a **second user** who
   can access the **NorthwindTest** Postgres connection.
2. **NorthwindTest** exists only on dev.
3. Names in this test: query `lifecycleQuery`, project
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

8. **Change the query's SQL.**
   - Under **NorthwindTest**, right-click `lifecycleQueryRenamed` and
     choose **Edit...**.
   - Change the query text to `select * from public.orders limit 100`.
   - Click **Save** on the ribbon.
   - Close the query editor.

9. **The owner reopens the project with the new SQL.**
   - Reload the browser tab.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleQueryProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleQueryProj` tile.
   - **Verify:** the table opens with 100 rows (not the 830 it was
     saved with).
   - **Verify:** no error dialog appears.
   - Right-click the left sidebar and select **Close All**.

10. **The recipient reopens the project (github-3550).**
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the second user's credentials.
    - Go to **Browse > Dashboards**.
    - Type `lifecycleQueryProj` into the search box.
    - Click the refresh icon.
    - Double-click the `lifecycleQueryProj` tile.
    - **Verify:** the table opens with 100 rows.
    - **Verify:** no error dialog appears.
    - Right-click the left sidebar and select **Close All**.
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the owner's credentials.

11. **Cleanup.**
    - Go to **Browse > Dashboards**.
    - Right-click the `lifecycleQueryProj` tile and choose **Delete
      Project**.
    - Click **DELETE**.
    - Wait until the dialog closes.
    - Under **NorthwindTest**, right-click `lifecycleQueryRenamed` and
      choose **Delete**.
    - **Verify:** the dialog says *Delete query
      "lifecycleQueryRenamed"?*.
    - Click **DELETE**.

## Expected results

- A query-based project reopens and runs the query again: after the
  SQL changes, the owner and the recipient see the new result.
- The recipient can open the shared project.
- After the query is renamed, the project still opens with its data.

## Automation notes

- Step 9 reloads the tab so that the owner's session does not reuse a
  query definition it loaded before the edit: after a script is edited,
  the same session can keep running the old body for about a minute.
  Whether a query edit has the same delay was not checked.
