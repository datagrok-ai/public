---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.op.share_with_recipient_open, projects.op.rename_project]
realizes: [views.projects]
realized_as:
  - projects-lifecycle-db-spec.ts
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-db.md
migration_date: 2026-05-04
related_bugs: []
---

# Projects — lifecycle of a database-based project

Two projects built from a database: one from a saved query and one
from a database table opened directly. Each is saved with Data sync,
reopened, shared with a second user and renamed. On reopening, the
data must be re-read from the database.

## Setup

1. Two accounts: the **owner** (test user) and a **second user** who
   can access the **NorthwindTest** Postgres connection.
2. Projects in this test: `lifecycleDbQuery` and `lifecycleDbTable`.

## Scenario

### Test 1 — saved query

1. **Run the query.**
   - Go to **Browse > Databases > Postgres > NorthwindTest**.
   - Double-click the query **PostgresAll**.
   - **Verify:** the result view opens with 830 rows.

2. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter `lifecycleDbQuery` as the name.
   - Leave **Data sync** ON.
   - Expand **CREATION SCRIPT**.
   - **Verify:** it shows a call of the **PostgresAll** query.
   - Click **OK**.
   - **Verify:** the balloon *Project "lifecycleDbQuery" uploaded.*
     appears.
   - In the **Share** dialog, click **CANCEL**.

3. **Reopen.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbQuery` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDbQuery` tile.
   - **Verify:** the table opens with 830 rows.
   - **Verify:** no error balloon appears.
   - Right-click the left sidebar and select **Close All**.

4. **Share.**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbQuery` into the search box.
   - Right-click the `lifecycleDbQuery` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.

5. **The recipient opens it.**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbQuery` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDbQuery` tile.
   - **Verify:** the `PostgresAll` table opens with 830 rows.
   - **Verify:** no error dialog appears.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

6. **Rename.**
   - Go to **Browse > Dashboards**.
   - Right-click the `lifecycleDbQuery` tile and choose **Rename...**.
   - Change the name to `lifecycleDbQueryRenamed`.
   - Click **OK**.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbQueryRenamed` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDbQueryRenamed` tile.
   - **Verify:** the table opens with 830 rows.
   - Right-click the left sidebar and select **Close All**.

7. **Cleanup.**
   - Right-click the `lifecycleDbQueryRenamed` tile and choose **Delete
     Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.

### Test 2 — database table

1. **Open the table.**
   - Go to **Browse > Databases > Postgres > NorthwindTest > Schemas >
     public**.
   - Right-click `orders` and choose **Get All**.
   - **Verify:** the `orders` view opens with 830 rows.

2. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter `lifecycleDbTable` as the name.
   - Leave **Data sync** ON.
   - Expand **CREATION SCRIPT**.
   - **Verify:** it shows `DbQuery(Dbtests:PostgresTest,
     "public.orders", …)`.
   - Click **OK**.
   - In the **Share** dialog, click **CANCEL**.

3. **Reopen.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbTable` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDbTable` tile.
   - **Verify:** the `orders` table opens with 830 rows.
   - Right-click the left sidebar and select **Close All**.

4. **Share.**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbTable` into the search box.
   - Right-click the `lifecycleDbTable` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.

5. **The recipient opens it.**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbTable` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDbTable` tile.
   - **Verify:** the `orders` table opens with 830 rows.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

6. **Cleanup.**
   - Right-click the `lifecycleDbTable` tile and choose **Delete
     Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.

## Expected results

- A project with a query or a database table reopens and reloads its
  data from the database.
- The recipient can open the shared project.
- Renaming the project does not break it.
