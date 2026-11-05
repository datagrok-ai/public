---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.op.share_with_recipient_open]
realizes: [views.projects, sharing.share-dialog]
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
from a database table opened directly. Each is saved with Data sync and
reopened; the reopened project must still carry Data sync and its
creation script. The table project is also shared: sharing a project
of a database table brings its connection along, so the recipient sees
the connection and opens the project with the data (intended, operator
2026-10-05; a project of a saved query does not bring its connection —
see `projects-lifecycle-query.md`).

## Setup

1. Two accounts: the **owner** (test user) and a **second user**.
2. **NorthwindTest** exists only on dev.
3. Before the test, the second user has no access to the
   **NorthwindTest** Postgres connection: signed in as the second user,
   **Browse > Databases > Postgres** does not list **NorthwindTest**. The
   connection is not shared by hand in this test.
4. Projects in this test: `lifecycleDbQuery` and `lifecycleDbTable`.

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
   - Click **SAVE** on the ribbon.
   - **Verify:** **Data sync** is ON for the table, and **CREATION
     SCRIPT** shows the call of **PostgresAll**.
   - Click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

4. **Cleanup.**
   - Right-click the `lifecycleDbQuery` tile and choose **Delete
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
   - Click **SAVE** on the ribbon.
   - **Verify:** **Data sync** is ON for `orders`, and **CREATION
     SCRIPT** shows the `DbQuery(...)` call.
   - Click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

4. **Share the project only.**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbTable` into the search box.
   - Right-click the `lifecycleDbTable` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.

5. **The recipient sees the connection and opens the project.**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Databases > Postgres**.
   - **Verify:** **NorthwindTest** is listed now: the share of the
     table project brought its connection along.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbTable` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDbTable` tile.
   - **Verify:** the `orders` table opens with 830 rows, and no error
     dialog appears.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

6. **Cleanup.**
   - In **Browse > Dashboards**, type `lifecycleDbTable` into the search
     box.
   - Right-click the `lifecycleDbTable` tile and choose **Delete
     Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.

## Expected results

- A project with a query or a database table reopens with its data and
  keeps Data sync and its creation script.
- Sharing a database-table project brings its connection along: the
  recipient, who could not see the connection before, sees it and opens
  the project with its data.

## Automation notes

- The tests cannot change the NorthwindTest data, so they do not prove
  that the reopened table was read from the database again rather than
  from a snapshot; the Save dialog of the reopened project is the
  closest observable check. `projects-lifecycle-query.md` proves the
  re-read by changing its own query's SQL.
- The help page `help/datagrok/concepts/project/dashboard.md` says that
  a recipient without the connection permission cannot run the query;
  that is the saved-query case (`projects-lifecycle-query.md`). On dev
  (2026-10-12) the second user, with no grant of its own on NorthwindTest,
  does not see the connection before the share of Test 2 and sees it
  after; the project opens with its 830 rows (intended, operator
  2026-10-05).
