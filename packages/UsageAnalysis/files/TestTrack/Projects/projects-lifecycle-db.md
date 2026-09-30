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
creation script. The table project is also shared. Sharing a project
does not share the database connection behind it: a recipient without
**View and use** on the connection cannot load the data and must be
told so. Once the connection is shared with them, the project opens
with data.

## Setup

1. Two accounts: the **owner** (test user) and a **second user**.
2. **NorthwindTest** exists only on dev.
3. At the start, the second user must not have access to the
   **NorthwindTest** Postgres connection: signed in as the second user,
   **Browse > Databases > Postgres** does not list **NorthwindTest**.
   Test 2 shares the connection with them.
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

5. **The recipient opens it without access to the connection.**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbTable` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDbTable` tile.
   - **Verify:** the **Data loading error** dialog opens and shows the
     `DbQuery(...)` line with `public.orders`; it offers **OPEN ANYWAY**
     and **CLOSE PROJECT**, and no **EDIT SCRIPT...**.
   - **Verify:** the page does not stay on a loading spinner, and no
     empty table opens without a message.
   - Click **CLOSE PROJECT**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

6. **Share the connection.**
   - Go to **Browse > Databases > Postgres**.
   - Right-click **NorthwindTest** and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.

7. **The recipient opens it with access.**
   - Sign out and sign in with the second user's credentials, as in
     step 5.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleDbTable` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleDbTable` tile.
   - **Verify:** the `orders` table opens with 830 rows.
   - **Verify:** no error dialog appears.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

8. **Cleanup.**
   - Go to **Browse > Databases > Postgres**, right-click
     **NorthwindTest** and choose **Share...**. Next to the second user,
     click the cross icon (tooltip *Remove*). Click **OK**.
   - In **Browse > Dashboards**, type `lifecycleDbTable` into the search
     box.
   - Right-click the `lifecycleDbTable` tile and choose **Delete
     Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.

## Expected results

- A project with a query or a database table reopens with its data and
  keeps Data sync and its creation script.
- Sharing a project does not share the database connection behind it.
- A recipient without access to the connection gets a clear error, not
  an endless spinner or an empty table.
- After the connection is shared, the recipient opens the shared table
  project with data.

## Automation notes

- The tests cannot change the NorthwindTest data, so they do not prove
  that the reopened table was read from the database again rather than
  from a snapshot; the Save dialog of the reopened project is the
  closest observable check. `projects-lifecycle-query.md` proves the
  re-read by changing its own query's SQL.
- The help page `help/datagrok/concepts/project/dashboard.md` says that
  a recipient without the connection permission cannot run the query;
  its planned illustration shows the page staying on the loading
  spinner. What the product shows today (Test 2, step 5) was not
  confirmed from source; record it on the first run.
- The cross icon with the tooltip *Remove* in the Share dialog is read
  from `permissions_browser.dart`.
