---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.op.share_with_recipient_open, projects.op.rename_project]
realizes: [views.projects]
realized_as:
  - projects-lifecycle-files-spec.ts
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-files.md
migration_date: 2026-05-04
related_bugs: []
---

# Projects — lifecycle of a file-based project

A project built from a file (`demog.csv`) is saved, shared with a
second user at two access levels, and renamed. Both the owner and the
recipient must be able to open it at every stage.

## Setup

1. Two accounts: the **owner** (test user) and a **second user**. The
   second user must not own the project.
2. The project in this test is `lifecycleFiles`.

## Scenario

1. **Open the file.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.

2. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter `lifecycleFiles` as the name.
   - Leave **Data sync** ON.
   - Click **OK**.
   - **Verify:** the balloon *Project "lifecycleFiles" uploaded.*
     appears.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

3. **Share for viewing.**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleFiles` into the search box.
   - Right-click the `lifecycleFiles` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Switch **Send notifications** off.
   - Click **OK**.
   - **Verify:** the balloon *Shared* appears.
   - Click the `lifecycleFiles` tile.
   - In the **Context Panel**, expand **Sharing**.
   - **Verify:** the second user is listed with the words *has special
     permissions*.

4. **The recipient opens it.**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleFiles` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleFiles` tile.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - **Verify:** no error balloon appears.
   - Click **SAVE** on the ribbon.
   - **Verify:** **Save original project** is disabled.
   - **Verify:** **Save a copy** is selected.
   - Click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

5. **Grant full access.**
   - Go to **Browse > Dashboards**.
   - Right-click the `lifecycleFiles` tile and choose **Share...**.
   - **Verify:** the second user is listed with **View and use**.
   - Click **View and use** next to the second user.
   - **Verify:** a privilege tree opens with **Full access**, **View and
     use**, **Write access**, **Edit**, **Delete** and **Share**.
   - Tick **Full access**.
   - Click outside the tree.
   - **Verify:** the second user's row says **Full access**.
   - Click **OK**.
   - **Verify:** the balloon *Shared* appears.

6. **Rename the project.**
   - Right-click the `lifecycleFiles` tile and choose **Rename...**.
   - Change the name to `lifecycleFilesRenamed`.
   - Click **OK**.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleFilesRenamed` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleFilesRenamed` tile.
   - **Verify:** the `demog` view opens.
   - Right-click the left sidebar and select **Close All**.

7. **The recipient opens the renamed project.**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleFiles` into the search box.
   - **Verify:** the tile `lifecycleFilesRenamed` is shown.
   - Double-click the `lifecycleFilesRenamed` tile.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - Click **SAVE** on the ribbon.
   - **Verify:** **Save original project** is enabled and selected.
   - Click **OK**.
   - **Verify:** the balloon *Project "lifecycleFilesRenamed" uploaded.*
     appears.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

8. **Cleanup.**
   - Go to **Browse > Dashboards**.
   - Right-click the `lifecycleFilesRenamed` tile and choose **Delete
     Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.

## Expected results

- A shared file-based project opens for the recipient.
- With **View and use** the recipient cannot overwrite the project;
  with **Full access** they can.
- Renaming the project does not break sharing.
