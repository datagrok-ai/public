---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.cp.share-spaces-datasync, projects.op.rename_external_dep, projects.op.share_with_recipient_open, projects.op.rename_project]
realizes: [views.projects, views.space]
realized_as:
  - projects-lifecycle-spaces-spec.ts
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-spaces.md
migration_date: 2026-05-04
related_bugs:
  - GROK-18345
---

# Projects — lifecycle of a Space-based project

A project built from a file stored in a Space is saved with Data
sync, shared with a second user, and then the Space and the project
are renamed. This reproduces GROK-18345: the recipient could not open
a shared, data-synced project whose table comes from a Space.

## Setup

1. Two accounts: the **owner** (test user) and a **second user**.
2. Names in this test: Space `lifecycleSpace`, project
   `lifecycleSpaceProj`.

## Scenario

1. **Create the Space.**
   - In **Browse**, right-click **Spaces** and choose **Create
     Space...**.
   - Replace *New Space* with `lifecycleSpace`.
   - Click **OK**.
   - Expand **Spaces**.
   - **Verify:** `lifecycleSpace` is listed.

2. **Copy a file into it.**
   - Go to **Browse > Files > Demo**.
   - Drag `demog.csv` onto the `lifecycleSpace` Space in the Browse
     tree.
   - **Verify:** the **Move entity** dialog opens with **Move**,
     **Link** and **Copy**, and **Link** is selected.
   - Select **Copy**.
   - Click **YES**.
   - Click the `lifecycleSpace` Space.
   - **Verify:** `demog.csv` is listed.
   - Click **Browse > Files > Demo**.
   - **Verify:** `demog.csv` is still listed.

3. **Open the file from the Space.**
   - Click the `lifecycleSpace` Space.
   - Double-click `demog.csv`.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - **Verify:** the URL ends with `/file/LifecycleSpace.Files/demog.csv`.

4. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter `lifecycleSpaceProj` as the name.
   - Leave **Data sync** ON.
   - Expand **CREATION SCRIPT**.
   - **Verify:** it shows `OpenFile("LifecycleSpace:Files/demog.csv")`.
   - Click **OK**.
   - **Verify:** the balloon *Project "lifecycleSpaceProj" uploaded.*
     appears.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

5. **Share the project only.**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleSpaceProj` into the search box.
   - Right-click the `lifecycleSpaceProj` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.

6. **The recipient opens it (GROK-18345).**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleSpaceProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleSpaceProj` tile.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - **Verify:** no error dialog appears.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

7. **Rename the Space.**
   - In **Browse > Spaces**, right-click `lifecycleSpace` and choose
     **Rename...**.
   - Change the name to `lifecycleSpaceRenamed`.
   - Click **OK**.

8. **Both users reopen the project.**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleSpaceProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleSpaceProj` tile.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleSpaceProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleSpaceProj` tile.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

9. **Rename the project.**
   - Go to **Browse > Dashboards**.
   - Right-click the `lifecycleSpaceProj` tile and choose **Rename...**.
   - Change the name to `lifecycleSpaceProjRenamed`.
   - Click **OK**.

10. **Both users open the renamed project.**
    - Go to **Browse > Dashboards**.
    - Type `lifecycleSpaceProjRenamed` into the search box.
    - Click the refresh icon.
    - Double-click the `lifecycleSpaceProjRenamed` tile.
    - **Verify:** the `demog` view opens with 5,850 rows.
    - Right-click the left sidebar and select **Close All**.
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the second user's credentials.
    - Go to **Browse > Dashboards**.
    - Type `lifecycleSpaceProjRenamed` into the search box.
    - Click the refresh icon.
    - Double-click the `lifecycleSpaceProjRenamed` tile.
    - **Verify:** the `demog` view opens with 5,850 rows.
    - Right-click the left sidebar and select **Close All**.
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the owner's credentials.

11. **Cleanup.**
    - Go to **Browse > Dashboards**.
    - Right-click the `lifecycleSpaceProjRenamed` tile and choose
      **Delete Project**.
    - Click **DELETE**.
    - Wait until the dialog closes.
    - In **Browse > Spaces**, right-click `lifecycleSpaceRenamed` and
      choose **Delete Space**.
    - **Verify:** the dialog says *Delete space "lifecycleSpaceRenamed"?
      This will delete space and its related data…*.
    - Click **DELETE**.

## Expected results

- A project built from a Space file reopens with Data sync.
- The recipient can open the shared project and sees the data.
- Renaming the Space and the project does not break the project.
