---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: [views.projects, views.space, sharing.share-dialog]
realized_as:
  - complex-move-spec.ts
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
related_bugs: []
---

# Complex — Move a project between Spaces, and access through a Space

Checks that a saved project can be moved into a Space, into its child
Space and from there into another Space, and that it still opens with
its data after each move. A project in a Space takes the Space's
permissions, and a child Space inherits them from its root: the second
user gets access to the project only by the Space being shared with
them — the project itself is never shared. Moving the project to a
Space that is not shared takes the access away.

## Setup

1. Two accounts: the **owner** (test user) and a **second user**.
2. Names in this test: project `moveTest`, Spaces `moveA` (with the
   child Space `moveAChild`) and `moveB`.

## Scenario

1. **Create two Spaces.**
   - In **Browse**, right-click **Spaces** and choose **Create
     Space...**.
   - In the **Create Space** dialog, replace *New Space* with `moveA`.
   - Click **OK**.
   - Right-click **Spaces** and choose **Create Space...**.
   - Replace *New Space* with `moveB`.
   - Click **OK**.
   - Expand **Spaces**.
   - **Verify:** `moveA` and `moveB` are listed.
   - Right-click `moveA` and choose **Create Child Space...**.
   - Replace *New Space* with `moveAChild`.
   - Click **OK**.

2. **Save the project.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.
   - Click **SAVE** on the ribbon.
   - Enter `moveTest` as the name.
   - Leave **Data sync** ON.
   - Click **OK**.
   - **Verify:** the balloon *Project "moveTest" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

3. **Move it to the first Space.**
   - Go to **Browse > Dashboards**.
   - Type `moveTest` into the search box.
   - Right-click the `moveTest` tile and choose **Move to Space...**.
   - In the **Move to space** dialog, set **Space** to `moveA`.
   - Click **OK**.
   - Click **Spaces > moveA**.
   - **Verify:** `moveTest` is listed in `moveA`.

4. **Open it from the Space.**
   - Double-click `moveTest` in `moveA`.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - **Verify:** no error balloon appears.
   - Right-click the left sidebar and select **Close All**.

5. **Share the Space, not the project.**
   - In **Browse > Spaces**, right-click `moveA` and choose
     **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.
   - In `moveA`, right-click `moveTest` and choose **Share...**.
   - **Verify:** the dialog shows *Inherited from* `moveA`, with the
     second user under it.
   - Click **CANCEL**.

6. **The second user opens the project.**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Spaces > moveA**.
   - **Verify:** `moveTest` is listed.
   - Double-click `moveTest`.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - **Verify:** no error dialog appears.
   - Click **SAVE** on the ribbon.
   - **Verify:** **Save original project** is disabled.
   - Click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

7. **Move it into the child Space.**
   - Go to **Browse > Spaces > moveA**.
   - Drag `moveTest` onto `moveAChild` in the Browse tree.
   - In the **Move entity** dialog, select **Move**.
   - Click **YES**.
   - Click **Spaces > moveA > moveAChild**.
   - **Verify:** `moveTest` is listed.

8. **The second user opens it in the child Space.**
   - Sign out and sign in with the second user's credentials, as in
     step 6.
   - Go to **Browse > Spaces > moveA > moveAChild**.
   - Double-click `moveTest`.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - Right-click the left sidebar and select **Close All**.
   - Sign out and sign in with the owner's credentials.

9. **Move it to the second Space.**
   - In `moveAChild`, right-click `moveTest` and choose **Move to
     Space...**.
   - Set **Space** to `moveB`.
   - Click **OK**.
   - Click **Spaces > moveB**.
   - **Verify:** `moveTest` is listed in `moveB`.
   - Click **Spaces > moveA > moveAChild**.
   - **Verify:** `moveTest` is not listed in `moveAChild`.

10. **Open it again.**
    - Click **Spaces > moveB**.
    - Double-click `moveTest`.
    - **Verify:** the `demog` view opens with 5,850 rows.
    - Right-click the left sidebar and select **Close All**.

11. **The second user no longer sees it.**
    - Sign out and sign in with the second user's credentials.
    - Go to **Browse > Dashboards**.
    - Type `moveTest` into the search box.
    - Click the refresh icon.
    - **Verify:** no `moveTest` tile is shown.
    - Go to **Browse > Spaces**.
    - **Verify:** `moveB` is not listed, and `moveA > moveAChild` does
      not list `moveTest`.
    - Sign out and sign in with the owner's credentials.

12. **Cleanup.**
    - In `moveB`, right-click `moveTest` and choose **Delete Project**.
    - Click **DELETE**.
    - Wait until the dialog closes.
    - Right-click `moveA` and choose **Delete Space**.
    - **Verify:** the dialog says *Delete space "moveA"? This will
      delete space and its related data…*.
    - Click **DELETE**.
    - Right-click `moveB` and choose **Delete Space**.
    - Click **DELETE**.

## Expected results

- A project can be moved into a Space and between Spaces.
- The project opens with its data after every move.
- A project in a shared Space opens for the users the Space is shared
  with, without being shared itself.
- A child Space passes the access on from its root.
- Moving the project to a Space that is not shared takes the access
  away.

## Automation notes

- **Move to Space...** lists only root Spaces (`project_meta.dart`,
  `getRootSpaces`), so step 7 drags the project onto the child Space;
  the drag of a project within the Browse tree was not run by hand
  while writing this file.
- The *Inherited from* line in the Share dialog is read from
  `permissions_browser.dart`.
- Deleting `moveA` in the cleanup is expected to delete its child Space
  `moveAChild` with it; if `moveAChild` is still listed afterwards,
  delete it the same way.
