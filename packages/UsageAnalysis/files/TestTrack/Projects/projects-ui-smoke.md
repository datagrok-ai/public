---
feature: projects
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [projects.cp.upload-save-reopen-golden, projects.op.share_with_recipient_open, projects.op.rename_project]
realizes: [views.projects, sharing.share-dialog]
realized_as:
  - projects-ui-smoke-spec.ts
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - save-project-dialog
  - data-sync-toggle
  - share-dialog-dismiss
  - browse-dashboards-search
  - browse-dashboards-tile-visibility
  - pcmdShareProject
  - share-dialog-recipients
  - share-dialog-permissions-editor
  - project-double-click-open
  - pcmdDeleteProject
  - delete-confirmation-dialog
  - pcmdRename
  - rename-dialog
  - pcmdSaveAsZip
  - pcmdCopy
  - pcmdCopyId
  - pcmdCopyGrokName
  - pcmdCopyMarkup
  - pcmdCopyUrl
  - pcmdAddToFavorites
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-ui-smoke.md
migration_date: 2026-05-04
related_bugs: []
---

# Projects — UI smoke

A short pass over the main project UI: the Save dialog, the project
tile in **Browse > Dashboards**, its right-click menu, and the Share,
Rename and Delete dialogs. It takes about 5 minutes and uses one
project with `demog.csv`.

## Setup

1. Log in as the test user.
2. The project in this test is `uiSmoke` (renamed to `uiSmokeRenamed`
   in step 6).
3. Recipient for sharing: any existing user other than yourself, for
   example `qa_playwright`.

## Scenario

1. **Open a table.**
   - Go to **Browse > Files > Demo**.
   - Right-click `demog.csv` and choose **Open**.
   - **Verify:** the `demog` view opens.

2. **Save the project.**
   - Click **SAVE** on the ribbon.
   - In the **Save project** dialog, enter `uiSmoke` as the name.
   - Leave **Data sync** ON.
   - Click **OK**.
   - **Verify:** the balloon *Project "uiSmoke" uploaded.* appears.

3. **Dismiss the Share dialog.**
   - **Verify:** the **Share** dialog opens.
   - Click **CANCEL**.
   - **Verify:** the dialog closes.

4. **Find the tile.**
   - Go to **Browse > Dashboards**.
   - Type `uiSmoke` into **Search projects by name or by #tags**.
   - **Verify:** the `uiSmoke` tile is shown.

5. **Share.**
   - Right-click the `uiSmoke` tile and choose **Share...**.
   - Type the recipient into **User, group, or email**.
   - Pick the recipient from the suggestion list.
   - Leave **View and use** selected.
   - Switch **Send notifications** off.
   - Click **OK**.
   - Click the `uiSmoke` tile.
   - In the **Context Panel**, expand **Sharing**.
   - **Verify:** the recipient is listed with the words *has special
     permissions*.

6. **Rename.**
   - Right-click the `uiSmoke` tile and choose **Rename...**.
   - Change the name to `uiSmokeRenamed`.
   - Click **OK**.
   - **Verify:** the tile shows `uiSmokeRenamed`.

7. **Copy ID.**
   - Right-click the `uiSmokeRenamed` tile and choose **Copy > ID**.
   - Paste the clipboard into any text field.
   - **Verify:** the pasted text is a UUID.

8. **Copy Grok name.**
   - Right-click the tile and choose **Copy > Grok name**.
   - Paste the clipboard into any text field.
   - **Verify:** the pasted text is your login, a colon and
     `UiSmokeRenamed`.

9. **Copy Markup.**
   - Right-click the tile and choose **Copy > Markup**.
   - Paste the clipboard into any text field.
   - **Verify:** the pasted text is `#{x.` + the Grok name +
     `."uiSmokeRenamed"}`.

10. **Copy URL.**
    - Right-click the tile and choose **Copy > URL**.
    - Paste the clipboard into any text field.
    - **Verify:** the pasted text is the server address, then `/p/`,
      your login and `.UiSmokeRenamed`.

11. **Add to favorites.**
    - Right-click the tile and choose **Add to favorites**.
    - Go to **Browse > My stuff > Favorites**.
    - **Verify:** `uiSmokeRenamed` is listed.
    - Right-click `uiSmokeRenamed` and choose **Remove from favorites**.

12. **Save as Zip.**
    - Go to **Browse > Dashboards**.
    - Right-click the `uiSmokeRenamed` tile and choose **Save as Zip**.
    - **Verify:** the browser downloads `uiSmokeRenamed.zip`.

13. **Reopen.**
    - Right-click the left sidebar and select **Close All**.
    - In **Browse > Dashboards**, double-click the `uiSmokeRenamed`
      tile.
    - **Verify:** the `demog` view opens.
    - **Verify:** the status bar shows **Rows: 5,850**.

14. **Delete.**
    - Right-click the left sidebar and select **Close All**.
    - In **Browse > Dashboards**, right-click the `uiSmokeRenamed` tile
      and choose **Delete Project**.
    - **Verify:** the dialog *Are you sure? Delete project
      "uiSmokeRenamed"?* opens.
    - Click **DELETE**.
    - Wait until the dialog closes.
    - Click the refresh icon next to the search box.
    - **Verify:** the `uiSmokeRenamed` tile is gone.

## Expected results

- The project is created, shared, renamed, reopened and deleted
  entirely through the UI.
- Every right-click menu item used in the test works as described.
- After deletion the project no longer appears in **Dashboards**.
