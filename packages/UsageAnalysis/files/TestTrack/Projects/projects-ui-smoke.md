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
Rename and Delete dialogs. The project's description is entered in the
**Save project** dialog and shown on its tile and in the **Context
Panel**; a tag is added in **Context Panel > Details > Tags**, and the
**Dashboards** search finds the project by `#tag`. The test uses one
main project and a second one without a description, both with
`demog.csv`.

## Setup

1. Log in as the test user.
2. Names in this test: projects `uiSmoke` (renamed to `uiSmokeRenamed`
   in step 9) and `uiSmokeOther`; tag `qatag1`.
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
   - Click the **Description** field and type `Demog snapshot for the
     UI smoke test`.
   - Leave **Data sync** ON.
   - Click **OK**.
   - **Verify:** the balloon *Project "uiSmoke" uploaded.* appears.

3. **Dismiss the Share dialog.**
   - **Verify:** the **Share** dialog opens.
   - Click **CANCEL**.
   - **Verify:** the dialog closes.

4. **Save a second project without a description.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Files > Demo** and double-click `demog.csv`.
   - Click **SAVE**, enter `uiSmokeOther` without a description, click
     **OK**, and **CANCEL** in **Share**.
   - Right-click the left sidebar and select **Close All**.

5. **Find the tile and its description.**
   - Go to **Browse > Dashboards**.
   - Type `uiSmoke` into **Search projects by name or by #tags**.
   - **Verify:** the `uiSmoke` tile is shown.
   - **Verify:** the `uiSmoke` tile shows *Demog snapshot for the UI
     smoke test* over its thumbnail; the `uiSmokeOther` tile shows no
     text there.
   - Click the `uiSmoke` tile.
   - **Verify:** the **Context Panel** shows the description under the
     project picture in **Details**.

6. **Add a tag.**
   - In **Details**, next to **Tags**, click the pencil icon (tooltip
     *Edit tags and metadata*).
   - Click the **Add tag** field, type `qatag1` and press Enter.
   - Press Esc.
   - **Verify:** **Tags** shows `#qatag1`.

7. **Find the project by its tag.**
   - Clear the search box and type `#qatag1`.
   - **Verify:** only the `uiSmoke` tile is shown.
   - Clear the search box and type `#qatag2`.
   - **Verify:** no tile is shown.
   - Clear the search box and type `uiSmoke`.

8. **Share.**
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

9. **Rename.**
   - Right-click the `uiSmoke` tile and choose **Rename...**.
   - Change the name to `uiSmokeRenamed`.
   - Click **OK**.
   - **Verify:** the tile shows `uiSmokeRenamed`.

10. **The tag and the description survive a reload.**
    - Reload the browser tab.
    - Go to **Browse > Dashboards**, type `#qatag1` into the search box.
    - Click the `uiSmokeRenamed` tile.
    - **Verify:** **Details** shows the description and `#qatag1`.

11. **Copy ID.**
    - Right-click the `uiSmokeRenamed` tile and choose **Copy > ID**.
    - Paste the clipboard into any text field.
    - **Verify:** the pasted text is a UUID.

12. **Copy Grok name.**
    - Right-click the tile and choose **Copy > Grok name**.
    - Paste the clipboard into any text field.
    - **Verify:** the pasted text is the project's namespace (for the
      admin account, `Admin`), a colon and `UiSmokeRenamed`.

13. **Copy Markup.**
    - Right-click the tile and choose **Copy > Markup**.
    - Paste the clipboard into any text field.
    - **Verify:** the pasted text is `#{x.` + the Grok name +
      `."uiSmokeRenamed"}`.

14. **Copy URL.**
    - Right-click the tile and choose **Copy > URL**.
    - Paste the clipboard into any text field.
    - **Verify:** the pasted text is the server address, then `/p/`,
      the project's namespace (the part of the Grok name from step 12
      before the colon) and `.UiSmokeRenamed`.

15. **Add to favorites.**
    - Right-click the tile and choose **Add to favorites**.
    - Click the tile.
    - **Verify:** the star next to the project name in the **Context
      Panel** header is filled.
    - Go to **Browse > My stuff > Favorites**.
    - **Verify:** `uiSmokeRenamed` is listed.
    - Right-click `uiSmokeRenamed` and choose **Remove from favorites**.

16. **Save as Zip.**
    - Go to **Browse > Dashboards**.
    - Right-click the `uiSmokeRenamed` tile and choose **Save as Zip**.
    - **Verify:** the browser downloads `uiSmokeRenamed.zip`.

17. **Reopen.**
    - Right-click the left sidebar and select **Close All**.
    - In **Browse > Dashboards**, double-click the `uiSmokeRenamed`
      tile.
    - **Verify:** the `demog` view opens.
    - **Verify:** the status bar shows **Rows: 5,850**.

18. **Change the description.**
    - Click **SAVE** on the ribbon.
    - **Verify:** the **Description** field holds *Demog snapshot for the
      UI smoke test*.
    - Replace it with `Updated description`.
    - Leave **Save original project** selected.
    - Click **OK**.
    - Right-click the left sidebar and select **Close All**.
    - In **Browse > Dashboards**, type `uiSmokeRenamed` into the search
      box and click the refresh icon.
    - **Verify:** the `uiSmokeRenamed` tile shows *Updated description*.

19. **Delete.**
    - Right-click the left sidebar and select **Close All**.
    - In **Browse > Dashboards**, right-click the `uiSmokeRenamed` tile
      and choose **Delete Project**.
    - **Verify:** the dialog *Are you sure? Delete project
      "uiSmokeRenamed"?* opens.
    - Click **DELETE**.
    - Wait until the dialog closes.
    - Click the refresh icon next to the search box.
    - **Verify:** the `uiSmokeRenamed` tile is gone.

## Cleanup

- In **Browse > Dashboards**, right-click the `uiSmokeOther` tile and
  choose **Delete Project**. Click **DELETE** and wait until the dialog
  closes. Do the same for `uiSmokeRenamed` (or `uiSmoke`) if the run
  stopped before step 19.

## Expected results

- The project is created, shared, renamed, reopened and deleted
  entirely through the UI.
- Every right-click menu item used in the test works as described.
- After deletion the project no longer appears in **Dashboards**.
- The description from the Save dialog is shown on the tile and in the
  **Context Panel**, and can be changed by saving again.
- A tag added in **Context Panel > Details** is kept, and `#tag` in the
  **Dashboards** search finds the tagged project only.

## Automation notes

- The description placeholder *Description*, the pencil tooltip, the
  **Add tag** field and the `#tag` chips are read from
  `project_meta.dart`; the tile shows the description in
  `grok-gallery-grid-item-desc`. That Esc leaves the tag editor and
  saves the tags is read from the same file.
