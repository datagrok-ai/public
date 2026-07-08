---
feature: projects
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: [upload-save-reopen-golden, share_with_recipient_open, rename_project]
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

Single short UI-driven flow that exercises the right-click context
menu, the Share / Delete / Rename dialogs, and the Browse >
Dashboards gallery. This is the only scenario in the Projects section
that drives those UI surfaces directly — every other scenario
substitutes the JS API for the same operations and relies on this
scenario for their UI coverage.

Target runtime: about 5 minutes. Uses one project with one source
(`demog.csv` from `System:DemoFiles`). The recipient for the Share
step is a placeholder, `<RECIPIENT_GROUP_USERNAME_TBD>`, resolved to
a real per-environment group account at automation time.

## Setup

1. Authenticate to Datagrok as the test user. The session must be UI
   (Playwright); JS API substitution is NOT allowed in this scenario —
   that's the whole point.
2. Project name: pick a unique value like
   `ui-smoke-${Date.now()}` to avoid collisions with parallel runs.
3. Recipient placeholder: `<RECIPIENT_GROUP_USERNAME_TBD>` (resolved
   at Automator stage). Do NOT bind to a literal user/group in this
   `.md`.
4. Cleanup contract: this scenario deletes its own project at
   Step 9-11 (terminal Delete via UI) — self-cleaning by design.

## Scenarios

### Main flow — open, save, share, reopen, rename, delete

1. **Open table from File share (UI).** Navigate
   `Browse` > `Files` > `Demo` and right-click `demog.csv` and select Open.
    Verify the table opens in the workspace
   (`grok.shell.tables.length > 0`).
2. **Save Project with Data Sync ON (UI).** Trigger
   Click the Save ribbon button — NOT the
   Ctrl+S keybind, per `feedback_no_ctrlS_for_layouts`. In the Save
   Project dialog: enter the unique project name from Setup, ensure
   the **Data Sync** toggle is **ON**, click **OK**. Verify the save
   completes (project appears in `Browse` > `Dashboards` OR the
   underlying `POST /projects` succeeds via network log).
3. **Cancel auto-opened Share dialog (UI).** After save, the platform
   opens a Share dialog automatically. Click **Cancel** (or close the
   dialog via the X button). Verify the dialog disappears and no
   share permissions were granted (no entries on the project's
   Sharing tab).
4. **Browse > Dashboards — assert tile visible (UI).** Navigate to
   `Browse` > `Dashboards`. Locate the saved project's tile in the
   gallery (search by name if needed). Assert the tile is visible
   and clickable.
5. **Right-click tile → Share menu item (UI — covers `pcmdShareProject`).**
   Right-click the project tile. From the context menu, click
   **Share**. Verify the Share dialog opens.
6. **Share dialog → recipient → access level → OK (UI — covers Share
   dialog modal).** In the Share dialog:
   - Type `<RECIPIENT_GROUP_USERNAME_TBD>` into the recipient
     picker (covers `share-dialog-recipients`).
   - Pick **View-and-Use** access level from the permissions editor
     dropdown (covers `share-dialog-permissions-editor`).
   - Click **OK**. Verify the dialog closes and the recipient appears
     on the project's Sharing tab (Context Panel verification).
7. **Reopen project — double-click tile (UI — covers
   `project-double-click-open`).** Close all current views
   (`grok.shell.closeAll()` is acceptable here as a state-reset). On
   `Browse` > `Dashboards`, double-click the project tile. Verify
   the project reopens and `grok.shell.tables.length > 0`.
8. **Verify table loaded.** Assert the `demog` table is loaded with
   non-zero row count.
9. **Right-click tile → Delete menu item (UI — covers
   `pcmdDeleteProject`).** Navigate back to `Browse` > `Dashboards`.
   Right-click the project tile. From the context menu, click
   **Delete**. Verify the DELETE confirmation dialog opens.
10. **DELETE confirmation → DELETE button (UI — covers Delete dialog).**
    In the confirmation dialog, click **DELETE** (the destructive-
    action button, not Cancel). Verify the dialog closes.
11. **Assert tile gone from gallery (UI).** Re-search Browse >
    Dashboards for the project name. Assert the tile is no longer
    present.

### Secondary coverage — context-menu items (smoke depth)

The following context-menu items are touched as part of the smoke run
(not the full main flow, but reachable from the right-click menu and
asserted to open / dispatch):

12. **Right-click tile → Rename (UI — covers `pcmdRename`,
    `rename-dialog`).** Before the Delete in Step 9, optionally
    perform a Rename: right-click tile > **Rename**, enter a new
    name (e.g. `<original-name>_renamed`), click **OK**. Verify the
    tile reflects the new name. (This step can be ordered before
    Step 9 in the spec implementation.)
13. **Right-click tile → Save as Zip (UI — covers `pcmdSaveAsZip`).**
    From the context menu, click **Save as Zip**. Verify the
    download dispatches (browser download triggered or a file blob
    is produced via Playwright download handler).
14. **Right-click tile → Copy submenu (UI — covers `pcmdCopy`,
    `pcmdCopyId`, `pcmdCopyGrokName`, `pcmdCopyMarkup`,
    `pcmdCopyUrl`).** Hover over the **Copy** submenu. Click each
    of the 4 sub-items in turn (Copy ID, Copy Grok Name, Copy
    Markup, Copy URL). For each, verify the clipboard receives the
    expected payload (use `navigator.clipboard.readText()` after
    each click).
15. **Right-click tile → Add to favorites (UI — covers
    `pcmdAddToFavorites`).** Click **Add to favorites**. Verify
    the tile appears under `Browse` > `Favorites` after the click
    (or `grok.dapi.favorites` reflects the addition).

### Expected results

After completing the main flow (Steps 1-11):
- The project went through the full UI lifecycle: created, saved,
  shared, reopened, deleted — all via UI driving, no JS API
  substitution.
- The Sharing tab confirms the recipient grant.
- The Browse > Dashboards gallery confirms the deletion (tile gone).

After completing the secondary coverage (Steps 12-15):
- Each context-menu item dispatches its expected action (rename
  applied, download triggered, clipboard populated, favorites
  updated).

## Notes

- **UI driving is the entire point.** JS API substitution is
  explicitly out of scope here — if any step were replaced by
  `grok.dapi.*` or `grok.shell.project = ...`, the scenario would
  lose its purpose.
- **UI coverage owned here.** This scenario witnesses the entire
  right-click context menu (11 items) plus the four dialog surfaces
  (Save, Share, Delete, Rename) plus auxiliary surfaces (auto-share
  dialog dismiss, Browse > Dashboards tile visibility / search). It is
  the canonical UI-coverage owner that the section's other (JS-API)
  scenarios delegate to.
- **Recipient placeholder convention.** Real recipient identities are
  intentionally not hard-coded in this `.md`. The placeholders
  (`<RECIPIENT_USERNAME_TBD>`, `<RECIPIENT_GROUP_USERNAME_TBD>`) are
  resolved to real per-environment accounts at automation time.
- **Self-cleaning.** Steps 9-11 delete the project. This scenario
  doesn't produce fixtures consumed by other scenarios.
