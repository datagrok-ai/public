---
feature: projects
sub_features_covered:
  - projects.upload
  - projects.api.save
  - projects.shell.open
  - projects.shell.share-via-context-menu
  - projects.api.search
  - projects.api.delete
target_layer: playwright
coverage_type: smoke
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

Single short UI-driven flow that exercises the right-click context menu,
Share / Delete / Rename dialogs, and Browse > Dashboards gallery. This is
the **only** scenario in the Projects chain that drives those UI surfaces;
every other scenario delegates UI coverage here per chain rev 4
`ui_coverage_plan.smoke_scenario` and substitutes JS API for the same
flows.

Target runtime: ~5 minutes. One project, one source (`demog.csv` from
`System:DemoFiles`). Recipient placeholder
`<RECIPIENT_GROUP_USERNAME_TBD>` for the Share step (resolved at
Automator stage to a per-environment group account).

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
   Step 9-11 (terminal Delete via UI). NO handoff to `deleting.md` —
   self-cleaning by design.

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

- **Origin: chain rev 3 ui-smoke consolidation (Olena's Option ε).**
  This scenario was authored per Phase A of the
  `projects migrate --force` cycle to address the F-UI-COVERAGE-01
  SCOPE_REDUCTION. Chain rev 3 emitted only DELEGATION (single
  smoke_scenario = upload-project.md, smoke_covers = 5 flows); F's
  audit found 12 of 16 extracted UI flows uncovered. Olena's
  decision: author a dedicated consolidated ui-smoke scenario,
  redirect 4 existing scenarios' `ui_coverage_delegated_to` here.
  Chain re-emit (rev 4) will register this scenario as the new
  `smoke_scenario` and update `delegated_scenarios` accordingly
  (Phase B).
- **Pyramid layer: ui-smoke (Rule 1).** Single short UI flow over
  one source; not matrix; not bug-focused; not manual. Rule 1
  applies. Per chain rev 4 plan, this scenario is THE smoke owner
  for the Projects section.
- **target_layer: playwright.** UI driving is the entire point of
  this scenario. JS API substitution is explicitly forbidden — if
  any step is replaced by `grok.dapi.*` or
  `grok.shell.project = ...`, the scenario loses its purpose.
- **UI coverage owned (chain rev 4 ui_coverage_responsibility, 20
  flows).** This scenario witnesses the entire right-click context
  menu (11 items per `references/projects.md` ## Context Menu
  Items table) plus the 4 verb-form dialog surfaces (Save, Share,
  Delete, Rename) plus auxiliary surfaces (auto-share dialog
  dismiss, Browse Dashboards tile visibility / search). 4 sibling
  scenarios delegate to this one:
  `upload-project.md`, `share-project.md`, `opening.md`,
  `deleting.md` — each updates `ui_coverage_delegated_to:
  projects-ui-smoke.md` per Phase A revisions.
- **Recipient placeholder convention.** Per Phase A authoring
  conventions, real recipient identities (`Olena Ahadzhanian`,
  `Test permission group`) are NOT in this `.md`. Placeholders
  (`<RECIPIENT_USERNAME_TBD>`, `<RECIPIENT_GROUP_USERNAME_TBD>`)
  are resolved at Automator stage from per-environment helper
  config.
- **Self-cleaning.** Steps 9-11 delete the project. This scenario
  is NOT in `deleting.md.depends_on` — its cleanup is internal.
  Does not produce fixtures consumed by other scenarios.
- **Reference selectors.** Selectors required by this scenario
  (`pcmdShareProject`, `pcmdDeleteProject`, `pcmdRename`, Share
  dialog widget, Delete confirmation dialog, Rename dialog,
  Browse > Dashboards per-tile selector) — per
  PROJECTS-SPLIT-COMPLETION-PLAN Wave 2D Step 1, these MUST be
  registered in `.claude/skills/grok-browser/references/projects.md`
  before Automator stage. Currently the Plan flags those as B14
  reference-additions awaiting approval. Phase A authoring of THIS
  `.md` does NOT depend on the selectors being registered yet —
  the `.md` describes UI flows by intent. Spec generation
  (Automator, Phase C) is the consumer of those selectors.
