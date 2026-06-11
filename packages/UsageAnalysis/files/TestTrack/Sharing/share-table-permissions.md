---
feature: sharing
sub_features_covered:
  - sharing.context-panel-pane
  - sharing.share-dialog
  - sharing.permissions-editor
  - sharing.advanced-editor
  - sharing.browse-shared-with-me
  - sharing.entity-types.table
  - sharing.notification
target_layer: manual-only
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/sharing/share-table-permissions.md
migration_date: '2026-06-11'
source_text_fixes:
  - trailing-space-in-share-dialog-title-block-a-step-2
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs: []
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-11-sharing-migrate-02
    timestamp: '2026-06-11T10:00:00Z'
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-11-sharing-migrate-02
    timestamp: '2026-06-11T10:30:00Z'
    failure_keys: []
    review_round: 1
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
---

# Sharing & Permissions — Table (stored TableInfo / dataframe)

## Setup

- `Owner` — the tester's own account, logged into the primary browser session
- `Recipient` — a second, **non-admin** account
- `Custom group` — a dedicated group containing the recipient, for the group-sharing block
- **Entity under test** — a **stored TableInfo owned by the owner** (e.g. a `demog` table saved as a dashboard)

> **Why manual-only:** Blocks D, E, and F require two independent browser sessions (owner + non-admin recipient). The atlas `manual_only[]` entry for `sharing.entity-types.table` confirms this constraint — two-actor authentication cannot be reduced to deterministic single-session UI automation.
>
> **Entity-specific notes:** TableInfo has no shareable dependencies, so Block B shows **no** cascade notice. The Advanced editor (Block C) includes a `ReadData` entity-specific column in addition to the common View/Edit/Delete/Share columns.

## Scenarios

### Block A — Open the Share dialog from the Sharing pane (owner)

1. As **owner**, set the table as the current object (Browse > My stuff, expand the created earlier dashboard, click the demog table so the Context Panel shows it). Expand the **Sharing** pane in the Context Panel.
   * Expected result: the Sharing pane lists the current grant — **"`owner` can view and use"** and **"You are the owner"** — and shows a **SHARE...** button.
2. Click **SHARE...**.
   * Expected result: the **Share demog** dialog opens with: a **"User, group, or email"** autocomplete input, a level dropdown defaulting to **View and use**, the existing owner grant row (with a **×**), an **Advanced editor...** link, and **OK** / **CANCEL**.

### Block B — Share-dialog controls (owner)

1. Type a few letters of the recipient into the **"User, group, or email"** input.
   * Expected result: an autocomplete list suggests matching users/groups; selecting one adds it as a pill with its own level dropdown and a **×** to remove.
2. Observe the dependent-entity notice and the notification controls.
   * Expected result: because a standalone table has no shareable dependencies, the dialog shows **no** "This entity will also be shared..." cascade notice. When adding a new recipient, a **Send notifications** checkbox and a message box ("Type in message here") are available.
3. Click **CANCEL**.
   * Expected result: the dialog closes; no grant is changed (Sharing pane still shows owner-only).

### Block C — Advanced editor per-permission matrix (owner)

1. Open the Share dialog again and click **Advanced editor...**.
   * Expected result: a permissions view opens (route `/permissions/<id>`) showing a **Group × Object** matrix. Columns are grouped under **Common** — **View**, **Edit**, **Delete**, **Share** — plus the entity-specific **ReadData** column for a TableInfo. The owner's existing rows show **View** (and **ReadData**) checked.
2. Inspect the bottom **"Type in user, role or group to add..."** row.
   * Expected result: a new group/user can be added and individual permission checkboxes toggled independently, with a **SAVE** action.
3. Close without saving.
   * Expected result: the matrix is dismissed with no permission changes applied.

### Block D — Share "View and use" to the recipient; recipient gains read + ReadData (two-actor)

1. As **owner**, open the Share dialog, add the `Recipient` name (or the **`Custom Group`**), leave the level at **View and use**, leave **Send notifications** as desired, and click **OK**.
   * Expected result: the dialog closes with no error; the Sharing pane now lists the recipient with **"can view and use"** in addition to the owner.
2. Log out from the `owner` session and log in as the `recipient`. In the **recipient** browser session, open **My stuff > Shared with me**.
   * Expected result: the shared table appears under **Shared with me** for the recipient.
3. As **recipient**, **open the table** and let its data load.
   * Expected result: the table opens and **its rows load** successfully.

### Block E — Negative: with "View and use", recipient cannot delete or re-share (two-actor)

1. As **recipient**, attempt to **delete** the table (Delete action / context menu).
   * Expected result: delete is **not permitted**; the table still exists.
2. As **recipient**, open the table's **Sharing** pane and attempt to **re-share** it (the way you would as owner).
   * Expected result: there is **no usable SHARE... affordance** (absent / disabled); the grant list is unchanged.

### Block F — Revoke access from the recipient; access lost (two-actor) — do this LAST

1. Log out from the `recipient` session and log in as the `owner`.
2. As **owner**, open the Share dialog and click the **×** next to the recipient's grant; click **OK** (alternatively, use the Advanced editor to uncheck the grants and SAVE).
   * Expected result: the recipient grant is removed; the Sharing pane shows owner-only again.
3. Log out from the `owner` session and log in as the `recipient`.
4. In the **recipient** session, refresh and look under **Shared with me** / try to open the table.
   * Expected result: the table no longer appears under Shared with me and the recipient can no longer open it or load its data (access revoked).

## Notes

- Atlas source: `core/client/xamgle/lib/src/commands/file/share_dataset.dart#L7` (Share dialog command registration); `core/client/xamgle/lib/src/features/users/permissions_browser.dart#L403` (Advanced editor / PermissionsView); `core/client/xamgle/lib/src/features/browse_panel/browse_panel_tree.dart#L96` (Browse > Shared with me node).
- TableInfo is a standalone entity (`source_classes.table.external_deps: []`); no transitive dependency co-sharing required.
- The `ReadData` permission column in the Advanced editor is TableInfo-specific (grants read access to the table's row data). Common columns are View/Edit/Delete/Share.
- Related help: `public/help/govern/access-control/access-control.md` §Permissions covers the TableInfo permission model.

---
{
  "order": 5,
  "datasets": []
}
