---
feature: sharing
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [cp-share-via-context-menu, cp-context-panel-share-button, cp-advanced-editor-matrix, cp-shared-with-me-browse-node]
realizes: []
realized_as:
  - share-layout-permissions-spec.ts
related_bugs: []
---


# Sharing & Permissions — Layout

Verifies the Share dialog and the Advanced editor permissions matrix for a saved view
layout: opening the Sharing pane from the Context Panel, sharing at "View and use", and
confirming that a recipient with that grant can see and apply the layout to a table —
but cannot rename, delete, or re-share it — until the owner revokes access.

## Setup

- `Owner` — the tester's own account, logged into the primary browser session
- `Recipient` — a second, non-admin account
- `Custom group` — a dedicated group containing the recipient, for the group-sharing block
- Entity under test — a ViewLayout owned by the owner. The owner must NOT have previously shared this layout (start from a clean grant list — owner only). Create a layout for the demog dataset by adding a scatterplot viewer.

> Two independent browser sessions (owner + non-admin recipient) are required for Blocks
> D, E, F — no single-session automation is possible.
> Source: `core/client/xamgle/lib/src/commands/file/share_dataset.dart#L7`.

## Scenarios

### Block A — Open the Share dialog from the Sharing pane (owner)

1. As **owner**, set the layout as the current object (go to Browse > Platform > Layouts and click a layout so the Context Panel shows it). Expand the **Sharing** pane in the Context Panel.
   * Expected result: the Sharing pane lists the current grant — **"`owner` can view and use"** and **"You are the owner"** — and shows a **SHARE...** button.
2. Click **SHARE...**.
   * Expected result: the **Share `layout name`** dialog opens with: a **"User, group, or email"** autocomplete input, a level dropdown defaulting to **View and use**, the existing owner grant row (with a **×**), an **Advanced editor...** link, and **OK** / **CANCEL**.

### Block B — Share-dialog controls (owner)

1. Type a few letters of the recipient into the **"User, group, or email"** input.
   * Expected result: an autocomplete list suggests matching users/groups; selecting one adds it as a pill with its own level dropdown and a **×** to remove.
2. Observe the dependent-entity notice and the notification controls.
   * Expected result: because a layout has no separately-shareable member entities, the dialog shows **no** "This entity will also be shared..." cascade notice. When adding a new recipient, a **Send notifications** checkbox and a message box ("Type in message here") are available.
3. Click **CANCEL**.
   * Expected result: the dialog closes; no grant is changed (Sharing pane still shows owner-only).

### Block C — Advanced editor per-permission matrix (owner)

1. Open the Share dialog again and click **Advanced editor...**.
   * Expected result: a permissions view opens (route `/permissions/<id>`) showing a **Group × Object** matrix. Columns are grouped under **Common** — **View**, **Edit**, **Delete**, **Share**. A ViewLayout uses the **standard set only** — there is **no** extra entity-specific use-permission column (applying a layout is a View-level use). The owner's existing row shows **View** checked.
2. Inspect the bottom **"Type in user, role or group to add..."** row.
   * Expected result: a new group/user can be added and individual permission checkboxes toggled independently, with a **SAVE** action. Close without saving.

### Block D — Share "View and use" to the recipient; recipient gains read + use (two-actor)

1. As **owner**, open the Share dialog, add the `Recipient` name (or the **`Custom Group`**), leave the level at **View and use**, leave **Send notifications** as desired, and click **OK**.
   * Expected result: the dialog closes with no error; the Sharing pane now lists the recipient with **"can view and use"** in addition to the owner.
2. Log out from the `owner` session and log in as the `recipient`. In the **recipient** browser session, open **My stuff > Shared with me**.
   * Expected result: the shared layout appears under **Shared with me** for the recipient, and shows up in the recipient's **Layouts**.
3. As **recipient**, open the demog table. Go to View > Layout > Open Gallery, locate and apply the shared layout.
   * Expected result: the layout applies successfully — its columns/viewers are restored on the recipient's table.

### Block E — Negative: with "View and use", recipient cannot edit, delete, or re-share (two-actor)

1. As **recipient**, attempt to **rename** the layout and try to **save**.
   * Expected result: editing is **not permitted**.
2. As **recipient**, attempt to **delete** the layout (Delete action / context menu).
   * Expected result: delete is **not permitted**; the layout still exists.
3. As **recipient**, open the layout's **Sharing** pane and attempt to **re-share** it (the way you would as owner).
   * Expected result: there is **no usable SHARE... affordance** (absent / disabled), or any attempt to grant is **rejected** — "View and use" does not include the **Share** permission. The grant list is unchanged.

### Block F — Revoke access from the recipient; access lost (two-actor) — do this LAST

1. Log out from the `recipient` session and log in as the `owner`.
2. As **owner**, open the Share dialog and click the **×** next to the recipient's grant; click **OK**.
   * Expected result: the recipient grant is removed; the Sharing pane shows owner-only again.
3. Log out from the `owner` session and log in as the `recipient`.
4. In the **recipient** session, refresh and look under **Shared with me**.
   * Expected result: the layout no longer appears under Shared with me / Layouts and the recipient can no longer view or apply it (access revoked).
