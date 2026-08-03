
# Sharing & Permissions — Script

## Setup

- `Owner` — the tester's own account, logged into the primary browser session
- `Recipient` — a second, **non-admin** account
- `Custom group` — a dedicated group containing the recipient, for the group-sharing block
- **Entity under test** — a **Script owned by the owner**, e.g. a JS script, created from Browse > Platform > Functions > Scripts, click New > JS Script, and insert:

```
//name: Template
//description: Hello world script
//language: javascript

grok.shell.info('hello world')
```

Click SAVE.

## Scenarios

### Block A — Open the Share dialog from the Sharing pane (owner)

1. As **owner**, set the script as the current object (Browse > Platform > Functions > Scripts and click the script, so the Context Panel shows it). Expand the **Sharing** pane in the Context Panel.
   * Expected result: the Sharing pane lists the current grant — **"`owner` can view and use"** and **"You are the owner"** — and shows a **SHARE...** button.
2. Click **SHARE...**.
   * Expected result: the **Share connection name** dialog opens with: a **"User, group, or email"** autocomplete input, a level dropdown defaulting to **View and use**, the existing owner grant row (with a **×**), an **Advanced editor...** link, and **OK** / **CANCEL**.

### Block B — Share-dialog controls (owner)

1. Type a few letters of the recipient into the **"User, group, or email"** input.
   * Expected result: an autocomplete list suggests matching users/groups; selecting one adds it as a pill with its own level dropdown and a **×** to remove.
2. Observe the dependent-entity notice and the notification controls.
   * Expected result: a script with no shareable dependencies shows no such notice. When adding a new recipient, a **Send notifications** checkbox and a message box ("Type in message here") are available.
3. Click **CANCEL**.
   * Expected result: the dialog closes; no grant is changed (Sharing pane still shows owner-only).

### Block C — Advanced editor per-permission matrix (owner)

1. Open the Share dialog again and click **Advanced editor...**.
   * Expected result: a permissions view opens (route `/permissions/<id>`) showing a **Group × Object** matrix. Columns are grouped under **Common** — **View**, **Edit**, **Delete**, **Share** — plus the entity-specific **Execute** column for a Script. The owner's existing rows show **View** (and **Execute**) checked.
2. Inspect the bottom **"Type in user, role or group to add..."** row.
   * Expected result: a new group/user can be added and individual permission checkboxes toggled independently, with a **SAVE** action.
3. Close without saving.

### Block D — Share "View and use" to the recipient; recipient gains read + execute (two-actor)

1. As **owner**, open the Share dialog, add `Recipient` name (or the **`Custom Group`**), leave the level at **View and use**, leave **Send notifications** as desired, and click **OK**.
   * Expected result: the dialog closes with no error; the Sharing pane now lists the recipient with **"can view and use"** in addition to the owner.
2. Log out from the `owner` session and log in as the `recipient`. In the **recipient** browser session, open **My stuff → Shared with me** (or search for the script).
   * Expected result: the shared script appears under **Shared with me** for the recipient.
3. As **recipient**, open the script and **run/execute** it (context menu > Run).
   * Expected result: the script runs successfully and returns its output.

### Block E — Negative: with "View and use", recipient cannot edit, delete, or re-share (two-actor)

1. As **recipient**, attempt to **edit** the script — change its source code and try to **save**.
   * Expected result: editing is **not permitted**; the script source is unchanged.
2. As **recipient**, attempt to **delete** the script (Delete action / context menu).
   * Expected result: delete is **not permitted**; the script still exists.
3. As **recipient**, open the script's **Sharing** pane and attempt to **re-share** it (the way you would as owner).
   * Expected result: there is **no usable SHARE... affordance**. The grant list is unchanged.

### Block F — Revoke access from the recipient; access lost (two-actor) — do this LAST

1. Log out from the `recipient` session and log in as the `owner`.
2. As **owner**, open the Share dialog and click the **×** next to the recipient's grant; click **OK** (alternatively, use the Advanced editor to uncheck the grants and SAVE).
   * Expected result: the recipient grant is removed; the Sharing pane shows owner-only again.
3. Log out from the `owner` session and log in as the `recipient`.
4. In the **recipient** session, refresh and look under **Shared with me** / try to run the script.
   * Expected result: the script no longer appears under Shared with me and the recipient can no longer run it (access revoked).
