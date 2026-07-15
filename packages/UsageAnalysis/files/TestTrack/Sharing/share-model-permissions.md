---
feature: sharing
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: []
realizes: []
realized_as:
  - share-model-permissions-spec.ts
related_bugs: []
---


# Sharing & Permissions — Model (predictive model)

Verifies the Share dialog and the Advanced editor permissions matrix for a predictive
model: opening the Sharing pane from the Context Panel, sharing at "View and use", and
confirming that a recipient with that grant can open and apply the model to a dataset —
but cannot edit, delete, or re-share it — until the owner revokes access. Also confirms
a standalone model shares only itself, with no cascade to other entities.

## Setup

- `Owner` — the tester's own account, logged into the primary browser session
- `Recipient` — a second, **non-admin** account
- `Custom group` — a dedicated group containing the recipient, for the group-sharing block
- **Entity under test** — a **predictive Model created by the owner** (e.g. a model like
  Predict SEX by HEIGHT, WEIGHT for the demog dataset).

## Scenarios

### Block A — Open the Share dialog from the Sharing pane (owner)

1. As **owner**, set the model as the current object (Browse > Platform > Predictive Models and click the model, so the Context Panel shows it). Expand the **Sharing** pane in the Context Panel.
   * Expected result: the Sharing pane lists the current grant — **"`owner` can view
     and use"** and **"You are the owner"** — and shows a **SHARE...** button.
2. Click **SHARE...**.
   * Expected result: the **Share `model name`** dialog opens with: a **"User, group, or
     email"** autocomplete input, a level dropdown defaulting to **View and use**, the
     existing owner grant row (with a **×**), an **Advanced editor...** link, and
     **OK** / **CANCEL**.

### Block B — Share-dialog controls (owner)

1. Type a few letters of the recipient into the **"User, group, or email"** input.
   * Expected result: an autocomplete list suggests matching users/groups; selecting
     one adds it as a pill with its own level dropdown and a **×** to remove.
2. Observe the dependent-entity notice and the notification controls.
   * Expected result: because a standalone model shares only itself, **no cascade
     notice** appears (no "This entity will also be shared..." line). When adding a new
     recipient, a **Send notifications** checkbox and a message box ("Type in message
     here") are available.
3. Click **CANCEL**.
   * Expected result: the dialog closes; no grant is changed (Sharing pane still shows
     owner-only).

### Block C — Advanced editor per-permission matrix (owner)

1. Open the Share dialog again and click **Advanced editor...**.
   * Expected result: a permissions view opens (route `/permissions/<id>`) showing a
     **Group × Object** matrix. Columns are grouped under **Common** — **View**,
     **Edit**, **Delete**, **Share** — the standard set for a Model. There is **no
     separate Execute column** for a Model unless present; **applying the model maps to
     the View-level use**. The owner's existing rows show **View** checked.
2. Inspect the bottom **"Type in user, role or group to add..."** row.
   * Expected result: a new group/user can be added and individual permission
     checkboxes toggled independently, with a **SAVE** action. Close without saving.

### Block D — Share "View and use" to the recipient; recipient gains read + use (two-actor)

1. As **owner**, open the Share dialog, add the `Recipient` name (or the **`Custom Group`**), leave the level at **View and use**, leave **Send notifications** as desired, and click **OK**.
   * Expected result: the dialog closes with no error; the Sharing pane now lists the recipient with **"can view and use"** in addition to the owner.
2. Log out from the `owner` session and log in as the `recipient`. In the **recipient** browser session, open **My stuff > Shared with me** (or search for the model in Browse > Platform > Predictive Models view).
   * Expected result: recipient has access to the shared model.
3. As **recipient**, open the demog dataset, right-click the model and **apply** it to a dataset.
   * Expected result: the model opens and can be **applied** to a dataset.

### Block E — Cascade: no outbound cascade for a standalone model (two-actor)

1. As **recipient**, confirm only the model itself was shared — **no additional
   dependent entity** is listed as auto-shared.
   * Expected result: **no outbound cascade for a standalone model (model shares only
     itself).** A model may depend on its training pipeline, but for sharing a
     standalone model shares only itself; a cascade notice appears only if it bundles
     dependencies.

### Block F — Negative: with "View and use", recipient cannot edit, delete, or re-share (two-actor)

1. As **recipient**, attempt to **edit** the model — rename / add description to it and try to **save**.
   * Expected result: editing is **not permitted**; the model is unchanged.
2. As **recipient**, attempt to **delete** the model (Delete action / context menu).
   * Expected result: delete is **not permitted**; the model still exists.
3. As **recipient**, open the model's **Sharing** pane and attempt to **re-share** it (the way you would as owner).
   * Expected result: there is **no usable SHARE... affordance** (absent / disabled), or any attempt to grant is **rejected**. The grant list is unchanged.

### Block G — Revoke access from the recipient; access lost (two-actor) — do this LAST

1. Log out from the `recipient` session and log in as the `owner`.
2. As **owner**, open the Share dialog and click the **×** next to the recipient's grant; click **OK**.
   * Expected result: the recipient grant is removed; the Sharing pane shows
     owner-only again.
3. Log out from the `owner` session and log in as the `recipient`.
4. In the **recipient** session, refresh and look under **Shared with me** / try to locate the model.
   * Expected result: the model no longer appears under Shared with me and the recipient can no longer open or apply it (access revoked).
