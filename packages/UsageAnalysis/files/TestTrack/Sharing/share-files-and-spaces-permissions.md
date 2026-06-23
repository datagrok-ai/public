
# Sharing & Permissions — File shares & Spaces

## Setup

- `Owner` — the tester's own account, logged into the primary browser session
- `Recipient` — a second, **non-admin** account
- **Files** — a **user file share** (e.g. My Files)
- **Spaces** — a **Space the tester owns** (e.g. a newly created space: Browse > Spaces; right-click the Spaces item and select Create Space...; Enter a unique name and click OK) containing a dataset.

## Scenarios

### Block A — Files are shared via the file-share connection (two-actor)

1. As **owner**, share the **file-share connection** behind the tester's share with
   `Recipient` at **View and use**.
   * Expected result: the connection is shared; **ListFiles** is part of view-and-use.
2. Log out from the `owner` session and log in as the `recipient`.
3. As **recipient**, browse the shared file share.
   * Expected result: the recipient can **list and open files** in the share, but
     cannot edit the connection or write files.

### Block B — File-share negatives, then revoke (two-actor)

1. Recipient attempts to **edit the connection / delete the share / re-share it**.
   * Expected result: all **denied** — view-and-use grants browse/read only.
2. Log out from the `recipient` session and log in as the `owner`.
3. As **owner**, open the Share dialog and click the **×** next to the recipient's grant; click **OK**.
   * Expected result: the recipient grant is removed; the Sharing pane shows
     owner-only again.

### Block C — Share a Space (namespace project) cascades to its contents (two-actor)

1. As **owner**, set the **Space** as the current object and open its **Sharing**
   surface; click **SHARE...** (a Space is a Project, so it uses the standard Share
   dialog).
   * Expected result *(confirm live)*: the **Share <space>** dialog opens with the
     **View and use / Full access** dropdown and the cascade notice **"This entity will
     also be shared to view and use:"** listing the space's contained entities.
2. As **owner**, add `recipient` name and click **OK**.
   * Expected result: the space's grant + cascaded view-and-use on its contents are
     created.
3. Log out from the `owner` session and log in as the `recipient`.
4. As **recipient**, open **Spaces** / **Shared with me**.
   * Expected result: the shared space appears for the recipient and its contained
     entities are **viewable/usable** (tables load, queries run, layouts apply) — the
     namespace-level project cascade. (Same mechanism as the Project/Dashboard case.)

### Block D — Space negatives, then revoke (two-actor)

While the recipient holds **View and use** on the space (before revoke):
1. As **recipient**, attempt to **edit** the space or its contents, **delete**, and
   **re-share**.
   * Expected result: all **denied** — view-and-use grants no Edit / Delete / Share on
     the space or its members.
2. Log out from the `recipient` session and log in as the `owner`.
3. As **owner**, **revoke** the space share last.
   * Expected result: the space and its cascaded contents disappear from the
     recipient's browse.
