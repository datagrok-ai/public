---
feature: sharing
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p0
realizes_atlas: [sharing.cp.cp-share-via-context-menu, sharing.cp.cp-context-panel-share-button, sharing.cp.cp-shared-with-me-browse-node]
realizes: [sharing.share-dialog, sharing.browse-shared-with-me, views.queries, views.space]
related_bugs:
  - GROK-19235
  - GROK-19401
  - GROK-18525
  - GROK-18210
  - GROK-21142
---

# Sharing — what the other user gets, and losing it

A query shared once through the Share dialog, seen from the other side: the second user finds it
under **Shared with me** and runs it, cannot pass it on, cannot put a shared query into a space
shared with them for viewing only, and loses the query when you take the share back. The same path
holds for every kind of entity; which privileges each kind gets is an API test.

## Setup

1. A second user account other than yours exists on the server (for example `bddsecond`). It has
   no notifications.
2. A query `ShareRecipient1` on the `Datagrok` connection (`System:Datagrok`), with the text
   `select 1 as x`, owned by you. The `Datagrok` connection is not shared with the second user for
   this test.
3. A query `ShareRecipient2` on the `Datagrok` connection, with the text `select 2 as y`, owned by
   you and shared with the second user at **View and use**.
4. A space `ShareSpace1`, owned by you, empty, and shared with the second user at **View and use**.

## Scenarios

### 1. You share the query to view and use (GROK-19235)

1. Go to **Browse > Platform > Functions > Queries** and type `ShareRecipient1` into
   **Search queries by name or by #tags**.
2. Right-click the `ShareRecipient1` card and choose **Share...**.
3. In **User, group, or email**, type the first letters of the second user's login and pick the
   user.
4. Leave the level at **View and use**. Clear **Send notifications**.
5. Click **OK**.
6. **Verify:** the dialog closes and the balloon *Shared* appears; no error.
7. Click the `ShareRecipient1` card and expand **Sharing** in the **Context Panel**.
8. **Verify:** the pane reads *You are the owner* and, under *On ShareRecipient1*, names the second
   user with *can view and use*.

### 2. The second user finds the query under Shared with me and runs it (GROK-19401, GROK-18525, GROK-18210)

1. Sign in as the second user.
2. In **Browse**, expand **My stuff > Shared with me**.
3. **Verify:** **Shared with me** lists your name.
4. Expand your name.
5. **Verify:** `ShareRecipient1` is listed under it.
6. Double-click `ShareRecipient1`.
7. **Verify:** the query runs and the `ShareRecipient1` table view opens with one column, `x`, and
   one row; no error.

### 3. The second user cannot share the query further

1. As the second user, in **Browse > My stuff > Shared with me > (your name)**, right-click
   `ShareRecipient1`.
2. **Verify:** the menu has **Share...**.
3. Choose **Share...**.
4. **Verify:** the **Share ShareRecipient1** dialog opens.
5. In **User, group, or email**, type the first letters of your name and pick yourself. Clear
   **Send notifications**.
6. Click **OK**.
7. **Verify:** the error balloon *You don't have a permission to share this object* appears and the
   dialog stays open.
8. Click **CANCEL**.

### 4. The second user cannot put a query into a space shared for viewing (GROK-21142)

1. As the second user, in **Browse > My stuff > Shared with me > (your name)**, find
   `ShareRecipient2` and `ShareSpace1`.
2. Drag `ShareRecipient2` onto `ShareSpace1`.
3. **Verify:** the **Move entity** dialog opens, reading *to ShareSpace1 project*, with a choice of
   **Move** and **Link** that reads **Link**, and the query `ShareRecipient2`.
4. Choose **Move** and click **YES**.
5. **Verify:** the error balloon *Unable to add entity ... to the project, insufficient privileges
   to share it* appears and the dialog stays open.
6. Click **CANCEL**.
7. Drag `ShareRecipient2` onto `ShareSpace1` again. In the **Move entity** dialog, leave **Link**
   and click **YES**.
8. **Verify:** an error balloon appears; the link is not made.
9. Double-click `ShareSpace1`.
10. **Verify:** the `ShareSpace1` view opens and does not list `ShareRecipient2`.

### 5. You take the share back

1. Sign in as yourself again.
2. Go to **Browse > Platform > Functions > Queries**, find `ShareRecipient1`, right-click its card
   and choose **Share...**.
3. **Verify:** the dialog lists the second user with **View and use**.
4. Hover over the second user's row and click its **×**.
5. **Verify:** the row disappears from the dialog.
6. Click **OK**.
7. **Verify:** the dialog closes and the balloon *Shared* appears.
8. Click the `ShareRecipient1` card and expand **Sharing** in the **Context Panel**.
9. **Verify:** the pane reads *You are the owner* and does not name the second user.

### 6. The second user no longer has the query

1. Sign in as the second user.
2. In **Browse**, expand **My stuff > Shared with me > (your name)**.
3. **Verify:** `ShareRecipient1` is not listed.
4. Sign in as yourself again.

## Cleanup

1. In **Browse > Platform > Functions > Queries**, delete `ShareRecipient1` and `ShareRecipient2`
   (right-click the card, **Delete**, confirm).
2. In **Browse > Spaces**, right-click `ShareSpace1`, choose **Delete Space** and confirm.

## Automation notes

- The feature is `@serial`: it switches the signed-in account on its page with `user signs in as
  the sharing user` / `user signs in as themselves again`.
- Setup items 3 and 4 are the second user's starting state, not the subject; they may be made
  through the API (the share of a space through the dialog is `spaces-share-and-edges`).
- Scenario 3, step 5: picking the running account in the user field is not a library step yet
  (there is one for the second user only).
- Scenario 4 is a known failure while GROK-21142 is open: on localhost (2026-10-09) **Move** was
  refused as written, but **Link** made the link without any message and `ShareSpace1` then listed
  `ShareRecipient2` for the second user. Remove the known-failure mark when a link by a user with
  **View and use** on the space is refused.
- Scenario 4, steps 3 and 7: the **Move**/**Link** choice is a plain select in the dialog; it opens
  on **Link** for the second user (localhost, 2026-10-12).
- That the second user may not edit, rename or delete the query is the server's answer and is an
  API test.
