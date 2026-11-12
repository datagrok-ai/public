---
feature: sharing
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
priority: p0
realizes_atlas: [sharing.cp.cp-share-via-context-menu, sharing.cp.cp-context-panel-share-button]
realizes: [sharing.share-dialog]
related_bugs:
  - GROK-20755
  - GROK-20758
---

# Sharing — the Share dialog and the Sharing pane

What you see when you are about to share something: the **Sharing** pane of the Context Panel, the
**Share** dialog it opens, and the **Resulting permissions** dialog. Nothing is shared here; the
share itself and what the other user then gets are in `sharing-recipient.md`.

## Setup

1. A second user account other than yours exists on the server (for example `bddsecond`).
2. A query `ShareDialog1` on the `Datagrok` connection (`System:Datagrok`), with the text
   `select 1 as x`, owned by you.
3. Go to **Browse > Platform > Functions > Queries** and type `ShareDialog1` into
   **Search queries by name or by #tags**.

## Scenarios

### 1. The Sharing pane names the owner and opens the dialog

1. Click the `ShareDialog1` card.
2. In the **Context Panel**, expand **Sharing**.
3. **Verify:** the pane reads *You are the owner* and does not name the second user. Below it are a
   **SHARE...** button and a calculator icon.
4. Click **SHARE...**.
5. **Verify:** the **Share ShareDialog1** dialog opens.
6. Click **CANCEL**.
7. **Verify:** the dialog closes; the **Sharing** pane still reads only *You are the owner*.

### 2. The dialog asks who and how much, and offers a message once someone is picked

1. Right-click the `ShareDialog1` card.
2. **Verify:** the menu has **Share...**.
3. Choose **Share...**.
4. **Verify:** the **Share ShareDialog1** dialog shows a **User, group, or email** field and an
   access level that reads **View and use**, and the **OK** and **CANCEL** buttons. Your own account
   is not listed. There is no message box and no **Send notifications** checkbox.
5. Type the first letters of the second user's login into **User, group, or email** and pick the
   user from the suggestions.
6. **Verify:** the user is shown in the field with an **×**. A **Type in message here** box and a
   **Send notifications** checkbox appear; the checkbox is checked and is followed by
   *(Will be sent to 1 user)*. The dialog shows nothing else.
7. Click **CANCEL**.
8. **Verify:** the dialog closes. Click the `ShareDialog1` card: its **Sharing** pane does not name
   the second user.

### 3. Your own account is not offered (GROK-20755)

1. Right-click the `ShareDialog1` card and choose **Share...**.
2. Type the first three letters of your own login into **User, group, or email**.
3. **Verify:** suggestions drop down, and none of them is your own account; no error appears.
4. Click **CANCEL**.

### 4. Escape and CANCEL close what they should (GROK-20758)

1. Right-click the `ShareDialog1` card and choose **Share...**.
2. Click the access level (**View and use**).
3. **Verify:** a list of levels opens: **Full access**, **View and use** (checked), **Edit**,
   **Delete**, **Share**.
4. Press **Escape**.
5. **Verify:** the list closes and the **Share ShareDialog1** dialog stays open.
6. Type one letter into **User, group, or email**.
7. **Verify:** a list of suggestions drops down.
8. Click **CANCEL**.
9. **Verify:** the dialog closes.

### 5. Resulting permissions shows what a user gets

1. Click the `ShareDialog1` card and expand **Sharing** in the **Context Panel**.
2. Click the calculator icon next to **SHARE...**.
3. **Verify:** the **Resulting permissions** dialog opens. It names `ShareDialog1` and has a
   **User, group or role** field and a **CANCEL** button.
4. Type the first letters of the second user's login into **User, group or role** and pick the
   user.
5. **Verify:** the dialog says *No access* for that user.
6. Click **CANCEL**.

## Cleanup

1. Right-click the `ShareDialog1` card, choose **Delete** and confirm.

## Automation notes

- Scenario 1: the calculator icon is the button named "Calculate resulting permissions for this
  entity".
- Scenario 2, step 4: your own grant is in the dialog as a hidden row (`grok-permissions-self`);
  the claim is that it is not shown, not that it is absent from the page.
- Scenario 3: reading the rows of the suggestion list is not a library step yet.
- Scenario 4: the list of levels is a popup outside the dialog; naming it as an element is not a
  library step yet, and picking a level from it is not automated.
- Scenario 5: the dialog has no OK button; **CANCEL** closes it.
- Which privileges each level holds for each kind of entity, and that the server refuses the rest,
  is an API test.
