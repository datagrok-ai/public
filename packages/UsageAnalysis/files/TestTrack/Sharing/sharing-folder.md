---
feature: sharing
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p1
realizes_atlas: [sharing.cp.cp-shared-with-me-browse-node]
realizes: [sharing.share-dialog, sharing.browse-shared-with-me]
related_bugs:
  - GROK-19264
  - GROK-12795
  - GROK-11266
  - GROK-10449
  - public#4268
---

# Sharing — a folder in your files

A single file has no share command; a folder has **Share folder...**. The first share of a folder
makes a connection for that folder, named `My files/<folder>`, and shares it. The second user then
finds the folder under **Shared with me**.

## Setup

1. A second user account other than yours exists on the server (for example `bddsecond`).
2. In **Browse > Files > My files**, a folder `ShareFolder1` that holds one file, `data1.csv`.
   `ShareFolder1` has never been shared, and **Browse > Files** has no `My files/ShareFolder1`
   node.

## Scenarios

### 1. Share folder opens and cancels cleanly (GROK-19264)

1. In **Browse**, expand **Files > My files** and right-click `ShareFolder1`.
2. **Verify:** the menu has **Share folder...**.
3. Choose **Share folder...**.
4. **Verify:** the **Share folder** dialog opens with a **User, group, or email** field, the access
   level **View and use**, and **OK** and **CANCEL**.
5. Click **CANCEL**.
6. **Verify:** the dialog closes, no error appears, and **Browse > Files** still has no
   `My files/ShareFolder1` node.

### 2. The folder is shared with the second user (GROK-12795, GROK-11266)

1. Right-click `ShareFolder1` and choose **Share folder...**.
2. In **User, group, or email**, type the first letters of the second user's login and pick the
   user. Clear **Send notifications**.
3. Click **OK**.
4. **Verify:** the dialog closes and the balloon *Shared* appears; no error.
5. **Verify:** **Browse > Files** now has a `My files/ShareFolder1` node.
6. Click `My files/ShareFolder1`.
7. **Verify:** the view lists `data1.csv`.
8. Right-click `ShareFolder1` under **My files** and choose **Share folder...** again.
9. **Verify:** the dialog now shows a **Copy link** icon at the top and lists the second user.
10. Click **CANCEL**.

### 3. The second user sees the folder (GROK-10449, public#4268)

1. Sign in as the second user.
2. In **Browse**, expand **My stuff > Shared with me > (your name)**.
3. **Verify:** `My files/ShareFolder1` is listed.
4. Click `My files/ShareFolder1`.
5. **Verify:** the view lists `data1.csv`; it does not ask for credentials.
6. Sign in as yourself again.

## Cleanup

1. In **Browse > Files**, right-click `My files/ShareFolder1`, choose **Delete...** and click
   **DELETE** in the dialog *Delete connection "My files/ShareFolder1"?*.
2. In **Browse > Files > My files**, right-click `ShareFolder1`, choose **Delete...** and click
   **DELETE** in the dialog *Delete folder "...ShareFolder1"?*.

## Automation notes

- The feature is `@serial`: it switches the signed-in account with `user signs in as the sharing
  user` / `user signs in as themselves again`.
- The connection the share makes is named `My files/ShareFolder1` in the tree and
  `<your namespace>:MyFilesShareFolder1` on the server. Deleting it through **Delete...** removed it
  and left no project behind (localhost, 2026-10-09), so the feature can sweep it by that name.
- Scenario 2, step 9: the **Copy link** icon (`icon-link`) is there only once the folder has been
  shared; what it copies is not checked.
- Scenario 3, step 5: the second user gets a *Credentials required* view (Login, Password, SAVE)
  for `My files/ShareFolder1` instead of the file list, while the owner's own view lists the file —
  public#4268 "Sharing: My files: error on opening the shared folder".
- What the second user may do with the files (read, not write) and that files outside the folder
  stay closed to them is an API test.
