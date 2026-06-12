
# Sharing — UI Smoke (Single-Actor)

## Setup

1. Log in to Datagrok as the test owner account (single browser session; no second-actor
   recipient required).
2. Ensure at least one shareable entity exists under the tester's account — use an
   existing project or upload the demog dataset as a TableInfo via `readDataframe`.
3. Open the Browse panel if not already visible (View > Browse).

## Scenarios

### Scenario 1: Open Share dialog via context menu (pcmdShare)

Steps:
1. In the Browse gallery, right-click an entity card owned by the tester (e.g. a project
   or TableInfo in Browse > Platform > Tables).
2. In the context menu, click the **Share...** item (DOM: `[name="div-Share..."]`,
   flow-id: `pcmdShare`).
3. Observe the Share dialog opens.

Expected:
- The Share dialog is visible and contains a PermissionsEditor widget with:
  - A user/group autocomplete input field.
  - An access-level dropdown ("View and use" / "Full access").
  - Existing grant rows for the current user (owner).
  - An "Advanced editor..." link.
  - A notification message textarea and "Send notifications" checkbox.
  - OK and CANCEL buttons.

# atlas entry derived from source: core/client/xamgle/lib/src/commands/file/share_dataset.dart#L7

### Scenario 2: Share dialog — add and cancel (flow-open-share-dialog)

Steps:
1. Open the Share dialog from the context menu as in Scenario 1.
2. Type a partial username into the autocomplete input (e.g. the first 3 characters of
   the tester's own username to get an autocomplete suggestion).
3. Select a suggestion from the autocomplete dropdown.
4. Confirm the access level is set to "View and use" in the dropdown.
5. Click CANCEL.

Expected:
- No permission change is applied.
- The dialog closes without error.
- The entity's permission state remains unchanged (verifiable via
  `grok.dapi.permissions.get(entity)`).

### Scenario 3: Context panel Sharing pane (flow-sharing-context-panel)

Steps:
1. Open an entity in the main view so it becomes the current object (e.g. open the
   demog table in a table view).
2. Open the Context Panel (View > Context Panel or via the right-side panel icon).
3. Click the **Sharing** pane tab in the Context Panel.
4. Observe the grant list rendered by `dapi.privileges.getHumanReadablePermissions()`.
5. Click the **SHARE...** button in the Sharing pane.

Expected:
- The Sharing pane is visible and lists at least the owner's grant as a human-readable
  string ("You are the owner" or "X has full access").
- Clicking SHARE... opens the Share dialog for the current entity (same modal as
  Scenario 1).

# atlas entry derived from source: core/shared/grok_shared/lib/src/http_client/privileges_client.dart#L51

### Scenario 4: Navigate to Advanced editor (flow-advanced-editor)

Steps:
1. Open the Share dialog for any owned entity as in Scenario 1.
2. Click the **Advanced editor...** link inside the Share dialog.
3. Observe the PermissionsView opens at the route `/permissions/<entityId>`.
4. Verify the permission matrix is visible with at least the Common permission columns
   (View, Edit, Delete, Share).
5. Click the browser Back button or navigate away without saving.

Expected:
- The PermissionsView route opens at `/permissions/<entityId>`.
- The permission matrix renders rows for each current grantee and the standard
  Common permission columns (View, Edit, Delete, Share).
- Navigating away without clicking SAVE does not alter permissions.

# atlas entry derived from source: core/client/xamgle/lib/src/features/users/permissions_browser.dart#L403
