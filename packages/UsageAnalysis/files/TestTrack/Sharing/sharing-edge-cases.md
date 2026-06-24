
# Sharing — Edge Cases (Permission Policy Enforcement)

## Setup

1. Log in as the test owner account (Browser session A — owner).
2. Ensure you own at least one shareable entity (e.g. a project created from the demog
   table, or a fresh script created via Browse > Platform > Functions > Scripts > New).
3. Confirm the tester account does NOT hold the `ShareWithEveryone` global permission
   unless otherwise noted in the scenario. Check via
   `dapi.privileges.checkGlobalPermission('ShareWithEveryone')`.

## Scenarios

### Scenario 1: SHARE_WITH_EVERYONE gate — blocked when permission absent

Steps:
1. Log in as a non-admin user who does NOT hold the `ShareWithEveryone` global
   permission (confirm via `GET /privileges/permissions/check/global/ShareWithEveryone`
   returns false).
2. Right-click an owned entity in Browse gallery → click **Share...** to open the
   Share dialog.
3. In the autocomplete input, type "All users" and select the **All users** built-in
   group from the dropdown.
4. Click OK to attempt sharing with the All users group.
5. Observe the dialog response.

Expected:
- The Share dialog either:
  (a) Blocks the "All users" selection with an inline warning ("You do not have
      permission to share with all users"), OR
  (b) Accepts the input but the server returns an error on POST
      `/privileges/permissions` and the dialog surfaces the error rather than
      silently applying a partial grant.
- Under no circumstances should the entity appear as viewable by anonymous or
  all-user sessions without the owner holding `ShareWithEveryone`.

# atlas entry derived from source: core/shared/grok_shared/lib/src/privileges.dart#L33

### Scenario 2: Owner always retains full access (server-side enforcement)

Steps:
1. As the owner, open the Share dialog for an owned entity → click **Advanced editor...**
   to open PermissionsView at `/permissions/<entityId>`.
2. In the permission matrix, locate the row(s) for the current owner (your own user or
   the group to which you belong as owner).
3. Attempt to uncheck ALL permission cells in the owner's row(s) — View, Edit, Delete,
   Share — so that the row has zero checked cells.
4. Click the **SAVE** ribbon button to commit the change via
   `dapi.privileges.savePermissions()` and `dapi.privileges.deletePermission()`.
5. After save, reload the Advanced editor for the same entity (navigate to
   `/permissions/<entityId>` again) and inspect the permission matrix.
6. Also attempt to access the entity from a fresh navigation to confirm the owner can
   still open it.

Expected:
- After saving, the owner retains at least View + Edit + Delete + Share permissions
  on the entity — either:
  (a) The server silently preserves the owner grant row even after the UI sent
      delete requests for those cells (the matrix reload shows the owner rows
      re-populated), OR
  (b) The UI blocks removal of the owner's grant row before the save request is
      issued (e.g. cells are disabled or locked for the owner's own entry).
- The owner can still open and modify the entity after the save attempt.
- The server does NOT respond with 200 OK for a delete on the owner's primary grant
  and leave the entity ownerless.

# atlas entry derived from source: core/server/datlas/lib/src/routers/privileges.dart#L4

### Scenario 3: Co-shared dependent entity notice in PermissionsEditor

Steps:
1. As the owner, ensure a Script entity is a member of at least one existing Project
   (i.e. the script is added to the project's entity list).
2. Right-click the **Script** entity directly (not the Project) in Browse > Platform >
   Functions > Scripts → click **Share...**.
3. Observe the PermissionsEditor widget in the Share dialog before adding any recipient.

Expected:
- The PermissionsEditor displays a notice or info panel listing the Projects that
  contain this Script — indicating those Projects will be co-affected when the Script
  is shared.
- The notice text references at least the name of the parent Project that embeds the
  Script.
- The dialog still allows the share to proceed after acknowledging the notice.

# atlas entry derived from source: core/client/xamgle/lib/src/features/users/permissions_browser.dart#L5
