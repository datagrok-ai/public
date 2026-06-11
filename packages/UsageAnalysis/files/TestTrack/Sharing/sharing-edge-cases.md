---
feature: sharing
sub_features_covered:
  - sharing.advanced-editor
  - sharing.server.privileges-router
  - sharing.permission-types
  - sharing.share-dialog
target_layer: manual-only
coverage_type: edge
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
---

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

## Notes

- target_layer rationale: Scenarios 1 and 2 require server-side policy enforcement
  responses that cannot be mocked in Playwright without a live server. Scenario 1
  requires a user account that genuinely lacks `ShareWithEveryone` — a real permission
  state, not simulatable in a single-session environment without multi-account setup.
  Scenario 2's expected result depends on server behavior at `/privileges/permissions`
  DELETE endpoint (whether it enforces owner-retention or silently ignores the delete).
  All three scenarios require live server interaction and, for Scenarios 1-2, a known
  account permission state — classified manual-only.
- Atlas manual_only[] note: None of the sub_features_covered here appear in manual_only[].
  The manual-only classification is driven by the server-side policy enforcement
  requirement, not the two-actor constraint that governs entity-type scenarios.
- Scenario 2 directly addresses the unresolved_ambiguity in sharing-all-users-and-owner.md
  (Block C step 1: "open the Advanced editor and remove every grant — or attempt to").
  The expected result in this scenario refines that ambiguity: mechanism is either
  silent-preserve by server OR UI lock — both are acceptable; ownerless entity is not.
- Deferrals: New-user invitation flow (atlas edge_cases[4], sharing.notification +
  sharing.share-dialog) deferred — requires a real email address without an existing
  Datagrok account and control over mailbox delivery; no email-server helper available.
  Project-dependency cascade edge case (atlas edge_cases[0], GROK-19403) deferred —
  sharing.entity-types.project is in manual_only[] (two-actor requirement) and the
  cascade behavior is a project-entity-type test, not a permission-policy test.
- See: core/docs/ENTITY_PERMISSIONS.md#Five Core Concepts
- See: core/docs/ENTITY_PERMISSIONS.md#Privileges (what you can do)
- See: public/help/govern/access-control/access-control.md#Authorization
- See: public/help/govern/access-control/access-control.md#Global Permissions
