---
feature: sharing
target_layer: apitest
coverage_type: regression
priority: p0
realizes_atlas: [cp-api-grant-check-revoke, cp-api-check-permission]
realizes: []
realized_as:
  - sharing-api-permissions-api-spec.ts
related_bugs: []
---


# Sharing — API Permissions (JS API Grant / Check / Revoke Lifecycle)

Verifies the JS API's `grok.dapi.permissions` methods — `grant()`, `get()`, `check()`, and
`revoke()` — directly, without going through the Share dialog UI: granting view or edit
access to a group, confirming the grant is visible and enforced, and confirming it is
gone again after revoke.

## Setup

1. Authenticate as the test owner account via `loginToDatagrok`.
2. Create a throwaway entity to use as the permission target throughout all scenarios.
   A Script entity is the lightest option (no external connection dependency):
   ```
   const script = await grok.dapi.scripts.save(DG.Script.create('test_sharing_script'));
   ```
   Record `script.id` for cleanup.
3. Obtain a test group to grant permissions to — use an existing non-admin group (e.g.
   "All users" for view-only tests, or a dedicated QA group). Fetch via:
   ```
   const group = await grok.dapi.groups.filter('All users').first();
   ```
4. Teardown: after all scenarios, revoke any remaining grants on the script entity and
   delete the script via `grok.dapi.scripts.delete(script)`.

## Scenarios

### Scenario 1: Full grant-get-revoke lifecycle

Steps:
1. Call `grok.dapi.permissions.grant(script, group, false)` to grant view permission
   (edit=false) to the test group.
2. Call `grok.dapi.permissions.get(script)` and inspect the returned Map.
3. Verify the test group appears in the returned Map's view array.
4. Verify the test group does NOT appear in the returned Map's edit array.
5. Call `grok.dapi.permissions.revoke(group, script)` to revoke all permissions.
6. Call `grok.dapi.permissions.get(script)` again.
7. Verify the test group no longer appears in either view or edit array.

Expected:
- After step 1: `grok.dapi.permissions.grant()` completes without error (backed by
  POST `/privileges/permissions`).
- After step 2-4: `grok.dapi.permissions.get()` returns a Map with `view` array
  containing the test group's id, and `edit` array not containing the test group.
- After step 5: `grok.dapi.permissions.revoke()` completes without error (backed by
  DELETE `/privileges/permissions/<id>`).
- After step 6-7: `grok.dapi.permissions.get()` returns a Map where the test group
  is absent from both `view` and `edit` arrays.

### Scenario 2: Grant full access and verify via get()

Steps:
1. Call `grok.dapi.permissions.grant(script, group, true)` to grant edit permission
   (edit=true) to the test group.
2. Call `grok.dapi.permissions.get(script)` and inspect the returned Map.
3. Verify the test group appears in the returned Map's edit array.
4. Call `grok.dapi.permissions.revoke(group, script)` to clean up.

Expected:
- After step 1: grant completes without error.
- After step 2-3: `get()` returns a Map with `edit` array containing the test group.
- After step 4: revoke completes without error; group is absent from `get()` results.

### Scenario 3: Permission check — view vs edit boundary

Steps:
1. Call `grok.dapi.permissions.grant(script, group, false)` to grant view-only access.
2. Call `grok.dapi.permissions.check(script, 'Edit')` for the current user's perspective
   on the test group's effective permissions.
3. Call `grok.dapi.permissions.check(script, 'View')` for the same entity.
4. Call `grok.dapi.permissions.revoke(group, script)`.
5. Call `grok.dapi.permissions.grant(script, group, true)` to grant edit access.
6. Call `grok.dapi.permissions.check(script, 'Edit')` again.
7. Call `grok.dapi.permissions.revoke(group, script)` to clean up.

Expected:
- After step 2 (view-only grant): `check(script, 'Edit')` returns false.
- After step 3 (view-only grant): `check(script, 'View')` returns true.
- After step 6 (edit grant): `check(script, 'Edit')` returns true.
- All check() calls are backed by GET
  `/privileges/permissions/check/<entityId>/<permission>` and return boolean.

### Scenario 4: Server router — batch permissions list endpoint

Steps:
1. Grant view permission to two different groups on the script entity using
   `grok.dapi.permissions.grant()` twice.
2. Call `grok.dapi.permissions.get(script)` and verify both groups appear in the
   view array.
3. Revoke both groups via two `grok.dapi.permissions.revoke()` calls.
4. Call `grok.dapi.permissions.get(script)` and verify both groups are absent.

Expected:
- POST `/privileges/permissions` (invoked by grant) processes each grant independently
  and returns success for both.
- GET `/privileges/permissions?entityId=<id>&all=true` (invoked by get) returns a
  consistent view reflecting both grants, then zero grants after revocation.
- DELETE `/privileges/permissions/<id>` (invoked by revoke) removes each grant row
  independently.
