# Users, groups and roles edge cases — what the features cannot say yet

Written for `group-edge-cases.feature`, `role-edge-cases.feature` and `user-edge-cases.feature`, with no library,
binding or code change. The features claim what the editors show: the only admin of a group refused in MANAGE, a
group nested into its own member refused, a name taken in the New Group dialog refused for a group and for a role,
a name of spaces only refused, CANCEL and Escape changing nothing, and no Disable... on the running account's own
card. What follows has no source for its expected behaviour, belongs to another test layer, or cannot run safely on
a shared stand.

## Still open

### Questions for owners

1. **Two groups with the same visible name.** The name guard (`groups_repository.dart` save, "A group or role with
   the same name already exists") compares the internal names. The New Group dialog sends no internal name and the
   server derives one from the name typed; a group made through the JS API (`DG.Group.create(name)`) has the name as
   its internal name. On localhost (2026-10-26) OK in the New Group dialog with the name of such a group closed the
   dialog with no message, and two groups showed the same name. Should the visible name be unique, whatever made the
   first group?
2. **Invite a Friend... with the email of an existing user.** No source says whether the dialog refuses it, invites
   again, or points to the user.

### API or server tests to write

1. **A user disabling their own account through the API.** Help (govern/access-control/users-and-groups.md,
   Disabling accounts): "you can't disable your own account". The client hides Disable... for the running account
   (`user-edge-cases.feature` claims that), but `UserManagementService.blockUser` checks only the Edit User
   permission: `POST /api/users/block` with the caller's own id answered 200 and ended the session on localhost. A
   datlas test in `users_groups/user_disable_test.dart`, beside a server guard.
2. **The last-admin guard on every route.** The server refuses a group with no admin when the admin unticks Admin
   or is removed through `/groups/{id}/children`, when the group is saved with relations and no admin, and when the
   member revokes its own Admin or leaves through `/groups/{personal}/parents` (`groups_repository.dart` save; all
   five reproduced 2026-10-26). No datlas test covers them; `users_groups/groups_test.dart` is their place. The
   feature claims the MANAGE part once.
3. **A group or role named like an existing group or role.** The same guard for a rename through Properties and
   for a role named after a group, on the server — `groups_test.dart`.
4. **core#561: a role holding a global permission deleted through `grok.dapi.entities.delete`.** The grant stays
   with no group, and the Global Permissions pane of every role shows "error". The test belongs on the server
   (deleting a group as an entity removes its global grants). Not a feature: the role must be deleted through the
   API, and while the defect is there the grant left behind breaks the pane of every role on the stand for every
   other feature until it is deleted by id. The library has no step for that cleanup. Deleting the role from the
   UI is covered by datlas `deleteGroupWithGlobalPermission` (GROK-20904).
5. **core#608: a group holding a domain grant cannot be deleted.** A server outcome — `groups_test.dart`.
6. **A user saved with an email another user has.** The server refuses it ("User with same e-mail already exists",
   `users_repository.dart`); `users_test.dart` covers only a taken login.

## Won't do

- **Saving a user to test the user editor's edge cases** (a taken email, a taken login through the dialog): every
  run would leave a user on the stand for good, since users cannot be deleted. The server guards are datlas
  `users_test.dart` (a taken login; a taken email is item 6 above).
- **MANAGE adding an email that matches no one.** The membership editor offers to create a user for an unknown
  email (`grok_group_meta.dart` `getGroupMembershipEditor`), and SAVE would make one that cannot be deleted.
