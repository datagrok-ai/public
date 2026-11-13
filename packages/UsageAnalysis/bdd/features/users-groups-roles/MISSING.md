# Users, groups and roles — what the features leave out

## Still open

- **An account can disable itself through the API.** Help: "you can't disable your own account"; only the client
  hides Disable... (`users-view` claims that), while `UserManagementService.blockUser` accepts the caller's own id.
  Needs a server guard and a ticket; claim when fixed: a datlas test in `user_disable_test.dart`.
- core#561 (a role deleted through `grok.dapi.entities.delete` leaves its global grant, and the Global Permissions
  pane of every role says "error") and core#608 (a group holding a domain grant cannot be deleted): server defects;
  claim when fixed: datlas `groups_test.dart` (`deleteGroupWithGlobalPermission` covers `groups.delete` only).

## Won't do

- A taken group or role name, the only admin unticking Admin or leaving, a group nested into its own member: ApiTests
  `Dapi: groups` ('a name a group or role already has is refused', 'the only admin cannot leave a group or lose
  Admin') and datlas `saveGroupWithCircularReference`; the UI claims one refusal (`groups-lifecycle`, Groups-04).
- A name of spaces in New Group, New Role or Properties...: the same `propertiesDialog` check as the emptied name
  (Groups-05, Roles-05).
- CANCEL in MANAGE and in Global Permissions: both dialogs write only on SAVE.
- A user saved with an email another user has: the server's refusal ("User with same e-mail already exists"), a
  datlas `users_test.dart` case to add (it has a taken login only), never a feature (users cannot be deleted).
- Saving a user from the editor (a taken email or login), MANAGE adding an unknown email: each leaves a user on the
  stand for good.
