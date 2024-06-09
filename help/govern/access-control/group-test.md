<!-- TITLE: Tests: Group management -->
<!-- SUBTITLE: -->

# Tests: group management

User [group](../users-and-groups#users) allow you to flexible configure user rights. [Users_deprecated/user.mdmd) can request membership
in [groups](../users-and-groups#users), start chat with other members, view [group](../users-and-groups#users) members and etc.

## Testing scenarios

1. Open ```Groups``` from ```Admin``` menu

1. Click on ```New group``` button

* "Create new [group](../users-and-groups#users)" dialog is open

1. Click ```OK``` (without filling fields)

* Warning "Group name is empty"
* New [group](../users-and-groups#users) was not created

1. Enter "test" in "Name" field and then click ```OK```

* New [group](../users-and-groups#users) was created with name "test"

1. Open "Edit" dialog for "test" [group](../users-and-groups#users). (From context menu or from "Actions"
   tab on [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel))

* Editing dialog is open. Here you can change name and description of [group](../users-and-groups#users)

1. Open "Edit members" dialog for "test" [group](../users-and-groups#users). (From context menu or from "
   Actions" tab on [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel))

* "Edit members" dialog is open
* In list of [group](../users-and-groups#users) members only current user (who created [group](../users-and-groups#users)) with admin
  role

1. Click on "admin" label near the username

* User role is switched to "user" and display accordingly

1. Click ```Save``` button

* Warning "[group](../users-and-groups#users) should have at least one admin"

1. Switch role back to "Admin"

1. Start typing username in search field

* Drop-down list shows the users found
* If found user is not in [group](../users-and-groups#users), you can add it (next to username will be button ```+``` )
* If found user is already a [group](../users-and-groups#users) member you can change its role or delete it from [group](../users-and-groups#users)

1. Start typing name of another [group](../users-and-groups#users) in search field

* You can add one [group](../users-and-groups#users) to another
* All users belonging to nested [group](../users-and-groups#users) get role and rights that are given to this
  [group](../users-and-groups#users) within the parent [group](../users-and-groups#users)
* Group Management as member of another [group](../users-and-groups#users) is similar to managing users in [groups](../users-and-groups#users)

1. Add one user and one [group](../users-and-groups#users) to "test" [group](../users-and-groups#users)

1. Change your role to "user", and roles of other members to "admin"

1. Click on delete user icon next to username ( ```X``` ) of another user in "Edit members" dialog, then
   click ```SAVE``` button

* User was not deleted
* Warning massage "Insufficient privileges"

1. Click on delete user icon next to self username ( ```X``` ) in "Edit members" dialog, then click ```SAVE``` button

* You left the group "test"

1. Request membership in "test" [group](../users-and-groups#users). (From context menu or from "Actions" tab
   on [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel))

* Balloon with information about successful request
* All admins and members of nested groups with admin role received a notification about request
* User's notification displays this request with current status (```?``` - pending request, ```x``` - rejected
  request, ```✓``` - request accepted)

1. Open "Edit memberships" dialog for "test" [group](../users-and-groups#users). (From context menu or from "Actions"
   tab on [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel), it is necessary to do from user with admin role in
   this group)

* "Edit memberships" dialog is open
* In this dialog you can request membership of current group in other groups
* Management in this dialog is similar to others, but carried out from the position of nested group

1. Open "Chat" for "test" [group](../users-and-groups#users). (From context menu or from "Actions" tab
   on [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel))

1. Delete "test" [group](../users-and-groups#users) (From context menu or from "Actions" tab
   on [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel))

See also:

* [Group](../users-and-groups#users)
* [User](../_deprecated/user.md)
