<!-- TITLE: Tests: Group management -->
<!-- SUBTITLE: -->

# Tests: group management

User [group](../govern/group.md) allow you to flexible configure user rights. [Users](user.md) can request membership
in [groups](../govern/group.md), start chat with other members, view [group](../govern/group.md) members and etc.

## Testing scenarios

1. Open ```Groups``` from ```Admin``` menu

1. Click on ```New group``` button

* "Create new [group](../govern/group.md)" dialog is open

1. Click ```OK``` (without filling fields)

* Warning "Group name is empty"
* New [group](../govern/group.md) was not created

1. Enter "test" in "Name" field and then click ```OK```

* New [group](../govern/group.md) was created with name "test"

1. Open "Edit" dialog for "test" [group](../govern/group.md). (From contex menu or from "Actions"
   tab on [Property Panel](../overview/navigation.md#properties))

* Editing dialog is open. Here you can change name and description of [group](../govern/group.md)

1. Open "Edit members" dialog for "test" [group](../govern/group.md). (From contex menu or from "
   Actions" tab on [Property Panel](../overview/navigation.md#properties))

* "Edit members" dialog is open
* In list of [group](../govern/group.md) members only current user (who created [group](../govern/group.md)) with admin
  role

1. Click on "admin" label near the username

* User role is switched to "user" and display accordingly

1. Click ```Save``` button

* Warning "[group](../govern/group.md) should have at least one admin"

1. Switch role back to "Admin"

1. Start typing username in search field

* Drop-down list shows the users found
* If found user is not in [group](group.md), you can add it (next to username will be button ```+``` )
* If found user is already a [group](group.md) member you can change its role or delete it from [group](group.md)

1. Start typing name of another [group](group.md) in search field

* You can add one [group](group.md) to another
* All users belonging to nested [group](group.md) get role and rights that are given to this
  [group](group.md) within the parent [group](group.md)
* Group Management as member of another [group](group.md) is similar to managing users in [groups](group.md)

1. Add one user and one [group](group.md) to "test" [group](group.md)

1. Change your role to "user", and roles of other members to "admin"

1. Click on delete user icon next to username ( ```X``` ) of another user in "Edit members" dialog, then
   click ```SAVE``` button

* User was not deleted
* Warning massage "Insufficient privileges"

1. Click on delete user icon next to self username ( ```X``` ) in "Edit members" dialog, then click ```SAVE``` button

* You left the group "test"

1. Request membership in "test" [group](group.md). (From contex menu or from "Actions" tab
   on [Property Panel](../overview/navigation.md#properties))

* Balloon with information about successful request
* All admins and members of nested groups with admin role received a notification about request
* User's notification displays this request with current status (```?``` - pending request, ```x``` - rejected
  request, ```✓``` - request accepted)

1. Open "Edit memberships" dialog for "test" [group](group.md). (From contex menu or from "Actions"
   tab on [Property Panel](../overview/navigation.md#properties), it is necessary to do from user with admin role in
   this group)

* "Edit memberships" dialog is open
* In this dialog you can request membership of current group in other groups
* Management in this dialog is similar to others, but carried out from the position of nested group

1. Open "Chat" for "test" [group](group.md). (From contex menu or from "Actions" tab
   on [Property Panel](../overview/navigation.md#properties))

1. Delete "test" [group](group.md) (From contex menu or from "Actions" tab
   on [Property Panel](../overview/navigation.md#properties))

See also:

* [Group](group.md)
* [User](user.md)
