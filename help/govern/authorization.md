<!-- TITLE: Authorization -->
<!-- SUBTITLE: -->

# Overview

Datagrok has built-in authorization system, based on user groups permissions. A permissions are an action that can be
applied to an entity type. Each [entity](../overview/objects.md) has list of permissions, which define it's visibility
and behaviour. Most of the entities have permissions to:

* View - User or group can see the entity
* Edit - User or group can edit entity attributes
* Delete - User or group can delete entity
* Share - User or group can edit entity permissions

Some entities can have additional permissions.

All permissions are grouped to "Can view" and "Can edit" groups, so, when user grants "Can edit"
permission on entity - actually all permissions are granted.

Permission can be set to a group, or to a user, using "personal" user group, which is created automatically. Groups can
be nested, so all group members get permission, that were set to a parent group.

See also:

* [Entities](../overview/objects.md)
* [Groups](../govern/group.md)
* [Authentication](authentication.md)
* [Security](security.md)
