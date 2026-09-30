---
title: "Roles"
sidebar_position: 1.5
description: Use roles to grant capabilities such as global permissions and access to functions and packages, separately from the groups that mirror your organization.
keywords:
  - roles
  - role-based access control
  - global permissions
  - create role
  - administrators role
---

A _role_ is a named set of capabilities you assign to people. Roles and
[groups](users-and-groups.md#groups) share the same membership mechanics, and technically either
can be granted any permission: a group can hold a global permission, and a space or connection can
be shared with a role. We recommend using them for different jobs:

| By convention | Groups | Roles |
|---|---|---|
| **Describe** | Your organization: departments, sites, teams, programs | What someone is allowed to do |
| **Carry** | Sharing: who can view and edit a space, query, connection, or dashboard | Capabilities: [global permissions](access-control.md#global-permissions) and access to functions and packages |
| **Typically come from** | Your identity provider, through [group synchronization](../../deploy/complete-setup/configure-auth.md#group-synchronization) | Datagrok administrators |
| **Examples** | Oncology, Chemistry, a partner program | Authors, Consumers, Developers |

If administrators follow this convention, a change in someone's organization never changes what
they can do, and a change in their job never changes what they can see. For example, Priya is a
chemist in the _Oncology_ group with the _Authors_ role. The group
determines what she can see: Oncology's spaces, dashboards, and database credentials. The role
determines what she can do: create and share dashboards and queries. If she moves to Immunology,
her group changes and she sees Immunology's data instead, and she can still create dashboards. If
she becomes her team's assay owner, an administrator adds the _Content Owners_ role, and she can
publish shared queries without getting access to any new data. See
[Running Datagrok in the enterprise](../enterprise-guide.md#groups-and-roles).

Roles are available in Datagrok 1.27.0 and later.

## How roles work

* **Members.** A role's members can be users, groups, or other roles. Permissions flow from parent
  to child, so a member inherits everything granted to the role. Assigning a role to a group, for
  example the Authors role to a Chemistry group, makes the role a parent of the group, and everyone
  in the group inherits it. On a group, the **Roles** tab does this for you.
* **Nesting.** Roles can nest like groups. Circular membership is rejected, and a group or role
  can't be its own parent or child.
* **Names.** Groups and roles share one namespace, so a role can't have the same name as a group.
* **Additive permissions.** A user's effective permissions are the union of everything granted to
  every group and role they hold. There is no deny: a role can't take away a capability that
  another group or role grants.
* **Sharing.** A role counts as a common group for sharing: the **Share With Everyone** permission
  is only needed to share with someone you have no groups _or roles_ in common with.

### The Administrators role

**Administrators** is a built-in role. Its members hold all global permissions, and groups created by
[group synchronization](../../deploy/complete-setup/configure-auth.md#group-synchronization) are
owned by it. Keep its membership very small, and delegate day-to-day membership work to
[group admins](users-and-groups.md#managing-groups) instead.

:::danger

Don't remove the last member of the Administrators role, or change its permissions, without a
working replacement. You can lose all administrative access to the platform.

:::

## What to grant to roles

* **Global permissions.** Creation rights such as **Create Dashboard**, **Create Script**, and
  **Create Data Query**; administrative rights such as **Publish Package** and **Edit Group**; and
  the **Browse** permissions that show or hide sections of the Browse panel. See the full list in
  [Global permissions](access-control.md#global-permissions).
* **Functions and packages.** Share a function, script, query, or package with a role to control
  who can run it. This is ordinary [sharing](access-control.md#permissions), with the role as the
  recipient; there is no separate package permission for roles.

Keep content sharing on groups. If a space or dashboard is shared with a role, access to it follows
a capability rather than a team, which is hard to audit later.

## Managing roles

To create a role, you need the **Create Role** global permission. To see the **Roles** section in
the Browse panel, you need **Browse Roles**. Both are held by the Administrators role by default.

Roles interact with the rest of the platform as follows:

* **Group synchronization never touches roles.** OpenID, Entra ID, and Google Workspace sync match,
  create, and remove groups only. Roles, including Administrators, are never matched, so sync doesn't
  grant or remove capabilities. (Keycloak synchronization is different: it maps Datagrok roles to
  Keycloak realm roles.)
* **Scheduled functions can run as a role.** When you
  [schedule a function](../../datagrok/concepts/functions/functions.md), you can choose a group or
  role to run it as, if you're a member.
* **Automation.** In the API, a role is a group with a role flag, and the group endpoints return
  and accept roles. The [`grok s groups`](../../develop/server-management.md#manage-users-and-groups)
  commands work on roles too: `groups save` with `"isRole": true` in the JSON creates a role, and
  `add-members` and `list-members` manage its members.

## See also

* [Users and groups](users-and-groups.md)
* [Access control](access-control.md)
* [Running Datagrok in the enterprise](../enterprise-guide.md)
* [Configure authentication](../../deploy/complete-setup/configure-auth.md)
