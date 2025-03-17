---
title: 'Users and groups'
format: 'mdx'
unlisted: false
sidebar_position: 1
---

```mdx-code-block
import NewUser from '../img/create-new-user-dialog.png';
import NewServiceUser from '../img/create-new-service-user-dialog.png';
import InviteUser from '../img/invite-a-friend-dialog.png';
```

## Users

A _user_ represents the identity of a person, and is used for [group](#groups) and role
management. To view or manage users, open the [Users View](https://public.datagrok.ai/users?) (**Sidebar > Browse (<FAIcon
icon="fa-solid fa-compass"/>) > Platform > Users**).

A user is an [entity](../../datagrok/concepts/objects.md). This means a common set of
operations apply to it, like getting its URL or using it as a parameter in the
[audit record](../audit/audit.md). Like with other entities, you can [search](../../datagrok/navigation/views/browse.md#entity-search) for users using its [parameters](../../datagrok/concepts/objects.md#parameters).

<details>
<summary>User entity parameters</summary>

| Field       | Description                                            |
|-------------|--------------------------------------------------------|
| ID          |                                                        |
| name        |                                                        |
| login       |                                                        |
| email       |                                                        |
| createdOn   |                                                        |
| updatedOn   |                                                        |
| group       | [Group](#groups) object: User personal security group |
| author      | [User](#users) object                                 |
| starredBy   | [User](#users) object                                 |
| commentedBy | [User](#users) object                                 |

</details>

### Adding users

To create a new user: 
1. On the **Sidebar**, click **Browse (<FAIcon icon="fa-solid fa-compass"/>) > Platform > Users**.
1. Click the **NEW** button and select the option you want:

   * **User...**: Select this option to create a user without sending an invitation to sign up. This option is typically used to create the initial admin users who can later invite other users to sign up via URL or email invitation, or add users using services like OAuth or OpenID.
     <details>
     <summary>Instructions</summary>
     1. In the **New User** dialog, enter the user's name, email, and login.
     <img src={NewUser} width="200px"/>
     1. Click **OK** to create a user.
     </details> 

   * **Service user...**: Choose this option to create a service user account.
     <details>
     <summary>Instructions</summary>
     1. In the **New Service User** dialog, enter the login for the service user.
     <img src={NewServiceUser} width="200px"/>
     1. Click **OK** to create the service user.   
     </details>

   * **Invite a Friend...**: Select this option to email an invitation to sign up to the platform.
     <details>
     <summary>Instructions</summary>
     1. In the **Invite User** dialog, enter the email address of the user you want to invite.
     <img src={InviteUser} width="200px"/>
     1. Click OK to send the invitation.
     </details>

To manage users, in the **Users View**, find the user you want, and right-click it to access available actions.

Only [Administrators](#group-types) with global permissions can add or remove users.

### Disabling accounts

To disable a user account, login as [administrator](#group-types), navigate
to **Browse | Platform | Users**, right-click on the user, and select "Block".
This will prevent user from logging in and using the platform, and this user 
will not count towards the license. 

All assets that the user has created will continue to be available in the system.
Administrators can [share](../../datagrok/navigation/basic-tasks/basic-tasks.md#share) 
them with others if necessary.  

Currently, there is no way to permanently delete a user. We are planning 
to implement it in the future versions.

## Groups

A _group_ is a named collection of users that share permissions. To view 
or manage groups, open the [Groups View](https://public.datagrok.ai/groups?) (**Sidebar > Browse (<FAIcon
icon="fa-solid fa-compass"/>) > Platform > Groups**).

A _group_ is an [entity](../../datagrok/concepts/objects.md), which means a common set of
operations apply to it, like getting its URL or using it as a parameter in the
[audit record](../audit/audit.md). Like with other entities, you can [search](../../datagrok/navigation/views/browse.md#entity-search) for groups using the group's [parameters](../../datagrok/concepts/objects.md#parameters).

<details>
<summary>Group entity parameters</summary>

You can use these fields to filter groups with smart search:

| Field       | Description                                        |
|-------------|----------------------------------------------------|
| ID          |                                                    |
| name        |                                                    |
| isPersonal  |                                                    |
| parents     |                                                    |
| members     |                                                    |
| createdOn   |                                                    |
| updatedOn   |                                                    |
| user        | [User](#users) object: User, if group is personal |

</details>

### Group structure

Every user automatically forms their own group, but groups can also consist
of multiple users and other groups. This grouping mechanism lets you
create hierarchies that mirror your organization’s structure or security
requirements. For example:

* **Group 1** has members **User 1**  and **Group 2**
* **Group 2** has members **User 2** and **User 3**
* Because **Group 2** is a member of **Group 1**, **User 2** and **User 3** are also members of **Group 1**.

![Group members and memberships](../img/groups.png)

  > Note: Child groups inherit [permissions](access-control.md#permissions) of the parent group.

Regardless of group membership, any user can do the following actions with respect to a group:

* Chat with group members (right-click > **Chat**)
* Request membership (right-click > **Request membership**)
* [Share entities](../../datagrok/navigation/basic-tasks/basic-tasks.md)

### Group types

Datagrok automatically creates several key groups upon deployment, each designed
with specific roles and permissions:

* **All users**: This group includes all users and groups and initially comes
  with a basic set of [permissions](access-control.md#permissions).
* **Administrators**:
     * During the deployment process, the Administrators group is created and
       granted all available permissions, ensuring complete control over the
       platform.
     * An 'admin' user and password is provided in the deployment script.
     * Immediately following the deployment, logging in as the admin user is the
       standard procedure to begin configuring and managing the Datagrok
       instance.
     
       :::danger
       
       Exercise caution when modifying the Administrators group. Modifying
       or deleting this group without a functional replacement may result in a
       loss of all administrative capabilities on the platform.
       
       :::

* **Developers**: Initially created as a child group under Administrators,
  this group inherits the permissions from its parent group.

Members of the Administrators group have global permissions, accessible via **Top Menu > Admin > Global Permissions...** 

The following operations require global permissions:
    * Creating a new user - `CreateUser`
    * Inviting a user - `InviteUser`
    * Editing a user - `EditUser`
    * Editing a group - `EditGroup`
    * Editing global permissions - `EditGlobalPermissions`
    * Editing server settings - `EditPluginsSettings`
    * Start Admin Session (disable all permissions check during current session) - `StartAdminSession`
    * Deploy or install a package - `PublishPackage`
    * Delete a comment in any chat - `DeleteComments`
    * Create or edit entity type (see [Sticky Meta](../catalog/sticky-meta.md)) - `SaveEntityType`
    * Modify any system pre-created data connection, such as `Datagrok`, `DatagrokAdmin` or `AppData` - `AdminSystemConnections`
    * Create anything - `CreateEntity`
    * Create a [script](../../compute/scripting/scripting.mdx) - `CreateScript`

### Managing groups

The following actions are available from the group's context menu (available on right-click):

* **Edit...**: Edit the group's name and description. Generate a URL link.
* **Edit members**: Add or remove group members
* **Edit memberships**: Add or remove the group to/from other groups
* **Delete**: Delete a group

#### Creating a group

Any user can create a group. To create a group:

1. Go to **Sidebar** > **Browse** (<FAIcon icon="fa-solid fa-compass"/>) > **Platform** > **Groups**. A **NEW GROUP...** button appears on the **Top Menu**.
1. On the **Top Menu**, click the **NEW GROUP...** button to open the **Create New Group** dialog.
1. In the dialog, fill out the group's name and, optionally, description.
1. Click **OK**.

Within a group, one or more members can be assigned as group admins. Group admins can:

* add and remove members
* approve and deny membership requests
* edit group name

If user is the only admin in group, they can't leave the group or revoke their admin privileges.

#### Adding users to a group

To add users to a group, you have two options:
1. **Manually add members** (right click the group and select **Add members...**)
1. Invite users to sign up via URL

##### Inviting users via URL

Group administrators can invite users to join a particular group via URL. This
feature is useful for onboarding multiple users simultaneously or for tracking
signups after specific events such as webinars or collaborative projects.

To create an invitation link:

1. Go to **Sidebar** > **Browse** (<FAIcon icon="fa-solid fa-compass"/>) > **Platform** > **Groups**. A **Groups View** opens.
1. In the **Groups View**, right-click the group and select **Edit…** A dialog opens.
1. Click the **Gear (<FAIcon icon="fa-solid fa-gear"/>) icon** and enter or generate a password in the **Password** field. By default, the dialog displays an autogenerated password.
1. Copy the password and use it to create an invitation link as follows:

      `<instance URL>?groupPassword=<password>`

      For example:

      `public.datagrok.ai/?groupPassword=w0TDE6RcpH8XO0ZIBSauVLos`
    
1. Copy the URL and distribute to recipients.

![Granting membership via URL](../img/group-membership.gif)

Only the group administrator can create an invitation link.



