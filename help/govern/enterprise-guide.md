---
title: "Running Datagrok in the enterprise"
sidebar_position: 0
description: Recommendations for a large Datagrok deployment - users and license tiers, data isolation, content versions, monitoring, and environments.
keywords:
  - enterprise governance
  - license tiers
  - roles and groups
  - data isolation
  - content versioning
  - environment management
  - dev test prod
  - observability
---

This page collects our recommendations for running Datagrok in a large organization:

1. [Users and license tiers](#users-and-license-tiers)
2. [Data isolation](#data-isolation)
3. [Content versions and ownership](#content-versions-and-ownership)
4. [Monitoring](#monitoring)
5. [Environments](#environments)

Settle users, data isolation, and environments before people start creating content. They are
harder to change later.

## Users and license tiers

### Groups and roles

Datagrok grants permissions to [groups](access-control/users-and-groups.md#groups) and
[roles](access-control/roles.md):

* **Groups** mirror your organization and control sharing: who can view or edit a space, a query, a
  connection, or a dashboard.
* **Roles** control capabilities: [global permissions](access-control/access-control.md#global-permissions)
  such as **Create Dashboard** or **Publish Package**, and access to functions and packages.

For example, Priya is a medicinal chemist. Your identity provider puts her in the _Oncology_ and
_PRAME program_ groups, so she can open the PRAME space and use the Oncology credential on the
warehouse connection. In Datagrok she has the _Authors_ role, so she can create and share
dashboards and queries.

* When she moves to Immunology, the identity provider changes her groups. At her next login she sees
  Immunology's data instead of PRAME's. Her role is unchanged, so she can still create dashboards.
* When she becomes her team's assay owner, an administrator adds the _Content Owners_ role. She can
  publish the team's shared queries, and her data access is unchanged.

Grant capabilities on roles and use groups for sharing. Groups and roles can be nested; keep it to
two or three levels so that it stays clear who can see what.

Every user has a personal group. Don't share with it: access given to one person is lost when they
change jobs.

### License tiers

Skip this section if your license doesn't define user tiers.

Create one group per tier in your identity provider, for example `DG-Tier-Full` and
`DG-Tier-Viewer`. In Datagrok these groups are used to count seats. Don't grant them permissions.

For each tier, track:

* **Provisioned**: members of the tier group in your identity provider
* **Active**: provisioned users who signed in during the last 30 days, from the sign-in records in
  the [audit trail](audit/audit.md)
* **Dormant**: provisioned users who aren't active. To free a seat, remove the person from the tier
  group and [block the user](access-control/users-and-groups.md#disabling-accounts). A blocked user
  can't sign in, and their content stays in Datagrok.

Login-time [group synchronization](../deploy/complete-setup/configure-auth.md#group-synchronization)
only updates a user's groups when they sign in. A user who stops signing in stays in the tier group,
so the count drifts upward. For accurate counts, run a scheduled job that reads the tier groups from
your directory and applies them with `grok s groups add-members` and `remove-members`. See
[Sync an AD group with Datagrok](../develop/server-management.md#sync-an-ad-group-with-datagrok).

### Group synchronization

Set up single sign-on and turn on
[group synchronization](../deploy/complete-setup/configure-auth.md#group-synchronization):

* Sync runs at sign-in, so a membership change takes effect at the user's next login. To apply
  changes sooner, use a scheduled job as described above.
* Sync doesn't remove groups or memberships that an administrator created, and it never adds anyone
  to a role.
* Sync matches groups by name. If you create a group with the same name as a group in your identity
  provider, sync adds the provider's members to it. Use a naming convention for groups you create in
  Datagrok.
* With Microsoft Entra ID, every group a user belongs to becomes a Datagrok group, and the Microsoft
  Graph permissions need tenant admin consent. Request consent early; in a large organization it
  usually needs a change request.

### Administrators

Keep the `Administrators` role small, and let group admins manage membership of their own groups.

Administrators see everything, so test access restrictions with an ordinary user account.

## Data isolation

### One space per project

We recommend a [space](../datagrok/concepts/project/space.md) for each discovery project, owned by
the project team's group, with access set from the project record. Create and update these spaces
from your project system rather than by hand. A scheduled script using the
[`grok s`](../develop/server-management.md) CLI can do this. It should:

* Create a space for each new project, and archive the spaces of closed projects. Deleting a space
  removes its content for everyone.
* Use the project's ID, not its name, to identify the space, so renaming a project changes only the
  label
* Give the owning group edit access and the reader group view access, based on the project's access
  columns, on every run
* Have a dry-run mode that lists the changes without making them

If a program has several projects and the program team needs access to all of them, add a program
space above the project spaces.

Set this up before people start saving content. Otherwise content stays in personal spaces and is
shared item by item.

### Organizing spaces

A space is a folder with permissions. Child spaces inherit the parent's permissions, and an item
moved into a space gets the space's permissions. To show a dashboard in several spaces, add a link
to it rather than a copy.

Because child spaces inherit permissions, every level in the space tree is an access boundary:

* **Add a level only where access differs.** A project or a partnership needs its own access. A
  therapeutic area or a modality usually doesn't.
* **Record classifications as properties.** With [sticky meta](catalog/sticky-meta.md), a space can
  have properties such as therapeutic area, target, modality, and phase, and people can search and
  filter spaces by them. If a project changes therapeutic area, edit the property. Moving the space
  would change who can see it.
* **Keep properties few.** Add a property only if people search for spaces by it, and don't use
  properties to control access.

### Where access is enforced

1. **In the database.** Use a separate
   [credential](access-control/access-control.md#credentials-management-system) per group on one
   connection, or an [OAuth connector](../access/databases/connectors/oauth-connectors.md) that runs
   queries as the signed-in user so the database's own row-level security applies. Use this option
   wherever the data source supports it.
2. **Row-level, in Datagrok** (Beta). With [domain schemas](../develop/how-to/db/domain-schemas.md),
   records such as studies, plates, or compound sets can be shared individually, and detail records
   follow the record they belong to.
3. **Column-level, in Datagrok** (Beta). A user sees a column only if one of their groups has access
   to it. Hidden columns are left out of query results, exports, and filters.

[Managed settings](access-control/managed-settings.md) and interface locks change what people see
in the UI. Don't use them to protect data.

### External collaborators

Give CRO and partner groups **View and use** access only, without permission to create entities or
to **Share With Everyone**. Add them with an
[invitation link](access-control/users-and-groups.md#inviting-users-via-url), and block them when
the collaboration ends.

## Content versions and ownership

Datagrok doesn't keep a version history for
[dashboards](../datagrok/concepts/project/dashboard.md#versioning) and projects. Saving over one
replaces it, including changes a colleague saved in the meantime. Until that changes:

* Give each shared dashboard one owner who edits it. Everyone else runs it or saves a copy.
* Save a copy before a large change.
* To keep a result reproducible, save with data sync off. The data is stored with the layout.
* Save layouts to the gallery, so you can reapply them after a bad edit.

For anything that needs change control, such as NCA, QSAR, curve fitting, a scoring model, or a
critical dashboard, use a package. Packages are [versioned](../develop/develop.md#version-control):
each publish creates a version, an administrator chooses the active version and can roll back, and
a version can be assigned to one group before everyone gets it. See
[publishing packages](../develop/how-to/packages/publish-packages.md) and the
[versioning policy](../develop/dev-process/versioning-policy.md).

## Monitoring

Datagrok records:

* An [audit trail](audit/audit.md) of actions on every object (create, edit, delete, share, run),
  and security events such as sign-ins, failed logins, impersonation, admin sessions, and settings
  changes
* Service health: `/admin/health` returns each service's status without a login, for load balancer
  and Kubernetes probes, and `grok s healthcheck` summarizes it (see
  [server health](../develop/server-management.md#server-health))
* Logs, which can be exported to [CloudWatch](../datagrok/solutions/teams/it/log-export-cloud-watch.md)
  or Google Cloud Logging

:::note Bleeding-edge build

Built-in alerts are available on the bleeding-edge build and not yet in a stable release. They cover
an unhealthy service, one error affecting several users, repeated failed logins for one account,
and an [unreachable data connection](../access/databases/monitor-connections.md).

:::

## Environments

We recommend four environments:

| Environment | Used for | Version | Access and configuration |
|---|---|---|---|
| **Dev / sandbox** | Package development, trying new features | Bleeding edge | Developers only; no production connections |
| **Test** | Checking upgrades before production | Next stable | Same configuration as production. Package database migrations can't be rolled back, so test them here first. |
| **Training** | Onboarding and demos | Stable | Broad access; synthetic data; no production credentials |
| **Production** | Daily work | Stable | Few administrators; approved package versions; configuration from your deployment chart |

Use separate credentials for each group in each environment, and masked or synthetic data outside
production.

### Moving content between environments

Keep these in files and apply them to each server:

* **Packages**: in git, published by CI
* **Server configuration**: [environment parameters](../deploy/configuration.md) and global
  permissions, next to your deployment chart values
* **Users, groups, connections, and sharing**: JSON applied with
  [`grok s`](../develop/server-management.md). Applying the same file twice changes nothing.
* **The default package list**: in your deployment chart

Content created in the UI (connections, queries, scripts, dashboards, spaces, layouts, tables,
files, and models) moves between servers as a bundle with `grok s pull` and `grok s push`.
Dependencies are included, items keep the same IDs on both servers, and you can preview what a push
will change. See [`pull`, `push`, and `migrate`](https://github.com/datagrok-ai/public/blob/master/tools/GROK_S.md).

Set these up on the target server first:

* **Credentials**, which are never copied
* **Users**. Otherwise their content is assigned to the account doing the push.
* **Packages** that own any of the pushed items
* **Platform file shares**, which come from the deployment

## See also

* [Users and groups](access-control/users-and-groups.md)
* [Roles](access-control/roles.md)
* [Access control](access-control/access-control.md)
* [Configure authentication](../deploy/complete-setup/configure-auth.md)
* [Spaces](../datagrok/concepts/project/space.md)
* [Audit](audit/audit.md)
* [Server management](../develop/server-management.md)
