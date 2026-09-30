---
title: "Running Datagrok in the enterprise"
sidebar_position: 0
description: Recommended defaults for governing a large Datagrok deployment - identity and license tiers, data isolation, content lifecycle, observability, and environments.
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

A Datagrok deployment at enterprise scale comes down to five decisions:

1. [**Who gets in**](#who-gets-in): identity, groups, roles, and, if your license has tiers, seat counting
2. [**What they can see**](#what-they-can-see): isolation between programs, partners, and data classes
3. [**What they can change**](#what-they-can-change): content ownership and lifecycle
4. [**How you know it's healthy**](#how-you-know-its-healthy): audit, usage, and operational signal
5. [**Where work happens**](#where-work-happens): dev, test, training, production, and what moves between them

Decisions 1, 2, and 5 are structural. They are cheap to make before users start creating content and
expensive to change afterward, so settle them first. Decisions 3 and 4 can evolve as usage grows.

This page gives a recommended default for each decision, says plainly where the platform enforces it
and where your admin process has to, and links to the reference pages for the mechanics.

## Who gets in

### Groups for structure, roles for capabilities

Datagrok grants permissions to [groups](access-control/users-and-groups.md#groups) and
[roles](access-control/roles.md), never to individuals directly. Use each for one job:

* **Groups** mirror your organization and carry sharing: who can see and edit a space, a query, a
  connection, or a dashboard.
* **Roles** carry capabilities: [global permissions](access-control/access-control.md#global-permissions)
  such as **Create Dashboard** or **Publish Package**, and access to functions and packages.

Both can be nested and assigned to each other, and children inherit their parent's permissions.
Keep nesting to two or three levels; beyond that, "who can see this?" stops having a checkable
answer.

A new entity is visible only to its author until it is shared, which is the safe default everything
else builds on. Every user also gets an automatic personal group. Don't share with personal groups:
access granted that way disappears when the person changes role, and nobody notices until someone
can't open their program's data.

### An example

Priya is a medicinal chemist in Oncology.

* **Her groups** come from your identity provider: _Oncology_ and _PRAME program_. They decide what
  she can see: the PRAME space and its dashboards, and the Oncology credential on the shared
  warehouse connection.
* **Her role** is assigned in Datagrok: _Authors_. It decides what she can do: create dashboards and
  queries and share them.

When she transfers to Immunology, your identity provider moves her groups. At her next login she
loses the PRAME data and gains Immunology's. Her role doesn't change, so she can still author, and
no Datagrok administrator is involved.

When she becomes the assay owner for her team, an administrator adds the _Content Owners_ role. She
can now publish the team's canonical queries, and she sees no new data.

An org change moves groups, a job change moves roles, and neither changes the other by accident.
Grant every capability on roles, and use groups only for sharing.

### License tiers

With an enterprise license, everyone has the same entitlement and you can skip this section.

If your license defines tiers of users, keep the tier assignment where you keep other entitlements:
in your identity provider, as one group per tier, for example `DG-Tier-Full` and `DG-Tier-Viewer`.
In Datagrok, the tier groups are for counting seats only. They carry no permissions, so assigning a
tier never changes what anyone can do.

Seat utilization is then three numbers per tier:

* **Provisioned**: members of the tier group in your identity provider
* **Active**: provisioned users who signed in during the last 30 days, from the sign-in records in
  the [audit trail](audit/audit.md)
* **Dormant**: the difference. To reclaim a seat, remove the person from the tier group and
  [block the user](access-control/users-and-groups.md#disabling-accounts). A blocked user can't sign
  in and stops counting toward the license, and their work stays in the system.

How the tier groups get into Datagrok matters for the count:

* **Login-time [group synchronization](../deploy/complete-setup/configure-auth.md#group-synchronization)**
  updates a person's groups only when they sign in. Someone who has never signed in isn't in the
  group yet, and someone who left the tier, or the company, stays in it because they never sign in
  again. That's fine for structure groups, but it makes the tier count drift upward.
* **A scheduled reconciliation job** reads the tier groups from your directory and applies them with
  `grok s groups add-members` and `remove-members` (see
  [Sync an AD group with Datagrok](../develop/server-management.md#sync-an-ad-group-with-datagrok)).
  Membership then matches the directory on every run, including people who haven't signed in, and
  the same job can block people who have left. We can provide this job for your directory.

### Identity and group sync

Set up single sign-on with your identity provider, and let
[group synchronization](../deploy/complete-setup/configure-auth.md#group-synchronization) maintain
your structure groups:

* Sync runs at sign-in, so a membership change takes effect at the user's next login. If you need
  changes to land sooner, or need accurate counts (see [License tiers](#license-tiers)), add a
  [scheduled reconciliation job](../develop/server-management.md#sync-an-ad-group-with-datagrok).
* Groups and memberships an administrator created by hand are never removed by sync, and roles are
  never matched. Sync cannot lock you out.
* Sync matches groups by name. If a group you created by hand has the same name as one your
  identity provider asserts, sync adopts it: it adds the provider's members and never removes the
  ones you added. Give hand-built groups a naming convention so that only happens on purpose.
* With Microsoft Entra ID, every group a user belongs to becomes a Datagrok group, and the Graph
  permissions it needs require tenant admin consent. In a large organization, that consent is
  usually a change request, so start it early.

### Administration

Keep the `Administrators` role very small. Delegate day-to-day membership work to group admins, who
can add and remove members and approve requests without touching global permissions.

Never use an administrator account to check whether a restriction works. Administrators see
everything, so verify every access boundary with a second, ordinary account.

## What they can see

### Spaces are the unit of sharing

A [space](../datagrok/concepts/project/space.md) is a container with permissions. Root spaces hold
child spaces, children inherit the root's privileges, and moving an entity into a space makes it
adopt that space's permissions. To show one canonical dashboard in several places, link it rather
than copying it: a link is a live, view-only reference to the original.

The most common failure state is an estate governed by hundreds of individual sharing decisions:
everything a user creates lands in their personal space until someone moves it. It's recoverable,
because entities move between spaces with drag and drop, but it's much cheaper to set up the space
structure first.

### Hierarchy is permissions; everything else is metadata

Because child spaces inherit their parent's privileges, every level you add to the space tree is an
access boundary, whether you meant it to be one or not.

* **Add a level only where the access boundary genuinely differs.** A program or a partnership is an
  access boundary. A therapeutic area or a modality usually isn't; it's a way of looking at the
  portfolio.
* **Put classifications on the space as metadata instead.** With
  [sticky meta](catalog/sticky-meta.md), a space can carry properties such as therapeutic area,
  target, modality, and phase, and people can search and filter spaces by them. If a program moves
  to a different therapeutic area, you edit the property. Moving the space instead would change who
  can see it.
* **Keep metadata small.** Add a property to a space only if people search for spaces by it. Other
  program details belong in your system of record, and nothing that controls access should be
  stored as metadata.

### Provision program spaces from your system of record

If your organization keeps programs or projects in a database, create spaces from it with a
scheduled job instead of by hand. The job should:

* Key each space on the program's stable ID, not its name, so a rename updates a label instead of
  orphaning a folder
* Grant the owning group edit rights and the reader group view rights, reconciling grants on every
  run from the live access flags
* Archive closed programs rather than deleting them, because deletion removes content from everyone
  irreversibly
* Run incrementally and idempotently, with a dry-run mode that lists what it would change

### Three layers of data isolation, strongest first

1. **The database enforces it.** Either bind several
   [credentials](access-control/access-control.md#credentials-management-system) to one connection,
   each for a group or a user, or pass each user's own identity through with an
   [OAuth connector](../access/databases/connectors/oauth-connectors.md) so the database's own
   row-level security applies. Either way, your database decides. Prefer this layer wherever the
   source supports it; it's the one that survives an audit.
2. **Row-level security inside the platform** (Beta). With
   [domain schemas](../develop/how-to/db/domain-schemas.md), registration records such as studies,
   plates, or compound sets can be shared row by row, and detail records inherit from the row they
   reference.
3. **Column-level security inside the platform** (Beta). A user sees a column only if one of their
   groups has access to it. Hidden columns never leave the server, so they're absent from results,
   exports, and filters.

[Managed settings](access-control/managed-settings.md) and interface locks are not a security
boundary. They shape what people see and keep teams consistent. Don't use them to protect anything
sensitive.

### External collaborators

Give CRO and collaborator groups narrow **View and use** access, no entity creation, and no
**Share With Everyone**. Onboard them by
[invitation link](access-control/users-and-groups.md#inviting-users-via-url) and offboard them by
blocking.

## What they can change

:::caution No version history for dashboards

Datagrok keeps no version history for [dashboards](../datagrok/concepts/project/dashboard.md#versioning)
and projects. Saving over one replaces the previous state, including changes a colleague saved in
the meantime.

:::

These practices cover most needs:

* **Save a copy before restructuring a dashboard**, named for its purpose. Share the copy once it's
  accepted, and retire the original.
* **Snapshot anything that must be reproducible** by saving with data sync off, which freezes the
  data together with the layout. A live dashboard can't be reproduced later.
* **Keep layouts in the gallery.** A layout is a separate object, so it survives a bad edit and can
  be reapplied to new data.
* **Give each canonical dashboard one owner.** Assay and program owners hold edit rights; everyone
  else runs it. This is the control that prevents forked copies from drifting apart.

### Put critical content in packages

A dashboard or query that many people rely on can be shipped in a package, and so can anything a
decision rests on, such as NCA, QSAR, curve fitting, or a scoring model. Packages are
[versioned](../develop/develop.md#version-control): every publish creates a version record,
several versions can be deployed at once, an administrator chooses the active one and can roll
back, and a specific version can be assigned to a specific group to pilot an upgrade. Publish from
git or from CI (see [publishing packages](../develop/how-to/packages/publish-packages.md) and the
[versioning policy](../develop/dev-process/versioning-policy.md)).

## How you know it's healthy

The platform records:

* **An [audit trail](audit/audit.md)** of every action on every object (created, edited, deleted,
  shared, executed), and a security trail of sign-ins, failed logins, impersonation, admin sessions,
  and settings changes
* **Health checks**: `/admin/health` lists each service's status and needs no login, for load
  balancer and Kubernetes probes; `grok s healthcheck` summarizes it for operators (see
  [server health](../develop/server-management.md#server-health))
* **Log export** to CloudWatch or Google Cloud Logging
  (see [log export](access-control/data-connection-credentials.md#for-logs-export-to-cloudwatch))

The next release adds built-in alerts, so you don't have to build them yourself:

* A service becomes unhealthy
* One error reaches several users at once, which usually means a bad package upgrade
* Repeated failed logins for one account
* A [data connection becomes unreachable](../access/databases/monitor-connections.md)

We're extending the same alerts to failing [query tests](../develop/how-to/tests/test-content.md)
and failed scheduled runs. If you need to know about a problem that isn't on this list, tell us.

## Where work happens

### Four servers, each with a stated job

| Environment | Purpose | Version | Posture |
|---|---|---|---|
| **Dev / sandbox** | Package development, developer builds, schema experiments, trying new features | Latest | Developers publish here; small population; no production connections |
| **Test / validation** | Release builds, mirroring production configuration | Next stable, before production | Where upgrades and schema migrations are verified. Not optional: package database migrations don't roll back |
| **Training** | Onboarding, tutorials, demo and synthetic data | Stable | Wide access and permissive; no production credentials |
| **Production** | The real thing | Stable | Very small admin group, approved package versions only, configuration from your deployment chart rather than the UI |

Use separate credentials per group per environment, and masked or synthetic data everywhere except
production.

### What moves between environments

**Rebuilt from files:**

* **Packages**: versions in git, published by CI, installed and upgraded by command
* **Server configuration**: [environment parameters](../deploy/configuration.md) and the global
  permission set, checked in beside your deployment chart values
* **Users, groups, connections, and sharing**: declared as JSON and applied idempotently with
  [`grok s`](../develop/server-management.md)
* **The default package set**: in your deployment chart, so every server starts identical

**Promoted as a bundle.** Content built in the UI (connections, queries, scripts, dashboards,
spaces, layouts, tables, files, and models) moves between servers as a bundle: a directory of JSON,
one file per entity, that you can review and move like any other code. Dependencies come along
automatically, and entities keep the same identity on both servers, so a repeated push changes
nothing. Preview what a push would change before it writes. See
[`grok s pull`, `push`, and `migrate`](https://github.com/datagrok-ai/public/blob/master/tools/GROK_S.md).

Four things belong to the target server rather than the bundle. Prepare them first:

* **Credentials** never travel. Set them on the target.
* **Users** aren't included. Create them first, or their content lands under the account doing the
  push.
* **Package-owned entities** need the package published on the target.
* **Platform shares** come from the deployment itself.

## See also

* [Users and groups](access-control/users-and-groups.md)
* [Roles](access-control/roles.md)
* [Access control](access-control/access-control.md)
* [Configure authentication](../deploy/complete-setup/configure-auth.md)
* [Spaces](../datagrok/concepts/project/space.md)
* [Audit](audit/audit.md)
* [Server management](../develop/server-management.md)
