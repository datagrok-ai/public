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

1. [**Who gets in**](#who-gets-in): identity, groups, roles, and, if your license has tiers, which seat each person holds
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

### Independent axes

We recommend modeling users along axes that never mix:

| Axis | Modeled as | Example | Source |
|---|---|---|---|
| **Structure** | Groups | Oncology, a site, Chemistry | Synced from your identity provider |
| **Role** | Roles | Authors, Consumers, Developers | Built in Datagrok |
| **License tier** (only if your license has tiers) | Roles | Full, Standard, Viewer | Built in Datagrok, exactly one per user |

Structure and role describe your organization. Keeping them separate means an organizational
change never silently changes what anyone is allowed to do. The same principle applies to every
capability: grant it on roles, and use groups only for sharing.

### If your license has tiers

With an enterprise license, everyone has the same entitlement and you can skip this section. If
your license defines tiers of users, model the tier as a third axis. License tier is a commercial
fact, and keeping it separate from structure and role means "who is entitled to what" stays
answerable as the organization changes shape underneath.

Datagrok has no built-in license-tier object. A tier is a convention built from roles, and it holds
only if two rules hold:

1. **Every user holds exactly one tier role.** Assign it per person, since it represents a seat.
   No tier role, no seat.
2. **Capabilities are granted on roles and nowhere else.** Structure groups carry content sharing
   only, never **Create Dashboard**, **Create Script**, **Publish Package**, or the **Browse**
   permissions.

:::caution Permissions are additive

A user's effective permissions are the union of every group and role they hold. There is no deny.
A "Viewer" tier role cannot take away a capability that another group or role already granted.
Rule 2 is what makes a tier mean something, and your admin process, not the product, holds that
line.

:::

With those rules in place, utilization becomes three numbers per tier: provisioned (tier-role
members), active in the last N days (from [Usage Analysis](audit/usage-analysis.md)), and dormant.
Reconcile monthly. To reclaim a seat, [block the user](access-control/users-and-groups.md#disabling-accounts):
a blocked user can't sign in and stops counting toward the license, and their work stays in the
system and stays shareable.

Write down the mapping from each tier to its exact permission set. It is the artifact that makes
tiers auditable.

### Identity and group sync

Set up single sign-on with your identity provider, and let
[group synchronization](../deploy/complete-setup/configure-auth.md#group-synchronization) maintain
the structure axis:

* Sync runs at sign-in, so a membership change takes effect at the user's next login. If you need
  changes to land sooner, or want to rebuild a server's users and groups from files, add a
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
  [sticky meta](catalog/sticky-meta.md), therapeutic area, target, modality, and phase become typed,
  searchable properties. Reclassifying a program is then an edit, not a move that silently
  inherits a different parent's permissions.
* **Treat every copied attribute as a cache.** Copy a field onto a space only if people need to find
  spaces by it before any data loads. Read everything else live from your system of record, and keep
  access-controlling fields out of metadata entirely.

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

Build your process around that:

* **Save a copy before a restructuring**, named for its purpose. Share the copy once it's accepted,
  and retire the original.
* **Snapshot anything that must be reproducible** by saving with data sync off, which freezes the
  data together with the layout. A live dashboard can't be reproduced later.
* **Keep layouts in the gallery.** A layout is a separate object, so it survives a bad edit and can
  be reapplied to new data.
* **Give each canonical dashboard one owner.** Assay and program owners hold edit rights; everyone
  else runs it. This is the control that prevents forked copies from drifting apart.

### Put decision-grade logic in packages

Anything a decision rests on, such as NCA, QSAR, curve fitting, or a scoring model, belongs in a
package or a script rather than a hand-edited dashboard. Packages are
[properly versioned](../develop/develop.md#version-control): every publish creates a version record,
several versions can be deployed at once, an administrator chooses the active one and can roll
back, and a specific version can be assigned to a specific group to pilot an upgrade. Publish from
git or from CI (see [publishing packages](../develop/how-to/packages/publish-packages.md) and the
[versioning policy](../develop/dev-process/versioning-policy.md)).

## How you know it's healthy

The platform already records:

* **An [audit trail](audit/audit.md)** of every action on every object (created, edited, deleted,
  shared, executed), structured and queryable
* **A security trail**: sign-ins and sign-outs, failed logins with reason, impersonation, admin
  sessions, key generation, and settings changes
* **[Usage Analysis](audit/usage-analysis.md)**: active and new users, package and function usage,
  per-project access frequency, and function timings
* **Health checks**: `/admin/health` lists each service's status and needs no login, for load
  balancer and Kubernetes probes; `grok s healthcheck` summarizes it for operators (see
  [server health](../develop/server-management.md#server-health))
* **Log export** to CloudWatch or Google Cloud Logging, by record type
  (see [log export](access-control/data-connection-credentials.md#for-logs-export-to-cloudwatch))
* **Per-group log verbosity**, so you can turn up detail for one team without flooding the rest

We recommend alerting on problems rather than watching dashboards of metrics:

* A service goes unhealthy or stops reporting
* One error reaches several users at once, which usually means a bad package upgrade
* A route or connection slows down or starts failing
* Repeated failed logins, unexpected admin sessions, or settings changes

Separately, review active versus provisioned users monthly. If your license has tiers, it's the
same number that drives your tier reconciliation.

For step-by-step guidance, see [Monitor data connections](../access/databases/monitor-connections.md)
and [Test dashboards and queries](../develop/how-to/tests/test-content.md).

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
