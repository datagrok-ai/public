---
title: "Manage an enterprise instance"
sidebar_label: "Manage enterprise"
sidebar_position: 0
description: A guide for administrators running Datagrok for a large organization, from identity and access to data connections, dashboards, monitoring, and upgrades.
keywords:
  - enterprise administration
  - groups vs roles
  - global permissions
  - monitoring data connections
  - health check
  - dashboard version control
  - dev to prod promotion
  - usage analysis
  - log export
---

This page follows an enterprise instance through its life: from the day IT
deploys it, through onboarding thousands of users and hundreds of data sources,
to monitoring it in production and upgrading it. Each section tells you what to
decide and links to the page that shows how.

| You want to...                                  | Go to                                                        |
|-------------------------------------------------|--------------------------------------------------------------|
| Connect your identity provider                  | [Identity](#identity)                                        |
| Understand groups vs. roles                     | [Groups and roles](#groups-and-roles)                        |
| Decide who can create and share what            | [Global permissions](#global-permissions)                    |
| Govern database access and credentials          | [Data connections](#data-connections)                        |
| Know when a connection breaks                   | [Monitor connections](#monitor-connections)                  |
| Keep dashboards safe from bad edits             | [Dashboards](#dashboards)                                    |
| Move content from dev to prod                   | [Promote content](#promote-content)                          |
| Feed Datagrok into your monitoring stack        | [Monitor the platform](#monitor-the-platform)                |
| Audit who did what                              | [Audit](#audit)                                              |
| Plan upgrades and backups                       | [Operate](#operate)                                          |

## Set up

A typical enterprise runs Datagrok on-premises or in its own cloud account, on
[Kubernetes or Docker Compose](../deploy/deploy.md). Plan for at least two
instances: a **dev/test** instance where developers publish packages and
analysts try new dashboards, and a **production** instance that users rely on.
Many organizations add a third, **validation**, instance for upgrades.

After the deployment, sign in as `admin` and complete the setup:

1. [Configure authentication](../deploy/complete-setup/configure-auth.md) and
   verify it in an incognito window before turning off login/password sign-in.
1. [Configure email](../deploy/complete-setup/configure-smtp.md) so that sharing,
   notifications, and password resets work.
1. [Install the packages](../deploy/complete-setup/install-packages.md) your
   users need.
1. Review the **Health** step of the setup wizard (`/settings/initial/health`).
   It lists the checks the server runs continuously (Core, Grok Connect,
   Jupyter, Credentials Server, Grok Spawner, and others). Fix anything that
   isn't ready before you invite users.
1. Set [defaults for everyone](access-control/managed-settings.md), such as
   formats, colors, and hidden menu sections, and lock the ones users
   shouldn't change.
1. Turn on the garbage collector in [server settings](access-control/admin-settings.md#data-retention).
   It's off by default, so audit records, logs, and old package versions are
   never deleted and the database keeps growing.

The server settings pages are listed in
[Server settings](access-control/admin-settings.md). Server-side parameters,
such as queue, Docker, and connector settings, are described in
[Server configuration](../deploy/configuration.md). You can set them in the UI
or pass them as `GROK_PARAMETERS` to keep them in your infrastructure code.

## Identity

Users sign in through your identity provider. Datagrok supports
[OpenID Connect](../deploy/complete-setup/configure-auth.md#openid-authentication)
(for example, Microsoft Entra ID or Okta),
[SAML](../deploy/complete-setup/configure-auth.md#saml-authentication),
[LDAP](../deploy/complete-setup/configure-auth.md#ldap-authentication), and
[Google IAP](../deploy/complete-setup/configure-auth.md#iap-authentication).
A user account is created on first sign-in.

Let the identity provider own group membership, and let Datagrok own what those
groups are allowed to do:

* With OpenID, Datagrok reads the `groups` claim at every sign-in and
  [synchronizes group membership](../deploy/complete-setup/configure-auth.md#group-synchronization).
  Missing groups are created automatically. Groups and memberships you added by
  hand in Datagrok are never removed by the sync.
* With SAML or LDAP, membership isn't synchronized automatically. Use the
  [`grok s` CLI](../develop/server-management.md#sync-an-ad-group-with-datagrok)
  to sync groups on a schedule from your directory.

For people who leave, [disable the account](access-control/users-and-groups.md#disabling-accounts).
Their sessions end immediately, they stop counting toward the license, and
their dashboards and queries stay available. The **Disable user** dialog lets
you move what they owned to a space.

For automation, create a
[service user](access-control/users-and-groups.md#adding-users) and give it a
[key pair](access-control/keypair-authentication.md#ci-and-automation) rather
than a person's credentials.

## Groups and roles

Datagrok has one mechanism for granting access: **permissions are granted to
groups**. Users, groups, and roles are all ways of putting people into groups.

| Concept            | What it is                                                           | Typical source                  | Example                  |
|--------------------|----------------------------------------------------------------------|---------------------------------|--------------------------|
| **Personal group** | Created automatically for every user. Sharing "with a user" grants to this group | Automatic              | `jdoe`                   |
| **Group**          | A set of users and other groups. Answers "*who* are these people?"   | Your identity provider          | `Oncology Discovery`     |
| **Role**           | A group marked as a role. Holds permissions and is assigned to groups. Answers "*what* may they do?" | Created by an administrator | `Data Steward`, `Dashboard Author` |

A role is technically a group with a role flag. It nests and inherits exactly
like a group does. The difference is how you use it:

* **Groups collect people.** Let them mirror your organization and come from
  your identity provider.
* **Roles collect permissions.** Grant global permissions and entity
  permissions to a role, then **assign the role to groups**. Every member of the
  group inherits what the role grants.

This separation keeps the identity provider in charge of *who* and Datagrok in
charge of *what*. Group synchronization never matches roles, so a group created
in the identity provider can't grant itself a Datagrok role. A role is assigned
by an administrator or by one of the role's own admins, the members marked
**Can assign**.

Permissions flow down the hierarchy: members of a child group get everything
granted to the parent. Assigning a role to a group makes that group a member of
the role.

To manage roles, go to **Browse** > **Platform** > **Roles**. To assign a role
to a group, open the group's editor and use the **Roles** tab. A role's
**Assigned to** tab shows which groups have it.

<details>
<summary>Built-in groups and roles</summary>

| Name               | Kind                                  | Purpose                                                               |
|--------------------|---------------------------------------|-----------------------------------------------------------------------|
| **All users**      | Group                                 | Contains every user and group. Its permissions apply to everyone      |
| **Administrators** | Role                                  | Holds all global permissions. Keep at least one working member         |
| **Admin**          | Personal group of the `admin` user    | Member of Administrators                                               |
| **System**         | Internal service account              | Used by the platform itself. Don't modify                              |

</details>

To learn more, see [Users and groups](access-control/users-and-groups.md) and
[Access control](access-control/access-control.md#authorization).

## Global permissions

Global permissions decide what people can do across the platform, as opposed to
what they can do with a particular dashboard or connection. To edit them, go to
**Settings** > **Global Permissions**, or select a group or role and use
**Global Permissions** on the **Context Panel**. See the
[full list](access-control/access-control.md#global-permissions).

A new instance is open by default: **All users** can create connections,
queries, scripts, dashboards, and spaces, invite users, and share with anyone.
This suits a small team. For an enterprise, review these defaults before the
rollout. A common setup:

| Role               | Global permissions                                                                 |
|--------------------|------------------------------------------------------------------------------------|
| All users          | **Create Dashboard**, **Create Space**, **Browse** permissions for what they use    |
| Dashboard Author   | adds **Create Data Query**, **Share With Everyone**                                |
| Data Steward       | adds **Create Database Connection**, **Create File Connection**, **Create Security Connection** |
| Developer          | adds **Create Script**                                                             |
| Administrators     | everything                                                                         |

The **Browse** permissions only hide sections of the **Browse** tree, which
simplifies the UI for business users. They don't restrict access. Access is
controlled by entity permissions.

## Content and spaces

Everything users create, such as dashboards, queries, scripts, and files, is an
entity with its own [permissions](access-control/access-control.md#permissions):
**View**, **Edit**, **Delete**, and **Share**, plus type-specific ones such as
**Execute** for queries. New content is private to its author until shared.

Organize shared content in [spaces](../datagrok/concepts/project/space.md), one
per team or project. Grant a team's group access to its space once, and
everything moved into the space inherits those permissions. Give each space an
owner group responsible for keeping it tidy.

## Data connections

Connections are where enterprise governance matters most, because they reach
your databases.

* **Who can create connections.** Restrict **Create Database Connection** to
  data stewards (see [Global permissions](#global-permissions)).
* **Who can use them.** Share a connection with **View and use** to let a group
  run its queries. A group can open a dashboard built on a query without being
  able to see or edit the connection behind it.
* **Who can write.** Write access (adding, changing, and removing rows) and
  schema changes are separate permissions on the connection. Grant them only to
  the groups that need them.
* **Whose credentials.** Credentials are
  [stored encrypted](access-control/access-control.md#credentials-storage),
  separately from the rest of the metadata, and can differ by group: for
  example, a read-only database account for **All users** and a read-write one
  for data stewards. To keep secrets in your vault, use
  [AWS Secrets Manager or GCP Secret Manager](access-control/data-connection-credentials.md).
  For databases that support it, use
  [OAuth](../access/databases/connectors/oauth-connectors.md) so that each user
  queries under their own identity.
* **How fast.** Turn on [result caching](../access/databases/databases.md#refreshing-and-caching)
  for slow queries whose data changes on a known schedule.

To learn more, see [Databases](../access/databases/databases.md).

### Monitor connections

Datagrok doesn't probe every connection on a schedule. Use these tools to know
when one breaks:

| Signal                          | Where                                                                                        |
|---------------------------------|----------------------------------------------------------------------------------------------|
| Test a connection now           | Right-click the connection > **Test connection**, or `grok s connections test <id>` from your monitoring job |
| Grok Connect (the connector service) is up The [health endpoints](../develop/server-management.md#server-health) or `grok s healthcheck` |
| Who ran which query, when, how long it took, and whether it failed | The query's **Activity** and **Usage** panes on the **Context Panel** ([audit](audit/audit.md#accessing-audit-logs)) |
| Failing queries across the platform | **Usage Analysis** > **Errors** and **Functions** ([Usage Analysis](audit/usage-analysis.md)) |
| Why a query is slow             | The **Debug** tab of the [query editor](../access/databases/databases.md#query-editor)       |

To actively watch critical connections, create a small test query for each and
[schedule it](../datagrok/concepts/functions/functions.md#scheduling), or run
`grok s connections test` from your existing monitoring system and alert on
failure. Failed runs are recorded as errors and appear in Usage Analysis and in
[exported logs](#monitor-the-platform).

## Dashboards

A [dashboard](../datagrok/concepts/project/dashboard.md) is a saved project:
tables, viewers, and layout, optionally re-running its queries every time it
opens ([Data sync](../datagrok/concepts/project/dashboard.md#data-sync)).

### Version control

Datagrok doesn't keep a version history for dashboards created in the UI.
Saving over the original replaces it for everyone. Set expectations with
authors and follow the
[versioning practices](../datagrok/concepts/project/dashboard.md#versioning):
save a copy before a risky change, save reproducible reports with Data sync off,
and keep important layouts in the gallery.

For dashboards that must be versioned, reviewed, and released like software,
ship them in a [package](../develop/develop.md#packages). A package holds
connections, queries, scripts, and dashboards in a Git repository and is
published as a numbered version. You can switch back to a previous version when
a release goes wrong (see [Publishing](../develop/develop.md#publishing)).

### Quality checks

There is no built-in check that opens every dashboard after a change. The
[Activity pane and the audit log](audit/audit.md) show when a dashboard was
edited and by whom, and **Usage Analysis** > **Errors** shows failures users hit
while opening it. When a query's output changes,
[this table](../datagrok/concepts/project/dashboard.md#when-the-source-changes)
shows what breaks.

For business-critical dashboards, test what they depend on:

* **Queries and scripts in packages** can declare test cases in their
  annotations. They run with the package tests in the
  [Test Manager](../develop/how-to/tests/test-packages.md#test-manager), from the
  command line with `grok test`, or in your CI
  (see [Testing functions](../develop/how-to/tests/add-package-tests.md#testing-functions)).
* **Package tests** can open a project and check its tables and viewers
  (see [Adding unit tests](../develop/how-to/tests/add-package-tests.md#adding-unit-tests)).
* **After an upgrade**, run the package tests on the validation instance before
  upgrading production (see [Operate](#operate)).

## Promote content

Content moves from dev to production in one of two ways:

* **As a package.** Developers keep connections, queries, scripts, and
  dashboards in Git and publish the same package version to each instance with
  `grok publish <instance> --release`. Connection settings that differ between
  instances, such as hosts and credentials, are
  [placeholders](../develop/develop.md#connections) filled in on each server.
  This is the recommended path for anything business-critical.
* **As entities.** For content authored in the UI, the `grok s` CLI can pull
  entities from one instance into a folder of files, show the difference, and
  push them to another instance, keeping their IDs so a repeated push updates
  rather than duplicates. Credentials never travel with the content. See the
  [Move content between instances](../develop/server-management.md#move-content-between-instances).

## Monitor the platform

Datagrok doesn't expose a Prometheus `/metrics` endpoint. Instead, it offers
health endpoints to probe, an in-app metrics dashboard, and a log stream to push
to your observability system.

| What                    | How                                                                                           | Use it for                                  |
|-------------------------|-----------------------------------------------------------------------------------------------|---------------------------------------------|
| Liveness                | `GET https://<host>/api/admin/health`. No sign-in needed. Returns every service with its status | Load balancer health checks, uptime probes |
| Detailed health         | `GET https://<host>/api/public/v1/healthcheck` with an API token. Returns `status: ok` or `degraded`. Always answers HTTP 200, so read `status` from the body | Synthetic monitoring   |
| Health from the shell   | `grok s healthcheck` ([Server health](../develop/server-management.md#server-health))        | Scripts and runbooks                        |
| Operational metrics     | **Usage Analysis** > **Metrics** (administrators): database health and size, storage and disk space, request latency, function-call queue, errors, sessions, slowest database statements | Capacity planning, triage |
| Usage                   | [Usage Analysis](audit/usage-analysis.md): active users, packages, functions, projects, errors | Adoption and value reporting              |
| Logs, alerts, heartbeat | **Settings** > **Logger** > **Log sync**: push to Amazon CloudWatch, Google Cloud Logging, or any OpenTelemetry (OTLP) collector ([setup](../datagrok/solutions/teams/it/log-export-cloud-watch.md)) | Central logging, alerting, SIEM |

The log stream includes audit events, errors, and alerts. The server opens an
alert when a health check fails, when requests slow down or fail across the
board, when one error hits many users, when an account keeps failing to sign
in, or when a user reports a problem, and resolves it when the condition
clears. It also sends a heartbeat every five minutes. Alert on a
missing heartbeat in your collector to detect an instance that stopped
reporting. Container-level metrics, such as CPU and memory, come from your
platform (Kubernetes, ECS, or Docker) as for any other service.

## Audit

Every server action and every explicitly run function is recorded with the
user, time, and parameters, in regular Postgres tables you can query directly.
Sign-ins, failed sign-ins, impersonation, admin sessions, and settings changes
form a separate security trail on **Usage Analysis** > **System Activity**. You
can raise or lower the level of detail per group, for example verbose logging
for a regulated team. See [Audit](audit/audit.md) and
[Data provenance](audit/data-provenance.md).

Users can report errors directly from the platform. Configure where reports go
in [Feedback](bug-reports.md#configuring-error-reporting-system).

## Operate

* **Upgrades.** Read the [release notes](../deploy/releases/releases.md) and the
  [compatible service versions](../deploy/releases/compatibility/compatibility.mdx).
  Upgrade the validation instance first, run your package tests and a smoke test
  of key dashboards, then upgrade production. The database schema migrates
  automatically on the first start of the new version, and migrations only run
  forward: rolling back means restoring the backup taken before the upgrade.
  See [Upgrade and roll back](../deploy/upgrade.md).
* **Backups.** Back up the Postgres database (metadata, audit, and encrypted
  credentials), the file storage, and the server keys together, so that they
  can be restored to the same point in time. A database restored without its
  keys loses every stored credential. See
  [Back up and restore](../deploy/complete-setup/backup.md).
* **Keys.** Rotate the [server keys](access-control/server-keys.md) that
  encrypt credentials according to your security policy.
* **Security posture.** Review [Security](../datagrok/solutions/teams/it/security.md)
  for vulnerability scanning and published VEX reports.

## See also

* [Access control](access-control/access-control.md)
* [Users and groups](access-control/users-and-groups.md)
* [Managed settings](access-control/managed-settings.md)
* [Server settings](access-control/admin-settings.md)
* [Server management with grok s](../develop/server-management.md)
* [Enterprise evaluation FAQ](../datagrok/solutions/teams/it/enterprise-evaluation-faq.md)
