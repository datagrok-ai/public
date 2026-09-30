---
title: "Server settings"
mdx:
  format: mdx
sidebar_position: 5
description: The platform-wide settings pages an administrator manages, what each controls, and where each is documented.
keywords:
  - server settings
  - admin settings
  - data retention
  - garbage collector
  - maintenance window
  - GROK_PARAMETERS
---

Datagrok settings come in two kinds:

* **Local** settings shape the user interface for one user, such as formats,
  colors, and hidden menus. Every user edits their own. Administrators can set
  defaults and lock them for groups (see [Managed settings](managed-settings.md)).
* **Server** settings apply to the whole instance, such as authentication,
  email, logging, and data retention. Only users with the
  **Edit Plugins Settings** [global permission](access-control.md#global-permissions)
  see and edit them.

To open them, on the **Sidebar**, click **Browse** (<FAIcon icon="fa-solid fa-compass"/>) > **Platform** > **Settings**.
Local settings are under **Local**, and server settings are under **Server**.

## Server pages

| Page                   | What it controls                                                                                              | Learn more                                                              |
|------------------------|---------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|
| **Admin**              | **Endpoints**: the public URL of the instance, used in links and invitations                                   | —                                                                        |
|                        | **E-mail**: SMTP, Mailgun, Amazon SES, or the Datagrok mailer                                                  | [Configure SMTP](../../deploy/complete-setup/configure-smtp.md)          |
|                        | **Error reporting**: where user error reports go, and whether errors are reported automatically               | [Feedback](../bug-reports.md#configuring-error-reporting-system)         |
|                        | **AI Providers**: the model provider behind AI features                                                         | [Configure AI providers](../../explore/ai/ai.md#configure-ai-providers)  |
|                        | **Garbage Collector** and **Maintenance**: data retention and the nightly cleanup window                        | [Data retention](#data-retention)                                        |
| **Users and Sessions** | Sign-up and sign-in rules, LDAP, OpenID, SAML, Google IAP, and group synchronization                           | [Configure authentication](../../deploy/complete-setup/configure-auth.md) |
| **Logger**             | Which events are logged for whom, and **Log sync** to Amazon CloudWatch, Google Cloud Logging, or OpenTelemetry | [Audit](../audit/audit.md), [Export logs](../../datagrok/solutions/teams/it/log-export-cloud-watch.md) |
| **Scripting**          | The script languages enabled on the instance                                                                   | [Scripting](../../compute/scripting/scripting.mdx)                       |
| **Cache**              | Server-side and client-side caching of function results and files, on or off for the whole instance            | [Cache function results](../../develop/how-to/functions/cache-function-results.md) |
| **Global Permissions** | What each group or role may do across the platform. Requires **Edit Global Permissions**                        | [Global permissions](access-control.md#global-permissions)               |

Two more nodes sit next to these pages under **Settings**:

* **Keys**: the server keys that encrypt stored credentials and sign sign-in
  tokens. Requires **Admin Keys**. See [Server keys](server-keys.md).
* **Repositories**: the package registries the server installs plugins from.
  Requires **Browse Plugins**.

## Data retention

The server deletes old log records, sessions, test runs, and orphaned data in a
nightly cleanup, the **garbage collector**. It is **off by default**, so on a
new instance nothing is deleted and the database keeps growing. Turn it on in
**Admin** > **Garbage Collector** > **Enabled**.

| Setting                                  | Default          |
|------------------------------------------|------------------|
| Audit records (the oldest record of each entity is kept) | 365 days |
| Sessions                                 | 365 days         |
| Errors, warnings, info, and usage logs   | 183 days         |
| Unused versions of locally published packages that are not current | 183 days |
| Test runs                                | 90 days          |
| Debug logs, HTTP request statistics, empty sessions, expired sign-in codes | 30 days |
| Orphaned tables and script runs          | Seven days       |

The cleanup starts on the **Maintenance** schedule (cron, UTC, default
`0 0 * * *`, midnight every day) and runs for at most the maintenance duration
(default three hours). The **Garbage Collector** entry of the health check
reports errors from the last run.

If your policy requires keeping audit records longer than the database does,
export them with **Log sync** and keep them in your log archive.

## Set settings in code

Server settings can be part of your infrastructure code instead of the UI.
On any server settings page, click `{}` to get the settings as JSON, then put
them in the `settings` section of `GROK_PARAMETERS`. See
[Server configuration](../../deploy/configuration.md#settings).

On instances managed by an external admin console, such as Datagrok SaaS,
some server pages are read-only and show **Managed by your administrator**.

See also:

* [Managed settings](managed-settings.md)
* [Access control](access-control.md)
* [Server configuration](../../deploy/configuration.md)
