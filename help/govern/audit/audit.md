---
title: "Audit"
mdx:
  format: mdx
sidebar_position: 2
toc_max_heading_level: 3
unlisted: false
description: How Datagrok records, stores, and exports audit logs of user and system activity across entities.
keywords:
  - audit log
  - activity log
  - track usage
  - event logging
  - cloudwatch export
  - audit record types
  - entity history
---

In Datagrok, users can perform a wide range of tasks, from self-service data analytics to using custom applications. Datagrok's audit system provides detailed insights into how
users interact with the platform. With this knowledge, you can:

* Improve UX
  * Remove unused menu items, etc.
* Understand your users better
  * Cluster usage by usage patterns
  * Correlate usage patterns with the organizational groups
* Track data history
  * Track origin of datasets
  * See history of modifications
  * See all actions performed on an [object](../../datagrok/concepts/objects.md)
* Analyze performance
  * See how long each operation takes
* Analyze function usage
  * Get a table of users, timestamps, inputs and outputs
* Track errors
  * See all errors across the platform
  * Identify actions that cause errors

## Audit log storage

Datagrok automatically records and stores all user actions performed on
[entities](../../datagrok/concepts/objects.md) in a structured format within a
[Postgres database](../../develop/under-the-hood/architecture.md#data-engine),
making it easy to query and analyze specific events. The audit data is stored in
the following tables:

* `events`
* `event_types`
* `event_parameters`
* `event_parameter_values`

Each event is associated with a fixed type and the user session that triggered it.

<details>
<summary> Audit record types</summary>

* query-created
* query-edited
* query-deleted
* query-start
* query-published
* query-transformations-edited
* connection-created
* connection-edited
* connection-deleted
* connection-published
* job-created
* job-edited
* job-deleted
* job-transformations-edited
* job-start
* script-created
* script-edited
* script-deleted
* script-start
* script-published
* predictive-model-created
* predictive-model-edited
* predictive-model-deleted
* predictive-model-start
* predictive-model-published
* action-start
* notebook-created
* notebook-edited
* notebook-deleted
* notebook-start
* notebook-published
* project-created
* project-edited
* project-deleted
* project-opened
* entity-shared
* entity-shared-silent
* user-invited
* comment-posted
* dialog-ok
* main-menu-item-click
* error
* tutorial-completed
* viewer-rendered
* log
* package-tested
* alert-opened
* alert-escalated
* alert-resolved
* server-started
* user-logged-in
* user-logged-out
* user-login-failed
* user-impersonated
* impersonation-failed
* admin-session-started
* admin-session-ended
* dev-key-generated
* settings-changed
* log-settings-changed

</details>

The last eleven are the platform's own security trail — who logged in and out, whose login
failed and why, impersonation, admin sessions, developer keys, settings changes, and server
starts. They are listed on the **System Activity** tab of [Usage Analysis](usage-analysis.md).

## Accessing audit logs

You can access audit logs in a number of ways:
* From the [Context Panel](../../datagrok/navigation/panels/panels.md#context-panel)
  * **Activity**: This [info pane](../../datagrok/navigation/panels/info-panels.md) shows server-side activity for the current entity over the past seven days
  * **Usage**: This info pane shows a high-level overview of server-side usage for the current entity in the last seven days
* From the [Console](../../datagrok/navigation/panels/panels.md#console)
* From a Datagrok app, [Usage Analysis](usage-analysis.md)
* In your own log backend, such as Amazon CloudWatch, Google Cloud Logging, or any
  OpenTelemetry collector (see [Export logs](#export-logs))

### Export logs

**Log sync** sends log events to external destinations as they are posted,
without waiting for them to be saved to the database. To add a destination:

1. On the **Sidebar**, click **Settings** (<FAIcon icon="fa-solid fa-gear"/>) > **Logger**, and expand **Log sync**.
1. Click **Add new sync block**.
1. Set **Cloud** to the destination and fill in its settings (see the table below).
1. Select the **Levels** to sync. Only levels that are posted can be synced. See [Changing settings for default logs](#changing-settings-for-default-logs).
1. Click **Save and apply**.

Add one block per destination. For example, send errors to one CloudWatch log
group and audit events to another.

<details>
<summary>Log sync settings</summary>

| Setting        | Applies to            | Description                                                                                          |
|----------------|-----------------------|------------------------------------------------------------------------------------------------------|
| **Enabled**    | All                   | Turns the destination off without deleting it                                                        |
| **Levels**     | All                   | Levels to sync: error, warning, info, debug, audit, usage                                            |
| **Event type** | All                   | Narrows the destination to a single audit or usage event type. Leave empty for all                   |
| **Params**     | All                   | Syncs only events whose parameters have the given values                                              |
| **Batch Size** | All                   | Maximum number of events sent in one batch                                                            |
| **Format**     | All                   | `json` (one object per event, queryable in the backend) or `text`                                     |
| **Cloud**      | All                   | **Amazon CloudWatch**, **Google Cloud Logging**, or **OpenTelemetry (OTLP)**                          |
| **Connection** | All                   | An AWS or GCP connection with write access. Leave empty to authenticate with the instance role       |
| **Log Group**, **Stream** | Amazon CloudWatch | Destination log group and stream. Without a connection, Datagrok creates them if they don't exist. With a connection, create them first |
| **Log Name**   | Google Cloud Logging  | Destination log                                                                                       |
| **Endpoint**   | OpenTelemetry         | OTLP/HTTP logs endpoint, for example `http://otel-collector:4318/v1/logs`                             |
| **Auth**       | OpenTelemetry         | `none`, `bearer` (with **Token**), `aws-sigv4`, or `gcp-oauth`                                        |
| **Headers**, **Compress**, **Max Payload** | OpenTelemetry | Extra HTTP headers, gzip compression, and the size at which a batch is split          |
| **Alerts**     | OpenTelemetry         | Sends alert and heartbeat records regardless of the selected levels. On by default                    |

</details>

Besides regular events, the server emits records that let your monitoring
system page you:

* `alert-opened`, `alert-escalated`, and `alert-resolved` when a problem starts, gets worse, or clears. Alerts
  cover a failed service health check, requests slowing down or failing across the board, one error hitting many
  users, an account that keeps failing to sign in, and a problem reported by a user.
* `alert-firing` every five minutes for each alert that is still open.
* `heartbeat` every five minutes. A missing heartbeat means the instance stopped reporting.

`alert-firing` and `heartbeat` aren't stored in the database. They exist only in Log sync destinations.

To learn more, see [Export logs to Amazon CloudWatch](../../datagrok/solutions/teams/it/log-export-cloud-watch.md).

:::note

On a platform-managed instance, the **Log sync** section is hidden. The deployment owns the destinations.

:::

### Retention

By default, Datagrok keeps all events. To delete old events automatically,
turn on the **Garbage Collector** in **Settings** > **Admin**. It removes audit
events after 365 days, errors, warnings, info, and usage events after 183 days,
and debug events after 30 days. You can change each period. To keep events
longer, [export them](#export-logs) to a backend with its own retention policy.

## Logging events

Datagrok logs both client-side and server-side activities. 

### Client-based actions

Thanks to Datagrok's [in-memory database](../../develop/under-the-hood/performance.md#in-memory-database), many
user activities related to exploratory data analysis, such as opening local
files, aggregating tables, or adding viewers, occur entirely on the client.
These actions generate internal named events that serve various purposes:

* Integration of different Datagrok components
* Ad-hoc customizations and extensions by plugins
* Audit logs

Datagrok logs a reasonable default set of client-based actions. For example, it
logs the opening of a file and its name, but not the content. Less important
actions, like "rows selected", aren't logged. To learn how to customize what
gets logs, see [Customize audit logging](#customize-audit-logging).

### Server-related actions

When an action modifies the server's state in any way, it's recorded in the
audit log in addition to the log file. You can see a seven-day history of such
activity on the **Context Panel** under **Activity**, and the
high-level overview of the entity's usage under **Usage**.

By default, each explicitly executed
[function](../../datagrok/concepts/functions/functions.md) is logged,
along with its parameter values.

## Customize audit logging

To customize which events are logged and how they are handled, you have these options:

1. Configure logging settings: For any user or group, add or remove logged
   events from the default set
1. Create custom logging events and make other customizations using [JS API](#customizing-logging-via-javascript-api)

### Changing settings for default logs

To access logging settings, on the **Sidebar**,  click **Settings (<FAIcon
icon="fa-solid fa-gear"/>) > Logger**. The **Settings - Logger** view opens.
From here, you can specify which audit events are logged for a user or group,
and set timeouts for logging events. 

To specify which audit events should be logged for a user or group:

1. Under **Log settings** add a user or group and the desired logging level (Standard, Verbose, Custom).
1. Optional. Adjust settings for the selected option by clicking the **Gear icon** next to the dropdown.

### Customizing logging via JavaScript API

To customize audit logs using Datagrok's JavaScript API:

1. Create a `Logger` object:

     ```javascript
     // Create logger
     let logger = new DG.Logger();
     // You can pass a callback that sets predefined parameters
     let logger = new DG.Logger((m) => m.params['persistent']= 'value');
     ```

1. To log events, use the `logger.log()` method:

     ```javascript
     // Default type is "log"
     logger.log('HELLO WORLD', {test: 'value', 'foo': 'bar'});
     // But you can specify another
     logger.log('HELLO WORLD', {test: 'value', 'foo': 'bar'}, 'my_log');
     ```

1. Errors are logged automatically, but you can log them explicitly and specify
   the stack trace for grouping errors:

     ```javascript
     // Specify stack trace to group errors
     logger.error('Error!', stackTrace);
     ```

See [this sample](https://public.datagrok.ai/js/samples/ui/ui-events)
for a demonstration of how to handle events.