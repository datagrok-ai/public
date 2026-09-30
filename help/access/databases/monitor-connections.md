---
title: "Monitor data connections"
description: Test a data connection on demand, or tag it so Datagrok tests it every five minutes and raises an alert when it fails.
keywords:
  - monitor data connections
  - connection health
  - test connection
  - connection alert
  - grok s connections test
---

This page shows how to check that a [data connection](databases.md) works, on demand or
continuously.

## Test a connection on demand

* **In the UI.** Open the connection's dialog and click **TEST**.
* **From the command line.** `grok s connections test <id-or-name>`, using the
  [`grok s` CLI](../../develop/server-management.md). It prints `ok` when the connection works, and
  the error otherwise.

A test checks that Datagrok can reach the source and sign in. It doesn't run a query, so it doesn't
detect a missing table or a revoked grant. For those, add
[tests to your queries](../../develop/how-to/tests/test-content.md).

## Monitor connections continuously

:::note Bleeding-edge build

Connection monitoring is available on the bleeding-edge build and not yet in a stable release.

:::

Add the `monitor` tag to a connection, and Datagrok tests it every five minutes, the same way the
**TEST** button does. After two failed tests in a row, Datagrok raises an alert with the connection
name and the type of failure: authentication, network, timeout, or driver. The alert closes after
the next successful test. Requiring two failures avoids alerts for a single dropped request.

Tag every connection that production dashboards use.

## See also

* [Databases](databases.md)
* [Test queries and dashboards](../../develop/how-to/tests/test-content.md)
* [Running Datagrok in the enterprise](../../govern/enterprise-guide.md)
