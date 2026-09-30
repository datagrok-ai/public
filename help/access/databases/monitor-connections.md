---
title: "Monitor data connections"
description: Test a data connection on demand, and let Datagrok check your critical connections continuously and alert when one becomes unreachable.
keywords:
  - monitor data connections
  - connection health
  - test connection
  - connection alert
  - grok s connections test
---

Dashboards and apps are only as healthy as the data sources behind them. Datagrok holds the
credentials for your [connections](databases.md), so it's also the place to find out when one stops
working.

## Test a connection on demand

* **In the UI.** Open the connection's dialog and click **TEST**.
* **From the command line.** `grok s connections test <id-or-name>`, using the
  [`grok s` CLI](../../develop/server-management.md). It prints `ok` when the connection works, and
  the error otherwise.

A test confirms that Datagrok can reach the source and sign in. It doesn't run a query, so it won't
notice a missing table or a revoked grant. To cover those, add
[tests to the queries](../../develop/how-to/tests/test-content.md) that matter.

## Monitor connections continuously

:::note Available in the next release

Continuous connection monitoring is available on the current bleeding-edge build and ships in the
next release.

:::

Tag a connection `monitor`, and Datagrok tests it every five minutes, the same way the **TEST**
button does. If two tests in a row fail, Datagrok raises an alert that names the connection and the
kind of failure: authentication, network, timeout, or driver. The alert closes by itself on the
first successful test, so a single network blip never raises one.

Tag every connection that a production dashboard depends on.

## See also

* [Databases](databases.md)
* [Test queries and dashboards](../../develop/how-to/tests/test-content.md)
* [Running Datagrok in the enterprise](../../govern/enterprise-guide.md)
