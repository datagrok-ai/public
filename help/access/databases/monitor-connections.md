---
title: "Monitor data connections"
description: Test data connections on demand and on a schedule, alert through your own tooling when a source fails, and use audit and usage data to spot slow sources.
keywords:
  - monitor data connections
  - connection health
  - test connection
  - canary query
  - scheduled query
  - grok s connections test
---

Dashboards and apps are only as healthy as the data sources behind them. This page shows how to
check a [connection](databases.md) on demand, how to check it automatically and alert when it
fails, and where to look for trends such as a source that is getting slower.

## Test a connection on demand

There are three ways to test a connection:

* **In the UI.** Open the connection's dialog and click **TEST**.
* **From the command line.** `grok s connections test <id-or-name>`, using the
  [`grok s` CLI](../../develop/server-management.md).
* **From code.** In a package or a script, `connection.test()` returns `"ok"` on success or the
  error message on failure. It doesn't throw.

:::note What a test checks

A connection test confirms that Datagrok can reach the source and sign in: it opens a connection
and closes it. It doesn't run a query, so it won't notice a missing table, a revoked grant, or a
schema change. To cover those, also run a small query.

:::

## Check connections automatically

A _canary_ is a small, scheduled check that fails loudly when a source stops working. We recommend
one canary for each connection your dashboards depend on, in two parts.

### Test from your own scheduler

Run the connection test from the scheduler your team already uses, such as cron, your CI system,
or your monitoring agent, so a failure reaches the people and channels that already handle
incidents. `grok s connections test` prints `ok` and exits with code 0 when the connection works.
When it fails, it prints the error as JSON and exits with code 1.

```bash
#!/usr/bin/env bash
# Test each critical connection; exit non-zero if any fails.
status=0
for conn in "Warehouse:Assays" "Warehouse:Compounds"; do
  if ! out=$(grok s connections test "$conn" --host prod 2>&1); then
    echo "Connection $conn failed:"
    echo "$out" | head -5
    status=1
  fi
done
exit $status
```

Use a service account for the `grok s` API key, with access only to the connections it tests.
If more than one connection shares a name, pass the connection's ID instead.

To check the data as well as the connection, have the same script run a small saved query.
`grok s functions run` prints the query's result as CSV and exits with code 1 if the query fails:

```bash
grok s functions run 'Warehouse:AssaysCanary()' --host prod > /dev/null || status=1
```

### Run a small query inside Datagrok

To check that the data itself is reachable, schedule a small query on the connection, for example
one that counts rows in a table your dashboards use. Scheduled queries run on the server with the
`#schedule` annotation (in SQL, `--schedule`), and `schedule.runAs` sets the group or role whose
permissions they use (see [Scheduling](../../datagrok/concepts/functions/functions.md#scheduling)).

```sql
--name: AssaysCanary
--connection: Warehouse:Assays
--schedule: */15 * * * *
--schedule.runAs: Monitoring
select count(*) from assay_results
```

Each scheduled run is recorded, and a failed run is recorded with its error. Datagrok doesn't
notify anyone when a scheduled run fails, so pair the query with an alert of your own, such as the
scheduler check above.

:::caution Keep canaries out of the cache

Query results are cached only if a query opts in with
[`meta.cache`](../../develop/how-to/functions/cache-function-results.md). A cached query keeps
answering after its source fails, until the cache expires. Don't enable caching on a canary query,
or it will report healthy while the database is down.

:::

## Watch for trends

Failures are only half the picture. A source that is getting slower affects every dashboard built
on it before anything fails outright.

* **Audit trail.** Every query run and every error is an [audit record](../../govern/audit/audit.md),
  as is every change to a connection.
* **Usage Analysis.** The [Usage Analysis](../../govern/audit/usage-analysis.md) app shows execution
  times for queries and functions, so you can see a source slowing down over weeks rather than
  hearing about it from users.
* **Log export.** Send errors to CloudWatch or Google Cloud Logging by record type, so connection
  failures appear where your team already looks (see
  [log export](../../govern/access-control/data-connection-credentials.md#for-logs-export-to-cloudwatch)).

## See also

* [Databases](databases.md)
* [Test dashboards and queries](../../develop/how-to/tests/test-content.md)
* [Running Datagrok in the enterprise](../../govern/enterprise-guide.md)
* [Server management with grok s](../../develop/server-management.md)
