---
title: "Test queries and dashboards"
description: Add test cases to the queries your dashboards depend on, so a broken query or data source is caught before users open the dashboard.
keywords:
  - test queries
  - test dashboards
  - query test annotation
  - quality engineering
  - upgrade verification
---

A dashboard is only as reliable as the queries behind it. If a query fails, or a column it returns is
renamed, the dashboards built on it break (see
[what happens when the source query changes](../../../datagrok/concepts/project/dashboard.md)).
Testing the queries catches most of those problems before users see them.

## Add a test to a query

Add a `--test:` line to the query's header. It calls the query with fixed parameters:

```sql
--name: AssayResults
--connection: Warehouse:Assays
--input: string compoundId
--test: AssayResults('REF-0001')
select * from assay_results where compound_id = @compoundId
```

The test fails if the query fails. Choose parameters whose results are stable, such as a reference
compound or a closed study. A query can have several `--test:` lines. The same annotation works for
scripts and other functions, and a test that evaluates to `false` also fails (see
[Testing functions](add-package-tests.md#testing-functions) for the full syntax).

## Run the tests

Tests on the queries in a package run with the package's tests: interactively in
[Test Manager](test-packages.md#test-manager), or from the command line with `grok test`, for
example before an upgrade reaches production.

:::note Coming next

We're extending the built-in alerts so query tests also run on a schedule and alert the query's
owner when one fails, the same way [connection monitoring](../../../access/databases/monitor-connections.md)
does.

:::

## See also

* [Add package tests](add-package-tests.md)
* [Test packages](test-packages.md)
* [Monitor data connections](../../../access/databases/monitor-connections.md)
* [Dashboards](../../../datagrok/concepts/project/dashboard.md)
* [Running Datagrok in the enterprise](../../../govern/enterprise-guide.md)
