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

When a query fails, or a column it returns is renamed, the dashboards that use it break (see
[what happens when the source query changes](../../../datagrok/concepts/project/dashboard.md)).
Adding tests to your important queries finds these problems before users do.

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

Tests on queries in a package run with the package's other tests: in
[Test Manager](test-packages.md#test-manager), or from the command line with `grok test`. Run them
on your test server before an upgrade reaches production.

Datagrok doesn't run query tests on a schedule.

## See also

* [Add package tests](add-package-tests.md)
* [Test packages](test-packages.md)
* [Monitor data connections](../../../access/databases/monitor-connections.md)
* [Dashboards](../../../datagrok/concepts/project/dashboard.md)
* [Running Datagrok in the enterprise](../../../govern/enterprise-guide.md)
