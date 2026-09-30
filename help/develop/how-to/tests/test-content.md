---
title: "Test dashboards and queries"
description: Write a smoke suite that opens your canonical dashboards and runs your canonical queries, and run it in CI and before every upgrade.
keywords:
  - test dashboards
  - test queries
  - smoke tests
  - quality engineering
  - upgrade verification
  - grok test
---

[Package tests](add-package-tests.md) are usually written by developers to test their own code. The
same framework also works for the content your organization depends on: the dashboards and queries
that people open every day. This page shows how to build a _smoke suite_ for that content, a short
set of tests that answers "does everything important still open and return sensible data?"

Run the suite:

* **Nightly** against a test server, to catch problems from changes in your data sources
* **Before every upgrade** reaches production, because package database migrations don't roll back
* **After a change** to a shared query or connection, since dashboards depend on them

:::note Why a suite is needed

Datagrok doesn't check a dashboard's dependencies for you. If a query a dashboard depends on is
deleted, the dashboard fails to open, and a renamed column can break its viewers (see
[what happens when the source query changes](../../../datagrok/concepts/project/dashboard.md)).
A smoke suite is how you find out before your users do.

:::

## Create a package for the suite

Keep the suite in its own package, owned by the team responsible for content quality. Create a
package, add a `src/tests` folder, and import your test files in `src/package-test.ts`, as described
in [Add package tests](add-package-tests.md).

## Test that a dashboard opens

A test can open a saved dashboard by name and check what it loaded: the tables, their columns, and
whether the row counts are plausible.

```js
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';

category('Canonical dashboards', () => {
  test('Assay QC opens with data', async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.open('Assay QC');

    const t = grok.shell.t;
    expect(t != null, true, 'No table loaded');
    expect(t.rowCount > 0, true, 'Table is empty');
    for (const name of ['compound_id', 'assay', 'ic50'])
      expect(t.columns.contains(name), true, `Missing column: ${name}`);
  });
});
```

## Test that a query returns what you expect

For each query your dashboards rely on, run it with fixed parameters and compare the result with
what you know to be true. Choose inputs whose results are stable, such as a reference compound or a
closed study.

```js
category('Canonical queries', () => {
  test('Assay results for the reference compound', async () => {
    const df = await grok.functions.call('Warehouse:AssayResults', {compoundId: 'REF-0001'});
    expect(df.rowCount >= 1, true, 'No rows for the reference compound');
    expect(df.columns.contains('ic50'), true, 'Missing column: ic50');
  });
});
```

For results that must match exactly, compare the whole table with a stored reference using
`expectTable`.

## Run the suite

From the package folder, run the suite against any server configured in your `grok` settings:

```bash
grok test --host test --csv
```

`grok test` builds the package, publishes it to that server in debug mode, and runs the tests in a
headless browser. `--csv` saves a report your CI can read, and `--category` runs a single category.
To run the same tests interactively, open **Test Manager** (see
[Test packages](test-packages.md#test-manager)).

Run the suite from your CI system on a schedule and before every upgrade. Datagrok doesn't run tests
on a schedule by itself.

## See also

* [Add package tests](add-package-tests.md)
* [Test packages](test-packages.md)
* [Monitor data connections](../../../access/databases/monitor-connections.md)
* [Dashboards](../../../datagrok/concepts/project/dashboard.md)
* [Running Datagrok in the enterprise](../../../govern/enterprise-guide.md)
