---
title: "Routing"
---

Sometimes it is useful to share a link with a colleague or business executive that opens an app in a specific state
depending on the URI and points to the important information that the one would like to share. Datagrok provides routing
support for package developers to achieve this.

## Routing tutorial

Routing is implemented by manipulating the `path` property of the `View`. Let's build a simple app that
opens tables, based on the URL path and applies selection according to the parameters passed in the URL.

We start off by getting the path and splitting it in segments:

```javascript
//name: Test
//input: string path {meta.url: true}
//input: string searchParam
//tags: app
export function test(url, searchParam) {
  const pathSegments = url.split('/').filter((s) => s != '');
  // etc.
}
```

Application gets the path in input parameter marked with `meta.url: true`, and all URL parameters mapped by names.

Depending on the `pathSegments` length, we choose to either start the application with the default behavior or handle
the URI parameters. Typically, when the app starts from the Apps page, the `pathSegments` would contain no elements because path equals `/`

Let's say the default behavior would be to add a couple of tables and set default paths for them.

```javascript
if (pathSegments.length == 0) {
  //Adding demog table view
  const demog = grok.data.testData('demog');
  const demogView = grok.shell.addTableView(demog);
  demogView.scatterPlot();
  demogView.path = '/Demog/All';
  grok.shell.v = demogView;

  //Adding random walk table view
  const wells = grok.data.testData('wells');
  const wellsView = grok.shell.addTableView(wells);
  wellsView.scatterPlot();
  wellsView.basePath = '/wells/All';
  grok.shell.v = wellsView;
}
```

Note that by setting `.path` we add path postfix after base path. In the case of demog table, the `pathname` would be `/apps/TestPackage/demog/All`. Link with such URI could be
opened by another user, but for this, to work we also need to handle the case, when the `pathSegments` contains segments `demog` and `All`

```javascript
const tableName = pathSegments[0];
const selectionLabel = pathSegments[1] ?? 'All';
const table = grok.data.testData(tableName as testData);
const tableView = grok.shell.addTableView(table);
tableView.path = `/${tableName}/${selectionLabel}`;
setSelection(tableView, selectionLabel);
grok.shell.v = tableView;
```

That's it! Now you've learned how to use routing to enhance your apps. For additional info see the full code and useful
links below.

## Function base path

By default, function base path looks like `/apps` of `/browse/apps` prefix with package name and then app name.
If there is an only app in the package, app name can be omitted;

For example:
```
/apps/TestPackage
/browse/apps/TestPackage
```
or in case of 2 and more apps in the package:
```
/apps/TestPackage/Test
/browse/apps/TestPackage/Test
```

You can customize package segment by adding meta section in package.json file:

```json
  "meta": {
    "url": "/some/custom/route"
  }
```

Now URL starts to look like this:

```
/apps/some/custom/route
/browse/apps/some/custom/route
```
or in case of 2 and more apps in the package:
```
/apps/some/custom/route/Test
/browse/apps/some/custom/route/Test
```
Package name based URL will also work, but will be immediately replaced by the new one, when user opens it.

Also, you can specify app URL alias with `meta.url` tag:

```javascript
//name: Test
//input: string path {meta.url: true}
//input: string searchParam
//meta.url: /application
//tags: app
export function test(url, searchParam) {
  const pathSegments = url.split('/').filter((s) => s != '');
  // etc.
}
```

In this case, URL becomes:
```
/apps/some/custom/route/application
/browse/apps/some/custom/route/application
```

If app URL contains the only slash, application becomes default for the package:
```
/apps/some/custom/route
/browse/apps/some/custom/route
```


## Tutorial code

Here's the full code used in the tutorial.

```javascript
/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

type testData = "wells" | "demog" | "biosensor" | "random walk" | "geo" | "molecules" | "dose-response";

//name: Test
//input: string path {meta.url: true}
//input: string filter
//tags: app
export function test(path: string, filter: string) {
  const pathSegments = url.split('/').filter((s) => s != '');

  if (pathSegments.length == 0) {  //Fresh app start
    //Adding demog table view
    const demog = grok.data.testData('demog');
    const demogView = grok.shell.addTableView(demog);
    demogView.scatterPlot();
    demogView.path = '/demog/All';
    grok.shell.v = demogView;

    //Adding random walk table view
    const wells = grok.data.testData('wells');
    const wellsView = grok.shell.addTableView(wells);
    wellsView.scatterPlot();
    wellsView.path = '/wells/All';
    grok.shell.v = wellsView;
  } else {  //Handle routing
    const tableName = pathSegments[0];
    const selectionLabel = pathSegments[1] ?? 'All';
    const table = grok.data.testData(tableName as testData);
    const tableView = grok.shell.addTableView(table);
    setSelection(tableView, tableName, selectionLabel, filter);
    grok.shell.v = tableView;
  }
}

function setSelection(tableView: DG.TableView, name: string, label: string, filter: string) {
  tableView.path = `/${name}/${label}`;
  if (label === 'All') {
    tableView.dataFrame.selection.init(_ => true);
    return;
  }
  const [colName, category] = filter.split('=');
  if (!colName || !category)
    throw new Error(`PathError: wrong filter format '${label}'. Should be 'colName=value'.`)
  tableView.dataFrame.rows.match(`${colName} = ${category}`).select();
}
```

See also:

* [How-to: Build an application](./build-an-app.md)
* [Routing](../../datagrok/navigation/routing.md)
