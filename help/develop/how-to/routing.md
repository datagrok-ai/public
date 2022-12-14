<!-- TITLE: Routing -->

# Routing

Sometimes it is useful to share a link with a colleague or business executive that opens an app in a specific state
depending on the URI and points to the important information that the one would like to share. Datagrok provides routing
support for package developers to achieve this.

## Routing tutorial

Routing is implemented by manipulating the `basePath` and `path` properties of the `View`. Let's build a simple app that
opens tables, based on the URL path and applies selection according to the parameters passed in the URL.

We start off by getting the path and splitting it in segments:

```javascript
//name: Test
//tags: app
export function test() {
  const pathSegments = window.location.pathname.split('/');
  // etc.
}
```

Depending on the `pathSegments` length, we choose to either start the application with the default behavior or handle
the URI parameters. Typically, when the app starts from the Apps page, the `pathSegments` would contain 3 elements:
`['', 'apps', 'Test']` because the `pathname` is `/apps/Test` by default for the app called `Test`.

Let's say the default behaviour would be to add a couple of tables and set default paths for them.

```javascript
if (pathSegments.length <= 3) {
  //Adding demog table view
  const demog = grok.data.testData('demog');
  const demogView = grok.shell.addTableView(demog);
  demogView.scatterPlot();
  demogView.basePath = '/demog';
  demogView.path = '/All';
  grok.shell.v = demogView;

  //Adding random walk table view
  const wells = grok.data.testData('wells');
  const wellsView = grok.shell.addTableView(wells);
  wellsView.scatterPlot();
  wellsView.basePath = '/wells';
  wellsView.path = '/All';
  grok.shell.v = wellsView;
}
```

Note that by setting `.basePath` we add a segment to the `pathname` and by setting `.path` we add another segment after
`.basePath`. In the case of demog table, the `pathname` would be `/apps/Test/demog/All`. Link with such URI could be
opened by another user, but for this, to work we also need to handle the case, when the `pathSegments` contains more
than 3 segments.

```javascript
const tableName = pathSegments[3];
const selectionLabel = pathSegments[4] ?? 'All';
const table = grok.data.testData(tableName as testData);
const tableView = grok.shell.addTableView(table);
tableView.basePath = `/${tableName}`;
setSelection(tableView, selectionLabel);
grok.shell.v = tableView;
```

Path segments can be used to provide intended behavior. In this case, the table name is retrieved from the fourth path
segment, which corresponds to the `.basePath` set previously, and the selection options can be retrieved from the fifth
path segment.

That's it! Now you've learned how to use routing to enhance your apps. For additional info see the full code and useful
links below.

## Tutorial code

Here's the full code used in the tutorial.

```javascript
/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

type testData = "wells" | "demog" | "biosensor" | "random walk" | "geo" | "molecules" | "dose-response";

//name: Test
//tags: app
export function test() {
  const pathSegments = window.location.pathname.split('/');

  if (pathSegments.length <= 3) {  //Fresh app start
    //Adding demog table view
    const demog = grok.data.testData('demog');
    const demogView = grok.shell.addTableView(demog);
    demogView.scatterPlot();
    demogView.basePath = '/demog';
    demogView.path = '/All';
    grok.shell.v = demogView;

    //Adding random walk table view
    const wells = grok.data.testData('wells');
    const wellsView = grok.shell.addTableView(wells);
    wellsView.scatterPlot();
    wellsView.basePath = '/wells';
    wellsView.path = '/All';
    grok.shell.v = wellsView;
  } else {  //Handle routing
    const tableName = pathSegments[3];
    const selectionLabel = pathSegments[4] ?? 'All';
    const table = grok.data.testData(tableName as testData);
    const tableView = grok.shell.addTableView(table);
    tableView.basePath = `/${tableName}`;
    setSelection(tableView, selectionLabel);
    grok.shell.v = tableView;
  }
}

function setSelection(tableView: DG.TableView, label: string) {
  tableView.path = `/${label}`;
  if (label === 'All') {
    tableView.dataFrame.selection.init(_ => true);
    return;
  }
  const [colName, category] = label.split('=');
  if (!colName || !category)
    throw new Error(`PathError: wrong path format '${label}'. Should be 'colName=value'.`)
  tableView.dataFrame.rows.match(`${colName} = ${category}`).select();
}
```

See also:

* [How-to: Build an application](./build-an-app.md)
* [Routing](../../datagrok/routing.md)
