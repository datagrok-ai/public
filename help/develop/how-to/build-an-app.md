<!-- TITLE: Building an application -->
<!-- SUBTITLE: -->

<!-- This is a developer's view on the Datagrok applications -->

# Datagrok Applications

Applications are built on top of the Datagrok platform. The Datagrok application serves a targeted
need around a particular problem or an area. For example, this may be a bioreactor modeling application,
an application for working with and understanding molecular sequence data, a molecular database browser,
a Covid-19 or weather info panel, etc.

Datagrok applications are developed in JavaScript / TypeScript using our rich [Datagrok JavaScript API](),
with also using a range of the [scripting]() languages we support, such as R or Python, by calling scripts,
and with connecting to third-party data sources and services.

Technically, applications are [functions](../overview/functions/function.md) tagged with the `#app` tag.
In this function, you take control over what the platform will let user do and see next after the app is
executed. Think of it as of the application's entry point, such as a `Main` function in C# console app.

## Entry point

The Datagrok platform is highly extensible. New functionality is delivered to a Datagrok instance as packages.
A Datagrok [package]() might contain zero, one, or more Datagrok applications. These come along with other
entities in the package, which these applications may be using, such as connections, viewers, scripts, etc.

Consider a simple example of a webpack-based package with just two trivial apps in one `src/package.js`:

```
import * as grok from 'datagrok-api/grok';

export let _package = new DG.Package();

//name: Test 1
//tags: app
export function test1() {
  grok.shell.info('Test 1');
}

//name: Test 2
//tags: app
export function test2() {
  grok.shell.info('Test 2');
}
```

To make this run on Datagrok, follow these `grok create` [steps](../develop/develop.md#getting-started)
to prepare this simple package and deploy it.

1. [Create a new package](../develop.md#getting-started)
2. Add an app to it via `grok add app <APP_NAME>` and see the minimal default structure for it, or alternatively just copy-paste the JS snippet from the above to `src/package.js`

After deploying this package, you'd find these 2 apps via `Functions | Apps` in the activity bar
situated on the left side of Datagrok's main window. Run both of the apps and notice the two
different tooltips popping up. You may also call these same entry points by an URL:
`https://public.datagrok.ai/apps/<PACKAGE_NAME>/TestUI1`, and a similar one for `TestUI2`.
Note that in case there is only one application `<APP>` defined in the package, the
corresponding URL will simply be `https://public.datagrok.ai/apps/<APP>`, omitting the
`<PACKAGE_NAME>` part.

This simple example finishes explaining the purpose of the entry points. Yet, trivial popups aren't
something one typically builds as an application. Let's look at a more UI-rich side of things.

## Main view

Most applications built on Datagrok start with a Datagrok's view. A view is a set of visualizations and
controls grouped together. Typically, the view is associated with a particular dataframe, in this case it's
called a table view. However, essentially a view can contain pretty much anything.

Imagine you are composing an application. You'd most likely start off with the root / main view,
add logical blocks to it either through simple [div-s](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/ui/sidebar.js),
or through [`splitH`/`splitV`](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/ui/layouts/splitters.js),
populate these blocks with visualizations and controls, maybe add a sidebar, add event handlers and so forth.
Our internal application [Usage Analysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis)
demonstrates such approach.

Another approach is found in a [Discovery](https://github.com/datagrok-ai/public/tree/master/packages/Discovery) application.
There we reuse a particular kind of view: a table view, and centralize the rest of the UI around it. In
this application you'd also find some useful techniques one needs sooner or later in many applications,
such as modifying the app's URI or hiding the side panels of the Datagrok's main UI.

Read more on creating custom views [here](./custom-views.md).

## Structuring code

Prehaps, one of the main things to know about webpack is that it allows you to structure JavaScript code
in a way similar to how you are used to it in enterprise-grade environments, such as Java or .NET.

Let's expand our previous example leveraging an `include` capability of webpack.

Create two separate files —

`src/test-app-01.js`:

```
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

export function makeTest01() {
  let t = grok.data.demo.demog();
  let view = grok.shell.addTableView(t);
  let acc = view.toolboxPage.accordion;
  acc.addPane('Demo 1', () => ui.divText(
    `Cells count: ${t.rowCount * t.columns.length}`), true, acc.panes[0]);
}
```

`src/test-app-02.js`:

```
import * as grok from 'datagrok-api/grok';

export function makeTest02() {
  grok.shell.info('Demo 2');
}
```

And modify the `src/package.js`:

```
import {TestUI1} from './test-app-01.js';
import {TestUI2} from './test-app-02.js';

export let _package = new DG.Package();

//name: Test 1
//tags: app
export function test01() {
    return makeTest01();
}

//name: Test 2
//tags: app
export function test02() {
    return makeTest02();
}
```

That would be an overkill to structure our super-trivial apps this way, though it is a highly desired
practice for anything real-life built for production use. Notice how we include parts of Datagrok API
for the corresponding features in `src/test-app-01.js`. You may find this `datagrok-api` is a
pre-defined location, as per `webpack.config.js`.

## Providing for IntelliSense

We recommend to restore the package dependencies before starting development with an IDE.
After the package is created, simply invoke `npm install` inside the package folder. This
will bring NPM modules defining the Datagrok API

An alternative way to bring IntelliSense capability for Datagrok classes is by cloning
the entire [public repo](https://github.com/datagrok-ai/public) and open this whole folder
in the IDE. This is suitable for the case when you develop a package which is
a [part](https://github.com/datagrok-ai/public/tree/master/packages) of our public repo.
The folder `js-api` in the root of the public repo is a source of our JS API, that's how
IDE discovers needed IntelliSense information.

# Application Development

Having all the above in mind, there is pretty much anything you can further do inside the application,
leveraging a full scale of platform capabilities. However, there are certain aspects of utmost interest
in almost any application:

* accessing and persisting data
* working with [dataframes]()
* performing computations
* interactive visualizations
* managing privileges
* application lifecycle

The following chapter guides through these key development topics. Take it as a birds-eye overview of the
applcation development area, grasp the major building blocks, and proceed to the articles and samples
referenced in the guide for the further details.

## Code samples

We provide a diversed set of code snippets of the API use, and sample packages with viewers, applications, and so forth.

* For short samples of using API, go to https://public.datagrok.ai/js and observe a "Samples" block, or alternatively —
access it via the "Help" button at the bottom of the activity bar on the left of the Datagrok's main window,
and then follow to `JavaScript API Samples`.

* The sources of these snippets are all located at https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples.

* Check all our demo codes at `https://github.com/datagrok-ai/public/tree/master/packages`. For instance, you may
find there a package [Viewers](https://github.com/datagrok-ai/public/tree/master/packages/Viewers) which showcases
a variery of viewers one can build on Datagrok, or a [Chem package](https://github.com/datagrok-ai/public/tree/master/packages/Chem),
which shows how in-browser computations with WebAssembly may be organized on the platform.

Developers who's been starting with the platform recently has shared their experience: the samples above became
a major source of knowledge for their daily work with the platform. We recommend you to recall these locations too,
and resort to them for resolving your daily technical questions in addition to posting questions and suggestions
on our [Community Forum](https://community.datagrok.ai/) at https://community.datagrok.ai/.

## Data Access

There's a variety of data sources which Datagrok can handle out of the box:

* Web Services (REST endpoints) specified via Swagger / OpenAPI ([link](../access/open-api.md))
* Arbitrary REST APIs (including these outside host's domain) ([link](./access-data.md#rest-endpoints))
* Access a file from the Datagrok server space (either shared or Home) ([link](./access-data.md#file-shares))
* Reading from databases: Datagrok supports more than 30 ([link](./access-data.md#parameters))
* Access a file from a network file share (Windows Shares, SAMBA and the like) ([link](./access-data.md#file-shares))

Connections to databases, respective database queries, and file shares connections live in Datagrok as entities.
These may be published to only you, the package author, or to any specific user group.
Main code work for accessing databases is concentrated around `grok.data.query` method ([a sample]()), whereas
accessing files from file shares and servers is made by `grok.functions.eval('OpenServerFile("...")')` ([a sample]()).

When it comes to connections, we strongly encourage you not to store their credentials directly inside your
application packages. The chapter ["Managing privileges"](#managing-privileges) discusses credentials management in detail.

## Dataframes

### Creating and working with dataframes

Dataframes are at the heart of Datagrok. DataFrame is a high-performance, easy to use tabular structure with
strongly-typed columns of different types. Supported types are: string, int, float, bool, DateTime, bigint and
[qualified numbers]().

Datagrok dataframes are highly optimized. They are implemented with the proprietary, unique technology allowing
to efficiently work with huge datasets in the browser. Essentially, it is a columnar in-memory database that was
engineered from scratch and optimized for the purpose of exploratory data analysis, interactive visualizations,
and machine learning.

Note that Datagrok dataframes operate entirely inside the browser, but not on our [compute server]().
However, it's possible to pass dataframes to scripts (in Python, R and others) which run on the server,
and [get dataframes in return](#computations). 

You get dataframes within your application in various ways. Dataframe may be a table rendered by a table view,
a new dataframe constructed from a set columns, a dataframe constructed from a file in a file share, a CSV file
uploaded to a browser, a dataframe return by a script, and so forth. You can add and delete rows and columns of
the dataframe, view it in a table view and let viewer be attached to it to render. Dataframes can also be calculated
on the flight for aggregations.

Dataframes are comprised of [columns](). Columns may be used as Datagrok functions arguments. For columns, it's
possible to get its underlying dataframe. In return, columns are comprised of cells, and it's possible to get cell's
underlying column. There is also a diversed [system of events](https://datagrok.ai/js-api/DataFrame) one can
subsrcribe on a dataframe.

Let's create a dataframe and check what we can do with it.

Try the below snippets in our interactive [JS playground](https://public.datagrok.ai/js), but
don't forget to request `import * as DG from 'datagrok-api/DG'` in case you're using this code
from a webpack package.

```
let df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.STRING, 'Person', ['Alice', 'Bob', 'Susan']),
  DG.Column.fromList(DG.TYPE.FLOAT, 'Height', [178, 190, 165])
]);
df.columns.addNew('City', 'string');
df.columns.byName('City').set(0, 'New York');
df.set('City', 1, 'Chicago');
df.set('City', 2, 'San Francisco');
df.rows.addNew(['John', 175, 'Dallas']);
df.rows.removeAt(0);
grok.shell.addTableView(df);
```

There we create a new dataframe from columns, add a column and a row, and delete a row.

A good overview of dataframes capabilities is available in our [API Samples](https://dev.datagrok.ai/js/samples/data-frame/).
Also check the API reference at [this link](https://datagrok.ai/js-api/DataFrame).

### Semantic annotation and metadata

### Aggregations and joining

## Persisting data

## Computations

## Visualizations

## Server API

## Managing privileges

### Groups and sharing

### Obtaining user info

### Managing credentials

We strongly advise you not to store service credentials directly inside your application packages
in any form, whether they reside as a part of a connection string, a parameter inside a connection file in 
Datagrok, or just login and password used in a POST request in your JS code.

Of course, it is suitable if one stores credentials for a database or a service from public domain
(for instance, of publicly available datasets like Chembl) right inside the package in a connection file.
For other types of credentials, there are suitable means in Datagrok described below.

There are two key facilities for managing your application credentials.

#### Pushing credentials by the Datagrok Server API

It's possible to programmatically push credentials to Datagrok and deliver them in to the connection
of interest. In such scenario, credentials are stored on a secured machine, and on demand,
e.g. through a deployment process, are delivered to Datagrok via a triggered bat/sh-script 
([link](./manage-credentials.md#database-connection-credentials)). For using it, you need
to provide an API developer's key, which is available in your user info pane in the Datagrok UI.
Access it with the user avatar button in the activity bar on the left side of Datagrok main window.

#### Storing credentials in the Datagrok Credentials Store

This is useful for credentials which are not part of Datagrok's OpenAPI web connection,
database or file share connection, but instead are used programmatically in JavaScript code of
the application to access 3-rd party REST APIs and the like ([link](./manage-credentials.md#package-credentials)).

In Datagrok, Credentials Store gives access to per-package key-value stores and accessible from package's code.
Credentials Store is a physically separate entity and may be placed on a separate machine,
which improves theft tolerance. It's also possible to deliver credentials to the Store programmatically
same way as discussed in p.1.

## UI and UX

Most of the UI capabilities Datagrok offers are described as samples in our [ApiSamples package]().
You may view them all conveniently by the following link: [https://public.datagrok.ai/js](https://public.datagrok.ai/js).

The main idea behind UI composition in Datagrok is straightforward: by minimal means, give to both non-experts
and professionals in JavaScript an ability to rapidly compose interfaces which look well by default. With this
goal in mind, along with a target of high performance and low latency of the UI, we've built our UI library 
from scratch instead of relying on one of the existing frameworks.

For most HTML elements such as buttons or dropdowns, our library is very lightweight, our classes such as
`ui.button` or `ui.choiceInput` are just tiny wrappers around regular DOM elements. In fact, you can always
access an underlying DOM object of our UI element with the `.root` property. There are also more advanced
controls not available in browsers out of the box, such as [panes](), [accordions]() and [dock views]().

With that being said, our customers do use other frameworks, such as React and Boostrap, in their applications
built on Datagrok. The final choice is up to the developer building on top of Datagrok.

This [JavaScript API help page](./develop/js-api.md) also gives a good idea of our UI/UX capabilities.

We are currently re-thinking our approach to UI composition with making it even simpler and more aligned
to some best practices of existing frameworks. For instance, we are going to provide for easy to understand
and easy to use controls catering for several typical application layouts, where the design has
to do with a choice between stretching and scrolling, fixed and dynamic sizing, and so forth.
Soon there comes an updated piece of documentation and samples on the subject.

## Working with packages

## Application lifecycle

## Debugging applications

Debugging applications follows same principles as to other packages, described [here](../develop.md#debugging).

As a note specifically for application developers, we recommend pointing your start URL in VS Code or other
IDE directly to the URL of your application: `https://<HOST_NAME>/apps/<APP_NAME>` for 1 app in the package
and `https://<HOST_NAME>/apps/<PACKAGE_NAME>/<APP_NAME>` for more than 1 app in the package. For VS Code,
this is done in `launch.json`.


See also:

  * [Grok JS development](develop.md)
  * [Developing grok applications](develop/develop.md#applications)
  * [Applications on Datagrok Public](https://public.datagrok.ai/apps)
  * [Development samples gallery](https://public.datagrok.ai/js)
