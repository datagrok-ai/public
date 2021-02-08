<!-- TITLE: Building an application -->
<!-- SUBTITLE: -->

<!-- This is a developer's view on the Datagrok applications -->

# Datagrok Applications

Applications are built on top of the Datagrok platform. The Datagrok application serves a targeted
need around a particular problem or an area. For example, this may be a bioreactor modeling application,
an application for working with and understanding molecular sequence data, a molecular database browser,
a Covid-19 or weather info panel, etc.

Datagrok applications are developed in JavaScript / TypeScript using our rich [Datagrok JavaScript API](),
with also using a range of the [scripting]() languages we support, such as R or Python, and connecting
to third-party data sources and services for data integration and orchestration.

Technically, applications are [functions](../overview/functions/function.md) tagged with the `#app` tag.
In this function, you take control over what the platform will let user do and see next, after the app is
executed. Think of it as of the application's entry point, such as a `Main` function in a C# console app.

## The Entry Point

The Datagrok platform is highly extensible. New functionality is delivered to a Datagrok instance as packages.
A Datagrok [package]() might contain zero, one, or more Datagrok applications. These come along with other
entities in the package, which the applications may be using, such as [connections](),
[viewers](), [scripts](), etc.

Consider a simple example of a webpack-based package with just one trivial app in a `src/package.js`:

```
import * as grok from 'datagrok-api/grok';

export let _package = new DG.Package();

//name: TestApp
//tags: app
export function test() {
  grok.shell.info('An app test');
}
```

To make this run on Datagrok, follow these `grok create` [steps](../develop/develop.md#getting-started)
to prepare this simple package and deploy it.

1. Install the [prerequisites](../develop.md#getting-started): 
   * A regular [Node.js](https://nodejs.org/en/download/) and [npm](https://docs.npmjs.com/about-npm) (comes with Node.js)
   * `npm install webpack -g` (`-g` will make `webpack` globally available, that's convenient) 
   * `npm install webpack-cli -g`
   * `npm install datagrok-tools -g`
2. [Create a new package](../develop.md#getting-started):
   * Make a new folder for the package
   * In this folder, call `grok create <PACKAGE_NAME> --ide=vscode`
   * The `--ide` key will create a setup for debugging the package with VS Code (currently only works on Windows)
   * If you run this for the first time, you'd be prompted to enter your Developer Keys for
     our Datagrok instances. Find this key in your [User Profile]() section in the Datagrok UI
3. Add an app to the package by `grok add app <APP_NAME>`, or just copy-paste the above JS snippet into `src/package.js`

After deploying this package to `https://public.datagrok.ai`, you'd find the `Test App` app via `Functions | Apps`
in the activity bar on the left side of Datagrok's main window. Run the app and notice the tooltip.
You may also call this entry point by an URL: `https://public.datagrok.ai/apps/TestApp`.
Note that in case there is only one application `<APP>` defined in the package, the corresponding URL
is simply `https://public.datagrok.ai/apps/<APP>`, but has the form of
`https://public.datagrok.ai/apps/<PACKAGE_NAME>/<APP_NAME>` for 2 and more apps in one package.

This simple example concludes on the entry point. Yet, trivial popups aren't
something one typically builds as an application. Let's look at a more UI-rich side of things.

## The Main View

Most applications built on Datagrok start with a Datagrok's view. A view is a set of visualizations and
controls grouped together. Typically, the view is associated with a particular [dataframe](#dataframes),
in this case it's called a [Table View](). However, essentially a view can contain pretty much anything.

Imagine you are composing an application. You likely start with the root / main view,
add logical blocks to it either through simple [div-s](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/ui/sidebar.js),
or through [`splitH`/`splitV`](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/ui/layouts/splitters.js),
populate these blocks with visualizations and controls, maybe add a sidebar, add event handlers and so forth.
Our internal application [Usage Analysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis)
demonstrates such an approach.

Another approach is found in a [Discovery](https://github.com/datagrok-ai/public/tree/master/packages/Discovery) application.
There we reuse a particular kind of view: a table view, and centralize the rest of the UI around it. In
this application you'd also find some useful techniques one needs in many applications,
such as modifying the app's URI or hiding the side panels of the Datagrok's main UI.

Read more on creating custom views [here](./custom-views.md).

## Working with IntelliSense

We recommend restoring the package dependencies before starting development with an IDE.
After the package is created, simply invoke `npm install` inside the package folder. This
will bring NPM modules defining the Datagrok API.

An alternative way to IntelliSense capability for Datagrok classes is by cloning
the entire [public repo](https://github.com/datagrok-ai/public) and opening its whole folder
in the IDE. This is suitable when you develop a package which is
a [part](https://github.com/datagrok-ai/public/tree/master/packages) of our public repo.
The folder `js-api` in the root of the public repo is a source of our JS API, that's how
IDE discovers needed IntelliSense information.

# Application Development

Starting with the above, there is pretty much anything you can further do inside the application,
leveraging a full scale of platform capabilities. However, there are certain aspects of interest
in almost any application:

* Accessing and persisting data
* Working with dataframes
* Performing computations
* Interactive visualizations
* Managing privileges
* Managing application lifecycle

The following chapter guides through these key development topics. Take it as a birds-eye overview of the
applcation development area, grasp the major building blocks, and proceed to the articles and samples
referenced in the guide for the further details.

## Code Samples

We provide a diversed set of code snippets of the API use, and sample packages with viewers, applications, and so forth.

* For short samples of using API, go to https://public.datagrok.ai/js and observe a "Samples" block, or alternatively —
access it via the "Help" button at the bottom of the activity bar on the left of the Datagrok's main window,
and then follow to `JavaScript API Samples`

* The sources of these snippets are all located at https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples

* Check all our demo codes at `https://github.com/datagrok-ai/public/tree/master/packages`. For instance, you may
find there a package [Viewers](https://github.com/datagrok-ai/public/tree/master/packages/Viewers) which showcases
a variery of viewers one can build on Datagrok, or a [Chem package](https://github.com/datagrok-ai/public/tree/master/packages/Chem),
which shows how in-browser computations with WebAssembly may be organized on the platform

As some experience was shared to us, the samples above became a major source of knowledge for newcomers'
daily work with the platform. We recommend you to recall these locations too, and resort to them
in your daily technical questions. That is, in addition to posting questions and suggestions
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
Main code works for accessing databases are concentrated around `grok.data.query` method ([a sample]()), whereas
accessing files from file shares and servers is made by `grok.functions.eval('OpenServerFile("...")')` ([a sample]()).

When it comes to connections, we strongly encourage you not to store their credentials directly inside your
application packages. The chapter ["Managing privileges"](#managing-privileges) discusses credentials management in detail.

### Accessing REST endpoints outside the host

Web services provide endpoints that you can programmatically connect to.
There are two main options for this: the first is to use
[OpenAPI/Swagger](https://swagger.io/docs/specification/about/) format supported
by Datagrok, the second one involves the
[use of the platform's server](https://dev.datagrok.ai/js/samples/dapi/fetch) to send a network request.
The details of Swagger-based connections are further explained in the
[dedicated article](../../access/open-api.md).

It's common for JavaScript applications to query external REST endpoints. Due to
[CORS](https://en.wikipedia.org/wiki/Cross-origin_resource_sharing),
the Web application won't be allowed to GET/POST with resources on different hosts.
To overcome this, you don't need to move data querying to the backend server yourself.
Datagrok provides for a proxy server to work with such REST end-points. The method `fetchProxy` has
an interface similar to the standard `fetch` API and described
[here](develop/how-to/access-data.md#rest-endpoints)
in detail.

*References:*

* [OpenAPI in Datagrok](access/open-api.md)
* [fetchProxy sample](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/dapi/fetch.js)
* [Access Data](develop/how-to/access-data.md)

## Dataframes

### Creating and working with Dataframes

Dataframes are at the heart of Datagrok. Dataframe is a high-performance, easy to use tabular structure with
strongly-typed columns of different types. Supported types are: string, int, float, bool, DateTime, bigint and
[qualified numbers]().

Datagrok dataframes are highly optimized. They are implemented with the proprietary, unique technology allowing
to efficiently work with huge datasets in the browser. Essentially, it is a columnar in-memory database that was
engineered from scratch and optimized for the purpose of exploratory data analysis, interactive visualizations,
and machine learning.

Note that Datagrok dataframes live and operate entirely inside the browser, but not on our [compute server]().
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
let col = df.columns.byName('City');
col.set(0, 'New York');
df.set('City', 1, 'Chicago');
df.set('City', 2, 'San Francisco');
df.rows.addNew(['John', 175, 'Dallas']);
df.rows.removeAt(0);
grok.shell.addTableView(df);
```

There we create a new dataframe from columns, add a column and a row, delete a row, and open the result
in a table view.

A good overview of dataframes capabilities is available in our [API Samples](https://dev.datagrok.ai/js/samples/data-frame/).
Also check the API reference at [this link](https://datagrok.ai/js-api/DataFrame).

### Semantic annotation and metadata

Most of the objects in Datagrok can be annotated with metadata (key-value pairs). The
metadata could be set manually; additionally, some of it gets assigned automatically.
Some keys affect the way an object (such as a column) interacts with the platform; other have no
effect at all, except that you can search objects by metadata.

There is a variery of metadata in Datagrok, discussed [here](discover/metadata.md).
Out of all metadata, column tags and semantic types are of particular interest in
application development and dataframes.

#### Column tags

Each column in a Datagrok dataframe has a `.tags` property, which is a JS `Map` you can read and write to.
You could use these `.tags` for any package-specific activity needed. For example, if you have a pair
of columns containing [smiles strings]() of molecules, and one column serves as a scaffold to which
the other column should be aligned to, it is sensible to put a `scaffold-col` tag on a column being aligned,
referencing a column with the scaffold. The Datagrok's [custom cell renderer]() would make use of this tag.

It's possible to [add and remove tags]() from the Datagrok UI. To edit column's metadata,
right-click on it and select "Properties..." (or press F2 in the grid).

*References:*

* [Tags](discover/tags.md)

#### Semantic types

Unlike in Excel, table columns in Datagrok are strongly-typed. In addition to that, a column
might also have semantic type associated with it. Semantic type identifies the meaning of the data.
For instance, a column could be of the data type "string", but the semantic type would be "country".

Semantic types are used in several ways. For example, OpenChemLib or RDKit [are used]()
for rendering string columns of a semantic type `'Molecule'`. Some viewers, such as [Map](),
use the semantic type to determine whether a viewer can be visualized against specific datasets.
Function parameters [could be annotated]() with the semantic type. This is  used for automatic
suggestions of applicable functions.

Semantic type is a special kind of a [column tag](tags.md#semantic-type).
It could be either detected automatically by column [semantic type detectors](), or set manually
in the JavaScript code with `.semType` property.

*References:*

* [Semantic Types](discover/semantic-types.md)

### Aggregations and joining

For those of you familiar with functional style in JavaScript (`.map`, `.filter`, etc),
with LINQ in C# or other kinds of fluent APIs, it should be straightforward to see
what happens to a dataframe in this snippet:

```
let agesByRace = grok.data.demo.demog()
  .groupBy(['race', 'sex'])
  .where('race = Asian and started after 2/2/2000')
  .avg('age')
  .avg('started')
  .add('med', 'age')
  .aggregate();

grok.shell.addTableView(agesByRace);
```

We first load the `demog` dataset, filter it, and for a compound grouping by age and sex
compute the average for `age` and `started`, and also a median for `age`.

The whole set of functions available for `.add` is located [here](transform/aggregation-functions.md).

*References:*

* [Aggregation functions](transform/aggregation-functions.md)

## Persisting Data

### User data storage

Often application settings or its inputs/outputs need to be shared between different applications,
different instances of the same application, or simply persisted for later use or autoload
on application start. This functionality is available in Datagrok as user data storage —
a Datagrok server storage which can be filled with new entries and from which these entries
can later be retrieved.

The data resides in the storage as a set of stores, each identified by a unique name,
with a key-value map placed in each store. There are several _asyncronous_ methods for storing
and retrieveing data from the user storage, such as `grok.dapi.userDataStorage.postValue`
for posting a single value, or `grok.dapi.userDataStorage.get` for getting the whole map.
Learn of all these methods [here](develop/user-data-storage.md), also check a complete example
in [API Samples](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/misc/user-data-storage.js).

As the name assumes, the storage is only seen to the user. However, it's possible to create
stores visible to all users. This is controlled by the last argument in the above methods,
a boolean `currentUser`, which is set to `true` by default.

*References:*

* User data storage [reference](develop/user-data-storage.md)
* User data storage [sample](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/misc/user-data-storage.js)

### Storing tables

## Computations

### Scripting

Beside application logic enabled on the client side (in the browser) with the Datagrok JavaScript API, 
Datagrok provides scripting, a mechanism for integration with languages for statistical/mathematical
computing.  Scripting combines fast interactive visualizations and other features of the Datagrok platform 
with thousands of statistical packages and visualizations available in  [R](https://www.r-project.org/about.html),
[Python](https://www.python.org), [Octave](https://octave.org/), [Julia](https://julialang.org),
or [JavaScript](https://www.javascript.com).

Let's look at a simple script. Create an empty Python script with `Functions | Scripts | New Python Script`
(find a `Functions` button on the left side of the Datagrok main window in the activity bar). Before
running, open some table, you may find a lot of them in the "Demo Files" section of the platform through
the `Data` pane opened by the "folder" button on the left side of the Datagrok's window. Then let's run
the below:

```python
#name: SimpleTestPython
#language: python
#input: dataframe table [Data table]
#output: int count [Number of cells in the table]
count = table.shape[0] * table.shape[1]
```

This script simply calculates a number of cells in the dataframe passed to the script.
Note how the script markup in the header comments specify this script as a function
with a signature of a name, typed inputs and outputs. This is a general notion of a function
you'd find across the entire platform.

Before reading more about Datagrok functions, let's check the output of the script.
To see it, open Datagrok's console. As in Doom, you may find it by hitting a `~`
anywhere inside Datagrok. There you can see the script's output, something like that:

```
> SimpleTestPython("demog")
  count: 64350
```

Stay in this console for a little longer and try executing this script again with the command:

```
SimpleTestPython("demog")
```

When the script is run, here is what happens under the hood:

* The dataframe and all other input parameters are serialized and sent to the Datagrok server
* The Datagrok server is hosting a so-called [compute virtual machine]() with a ready-to-execute
  instance of a [Jupyter Kernel]() for Python, as well as for [other]() supported scripting languages
* When the request to execute a script arrives to the server along with its parameters,
  CVM loads scripts' code from its storage and runs it with these parameters
* The output values are assigned anywhere along the script
* The script execution is fully stateless and isolated

You can even return graphics from the script! Check it with [this excercise on Scripting Viewers]().

### Datagrok functions

### Calling scripts from JavaScript code

## Visualizations

Datagrok provides for rich data visualization with more than 25+ viewers out of the box, including
[Scatter Plot](visualize/viewers/scatter-plot.md), [Histogram](visualize/viewers/histogram.md),
Line Chart, Bar Chart, Pie Chart, Trellis, Matrix Plot, 3D Scatter, Density Plot,
PC Plot, Word Cloud, Network, Box Plot, Tree Map, Heat Map, Statistics, Correlation, Calendar, Table Grid,
Markup, Tiles, Form, Map, Shape, Chord, and Tree. These viewers were crafted for web from scratch, and are
purpose-fit for life science applcations. For example, some viewers [provide for](visualize/viewers/box-plot.md)
built-in statistical hypothesis testing and Multivariate Analysis. Based on our unique in-browser
[data store](#creating-and-working-with-dataframes), they cater for high performance data exploration
on datasets of 10 million rows.

Viewers may be connected to other visual entities of your applciation, providing for interactive
filtering, selection, current row, mouse-over rows and other events. All existing viewers
are available both ad-hoc while exploring data and programmatically in applications' code.
It's possible to expand the platform with new viewers shared through versioned packages.
Viewers may be created in both native way in JavaScript code with Datagrok JS API, and via
R, Python, Julia-based visualizations. 

Viewers are either supplied with their own dataframes, or (more typical) are linked to
the dataframes objects which are re-used through the application's code.

*References:*

* [Viewers](visualize/viewers.md)

## Managing privileges

Datagrok has an enterprise-grade support for priveleges management and secure deployment.
In the context of Datagrok application development, this covers both fine-grained access control
to all Datagrok entities and facilities, but also security patterns for storing credentials
and authorized access to Datagrok using popular authentication protocols.

### Authorization, groups and sharing

Datagrok has a flexible mechanism for grouping users together. A user can belong to more than one group.
A group can be included in another group, which is useful for both reflecting organization hierarchy
and implementing  role-based [security](security.md). 

Many types of objects within the Datagrok platform can be shared with other users or
[groups](../govern/group.md). Such shareable objects are called [entities](../overview/objects.md).
When an object is shared, you are essentially granting a [privilege](../govern/authorization.md)
(typically, 'view' or 'edit') to a grantee. See the [Security](../govern/security.md) article
for details.

When launching a new application on Datagrok, we recommend considering a group which members
should have access to this application. Often this is a new group which you'd create for 
this specific application. With the application along come the queries and connections used in it.
This means grating access to an applicaiton implies granting access to the whole package.

Access to either whole package with all its content, or more fine-grained items such as applications
or connections, is managed per-user and per-group. Therefore, letting a user group access your applcaiton
usually means sharing the corresponding package to this group. When granting access to a package for
a particular group, access to all its items is also shared to this group.


*References:*

* [Sharing](/colaborate/sharing.md)
* [Groups](/govern/groups.md)
* [Security](/govern/security.md)
* [Authorization](/govern/authorization.md)

### Authentication

*References:*

* [Authentification](/govern/authentication.md)

### Obtaining groups and users info

It is possible to access groups info through Datagrok JavaScript API. The namespace `grok.dapi.groups`
provides for it. Find code snippets for this topic in
[/dapi/groups](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples/scripts/dapi)
part of our JS API samples. 

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

*References:*

* [Managing credentials](./develop/how-to/manage-credentials.md)

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

## Subscribing to events

We are exposing events coming out of the platform as a stream via the [Rx.JS](rxjs.dev) library
that makes it easy to compose asynchronous or callback-based code. The API makes easy to subscribe
to either global, or instance-related events:

```javascript
let v = grok.shell.newView('Range Slider');
let rangeSlider = ui.rangeSlider(0, 10, 2, 5);
rangeSlider.onValuesChanged.subscribe(
  (_) => grok.shell.info(`Range Slider changed to [${rangeSlider.min}, ${rangeSlider.max}]`));
v.append(rangeSlider);
```

Paste it to the [JavaScript fiddle](https://public.datagrok.ai/js), drag the range slider
to see the tooltip with current bounds values.

There are other styles of events subscription offered in Datagrok APIs. They are introduced for
brevity in places where it's natural. For instance, here we create a dialog and subscribe
to pressing an `OK` button:

```javascript
ui.dialog('Windows')
  .add(ui.span(['A message']))
  .onOK(() => grok.shell.info('OK!'))
  .show();
```

To figure out what events are coming out of the platform, use the Inspector tool.
Open it (`Alt+I`), go to the "Client Log" tab, and perform the action that you want
to intercept. In the panel, you will see one or more of the events, click on them
to inspect event parameters. To simplify the development process, we also generate
JavaScript code for handling this particular event, copy-paste it from the
property panel into your code if needed. 

Remember that the underlying UI controls in Datagrok are still DOM elements,
so you can always subscribe to the regular events using `.addEventListener`:

```javascript
let v = grok.shell.newView('Events');
let button = ui.bigButton("Run");
v.append(button);
button.addEventListener("click", () => alert("Hello!"));
```

This is useful in case you decouple interface specification and the event handlers
implementation into two separate files. However, we are homogenizing the current state
of UI events subcription so that it's always possible to subcribe using rxjs
`.subscribe` method on an event object.


Read more about Datagrok events [here](develop/js-api.md#Events).

*References:*

* [Events](develop/js-api.md#Events)
* [Global events](https://public.datagrok.ai/js/samples/events/global-events)
* [DataFrame events](https://public.datagrok.ai/js/samples/data-frame/events)

## Working with packages

### Structuring code

Prehaps, one of the main things to know about webpack is that it allows you to structure JavaScript code
in a way similar to how you are used to it in enterprise-grade environments, such as Java or .NET.

Let's expand our previous example leveraging an `import` capability of webpack and modern JavaScript.

Create two separate files:

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
