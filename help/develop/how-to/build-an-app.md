<!-- TITLE: Building an application -->
<!-- SUBTITLE: -->

<!-- This is a developer's view on the Datagrok applications -->

# Grok Application

Applications are built on top of the Datagrok platform. The Datagrok application serves a targeted
need around a particular problem or an area. For example, this may be a bioreactor modeling application,
an application for working with and understanding molecular sequence data, a molecular database browser,
a Covid-19 or weather info panel, etc.

Datagrok applications are developed in JavaScript / TypeScript using our rich [Datagrok JavaScript API](),
with also using any of the [scripting]() languages we support, such as R or Python, via calling scripts, and
connecting to third-party data sources and services.

Technically, applications are [functions](../overview/functions/function.md) tagged with the `#app` tag.
In this function, you take control over what the platform will let user do and see next once the app is
executed. Think of it as of the application's entry point (`Main` function and the like).

## Entry point

A Datagrok [package]() might contain zero, one, or more Datagrok applications. These come along any other
entities in the package, which these applications may be using, such as connections, viewers, scripts, etc.

Let's consider a very simple example of a webpack-based package with just two trivial apps in one `src/package.js`:

```
import * as grok from 'datagrok-api/grok';

export let _package = new DG.Package();

//name: TestUI1
//tags: app
export function testui1() {
  grok.shell.info('Test 1');
}

//name: TestUI2
//tags: app
export function testui2() {
  grok.shell.info('Test 2');
}
```

To make this run on Datagrok, follow these `grok create` [steps](../develop/develop.md#getting-started)
to prepare such a simple package and deploy it.

In very short:

1. With the first use of `grok tools` after you install it, you'll be prompted to enter your API developer's token to allow pushing packages to Datagrok; this key will be stored inside `.grok` settings file on your drive
2. You'd choose to create a package via `grok create <PACKAGE_NAME>`; this will initialize a standard directory structure for your package
3. Then you may choose to add an app to it via `grok add app <APP_NAME>` and see the minimal default structure for it, or alternatively just copy-paste the JS snipped from the above to `package.js`

After deploying this package, you'd find these 2 apps via `Functions | Apps` in the activity bar
situated on the left side of Datagrok's main window. Run both of the apps and notice the two different
tooltips popping up. You'd also call these same entry points by an URL:
`https://public.datagrok.ai/apps/<PACKAGE_NAME>/TestUI1` and a similar one for `TestUI2`.
Note that, in case there is only one application `<APP>` defined in the package, the
corresponding URL will be simply `https://public.datagrok.ai/apps/<APP>`, omitting the
`<PACKAGE_NAME>` part.

Read more on creating custom views [here](./custom-views.md).

This simple example finishes explaining the purpose of the entry points. Yet, trivial popups aren't
something one typically builds as an application. Let's look at a more UI-rich side of things.

## Main view

Most applications built on Datagrok start with a view. A view is a set of visualizations grouped together.
Typically, it is associated with a particular dataframe, then it is called a table view. However, essentially
a view can contain pretty much anything.

Imagine you are composing a coherent application.
You'd most likely start off with the root / main view, add logical blocks to it either through simple
[div-s](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/ui/sidebar.js), or through
[`splitH`/`splitV`](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/ui/layouts/splitters.js),
populate these blocks with visualizations and controls, maybe add some sidebard, add event handlers and so forth.
Our internal application [Usage Analysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis)
demonstrates such approach.

Another approach is found in a [Discovery](https://github.com/datagrok-ai/public/tree/master/packages/Discovery) application.
There we reuse a particular kind of view: a table view, and centralize the rest of the UI around it. In
this application you'd also find some useful techniques one needs sooner or later in many applications,
such as modifying the app's URI or hiding the side panels of the Datagrok's main UI.

With that being said, there is pretty much anything you can further do inside the application, leveraging
a full scale of platform capabilities. However, there are certain things of utmost interest in almost any
application:

* accessing and persisting data
* working with [dataframes]()
* performing computations
* interactive visualizations
* managing privileges
* application lifecycle

This is what we are going to cover further in this how-to.

# Application Development

This chapter serves as a guide to the detailed articles on key topics around applications building.

## Code samples

We provide a diversed set of code snippets of the API use, and sample packages with viewers, applications, etc.

* For short samples of using API, go to https://public.datagrok.ai/js and observe "Samples" block, or alternatively access it via a "Help" button at the bottom of the activity bar on the left of the Datagrok's main window (then follow to `JavaScript API Samples`).

* The sources of these snippets are all located at https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples.

* Observe the entirety of our demo codes at `https://github.com/datagrok-ai/public/tree/master/packages`. For instance, you may find there a package [Viewers](https://github.com/datagrok-ai/public/tree/master/packages/Viewers) which showcases a variery [Chem package](https://github.com/datagrok-ai/public/tree/master/packages/Chem) 

These developers who's been starting with the platform recently are telling us the samples above became a major
source of knowledge for their daily work with the platform. We recommend you to recall these locations too,
and resort to them for resolving your daily technical questions in addition to posting questions and suggestions
on our [Community Forum](https://community.datagrok.ai/) at https://community.datagrok.ai/.

## Data Access

There's a variety of data sources which Datagrok can handle out of the box:

* Web Services (REST endpoints) specified via Swagger / OpenAPI ([link]())
* Arbitrary REST APIs (including these outside host's domain) ([link]())
* Access a file from the Datagrok server space (either shared or Home) ([link]())
* Reading from databases: Datagrok supports more than 30 ([link]())
* Access a file from a file share (Windows Network, SAMBA and the like) ([link]())

Connections to databases, respective database queries, and file shares live in Datagrok as entities.
These may be published to only you, the package author, or to any specific user group.
Main code work for accessing databases is concentrated around `grok.data.query` method, whereas
accessing files from file shares and servers is made by `grok.functions.eval(`OpenServerFile("..."))`. 

## Persisting data

## Computations

## Visualizations

## Server API

## Managing privileges

## UI and UX

## Working with packages

### Package structure

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