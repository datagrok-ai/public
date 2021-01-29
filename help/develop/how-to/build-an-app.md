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

## Structuring code





# Application Development

With the above being said, there is pretty much anything you can further do inside the application, leveraging
a full scale of platform capabilities. However, there are certain aspects of utmost interest in almost any
application:

* accessing and persisting data
* working with [dataframes]()
* performing computations
* interactive visualizations
* managing privileges
* application lifecycle

The following chapter guides through these key development topics. Take it as a birds-eye overview of the
applcation development area, and proceed to the articles and samples referenced in the guide for the
further details.

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

* Web Services (REST endpoints) specified via Swagger / OpenAPI ([link](../access/open-api.md))
* Arbitrary REST APIs (including these outside host's domain) ([link](./access-data.md#rest-endpoints))
* Access a file from the Datagrok server space (either shared or Home) ([link](./access-data.md#file-shares))
* Reading from databases: Datagrok supports more than 30 ([link](./access-data.md#parameters))
* Access a file from a network file share (Windows Shares, SAMBA and the like) ([link](./access-data.md#file-shares))

Connections to databases, respective database queries, and file shares live in Datagrok as entities.
These may be published to only you, the package author, or to any specific user group.
Main code work for accessing databases is concentrated around `grok.data.query` method, whereas
accessing files from file shares and servers is made by `grok.functions.eval('OpenServerFile("...")')`.

We strongly encourage you not to store credentials directly inside your application packages.
The chapter ["Managing privileges"](#managing-privileges) discusses credentials management in detail.

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
Datagrok, or just login and password used in a POST request in JS code.

Of course, it is suitable if one stores credentials for a database or a service from public domain
(for instance, of publicly available datasets like Chembl) right inside the package in a connection file.
For other types of credentials, there are suitable means in Datagrok.

There are two key facilities for managing your application credentials.

1. *Pushing credentials via REST API using a developer's key* ([link](./manage-credentials.md#database-connection-credentials)). It's possible to programmatically push credentials to Datagrok and deliver them in to the database connection of interest. In such scenario, credentials are stored on a secured machine, and on demand, e.g. through a deployment process, are delivered to Datagrok via a triggered bat/sh-script.

2. *Storing credentials in Datagrok's Credentials Store* ([link](./manage-credentials.md#package-credentials)). This is useful for credentials which are not part of Datagrok's OpenAPI web connection, database or file share connection, but instead are used programmatically in JavaScript code of the application to access 3-rd party REST APIs and the like. In Datagrok, Credentials Store gives access to per-package key-value stores available per package and accessible from package's code.  
Credentials Store is a physically separate entity and may be placed on a separate machine, which improves theft tolerance. It's also possible to deliver credentials to the Store programmatically same way as described in p.1.

## UI and UX

Most of the UI capabilities Datagrok offers are described as samples in our [ApiSamples package]().
You may view them all conveniently by the following link: https://public.datagrok.ai/js.

The main idea behind UI composition in Datagrok is very simple: by minimal means give to both non-experts
and professionals in JavaScript an ability to rapidly compose interfaces which look well by default. With this
goal in mind, along with a target of high performance and low latency of the UI, we've built our UI library 
from scratch instead of relying on one of the existing frameworks.

For most HTML elements such as buttons or dropdowns, our library is very lightweight, our classes such as
`ui.button` or `ui.choiceInput` are just tiny wrappers around regular DOM elements. In fact, you can always
access an underlying DOM object of our UI element with `.root` property. There are also more advanced
controls not available in browsers out of the box, such as [panes](), [accordions]() and [dock views]().

With that being said, our customers do use other frameworks, such as React and Boostrap, in their applications
built on Datagrok. The final choice is up to the applied developer.

This [JavaScript API help page](./develop/js-api.md) also gives a good idea of our UI/UX capabilities.

We are currently re-thinking our approach to UI composition to make it even simpler and more aligned
to best practices of existing frameworks. For instance, we are going to provide for easy to understand
and easy to use controls carrying of several typical application layouts, where the design has
to do with a choice between stretching and scrolling, fixed and dynamic sizing. Soon there comes
an updated piece of documentation and samples on the subject.

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