<!-- TITLE: Building an application -->
<!-- SUBTITLE: -->

<!-- This is a developer's view on the Datagrok applications -->

# A Grok Application

Applications are built on top of the Datagrok platform. The Datagrok application serves a targeted
need around a particular problem or an area. For example, this may be a bioreactor modeling application,
an application for working with and understanding molecular sequence data, a molecular database browser,
a Covid-19 or weather info panel, etc.

Datagrok applications are developed in JavaScript / TypeScript using our rich [Datagrok JavaScript API](),
with also using any of the [scripting]() languages we support, such as R or Python, via calling scripts.

Technically, applications are [functions](../overview/functions/function.md) tagged with the `#app` tag.
In this function, you take control over what the platform will let user do and see next once the app is
executed. Think of it as of the application's entry point (`Main` function and the like).

## Entry points

A Datagrok [package]() might contain zero, one, or more Datagrok applications. These come along any other
entities in the package, which these applications may be using, such as connections, viewers, scripts, etc.

Let's consider a very simple example of a webpack-based package with two apps in one `package.js`:

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

To make this run on Datagrok, follow these `grok create app` [steps](../develop/develop.md#getting-started)
to prepare this simple package and deploy it.

In very short:

1. You'll be prompted to enter your API developer's token to allow pushing packages to Datagrok; this key will be stored inside `.grok` settings file on your drive
2. You'd choose to create a package via `grok create <PACKAGE_NAME>`; this will initialize a standard directory structure for your package
3. Then you may choose to add an app to it via `grok add app <APP_NAME>` and see the minimal default structure for it, or alternatively just copy-paste the JS snipped from the above to `package.js`

After deploying this package, you'd find these 2 apps via `Functions | Apps` in the activity bar
situated on the left side of Datagrok's main window. Run both of the apps and notice the two different
tooltips popping up. You'd also call these same entry points by an URL:
`https://public.datagrok.ai/apps/<PACKAGE_NAME>/TestUI1` and a similar one for `TestUI2`.
Note that, in case there is only one application `<APP>` defined in the package, the
corresponding URL will be simply `https://public.datagrok.ai/apps/<APP>`, omitting the
`<PACKAGE_NAME>` part.

This simple example finishes explaining the purpose of the entry points. Yet, tooltips aren't
something one typically builds as an application. Let's proceed to a more UI-rich kind of things.

## The main view

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

## Computations

## Visualizations

## Server API

## Managing privileges

## UI and UX

## Working with packages

### Package structure

# Debugging applications

See also:

  * [Grok JS development](develop.md)
  * [Developing grok applications](develop/develop.md#applications)
  * [Applications on Datagrok Public](https://public.datagrok.ai/apps)
  * [Development samples gallery](https://public.datagrok.ai/js)