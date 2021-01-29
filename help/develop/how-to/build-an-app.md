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

To make this run on Datagrok, follow `grok create app` [steps](../develop/develop.md#getting-started)
to prepare this simple package and deploy it.

After deploying this package, you'd find these 2 apps via `Functions | Apps` in the activity bar
situated on the left side of Datagrok's main window. Run both of the apps and notice the two different
tooltips popping up. You'd also call these same entry points by an URL:
`https://public.datagrok.ai/apps/<PACKAGE_NAME>/TestUI1` and a similar one for `TestUI2`.
Note that, in case there is only one application `<APP>` defined in the package, the
corresponding URL will be simply `https://public.datagrok.ai/apps/<APP>`, omitting the
`<PACKAGE_NAME>` part.

This simple example finishes explaining the purpose of the entry points. Yet, tooltips aren't
something one typically builds as an application. So, let's proceed to a more UI-rich kind of things.

## The main view

# Application Development

This chapter serves as a guide to the detailed articles on key topics around applications building.

## Data Access

## Computations

## Visualizations

## Server API

## Managing privileges

## UI and UX

## Working with packages

# Debugging applications

See also:

  * [Grok JS development](develop.md)
  * [Developing grok applications](develop/develop.md#applications)
  * [Applications on Datagrok Public](https://public.datagrok.ai/apps)
  * [Development samples gallery](https://public.datagrok.ai/js)