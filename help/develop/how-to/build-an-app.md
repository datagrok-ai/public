<!-- TITLE: Building an application -->
<!-- SUBTITLE: -->

<!-- This is a developer's view on the Datagrok applications -->

# A Grok Application

Applications are built on top of the Datagrok platform. The Datagrok application typically serves a targeted
need around a particular problem or area. For example, this may be a bioreactor modeling application,
an application for working with and understanding molecular sequence data, molecular database browser,
a Covid-19 or weather info panel.

Datagrok applications are built in JavaScript or TypeScript using our rich [Datagrok JavaScript API]() 

Technically, applications are [functions](../overview/functions/function.md) tagged with the `#app` tag.
In this function, you take control over what the platform will let user do and see next once the app is
executed. Think of it as of the application's entry point (`Main` function and the like).

# Entry points

A Datagrok package might contain zero, one, or more Datagrok applications. These come along any other
entities, which these applications may be using, such as connections, viewers, scripts, etc.

Let's consider a very simple example of a webpack-based package with two apps in one `package.js`:

```
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

Follow `grok create app` [steps](../develop/develop.md#getting-started) to prepare this simple package.

# Application Development

# Debugging applications

See also:

  * [Grok JS development](develop.md)
  * [Developing grok applications](develop/develop.md#applications)
  * [Applications on Datagrok Public](https://public.datagrok.ai/apps)
  * [Development samples gallery](https://public.datagrok.ai/js)