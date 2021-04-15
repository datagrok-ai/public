<!-- TITLE: Grok Application -->
<!-- SUBTITLE: -->

<!-- This is a user-centric view on the Datagrok applications -->

# Grok Application

A Datagrok application is built on top of the Datagrok platform.

It typically serves a targeted need around a particular problem or area. For example:

* customer-specific applications, such as a bioreactor modeling application or an application for working with and understanding molecular sequence data
* interactive applications and dashboards, such as a molecular database browser (see our [Enamine Store](https://github.com/datagrok-ai/public/tree/master/packages/EnamineStore) and [ChemblBrowser](https://github.com/datagrok-ai/public/tree/master/packages/ChemblBrowser)), a Covid-19 or a weather data info panel, and so forth
* integration applications, where data from several sources is interpreted and visualized in a coherent way

Datagrok applications are developed in JavaScript / TypeScript using our rich
[Datagrok JavaScript API](develop/js-api.md), with parts in scripting languages
we support, such as R or Python, by calling [Datagrok scripts](scripting.md).

Application is typically comprised of many [UI features](develop/ui.md) of the platform,
including [input controls](develop/ui.md), [accordion panes](develop/ui.md#accordions),
[grid viewers](develop/ui.md#grid), [scatter plots](visualize/viewers/scatter-plot.md),
[line charts](visualize/viewers/line-chart.md), and many other [viewers](visualize/viewers);
functional features, such as these implemented in Python and R [scripts](develop/scripting.md),
accessing external web services in JavaScript via
[OpenAPI specs]() and [REST API](), and many others.

One [Datagrok package](../develop/develop.md#packages) may contain zero, one or several applications.

# Launching applications

To open the application launcher and see all [available](https://public.datagrok.ai/apps) applications,
follow to `Functions | Apps` from the Datagrok sidebar, or follow [this link](https://public.datagrok.ai/apps).
To launch a particular app automatically, open the following URL: `https://public.datagrok.ai/apps/<APP_NAME>`,
if there is only one app in the package, or `https://public.datagrok.ai/apps/<PACKAGE_NAME>/<APP_NAME>`,
if there are several. You'd learn the right URL first time you run the app from within the application launcher; we recommend you to bookmark it for later use.

# Developing applications

There's a handful of concepts and patterns in Datagrok especially useful
for building applications. The application building story is curated in
[this guide](develop/how-to/build-an-app.md), which we highly recommend prior
to starting your first Datagrok application development.

See also:

  * [Grok JS development](develop.md)
  * [Developing grok applications](develop/develop.md#applications)
  * [Applications on Datagrok Public](https://public.datagrok.ai/apps)
  * [Development samples gallery](https://public.datagrok.ai/js)
  * [A guide: How to build an app](develop/how-to/build-an-app.md)