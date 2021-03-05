<!-- TITLE: Grok Application -->
<!-- SUBTITLE: -->

<!-- This is a user-centric view on the Datagrok applications -->

# Grok Application

A Datagrok application is built on top of the Datagrok platform.

The Datagrok application typically serves a targeted need around a particular problem or area. For example:

* customer-specific applications, such as a bioreactor modeling application or an application for working with and understanding molecular sequence data
* interactive applications and dashboards, such as molecular database browser (see our [Enamine Store application](https://github.com/datagrok-ai/public/tree/master/packages/EnamineStore)), a Covid-19 or weather data info panel, etc.
* integration applications, where data from several sources is interpreted and visualised in a coherent way

Datagrok applications are developed in JavaScript / TypeScript using our rich [Datagrok JavaScript API](), with using
any of the [scripting]() languages we support, such as R or Python, via calling scripts.

Application is typically comprised of many UI features of the platform, including [input controls](),
[accordion panes](), [grid viewers](), scatter plots, line charts and many other [viewers](../viewers/scatter-plot.md),
functional features, such as calling Python and R [scripts](../develop/scripting.md), accessing external web services
via [OpenAPI specs]() and [REST API](), and many others.

One [Datagrok package](../develop/develop.md) may contain zero, one or several applications.

# Launching applications

To open the application launcher and see all [available](https://public.datagrok.ai/apps) applications,
follow to `Functions | Apps` from the Datagrok sidebar, or follow [this link](https://public.datagrok.ai/apps).
To launch a particular app automatically, open the following URL: `https://public.datagrok.ai/apps/<APP_NAME>`,
if there is only one app in the package, or `https://public.datagrok.ai/apps/<PACKAGE_NAME>/<APP_NAME>`,
if there are several. You'd learn the right URL first time you run the app from within the application launcher,
so we recommend you to bookmark it for later use.

See also:

  * [Grok JS development](develop.md)
  * [Developing grok applications](develop/develop.md#applications)
  * [Applications on Datagrok Public](https://public.datagrok.ai/apps)
  * [Development samples gallery](https://public.datagrok.ai/js)