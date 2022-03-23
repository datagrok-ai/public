<!-- TITLE: JavaScript Development Overview -->
<!-- ORDER: 0 -->

# JavaScript development overview

This document provides a high-level overview of JavaScript development on the Datagrok platform.
Using [our feature-rich JavaScript API], you can develop packages to perform various tasks on Datagrok, namely:

* [Manipulate data]
* [Add views] and [viewers]
* [Develop custom viewers]
* [Register functions]
* [Train and apply predictive models]
* [Build custom apps]

We do not impose any requirements on the UI frameworks or technologies used for the JavaScript plugins. We encourage
you, however, to keep it simple.

## Your workflow

Compared to traditional web applications, the packages you create run _inside the Datagrok platform_
, which makes the development and debugging experience different. Essentially, you develop packages locally, and
Datagrok runs them remotely.

Datagrok provides two options to run the JavaScript code:

* Use the built-in JavaScript editor: On the [Datagrok Public site], go to **Functions** > **
  Scripts** > **New JavaScript Script**. The [Scripting section](../../compute/scripting.md) provides more details on using
  Datagrok's JavaScript editor.
* For reusable functions, viewers, and applications, you can [develop packages].

While the built-in JavaScript editor is useful to quickly try out the Datagrok features, _packages_
are what you will focus on most of the time developing for Datagrok.

The following sections provide more details on packages:

* [Packages], covers the creation and general aspects of building a Datagrok package
* [Package function types], explains what type of functions you typically create in a package
* [Debugging], explains how to configure Visual Studio Code and JetBrains WebStorm to debug packages
* [Publishing], discusses how to publish a package to Datagrok

To streamline your experience, we've established the following workflow:

1. Your package is uploaded to the remote server at startup.
2. The server runs the package.

You can still set breakpoints and debug your package as usual.

You should deploy packages [in development mode] so that only you can see and use this package in Datagrok. This ensures
that multiple people can simultaneously work on the same package.

## Environments

Datagrok comes with two environments:

* [Public][datagrok-production-environment]. Use this environment to publish your package for all Datagrok users.
* [Dev][datagrok-development-environment]. Use this environment to develop a package for Datagrok.

To configure the Datagrok servers that you use, refer to the [Datagrok configuration section](_datagrok-configuration.md). After
you configure the development server, use the `dev` scripts in your `package.json` to run the package on the
**Dev** environment.

### Datagrok inspector

To simplify development, Datagrok provides an Inspector, a tool that helps you peek under the hood of the platform and
see which events get fired and when, how views and viewers are serialized, what widgets are registered by the platform,
and other important aspects of how Datagrok works.

![](./datagrok-inspector.png)

To run the Inspector, on the [Datagrok site](https://dev.datagrok.ai), select **Alt** + **I**.

## Documentation

Here are a few sources that you might want to run through to learn more about Datagrok:

* [JavaScript code samples]
  Datagrok provides an interactive tool for browsing, editing, and running JavaScript samples. Samples are grouped by
  domain, such as data manipulation, visualization, or cheminformatics.
* [JavaScript API] and [reference]
  Use the JavaScript API to create packages for Datagrok. The JavaScript code samples also use the API.
* [Datagrok community]
* [Slack space]

## What's next?

* [Packages](./_packages.md)

[Manipulate data]: https://datagrok.ai/help/develop/js-api#data-manipulation

[Add views]: https://datagrok.ai/help/develop/js-api#views

[viewers]: https://datagrok.ai/help/develop/how-to/manipulate-viewers

[Develop custom viewers]: https://datagrok.ai/help/develop/how-to/develop-custom-viewer

[Register functions]: https://datagrok.ai/help/develop/js-api#registering-functions

[Train and apply predictive models]: https://datagrok.ai/help/learn/predictive-modeling

[Build custom apps]: https://datagrok.ai/help/develop/package-function-types

[datagrok-production-environment]: https://public.datagrok.ai/

[datagrok-development-environment]: https://dev.datagrok.ai/

[JavaScript code samples]: https://public.datagrok.ai/js

[JavaScript API]: https://datagrok.ai/help/develop/js-api

[reference]: https://datagrok.ai/js-api/

[Datagrok community]: https://community.datagrok.ai/

[Slack space]: https://datagrok.slack.com