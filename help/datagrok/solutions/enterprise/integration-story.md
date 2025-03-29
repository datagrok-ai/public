---
title: "Integration"
---

Datagrok was designed to be as extensible and customizable as possible. The platform exposes a number of integration
points, allowing you to customize it according to your organization's needs, and integrate with existing systems.

In short, if anything can be integrated or automated in some way, chances are we have a mechanism for doing that
already; if not, let us know and we'll come up with a solution!

## Custom scripts

[Scripting](../../../compute/scripting/scripting.mdx) is an integration mechanism with different languages, mostly used for statistical
computing and machine learning. Scripting combines fast interactive visualizations and other features of the Datagrok
platform with thousands of statistical packages and visualizations available in
[R](https://www.r-project.org/about.html), [Python](https://www.python.org)
, [Octave](https://octave.org/),
[Julia](https://julialang.org), or [JavaScript](https://www.javascript.com).

Note that other languages, such as Java, C#, or Node.js can be integrated in a similar way.

## Custom data connectors

Out of the box, the platform comes with the data connectors
for [30+ popular databases](../../../access/databases/connectors/connectors.md), and the list is constantly growing. In addition to
that, it is possible to develop your own data connectors, and seamlessly integrate them into the platform.

In order to do that, a "Grok Connect" REST endpoint that implements a few methods has to be registered with the
platform. The methods are:

* getConnectors - returns all connectors that the endpoint supports (one for database type)
* getSchema(connection) - if applicable, returns database schema for the given connection
* testConnection(connection) - tests the connection
* execute(query) - executes the specified query
* queryTable(structuredQuery) - executes a structured query

At startup, the server asks each registered endpoint for the list of supported connectors, and creates a global list of
supported connectors. The client asks the server for the available connectors, and populates the UI accordingly. Later
on, when client makes a request to query a database, this request gets accepted by a server, and then routed to the
corresponding database connector.

All Datagrok connectors are open-sourced under the MIT license and reside in the
[Datagrok public repository](https://github.com/datagrok-ai/public/tree/master/connectors). A
command-line ["GrokConnectTest"](https://github.com/datagrok-ai/public/tree/master/connectors/grok_connect/src/test/java/grok_connect)
application could be useful for testing and debug purposes.

## Openapi

Datagrok integrates with [OpenAPI](../../../access/open-api.md) really well, automatically mapping
OpenAPI's [paths](https://swagger.io/docs/specification/basic-structure/) to
Grok's [functions](../../concepts/functions/functions.md). This has many benefits:

* Ability to easily call that web method from:
  * [Console](../../navigation/panels/panels.md#console)
  * Event handler
  * [Info pane](../../../datagrok/navigation/panels/info-panels.md)
* Audit
* Data lineage

## JavaScript API

JavaScript is the lingua franca of the web, and naturally it is a first-class language in the Datagrok ecosystem. We
expose a [JavaScript API](../../../develop/packages/js-api.md) that allows you to control most of the platform, including data
manipulation, handling platform events, creating custom viewers, controlling window docking, customizing the platform,
etc.

## Web API

Datagrok's server provides a REST API that lets you programmatically invoke server-side methods.
(OpenAPI is work in progress)

## Grok SDK

A number of command-line utilities for server management. (Work in progress)

## Client-side settings

Client-side settings are specific to the user, and are controlled by the user (
unless the organization's IT policy overrides it).

## Server-side settings

[Server settings](../../../deploy/configuration.md) can only be set at the platform start. Their change requires server restart.
