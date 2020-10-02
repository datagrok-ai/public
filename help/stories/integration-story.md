<!-- TITLE: Integration -->
<!-- SUBTITLE: -->

# Overview

Datagrok was designed to be as extensible and customizable as possible. The platform exposes a number of 
integration points, allowing you to customize it according to your organization's needs, and integrate
with existing systems.

In short, if anything can be integrated or automated in some way, chances are we have a mechanism 
for doing that already; if not, let us know and we'll come up with a solution! 

## Custom Scripts

[Scripting](../develop/scripting.md) is an integration mechanism with different languages,
mostly used for statistical computing and machine learning. 
Scripting combines fast interactive visualizations and other features of the Datagrok platform 
with thousands of statistical packages and visualizations available in 
[R](https://www.r-project.org/about.html), [Python](https://www.python.org),
[Julia](https://julialang.org), or [JavaScript](https://www.javascript.com).

Note that other languages, such as Java, C#, or Node.js can be integrated in a similar way.

## Custom Data Connectors

Out-of-the box, the platform comes with the data connectors  
for [30+ popular databases](../access/data-connection.md#connectors), and the list is constantly growing.
In addition to that, it is possible to develop your own data connectors, and seamlessly integrate them
into the platform.

In order to do that, a "Grok Connect" REST endpoint that implements a few methods has to be registered with the platform. 
The methods are:
* getConnectors - returns all connectors that the endpoint supports (one for database type)
* getSchema(connection) - if applicable, returns database schema for the given connection
* testConnection(connection) - tests the connection
* execute(query) - executes the specified query
* queryTable(structuredQuery) - executes a structured query

At startup, the server asks each registered endpoint for the list of supported connectors, and 
creates a global list of supported connectors. The client asks the server for the available connectors,
and populates the UI accordingly. Later on, when client makes a request to query a database, this request
gets accepted by a server, and then routed to the corresponding database connector.

All Datagrok connectors are open-sourced under the MIT license and reside in the 
[Datagrok public repository](https://github.com/datagrok-ai/public/tree/master/connectors).
A command-line ["GrokConnectTest"](https://github.com/datagrok-ai/public/tree/master/connectors/grok_connect/src/test/java/grok_connect) 
application could be useful for testing and debug purposes. 

## OpenAPI

Datagrok integrates with [OpenAPI](../access/open-api.md) really well, automatically mapping
OpenAPI's [paths](https://swagger.io/docs/specification/basic-structure/) to Grok's [functions](../overview/functions/function.md).
This has many benefits:
* Ability to easily call that web method from:
  * [Console](../overview/console.md)
  * Event handler
  * [Info panel](../discover/info-panels.md)
* Audit
* Data lineage

## Grok JavaScript API

JavaScript is the lingua franca of the web, and naturally it is a first-class language in the Grok ecosystem. 
We expose a [JavaScript API](../develop/js-api.md) that allows you to control most of the platform, including
data manipulation, handling platform events, creating custom viewers, controlling window docking,
customizing the platform, etc.  

## Grok Web API

Datagrok's server provides a REST API that lets you programmatically invoke server-side methods. 
(OpenAPI is work in progress)

## Grok SDK

A number of command-line utilities for server management. (Work in progress)

## Client-side settings

[Client-side settings](../overview/settings.md) are specific to the user, and are controlled by the user (unless 
organizations's IT policy overrides it).

## Server-side settings

[Server settings](../overview/settings-server.md) are controlled by platform's [administrators](../govern/security.md).
