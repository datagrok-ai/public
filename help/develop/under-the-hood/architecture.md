---
title: "Architecture"
---

<!--Incorrect links
This group of pages describes the internal mechanics of Datagrok. It is meant for more [advanced] users who would like to better understand how we [do X] to achieve excellent performance and provide features such as [1],[2],[3],[4],[5], and [6].
* Performance
* Architecture
* Grok script-->

## Goals

Datagrok is designed with the following main objectives:

* **Performance**. This cannot be overstated. A modern system should be capable of working with tens of millions of
  datapoints, and every bit of performance counts. Memory consumption, latency, and network traffic fall in the same
  category.
* **Security**. We should be in complete control of the data.
* **Usability**. A system should be a joy to work with.
* **Modularity**. We should safely turn different features on or off without affecting the rest of the platform.
* **Extensibility**. We should be able to add new functionality easily.
* **Interoperability**. Other systems should be able to talk to us via open protocols, as well as seamlessly integrate
  with other platforms, environments, and programming languages.
* **Maintainability**. In the long run, it is indispensable to have the ability to navigate in the code easily, make
  changes, and reason about it. This is only possible when working with a clean, minimalistic, coherent codebase.

These goals heavily influenced the design. The platform appears quite traditional on the surface, as it consists of two
parts - client and server. However, the way these parts work and connect is quite unusual.

## Design philosophy

A traditional approach would be to pick popular libraries — a UI framework, a
visualization library, a serialization library — and integrate them. While
predictable, this approach produces mediocre results because of the inherent
compromises. You cannot build a fighter jet from random Lego parts and duct tape.

Our firm belief, backed by 20 years of experience in developing analytical solutions and data analysis platforms, is
that we need a clean break for such an ambitious project. Before building high-level features, we need to dive as deep
as possible, build as an efficient foundation as possible, own every single detail of the implementation, and only then
build our way up.

### Data engine

First of all, we built a fit-for-purpose [data engine](scaling.md#in-memory-database), which is essentially an
in-memory columnar database. It is the heart of the system and is highly optimized for the routines typical for
exploratory data analysis, i.e., fast sequential access for column data, efficient data storage, data transfer, data
compression, filtering, aggregations, joins, caches, descriptive statistics, etc. It has an expressive and clean API.

### Viewers

On top of the data engine, we built several [high-performance viewers](scaling.md#viewers)
that make heavy use of the data engine. The ability to use that engine unlocks unique possibilities for viewers, such as
easy data aggregations or passing data to web workers for multithreaded rendering. All viewers access the same data (no
copies are made), share specific statistics and cached calculations, have the same look and feel and usage patterns, and
cooperate in particular tasks.

### Server

Also, on top of the data engine, we built a [server](infrastructure.md#1-core-components). It is used for many purposes: data
retrieval from different data sources, security, user management, metadata storage, multiple object repositories, user
collaboration, running arbitrary scripts written in different languages, building and applying predictive models, and a
lot more. The server utilizes extensible plugin architecture and currently has over 20 plugins. The server enables the
client to do what is needed most efficiently, minimizing memory consumption and network traffic. It helps that both
client and server are written in Dart and share a lot of common libraries.

### Application client

Next, we built a [data analytics application client](infrastructure.md#1-core-components) that integrates our data engine,
viewers, and server-based capabilities. The application can run entirely autonomously in the browser without making
requests to the server. In that mode, it lets users import data of non-trivial sizes (tens of millions of rows) from
local files in many popular data formats, clean and transform it using some tools and routines, interactively visualize
it, and export the data. Using the server-based capabilities, it gets to another level by centrally managing data
connections, enabling data governance and provenance, enabling scripting, providing data modeling capabilities, and
dozens of other features outlined in our [knowledge base](../../datagrok/datagrok.md).

In the results, Datagrok is an application that enables users to ingest, transform, analyze, visualize, model, and share
data by combining powerful client and server capabilities. The application is highly flexible and extensible and allows
creating new data pipelines, models dynamically and building viewers, plugins, and applications on top of it.

### Compute Engine

The compute engine runs server-side computations. You can
[invoke scripts](../../compute/scripting/scripting.mdx) written in R, Python, and
other languages and expose them to the Datagrok ecosystem as functions. It supports
multiple backends for [predictive modeling](../../learn/learn.md) with a consistent
interface. [Jupyter Notebook](../../compute/jupyter-notebook.md) integration lets
you create and share documents with live code, equations, visualizations, and
narrative text. To scale, spin up multiple instances of a component — the load
balancer dispatches computations automatically.

## Useful links

* [Enterprise evaluation FAQ](../../datagrok/solutions/teams/it/enterprise-evaluation-faq.md)
* [Infrastructure](infrastructure.md) with detailed description of every component
* [Deployment](../../deploy/deploy.md) instruction to install and try Datagrok
