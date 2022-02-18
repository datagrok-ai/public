<!-- TITLE: Architecture -->
<!-- SUBTITLE: -->

# Architecture

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

A traditional approach to building such a system would be to commit to widely used technology such as JavaScript, select
a bunch of widely used libraries and frameworks, including a UI framework, a visualization library, a serialization
library, and a dozen of other ones, and plumb it all together. While the result would be somewhat predictable, it will
also inevitably be mediocre due to all the compromises this path implies. You cannot build a fighter jet with a steam
engine, random Lego parts, and a lot of duct tape.

Our firm belief, backed by 20 years of experience in developing analytical solutions and data analysis platforms, is
that we need a clean break for such an ambitious project. Before building high-level features, we need to dive as deep
as possible, build as an efficient foundation as possible, own every single detail of the implementation, and only then
build our way up.

### Data engine

First of all, we built a fit-for-purpose [data engine](infrastructure.md#in-memory-database), which is essentially an
in-memory columnar database. It is the heart of the system and is highly optimized for the routines typical for
exploratory data analysis, i.e., fast sequential access for column data, efficient data storage, data transfer, data
compression, filtering, aggregations, joins, caches, descriptive statistics, etc. It has an expressive and clean API.

### Viewers

On top of the data engine, we built several [high-performance viewers](infrastructure.md#viewers)
that make heavy use of the data engine. The ability to use that engine unlocks unique possibilities for viewers, such as
easy data aggregations or passing data to web workers for multithreaded rendering. All viewers access the same data (no
copies are made), share specific statistics and cached calculations, have the same look and feel and usage patterns, and
cooperate in particular tasks.

### Server

Also, on top of the data engine, we built a [server](infrastructure.md#datlas). It is used for many purposes: data
retrieval from different data sources, security, user management, metadata storage, multiple object repositories, user
collaboration, running arbitrary scripts written in different languages, building and applying predictive models, and a
lot more. The server utilizes extensible plugin architecture and currently has over 20 plugins. The server enables the
client to do what is needed most efficiently, minimizing memory consumption and network traffic. It helps that both
client and server are written in Dart and share a lot of common libraries.

### Application client

Next, we built a [data analytics application client](infrastructure.md#web-application) that integrates our data engine,
viewers, and server-based capabilities. The application can run entirely autonomously in the browser without making
requests to the server. In that mode, it lets users import data of non-trivial sizes (tens of millions of rows) from
local files in many popular data formats, clean and transform it using some tools and routines, interactively visualize
it, and export the data. Using the server-based capabilities, it gets to another level by centrally managing data
connections, enabling data governance and provenance, enabling scripting, providing data modeling capabilities, and
dozens of other features outlined in our [knowledge base](../../home.md).

In the results, Datagrok is an application that enables users to ingest, transform, analyze, visualize, model, and share
data by combining powerful client and server capabilities. The application is highly flexible and extensible and allows
creating new data pipelines, models dynamically and building viewers, plugins, and applications on top of it.

### Compute Engine

To extend already massive functionality, we created a compute engine for performing on-server computations. Compute
engine lets the user easily [invoke scripts](../../compute/scripting.md) written in languages such as R or Python and
expose them to the Datagrok ecosystem as functions. Also, it utilizes a few different backends
for [predictive modeling](../../learn/predictive-modeling.md) while providing users the same look and feel. To create
and share documents that contain live code, equations, visualizations, and narrative text, we added
[Jupyter Notebook](../../compute/jupyter-notebook.md) features to compute server. Datagrok allows efficiently creating,
editing, importing, linking, and applying Notebooks into tables. To scale computations, spin out multiple instances of
some component, and our load balancer will take care of dispatching computations.

## Useful links

* [Enterprise evaluation FAQ](enterprise-evaluation-faq.md)
* [Infrastructure](infrastructure.md) with detailed description of every component
* [Deployment](deploy.md) instruction to install and try Datagrok
