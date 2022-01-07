<!-- TITLE: Architecture Design -->
<!-- SUBTITLE: -->

# Architecture Design

## Goals

Datagrok was designed with the following main objectives:

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
  changes, and reason about it. This is only possible when working with a clean, minimalistic, and coherent codebase.

These goals heavily influenced the design. On the surface, the platform appears quite traditional, as it consists of two
parts - client and server. However, the way these parts work and connect is quite unusual.

### Design philosophy

A traditional approach to building such a system would be to commit to a widely used technology such as JavaScript,
select a bunch of widely used libraries and frameworks, including a UI framework, a visualization library, a
serialization library, and a dozen of other ones, and plumb it all together. While the result would be somewhat
predictable, it will also inevitably be mediocre due to all the compromises this path implies. You can’t build a fighter
jet with a steam engine, random Lego parts, and a lot of duct tape.

Our firm belief, backed by 20 years of experience in developing analytical solutions and data analysis platforms, that
for such an ambitious project, we need a clean break. Before building high-level features, we need to dive as deep as
possible, build as efficient foundation as possible, own every single detail of the implementation, and only then build
our way up.

First of all, we built a fit-for-purpose [data engine](#in-memory-database), which is essentially an in-memory columnar
database. It is the heart of the system and is highly optimized for the routines typical for exploratory data analysis,
i.e., fast sequential access for column data, efficient data storage, data transfer, data compression, filtering,
aggregations, joins, caches, descriptive statistics, etc. It has an expressive and clean API.

On top of the data engine, we built several
[high-performance viewers](#viewers)
that make heavy use of the data engine. The ability to use that engine unlocks unique possibilities for viewers, such as
easy data aggregations or passing data to web workers for multithreaded rendering. All viewers access the same data (no
copies are made), share specific statistics and cached calculations, have the same look and feel and usage patterns, and
cooperate in particular tasks.

Also, on top of the data engine, we built a [server](architecture-details.md#datlas). It is used for many purposes: data
retrieval from different data sources, security, user management, metadata storage, multiple object repositories, user
collaboration, running arbitrary scripts written in different languages, building and applying predictive models, and a
lot more. The server utilizes extensible plugin architecture and currently has over 20 plugins. The server enables the
client to do what is needed most efficiently, minimizing memory consumption and network traffic. It helps that both
client and server are written in Dart and share a lot of common libraries.

Next, we built a [data analytics application](#application) that puts together our data engine, viewers, and
server-based capabilities. The application is capable of running entirely autonomously in the browser without making
requests to the server. In that mode, it lets users import data of non-trivial sizes (tens of millions of rows) from
local files in many popular data formats, clean and transform it using some tools and routines, interactively visualize
it, and export the data. By using the server-based capabilities, it gets to another level by centrally managing data
connections, enabling data governance and data provenance, enabling scripting, providing data modeling capabilities, and
dozens of other features outlined in our [knowledge base](../../home.md).

## In-memory database

At the heart of the platform is the proprietary, unique technology that lets us efficiently work with massive datasets
in the browser. Essentially, it is a columnar in-memory database engineered from scratch and optimized for exploratory
data analysis, interactive visualizations, and machine learning.

* Efficient in-memory storage
    * Column-based data store
    * High-performance bitsets
    * Versionable tables, columns, and bitsets
    * Adaptive bit storage for ints
    * Out-of-box support for big integers
    * All strings are categories
    * Backed by typed arrays
        * Cache locality
        * SIMD instructions
        * Lightweight cloning (which enables multithreading via isolates)
* Custom binary serialization with compression - fast persistence
    * Adaptive bit storage
    * Automatic pattern recognition
    * RLE, Huffman
    * FCP encoding for floating-point data
    * Designed to be extended with the domain-specific compression algorithms
* Support for user-defined types and extending existing types
    * Parsing
    * Formatting
    * Aggregation
    * Comparison and sorting
    * Persistence
* Built-in intelligent CSV parser
    * Automatically handles delimiters and comments
    * Adaptive reading: makes one speedy pass to find out some features about the file (whether it has quotes,
      multi-lines, etc.) and chooses the strategy accordingly
    * Adjusts strategy for performance/memory as it reads
    * Smart parsing of dates
    * Multithreaded parsing
        * Does not block the main UI thread
        * 10x faster compared to the competition (coming soon)
    * Preview (isolates)
* Flexible sorting
    * Simple API for sorting by simple/multiple columns
    * Default natural sorting for strings (“study2” comes before “study10”)
    * Built-in custom category sorting (ex: Mon, Tue, Wed, Thurs, Fri, Sat, Sun)
    * Ability to pass custom comparison function that will also be used by the rest of the engine (
      sorting, aggregations, grouping, etc.)
* Built-in high-performance descriptive statistics
    * Counts, min, max, sum, avg, med, avg, stdev, q1, q2, q3
    * All stats calculated in one pass where possible
    * Auto-cacheable (client code should not worry about calculating it twice)
* Joins
    * Joined table have metadata that helps to link it back to the source tables easily
* Aggregations
    * Uses the same set of high-perf statistical routines
    * Fluent API:
    * t.groupBy(["race", "sex"]).avg("height").count("subjId").aggregate();
    * Pivoting
    * Supports user-defined aggregations
    * Free-text, SQL-like queries
* Metadata on column and data frame levels (units, quality, auto-formatting, etc.)
* Change notifications
    * A custom eventing mechanism used across the whole platform allows for easy listening to, aggregation, filtering,
      routing, and logging events.

## Viewers

Just like the in-memory database, our [viewers](../../visualize/viewers.md)
were built from scratch to be able to work with millions of data points at once interactively. All of them make heavy
use of the in-memory database. The ability to use that engine unlocks unique possibilities for viewers, such as easy
data aggregations or passing data to web workers for multithreaded rendering. All viewers access the same data, so no
copies are made, they all share certain statistics and cached calculations, have the same look and feel and usage
patterns, and cooperate on certain tasks.

* Fast, slick, relevant.
* Engineered to take full advantage of DDT
    * Uses DDT’s data frames - super-fast and no additional memory overhead
    * Uses the same cached descriptive statistics, sorted orders, etc
    * Many viewers use lightweight, calculated on-the-fly dataframes as an aggregate data source
    * Picking up column metadata (formats, etc.)
    * Fast, extensible, annotated aggregation functions that work across all viewers
* High-performance rendering
    * Choosing the best option for rendering (HTML / canvas / SVG / WebGL) based on the viewer’s distinctive features,
      without compromising performance. The stretch goal is for all viewers to be able to visualize a billion rows (
      certain viewers will resort to auto-sampling in order to still be interactive during the data exploration stage).
      Many viewers utilize hybrid rendering systems, i.e., SVG for high-level controls and canvas for performance and
      memory consumption reasons
    * Immediate-mode canvas rendering
    * Renders millions of primitives quickly
        * Adaptive marker rendering - switches between drawing directly on the canvas, rendering from the cache, or
          rendering into an array of bytes. This is transparent to viewers’ code.
    * Multithreaded rendering
    * WebGL-accelerated rendering with custom shaders (coming soon)
    * Adaptive rendering behavior - the system keeps track of how long it took each viewer to render and optimizes
      accordingly - for instance, “fast” viewers are rendered first, and “slow” viewers are not re-rendered while a
      slider is being dragged.
* Interactivity and synchronization
    * Current row, mouse-over row, current column, mouse-over column, mouse-over row group
* Viewers as first-class citizens
    * Register, query, instantiate, attach to a data source, add to view, use as a tooltip for row groups, render
      viewers dynamically. Usage example: a full-screen mode that applies to all viewers.
    * Viewer descriptors: name, tooltip, best size/position, type of accepted data
* Properties infrastructure
    * Persistence
    * Consistent names enforced by conventions
    * Discoverability
    * Easy UI bindings (property grid, menu, dialogs)
    * out-of-box automatic validation
    * Minimum overhead (convention over configuration)
    * Automatic, seamless code generation at build time (small code size)
    * Change notification
    * Support for categories
    * Applies to both visual and non-visual objects (such as tables or columns)
    * The same infrastructure is used for editing, property notification, and serialization
    * Minimizing serialized side by comparing to the default value
    * Automatic JSON and binary serialization
* Clean, functional decomposition
    * Clean separation between settings and viewer fields
    * Each viewer consists of three main classes - core, look, and meta. Each of them can be accessed dynamically, which
      allows operating on categories of entities. This dramatically increases code reuse and allows for complex
      customizations to be implemented very easily.
* Share common base, utilize the same tricks, same naming conventions
    * Dense, straight-to-business, and easy-to-understand code
    * Data-aware axes across all viewers (adaptive resolution for time series)
    * One-liners for zooming/selecting (including selection rectangle)
    * Generalized support for viewer extensions (micro-plugins)
* A number of convenience helpers (Rect, Color)
* Harmonized UX across all viewers. Hand-crafted widgets with clean, expressive API
    * Tooltips
    * Notifications
    * Popup menus
    * Main menus
    * In-viewer interactive help system
    * Progress indicators
    * Color pickers
    * Sliders
    * Column selectors (lightweight, searchable, sortable, draggable, extensible, customizable)
* Drag-n-drop
    * Drag files right into the browser (import)
    * Drag columns
        * from anywhere (column manager, grid, column selectors)
        * to anywhere (dialogs, column selectors, viewer settings)
        * very slick and interactive
    * Drag rows (extract rows, etc.)
    * Drag any objects (users, scripts, tables, statistics, viewers, etc.)
* Composable rendering
    * Many controls, such as histogram, are capable of rendering themselves on a canvas that another viewer owns. That
      allows for a lightweight, memory-efficient rendering of complex scenes (ex:
      histograms on a line chart)
* Event bus for common viewer events for decoupling and easier event handling
* Filters
    * Collaborative - each filter has a say in determining whether the row passes the filter
    * Interpretable - a string description of what is filtered in or out
* Easy orchestration - making viewers work together as a team
    * Complementing each other
    * Passing information between viewers (current, mouse-over record)
    * Filter-viewer relationship
* Dialogs
    * Fluent API for easy programmatic construction of dialogs
    * Standard validation across all UI
    * Out-of-box, opt-in persistence
    * Out-of-box, opt-in logging
* Standard (but extensible) editors for data types used across dialogs/grids
* Flex tooltips
    * Row tooltips: Ability to select columns to show on a tooltip (including row viewers)
    * Row group tooltip: Use any viewer as a tooltip (for example, when the mouse is over a particular histogram bin,
      the tooltip contains a scatterplot with values that fall in that bin)
    * Object tooltip (users, etc.)

## Application

We built an application that enables users to ingest, transform, analyze, visualize, model, and share data by combining
powerful client and server capabilities. The application is highly flexible and extensible and allows creating new data
pipelines, models dynamically and build viewers, plugins, and applications on top of it.

![](architecture-diagram1.png)

Datagrok installation consists of two virtual machines:

* [Datagrok Virtual Machine](#datagrok-virtual-machine)
* [Compute Virtual Machine](#compute-virtual-machine)

Also, it needs [database](#database) and [persistent file storage](architecture-details.md#storage). Both of the virtual
machines can be deployed as Docker containers in [AWS EC2](deploy-amazon-ec2.md)
,
[AWS ESC](deploy-amazon-ecs.md) or [regular host machine](deploy-regular.md).

See [Architecture Details](architecture-details.md) for more details.

## Datagrok virtual machine

This machine is the heart of the platform and is required for all activities.

See [Datagrok virtual machine](architecture-details.md#datagrok-virtual-machine) for more details.

### Datlas

Grok server (also referred to as “Datlas”) is a Dart stand-alone application that creates REST endpoints that the client
consumes. The same API can also be used by other parties (for instance, by the IT department).

* Implemented in Dart, shares much code with the client
* Exposes REST API
* Extendable via server plugins
* Exposes some services
* Works on either Windows or Linux
* Built-in support for server configuration and deployment schemas
* Secured data transfer system, based on asymmetric keys
* Customizable storage
* Supports SSL

#### ORM

To efficiently work with the database, we have built a custom object-relational mapper library (ORM)
. It takes advantage of our coding conventions, is fit and tuned for our goals, and lets us do the following:

* Rapidly develop new classes and bind them to database tables
* Use either SQL or API for working with DB entities
* Avoid all boilerplate code for CRUD operations
* Use the same models on both server-side and client-side
* Work with property metadata and use it to optimize queries in run-time
* Free-text filtering (SQL-like language)
* ORM can extract object fields from a database without mapping them to objects

### Full-text search

Datagrok’s full-text search (press Alt+Q to search) is backed by the ElasticSearch. It searches in Wiki, forums, and
datasets.

## Compute virtual machine

Compute Virtual Machine is used for performing on-server computations. It is used for scripting, training and applying
predictive models, and cheminformatics. You might not need it if your use cases do not involve the below-mentioned
tasks.

To scale computations, you might want to spin out multiple compute instances, and our load balancer will take care of
dispatching computations.

See [Compute virtual machine](architecture-details.md#compute-virtual-machine) for more details.

### Scripting

[Scripting](../../compute/scripting.md) lets you easily invoke scripts written in languages such as R or Python, and
expose them to the Datagrok ecosystem as functions.

Scripting works by sending code and data to the Jupyter kernel gateway, one of the ComputeVM components.

Available languages are: Python, R, JS, Octave, Julia.

### Predictive modeling

Datagrok utilizes a few different backends for [predictive modeling](../../learn/predictive-modeling.md), while
providing users the same look and feel. Training and applying models for either of the backends are performed on the
Compute Virtual Machine.

* [H20](https://h2o.ai)
* [Caret](http://topepo.github.io/caret/index.html): Classification And Regression Training
* [ChemProp](https://github.com/swansonk14/chemprop): Chemical Property Prediction with Graph Convolutional Networks

### Jupyter Notebook

[Jupyter Notebook](../../compute/jupyter-notebook.md) is an open-source web application that allows you to create and
share documents that contain live code, equations, visualizations, and narrative text. Datagrok allows creating, edit,
import, link and apply Notebooks into tables.

## Database

Metadata associated with users, datasets, algorithms, predictive models, etc. is kept in a Postgres database. Having the
data stored in a relational database brings well-defined relations and data consistency, enables us to work efficiently
with complex queries, and adds transactional support.

Postgres is free and can easily be deployed locally for development. Also, Postgres data protocol is prevalent and is
used in many big data solutions, as well as in cloud databases. If necessary, we can switch to a scalable solution like
Aurora or CocroachDB without changing much code.

See [Database](architecture-details.md#database) for more details.

## Deployment

Datagrok is a web application, which means no deployment efforts per user once the server is set up. Administration
tasks could be performed via the web interface as well.

Enterprises typically prefer on-premise deployment for multiple reasons, such as security, ability to easily access
internal data, and other features such as integration with the enterprise
[authentication](../../govern/authentication.md) methods. One of the most convenient ways of doing that is deploying the
platform on the AWS in the Virtual Private Cloud. It is as easy as spinning out Amazon’s EC2 machines with our docker
images inside your virtual private network. See [Datagrok deploy on AWS EC2](deploy-amazon-ec2.md),
[Datagrok deploy on AWS ECS](deploy-amazon-ecs.md),
[Running Datagrok docker container](deploy-regular.md) for details.

See [Enterprise evaluation FAQ](enterprise-evaluation-faq.md) for more details.
