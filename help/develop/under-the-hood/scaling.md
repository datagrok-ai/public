---
title: "Scaling"
---

Datagrok was designed to connect millions of users with millions of datasets, data sources, algorithms, scripts, and
applications. It should be able to handle such volumes, and it should do it fast and without cluttering the UI.

Below are some scaling-related approaches that we took. Some of them are related to the infrastructure, others to the
application, and some to the ecosystem we are building.

## Infrastructure

### Database

By keeping metadata separately from data and utilizing scalable storage for both of them (Postgres for metadata, and S3
for data) Datagrok is capable of storing many datasets, both in terms of numbers and volume.

By utilizing our proprietary [in-memory data engine](../under-the-hood/performance.md#in-memory-database)
on both the client and server sides, we can transfer datasets between tiers a lot faster than other systems.

The recommended metadata engine for enterprises is RDS.

### Application server

Our application server is Dart-based and uses asynchronous coding techniques similar to Node.js. It serves a lot of
clients simultaneously
(see [stress testing](../../datagrok/solutions/enterprise/stress-testing-results.md)). Most scientific computations are off-loaded to the
special ["compute" machines](#computations) with autoscaling enabled. Nevertheless, some tasks performed on the app
server are still computationally-intensive (for example, parsing CSV files), so our app server is also multithreaded in
addition to running asynchronous code. It takes advantage of the modern multi-core architecture.

As can be seen from the [stress testing](../../datagrok/solutions/enterprise/stress-testing-results.md), app server scaling won't be needed for most of the
workloads in the enterprise. However, this still could be done if such a need arises.

### Data storage

The recommended data storage engine for enterprises is S3, which is inherently scalable.

### Computations

CVM is used for scientific computations. For high loads, it should be deployed with the auto-scaling mode enabled.

## Application

The web-based application can interactively work with datasets consisting of tens of millions of rows or columns,
entirely on the client, in the browser. To do that, a radical break from the traditional web-based applications was
required.

### Big datasets

In order to work with big datasets right in the browser, we developed our
proprietary [in-memory database](#in-memory-database).

#### In-memory database

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

### Visualizations

While some popular charting libraries like D3 are great and do produce nice-looking results, we are hitting their limits
quite soon once we start working with sizable datasets. First of all, it is not feasible to keep more than 10,000
objects in the browser's DOM tree. Second, since they were not designed for working with big datasets, even the amount
of memory they allocate is sizable (usually, the data needs to be in JSON format, which is suboptimal).

To address that, we had to bite the bullet and implement several
[high-performance visualizations](#viewers)
from scratch. All performance-critical viewers use immediate-mode canvas-based rendering, thus solving the DOM issue.
Moreover, they use our high-performance data engine. The data is not only efficiently stored, but viewers allocate no
additional memory since they work directly with the data engine. Of course, that required a lot of work and meticulous
engineering, but the result is definitely worth it.

#### Viewers

Just like the in-memory database, our [viewers](../../visualize/viewers/viewers.md)
were built from scratch to be able to work with millions of data points at once interactively. All of them make heavy
use of the in-memory database. The ability to use that engine unlocks unique possibilities for viewers, such as easy
data aggregations or passing data to web workers for multithreaded rendering. All viewers access the same data, so no
copies are made, they all share certain statistics and cached calculations, have the same look and feel and usage
patterns, and cooperate on certain tasks.

* Fast, slick, relevant.
* Engineered to take full advantage of DDT
  * Uses DDT's data frames - super-fast and no additional memory overhead
  * Uses the same cached descriptive statistics, sorted orders, etc
  * Many viewers use lightweight, calculated on-the-fly dataframes as an aggregate data source
  * Picking up column metadata (formats, etc.)
  * Fast, extensible, annotated aggregation functions that work across all viewers
* High-performance rendering
  * Choosing the best option for rendering (HTML / canvas / SVG / WebGL) based on the viewer's distinctive features,
    without compromising performance. The stretch goal is for all viewers to be able to visualize a billion rows (
    certain viewers will resort to auto-sampling in order to still be interactive during the data exploration stage).
    Many viewers utilize hybrid rendering systems, i.e., SVG for high-level controls and canvas for performance and
    memory consumption reasons
  * Immediate-mode canvas rendering
  * Renders millions of primitives quickly
    * Adaptive marker rendering - switches between drawing directly on the canvas, rendering from the cache, or
      rendering into an array of bytes. This is transparent to viewers' code.
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

### Extensibility

As the program grows, so does the complexity. Eventually, even a program with the well-thought plugin-based architecture
becomes unwieldy, and users get lost in the hundreds of different menus and thousands of options.

We had to come up with a solution that is not only highly extensible but also has a built-in mechanism for restricting
the feature creep and keeping the UI clean. We addressed these seemingly self-contradictory requirements by designing a
solution that combines first-class support for [functions](../../datagrok/concepts/functions/functions.md)
with the [data augmentation](../../explore/data-augmentation/data-augmentation.md) mechanism.

In the end, we came up with a solution that not only satisfied the initial requirements, but provides a solid foundation
for evolving the ecosystem on a global scale:

* First-class functions, reflectable, with metadata on parameters
* Support for [multiple languages](../../compute/scripting/scripting.mdx)
* Multiple dynamic backends for functions
* Everything is a function (db queries, web services, predictive models, etc)
* Flexible way for [packaging](../../develop/develop.md#packages), deploying, and targeting functions
* Dynamic, asynchronous loading of functions
* Applicable actions are [suggested based on the current context](../../explore/data-augmentation/data-augmentation.md). Out of thousands
  of functions available in the repository, only a handful of relevant ones are being suggested to the user.

See also:

* [Performance](performance.md)
