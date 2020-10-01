<!-- TITLE: Scaling -->
<!-- SUBTITLE: -->

# Scaling

Datagrok was designed with the goal to connect millions of users with millions
of datasets, data sources, algorithms, scripts, and applications. Not only it should be
able to handle such volumes, but it should do it fast, and without cluttering the UI.

Below are some of the scaling-related approaches that we took.
Some of them are related to the infrastructure, others to the application, 
and some are related to the ecosystem we are building. 

## Infrastructure

### Database

By keeping metadata separately from data, and utilizing scalable storage for both of 
them (Postgres for metadata, and S3 for data) Datagrok is capable of storing many
datasets, both in terms of numbers and volume.

By utilizing our proprietary [in-memory data engine](../performance.md#in-memory-database)
on both client and server side, we can transfer datasets between tiers a lot faster
than other systems.

The recommended metadata engine for enterprises is RDS. 

### Application Server

Our application server is written in Dart, and uses asynchronous coding techniques
similar to Node.js. It serves a lot of clients simultaneously 
(see [stress testing](stress-testing-results.md)). Most scientific
computations are off-loaded to the special ["compute" machines](#computations) with 
autoscaling enabled. Nevertheless, some tasks performed on the app server are still
computationally-intensive (example is parsing CSV files), so in addition to 
running asynchronous code our app server is also multithreaded, and takes advantage
of the modern multi-core architecture.

As can be seen from the [stress testing](stress-testing-results.md), app server scaling 
won't be needed for most of the workloads in the enterprise. However, this still could
be done if such a need arises.    

### Data Storage

The recommended data storage engine for enterprises is S3, which is inherently
scalable. 

### Computations

CVM is used for scientific computations. For high loads, it should be deployed with
the auto-scaling mode enabled. 
See [Compute Virtual Machine](architecture.md#compute-virtual-machine) for details.

## Application

The web-based application  is capable of interactively working with datasets consisting
of tens of millions of rows, or millions of columns, completely on the client, 
in the browser. In order to do that, a radical break from the traditional web-based
applications was required.

### Big Datasets

In order to work with big datasets right in the browser, we developed our
proprietary [in-memory database](architecture.md#in-memory-database).

### Visualizations

While some popular charting libraries like D3 are great and do produce nice-looking results, 
we are hitting their limits quite soon once we start working with sizable datasets. First of all,
it is not feasible to keep more than 10,000 objects in the browser's DOM tree. Second,
since they were not designed for working with big datasets, even the amount of memory they
allocate is sizable (usually, the data needs to be in some sort of JSON format, which is 
suboptimal).

To address that, we had to bite the bullet and implement a number of 
[high-performance visualizations](architecture.md#viewers) 
from scratch. All performance-critical viewers use immediate-mode canvas-based
rendering, thus solving the DOM issue. Moreover, they use
our high-performance data engine. Not only the data is efficiently stored,
but viewers allocate no additional memory since they work directly with the data engine.
Of course, all of that required a lot of work and very careful engineering, but the
result is definitely worth it.   

### Extensibility

As the program grows, so does the complexity. Eventually, even a program with the 
well-thought plugin-based architecture becomes unwieldy, and users get lost in the 
hundreds of different menus and thousands of options.

We had to come up with a solution that is not only highly extensible, but also
has a built-in mechanism for restricting the feature creep, and keeping the UI clean.
We addressed these seemingly self-contradictory requirements by designing a solution 
that combines first-class support for [functions](../../overview/functions/function.md)
with the [data augmentation](../../discover/data-augmentation.md) mechanism.

In the end, we came up with a solution that not only satisfied the initial requirements,
but provides a solid foundation for evolving the ecosystem on a global scale:

* First-class functions, reflectable, with metadata on parameters
* Support for [multiple languages](../../compute/scripting.md)
* Multiple dynamic backends for functions
* Everything is a function (db queries, web services, predictive models, etc) 
* Flexible way for [packaging](../../develop/develop.md#packages), deploying, and targeting functions
* Dynamic, asynchronous loading of functions
* Applicable actions are [suggested based on the current context](../../discover/data-augmentation.md). 
  Out of thousands functions available in the repository, only a handful relevant ones are
  being suggested to the user.   

See also:

  * [Performance](../performance.md) 
