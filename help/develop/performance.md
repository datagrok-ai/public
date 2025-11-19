---
title: "Performance"
sidebar_position: 0
---

# Performance

At Datagrok, we don't write slow software. Efficient programming is a craft that needs 
to be developed. Nothing beats deliberate practice, but here are some of the approaches you
might find useful:

## General recommendations

Some general recommendations on writing high-performance JavaScript code:

* [Writing Fast, Memory-Efficient JavaScript](https://www.smashingmagazine.com/2012/11/writing-fast-memory-efficient-javascript/)
* [JavaScript Performance](https://developer.mozilla.org/en-US/docs/Learn/Performance/javascript_performance)
* [Chrome DevTools Documentation](https://developer.chrome.com/docs/devtools/)
* [Web development best practices](https://web.dev/fast/)

Below, we will discuss some performance-related topics important to building solutions with Datagrok. However, nothing
beats common sense, benchmarks, and eventually developing an intuition of how fast or slow a particular method would
work, and why. Here are some universal recommendations:

* Know how long each operation (network call, memory access, etc) should take
* Master the tooling (Chrome Profiler, Network tab, etc)
* Maintain a suite of benchmarks
* Look into problems, and check your assumptions using the Chrome Profiler

## Dataframe

**DO NOT** use row-based access for iterating over rows when [performance](./advanced/performance-tips.md)
matters (pretty much anytime when the size of the dataset is not known in advance). Each call to `row(i)` creates
a `Row` object that is unnecessary, since the underlying storage is columnar. Only use it for passing a reference to a
particular row. Prefer using `column.get(i)` methods, instead.

## Iterables and arrays

**PREFER** using [typed arrays](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays)
instead of the standard JS lists when working with numbers or bitsets, if possible.

When working with iterables, **DO NOT** create arrays for the sole purpose of iterating over their elements, this is
wasteful. Consider using [wu](https://github.com/fitzgen/wu.js/) instead, it is already included with Datagrok. For
instance, instead of using `Array.from(values).length` use `wu(values).length`.

**DO NOT** use temporary arrays for implementing a generic functionality as a one-liner. Instead, consider developing or
finding a utility function for that. The following code is wasteful (unnecessary allocation)
, slow (use of lambdas) and inefficient (likely a typed array would perform better).
`new Array(nItems).fill(0).map((_, i) => begin + i)`. A better way would be a specialized utility function if a real
array is needed, or [wu.count](https://fitzgen.github.io/wu.js/#count) if you only need to iterate over indexes.

## Common root causes of performance problems

* **Using the suboptimal algorithm**. Think whether a different approach would
  be faster.
* **Doing unnecessary computations upfront**. See if the result
  is needed right now, or perhaps it could be computed later. Example: column tooltips
  are calculated dynamically right when a user needs to see it, but not earlier.
* **Computing what already has been computed**. See if the input has changed,
  perhaps there is no need for recalculation.
* **Recalculating too often**. In case of applications reacting to the streams
  of events, consider using event debouncing to allow multiple events to
  be fired, and then recalculating only once after that.
* **Using DataFrame's API for number crunching**. For maximum performance, consider
  working with raw data instead.
* **Inefficiently using memory**. For number crunching, one of the most important
  things that influences the performance is cache locality. Try to arrange data
  (usually in raw memory buffers) in such a way that your algorithm would access
  data sequentially. Minimize the memory footprint.
* **Using the wrong containers for the job**. Typical error is creating a list
  of objects for the sole purpose to find out whether another object is in the
  list. Consider using Map.
* **Creating too many objects**. Each object comes at a cost - this includes
  allocation, memory consumption, and garbage collection. Think whether you
  need it, especially within inner loops, or as part of a commonly used structure.
* **Chatty client-server interactions**. Calling a web service is expensive, and
  each call introduce additional time penalty. Consider minimizing the number of
  calls, and the amount of data transferred. Use `Network` tab in the Chrome Dev Tools
  to see what's happening.
* **Not caching results of data queries**. If the underlying data has a known
  change cadence (for instance, it's not changing at all, or ETL is run overnight),
  consider using Datagrok's built-in query caching mechanism.

## Visualizations

Here are some approaches we can use to make our [viewers](../visualize/viewers/viewers.md) fast: 

* **Raw performance improvements**.
  Use BitSet efficiently, switch to raw arrays, minimize the number of arithmetic
  operations, do micro optimizations, get rid of lambda functions, extract conditional logic outside loops, use “for”
  loop instead of iterators, etc. Profiler is our best friend here :)
* **Caching**.
  We have versionable columns that let us cache results (such as min max histogram etc) and reuse between viewers. This
  starts to matter on big datasets. LruMap is also convenient.
* **Render only what's visible**.
  Render only after the user scrolls to it. Not applicable that often but was very useful for filter group
* **Adaptive debouncing**
  on big datasets (not doing computations immediately on user input such as dragging the slider).
* **Adaptive settings**.
  Example: turning off distributions on filters when we have more than 10M rows. Or turn off
  “zoom to filter” in scatter plots. You can still turn them on.
* **Adaptive initial data configuration**. 
  If we have a huge dataset, it makes no sense no show all 5 numerical columns
  in a line chart, one is enough - a 5x speedup right here
* **Adaptive interactivity**. 
  As the data size grows, disable mouse-over effects
* **Off-loading computations to workers**.
  Currently we do it in the core for CSV importing already so we can do it for
  viewers as well (remains to be seen if we want to complicate the code)
* **WebGPU rendering and computations**. 
  Super fast but hard to implement, does not work everywhere, and not trivial to
  integrate. Perhaps could be used as a drop-in replacement for existing functions.