<!-- TITLE: Parallel computing -->

# Parallel computing in a browser

In this section we bring the best practices of using Web Workers in JavaScript/TypeScript code
of your Datagrok packages. We recommend them as they help provide a more low-latency,
high-performance experience for your end-users.

This article does not require you to have any background in Web Workers or parallel computing.
After reading it, you will be able to use (with a good degree of convenience!) the power of
parallel computation in a browser on Datagrok platform.

Threading is organized in modern web browsers through Web Workers. They are part of ECMAScript
standard, though details may vary from browser to browser. We rely on Google Chrome as a
reference browser.

Though the basic techniques of using Web Workers covered in this article are not specific to
Datagrok, there are important aspects of using Web Workers as part of a webpack-based and/or
TypeScript packages. In addition, we are providing convenience parallel computation libraries
built on top of Web Workers as part of the Datagrok distribution.

## Using a single Web Worker for a heavy task to unblock the UI

When a piece of JavaScript/TypeScript code takes longer than a few hundred milliseconds, it is
worth running it in a separate thread. In this case such computation will not block the main UI
thread. This is pivotal for providing uninterrupted user experience.

We start with a worker itself. Let's create one which computes primes within a given range.
We won't get far into making the computation itself very efficient, and use instead a rather
naive version. The purpose is to bring to the scene a function taking considerable time to run.

```
```

## Using a thread pool to optimize execution of independent tasks

## Parallelizing a single task using multiple Workers

# Parallel computing on a server