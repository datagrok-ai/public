# Datagrok package repository

This is a public repository for the API, tools,
and [packages](https://datagrok.ai/help/develop/develop#packages)
available for [Datagrok™](https://datagrok.ai), a next-generation web-based data analytics platform.
The platform is very extensible, and almost anything could be implemented as a package:

* Support for scientific domains, such as [cheminformatics](packages/Chem/README.md)
* Applications, such as [Clinical Case](packages/ClinicalCase/README.md)
  or [Peptides](packages/Peptides/README.md)
* Connectors
  to [OpenAPI web services](https://github.com/datagrok-ai/public/tree/master/packages/Swaggers)
* Visualizations, such as [Leaflet](packages/Leaflet/README.md)
* Importing and previewing files, such as
  [SQLite](packages/SQLite),
  [PDF](packages/FileEditors/README.md), or
  [CIF](packages/NglViewer/README.md)
* Scientific methods implemented in R, Python, or Julia
* File metadata extractors, such as [Tika](packages/Tika/README.md)
* Custom predictive models that work with the
  built-in [predictive modeling](https://datagrok.ai/help/learn)
  , such as [TensorFlow.js](packages/TensorFlow.js/README.md)
* Platform enhancements, such as [PowerPack](packages/PowerPack/README.md)
  or [UsageAnalysis](packages/UsageAnalysis)
* ... and other types of extensions documented [here](https://datagrok.ai/help/develop/packages/extensions).

These open-source packages are free to use by anyone, although for
the [public environment](https://public.datagrok.ai)
there are some restrictions related to the server computational capacities. Organizations that
deploy Datagrok
[on their premises](https://datagrok.ai/help/develop/under-the-hood/architecture) also can access public packages.
In addition to that, enterprises typically establish their own private repositories that contain
proprietary extensions.

For developers: check out [getting started](https://datagrok.ai/help/develop)
and [contributor's guide](CONTRIB.md).

## Building

The repository is one pnpm workspace with a Turborepo task graph: `grok setup` once (or `pnpm install` at the root),
then `grok build` in any package (or `grok build --all`). Details, and what a package declares, are in
[packages/BUILD.MD](packages/BUILD.MD).

## Academia

Datagrok grants free license to academic institutions to use it in any context, either research or
educational. Moreover, publishing scientific methods as Datagrok packages provides a number of
unique benefits that are specifically important to academia:

* [Reproducible and scalable computations](https://datagrok.ai/help/compute)
* Making your research globally available by
  using [data augmentation](https://datagrok.ai/help/explore/data-augmentation) capabilities. The platform
  proactively suggests contextual actions and enriches the current object
  using [functions](https://datagrok.ai/help/datagrok/concepts/functions)
  implemented in [R, Python, Julia, Matlab, or other language](https://datagrok.ai/help/compute/scripting). In other
  words, Datagrok not only can run a function, but also suggests _what_ could be derived from your
  dataset. This cross-pollination of knowledge could be transformative within and across a broad
  range of scientific disciplines.

For academic collaborations, please email `info@datagrok.ai`.

## Ideas for contributions

If you want to get familiar with the platform, here are some ideas. Pick whatever interests you, and
reach out to Andrew (askalkin@datagrok.ai) or post on
our [community forum](https://community.datagrok.ai/).

* Visualizations
  * Gantt chart
  * Port visjs-based [network diagram](https://datagrok.ai/help/visualize/viewers/network-diagram)
    from Dart to JavaScript
  * WebGL-based rendering of the 2D scatter plot to work with 10M+ points
  * [Event drops](https://github.com/marmelab/EventDrops)
* Scientific methods
  * Statistical hypothesis testing
  * Bayesian statistics
  * Computer vision
  * [NLP](packages/NLP)
* File editors and viewers
* File metadata extractors (see Apache Tika)
* [WASM-based support for digital signal processing](packages/DSP)
* Domain-specific algorithms
* Connectors to web services and open datasets
* [Bioinformatics](packages/Bio)
* Telecom
* Fintech

## See also

* [Datagrok home](https://datagrok.ai/)
* [JavaScript development](https://datagrok.ai/help/develop)
* [Community forum](https://community.datagrok.ai/)
