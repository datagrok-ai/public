---
title: "Datagrok package repository"
---

This is a public repository for the API, tools, and [packages](https://datagrok.ai/help/develop/develop#packages)
available for [Datagrok™](https://datagrok.ai), a next-generation web-based data analytics platform. The platform is
very extensible, and almost any functionality can be implemented as a package:

* Support for scientific domains, such
  as [cheminformatics](https://github.com/datagrok-ai/public/tree/master/packages/Chem/README.md)
* Applications, such
  as [Clinical Case](https://github.com/datagrok-ai/public/tree/master/packages/ClinicalCase/README.md)
  or [Peptides](https://github.com/datagrok-ai/public/tree/master/packages/Peptides/README.md)
* Connectors to [OpenAPI web services](https://github.com/datagrok-ai/public/tree/master/packages/Samples/swaggers)
* Visualizations, such as [Charts](https://github.com/datagrok-ai/public/blob/master/packages/Charts/README.md), [Biostructure Viewer](https://github.com/datagrok-ai/public/blob/master/packages/BiostructureViewer/README.md), or [GIS](https://github.com/datagrok-ai/public/blob/master/packages/GIS/README.md)
* Importing and previewing files, such as
  [SQLite](https://github.com/datagrok-ai/public/tree/master/packages/SQLite),
  [PDF](https://github.com/datagrok-ai/public/tree/master/packages/FileEditors/README.md), or
  [JDX](https://github.com/datagrok-ai/public/blob/master/packages/nmrium/README.md)
* Scientific methods implemented in R, Python, or Julia
* Custom predictive models that work with the built-in [predictive modeling](../learn/learn.md)
  , such as [EDA](https://github.com/datagrok-ai/public/blob/master/packages/EDA/README.md)
* Platform enhancements, such
  as [PowerPack](https://github.com/datagrok-ai/public/tree/master/packages/PowerPack/README.md)
  or [UsageAnalysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis)
* ... and other types of extensions documented [here](../develop/packages/extensions.md).

These open-source packages are free to use by anyone, although for the [public environment](https://public.datagrok.ai)
there are some restrictions related to the server computational capacities. Organizations that deploy Datagrok
[on their premises](../deploy/deploy.md) also can access public packages. In addition to that, enterprises
typically establish their own private repositories that contain proprietary extensions.

For developers: check out [getting started](../develop/develop.md)
and [contributor's guide](https://github.com/datagrok-ai/public/blob/master/CONTRIB.md).

## Academia

Datagrok grants free license to academic institutions to use it in any context, either research or educational.
Moreover, publishing scientific methods as Datagrok packages provides a number of unique benefits that are specifically
important to academia:

* [Reproducible and scalable computations](../compute/compute.md)
* Making your research globally available by using [data augmentation](../explore/data-augmentation/data-augmentation.md) capabilities.
  The platform proactively suggests contextual actions and enriches the current object
  using [functions](../datagrok/concepts/functions/functions.md)
  implemented in [R, Python, Julia, Matlab, or other language](../compute/scripting/scripting.mdx). In other words, Datagrok not
  only can run a function, but also suggests _what_ could be derived from your dataset. This cross-pollination of
  knowledge could be transformative within and across a broad range of scientific disciplines.

For academic collaborations, please email `info@datagrok.ai`.

## See also

* [Datagrok home](https://datagrok.ai/)
* [JavaScript development](../develop/develop.md)
* [Community forum](https://community.datagrok.ai/)
