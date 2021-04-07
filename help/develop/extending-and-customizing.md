!-- TITLE: Datagrok Extensions -->
<!-- SUBTITLE: -->

# Extending and customizing Datagrok

Datagrok is built highly extensible, composable and customizable. Many parts of the Datagrok platform can be
enhanced using our [JavaScript API](). The enhancements are structured and delivered to the platform using
[Datagrok packages](js-api.md). Many features of the platform, such as a [Timelines]() or [Sunburst]() viewers,
or a [Cheminformatics package](), are already built as extensions. It is straightforward to create your own ones,
using the [existing extensions](https://github.com/datagrok-ai/public/tree/master/packages) as examples, and
following our [guides](develop.md), [API samples](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples), [how-to's](how-to/develop-custom-viewer.md), and [exercises](exercises.md).

Both the existing extensions and your own-made ones are intended highly customizable. For instance,
look at some properties of the [grid viewer](), where you can set which tooltips to display, which grid lines
to present, and even add [conditional color coding](). Datagrok provides for introducing your own custom
properties specific to your extensions.

## What can be extended

With using our [JavaScript API](), you can create your own:

* [functions](), which may be written in any [scripting language we support](), and later be reused
  in various contexts, including other functions
* [visualizations (viewers)]() — to view tabular data in new ways, in addition to our 30+ viewers
* [file viewers]() — to support new data formats in addition to many we already recognize,
* [cell renderers]() — to visualize certain semantic types, such as [molecules]() or [nucleotide sequences](),
  in their native renders, inside contexts such as a grid cell, a tooltip, or an axis label in a viewer
* [semantic type detectors]() — to attach semantic types to columns of particular data types to later re-use this knowledge
* Web [applications]() focused on specific tasks, such as an interactive dashboard or a data set browser, as [the one for Chembl]()
* [menus](), which may be embed into virtually any context inside the Datagrok UI, such as a top menu or a context menu
  of a viewer
* [info panels](), which augment datasets with all possible kinds of information based on the original dataset contents
* [connections](), to add new public or in-house data sources to the Datagrok instance, such as [Chembl]() or [ENA]()

## What can be customized

Using packages 

<!-- Grid properties, etc. TBD -->