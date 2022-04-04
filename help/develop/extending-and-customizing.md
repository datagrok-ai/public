<!-- TITLE: &#8204;Extending Datagrok -->
<!-- SUBTITLE: -->

# Extending and customizing Datagrok

Datagrok is built highly extensible, composable and customizable. Many parts of the Datagrok platform can be enhanced by
plugins using our [JavaScript API](js-api.md). The plugins are structured and delivered to the platform
using [Datagrok packages](develop.md#packages). Many features of the platform, such as a
[Timelines](https://github.com/datagrok-ai/public/tree/master/packages/Viewers) or
[Sunburst](https://github.com/datagrok-ai/public/tree/master/packages/Sunburst) viewers, as well
as [cheminformatics package](https://github.com/datagrok-ai/public/tree/master/packages/Chem), are already built as
plugins. It is straightforward to create your own ones, using
the [existing packages](https://github.com/datagrok-ai/public/tree/master/packages) as examples, and following
our [guides](develop.md), [API samples](https://public.datagrok.ai/js)
, [how-to's](how-to/develop-custom-viewer.md), and [exercises](exercises/exercises.md).

Both the existing extensions and your own-made ones are intended highly customizable. For instance, look at some
properties of the [grid viewer](../visualize/viewers/grid.md), where you can set which tooltips to display, which grid
lines to present, and even add [conditional color coding](../visualize/viewers/grid.md#color-coding)! Datagrok also
provides for introducing your own custom properties specific to your plugins.

## What can be extended

With using our [JavaScript API](js-api.md), you can create your own:

* [functions](../overview/functions/function.md), which may be written in any
  [scripting language we support](../compute/scripting.md), and later be reused in various contexts, including other
  functions, or called directly from Datagrok UI or the [console](../overview/navigation.md#console)
* [visualizations](../visualize/viewers.md) — to visualize data in new ways, in addition to our 30+ core viewers
* [file viewers](how-to/custom-file-viewers.md) — to support new data formats in addition to many we already recognize
* [cell renderers](../visualize/viewers/grid.md#custom-cell-renderers) — to visualize certain semantic types, such
  as [molecules](https://github.com/datagrok-ai/public/blob/master/packages/Chem/src/rdkit-api.ts)
  or [nucleotide sequences](https://github.com/datagrok-ai/public/tree/master/packages/Sequence/web-logo-viewer)
  , in their native-looking renders, inside contexts such as a grid cell, a tooltip, or an axis label in a viewer
* [semantic type detectors](how-to/define-semantic-type-detectors.md) — to attach semantic types to columns of
  particular data types to later re-use this knowledge
* Web [applications](how-to/build-an-app.md) focused on specific tasks, such as an interactive dashboard or a data set
  browser, as [the one for Chembl](https://github.com/datagrok-ai/public/tree/master/packages/ChemblBrowser)
* menus, which may be embedded into virtually any context inside the Datagrok UI, such as a
  [top menu](https://public.datagrok.ai/js/samples/ui/menu) or
  a [context menu](https://public.datagrok.ai/js/samples/events/viewer-events) of a viewer
* [info panels](how-to/add-info-panel.md), to augment data with additional information retrieved or calculated
  on-the-fly
* [connections](../access/data-connection.md), to add new public or in-house data sources, such
  as [Chembl](https://www.ebi.ac.uk/chembl/) or [ENA](https://www.ebi.ac.uk/ena/browser/),
  custom [filters](https://github.com/datagrok-ai/public/blob/master/packages/Widgets/src/filters/radio-button-filter.ts)
  , which allow adding a filtering mask to an active dataset; in addition, viewers themselves may act as filters, and
  filtering through one viewer shall reflect the state in all the other active viewers
* [accordion sections](ui.md#accordions) — accordion is an area on the left of the Datagrok UI, useful for additional
  custom functionality

You could spot other aspects in Datagrok that allow for customization while getting introduced to the platform.

## What can be customized

Beside a regular CSS-based customization, here are some of the things which you can customize both programmatically and
through the UI:

* Every viewer exposes a diverse set of [properties](../overview/navigation.md#properties), accessible on the right side
  in the Datagrok's
  `Property Panel`. Often you need to modify these properties programmatically. The mapping is straightforward: take the
  property's name, say, `"Col Header Height"`, modify it to the camel case: `colHeaderHeight`, and apply the property to
  your viewer via `.setOptions`. Say, for the grid `grid` call `grid.setOptions({ colHeaderHeight: 50 });`. Read more
  about this
  [here](../develop/how-to/develop-custom-viewer.md)

* You can control the platform's shell aspects programmatically.
  Study [this sample](https://public.datagrok.ai/js/samples/shell/ui-parts)
  to learn how to hide the Datagrok's sidebar or the toolbox, which is often useful in
  the [Datagrok applications](how-to/build-an-app.md)

* Use features like [conditional color coding](https://dev.datagrok.ai/js/samples/grid/color-coding-conditional)
  for making your views appear indicative to the data they exhibit

* Use the [dock system](how-to/manipulate-viewers.md#docking-viewers)
  and [custom views](how-to/custom-views.md)
  to provide for completely custom layouts composed inside the Datagrok platform, tailored to your specific use cases.
  Learn more about Datagrok's composable UI [here](ui.md)

### What's next

Learn more about the platform by watching our educational content and hands-on practicing:

* [Getting Started Video Walkthrough](getting-started.md#6-videos)
* [Development Exercises](exercises/exercises.md) to practice developing on Datagrok in a good sequence and through an
  interesting story
* [Community forum](https://community.datagrok.ai/) to discuss anything

See also:

* [JS API](js-api.md)
* [JS API Samples](https://public.datagrok.ai/js) (expand "JavaScript API Examples")
* [Public package repository](https://github.com/datagrok-ai/public)
