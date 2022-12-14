<!-- TITLE: &#8204;Extending Datagrok -->
<!-- SUBTITLE: -->

# Extending and customizing Datagrok

Datagrok is built highly extensible, composable and customizable. Many parts of the Datagrok platform can be enhanced by
plugins using our [JavaScript API](../js-api.md). The plugins are structured and delivered to the platform
using [Datagrok packages](../develop.md#packages). Many features of the platform, such as a
[Timelines](https://github.com/datagrok-ai/public/tree/master/packages/Viewers) or
[cheminformatics package](https://github.com/datagrok-ai/public/tree/master/packages/Chem), are already built as
plugins. It is straightforward to create your own ones, using
the [existing packages](https://github.com/datagrok-ai/public/tree/master/packages) as examples, and following
our [guides](../develop.md), [API samples](https://public.datagrok.ai/js)
, [how-to's](../how-to/develop-custom-viewer.md), and [exercises](../exercises/exercises.md).

## What can be extended

With using our [JavaScript API](../js-api.md), you can create your own:

* [functions](../../datagrok/functions/function.md), which may be written in any
  [scripting language we support](../../compute/scripting.md), and later be reused in various contexts, including other
  functions, or called directly from Datagrok UI or the [console](../../datagrok/navigation.md#console)
* [visualizations](../../visualize/viewers.md) — to visualize data in new ways, in addition to our 30+ core viewers
* [file viewers](../how-to/create-custom-file-viewers.md) — to support new data formats in addition to many we already recognize
* [cell renderers](../function-roles.md#cell-renderers) — to visualize certain semantic types, such
  as [molecules
  or nucleotide sequences](https://github.com/datagrok-ai/public/blob/master/libraries/chem-meta/src/rdkit-api.ts)
  , in their native-looking renders, inside contexts such as a grid cell, a tooltip, or an axis label in a viewer
* [semantic type detectors](../how-to/define-semantic-type-detectors.md) — to attach semantic types to columns of
  particular data types to later re-use this knowledge
* Web [applications](../how-to/build-an-app.md) focused on specific tasks, such as an interactive dashboard or a data set
  browser<!--, as [the one for Chembl](https://github.com/datagrok-ai/public/tree/master/packages/ChemblBrowser)-->
* menus, which may be embedded into virtually any context inside the Datagrok UI, such as a
  [top menu](https://public.datagrok.ai/js/samples/ui/menu) or
  a [context menu](https://public.datagrok.ai/js/samples/events/viewer-events) of a viewer
* [info panels](../how-to/add-info-panel.md), to augment data with additional information retrieved or calculated
  on-the-fly
* [connections](../../access/data-connection.md), to add new public or in-house data sources, such
  as [Chembl](https://www.ebi.ac.uk/chembl/) or [ENA](https://www.ebi.ac.uk/ena/browser/),
  custom [filters](https://github.com/datagrok-ai/public/blob/master/packages/Widgets/src/filters/radio-button-filter.ts)
  , which allow adding a filtering mask to an active dataset; in addition, viewers themselves may act as filters, and
  filtering through one viewer shall reflect the state in all the other active viewers
* [accordion sections](../advanced/ui.md#accordions) — accordion is an area on the left of the Datagrok UI, useful for additional
  custom functionality

### What's next

Learn more about the platform by watching our educational content and hands-on practicing:

* [Getting Started Video Walkthrough](../getting-started.md#6-videos)
* [Development Exercises](../exercises/exercises.md) to practice developing on Datagrok in a good sequence and through an
  interesting story
* [Community forum](https://community.datagrok.ai/) to discuss anything

See also:

* [JS API](../js-api.md)
* [JS API Samples](https://public.datagrok.ai/js) (expand "JavaScript API Examples")
* [Public package repository](https://github.com/datagrok-ai/public)
