<!-- TITLE: Function roles-->

# Function roles

A package can contain a variety of functions, so it will be appropriate to give an overview of the most common ones.
Typically, each function type has a special tag denoting what the function does, for example:

* `#app` for [applications](#applications)
* `#dashboard` for [dashboards](#dashboards)
* `#panel` for [info panels](#info-panels)
* `#init` and `#autostart` for [pre-run functions](#pre-run-functions)
* `#semTypeDetector` for [semantic types detectors](#semantic-type-detectors)
* `#cellRenderer` for custom [cell renderers](#cell-renderers)
* `#fileViewer` and `#fileExporter` for [file viewers](#file-viewers)
  and [exporters](#file-exporters)
* `#packageSettingsEditor` for [settings editors](#settings-editors)

You can use these tags to search for certain functions either from the platform's interface
([https://public.datagrok.ai/functions?q](https://public.datagrok.ai/functions?q)) or from within your code:

```js
const applications = DG.Func.find({tags: [DG.FUNC_TYPES.APP]});
```

**TIP** To disable all package functions (for debug purposes), use the
`initPackageFunctions=false` flag in the start URL, such as
`https://public.datagrok.ai?initPackageFunctions=false`.

## Applications

Applications are [functions](../overview/functions/function.md) tagged with the `#app` tag.
Use `datagrok-tools` to get a template:

```shell
cd <package-name>
grok add app <name>
```

*Details:* [How to build an application](how-to/build-an-app.md)

## Info panels

Functions tagged as `#panel` extend the property panel with additional content for the current object.
Use `datagrok-tools` to get a template:

```shell
cd <package-name>
grok add function panel <name>
grok add script panel <language> <name>
```

*Details:* [How to add an info panel](how-to/add-info-panel.md)

## Pre-run functions

The purpose of pre-run functions is to prepare the main package code for execution. This includes fetching specific
pieces of data, subscribing to global events, changing the user interface right after the platform starts, connecting to
external services with refined configuration parameters, and so on.

There are two types of functions serving this purpose: `init` and `autostart`. The function tagged with `init` gets
invoked when the containing package is initialized. This typically happens the first time any of the functions in the
package is called. It is guaranteed that this function gets invoked *once* before the execution of the rest of the code
and will not be re-executed on subsequent calls. The `autostart` functions are similar to the first type, but differ
from it in a few aspects. Firstly, these functions are called at the platform startup, not necessarily when some package
function is invoked. Moreover, if you decide to call a regular function from your package, there is no guarantee that
its code will wait until the `autostart` completes. Another caveat is that the whole package will get initialized along
with `autostart`, so use this type of functions wisely. If possible, stick to the `init` tag while developing your
programs.

Use `datagrok-tools` to get a template:

```shell
cd <package-name>
grok add function init <packageName>Init
```

## Semantic type detectors

Functions that define semantic types have the `#semTypeDetector` tag.
Use `datagrok-tools` to get a template:

```shell
cd <package-name>
grok add detector <semantic-type-name>
```

*Details:* [How to define semantic type detectors](how-to/define-semantic-type-detectors.md)

## Cell renderers

Cell renderers allow customizing the appearance of cells in the [grid](../visualize/viewers/grid.md).
These functions are annotated with two special tags: `cellRenderer` and `cellRenderer-<type>`.

## File viewers

File viewers are used in Datagrok's [file share browser](../access/file-shares.md). The platform provides a way to
define custom viewers (or editors) in addition to the native ones. These functions work on files with a specific
extension, which is derived from the `fileViewer-<extension>` tag.

*Details:* [How to develop custom file viewers](how-to/custom-file-viewers.md)

## File exporters

A file exporter is a function used for loading data from the platform. It is annotated with the `#fileExporter` tag.
Exporters reside in the platform's top menu "export" section.

*Details:* [How to create file exporters](how-to/file-exporters.md)

## Settings editors

Settings editors work with [package properties](#package-settings) and define how they will be displayed in
the `Settings` pane of the property panel. An editor function should return a widget (`DG.Widget`) and be tagged as
`#packageSettingsEditor`.

*Details:* [How to write custom package settings editors](how-to/custom-package-settings-editors.md)
