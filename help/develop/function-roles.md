---
title: "Function roles"
---

A package can contain functions that get discovered and integrated with the platform at runtime.
Typically, each function type has a special tag denoting what the function does, for example:

* `#app` for [applications](#applications)
* `#dashboard` for [dashboards](../visualize/dashboard.md)
* `#panel` for [info panels](#info-panels)
* `#init` for [package initialization](#package-initialization)
* `#autostart` for [automatic execution at platform startup](#autostart)
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

Applications are [functions](../datagrok/functions/functions.md) tagged with the `#app` tag. Use `datagrok-tools` to get
a template:

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

## Package initialization

`init` function gets invoked before the first time any of the package functions is invoked. This is a good place to
initialize some common structures
(load WebAssembly, fetch files) that some of exposed functions use. It gets invoked at most once.

Use `datagrok-tools` to get a template:

```shell
cd <package-name>
grok add function init <packageName>Init
```

See also [autostart](#autostart).

## Autostart

`autostart` function get invoked at platform startup. It starts immediately if it is tagged
as `meta.autostartImmediate: true`, otherwise it starts three seconds later.

Use the `meta.autostartImmediate: true` mode is good to subscribe to global events, change some default settings, or
change the UI right at the platform start. The default mode is good for preloading popular big packages (for instance,
Chem package preloads to eliminate the delays when a user opens a file with molecules). Keep in mind that some packages
are big, and unless the `autostart` function resides in the `detectors.js` file, the whole content of the package gets
loaded, so use this option wisely.

The `autostart` function can reside in the `detectors.js` file, in which case full package content is not loaded. You
might want to use it to quickly add some menu items or start listening to some events, and then you can load the full
package only when a user launches that functionality.

If the `autostart` function is defined in the package where the [`init`](#package-initialization)
is also defined, the `init` function gets executed first. No one awaits on the `autostart` function. You might have
zero, one, or more `autostart` functions in a package.

Example from the [PowerPack](https://github.com/datagrok-ai/public/tree/master/packages/PowerPack)
plugin that introduces a "Welcome" view at startup:

```js
//name: welcomeView
//tags: autostart
//meta.autostartImmediate: true
export function _welcomeView(): void {
  if (_properties['showWelcomeView'])
    welcomeView();
}
```

## Semantic type detectors

Functions that define semantic types have the `#semTypeDetector` tag. Use `datagrok-tools` to get a template:

```shell
cd <package-name>
grok add detector <semantic-type-name>
```

*Details:* [How to define semantic type detectors](how-to/define-semantic-type-detectors.md)

## Cell renderers

Cell renderers allow customizing the appearance of cells in the [grid](../visualize/viewers/grid.md). These functions
are annotated with two special tags: `cellRenderer` and `cellRenderer-<type>`.

*Details:* [How to develop custom cell renderers](how-to/custom-cell-renderers.md)

## File viewers

File viewers are used in Datagrok's [file share browser](../access/connect-a-file-share.md). The platform provides a way
to define custom viewers (or editors) in addition to the native ones. These functions work on files with a specific
extension, which is derived from the `fileViewer-<extension>` tag.

*Details:* [How to develop custom file viewers](how-to/create-custom-file-viewers.md)

## File exporters

A file exporter is a function used for loading data from the platform. It is annotated with the `#fileExporter` tag.
Exporters reside in the platform's top menu "export" section.

*Details:* [How to create file exporters](how-to/file-exporters.md)

## Settings editors

Settings editors work with [package properties](develop.md#package-settings) and define how they will be displayed in
the `Settings` pane of the property panel. An editor function should return a widget (`DG.Widget`) and be tagged as
`#packageSettingsEditor`.

*Details:* [How to write custom package settings editors](how-to/custom-package-settings-editors.md)
