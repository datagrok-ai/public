<!-- TITLE: Function Types -->
<!-- ORDER: 4 -->

# Function types

A package can contain various functions, and each function is annotated with a tag that defines what this function does:

* `#app` for [applications](#applications)
* `#dashboard` for [dashboards](../../datagrok/concepts/project/dashboard.md)
* `#panel` for [info panels](#info-panels)
* `#init` and `#autostart` for [pre-run functions](#pre-run-functions)
* `#semTypeDetector` for [semantic types detectors](#semantic-type-detectors)
* `#cellRenderer` for custom [cell renderers](#cell-renderers)
* `#fileViewer` for [file viewers](#file-viewers)
* `#fileExporter` for [exporters](#file-exporters)
* `#packageSettingsEditor` for [settings editors](#package-settings-editors)

You can use the tags to search for certain functions from [Datagrok > Functions] or from within your code:

```js
const applications = DG.Func.find({ tags: [DG.FUNC_TYPES.APP] });
```

## Applications

Applications are [functions](../../datagrok/functions/function.md) tagged with `app`:

```js
//name: Test Application
//tags: app
export function app() {
  grok.shell.info('This is Test Application!');
}
```

**See more**: [BioSignals application], [ClinicalCase application]

On Datagrok, you can view applications in the [application launcher]. To launch a particular app automatically, go
to `https://public.datagrok.ai/apps/<APP_NAME>`.

A package can contain more than one application. To create a template for an `#app` function, from your package
directory, run:

```shell
grok add app <application-name>
```

For more information on developing an application, refer to
the [How to build an application section](../how-to/build-an-app.md)

## Dashboards

To be continued...

## Info panels

To be continued...

## Pre-run functions

The purpose of pre-run functions is to prepare the main package for execution. This includes fetching specific pieces of
data from an API, subscribing to global events, changing the user interface after the platform starts, connecting to
external services with refined configuration parameters, and so on.

There are two types of functions that can prepare the package for execution &mdash; `init`
and `autostart`. A function tagged with `#init` gets invoked when the package is initialized. This happens the first
time any of the functions in the package is called. It is guaranteed that the `#init` function is invoked _once_ before
the execution of the rest of the code and will not be re-executed on subsequent calls.

An `autostart` function is similar to the the `init` function, but is called at the platform startup, not necessarily
when a package function is invoked. Another difference between `#autostart`
and `#init` is that when you call a regular function from the package, there's no guarantee that its code will wait
until the `autostart` completes. Another caveat is that the whole package will get initialized along with `autostart`,
so use this type of functions wisely. If possible, stick to the `init` tag while developing your programs.

To generate an `init` function, from your package directory, run:

```shell
grok add function init <packageName>Init
```

## Semantic type detectors

To be continued...

## Cell renderers

To be continued...

## File viewers

File viewers are used in Datagrok's [file share browser](../../access/files/files.md). The platform provides a way to
define custom viewers or editors in addition to the native ones. These functions work with files with a specific
extension, which is derived from the `fileViewer-<extension>` tag.

*See more:* [How to Develop Custom File Viewers](../how-to/create-custom-file-viewers.md)

## File exporters

A file exporter is a function used for loading data from the platform. It is annotated with the `#fileExporter` tag.
Exporters reside in the platform's top menu "export" section.

*See more:* [How to Create File Exporters](../how-to/file-exporters.md)

## Package settings editors

To be continued...

[Datagrok > Functions]: https://public.datagrok.ai/functions?q

[Datagrok GitHub]: https://github.com/datagrok-ai/public/tree/master/packages

[application launcher]: https://public.datagrok.ai/apps

[BioSignals application]: https://github.com/datagrok-ai/public/tree/master/packages/BioSignals

[ClinicalCase application]: https://github.com/datagrok-ai/public/tree/master/packages/ClinicalCase

[the direct link]: https://public.datagrok.ai/apps
