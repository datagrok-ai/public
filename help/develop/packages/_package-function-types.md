<!-- TITLE: Function Types -->
<!-- ORDER: 4 -->

# Function Types

A package can contain various functions, and each function is annotated with a tag that defines what this function does:

* `#app` for [applications](#applications)
* `#dashboard` for [dashboards](#dashboards)
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

Applications are [functions](../../overview/functions/function.md) tagged with `app`:

```js
class EnamineStorePackage extends DG.Package {

  //tags: app
  //name: Enamine Store
  startApp(context) {
  }

  enamineStore(smiles) {
  }

  createSearchPanel(panelName, smiles) {
  }

  static dataToTable(data, name) {
  }

  static searchModeToCommand(name) {
  }
}
```

**See more**: [Enamine Store application]

On Datagrok, you can view the applications in the [application launcher]. To launch a particular app automatically, go 
to `https://public.datagrok.ai/apps/<APP_NAME>`.

To create a template for an `#app` function, from your package directory, run:

```shell
grok add app <application-name>
```

For more information on developing an application, refer to the [How to build an application section](../how-to/build-an-app.md)

## Dashboards

To be continued...

## Info Panels

To be continued...

## Pre-run Functions

The purpose of pre-run functions is to prepare the main package for execution. This includes fetching specific pieces of 
data from an API, subscribing to global events, changing the user interface after the platform starts, connecting 
to external services with refined configuration parameters, and so on.

There are two types of functions that can prepare the package for execution &mdash; `init` and `autostart`. A function 
tagged with `#init` gets invoked when the package is initialized.

The initialization phase happens the first time any functions in the package is called. Datagrok guarantees that the 
`#init` function is invoked _once_ before the rest of the code. 

An `autostart` function is similar to the the `init` function, but is called at the platform startup, not necessarily 
when a package function is invoked. Hence, an `#autostart` function is called earlier than an `#init` function.

Another difference between `#autostart` and `#init` is that when you call a regular function from the package, there's 
no guarantee that its code will wait until the `autostart` completes. Another caveat is that the whole package will get 
initialized along with `autostart`, so use this type of functions wisely. If possible, stick to the `init` tag while 
developing your programs.

To generate an `init` function, from your package directory, run:

```shell
grok add function init <packageName>Init
```

## Semantic Type Detectors

To be continued...

## Cell Renderers

To be continued...

## File Viewers

File viewers are used in Datagrok's [file share browser](../../access/file-shares.md).
The platform provides a way to define custom viewers or editors in addition to the native ones.
These functions work with files with a specific extension, which is derived from the `fileViewer-<extension>` tag.

*See more:* [How to Develop Custom File Viewers](../how-to/custom-file-viewers.md)

## File Exporters

A file exporter is a function used for loading data from the platform. It is annotated with the `#fileExporter` tag. 
Exporters reside in the platform's top menu "export" section.

*See more:* [How to Create File Exporters](../how-to/file-exporters.md)

## Package Settings Editors

To be continued...

[Datagrok Public > Functions]: https://public.datagrok.ai/functions?q
[Datagrok GitHub]: https://github.com/datagrok-ai/public/tree/master/packages
[Enamine Store application]: https://github.com/datagrok-ai/public/tree/master/packages/EnamineStore
[the direct link]: https://public.datagrok.ai/apps