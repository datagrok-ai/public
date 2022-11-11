<!-- TITLE: Work with package files -->

# How to work with package files

There are multiple ways to access your data in the Datagrok platform. To learn more about them,
please refer to the [main article](access-data.md) covering this topic. Here we will focus on
files that are part of the package. This option is well suited for public data used for
demonstration purposes, as the data is distributed directly with the package.

## Package structure conventions

In a package, the following directories are reserved for data: `files` and `tables`. As
the names suggest, the former is more generic, while the latter is limited to tabular data.
Note that subdirectories are allowed in both cases:

```
files/
  ...
tables/
  ...

package.json
```

These conventions are not set in stone, so you can create other folders for your data, using
the special directories Datagrok is aware of simply allows you to implement more flexible and
maintainable solutions. Let's take dataframes stored under the `tables` folder as an example.
Their import options from an additional file `<table_name>_csv_options.json`, if present,
are taken into account.

## Accessing files in package code

There are several ways to read files that reside within your package. The first is to download
them via HTTP. The package root for client-side can be found with the `webRoot` property of class
[Package](https://datagrok.ai/js-api/classes/dg.Package). Here are some examples:

```js
// `_package` is defined in package.js
export const _package = new DG.Package();

// Get `test.csv` from the `tables` folder and open a table view for it
grok.data.loadTable(`${_package.webRoot}tables/test.csv`)
  .then((t) => grok.shell.addTableView(t));

// Set a background image from the `images` subdirectory
const root = document.createElement('div');
root.style.backgroundImage = `url(${_package.webRoot}images/night-sky.png)`;
```

These methods work with any URL, be it a link to an external resource or a file from a package.
Data under `files`, however, can be accessed via standard file methods exposed in Datagrok JS
API (check out the [FileSource](https://datagrok.ai/js-api/classes/dg.FileSource) class):

```js
// `_package` is defined in package.js
export const _package = new DG.Package();

async function test() {
  // List files in `files/templates/` recursively
  const files = await _package.files.list('templates/', true);

  // Read a dataframe, json, or binary data
  const f1 = DG.DataFrame.fromCsv(await _package.files.readAsText('df.csv'));
  const f2 = JSON.parse(await _package.files.readAsText('template.json'));
  const f3 = await _package.files.readAsBytes('test.dat');
  const f4 = (await _package.files.readBinaryDataFrames('project.d42'))[0];
}
```

The package files are located in a storage specified in `GROK_PARAMETERS` (see
[configuration details](../admin/configuration.md)), for example, in an S3 bucket.
Users can browse these files from `Files | App Data | <Package>`:

![Browsing files](./app-data.gif "Find package data in the file browser")

See also:

* [JavaScript development](../develop.md)
* [Package structure](../develop.md#package-structure)
* [How to access data](access-data.md)
