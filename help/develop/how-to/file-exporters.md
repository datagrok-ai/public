<!-- TITLE: Create file exporters -->

# File exporters

Part of Datagrok's functionality is built for easy file management. File exporters, along
with [file viewers](custom-file-viewers.md), provide an example of such features. A file exporter is a function used for
loading data from the platform. Once registered, it appears at the file's "
export" menu:

![Save as SDF](file-exporter.gif "Save as SDF")

To write an exporter function, add a `fileExporter` tag to its annotation. Make sure to specify the description: it will
be used as the menu entry, e.g.,
`As ${fileExtension}`. Unlike file viewers, exporters don't need an extension-specific tag, you start with an open
table, modify it, and let the user download the converted version. Let's have a look at a function from the
[Chem](https://github.com/datagrok-ai/public/blob/73356b9c34e28fcd2278a8f60137c1c90684c8f3/packages/Chem/package.js)
package that exports a dataframe in a special file format for chemical data:

```js
//name: saveAsSdf
//description: Save as SDF
//tags: fileExporter
saveAsSdf() {
  let table = grok.shell.t;
  let structureColumn = table.columns.bySemType('Molecule');
  if (structureColumn == null)
    return;

  let result = '';

  for (let i = 0; i < table.rowCount; i++) {
    try {
      let mol = new OCL.Molecule.fromSmiles(structureColumn.get(i));
      result += `\n${mol.toMolfile()}\n`;

      for (let col of table.columns)
        if (col !== structureColumn) {
          result += `>  <${col.name}>\n${col.get(i)}\n\n`;
        }

      result += '$$$$'
    }
    catch (error) {
      console.error(error);
    }
  }

  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(result));
  element.setAttribute('download', table.name + '.sdf');
  element.click();
}
```

In this function, we obtain a currently open table with `grok.shell.t`, work on the output format, create an anchor
element with a filename in the `download`
attribute, and lastly fire the element's click event. Note that we don't return anything from the exporter. After the
package publication, the registered function will get attached to the file export menu at the platform's startup.

See also:

* [JavaScript development](../develop.md)
* [How to develop custom file viewers](custom-file-viewers.md)
