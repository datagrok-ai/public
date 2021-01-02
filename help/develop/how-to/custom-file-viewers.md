<!-- TITLE: Develop Custom File Viewers -->

# File Viewers

Datagrok provides a way to define custom file viewers that are used by the 
[file share browser](../../access/file-shares.md).
This could be done by define a function annotated in a special way. It should take a single 
argument of type `file`, return a `view`, and have at least two tags: `fileViewer`
and `fileViewer-<extension>` (specify the extension here). This is it!

The following example defines a file viewer for .mol, .sdf, and .cif files by visualizing them
with the NglViewer

```js
//tags: fileViewer, fileViewer-mol, fileViewer-sdf, fileViewer-cif
//input: file file
//output: view v
nglStructureViewer(file) {
  let view = DG.View.create();
  var host = ui.div([], 'd4-ngl-viewer');
  var stage = new NGL.Stage(host);

  file
    .readAsBytes()
    .then(bytes => stage.loadFile(new Blob([bytes])));

  view.append(host);
  return view;
}

```

And here is the result:

![](../../access/file-shares-file-viewers.gif)

See also:
* [File shares](../../access/file-shares.md)