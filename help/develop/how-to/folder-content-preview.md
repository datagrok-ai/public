<!-- TITLE: Develop custom file handlers -->

# Folder content preview

To provide custom folder content preview, register a function tagged as `folderViewer` that takes two
parameters `folder` and `files`, inspects them and returns a widget if a custom preview could be provided, or null
otherwise.

The following function adds the 'START' button if one of the files in that folder is named "
demog.csv":

```js
//tags: folderViewer
//input: file folder
//input: list<file> files
//output: widget
export function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]): DG.Widget | undefined {
  if (files.some((f) => f.fileName.toLowerCase() == 'demog.csv'))
    return DG.Widget.fromRoot(ui.div([ui.button('START', () => grok.shell.info('Foo'))]));
}
```

This is what you would see when you open a folder that contains "demog.csv":

![folder-content-preview](folder-content-preview.png)

See also

* [File handlers](file-handlers.md)
* [File exporters](file-exporters.md)
