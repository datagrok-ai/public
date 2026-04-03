---
name: create-file-viewer
description: Create a custom file viewer for the Datagrok file share browser
when-to-use: When user asks to create a file viewer, handle a new file type, or add file preview
effort: medium
argument-hint: "[extensions] [package-path]"
---

# Create File Viewer

Create a custom file viewer that is used by the Datagrok file share browser to display files with specific extensions.

## Usage

```
/create-file-viewer [extensions] [package-path]
```

## Instructions

When this skill is invoked, help the user create a custom file viewer for one or more file extensions.
Use the /ui skill for building the UI if necessary.

### Step 1: Create the viewer function

Add a function to `package.ts` (or a separate file imported by it) with the correct annotations. The function must:

- Accept a single `file` input
- Return a `view`
- Have the `fileViewer` tag plus a `fileViewer-<ext>` tag for each supported extension

```ts
//name: myFileViewer
//tags: fileViewer, fileViewer-xyz, fileViewer-abc
//input: file file
//output: view v
export function myFileViewer(file: DG.FileInfo) {
  let view = DG.View.create();

  // Read file content and render it
  file.readAsString().then((content) => {
    let host = ui.div([]);
    // Process and display content
    host.innerText = content;
    view.append(host);
  });

  return view;
}
```

### Step 2: Handle file content

Choose the appropriate method to read file content:

- `file.readAsString()` - for text-based files (returns `Promise<string>`)
- `file.readAsBytes()` - for binary files (returns `Promise<Uint8Array>`)

Example for binary files:

```ts
//tags: fileViewer
//input: file file
//output: view v
//meta.fileViewer: mol, sdf, cif
export function structureViewer(file: DG.FileInfo) {
  let view = DG.View.create();
  let host = ui.div([], 'd4-ngl-viewer');

  file
    .readAsBytes()
    .then(bytes => {
      // Process binary content
      let blob = new Blob([bytes]);
      // Render blob...
    });

  view.append(host);
  return view;
}
```

### Key points

- The function runs when a user clicks on a file with a matching extension in the file share browser
- The view is displayed in the file preview panel
- Return the view synchronously; load file content asynchronously inside it
- `meta.fileViewer` can take multiple extensions

## For extensions that are efficiently DataFrames

- Return a `DG.TableView` directly — do not wrap a grid in a custom `DG.View`
- Transfer all source metadata to tags (table-level and column-level) so it survives
  round-trips through the platform

```typescript
// Good — file viewer for a single-table format
static async previewFoo(file: DG.FileInfo): Promise<DG.View> {
  const bytes = await file.readAsBytes();
  const data = await parseFoo(bytes);
  const df = toDataFrame(data);
  const view = DG.TableView.create(df, false);
  view.name = file.name;
  return view;
}
```

- Set `df.setTag('source.format', '...')` to record where the data came from
- Preserve column-level metadata (description, format, categories) as tags
- Use a namespace prefix for format-specific tags (e.g., `minitab.version`, `prism.sheetId`)
- Use the standard `description` tag for column descriptions

## For extensions that contain multiple dataframes

- Use custom `DG.View`
- Show global file metadata on top
- 'Open all' big button on the ribbon
- Tab control with pane per dataframe 

## Behavior

1. Ask the user which file extension(s) they want to support if not specified
2. Determine whether the file is text-based or binary to choose `readAsString` vs `readAsBytes`
3. Create the viewer function with correct tags for all specified extensions
4. Add rendering logic appropriate for the file type

