---
name: work-with-package-files
description: Bundle and read static data files from a Datagrok package — webRoot URLs for `tables/` assets, `_package.files` (FilesDataSource) for `files/` content
---

# work-with-package-files

## When to use

You have demo data, templates, lookup tables, or images you want to ship
inside a Datagrok package — read them from package code without an
external URL or DB connection. Triggers: "load `tables/test.csv` into a
view", "read a JSON template from `files/`", "set a background image
from package assets", "list files under `files/templates/` recursively",
"read a `.d42` dataframe archive at startup."

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `datagrok-api` available — `import * as grok from 'datagrok-api/grok'`,
  `import * as DG from 'datagrok-api/dg'`.
- `_package` symbol exported from `src/package.ts` (every shipping
  package does this; the article's `package.js` mention is a doc
  artifact, drift `DG-FACT-DRIFT-066`).

## Steps

1. **Place data under one of the reserved directories.**
   The package root reserves `files/` for any data (binary, JSON, CSV,
   templates, images) and `tables/` specifically for tabular data
   (`DG-FACT-165`). Subdirectories are allowed in both. Other folder
   names work too — they're served via `webRoot` — but only `files/` is
   wired to the `_package.files` AppData FilesDataSource.
   ```text
   <package>/
     files/                  # → _package.files (AppData), also webRoot-served
       templates/foo.json
       df.csv
     tables/                 # tabular data, webRoot-served
       test.csv
     images/                 # arbitrary; webRoot-served
       night-sky.png
     package.json
     src/package.ts
   ```
   Expected: files are committed in the package source tree; `grok publish`
   uploads `files/` contents to the storage backend named in
   `GROK_PARAMETERS` (e.g., S3) and surfaces them under
   `Files | App Data | <Package>` in the UI (`DG-FACT-168`).

2. **Export `_package` once at the entrypoint.**
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   export * from './package.g';
   export const _package = new DG.Package();
   ```
   Expected: every other module in the package can `import {_package}`
   and call `_package.webRoot` / `_package.files.*` (`DG-FACT-166`).
   Subclassing (`class XPackage extends DG.Package`) is supported when
   you need a custom `init()`.

3. **Use `_package.webRoot` for HTTP-served assets (`tables/`, `images/`, READMEs).**
   `webRoot` is a string getter that returns the URL prefix under which
   the package's source tree is served — the trailing slash IS included,
   so suffix without a leading slash (`DG-FACT-167`).
   ```typescript
   // Load a CSV from `tables/` into a TableView
   grok.data.loadTable(`${_package.webRoot}tables/test.csv`)
     .then((t) => grok.shell.addTableView(t));   // DG-FACT-171

   // Background image from `images/`
   const root = document.createElement('div');
   root.style.backgroundImage = `url(${_package.webRoot}images/night-sky.png)`;

   // README in the help pane
   grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
   ```
   Expected: HTTP GET against the package web root resolves; in
   browser DevTools the request URL is `<webRoot>tables/test.csv`.
   Pattern matches `packages/HitTriage/src/app/pep-triage-views/info-view.ts:22,32`,
   `packages/Notebooks/src/package.js:257`.

4. **Use `_package.files` for content under `files/`.**
   `_package.files` is a fresh `DG.FilesDataSource` rooted at
   `System:AppData/<PackageName>/` — paths passed to its methods are
   relative to the package's source-tree `files/` subtree, NOT to the
   package root (`DG-FACT-168`). The class is `FilesDataSource`;
   `FileSource` (still linked in the article) is a deprecated alias —
   prefer `DG.FilesDataSource` in skill code (`DG-FACT-172`, drift
   `DG-FACT-DRIFT-068`).
   ```typescript
   // List files under files/templates/ recursively — DG-FACT-169
   const entries = await _package.files.list('templates/', true);
   //   FilesDataSource.list(path, recursive=false, searchPattern=null) → FileInfo[]

   // Read text / CSV — DG.DataFrame.fromCsv on the string — DG-FACT-170
   const df = DG.DataFrame.fromCsv(await _package.files.readAsText('df.csv'));

   // Read JSON — no readAsJson; parse the string yourself
   const cfg = JSON.parse(await _package.files.readAsText('templates/foo.json'));

   // Read raw bytes for a binary file
   const bytes: Uint8Array = await _package.files.readAsBytes('test.dat');

   // Read a .d42 dataframe archive (one or more frames; index [0] for the first)
   const archived = (await _package.files.readBinaryDataFrames('project.d42'))[0];
   ```
   Expected: all four reads succeed for files committed under `files/`
   in the package source. Pattern matches
   `packages/HitTriage/src/app/utils.ts:162-168` (list + JSON read),
   `packages/Bio/src/tests/utils.ts:11` (CSV read),
   `packages/DrugBank/src/package.ts:18` (`.d42` at init).

5. **Optional — load a static dataset once at package init.**
   For datasets read on every call (lookup tables, structure libraries),
   read them from `_package.files` inside an `@grok.decorators.init()`
   method on `PackageFunctions` and cache in module scope. Mirrors
   `packages/DrugBank/src/package.ts:15-21`.
   ```typescript
   let dbdf: DG.DataFrame;
   export class PackageFunctions {
     @grok.decorators.init()
     static async initData(): Promise<void> {
       dbdf = (await _package.files.readBinaryDataFrames('drugbank.d42'))[0];
     }
   }
   ```
   Expected: `dbdf` is populated before any function in the package
   runs; the `.d42` file is read from `System:AppData/<Package>/`.

## Common failure modes

- **Path passed to `_package.files.readAsText` is treated as relative
  to package root.** It isn't — it's relative to the package's source
  `files/` subtree (`DG-FACT-168`). Don't pass `'files/df.csv'`; pass
  `'df.csv'`. Verify by browsing `Files | App Data | <Package>` — the
  same path must resolve there.
- **Asset lives outside `files/` but you used `_package.files`.** Only
  `files/` is wired to the FilesDataSource. For `tables/`, `images/`,
  or `README.md`, use `_package.webRoot + 'path'` and an HTTP-style
  loader (`grok.data.loadTable`, `<img>`, `showHelp`) — `DG-FACT-167`,
  `DG-FACT-171`.
- **`webRoot` URL has a double or missing slash.** `webRoot` already
  ends in `/` (`DG-FACT-167`); the suffix must NOT start with `/`. Use
  `${_package.webRoot}tables/x.csv`, not `/tables/x.csv`.
- **Binary file mojibake.** You called `readAsText` on a binary
  payload. Use `readAsBytes` for raw bytes or `readBinaryDataFrames`
  for `.d42` archives (`DG-FACT-170`).
- **Followed `<table>_csv_options.json` sidecar from the article and
  it had no effect.** No public-repo package ships a `*_csv_options.json`
  file — the convention may be stale (drift `DG-FACT-DRIFT-067`). For
  explicit CSV options, pass them per-call to
  `_package.files.readCsv(name, options)` or `grok.data.parseCsv`.
- **Code references `FileSource` and TypeScript autocompletes it.**
  `DG.FileSource` is a deprecated alias for `DG.FilesDataSource`
  (`DG-FACT-172`, drift `DG-FACT-DRIFT-068`). New code should use
  `DG.FilesDataSource`.

## Verification

- `grok publish <host>` exits `0`.
- In the Datagrok UI, `Files | App Data | <Package>` lists everything
  that was committed under the source `files/` directory.
- A test or `init()` that calls `_package.files.list('', true)`
  returns a non-empty `FileInfo[]`; `readAsText` / `readAsBytes` on a
  known entry returns its content.
- A test that loads `${_package.webRoot}tables/<x>.csv` via
  `grok.data.loadTable` returns a `DG.DataFrame` with the expected
  columns.

## See also

- Source articles:
  - `help/develop/how-to/packages/work-with-package-files.md`
  - `help/develop/develop.md#package-structure` (linked from the article)
  - `help/develop/how-to/db/access-data.md` (linked from the article)
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-165` … `DG-FACT-172` and drifts
  `DG-FACT-DRIFT-066`..`DG-FACT-DRIFT-068`.
- Reference packages:
  - `packages/HitTriage/src/app/utils.ts:159-178` — `_package.files.list`
    on a folder, then `readAsText` + `JSON.parse` per entry.
  - `packages/Bio/src/tests/utils.ts:10-19` — `readAsText` + `DataFrame.fromCsv`.
  - `packages/DrugBank/src/package.ts:8-21` — `readBinaryDataFrames` of
    a `.d42` archive inside `@grok.decorators.init()`.
  - `packages/HitTriage/src/app/pep-triage-views/info-view.ts:22,32` —
    `_package.webRoot + 'README.md'` / `+ 'images/...'`.
- Related skills:
  - `file-handlers` (register a parser for a file extension users open).
  - `create-custom-file-viewers` (preview a file as a custom view).
  - `access-data` (parent topic — all the ways to bring data in).
