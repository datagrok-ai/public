---
name: work-with-package-files
version: 0.1.0
description: |
  Bundle and read static data (demo CSVs, lookup tables, JSON
  templates, images, `.d42` archives) from inside a Datagrok package
  without an external URL or DB query. For plugin authors who want to
  ship reference assets alongside their TypeScript and read them at
  runtime. Produces a package layout under `files/` and/or `tables/`
  plus the calls (`_package.webRoot`, `_package.files.*`) that resolve
  them.
  Use when asked to "ship a reference CSV with my plugin", "load a
  bundled lookup table at startup", or "set a background image from
  packaged assets".
triggers:
  - bundle a reference csv with my plugin
  - ship sample data alongside the package
  - load a lookup table at startup
  - read a json template from disk
  - background image from packaged assets
  - open a packaged d42 archive
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# work-with-package-files

## When to use

You have demo data, templates, lookup tables, images, or `.d42`
archives you want to ship inside a Datagrok package — read them from
package code without an external URL or DB connection.

## Prerequisites

- `import * as grok from 'datagrok-api/grok'`,
  `import * as DG from 'datagrok-api/dg'`.
- `_package` exported from `src/package.ts` (the article's `package.js`
  mention is a doc artifact — drift `DG-FACT-DRIFT-066`).

## Steps

1. **Place data under one of the reserved directories.**
   Only `files/` is wired to `_package.files` (AppData
   FilesDataSource); `tables/` / `images/` and anything else are
   served via `webRoot` (see `DG-FACT-165`, `DG-FACT-168`).
   ```text
   <package>/
     files/templates/foo.json
     files/df.csv
     tables/test.csv
     images/night-sky.png
     package.json
     src/package.ts
   ```
   `grok publish` uploads `files/` to AppData and surfaces it under
   `Files | App Data | <Package>`.

2. **Export `_package` once at the entrypoint.**
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   export * from './package.g';
   export const _package = new DG.Package();
   ```
   Other modules `import {_package}` and call `_package.webRoot` /
   `_package.files.*` (see `DG-FACT-166`). Subclass `DG.Package` for a
   custom `init()`.

3. **Use `_package.webRoot` for HTTP-served assets** (`tables/`,
   `images/`, READMEs). The getter returns a URL prefix WITH the
   trailing slash — never lead the suffix with `/` (see `DG-FACT-167`).
   ```typescript
   // CSV from tables/ → TableView
   grok.data.loadTable(`${_package.webRoot}tables/test.csv`)
     .then((t) => grok.shell.addTableView(t));   // DG-FACT-171

   // Background image from images/
   document.body.style.backgroundImage =
     `url(${_package.webRoot}images/night-sky.png)`;

   // README in the help pane
   grok.shell.windows.help.showHelp(_package.webRoot + 'README.md');
   ```
   DevTools shows the request URL as `<webRoot>tables/test.csv`. Refs:
   `packages/HitTriage/.../info-view.ts:22,32`.

4. **Use `_package.files` for content under `files/`.**
   `_package.files` is a `DG.FilesDataSource` rooted at
   `System:AppData/<PackageName>/`; paths are relative to the `files/`
   subtree, NOT the package root (see `DG-FACT-168`, `DG-FACT-172`).
   ```typescript
   // List recursively — list(path, recursive=false, searchPattern=null) (DG-FACT-169)
   const entries = await _package.files.list('templates/', true);

   // CSV → DataFrame via DG.DataFrame.fromCsv on the string (DG-FACT-170)
   const df = DG.DataFrame.fromCsv(await _package.files.readAsText('df.csv'));

   // JSON — no readAsJson; parse the string yourself
   const cfg = JSON.parse(await _package.files.readAsText('templates/foo.json'));

   // Raw bytes
   const bytes: Uint8Array = await _package.files.readAsBytes('test.dat');

   // .d42 archive — array of one or more frames; pick [0] for the first
   const archived = (await _package.files.readBinaryDataFrames('project.d42'))[0];
   ```
   Refs: `packages/HitTriage/src/app/utils.ts:162-168` (list +
   readAsText + JSON.parse); `packages/DrugBank/src/package.ts:18`
   (`.d42` at init).

5. **Optional — load a static dataset once at package init.**
   Read inside `@grok.decorators.init()` and cache in module scope.
   Mirrors `packages/DrugBank/src/package.ts:10-21`.
   ```typescript
   let dbdf: DG.DataFrame;
   export class PackageFunctions {
     @grok.decorators.init()
     static async initData(): Promise<void> {
       dbdf = (await _package.files
         .readBinaryDataFrames('drugbank-open-structures.d42'))[0];
     }
   }
   ```

## Common failure modes

- `_package.files.readAsText('files/df.csv')` 404s — path is relative
  to the `files/` subtree (`DG-FACT-168`); pass `'df.csv'`.
- Asset outside `files/` — `_package.files` doesn't cover it; use
  `${_package.webRoot}<suffix>` with no leading slash (`DG-FACT-167`,
  `DG-FACT-171`).
- Binary mojibake on `readAsText` — use `readAsBytes` /
  `readBinaryDataFrames` (`DG-FACT-170`).
- `<table>_csv_options.json` sidecar has no effect — convention is
  stale (`DG-FACT-DRIFT-067`); pass options per-call.
- `DG.FileSource` is deprecated — use `DG.FilesDataSource`
  (`DG-FACT-172`).

## See also

- Source articles:
  - `help/develop/how-to/packages/work-with-package-files.md`
  - `help/develop/develop.md#package-structure` (linked from the article)
  - `help/develop/how-to/db/access-data.md` (linked from the article)
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-165` … `DG-FACT-172` and drifts `DG-FACT-DRIFT-066` …
  `DG-FACT-DRIFT-068`.
- Reference packages:
  `packages/HitTriage/src/app/utils.ts:159-178`,
  `packages/Bio/src/tests/utils.ts:10-19`,
  `packages/DrugBank/src/package.ts:10-21`,
  `packages/HitTriage/src/app/pep-triage-views/info-view.ts:22,32`.
- Related skills: `file-handlers` (register a parser for a file
  extension users open), `create-custom-file-viewers` (preview a file
  as a custom view), `access-data` (parent topic).
