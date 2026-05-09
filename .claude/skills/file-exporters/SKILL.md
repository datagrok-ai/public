---
name: file-exporters
description: Register a file-exporter function so a Datagrok package adds a "Save as <FORMAT>" entry to the platform's Export menu
---

# file-exporters

## When to use

Your package owns a target file format (`.sdf`, `.fasta`, `.parquet`, an
in-house text dump, …) and you want users to download the currently open
table in that format from the top **Export** menu. Triggers: "save the
table as SDF", "let me export to Parquet", "add an export to our `.foo`
format". For *importing* files, use a file-handler / file-viewer skill
instead — exporters only push data out.

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `datagrok-api` available (`import * as grok / DG`); the article's
  snippet omits imports and won't compile as written
  (drift `DG-FACT-DRIFT-028`).
- Familiarity with `DG.DataFrame` and `grok.shell.t` (`DG-FACT-081`).

## Steps

1. **Pick a registration form (decorator vs. header comments).**
   The platform discovers exporters via the function role `fileExporter`
   (`DG-FACT-079`). Two equivalent surfaces:
   - **Decorator (canonical, used by every recent package):** put a
     `static` method on `PackageFunctions` decorated with
     `@grok.decorators.fileExporter({description: '<menu label>'})`.
   - **Header comments (article form):** annotate an exported function
     with `//meta.role: fileExporter` and `//description: <menu label>`.
     This form survives only in generated `package.g.ts` files
     (drift `DG-FACT-DRIFT-028`).
   Either way the regenerated `package.g.ts` wrapper exposes
   `static <method>(): void | Promise<void>` — the function takes NO
   inputs and returns NO outputs (`DG-FACT-079`).

2. **Pick the menu-label convention (`description` is required).**
   Use `'As <FORMAT>...'` (with the trailing ellipsis) when the exporter
   pops a follow-up dialog before downloading; use
   `'Save as <FORMAT>'` when it downloads directly. Examples:
   `Chem` → `'As SDF...'`, `Bio` → `'As FASTA...'`,
   `Arrow` → `'Save as Parquet'` / `'Save as Feather'`
   (`DG-FACT-080`). The string is a HARD requirement: `grok check`
   errors with "File exporters should have a description parameter"
   when it is empty (drift `DG-FACT-DRIFT-030`).

3. **Add the exporter (direct-download form).**
   Read `grok.shell.t`, guard for `null`, serialize, and hand the bytes
   to `DG.Utils.download(filename, content, contentType?)` — the
   canonical helper (`DG-FACT-081`). Use this form for binary payloads
   (Parquet, Feather, XLSX): the article's `data:text/plain;...,encodeURIComponent(...)`
   pattern silently corrupts non-text content beyond ~2 MB on most
   browsers (drift `DG-FACT-DRIFT-029`).
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG   from 'datagrok-api/dg';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.fileExporter({description: 'Save as Parquet'})
     static saveAsParquet(): void {
       const table = grok.shell.t;
       if (table == null) return;            // no open table -> no-op
       const bytes = toParquetBytes(table);  // your serializer
       DG.Utils.download(`${table.name}.parquet`, bytes);
     }
   }
   ```
   Expected: build succeeds; `src/package.g.ts` gains a wrapper with
   `//description: Save as Parquet` + `//meta.role: fileExporter`
   (compare `packages/Chem/src/package.g.ts:166-170`).

4. **Or use the dialog form when the user must choose options first.**
   Declare `async` and delegate to a dialog helper that builds the
   payload and ultimately calls `DG.Utils.download(...)`. This is the
   `Chem.saveAsSdf` / `Bio.saveAsFasta` pattern
   (`DG-FACT-080`,`DG-FACT-DRIFT-027`).
   ```typescript
   @grok.decorators.fileExporter({description: 'As SDF...'})
   static async saveAsSdf(): Promise<void> {
     saveAsSdfDialog();   // your DG.Dialog that ends in DG.Utils.download
   }
   ```
   Expected: clicking the menu opens your dialog; the download fires
   only after the user clicks OK.

5. **Do NOT add an extension annotation.**
   Unlike file viewers / handlers, exporters have NO
   `meta.fileExporter: <ext>` companion key and no `FileExporterOptions`
   type (`DG-FACT-082`). The decorator's options bag is the generic
   `FunctionOptions` — `description`, `name`, `tags`, `meta.icon` are
   carried through, but extension routing is decided at runtime via the
   filename you pass to `DG.Utils.download(...)`.

6. **Build, publish, and let the platform auto-register.**
   ```bash
   npm install
   grok check                     # exits 0; description must be non-empty
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the function role
   `fileExporter` is the registration; the platform attaches the
   function to the global **Export** menu at startup (`DG-FACT-079`).

## Common failure modes

- **Menu entry never appears.** The function role is wrong or missing.
  Inspect the regenerated `src/package.g.ts`: it MUST contain
  `//meta.role: fileExporter` (`DG-FACT-079`). The role token is
  case-sensitive — `fileexporter` / `FileExporter` won't register.
- **`grok check` fails with "File exporters should have a description
  parameter".** The decorator was called as
  `@grok.decorators.fileExporter({})` or `description: ''`. The CLI
  linter enforces a non-empty `description` (drift `DG-FACT-DRIFT-030`,
  `DG-FACT-080`).
- **Download is empty / contains "[object Object]" / corrupts a binary
  file.** You wrote the article's hand-rolled
  `data:text/plain;charset=utf-8,encodeURIComponent(...)` form. Use
  `DG.Utils.download(filename, bytes)` instead — pass a `Uint8Array` for
  binary, a `string` for text (drift `DG-FACT-DRIFT-029`,
  `DG-FACT-081`).
- **Runtime `Cannot read properties of null (reading 'columns')`.**
  `grok.shell.t` returned `null` because no table was open. Add the
  `if (table == null) return;` guard at the top of every exporter
  (`DG-FACT-081`).
- **Article snippet copy-pasted, won't compile.** The article shows a
  JS-class method body with `var`, `let table = grok.shell.t;` and no
  imports, and links to a code path
  (`packages/Chem/package.js@73356b9c…`) that no longer exists
  (drifts `DG-FACT-DRIFT-027`,`DG-FACT-DRIFT-028`). Translate to the
  decorator form on `PackageFunctions` with explicit
  `import * as grok from 'datagrok-api/grok';` and
  `import * as DG from 'datagrok-api/dg';`.
- **Trying to scope the exporter to one file extension via metadata.**
  There is no `meta.fileExporter: <ext>` key (`DG-FACT-082`). Every
  registered exporter shows up in the global Export menu — gate by
  inspecting `grok.shell.t` inside the function body (e.g. require a
  `Molecule` semType column) and `return` early when the data shape is
  wrong, mirroring `Chem.saveAsSdf`'s structure-column check.

## Verification

- `grok check` exits `0`; `grok publish <host>` exits `0`.
- The regenerated `src/package.g.ts` contains an exported wrapper with
  `//description: <your label>` AND `//meta.role: fileExporter`.
- In Datagrok, open any table → top menu **Export** lists your entry
  with the exact `description` text. Click it: the browser downloads a
  file whose name matches the `filename` argument you passed to
  `DG.Utils.download(...)`.
- Close the table and re-trigger the menu entry: the function returns
  silently (no exception, no empty file) thanks to the `null` guard.

## See also

- Source articles:
  - `help/develop/how-to/files/file-exporters.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-079` … `DG-FACT-082` and drifts
  `DG-FACT-DRIFT-027..030`.
- Reference packages:
  - `packages/Arrow/src/package.ts:64-76` — decorator + `DG.Utils.download`
    for binary payloads (Parquet, Feather).
  - `packages/Chem/src/package.ts:526-531` — async exporter that
    delegates to a dialog (`As SDF...`).
  - `packages/Bio/src/package.ts:1404-1407` — sync exporter that
    delegates to a dialog (`As FASTA...`).
- Related skills:
  - `create-custom-file-viewers` (sibling — preview a file the user
    clicked, vs. push the open table out as a file).
