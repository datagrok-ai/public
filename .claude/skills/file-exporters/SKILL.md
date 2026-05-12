---
name: file-exporters
version: 0.1.0
description: |
  Add a "Save as <FORMAT>" entry to Datagrok's top Export menu so users
  can download the currently open table in a format your package owns
  (`.sdf`, `.fasta`, `.parquet`, an in-house dump). For plugin authors
  whose users keep asking to dump a table to a domain-specific file.
  Produces a `static` method on `PackageFunctions` decorated with
  `@grok.decorators.fileExporter`, plus the regenerated `package.g.ts`
  wrapper the platform auto-registers at startup.
  Use when asked to "let users download this table in our format",
  "hook into the platform's save/download menu", or "add a save-as
  option for a domain-specific format".
triggers:
  - hook into the save menu
  - download the open table as a custom format
  - add a save-as option to the platform
  - emit a domain-specific dump from a table
  - register a downloader for our format
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
harness-authored: true
---

# file-exporters

## When to use

Your package owns a target file format and you want users to download
the currently open table in that format from the top **Export** menu.
For *importing* files (custom previewers, drop-in handlers), use a
`create-custom-file-viewers` or file-handler skill — exporters only push
data out.

## Prerequisites

- A package scaffold (`grok create <PackageName>`); all paths below are
  relative to the package root.
- `datagrok-api` available with `import * as grok from 'datagrok-api/grok'`
  and `import * as DG from 'datagrok-api/dg'` — the article snippet omits
  imports and won't compile as written.
- Familiarity with `DG.DataFrame` and `grok.shell.t` (knowledge
  `DG-FACT-081`).

## Steps

1. **Add a `static` method on `PackageFunctions` and decorate it.**
   The platform discovers exporters via the function role `fileExporter`
   — exactly that camelCase token (knowledge `DG-FACT-079`). The
   decorator's options bag is the generic `FunctionOptions`; the
   `description` field becomes the visible menu label and is REQUIRED —
   `grok check` exits non-zero when it is missing or empty (knowledge
   `DG-FACT-080`). See the "Label convention" note below for phrasing.
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.fileExporter({description: 'Save as <FORMAT>'})
     static saveAsFormat(): void {
       /* body — implemented in step 2 */
     }
   }
   ```
   Expected: after `npm run build` (or any tool that regenerates
   `package.g.ts`), `src/package.g.ts` contains an exported wrapper
   carrying `//description: Save as <FORMAT>` and `//meta.role:
   fileExporter` — same shape as `packages/Arrow/src/package.g.ts:60-64`.

2. **Implement the body: read the table, serialize, hand off to `DG.Utils.download`.**
   `grok.shell.t` returns the currently open `DG.DataFrame` or `null` —
   guard for `null` before touching `.columns` or `.rowCount`. Then call
   `DG.Utils.download(filename, content, contentType?)`
   (`js-api/src/utils.ts:230-237`) — the canonical helper that wraps
   `Blob` + `URL.createObjectURL` + an anchor click. `content: BlobPart`
   accepts `string`, `Uint8Array`, `ArrayBuffer`, or `Blob`; binary
   payloads (Parquet, Feather, XLSX) need a `Uint8Array` (knowledge
   `DG-FACT-081`). Under modern TS lib types, a bare `Uint8Array` is
   typed `Uint8Array<ArrayBufferLike>` and is NOT assignable to
   `BlobPart` (TS2345) — narrow the generic at the call site with
   `as Uint8Array<ArrayBuffer>`, in the same expression that
   null-coalesces an empty buffer (knowledge `DG-FACT-426`,
   `packages/Arrow/src/package.ts:68,75`).
   ```typescript
   @grok.decorators.fileExporter({description: 'Save as Parquet'})
   static saveAsParquet(): void {
     const table = grok.shell.t;
     if (table == null) return;                  // no open table -> no-op
     const bytes = toParquetBytes(table);        // your serializer
     DG.Utils.download(
       `${table.name}.parquet`,
       (bytes ?? new Uint8Array(0)) as Uint8Array<ArrayBuffer>,
     );
   }
   ```
   Expected: clicking the menu entry triggers a browser download of
   `<table-name>.parquet` whose payload deserializes correctly in the
   target tool.

3. **(Optional) Use a dialog form when the user must choose options first.**
   When the exporter needs follow-up input (which columns to include,
   width/precision, filename), declare `async` and delegate to a
   `DG.Dialog` helper that builds the payload and ultimately calls
   `DG.Utils.download(...)` itself. This is the `Chem.saveAsSdf` /
   `Bio.saveAsFasta` pattern (`packages/Chem/src/package.ts:526-531`,
   `packages/Bio/src/package.ts:1404-1407`).
   ```typescript
   @grok.decorators.fileExporter({description: 'As SDF...'})
   static async saveAsSdf(): Promise<void> {
     saveAsSdfDialog();   // your DG.Dialog ends in DG.Utils.download
   }
   ```
   Expected: clicking the menu opens your dialog; the download fires
   only after the user clicks OK.

4. **Build, publish, let the platform auto-register.**
   ```bash
   npm install
   grok check                     # exits 0; description must be non-empty
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the `fileExporter` role
   is the registration; the platform attaches the function to the
   global **Export** menu at package startup (knowledge `DG-FACT-079`).

## Notes

- **Label convention.** Pick the `description` phrasing to match
  runtime behavior (knowledge `DG-FACT-080`): use `'As <FORMAT>...'`
  with a trailing ellipsis when a dialog opens first (Chem
  `'As SDF...'`, Bio `'As FASTA...'`); use `'Save as <FORMAT>'` when
  the click downloads directly (Arrow `'Save as Parquet'`,
  `'Save as Feather'`).
- **No extension annotation, no per-format scoping.** Unlike file
  viewers/handlers, exporters have NO `meta.fileExporter: <ext>` key
  and no `FileExporterOptions` type — the options bag is the generic
  `FunctionOptions` (`js-api/src/decorators/functions.ts:350-356`,
  knowledge `DG-FACT-082`). Every registered exporter appears in the
  global Export menu. If your format only makes sense for a particular
  data shape (e.g. a `Molecule` semType column), gate inside the body
  and `return` early — mirror `Chem.saveAsSdf`'s structure-column null
  check.

## Common failure modes

- **Menu entry never appears.** The function role is wrong or missing.
  Inspect the regenerated `src/package.g.ts`: it MUST contain
  `//meta.role: fileExporter` (camelCase, knowledge `DG-FACT-079`).
  `fileexporter` / `FileExporter` won't register.
- **`grok check` fails with "File exporters should have a description
  parameter".** The decorator was called as
  `@grok.decorators.fileExporter({})` or `description: ''`. The CLI
  linter at `tools/bin/commands/check.ts:358-359` enforces a non-empty
  `description` (knowledge `DG-FACT-080`).
- **Download is empty / shows `[object Object]` / corrupts a binary
  file.** You wrote a hand-rolled
  `<a href="data:text/plain;...,encodeURIComponent(...)">` form. That
  pattern silently corrupts non-text content beyond ~2 MB on most
  browsers. Use `DG.Utils.download(filename, bytes)` — pass a
  `Uint8Array` for binary, a `string` for text (knowledge `DG-FACT-081`).
- **Runtime `Cannot read properties of null (reading 'columns')`.**
  `grok.shell.t` returned `null` because no table was open. Add the
  `if (table == null) return;` guard at the top of every exporter.
- **Trying to scope the exporter to one extension via metadata.**
  There is no `meta.fileExporter: <ext>` key (knowledge `DG-FACT-082`).
  Every registered exporter appears in the global Export menu — gate by
  inspecting `grok.shell.t` inside the body and returning early when
  the data shape is wrong.
- **`grok check && webpack` fails with TS2345 "Argument of type
  'Uint8Array<ArrayBufferLike>' is not assignable to parameter of type
  'BlobPart'".** Modern TS `lib.dom.d.ts` types `Uint8Array` with an
  `ArrayBufferLike` slot that includes `SharedArrayBuffer`, which is
  not a valid `BlobPart`. Narrow at the call site:
  `(bytes ?? new Uint8Array(0)) as Uint8Array<ArrayBuffer>` — the
  pattern both Arrow exporters use (knowledge `DG-FACT-426`).

## Verification

- `grok check` exits `0`; `grok publish <host>` exits `0`.
- `src/package.g.ts` contains an exported wrapper with
  `//description: <your label>` AND `//meta.role: fileExporter`
  (compare `packages/Arrow/src/package.g.ts:60-70`).
- In Datagrok, open any table → top menu **Export** lists your entry
  with the exact `description` text. Click it: the browser downloads a
  file whose name matches the `filename` argument you passed to
  `DG.Utils.download(...)`.
- Close the table and re-trigger the menu entry: the function returns
  silently (no exception, no empty download) thanks to the `null` guard.

## See also

- Source articles: `help/develop/how-to/files/file-exporters.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-079` (role token + signature), `DG-FACT-080` (description
  required + label convention), `DG-FACT-081` (`grok.shell.t` +
  `DG.Utils.download`), `DG-FACT-082` (no extension metadata),
  `DG-FACT-426` (`Uint8Array<ArrayBuffer>` cast for `BlobPart`).
- Reference packages:
  - `packages/Arrow/src/package.ts:64-76` — decorator + `DG.Utils.download`
    for binary payloads.
  - `packages/Chem/src/package.ts:526-531` — async exporter delegating
    to a dialog (`As SDF...`).
  - `packages/Bio/src/package.ts:1404-1407` — sync exporter delegating
    to a dialog (`As FASTA...`).
- Related skills: `create-custom-file-viewers` (inverse direction —
  preview a file the user clicked).
