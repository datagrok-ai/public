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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# file-exporters

## When to use

Your package owns a target file format and you want users to download
the currently open table in that format from the top **Export** menu.
For *importing* files (custom previewers, drop-in handlers), use a
`create-custom-file-viewers` or file-handler skill — exporters only push
data out.

## Prerequisites

- Familiarity with `DG.DataFrame` and `grok.shell.t` (knowledge
  `DG-FACT-081`).

## Steps

1. **Add a `static` method on `PackageFunctions` and decorate it.**
   Role token is `fileExporter` (camelCase — `DG-FACT-079`).
   `description` is REQUIRED and becomes the menu label — empty fails
   `grok check` (`DG-FACT-080`).
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

2. **Implement the body: read the table, serialize, hand off to `DG.Utils.download`.**
   `grok.shell.t` returns the currently open `DG.DataFrame` or `null`.
   `DG.Utils.download(filename, content, contentType?)` accepts
   `string`, `Uint8Array`, `ArrayBuffer`, or `Blob` (`DG-FACT-081`).
   For binary payloads under modern TS, narrow the type with
   `as Uint8Array<ArrayBuffer>` (`DG-FACT-426`).
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

3. **(Optional) Use a dialog form when the user must choose options first.**
   Declare `async` and delegate to a `DG.Dialog` helper that ultimately
   calls `DG.Utils.download(...)`:
   ```typescript
   @grok.decorators.fileExporter({description: 'As SDF...'})
   static async saveAsSdf(): Promise<void> {
     saveAsSdfDialog();   // your DG.Dialog ends in DG.Utils.download
   }
   ```

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

- **Label convention.** `'As <FORMAT>...'` (with ellipsis) when a
  dialog opens first; `'Save as <FORMAT>'` when the click downloads
  directly (`DG-FACT-080`).
- **No per-format scoping.** No `meta.fileExporter: <ext>` key
  (`DG-FACT-082`). Gate inside the body if your format only makes
  sense for certain data shapes.

## Common failure modes

- **Menu entry never appears.** `//meta.role: fileExporter` missing
  in `package.g.ts` (`DG-FACT-079`).
- **`grok check`: "File exporters should have a description
  parameter".** `description` empty or missing (`DG-FACT-080`).
- **Binary download corrupts.** Use `DG.Utils.download(filename,
  bytes)`, not a hand-rolled `data:` URL (`DG-FACT-081`).
- **Runtime null on `.columns`.** Guard `if (table == null) return;`.
- **TS2345 `Uint8Array<ArrayBufferLike>` not assignable to
  `BlobPart`.** Cast: `(bytes ?? new Uint8Array(0)) as
  Uint8Array<ArrayBuffer>` (`DG-FACT-426`).

## See also

- Source: `help/develop/how-to/files/file-exporters.md`.
- Knowledge: `DG-FACT-079`–`082`, `DG-FACT-426`.
- Related skills: `create-custom-file-viewers`.
