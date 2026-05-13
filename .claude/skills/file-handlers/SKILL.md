---
name: file-handlers
version: 0.1.0
description: |
  Claim a file extension (e.g. `.fasta`, `.sdf`, `.bam`) so that when a
  user opens such a file in Datagrok, your package parses it and the
  platform opens the resulting `DataFrame[]` as TableViews — instead of
  the generic text/binary fallback. For plugin authors whose users keep
  dragging in domain-specific files they want imported as tables.
  Produces a `static` method on `PackageFunctions` decorated with
  `@grok.decorators.fileHandler({ext})`, plus the regenerated
  `package.g.ts` wrapper the platform auto-registers at startup.
  Use when asked to "import a domain file as a table", "open our `.sdf`
  files as DataFrames", or "make Datagrok parse a custom format on
  drop".
triggers:
  - import a domain file as a table
  - parse a custom format on drop
  - open our sdf as dataframes
  - claim a file extension for import
  - turn a fasta drop into tables
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# file-handlers

## When to use

Your package owns a tabular-ish file format and you want a user dropping
a `.fasta` / `.sdf` / in-house file into Datagrok to get one or more
`DataFrame` tables opened automatically. A handler IMPORTS — not a
preview (`create-custom-file-viewers`), not a writer (`file-exporters`).
The three "open this file" roles are distinct (knowledge `DG-FACT-086`).

## Prerequisites

- A working parser that returns `DG.DataFrame[]`.

## Steps

1. **Add a `static` method on `PackageFunctions` and decorate it** with `@grok.decorators.fileHandler({ext})` — role `fileHandler` (camelCase), required `ext`, signature `(content: string | Uint8Array) => DG.DataFrame[]` (`DG-FACT-083`, `DG-FACT-085`).
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.fileHandler({
       description: 'Opens FASTA file',
       ext: 'fasta, fna, ffn, faa, frn, fa, fst',
     })
     static importFasta(fileContent: string): DG.DataFrame[] {
       return new FastaFileHandler(fileContent).importFasta();
     }
   }
   ```
   Expected: `src/package.g.ts` carries `//meta.role: fileHandler` and `//meta.ext: <your-list>`.

2. **Declare every extension you claim in `meta.ext`** — single comma-separated string, whitespace tolerated (`DG-FACT-084`). Repeating `//meta.ext:` keeps only the last entry; register one decorator per extension when bodies differ.

3. **Pick the input type to match your payload.** Text is the default; for binary, add `@grok.decorators.param({type: 'list'})` and type the parameter `Uint8Array` (`DG-FACT-087`).
   ```typescript
   @grok.decorators.fileHandler({description: 'Opens SDF file', ext: 'sdf,mol'})
   static importSdf(
     @grok.decorators.param({type: 'list'}) bytes: Uint8Array
   ): DG.DataFrame[] | void {
     try { return _importSdf(Uint8Array.from(bytes)); }
     catch (e: any) {
       grok.shell.warning('file is not supported or malformed');
       grok.shell.error(e);
     }
   }
   ```
   Expected: emitted wrapper has `//input: list bytes`.

4. **Return `DG.DataFrame[]`; the platform opens one TableView per element.** Returning `[]` is legal when a sibling `fileViewer` provides the visualization (`DG-FACT-086`).

5. **(Optional) Disambiguate ambiguous extensions with `fileViewerCheck`** (`DG-FACT-085`). The probe must itself be a registered function.
   ```typescript
   @grok.decorators.fileHandler({fileViewerCheck: 'Nmrium:checkNmriumJdx', ext: 'jdx'})
   static async jdxFileHandler(bytes: string) { return await addNmriumView('jdx', bytes); }

   @grok.decorators.func()
   static checkNmriumJdx(content: string): boolean { return content.includes('NMR SPECTRUM'); }
   ```

6. **Build, publish, let the platform auto-register.**
   ```bash
   npm install
   grok check
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` — role + `meta.ext` IS the registration (`DG-FACT-083`).

## Common failure modes

- **Nothing happens on file open** — regenerated `package.g.ts` is missing `//meta.role: fileHandler` or `//meta.ext` (`DG-FACT-083`).
- **Build error: `Property 'ext' is missing`** — `config` and its `ext` field are required (`DG-FACT-085`).
- **Binary parsed as text (garbled)** — add `@grok.decorators.param({type: 'list'})` + `Uint8Array` (`DG-FACT-087`).
- **Only one of two `meta.ext:` lines applies** — repeats don't merge; comma-separate or split decorators (`DG-FACT-084`).
- **Two handlers fight over the same extension** — add `fileViewerCheck` to each (`DG-FACT-085`).

## See also

- Source articles:
  - `help/develop/how-to/files/file-handlers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-083` (role + signature), `DG-FACT-084` (`meta.ext`
  comma-list), `DG-FACT-085` (`FileHandlerOptions` + `fileViewerCheck`),
  `DG-FACT-086` (handler vs viewer vs exporter), `DG-FACT-087` (codegen
  contract + binary input override).
- Reference packages:
  - `packages/Bio/src/package.ts:1089-1099` — text handler, comma-list.
  - `packages/Chem/src/package.ts:1908-1920` — binary handler with
    `@grok.decorators.param({type: 'list'})`.
  - `packages/nmrium/src/package.ts:20-46, 70-73` — per-extension
    decorators + `fileViewerCheck`.
- Related skills:
  - `create-custom-file-viewers` — visualize a file without importing
    it as a table.
  - `file-exporters` — write the active table out via the Export menu.
