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
harness-authored: true
---

# file-handlers

## When to use

Your package owns a tabular-ish file format and you want a user dropping
a `.fasta` / `.sdf` / in-house file into Datagrok to get one or more
`DataFrame` tables opened automatically. A handler IMPORTS — not a
preview (`create-custom-file-viewers`), not a writer (`file-exporters`).
The three "open this file" roles are distinct (knowledge `DG-FACT-086`).

## Prerequisites

- A package scaffold (`grok create <Name>`); paths are relative to the
  package root.
- `datagrok-api` imported as `import * as grok from 'datagrok-api/grok'`
  and `import * as DG from 'datagrok-api/dg'`.
- A working parser that returns `DG.DataFrame[]`.

## Steps

1. **Add a `static` method on `PackageFunctions` and decorate it.**
   The platform discovers handlers via the function role `fileHandler`
   — exactly that camelCase token in the emitted `package.g.ts`
   (knowledge `DG-FACT-083`). The decorator's options bag is
   `FileHandlerOptions { ext: string; fileViewerCheck?: string }` and
   `ext` is REQUIRED — `config` itself is non-optional, unlike
   `fileViewer(config?)` (knowledge `DG-FACT-085`). The handler signature
   is `(content: string | Uint8Array) => DG.DataFrame[]`.
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
   Expected: after `npm run build`, `src/package.g.ts` carries
   `//meta.role: fileHandler` and `//meta.ext: <your-list>` — same
   shape as `packages/Bio/src/package.g.ts:486-494`.

2. **Declare every extension you claim in `meta.ext`.**
   `ext` is a single comma-separated string; whitespace around commas
   is tolerated (knowledge `DG-FACT-084`). One handler entry, many
   extensions — `'fasta, fna, ffn, faa, frn, fa, fst'`
   (`packages/Bio/src/package.ts:1089-1094`), `'sdf,mol'`
   (`packages/Chem/src/package.ts:1908-1911`). Repeating `//meta.ext:`
   does NOT additively merge — the parser keeps only the last. If two
   extensions need different bodies, register two decorators — the
   per-extension dispatch pattern at
   `packages/nmrium/src/package.ts:20-46`.

3. **Pick the input type to match your payload — text or binary.**
   The default codegen contract is text:
   `inputs: [{name: 'content', type: 'string'}], outputs: [{name: 'tables', type: 'list'}]`
   (knowledge `DG-FACT-087`). For binary formats override the input
   shape with a `@grok.decorators.param({type: 'list'})` annotation on
   the parameter and type it as `Uint8Array` — the Chem SDF pattern at
   `packages/Chem/src/package.ts:1908-1920`. The regenerated wrapper
   then carries `//input: list bytes` (compare
   `packages/Chem/src/package.g.ts:793-799`).
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
   Expected: the emitted wrapper has `//input: list bytes` (NOT
   `//input: string content`) and your parser receives raw bytes.

4. **Return `DG.DataFrame[]`; the platform opens one TableView per
   element.** Empty `[]` is legal — nmrium returns `[]` after
   side-effect-opening its own view (`packages/nmrium/src/package.ts:12-17,25-46`).
   Use that pattern only when an accompanying `fileViewer` provides the
   visualization (knowledge `DG-FACT-086`).

5. **(Optional) Disambiguate ambiguous extensions with `fileViewerCheck`.**
   When the same extension is shared with another tool (e.g. `.jdx`,
   `.dx` JCAMP-DX), gate your handler on a content probe by setting
   `fileViewerCheck: '<Package>:<Function>'` — the platform calls the
   named function with file content (`string` for text,
   `Uint8Array` for binary) and routes only on truthy return
   (knowledge `DG-FACT-085`). The probe must itself be a Datagrok
   function — decorate it with plain `@grok.decorators.func()`. Pattern
   at `packages/nmrium/src/package.ts:20-40, 70-73`.
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
   No explicit `register(...)` call is needed — the `fileHandler` role
   plus `meta.ext` IS the registration; the platform invokes the
   function whenever a user opens a matching file (knowledge
   `DG-FACT-083`).

## Common failure modes

- **Nothing happens on file open.** The regenerated `src/package.g.ts`
  is missing `//meta.role: fileHandler` or `//meta.ext: <ext>`. Both
  lines MUST be present; the role token is camelCase as emitted
  (knowledge `DG-FACT-083`). Re-run `npm install` if it looks stale.
- **Build error: `Property 'ext' is missing` / `Expected 1 argument, but got 0`.**
  You wrote `@grok.decorators.fileHandler()` or `({})`. `config` and
  its `ext` field are both REQUIRED (knowledge `DG-FACT-085`).
- **Handler runs but the parsed binary is garbled.** You took the
  default text input (`content: string`) for a binary format. Add
  `@grok.decorators.param({type: 'list'})` to the parameter and type
  it `Uint8Array` (knowledge `DG-FACT-087`). The wrapper should then
  show `//input: list bytes`.
- **Only one of two `meta.ext:` lines applies.** Repeating the
  annotation does NOT additively merge — the parser keeps only the
  last entry (knowledge `DG-FACT-084`). Comma-separate into one entry
  or register one decorator per extension when bodies differ.
- **Two handlers fight over `.jdx`.** Both claim the same extension
  with no `fileViewerCheck`; the platform picks one non-deterministically.
  Add `fileViewerCheck: '<Pkg>:<Fn>'` to each (knowledge `DG-FACT-085`).

## Verification

- `grok check` exits `0`; `grok publish <host>` exits `0`.
- `src/package.g.ts` contains the exported wrapper with
  `//meta.role: fileHandler` and `//meta.ext: <your-list>` (compare
  `packages/Bio/src/package.g.ts:486-503` or
  `packages/Chem/src/package.g.ts:793-799`).
- In Datagrok, drag-drop or open a file with one of the registered
  extensions: a TableView appears for each returned `DG.DataFrame` in
  the order your handler returned them.
- Open an unrelated file (`.csv`, `.txt`): your handler is NOT invoked
  — the generic importer takes over.

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
