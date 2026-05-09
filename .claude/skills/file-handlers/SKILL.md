---
name: file-handlers
description: Register a file-handler function so a Datagrok package imports a custom file format and opens its contents as one or more tables
---

# file-handlers

## When to use

Your package owns one or more file extensions (`.fasta`, `.sdf`, `.bam`,
`.parquet`, an in-house text dump, …) and you want the platform to invoke
your parser whenever a user opens a file with that extension, returning
`DG.DataFrame[]` for the platform to display as TableViews. Triggers:
"open `.fasta` as a table", "import our custom binary format", "register
a parser for `.foo`". For *previewing* a file without importing as a
table (3D structures, NMR spectra), use `create-custom-file-viewers`
instead. For *exporting* the open table back out, use `file-exporters`.

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `datagrok-api` available (`import * as grok / DG`); the article omits
  imports and won't compile as written (drift `DG-FACT-DRIFT-032`).
- Familiarity with `DG.DataFrame` construction from
  `grok.data.parseCsv(...)` or your own parser.

## Steps

1. **Pick a registration form (decorator vs. header comments).**
   The platform discovers handlers via the function role `fileHandler`
   (`DG-FACT-083`). The role token is camelCase in every production
   package and in the article; the kebab form `file-handler` exists in
   `js-api/src/const.ts` but is not what gets emitted (drift
   `DG-FACT-DRIFT-031`). Two equivalent surfaces:
   - **Decorator (canonical, used by every recent package):** put a
     `static` method on `PackageFunctions` decorated with
     `@grok.decorators.fileHandler({ext: '<csv-of-extensions>', description: '<label>'})`.
   - **Header comments (article form):** annotate an exported function
     with `//meta.role: fileHandler` and `//meta.ext: <csv-of-extensions>`.
     This form survives only in generated `package.g.ts` files
     (drift `DG-FACT-DRIFT-032`).

2. **Declare the extensions in `meta.ext` as a comma-separated list.**
   One handler entry can claim many extensions:
   `meta.ext: fasta, fna, ffn, faa, frn, fa, fst` (Bio FASTA),
   `meta.ext: sdf,mol` (Chem SDF), `meta.ext: bam, bai` (Bio BAM)
   (`DG-FACT-084`). Whitespace around commas is tolerated. The `ext`
   field is REQUIRED on `FileHandlerOptions` — there is no default
   (`DG-FACT-085`). When per-extension bodies differ, write one
   decorator per extension (the nmrium `.jdx` / `.dx` / `.nmrium`
   pattern at `packages/nmrium/src/package.ts:20-46`).

3. **Pick the input type — `string` for text, `Uint8Array` for binary.**
   The codegen default is `inputs: [{name: 'content', type: 'string'}]`
   (`DG-FACT-087`). For binary formats override with `@grok.decorators.param({type: 'list'})`
   on a `Uint8Array` parameter — silently using the string default on a
   binary file corrupts the payload (drift `DG-FACT-DRIFT-033`).
   ```typescript
   // src/package.ts — text format (FASTA)
   import * as grok from 'datagrok-api/grok';
   import * as DG   from 'datagrok-api/dg';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.fileHandler({
       ext: 'fasta, fna, ffn, faa, frn, fa, fst',
       description: 'Opens FASTA file',
     })
     static importFasta(content: string): DG.DataFrame[] {
       return parseFasta(content);   // your parser → DG.DataFrame[]
     }
   }
   ```
   Expected: build succeeds; `src/package.g.ts` gains a wrapper with
   `//input: string content`, `//output: list<dataframe> result`,
   `//meta.role: fileHandler`, `//meta.ext: fasta, fna, ...`
   (compare `packages/Bio/src/package.g.ts:486-503`).

4. **Override to bytes for binary formats.**
   ```typescript
   @grok.decorators.fileHandler({description: 'Opens SDF file', ext: 'sdf,mol'})
   static importSdf(
     @grok.decorators.param({type: 'list'}) bytes: Uint8Array
   ): DG.DataFrame[] | void {
     try { return _importSdf(Uint8Array.from(bytes)); }
     catch (e: any) { grok.shell.warning('file is not supported or malformed'); grok.shell.error(e); }
   }
   ```
   Expected: `package.g.ts` emits `//input: list bytes` (NOT
   `//input: string content`) — confirms the override took effect
   (compare `packages/Chem/src/package.g.ts:793-799`,
   `DG-FACT-087`,`DG-FACT-DRIFT-033`).

5. **Disambiguate ambiguous extensions with `fileViewerCheck`.**
   Several packages legitimately claim the same extension (`.json`,
   `.txt`, `.xml`, `.dx`). Set `fileViewerCheck: '<Pkg>:<Fn>'` on the
   options bag — the platform calls the named check function with the
   file content and routes to this handler only on truthy return
   (`DG-FACT-085`, drift `DG-FACT-DRIFT-034`).
   ```typescript
   @grok.decorators.fileHandler({
     fileViewerCheck: 'Nmrium:checkNmriumJdx',
     ext: 'jdx',
     outputs: [{name: 'tables', type: 'list'}],
   })
   static async jdxFileHandler(bytes: string) {
     return await addNmriumView('jdx', bytes);   // returns []; view added as side effect
   }
   ```
   Expected: opening a `.jdx` file that fails `checkNmriumJdx` is
   handled by another package; one that passes is routed here.

6. **Returning `[]` is legal — handler may open a custom view as a side effect.**
   The "side-effect view" pattern is canonical: call
   `grok.shell.addView(...)` from inside the handler and return `[]`.
   nmrium does this so the platform opens no extra empty TableView on
   top of the custom NMR view (`DG-FACT-087`, drift `DG-FACT-DRIFT-035`,
   `packages/nmrium/src/package.ts:25-46`).

7. **Build, publish, and let the platform auto-register.**
   ```bash
   npm install
   grok check                     # exits 0
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the function role
   `fileHandler` IS the registration; opening any matching file in
   Datagrok dispatches into your function (`DG-FACT-083`).

## Common failure modes

- **Handler never fires; the platform falls back to the built-in
  importer.** The role token is wrong or missing. Inspect
  `src/package.g.ts`: it MUST contain `//meta.role: fileHandler`
  (camelCase) and a `//meta.ext: <ext>` line for each claimed extension
  (`DG-FACT-083`,`DG-FACT-084`). The token is case-sensitive —
  `filehandler` / `FileHandler` won't register; the kebab variant
  `file-handler` is not what packages emit (drift `DG-FACT-DRIFT-031`).
- **Binary file is mangled to mojibake.** You took the string default
  on a binary format. Override with
  `@grok.decorators.param({type: 'list'}) bytes: Uint8Array` (decorator
  form) or `//input: list bytes` (header form), and confirm
  `src/package.g.ts` shows `//input: list bytes` (drift
  `DG-FACT-DRIFT-033`,`DG-FACT-087`).
- **Two packages fight over the same extension.** Both register a
  handler for `.json` / `.txt` / `.dx` and the platform picks
  unpredictably. Add `fileViewerCheck: '<Pkg>:<Fn>'` to your
  `FileHandlerOptions` so the platform routes by content (drift
  `DG-FACT-DRIFT-034`,`DG-FACT-085`).
- **TypeScript error: `Property 'ext' is missing`.** `ext` is REQUIRED
  on `FileHandlerOptions` (no `?`); the decorator's `config` argument
  is itself non-optional unlike `fileViewer(config?)` (`DG-FACT-085`).
  Always pass at least `{ext: '<csv>'}`.
- **Article snippet copy-pasted, won't compile.** The article shows a
  bare JS function with header comments and no imports, and uses
  `return tables;` against an undefined variable (drift
  `DG-FACT-DRIFT-032`,`DG-FACT-DRIFT-035`). Translate to the decorator
  form on `PackageFunctions` with explicit
  `import * as grok from 'datagrok-api/grok';` and
  `import * as DG from 'datagrok-api/dg';`, and build a real
  `DG.DataFrame[]` from your parser before returning.

## Verification

- `grok check` exits `0`; `grok publish <host>` exits `0`.
- The regenerated `src/package.g.ts` contains, for each claimed
  extension, a wrapper with `//meta.role: fileHandler`,
  `//meta.ext: <your csv>`, and the right input line — `//input: string content`
  (text) or `//input: list bytes` (binary).
- In Datagrok, drag-drop or open a file with one of the registered
  extensions: the platform invokes your handler, and the returned
  `DG.DataFrame[]` opens as one TableView per dataframe. (For the
  side-effect-view / `return []` pattern, the custom view appears
  instead of a TableView.)
- Open a file whose extension is shared with another package and the
  `fileViewerCheck` returns falsy: another package handles it; truthy
  returns route here.

## See also

- Source articles:
  - `help/develop/how-to/files/file-handlers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-083` … `DG-FACT-087` and drifts
  `DG-FACT-DRIFT-031..035`.
- Reference packages:
  - `packages/Chem/src/package.ts:1908-1934` — decorator + bytes
    override for `.sdf`, `.mol`, `.smi`.
  - `packages/Bio/src/package.ts:1089-1111` — multi-extension `ext`
    list for `.fasta` family; bytes form for `.bam`, `.bai`.
  - `packages/nmrium/src/package.ts:19-46` — `fileViewerCheck` to
    disambiguate `.jdx`/`.dx` and the `return []` side-effect-view
    pattern.
- Related skills:
  - `create-custom-file-viewers` (sibling — preview the file as a
    custom view without importing as a table).
  - `file-exporters` (inverse — push the open table out as a file).
