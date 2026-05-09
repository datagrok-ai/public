---
name: folder-content-preview
description: Register a folderViewer function so a Datagrok package shows a custom widget (or view) when the user opens a matching folder in the file-share browser
---

# folder-content-preview

## When to use

Your package owns a *folder* shape — an SDTM study tree, a plates layout
directory, an HTRF run, a clinical-trial dump — and you want Datagrok's
file-share browser to render a custom panel/launcher when the user
clicks the folder, instead of the default file listing. Triggers:
"detect SDTM folders and show a 'Run study' button", "preview plate
folders as a plate map", "show a project-specific dashboard for any
folder containing `manifest.json`". For previewing a single *file* (not
the folder), use `create-custom-file-viewers`. For *importing* matched
files as tables, use `file-handlers`.

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `datagrok-api` available (`import * as grok / ui / DG`); the article's
  snippet omits imports and won't compile as written (drift `DG-FACT-DRIFT-036`).
- Familiarity with `DG.FileInfo`, `DG.Widget.fromRoot(...)`
  (`DG-FACT-092`,`DG-FACT-093`).

## Steps

1. **Pick a registration form (decorator vs. header comments).**
   The platform discovers folder previews via the function role
   `folderViewer` (camelCase — `DG-FACT-088`). Two equivalent surfaces:
   - **Decorator (canonical, used by every recent package):** put a
     `static` method on `PackageFunctions` decorated with
     `@grok.decorators.folderViewer({...})`. The `config` arg is
     OPTIONAL — there is no specialized `FolderViewerOptions`
     interface (`DG-FACT-090`).
   - **Header comments (article form):** annotate an exported function
     with `//meta.role: folderViewer`, `//input: file folder`,
     `//input: list<file> files`, `//output: widget` (`DG-FACT-088`,
     `DG-FACT-089`). This form survives only in generated `package.g.ts`
     (drift `DG-FACT-DRIFT-036`).

2. **Match the canonical signature.**
   The function-role descriptor declares
   `folderViewer(folder: File, files: list<file>): Widget`
   (`DG-FACT-089`). Codegen defaults to
   `inputs: [{name:'folder',type:'file'}, {name:'files',type:'list<file>'}]`
   and `outputs: [{name:'result',type:'widget'}]`
   (`tools/bin/utils/func-generation.ts:420-426`). Use
   `folder: DG.FileInfo` and `files: DG.FileInfo[]` in TypeScript
   (`DG-FACT-093`).

3. **Detect the folder shape and return `undefined` to opt out.**
   Decide cheaply (filename scan or `folder.name`) whether your custom
   preview applies. Returning `undefined` (or omitting the return)
   signals "no preview" — the platform falls back to the default
   listing (`DG-FACT-091`). Inspect either input: scan `files` for a
   sentinel filename (ClinicalCase pattern), or inspect `folder.name`
   (Plates pattern).
   ```typescript
   // src/package.ts — sentinel-file pattern (SDTM Demographics)
   import * as grok from 'datagrok-api/grok';
   import * as ui   from 'datagrok-api/ui';
   import * as DG   from 'datagrok-api/dg';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.folderViewer({})
     static async clinicalCaseFolderLauncher(
       folder: DG.FileInfo,
       files: DG.FileInfo[],
     ): Promise<DG.Widget | undefined> {
       if (!files.some((f) => f.fileName.toLowerCase() === 'dm.csv'))
         return undefined;                         // opt out — fall back
       const csv = await grok.dapi.files.readAsText(`${folder.fullPath}/dm.csv`);
       const studyId = DG.DataFrame.fromCsv(csv).get('STUDYID', 0) ?? 'unknown';
       return DG.Widget.fromRoot(ui.div([
         ui.divText(`SDTM study: ${studyId}`),
         ui.button('Run study', () => grok.shell.info(`Launching ${studyId}`)),
       ]));
     }
   }
   ```
   Expected: build succeeds; `src/package.g.ts` gains a wrapper with
   `//meta.role: folderViewer`, `//input: file folder`,
   `//input: list<file> files`, `//output: widget result`
   (compare `packages/ClinicalCase/src/package.g.ts:37-43`).

4. **Build the widget — wrap any `HTMLElement` with `DG.Widget.fromRoot`.**
   `DG.Widget.fromRoot(root: HTMLElement): Widget` is the canonical
   adapter (`DG-FACT-092`, `js-api/src/widgets/base.ts:374-376`). Build
   with `ui.div`, `ui.panel`, `ui.button` and pass the root element.
   The article's stub (`ui.button('START', () => grok.shell.info('Foo'))`)
   compiles but is not useful in production — read the matched files
   and wire a real launcher (drift `DG-FACT-DRIFT-037`).

5. **Return a `DG.ViewBase` instead of a widget — opt-in via `outputs: dynamic`.**
   The default codegen output is `widget`. To return `DG.ViewBase`
   (e.g. an entire view, not just a panel), override the output type
   with `outputs: [{name: 'result', type: 'dynamic'}]` on the
   decorator (`DG-FACT-094`); the regenerated wrapper then emits
   `//output: dynamic result`. Without this override, returning a
   `ViewBase` is a type mismatch at codegen (drift `DG-FACT-DRIFT-038`).
   ```typescript
   @grok.decorators.folderViewer({outputs: [{name: 'result', type: 'dynamic'}]})
   static async platesFolderPreview(
     folder: DG.FileInfo, files: DG.FileInfo[],
   ): Promise<DG.Widget | DG.ViewBase | undefined> {
     if (!folder.name?.toLowerCase().includes('plate')) return undefined;
     return getPlatesFolderPreview(files);                  // your view builder
   }
   ```
   Expected: `package.g.ts` shows `//output: dynamic result` (compare
   `packages/Plates/src/package.g.ts:17-23`).

6. **Build, publish, and let the platform auto-register.**
   ```bash
   npm install
   grok check                     # exits 0
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the function role
   `folderViewer` is the registration; opening any folder dispatches
   into every registered viewer until one returns a non-`undefined`
   result (`DG-FACT-088`, `DG-FACT-091`).

## Common failure modes

- **Preview never appears; the default folder listing shows.** Header
  is missing `//meta.role: folderViewer` (camelCase). Inspect the
  regenerated `src/package.g.ts`: it MUST contain
  `//meta.role: folderViewer` AND the `folder` + `files` inputs
  (`DG-FACT-088`, `DG-FACT-089`). The token is case-sensitive —
  `folder-viewer` / `FolderViewer` won't register.
- **Preview applies to every folder, even unrelated ones.** Your
  detector matched too broadly. Always start the body with a sentinel
  check (`files.some(f => f.fileName === '<expected>')` or a folder
  name predicate) and `return undefined` on miss (`DG-FACT-091`).
  Multiple packages each register a `folderViewer` and the platform
  asks each — sloppy detection collides with siblings.
- **Article snippet copy-pasted, won't compile.** The article uses the
  bare-function/header form with no imports and a `: DG.Widget | undefined`
  return on a non-imported `DG` (drift `DG-FACT-DRIFT-036`). Translate
  to the decorator form on `PackageFunctions` with explicit
  `import * as grok from 'datagrok-api/grok';`,
  `import * as ui from 'datagrok-api/ui';`,
  `import * as DG from 'datagrok-api/dg';`.
- **Type error: `ViewBase is not assignable to Widget`.** You returned
  a `DG.ViewBase` without overriding the output type. Add
  `outputs: [{name: 'result', type: 'dynamic'}]` to the decorator
  config (`DG-FACT-094`, drift `DG-FACT-DRIFT-038`); the codegen will
  then emit `//output: dynamic result`.
- **Stub demo widget shipped to production.** The article's `'START'`
  button is a placeholder. Production code (`packages/ClinicalCase`)
  reads the sentinel file, parses it, and calls a real app entry
  (`grok.functions.call('<Pkg>:<App>')`) — not `grok.shell.info('Foo')`
  (drift `DG-FACT-DRIFT-037`).

## Verification

- `grok check` exits `0`; `grok publish <host>` exits `0`.
- The regenerated `src/package.g.ts` contains a wrapper with
  `//meta.role: folderViewer`, `//input: file folder`,
  `//input: list<file> files`, and `//output: widget result` (or
  `//output: dynamic result` when you overrode `outputs`).
- In Datagrok, **Browse → Files**, navigate to a folder that satisfies
  your detector — your widget renders in place of the default listing.
- Navigate to an unrelated folder — the default listing shows
  (i.e. your function returned `undefined` and the platform fell back).

## See also

- Source articles:
  - `help/develop/how-to/files/folder-content-preview.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-088` … `DG-FACT-094` and drifts `DG-FACT-DRIFT-036..038`.
- Reference packages:
  - `packages/ClinicalCase/src/package.ts:200-228` — async, sentinel-file
    detection (`dm.csv`), reads the file, wires `grok.functions.call`.
  - `packages/Plates/src/package.ts:48-54` — `outputs: dynamic`
    override returning `DG.ViewBase | undefined`; folder-name detector.
  - `packages/ClinicalCase/src/package.g.ts:37-43` — auto-emitted
    wrapper showing the canonical role + input/output annotations.
- Related skills:
  - `create-custom-file-viewers` (sibling — preview a single *file*,
    not a folder).
  - `file-handlers` (sibling — import the matched files as
    `DG.DataFrame[]` instead of rendering a widget).
