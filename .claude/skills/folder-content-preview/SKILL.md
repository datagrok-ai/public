---
name: folder-content-preview
version: 0.1.0
description: |
  Register a function that detects a known folder shape (sentinel
  filename or folder-name pattern) and returns a custom widget or view
  in place of the default file listing when the user opens that folder
  in Datagrok's Browse → Files tree. For package authors who own a
  recognizable directory shape — an SDTM study tree, a plates layout
  dump, a docking-run output — and want to surface a launcher,
  thumbnail, or summary panel inline.
  Use when asked to "show a custom panel when a folder opens", "render
  a launcher for a known directory shape", or "give a project-specific
  thumbnail in the file tree".
triggers:
  - custom panel for a folder
  - launcher for a directory shape
  - thumbnail in the file tree
  - detect a study folder and show a button
  - browse-tree preview for a directory
  - inline summary when opening a folder
allowed-tools:
  - Read
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# folder-content-preview

## When to use

Your package owns a *folder* shape — an SDTM study tree, a plates
layout, an HTRF run, a docking-output dump — and you want Datagrok's
file-share browser to render a custom panel/launcher when the user
clicks the folder, instead of the default file listing. For previewing
a single *file*, use `create-custom-file-viewers`. For *importing*
matched files as `DG.DataFrame[]`, use `file-handlers`.

## Prerequisites

- Familiarity with `DG.FileInfo` and `DG.Widget.fromRoot(...)`
  (knowledge: `DG-FACT-092`, `DG-FACT-093`).

## Steps

1. **Add the decorator + canonical signature to `src/package.ts`.**
   Open `src/package.ts` and add a `static async` method on
   `PackageFunctions` decorated with
   `@grok.decorators.folderViewer({'name': '<fnName>'})`. The role
   string is `folderViewer` — camelCase, exactly — because the
   platform discovers folder previews via that function role
   (`DG-FACT-088`); `folder-viewer` / `FolderViewer` won't register.
   The config arg is optional (`DG-FACT-090`); pass `{'name': ...}`
   to match the article and pin the registered function name.
   Match the canonical signature
   `(folder: DG.FileInfo, files: DG.FileInfo[]) => Promise<DG.Widget | undefined>`
   (`DG-FACT-089`, `DG-FACT-091`, `DG-FACT-093`,
   `js-api/src/const.ts:593-597`). Add the three `datagrok-api`
   imports at the top — the article's snippet omits them.
   Do *not* hand-author the header-comments form
   (`//meta.role: folderViewer`, `//input: ...`, `//output: ...`) in
   `src/package.g.ts` — that wrapper is auto-generated from your
   decorator (compare `packages/ClinicalCase/src/package.g.ts:37-43`).
   ```typescript
   // src/package.ts — skeleton
   import * as grok from 'datagrok-api/grok';
   import * as ui   from 'datagrok-api/ui';
   import * as DG   from 'datagrok-api/dg';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.folderViewer({'name': 'clinicalCaseFolderLauncher'})
     static async clinicalCaseFolderLauncher(
       folder: DG.FileInfo,
       files: DG.FileInfo[],
     ): Promise<DG.Widget | undefined> {
       return undefined; // detector body added in step 2
     }
   }
   ```
   Expected: file saves; the decorator import resolves (no red squiggle on
   `grok.decorators.folderViewer`).

2. **Fill in the detector and widget return.**
   Inside the method, detect the folder shape cheaply — scan `files`
   for a sentinel filename (ClinicalCase pattern) or test `folder.name`
   (Plates pattern) — and `return undefined` on miss so the platform
   falls back to the default listing (`DG-FACT-091`). On match, wrap
   any `HTMLElement` (built with `ui.div` / `ui.button`) with
   `DG.Widget.fromRoot(...)` (`DG-FACT-092`,
   `js-api/src/widgets/base.ts:374-376`).
   ```typescript
   // src/package.ts — sentinel-file pattern (SDTM Demographics)
   @grok.decorators.folderViewer({'name': 'clinicalCaseFolderLauncher'})
   static async clinicalCaseFolderLauncher(
     folder: DG.FileInfo,
     files: DG.FileInfo[],
   ): Promise<DG.Widget | undefined> {
     if (!files.some((f) => f.fileName.toLowerCase() === 'dm.csv'))
       return undefined;                         // opt out → fallback
     const csv = await grok.dapi.files.readAsText(`${folder.fullPath}/dm.csv`);
     const studyId = DG.DataFrame.fromCsv(csv).get('STUDYID', 0) ?? 'unknown';
     return DG.Widget.fromRoot(ui.div([
       ui.divText(`SDTM study: ${studyId}`),
       ui.button('Run study', () => grok.shell.info(`Launching ${studyId}`)),
     ]));
   }
   ```
   Note: the article stops at a placeholder `'START'` button calling
   `grok.shell.info('Folder contains SDTM data')`. The `readAsText` +
   `STUDYID` parsing above is an addition grounded in the production
   form at `packages/ClinicalCase/src/package.ts:200-228`.
   Expected: `src/package.g.ts` (regenerated next build) gains a wrapper
   with `//meta.role: folderViewer`, `//input: file folder`,
   `//input: list<file> files`, `//output: widget result`.

3. **Optional — widen the output to `DG.ViewBase`.**
   Codegen defaults to `outputs: [{name:'result', type:'widget'}]`
   (`tools/bin/utils/func-generation.ts:420-426`). To return a full
   view instead of a widget, override with
   `outputs: [{'name': 'result', 'type': 'dynamic'}]` on the decorator
   (`DG-FACT-094`); the regenerated wrapper then emits
   `//output: dynamic result`. Grounded in
   `packages/Plates/src/package.ts:48-54`.
   ```typescript
   @grok.decorators.folderViewer({
     outputs: [{'name': 'result', 'type': 'dynamic'}],
   })
   static async platesFolderPreview(
     folder: DG.FileInfo, files: DG.FileInfo[],
   ): Promise<DG.Widget | DG.ViewBase | undefined> {
     if (!folder.name?.toLowerCase().includes('plate')) return undefined;
     return getPlatesFolderPreview(files);
   }
   ```
   Expected: `package.g.ts` shows `//output: dynamic result` (compare
   `packages/Plates/src/package.g.ts:17-22`).

4. **Install deps and run `grok check`.**
   ```bash
   npm install
   grok check
   ```
   Expected: `grok check` exits `0`; `src/package.g.ts` is
   regenerated and contains the auto-emitted wrapper with the
   `//meta.role: folderViewer` header plus the `folder` / `files`
   inputs and `widget` (or `dynamic`) output.

5. **Publish; the platform auto-registers.**
   ```bash
   grok publish <host>           # add --release once stable
   ```
   Expected: `grok publish` exits `0`. No explicit `register(...)`
   call is needed — the function role `folderViewer` IS the
   registration. Opening any folder dispatches into every registered
   viewer until one returns non-`undefined` (`DG-FACT-088`,
   `DG-FACT-091`).

## Common failure modes

- **Preview never appears; default listing shows.** Role string is
  wrong-case. It MUST be `folderViewer` (camelCase) — `folder-viewer`
  / `FolderViewer` won't register. Inspect `src/package.g.ts`: it
  must contain `//meta.role: folderViewer` AND both `folder` +
  `files` inputs (`DG-FACT-088`, `DG-FACT-089`).
- **Preview applies to every folder, even unrelated ones.** Detector
  matched too broadly. Always start with a sentinel check
  (`files.some(f => f.fileName === '<expected>')` or a folder-name
  predicate) and `return undefined` on miss (`DG-FACT-091`). Multiple
  packages may each register a `folderViewer`; sloppy detection
  collides with siblings.
- **Article snippet copy-pasted, won't compile.** The article uses the
  decorator form but omits imports and references an unimported `DG`.
  Add explicit `import * as grok from 'datagrok-api/grok';` /
  `* as ui from 'datagrok-api/ui';` / `* as DG from 'datagrok-api/dg';`
  at the top of `src/package.ts`.
- **Type error: `ViewBase is not assignable to Widget`.** You returned
  a `DG.ViewBase` without overriding the output type. Add
  `outputs: [{name:'result', type:'dynamic'}]` to the decorator config
  (`DG-FACT-094`); codegen will then emit `//output: dynamic result`.
- **Stub demo widget shipped to production.** The article's `'START'`
  button is a placeholder calling `grok.shell.info('Folder contains
  SDTM data')`. Production code (`packages/ClinicalCase`) reads the
  sentinel file, parses it, and calls a real app entry via
  `grok.functions.call('<Pkg>:<App>')`.

## See also

- Source articles:
  - `help/develop/how-to/files/folder-content-preview.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-088` … `DG-FACT-094`.
- Reference packages:
  - `packages/ClinicalCase/src/package.ts:200-228` — async, sentinel-file
    detection (`dm.csv`), reads the file, wires `grok.functions.call`.
  - `packages/Plates/src/package.ts:48-54` — `outputs: dynamic`
    override returning `DG.Widget | DG.ViewBase | undefined`;
    folder-name detector.
  - `packages/ClinicalCase/src/package.g.ts:37-43` — auto-emitted
    wrapper showing the canonical role + input/output annotations.
- Related skills:
  - `create-custom-file-viewers` (sibling — preview a single *file*).
  - `file-handlers` (sibling — import matched files as
    `DG.DataFrame[]` instead of rendering a widget).
