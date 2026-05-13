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
   Role string is `folderViewer` (camelCase, exact — `DG-FACT-088`).
   Signature: `(folder: DG.FileInfo, files: DG.FileInfo[]) =>
   Promise<DG.Widget | undefined>` (`DG-FACT-089`, `091`, `093`).
   Don't hand-author the `package.g.ts` wrapper — it's auto-generated.
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

2. **Fill in the detector and widget return.**
   Detect cheaply — sentinel filename (ClinicalCase) or `folder.name`
   pattern (Plates) — and `return undefined` on miss so the platform
   falls back to the default listing (`DG-FACT-091`). Wrap any
   `HTMLElement` with `DG.Widget.fromRoot(...)` (`DG-FACT-092`).
   ```typescript
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

3. **Optional — widen the output to `DG.ViewBase`.**
   Override the default `widget` output with
   `outputs: [{'name': 'result', 'type': 'dynamic'}]` on the decorator
   to return a full view (`DG-FACT-094`):
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

- **Preview never appears; default listing shows.** Role string
  wrong-case — must be `folderViewer` (`DG-FACT-088`, `089`).
- **Preview applies to every folder.** Detector too broad. Use a
  sentinel check, return `undefined` on miss (`DG-FACT-091`).
- **Article snippet won't compile.** Add explicit imports for
  `grok`, `ui`, `DG`.
- **Type error: `ViewBase is not assignable to Widget`.** Add
  `outputs: [{name:'result', type:'dynamic'}]` (`DG-FACT-094`).

## See also

- Source: `help/develop/how-to/files/folder-content-preview.md`.
- Knowledge: `DG-FACT-088`–`094`.
- Related skills: `create-custom-file-viewers`, `file-handlers`.
