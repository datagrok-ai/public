---
name: create-custom-file-viewers
version: 0.1.0
description: |
  Register a custom preview for a file extension so that when a user
  clicks a file of that type in Datagrok's file-share browser, your
  package renders it with your own UI instead of the generic text/binary
  view. For plugin authors who own a domain file format (`.pml`, `.cif`,
  in-house `.plate` text) and want a tailored visualization tab to open
  on click. Produces a `fileViewer` function (decorator or header form)
  wired into `src/package.ts`.
  Use when asked to "preview my files in Datagrok", "show a 3D structure
  when a `.cif` is opened", or "hook a custom preview for our file
  format".
triggers:
  - preview a file format in datagrok
  - render a 3d structure on file click
  - show custom preview when a file opens
  - hook a file extension to a preview tab
  - file-share browser custom preview
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# create-custom-file-viewers

## When to use

Your package owns a file format and you want the file-share browser to
render a tailored preview on click — a 3D structure for `.cif`/`.mol`,
a plate layout for plate-shaped `.txt`, a PDF/EPS render, an audio
player. You are NOT importing the file as a `DataFrame` (that's
`file-handlers`) and NOT extending right-click menus (`context-actions`).

## Prerequisites

- Familiarity with `DG.View` and `DG.FileInfo` (`DG-FACT-078`).

## Steps

1. **Register the function as a file viewer** via role `fileViewer` (`DG-FACT-073`). Decorator form on a `PackageFunctions` static method is canonical; the `fileViewer` key is required (`DG-FACT-076`).

   <details><summary>Variant: header-comment form (older packages)</summary>

   ```typescript
   //name: previewEps
   //tags: fileViewer
   //meta.role: fileViewer
   //meta.fileViewer: eps
   //input: file file
   //output: view v
   export function previewEps(file: DG.FileInfo): DG.View { /* … */ }
   ```
   Both forms produce the same wrapper; prefer the decorator in new code.
   </details>

2. **Write the viewer body (synchronous).** Return the view immediately, fill it after `file.readAsBytes().then(...)` so the tab opens before bytes load (`DG-FACT-077`). Set `view.name = file.name` for the tab caption (`DG-FACT-078`).
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as ui   from 'datagrok-api/ui';
   import * as DG   from 'datagrok-api/dg';
   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.fileViewer({fileViewer: 'eps'})
     static previewEps(file: DG.FileInfo): DG.View {
       const view = DG.View.create();
       view.name = file.name;
       const host = ui.div([], 'eps-viewer-host');
       view.append(host);
       file.readAsBytes().then((bytes) => host.appendChild(renderEps(bytes)));
       return view;
     }
   }
   ```
   Expected: build succeeds; `src/package.g.ts` carries `//meta.role: fileViewer` + `//meta.fileViewer: eps`.

   <details><summary>Variant: async (when construction must await content)</summary>

   Use INSTEAD of sync when the view cannot be assembled until bytes are parsed. Declare `async`, return `Promise<DG.View>` (`DG-FACT-077`). One form per viewer.
   ```typescript
   @grok.decorators.fileViewer({fileViewer: 'json'})
   static async escherViewer(file: DG.FileInfo): Promise<DG.View> {
     const view = DG.View.create('d4-escher-container');
     view.name = file.name;
     const mapData = JSON.parse(await file.readAsString());
     mountEscher(view.root, mapData);
     return view;
   }
   ```
   </details>

3. **Multiple extensions: comma-separate them** (whitespace tolerated, see `DG-FACT-074`).
   ```typescript
   @grok.decorators.fileViewer({fileViewer: 'mol,sdf,cif'})
   static viewStructure(file: DG.FileInfo): DG.View { /* … */ }
   ```

4. **(Optional) Disambiguate viewers on the same extension with `fileViewerCheck`** (`DG-FACT-075`). The probe must itself be a registered function (`@grok.decorators.func()`).
   ```typescript
   @grok.decorators.fileViewer({fileViewer: 'txt',
                                fileViewerCheck: 'Plates:checkFileIsPlate'})
   static previewPlate(file: DG.FileInfo): DG.View { /* … */ }

   @grok.decorators.func()
   static checkFileIsPlate(content: string): boolean { /* … */ }
   ```

5. **Build, publish, and let the platform auto-register.**
   ```bash
   npm install
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the `fileViewer` role
   plus `meta.fileViewer` IS the registration (`DG-FACT-073`).

## Common failure modes

- **Viewer never opens** — regenerated `package.g.ts` is missing `//meta.role: fileViewer` or `//meta.fileViewer: <ext>` (`DG-FACT-073`, `DG-FACT-074`). Re-run `npm install` if stale.
- **Build error: `Property 'fileViewer' is missing`** — empty decorator args. `fileViewer` key is required (`DG-FACT-076`).
- **Article snippet won't compile** — `file` needs `: DG.FileInfo` and the return needs `: DG.View`.
- **Two viewers fight over the same extension** — add `fileViewerCheck` to each (`DG-FACT-075`).
- **Tab caption shows function name** — missing `view.name = file.name` (`DG-FACT-078`).

## See also

- Source articles:
  - `help/develop/how-to/files/create-custom-file-viewers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-073` (role + signature), `DG-FACT-074` (extension list),
  `DG-FACT-075` (`fileViewerCheck`), `DG-FACT-076` (decorator form),
  `DG-FACT-077` (sync vs async), `DG-FACT-078` (`DG.View.create` +
  `FileInfo` content access).
- Reference packages:
  - `packages/EpsViewer/src/package.ts:122-146` — sync, decorator.
  - `packages/Plates/src/package.ts:56-66` — `fileViewerCheck`.
  - `packages/MetabolicGraph/src/package.ts:76-101` — async, header.
  - `packages/Media/src/package.ts:70-88` — multi-extension lists.
  - `packages/FileEditors/src/package.g.ts:11-15` — async wrapper.
- Related skills:
  - `file-handlers` — import the file as DataFrames (vs. preview).
  - `file-exporters` — write the active table out via Export menu.
  - `folder-content-preview` — preview a whole folder, not one file.
