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

1. **Register the function as a file viewer.**
   The platform discovers viewers via the function role `fileViewer`
   (`DG-FACT-073`). Use a `static` method on `PackageFunctions`
   decorated with `@grok.decorators.fileViewer({fileViewer: '<ext>'})`
   — the canonical form in recent packages. The `fileViewer` key is
   REQUIRED (`DG-FACT-076`). The regenerated `src/package.g.ts` wrapper
   will have signature `(file: DG.FileInfo) => DG.View | Promise<DG.View>`.

   <details><summary>Variant: header-comment form (older packages)</summary>

   The article uses header comments instead of a decorator:
   ```typescript
   //name: previewEps
   //tags: fileViewer
   //meta.role: fileViewer
   //meta.fileViewer: eps
   //input: file file
   //output: view v
   export function previewEps(file: DG.FileInfo): DG.View { /* … */ }
   ```
   Both forms produce the same wrapper. Prefer the decorator in new code
   (`DG-FACT-073`, `DG-FACT-074`).
   </details>

2. **Write the viewer body (synchronous).**
   Mint the view, attach an empty host, return immediately, then fill
   the host inside `file.readAsBytes().then(...)` so the tab appears
   before the bytes load (`DG-FACT-077`). Set `view.name = file.name`
   so the tab caption matches the file (`DG-FACT-078`).
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
   Compare `packages/EpsViewer/src/package.ts:122-146`.
   Expected: build succeeds; `src/package.g.ts` carries
   `//meta.role: fileViewer` + `//meta.fileViewer: eps` (compare
   `packages/EpsViewer/src/package.g.ts:4-10`).

   <details><summary>Variant: async (when construction must await content)</summary>

   Use INSTEAD of the sync form when the view cannot be assembled until
   the bytes are parsed (load JSON, build a model, then mount). Declare
   `async`, return `Promise<DG.View>` (`DG-FACT-077`). Pick one form per
   viewer — don't write both.
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
   Refs: `packages/MetabolicGraph/src/package.ts:82-101`,
   `packages/FileEditors/src/package.g.ts:11-15` (async wrappers are
   regenerated with `Promise<any>` — both `Promise<DG.View>` and
   `Promise<any>` are accepted, `DG-FACT-076`).
   </details>

3. **Multiple extensions: comma-separate them.**
   `'mol,sdf,cif'` and `'mp3, wav, flac, ogg'` both work — whitespace
   around commas is tolerated (`DG-FACT-074`). One function registers
   for every listed extension.
   ```typescript
   @grok.decorators.fileViewer({fileViewer: 'mol,sdf,cif'})
   static viewStructure(file: DG.FileInfo): DG.View { /* … */ }
   ```
   Compare `packages/Media/src/package.ts:73,83` and
   `packages/BiostructureViewer/src/package.g.ts:89,97,105,113`.

4. **(Optional) Disambiguate two viewers on the same extension with
   `fileViewerCheck`.**
   When two viewers register the same extension (e.g. `.txt` *only* if
   it parses as a plate file), each declares
   `fileViewerCheck: '<Package>:<Function>'` — a `bool`-returning probe
   the platform calls with file content (`string` for text,
   `Uint8Array` for binary), routing only on truthy return
   (`DG-FACT-075`). The article omits this; see `packages/Plates`.
   The probe is itself a Datagrok function and needs its own
   registration — use `@grok.decorators.func()` (parameter and return
   types are inferred from the TypeScript signature).
   ```typescript
   @grok.decorators.fileViewer({fileViewer: 'txt',
                                fileViewerCheck: 'Plates:checkFileIsPlate'})
   static previewPlate(file: DG.FileInfo): DG.View { /* … */ }

   @grok.decorators.func()
   static checkFileIsPlate(content: string): boolean { /* … */ }
   ```
   Compare `packages/Plates/src/package.ts:56-66` (viewer) and
   `:225-230` (check).

5. **Build, publish, and let the platform auto-register.**
   ```bash
   npm install
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the `fileViewer` role
   plus `meta.fileViewer` IS the registration (`DG-FACT-073`).

## Common failure modes

- **Viewer never opens for the file.** Regenerated `src/package.g.ts`
  is missing `//meta.role: fileViewer` or `//meta.fileViewer: <ext>`.
  Inspect — both lines MUST be present (`DG-FACT-073`, `DG-FACT-074`).
  Re-run `npm install` if the wrapper looks stale.
- **Build error: `Property 'fileViewer' is missing` / not assignable to
  `string`.** You wrote `@grok.decorators.fileViewer()` or `({})`. The
  `fileViewer` extension-list key is REQUIRED in the decorator form
  (`DG-FACT-076`).
- **Article snippet ported to TypeScript won't compile — `file`
  implicitly `any`.** The article's JS snippet types the parameter as
  bare `file`; the TypeScript wrapper expects `file: DG.FileInfo`. Add
  the type annotation (and `: DG.View` on the return) when porting.
- **Two viewers fight over the same extension.** Both registered for
  `.txt` with no `fileViewerCheck`; the platform picks one
  non-deterministically. Add `fileViewerCheck` to each viewer that
  should be content-gated (`DG-FACT-075`).
- **Tab caption shows the function name, not the file.** Missing
  `view.name = file.name;` after `DG.View.create()` (`DG-FACT-078`).

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
