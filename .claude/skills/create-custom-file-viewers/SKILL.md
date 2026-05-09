---
name: create-custom-file-viewers
description: Register a custom file viewer so Datagrok previews specific file extensions with your own UI in the file-share browser
---

# create-custom-file-viewers

## When to use

Your package owns a file format (`.eps`, `.mol`, an in-house `.plate`-style
text file, …) and you want Datagrok's file-share browser to render a
custom preview when the user clicks one. Triggers: "preview my files",
"show a 3D model when someone opens a `.cif`", "render our plate layout
for `.txt` files that look like plates". For ID detection in *cells*
(not file previews), use `register-identifiers` instead.

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the package root.
- `datagrok-api` available (`import * as grok / ui / DG`); the article's
  snippet omits imports and won't compile as written (drift
  `DG-FACT-DRIFT-026`).
- Familiarity with `DG.View` and `DG.FileInfo` (`DG-FACT-078`).

## Steps

1. **Pick a registration form (decorator vs. header comments).**
   The platform discovers viewers via the function role `fileViewer`
   (`DG-FACT-073`). Two equivalent surfaces:
   - **Decorator (canonical, used by every recent package):** put a
     `static` method on `PackageFunctions` decorated with
     `@grok.decorators.fileViewer({fileViewer: '<ext>'})`. `fileViewer`
     is REQUIRED in this form (`DG-FACT-076`).
   - **Header comments (article form):** annotate an exported function
     with `//meta.role: fileViewer`, `//meta.fileViewer: <ext>[,<ext>…]`,
     `//input: file file`, `//output: view v` (`DG-FACT-073`,
     `DG-FACT-074`).
   Either way the `package.g.ts` wrapper is regenerated to expose
   `static <method>(file: DG.FileInfo): DG.View | Promise<DG.View>`.

2. **Add the viewer function (sync form).**
   Mint a view, attach an empty host, return immediately, and finish
   filling the host inside `file.readAsBytes().then(...)` so the tab
   appears before the bytes arrive (`DG-FACT-077`,`DG-FACT-078`). Set
   `view.name = file.name` so the tab caption matches the file.
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
       file.readAsBytes().then((bytes) => {
         host.appendChild(renderEps(bytes));   // your renderer
       });
       return view;
     }
   }
   ```
   Expected: build succeeds; `src/package.g.ts` gains a wrapper with
   `//meta.role: fileViewer` + `//meta.fileViewer: eps`
   (compare `packages/EpsViewer/src/package.g.ts:6-10`).

3. **Or use the async form when content must be awaited first.**
   Declare `async` and return `view` after the `await`. Use this when
   construction depends on the file's content (e.g. parse to JSON, build
   a typed object, then mount) — `FileEditors.previewPdf` and
   `MetabolicGraph.escherFileViewer` are the canonical references
   (`DG-FACT-077`).
   ```typescript
   @grok.decorators.fileViewer({fileViewer: 'pdf'})
   static async previewPdf(file: DG.FileInfo): Promise<DG.View> {
     const bytes = await file.readAsBytes();
     const view  = DG.View.create();
     view.name   = file.name;
     view.append(renderPdf(bytes));     // your renderer
     return view;
   }
   ```
   Expected: wrapper signature is
   `static previewPdf(file: DG.FileInfo): Promise<any>` — both `Promise<DG.View>`
   and `Promise<any>` are accepted (`DG-FACT-076`).

4. **Multiple extensions: comma-separate them.**
   `'mol,sdf,cif'` and `'mp3, wav, flac, ogg'` both work — whitespace is
   tolerated (`DG-FACT-074`). One viewer registers for every listed
   extension.
   ```typescript
   @grok.decorators.fileViewer({fileViewer: 'mol,sdf,cif'})
   static viewStructure(file: DG.FileInfo): DG.View { /* … */ }
   ```

5. **(Optional) Disambiguate two viewers on the same extension with `fileViewerCheck`.**
   When two viewers in your package register the same extension (e.g.
   "`.txt` *only* if it's a plate file"), each one declares a
   `fileViewerCheck: '<Pkg>:<Fn>'` pointing to a `bool`-returning
   function that the platform calls with the file's content
   (`string` for text, `Uint8Array` for binary) (`DG-FACT-075`). The
   article doesn't mention this — see drift `DG-FACT-DRIFT-025`.
   ```typescript
   @grok.decorators.fileViewer({fileViewer: 'txt',
                                fileViewerCheck: 'Plates:checkFileIsPlate'})
   static previewPlate(file: DG.FileInfo): DG.View { /* … */ }

   @grok.decorators.func({})
   static checkFileIsPlate(content: string): boolean { /* … */ }
   ```
   Expected: only the viewer whose `fileViewerCheck` returns `true` is
   offered for that file (compare `packages/Plates/src/package.ts:56,176`).

6. **Build, publish, and let the platform auto-register.**
   ```bash
   npm install
   grok publish <host>            # add --release once stable
   ```
   No explicit `register(...)` call is needed — the function role
   `fileViewer` plus `meta.fileViewer` is the registration
   (`DG-FACT-073`).

## Common failure modes

- **Viewer never opens for the file.** Header is missing
  `//meta.role: fileViewer` or the decorator's `fileViewer:` key
  (REQUIRED — `DG-FACT-076`). Inspect the regenerated
  `src/package.g.ts`: it MUST contain `//meta.role: fileViewer` AND
  `//meta.fileViewer: <ext>` (`DG-FACT-073`,`DG-FACT-074`).
- **Build error: `fileViewer` is not assignable to type `string`.**
  You wrote `@grok.decorators.fileViewer({})` or
  `@grok.decorators.fileViewer()`. The `fileViewer` extension list is
  required in the decorator form, even when only the role matters
  (`DG-FACT-076`).
- **Article snippet copy-pasted, won't compile.** The article uses `var`
  and references `ui`, `DG`, `NGL` without imports
  (drift `DG-FACT-DRIFT-026`). Add explicit
  `import * as DG from 'datagrok-api/dg';` (and `ui`, `grok` as needed)
  and use `const`/`let`.
- **Two viewers fight over the same extension.** Both registered for
  `.txt` with no `fileViewerCheck`; the platform picks one
  non-deterministically. Add `meta.fileViewerCheck` /
  `fileViewerCheck:` to *each* viewer that should be content-gated
  (`DG-FACT-075`).
- **Reading docs that mention a `fileViewer-<ext>` tag.** Outdated
  JSDoc form — current authoring uses the structured `meta.fileViewer`
  annotation only (drift `DG-FACT-DRIFT-024`). Don't try to register
  via `tags`.
- **Tab caption shows the function name instead of the file name.**
  You forgot `view.name = file.name;` after `DG.View.create()`
  (`DG-FACT-078`).

## Verification

- `grok publish <host>` exits `0`; the regenerated `src/package.g.ts`
  contains an exported wrapper with `//meta.role: fileViewer` and
  `//meta.fileViewer: <ext>`.
- In Datagrok, open **Browse → Files** (or any file share), select a
  file with the registered extension — your view opens in a new tab
  whose caption equals `file.name`.
- Toggle the tab off and re-select the file; the tab re-opens (i.e.
  the view is rebuilt per click, not cached).
- For multi-viewer extensions: rename the test file to a content shape
  that only ONE `fileViewerCheck` accepts; only that viewer should
  appear.

## See also

- Source articles:
  - `help/develop/how-to/files/create-custom-file-viewers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-073` … `DG-FACT-078` and drifts `DG-FACT-DRIFT-024..026`.
- Reference packages:
  - `packages/EpsViewer/src/package.ts:122-146` — sync, single extension.
  - `packages/Plates/src/package.ts:56-66,176-226` — decorator +
    `fileViewerCheck` for content disambiguation.
  - `packages/FileEditors/src/package.g.ts:11-15` — async wrapper signature.
- Related skills:
  - `register-identifiers` (sibling — pattern detection for IDs *inside*
    files / cells, not file previews).
