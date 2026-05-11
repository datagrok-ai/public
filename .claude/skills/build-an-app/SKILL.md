---
name: build-an-app
version: 0.2.1
description: |
  Register a TypeScript function in a Datagrok package as a launchable
  application — entry point appears under Functions | Apps, owns a URL
  the team can share, and (optionally) binds URL path segments to
  arguments for deep links. For plugin authors who already have a view
  builder and need to wire it to the platform's app surface.
  Use when asked to "make a function launchable from the Apps page",
  "register an entry point a colleague can open via URL", or "wire a
  static method in src/package.ts as the launched function".
triggers:
  - register entry point in package.ts
  - mark a function as launchable
  - show up in functions catalog
  - bind url path segment to argument
  - add side navigation tree
  - deploy a launchable workflow
allowed-tools:
  - Read
  - Edit
  - Bash
harness-authored: true
---

# build-an-app

## When to use

You have a package and a view-rendering function, and need the platform
to expose it as a double-clickable entry under `Functions | Apps` with a
shareable URL. You may also want a side-navigation tree and a URL path
segment bound to a function argument for deep links.

## Prerequisites

- A package scaffold with `src/package.ts` declaring
  `export const _package = new DG.Package();` (related skill:
  `create-package`).
- `datagrok-tools` global (`npm i -g datagrok-tools`) and a Datagrok
  server alias in `~/.grok/config.yaml`.
- Canonical imports at the top of `src/package.ts`: `grok` from
  `datagrok-api/grok`, `DG` from `datagrok-api/dg` (lowercase), `ui`
  from `datagrok-api/ui`.

## Steps

1. **Declare the entry point with `@grok.decorators.app`.** Place it
   on a `static` method of a `PackageFunctions` class in
   `src/package.ts`. The decorator generates package metadata at build
   time — no runtime `grok.functions.register(...)` call is required.
   The article example returns `void`; reference packages return a
   `DG.ViewBase` so the platform can host it, which the snippet below
   follows. Knowledge: `DG-FACT-330`.
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG   from 'datagrok-api/dg';
   import * as ui   from 'datagrok-api/ui';

   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.app({
       name: 'Demo',
       description: 'Interactive demo',
       icon: 'images/icons/demoapp-icon.png',
     })
     static demoApp(): DG.ViewBase {
       const view = DG.View.create();
       view.name = 'Demo';
       view.root.appendChild(ui.divV([ui.h1('Welcome')]));
       return view;
     }
   }
   ```
   Expected: after `npm run build`, `src/package.g.ts` emits a wrapper
   for `demoApp` carrying `//meta.role: app`; the Apps page lists it.
   Source: `help/develop/how-to/apps/build-an-app.md:48-69`. Reference:
   `packages/Chem/src/package.ts:2919-2925`,
   `packages/ClinicalCase/src/package.ts:73-80`.

2. **(Optional) Bind a URL path segment to a function argument.** Deep
   links like `/apps/<AppName>/cohort-42` require exactly one string
   parameter annotated `metaUrl: true, optional: true`. Knowledge:
   `DG-FACT-332`.
   ```typescript
   @grok.decorators.app({name: 'Demo', description: 'Interactive demo'})
   static demoApp(
     @grok.decorators.param({options: {metaUrl: true, optional: true}})
     path?: string,
   ): DG.ViewBase {
     grok.shell.info(`Launched with path: ${path ?? '(none)'}`);
     return DG.View.create();
   }
   ```
   Expected: opening `/apps/<AppName>/foo/bar` (single-app package) or
   `/apps/<PackageName>/<AppName>/foo/bar` (multi-app) invokes
   `demoApp` with `path = 'foo/bar'`. Source:
   `help/develop/how-to/apps/build-an-app.md:62-67`. Reference:
   `packages/Chem/src/package.ts:2923-2924`,
   `packages/UsageAnalysis/src/package.ts:104-105`.

3. **(Optional) Attach a navigation tree with
   `@grok.decorators.appTreeBrowser`.** Add a sibling `static async`
   method that receives `treeNode: DG.TreeViewGroup` and returns
   `Promise<void>`. The platform calls it when the named app launches.
   `app:` must match the app `name` exactly (case-sensitive).
   Knowledge: `DG-FACT-331`.
   ```typescript
   @grok.decorators.appTreeBrowser({app: 'Demo'})
   static async demoAppTreeBrowser(
     treeNode: DG.TreeViewGroup,
   ): Promise<void> {
     treeNode.item('Home')
       .onSelected.subscribe((_) => grok.shell.info('Home selected'));
     const dataGroup = treeNode.group('Data');
     dataGroup.item('View 1')
       .onSelected.subscribe((_) => grok.shell.info('View 1'));
   }
   ```
   Expected: opening Demo shows the side panel; `package.g.ts`
   contains `//meta.role: appTreeBrowser` + `//meta.app: Demo`. Source:
   `help/develop/how-to/apps/build-an-app.md:110-133`. Reference:
   `packages/UsageAnalysis/src/package.ts:348-389`.

4. **Build the package.**
   ```bash
   npm install   # first-time only
   npm run build
   ```
   Expected: build exits 0, `src/package.g.ts` regenerates with
   `//meta.role: app` for the app function and (if step 2 applied)
   the `//input: string path { meta.url: ...; }` line.

5. **Publish in debug, then release.** Knowledge: `DG-FACT-334`,
   `DG-FACT-333`. Source:
   `help/develop/how-to/apps/build-an-app.md:21-28,96-100,851-857`.
   ```bash
   grok publish <host>            # debug: visible only to publisher
   grok publish <host> --release  # release: replaces installed version
   ```
   Expected: command exits 0. The app is reachable at
   `https://<host>/apps/<AppName>` when the package contains exactly
   one app, or at `https://<host>/apps/<PackageName>/<AppName>` when
   two or more apps live in the package.

## Common failure modes

- **App is missing from `Functions | Apps`.** Codegen wasn't
  regenerated, or the decorator is on a non-`static` method. Inspect
  `src/package.g.ts` for `//meta.role: app` on a wrapper; re-run
  `npm run build` — `.g.ts` is generated, not hand-edited.
- **`path` argument is always `undefined`.** The `@grok.decorators.param`
  options are missing `metaUrl: true`, the parameter has the wrong
  type (must be `string`), or two parameters declare it. Exactly one
  input must carry `{metaUrl: true}` (`DG-FACT-332`).
- **Import resolution fails on Linux CI but works on macOS.** Module
  path uses uppercase `datagrok-api/DG` as the article shows — the
  published npm path is lowercase. Always import
  `from 'datagrok-api/dg'` (`DG-FACT-DRIFT-BUILD-APP-003`).
- **`appTreeBrowser` never fires.** Its `app:` value doesn't match the
  registered app `name` (case-sensitive), or the decorator was attached
  to a non-`static` method.
- **Release pushes to everyone unintentionally.** `grok publish`
  defaults to debug (publisher-only); `--release` replaces the
  installed version for all users (`DG-FACT-334`). Verify the flag
  before publishing.

## Verification

- `src/package.g.ts` contains `//meta.role: app` for the app function
  and (if step 3) `//meta.role: appTreeBrowser` + `//meta.app: <AppName>`.
- After `grok publish <host>`, opening
  `https://<host>/apps/<AppName>` (single-app package) or
  `https://<host>/apps/<PackageName>/<AppName>` (multi-app package)
  renders the view from step 1.
- If step 2 applied, appending a tail to that URL
  (`.../<AppName>/foo/bar`) surfaces `foo/bar` in the `path` argument
  (verify via `grok.shell.info` or a debug breakpoint).
- If step 3 applied, the navigation tree appears alongside the view
  and item `onSelected` handlers fire.

## See also

- Source: `help/develop/how-to/apps/build-an-app.md` (mirror at
  `docs/_internal/articles-mirror/how-to/apps/build-an-app.md`).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-330` (`@grok.decorators.app`), `DG-FACT-331`
  (`appTreeBrowser`), `DG-FACT-332` (`metaUrl` param + alias),
  `DG-FACT-333` (URL form), `DG-FACT-334` (publish modes),
  `DG-FACT-DRIFT-BUILD-APP-003` (lowercase `dg` import),
  `DG-FACT-DRIFT-BUILD-APP-004` (canonical imports).
- Related skills:
  - `create-package` — produces the `src/package.ts` scaffold this
    skill mutates.
  - `routing` — deep-link state restoration on top of the
    `metaUrl: true` parameter introduced in step 2.
  - `add-info-panel` — sibling decorator-based registration pattern.
