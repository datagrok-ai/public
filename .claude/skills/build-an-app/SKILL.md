---
name: build-an-app
version: 0.3.1
description: |
  Register a TypeScript function inside a Datagrok package as a launchable
  entry — it surfaces under `Functions | Apps`, owns a shareable URL, and
  (optionally) binds a URL path segment to a string argument for deep
  links plus a side-navigation tree. For plugin authors who already have
  a view-building function and need to wire it to the platform's app surface.
  Use when asked to "make a function launchable from the Apps page",
  "register an entry point a colleague can open via URL", or "wire a
  static method in src/package.ts as the launched function".
triggers:
  - register entry point in package.ts
  - mark a function as launchable
  - show up in functions catalog
  - bind url path segment to argument
  - add side navigation tree
  - publish a launchable workflow
allowed-tools:
  - Read
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# build-an-app

## When to use

You have a Datagrok package and a function that builds a view, and need
the platform to expose it as a double-clickable entry under
`Functions | Apps` with a URL teammates can open. You may also want a
side-navigation tree or a URL path segment bound to an argument so deep
links restore state.

## Prerequisites

- A package scaffold with `src/package.ts` declaring
  `export const _package = new DG.Package();` — produce one with the
  `create-package` skill if absent.

## Steps

1. **Declare the entry point with `@grok.decorators.app`.** Place it on
   a `static` method of a `PackageFunctions` class in `src/package.ts`.
   The decorator emits package metadata at build time — no runtime
   `grok.functions.register(...)` call needed. The article returns
   `void`; production packages return `DG.ViewBase` (or `Promise` of
   one) so the platform can host the result. Knowledge: `DG-FACT-330`.
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
   for `demoApp` carrying `//meta.role: app`; the entry appears under
   `Functions | Apps`. Source:
   `help/develop/how-to/apps/build-an-app.md:48-69`. Reference:
   `packages/Chem/src/package.ts:2919-2925`,
   `packages/ClinicalCase/src/package.ts:73-80`.

2. **(Optional) Bind a URL path segment to a function argument.**
   Annotate exactly one `string` parameter with `metaUrl: true,
   optional: true` to capture everything after the app name in the URL.
   The article uses camelCase `metaUrl` exclusively
   (`help/develop/how-to/apps/build-an-app.md:64`); production packages
   follow the same form. The JS API source also accepts a dotted
   `'meta.url'` key (`js-api/src/decorators/functions.ts:101-102`), but
   camelCase is the canonical, documented variant
   (`DG-FACT-DRIFT-BUILD-APP-002`). Knowledge: `DG-FACT-332`.
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
   `packages/Chem/src/package.ts:2923-2924`.

3. **(Optional) Attach a navigation tree with
   `@grok.decorators.appTreeBrowser`.** Add a sibling `static async`
   method that receives `treeNode: DG.TreeViewGroup` and returns
   `Promise<void>`. The decorator's `app:` value must match the app's
   `name` exactly (case-sensitive). Knowledge: `DG-FACT-331`.
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
   Expected: launching Demo renders the side panel; `package.g.ts`
   contains `//meta.role: appTreeBrowser` and `//meta.app: Demo`.
   Source: `help/develop/how-to/apps/build-an-app.md:110-133`.
   Reference: `packages/ClinicalCase/src/package.ts:83-89`,
   `packages/UsageAnalysis/src/package.ts:348-389`.

4. **Build the package to regenerate `src/package.g.ts`.**
   ```bash
   npm install      # first time only
   npm run build
   ```
   Expected: build exits 0; `src/package.g.ts` regenerates with
   `//meta.role: app` for the entry function and (if step 2 applied)
   an `//input: string path { meta.url: ...; }` line on the wrapper.

5. **Publish (debug for self-test, `--release` to promote).**
   `grok publish` defaults to debug (per-user, visible only to
   publisher); add `--release` to replace the installed version for all
   users with proper privileges. Run the debug form first to validate,
   then re-run with `--release` once verified. Knowledge:
   `DG-FACT-333`, `DG-FACT-334`.
   ```bash
   grok publish <host>             # debug: visible only to publisher
   grok publish <host> --release   # release: replaces installed version
   ```
   Expected: command exits 0. The entry is reachable at
   `https://<host>/apps/<AppName>` when the package contains exactly
   one app, or at `https://<host>/apps/<PackageName>/<AppName>` when
   two or more apps live in the package. Source:
   `help/develop/how-to/apps/build-an-app.md:21-28,96-100,851-857`.

## Common failure modes

- **Entry missing from `Functions | Apps`.** Codegen wasn't refreshed,
  or the decorator sits on a non-`static` method. Inspect
  `src/package.g.ts` for `//meta.role: app` on a wrapper; if absent,
  re-run `npm run build` — `.g.ts` is generated, not hand-edited.
- **`path` argument always `undefined`.** The `@grok.decorators.param`
  options are missing `metaUrl: true`, the parameter is not typed
  `string`, or two parameters declare it. Exactly one input must carry
  `{metaUrl: true}` (`DG-FACT-332`).
- **Import resolution fails on Linux CI but works on macOS.** Module
  path uses uppercase `datagrok-api/DG` as the article shows — the
  published npm path is lowercase. Always import
  `from 'datagrok-api/dg'` (`DG-FACT-DRIFT-BUILD-APP-003`).
- **`appTreeBrowser` never fires.** The decorator's `app:` value
  doesn't match the registered app `name` (case-sensitive), or the
  decorator sits on a non-`static` method.
- **Release reaches everyone unintentionally.** `grok publish` alone
  is publisher-only; `--release` replaces the installed version for
  all users. Confirm the flag before publishing (`DG-FACT-334`).
- **`view.boxPlot({...})` build fails with TS2561 "Did you mean to
  write 'category1'?".** `IBoxPlotSettings` has no `category` key —
  the categorical-axis key is `category1` (optional second is
  `category2`), or use the array form `categoryColumnNames: ['<col>']`;
  `value` / `valueColumnName` are both valid. Datagrok viewer settings
  use verbose property names overall (`xColumnName`, `valueColumnName`,
  `colorColumnName`); `IScatterPlotSettings` is the laxer exception
  (both `x` and `xColumnName` exist). Reference:
  `packages/ApiSamples/scripts/demo/stock-broker.js:44-47`
  (`DG-FACT-339`).

## See also

- Source: `help/develop/how-to/apps/build-an-app.md`
  (mirror: `docs/_internal/articles-mirror/how-to/apps/build-an-app.md`).
- Knowledge (`docs/_internal/knowledge/knowledge-graph.md`):
  `DG-FACT-330` (`@grok.decorators.app`), `DG-FACT-331`
  (`appTreeBrowser`), `DG-FACT-332` (`metaUrl` parameter),
  `DG-FACT-333` (URL form), `DG-FACT-334` (publish modes),
  `DG-FACT-339` (viewer settings naming — `category1`/`valueColumnName`
  not `category`/`value`),
  `DG-FACT-DRIFT-BUILD-APP-002` (`metaUrl` vs dotted `'meta.url'`),
  `DG-FACT-DRIFT-BUILD-APP-003` (lowercase `dg` import),
  `DG-FACT-DRIFT-BUILD-APP-004` (canonical imports).
- Related skills: `create-package` (produces the scaffold this skill
  mutates), `routing` (deep-link state on top of step 2's `metaUrl`),
  `add-info-panel` (sibling decorator-based registration pattern).
