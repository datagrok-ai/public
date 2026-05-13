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
   a `static` method of `PackageFunctions` in `src/package.ts`. Return
   `DG.ViewBase` (or `Promise<DG.ViewBase>`) so the platform can host
   the result (`DG-FACT-330`).
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

2. **(Optional) Bind a URL path segment to a function argument.**
   Annotate exactly one `string` parameter with
   `{metaUrl: true, optional: true}` (camelCase is canonical —
   `DG-FACT-332`, `DRIFT-BUILD-APP-002`).
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
   Opening `/apps/<AppName>/foo/bar` invokes `demoApp` with
   `path = 'foo/bar'`.

3. **(Optional) Attach a navigation tree with
   `@grok.decorators.appTreeBrowser`.** Decorator's `app:` value must
   match the app's `name` exactly (case-sensitive — `DG-FACT-331`).
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

4. **Build the package to regenerate `src/package.g.ts`.**
   ```bash
   npm install      # first time only
   npm run build
   ```

5. **Publish (debug for self-test, `--release` to promote).**
   `grok publish` defaults to debug (publisher-only); `--release`
   replaces the installed version for all users (`DG-FACT-333`,
   `DG-FACT-334`).
   ```bash
   grok publish <host>             # debug: visible only to publisher
   grok publish <host> --release   # release: replaces installed version
   ```
   URL: `https://<host>/apps/<AppName>` (single-app package) or
   `https://<host>/apps/<PackageName>/<AppName>` (multi-app).

## Common failure modes

- **Entry missing from `Functions | Apps`.** Codegen not refreshed or
  decorator on a non-`static` method. Verify `//meta.role: app` in
  `src/package.g.ts`; re-run `npm run build`.
- **`path` argument always `undefined`.** Missing `metaUrl: true`,
  wrong type, or two params declare it (`DG-FACT-332`).
- **Import fails on Linux CI.** Always import from lowercase
  `'datagrok-api/dg'` (`DRIFT-BUILD-APP-003`).
- **`appTreeBrowser` never fires.** `app:` value doesn't match the
  app's `name` (case-sensitive), or method isn't `static`.
- **Release reaches everyone unintentionally.** Confirm the
  `--release` flag before publishing (`DG-FACT-334`).
- **TS2561 "Did you mean to write 'category1'?".** `IBoxPlotSettings`
  has no `category` — use `category1`/`category2` or
  `categoryColumnNames: ['<col>']` (`DG-FACT-339`).

## See also

- Source: `help/develop/how-to/apps/build-an-app.md`.
- Knowledge: `DG-FACT-330`–`334`, `339`,
  `DRIFT-BUILD-APP-002`/`003`/`004`.
- Related skills: `create-package`, `routing`, `add-info-panel`.
