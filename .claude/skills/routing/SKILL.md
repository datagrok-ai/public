---
name: routing
version: 0.1.0
description: |
  Make a Datagrok app respond to URL path segments so a link teammates
  paste into a browser reopens the app on exactly the view, table, and
  selection the sender was looking at, and so the breadcrumb and address
  bar follow the user's in-app navigation. For plugin authors whose app
  already launches but whose URL never moves past the home page.
  Use when asked to "share a link that opens at a specific state",
  "restore the exact selection a colleague was looking at", or "make
  every in-app navigation step its own shareable URL".
triggers:
  - share a link that reopens specific state
  - deep link into a selection
  - paste url to restore view
  - parse url path segments in app
  - drive breadcrumb and address bar from code
  - bind url suffix to app argument
allowed-tools:
  - Read
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# routing

## When to use

Your Datagrok app launches but the URL is frozen at the app root —
nothing teammates paste reopens the same table, viewer, or selection.
You want the address bar and breadcrumb to follow the user's in-app
navigation, and a pasted URL to rebuild that state on launch.

## Prerequisites

- App entry point already registered with `@grok.decorators.app` on a
  `static` method of `PackageFunctions` in `src/package.ts` — produce
  one with `build-an-app` if absent.
- `datagrok-tools` (`npm install -g datagrok-tools`) and a configured
  server alias for the publish step.

## Steps

1. **Bind the URL path suffix into a string parameter.** Annotate one
   `string` parameter with `metaUrl: true, optional: true` (see
   `DG-FACT-332`).
   ```typescript
   @grok.decorators.app({name: 'Cohorts'})
   static cohortsApp(
     @grok.decorators.param({options: {metaUrl: true, optional: true}})
     path?: string,
   ): DG.ViewBase {
     /* see step 2 */
     return DG.View.create();
   }
   ```
   `/apps/Cohorts/foo/bar` calls with `path = 'foo/bar'`. Refs:
   `packages/Tutorials/src/package.ts:148-150`.

2. **Split the path and branch on segment count.** Empty ⇒ fresh
   launch; non-empty ⇒ restore. `grok.data.testData(...)` is sync (see
   `DG-FACT-391`, `DG-FACT-392`).
   ```typescript
   const pathSegments = (path ?? '').split('/').filter((s) => s !== '');
   if (pathSegments.length === 0) openDefaultViews();
   else                          restoreFromPath(pathSegments);
   ```

3. **Set `view.path` whenever app state changes.** Assigning rewrites
   the address bar AND the breadcrumb in one call — there is no
   separate "push history". The full URL becomes
   `<function-base-path>/<view.path>`; never hardcode the base (see
   `DG-FACT-393`, `DG-FACT-DRIFT-ROUTING-003`).
   ```typescript
   function openDefaultViews(): void {
     const demog = grok.data.testData('demog');
     const view = grok.shell.addTableView(demog);
     view.addViewer(DG.Viewer.scatterPlot(demog));
     view.path = '/demog/All';     // ← drives URL + breadcrumb
     grok.shell.v = view;
   }
   ```
   Address bar and breadcrumb mirror `view.path`. Ref:
   `js-api/src/views/view.ts:169-171`.

4. **Rebuild state from path segments.** Always validate before
   feeding into a typed API — the `as DG.DemoDatasetName` cast at
   `routing.md:58` is a TS-only silencer. Reassign `view.path` once
   the view is on screen so the URL canonicalizes (a missing trailing
   label fills in to `/<table>/All`).
   ```typescript
   const DEMO: readonly DG.DemoDatasetName[] =
     ['wells', 'demog', 'biosensor', 'random walk',
      'geo', 'molecules', 'dose-response'];

   function restoreFromPath(segs: string[]): void {
     const name = segs[0];
     if (!(DEMO as readonly string[]).includes(name)) {
       grok.shell.warning(`Unknown dataset '${name}'`); return;
     }
     const label = segs[1] ?? 'All';
     const view = grok.shell.addTableView(
       grok.data.testData(name as DG.DemoDatasetName));
     applySelection(view, label);
     view.path = `/${name}/${label}`;
     grok.shell.v = view;
   }
   ```
5. **(Optional) Override the package's URL segment.** Add `meta.url`
   in `package.json` to rename the package-name segment (see
   `DG-FACT-389`, `DG-FACT-390`).
   ```json
   { "name": "@datagrok/cohorts", "meta": { "url": "/cohorts" } }
   ```
   Both `/apps/cohorts/...` and `/browse/apps/cohorts/...` serve the
   app; the old `/apps/Cohorts/...` redirects to the alias.

6. **Build and publish.**
   ```bash
   npm run build
   grok publish <host>            # debug, publisher-only
   grok publish <host> --release  # release for all users
   ```

## Common failure modes

- `path` always `undefined` — parameter not `string`, not flagged
  `metaUrl: true`, or two inputs claim it (`DG-FACT-332`).
- URL frozen at `/apps/<AppName>` — app never assigns `view.path`
  (`DG-FACT-393`); set on every state change.
- `TS2345 string is not assignable to DemoDatasetName` — narrow with
  an `includes` guard before calling `testData` (`DG-FACT-391`).
- Legacy URLs redirect after `meta.url` alias — expected
  (`DG-FACT-390`); update bookmarks.
- Multi-app packages prepend `<PackageName>/<AppName>` to `view.path`
  (`DG-FACT-DRIFT-ROUTING-003`); don't hardcode the base.

## See also

- Source: `help/develop/how-to/apps/routing.md`
  (mirror: `docs/_internal/articles-mirror/how-to/apps/routing.md`).
- Knowledge (`docs/_internal/knowledge/knowledge-graph.md`):
  `DG-FACT-332`, `DG-FACT-389` – `DG-FACT-393`,
  `DG-FACT-DRIFT-ROUTING-001` – `DG-FACT-DRIFT-ROUTING-003`.
- Related skills: `build-an-app` (registers the entry point this skill
  layers on), `create-package` (produces the scaffold).
