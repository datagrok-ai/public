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
   `string` parameter with `metaUrl: true, optional: true`. The platform
   passes everything after the app's base path into it — empty (or
   `undefined`) on a fresh launch from `Apps`. Knowledge: `DG-FACT-332`.
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
   Expected: opening `/apps/Cohorts/foo/bar` invokes `cohortsApp` with
   `path = 'foo/bar'`. Source:
   `help/develop/how-to/apps/routing.md:16-25,141-145`. Reference:
   `packages/Chem/src/package.ts:2923-2924`,
   `packages/Tutorials/src/package.ts:148-150`.

2. **Split the path and branch on segment count.** Empty path → "fresh
   launch, build defaults"; non-empty → "rebuild state encoded in the
   URL". `grok.data.testData(...)` is synchronous (no `await`);
   `DG.DemoDatasetName` allows only seven literal values — `wells`,
   `demog`, `biosensor`, `random walk`, `geo`, `molecules`,
   `dose-response` — so any URL-derived segment fed to it needs a
   runtime guard, not just the article's TS cast. Knowledge:
   `DG-FACT-391`, `DG-FACT-392`, `DG-FACT-DRIFT-ROUTING-001`.
   ```typescript
   const pathSegments = (path ?? '').split('/').filter((s) => s !== '');
   if (pathSegments.length === 0) openDefaultViews();
   else                          restoreFromPath(pathSegments);
   ```
   Expected: launching from `Apps` → `pathSegments` is `[]`; pasting
   `/apps/Cohorts/demog/All` → `pathSegments` is `['demog', 'All']`.
   Source: `help/develop/how-to/apps/routing.md:14-25,29-31,55-63`.
   Reference: `packages/Tutorials/src/package.ts:151-152`.

3. **Set `view.path` whenever app state changes.** Assigning to
   `view.path` rewrites the browser address bar AND the in-app
   breadcrumb in one call — no separate "push history" exists. The
   full URL becomes `<function-base-path>/<view.path>`; never hardcode
   the base path because it varies between one-app and multi-app
   packages and is overridable (step 5). Knowledge: `DG-FACT-393`,
   `DG-FACT-DRIFT-ROUTING-003`.
   ```typescript
   function openDefaultViews(): void {
     const demog = grok.data.testData('demog');
     const view = grok.shell.addTableView(demog);
     view.addViewer(DG.Viewer.scatterPlot(demog));
     view.path = '/demog/All';     // ← drives URL + breadcrumb
     grok.shell.v = view;
   }
   ```
   Expected: address bar shows `…/apps/<Package>/Cohorts/demog/All`
   (or `…/apps/Cohorts/demog/All` for one-app packages); breadcrumb
   mirrors the same segments. Source:
   `help/develop/how-to/apps/routing.md:32-50`. Reference:
   `js-api/src/views/view.ts:169-171`,
   `packages/UsageAnalysis/src/package.ts:461`.

4. **Rebuild state from path segments.** Validate each segment before
   feeding it to a typed API; the `as DG.DemoDatasetName` cast at
   `routing.md:58` silences the TS checker but does nothing at runtime.
   Reassign `view.path` once the view is on screen so the URL
   canonicalizes (missing trailing label fills in to `/<table>/All`).
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
   Expected: `/apps/Cohorts/demog/All` rebuilds the demog view and
   selection; `/apps/Cohorts/banana` triggers the toast warning
   instead of a TS-silent crash. Source:
   `help/develop/how-to/apps/routing.md:55-63,162-182`.

5. **(Optional) Override the package's URL segment.** Add `meta.url`
   in `package.json` to replace the package-name segment in every app
   URL the package hosts. The old `/apps/<PackageName>` keeps
   resolving but redirects on landing — share the alias instead. A
   matching `meta.url` on a function overrides the app-name segment
   for that one function. Knowledge: `DG-FACT-389`, `DG-FACT-390`.
   ```json
   { "name": "@datagrok/cohorts", "meta": { "url": "/cohorts" } }
   ```
   Expected: `/apps/cohorts/...` and `/browse/apps/cohorts/...` both
   serve the app; `/apps/Cohorts/...` rewrites to the alias on first
   render. Source:
   `help/develop/how-to/apps/routing.md:84-103,119-129`. Reference:
   `packages/UsageAnalysis/package.json:71-73`.

6. **Build and publish.**
   ```bash
   npm run build
   grok publish <host>            # debug, publisher-only
   grok publish <host> --release  # release for all users
   ```
   Expected: build exits 0; entry resolves at both `/apps/<AppName>`
   and `/browse/apps/<AppName>` (`DG-FACT-389`).

## Common failure modes

- **`path` is always `undefined`.** Parameter isn't typed `string`,
  isn't marked `metaUrl: true`, or two parameters claim it. Exactly
  one input must carry `{metaUrl: true}` (`DG-FACT-332`).
- **URL stays at `/apps/<AppName>` no matter what the user clicks.**
  The app never assigns to `view.path` after mutating state.
  `view.path` is the ONLY mutator that moves the address bar and
  breadcrumb in one call (`DG-FACT-393`); set it on every state change.
- **TS2345 `string is not assignable to DemoDatasetName`.** Raw URL
  segment fed into `grok.data.testData(...)`. Narrow with an `includes`
  guard against the seven literal values, or cast plus catch the
  runtime failure (`DG-FACT-391`).
- **Legacy URLs keep redirecting after a `meta.url` alias is added.**
  Expected — the platform rewrites `/apps/<PackageName>` to the alias
  on first render (`DG-FACT-390`). Update shared/bookmarked links.
- **Tutorial example assumes a single-app package.** The article's
  `/apps/TestPackage/demog/All` skips the app-name segment; in multi-app
  packages the form is `/apps/<PackageName>/<AppName>/demog/All`.
  `view.path` itself never changes — it's always appended to whatever
  base path the launcher resolved (`DG-FACT-DRIFT-ROUTING-003`).

## See also

- Source: `help/develop/how-to/apps/routing.md`
  (mirror: `docs/_internal/articles-mirror/how-to/apps/routing.md`).
- Knowledge (`docs/_internal/knowledge/knowledge-graph.md`):
  `DG-FACT-332`, `DG-FACT-389` – `DG-FACT-393`,
  `DG-FACT-DRIFT-ROUTING-001` – `DG-FACT-DRIFT-ROUTING-003`.
- Related skills: `build-an-app` (registers the entry point this skill
  layers on), `create-package` (produces the scaffold).
