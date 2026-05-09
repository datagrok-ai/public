---
name: routing
description: Wire a Datagrok app to URL paths so deep-links open it in a specific state
---

# routing

## When to use

Your package has an app and you want shareable URLs that restore a chosen
view, table, selection, or filter when a colleague opens the link. Triggers:
"deep-link my app", "URL state", "open this view via path", "share app link".

## Prerequisites

- A package with at least one function annotated `//meta.role: app`
  (see knowledge `DG-FACT-010`).
- `datagrok-tools` installed and the package builds locally (`webpack`).
- Datagrok server reachable; `grok publish` works for this package.

## Steps

1. **Declare the app function with a URL-bound `path` parameter.**
   Add the input header so Datagrok injects the URL path into your function.
   Knowledge: `DG-FACT-001`, `DG-FACT-010`.
   ```typescript
   //name: MyApp
   //input: string path {meta.url: true; optional: true}
   //input: string filter {optional: true}
   //meta.role: app
   //output: view result
   export function myApp(path?: string, filter?: string) {
     const segments = (path ?? '').split('/').filter((s) => s !== '');
     // dispatch on segments below
   }
   ```
   Expected: rebuilt `package.g.ts` contains the `meta.url: true` annotation
   on the `path` input. Query parameters (e.g. `?filter=...`) bind to other
   named inputs by name.

2. **Render default state when the path is empty.**
   `segments.length === 0` means the user opened the app fresh from the
   Apps page. Build the default view and set its `path`. Knowledge:
   `DG-FACT-002`, `DG-FACT-011`, `DG-FACT-012`.
   ```typescript
   if (segments.length === 0) {
     const demog = grok.data.testData('demog');
     const view = grok.shell.addTableView(demog);
     view.addViewer(DG.Viewer.scatterPlot(demog));
     view.path = '/demog/All';   // appended after the function base path
     grok.shell.v = view;
     return view;
   }
   ```
   Expected: opening `/apps/<PackageName>` shows the default table view, and
   the address bar updates to `/apps/<PackageName>/demog/All`.

3. **Dispatch on path segments to restore deep-linked state.**
   When segments are present, parse them and recreate the same view the
   sharer saw. Knowledge: `DG-FACT-013` (allowed `testData` names).
   ```typescript
   const [tableName, label = 'All'] = segments;
   const table = grok.data.testData(tableName as DG.DemoDatasetName);
   const view = grok.shell.addTableView(table);
   view.path = `/${tableName}/${label}`;
   if (label !== 'All' && filter) {
     const [col, value] = filter.split('=');
     if (!col || !value) throw new Error(`bad filter '${filter}'`);
     view.dataFrame.rows.match(`${col} = ${value}`).select();
   }
   grok.shell.v = view;
   ```
   Expected: opening `/apps/<PackageName>/demog/Subset?filter=age=30` opens
   the demog table with the matching rows selected.

4. **(Optional) Override the package URL segment.**
   To use a custom URL instead of the package name, edit `package.json`.
   Knowledge: `DG-FACT-005`.
   ```json
   {
     "meta": { "url": "/some/custom/route" }
   }
   ```
   Expected: app URL becomes `/apps/some/custom/route[/<AppName>]`. The
   package-name URL still works but auto-redirects (`DG-FACT-008`).

5. **(Optional) Override the app URL alias or make the app the package default.**
   Add `//meta.url:` at the function level. Use `//meta.url: /` to make this
   app the package default (no app segment in URL). Knowledge: `DG-FACT-006`,
   `DG-FACT-007`.
   ```typescript
   //name: MyApp
   //meta.role: app
   //meta.url: /todo
   //input: string path {meta.url: true; optional: true}
   ```
   Expected: app opens at `/apps/<PackageName>/todo` (or at the package
   root when alias is `/`).

6. **Build and publish.**
   ```bash
   webpack && grok publish
   ```
   Expected: package upload finishes without errors; the app is listed at
   `/apps/<PackageName>`.

## Common failure modes

- **`path` arrives `undefined` / empty even with deep link.** The input
  is not annotated `{meta.url: true}` (most common after copying from a
  non-app function). Fix: re-add the brace annotation, run `grok link`
  + `webpack` so `package.g.ts` regenerates.
- **`view.path` change doesn't update the address bar.** You set
  `view.basePath` instead — that getter is `@deprecated use path instead`
  (`DG-FACT-003`, drift `DG-FACT-DRIFT-001`). Fix: assign to `view.path`.
- **`Cannot read property 'split' of undefined`.** The function signature
  drifts from the header annotation (article tutorial bug,
  `DG-FACT-DRIFT-002`). Fix: keep the parameter name and the variable
  you call `.split` on identical, e.g. always `path`.
- **`testData('foo')` throws.** Only the seven names in
  `DG-FACT-DRIFT` knowledge entry `DG-FACT-013` are valid
  (`wells`, `demog`, `biosensor`, `random walk`, `geo`, `molecules`,
  `dose-response`). Fix: validate `tableName` against the allowed set
  before casting.

## Verification

- Open `/apps/<PackageName>` — default view loads, URL becomes
  `/apps/<PackageName>/<defaultPath>`.
- Copy the URL, open it in a new tab — the same state is restored.
- Edit the URL by hand to a different segment, hit Enter — the matching
  branch in your dispatcher runs and the view reflects it.

## See also

- Source articles:
  - `help/develop/how-to/apps/routing.md`
  - `help/develop/how-to/apps/build-an-app.md` (referenced from routing.md)
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts `DG-FACT-001`
    through `DG-FACT-013` and drifts `DG-FACT-DRIFT-001..004`.
- Related skills:
  - `build-an-app` (creates the app function this skill wires URL state into).
