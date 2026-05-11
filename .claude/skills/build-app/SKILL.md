---
name: build-app
description: Build a Datagrok application — register an app entry point, capture URL path, add a navigation tree, and deploy
harness-authored: true
---

# build-app

## When to use

Your package needs an **app** — a function with `meta.role: app` that
appears under `Functions | Apps`, owns a URL at `/apps/<PKG>` (or
`/apps/<PKG>/<APP>`), and renders a view. Triggers: "create a
Datagrok app", "wire an entry point", "add a side-navigation tree",
"make this function launchable from `/apps/...`".

## Prerequisites

- `datagrok-tools` global (`npm i -g datagrok-tools` — provides `grok`);
  article also asks for `npm i -g webpack webpack-cli`.
- Datagrok server alias in `~/.grok/config.yaml` (`grok config add`).
- A package scaffold (fresh from step 1, or existing).
- `DG.View` / `DG.ViewBase` and `ui.div*` familiarity for view bodies
  (related skill: `create-custom-view`).

## Steps

1. **Create the package (skip if you already have one).** `--ide=vscode`
   wires VS Code debug config — Windows-only per the article; harmless to
   omit elsewhere.
   ```bash
   grok create MyPackage --ide=vscode
   cd MyPackage
   npm install
   ```
   Expected: a `src/package.ts` with `export const _package = new DG.Package();`,
   a `package.json` with `datagrok-api` in `dependencies`, and (if `--ide`)
   `.vscode/launch.json`.

2. **Add an app function** — generates a stub `@grok.decorators.app`
   method on `PackageFunctions`.
   ```bash
   grok add app MyApp
   ```
   Expected: `src/package.ts` gains a static method with
   `@grok.decorators.app({name: 'MyApp', ...})`; on next build,
   `src/package.g.ts` emits a wrapper with `//meta.role: app`
   (`DG-FACT-010`). Compare
   `packages/UsageAnalysis/src/package.g.ts:36`, `:51`, `:61`, `:73`.

3. **Author the entry point (decorator form, canonical).** Declare
   exactly one URL-path parameter with `{meta.url: true, optional: true}`
   so the platform binds the path tail to that argument
   (`DG-FACT-001`). Query-string params bind by name.
   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG   from 'datagrok-api/dg';
   import * as ui   from 'datagrok-api/ui';

   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.app({
       name: 'MyApp',
       description: 'What the app does (shown in Functions | Apps)',
       icon: 'images/icons/myapp-icon.png',
     })
     static myApp(
       @grok.decorators.param({options: {'meta.url': true, 'optional': true}}) path?: string,
     ): DG.ViewBase {
       const view = DG.View.create();
       view.name = 'MyApp';
       view.root.appendChild(ui.divV([ui.h1('Welcome'), ui.divText(`path: ${path ?? ''}`)]));
       return view;
     }
   }
   ```
   Pick ONE output contract:
   - **Return `DG.ViewBase`** + codegen emits `//output: view result` —
     host places the view (Browse preview tile, main shell, etc.). Do
     NOT also call `grok.shell.addView` / `addPreview` / set
     `grok.shell.v` inside the body — that double-places the view.
   - **Return `void`** and place views yourself with
     `grok.shell.addView(...)` / `grok.shell.addTableView(...)`
     (`DG-FACT-011`, `DG-FACT-012`). Use when the app opens several
     table views or owns its own dock layout.

   Expected after `npm run build`: `src/package.g.ts` contains
   `//input: string path { meta.url: true; optional: true }`,
   `//output: view result` (for the view-returning form), and
   `//meta.role: app`. Compare
   `packages/UsageAnalysis/src/package.g.ts:28-41`.

4. **(Optional) Add a navigation tree** — link an
   `@grok.decorators.appTreeBrowser` function to the app via its
   `app:` option. The platform calls it with a root `treeNode` and
   renders the result alongside the app view.
   ```typescript
   @grok.decorators.appTreeBrowser({app: 'MyApp'})
   static async myAppTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
     treeNode.item('Home').onSelected.subscribe((_) => grok.shell.info('Home'));
     const data = treeNode.group('Data');
     data.item('View 1').onSelected.subscribe((_) => grok.shell.info('View 1'));
     data.item('View 2').onSelected.subscribe((_) => grok.shell.info('View 2'));
   }
   ```
   Expected: codegen emits `//meta.role: appTreeBrowser` +
   `//meta.app: MyApp` in `package.g.ts`. Compare
   `packages/UsageAnalysis/src/package.ts:348-389`,
   `packages/UsageAnalysis/src/package.g.ts:80-92`.

5. **(Optional) Override the URL alias.** Default app URL is
   `/apps/<PackageName>` for one app and
   `/apps/<PackageName>/<AppName>` for two-plus (`DG-FACT-004`).
   - Per-app override — function-level annotation `//meta.url: /<alias>`
     or `'url': '/<alias>'` in the decorator options (`DG-FACT-006`).
   - Package-wide override — `"meta": {"url": "/<segment>"}` in
     `package.json` (`DG-FACT-005`; see
     `packages/UsageAnalysis/package.json`).
   - Set an app's `meta.url` to exactly `/` to make it the package
     default — the URL collapses to just the package URL with no app
     segment (`DG-FACT-007`).

6. **Build and publish.** `grok publish` deploys in **debug** mode —
   only you see the new version; everyone else keeps the prior
   release. Add `--release` to replace the installed version for all
   users.
   ```bash
   npm run build
   grok publish <host>            # debug: visible only to you
   grok publish <host> --release  # release: replaces installed version
   ```
   Expected: publish exits 0; URL `https://<host>/apps/<PackageName>`
   (or `.../<PackageName>/MyApp`) opens the app.

## Common failure modes

- **App doesn't appear in `Functions | Apps`.** Codegen not
  regenerated, or decorator missing. Inspect `src/package.g.ts` for
  `//meta.role: app` on a wrapper for your function (`DG-FACT-010`).
  Re-run `npm run build` — `.g.ts` is generated, not hand-edited.
- **`path` argument is always `undefined`.** `meta.url: true` missing
  on the param, or two params declare it. Exactly one input must
  carry `{meta.url: true}` (`DG-FACT-001`); query string binds by name.
- **View shows up twice.** Function declares `//output: view result`
  AND calls `grok.shell.addView`/`addPreview`/sets `grok.shell.v`.
  Pick one shape: return the view OR return void and place it
  yourself, never both.
- **`appTreeBrowser` never fires.** Decorator's `app:` value doesn't
  match the registered app `name` (case-sensitive).
- **URL `/apps/<PackageName>` returns 404.** Package has 2+ apps —
  one-app form doesn't resolve. Use `/apps/<PackageName>/<AppName>`,
  or set `meta.url: '/'` on the default app (`DG-FACT-004`,
  `DG-FACT-007`).
- **`grok publish` fails with `unknown host`.** No alias in
  `~/.grok/config.yaml`. Run `grok config add`, or omit the host arg.

## Verification

- `npm run build` exits 0; `grok publish <host>` exits 0.
- `src/package.g.ts` contains a wrapper for the app function with
  `//meta.role: app` and (for view-returning apps) `//output: view result`.
- Open `https://<host>/apps/<PackageName>` (or
  `.../<PackageName>/<AppName>`): the view from step 3 renders.
- Append a tail to the URL (`.../<AppName>/foo/bar`) and observe the
  `path` argument receives `foo/bar` in the view body.
- (If step 4) The side navigation tree appears next to the app view
  and `onSelected` handlers fire.

## See also

- Source articles:
  - `help/develop/how-to/apps/build-an-app.md` (mirror in
    `docs/_internal/articles-mirror/how-to/apps/build-an-app.md`)
  - `help/develop/how-to/apps/routing.md` — companion article on URL
    routing, view paths, package/app URL overrides.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-001` (`meta.url: true` URL path param),
  `DG-FACT-004` (default `/apps/<pkg>[/<app>]` prefix),
  `DG-FACT-005` (package-level `meta.url` in `package.json`),
  `DG-FACT-006` (function-level `//meta.url` override),
  `DG-FACT-007` (`meta.url: '/'` collapses to package URL),
  `DG-FACT-010` (`meta.role: app` registers the function),
  `DG-FACT-011` (`grok.shell.v = view`),
  `DG-FACT-012` (`grok.shell.addTableView` returns `TableView`).
- Reference packages:
  - `packages/UsageAnalysis/src/package.ts:99-116` — `usageAnalysisApp`
    with `meta.url: '/'`, `browsePath`, `meta.url: true` path param;
    codegen at `packages/UsageAnalysis/src/package.g.ts:28-41`.
  - `packages/UsageAnalysis/src/package.ts:348-389` — two
    `@grok.decorators.appTreeBrowser` examples bound to different apps.
  - `packages/ClinicalCase/src/package.ts:73-80` — minimal
    `@grok.decorators.app` delegating to a view-builder helper.
- Related skills: `create-package` (prerequisite),
  `create-custom-view` (the view body), `add-info-panel` (sibling
  role-based codegen pattern).
