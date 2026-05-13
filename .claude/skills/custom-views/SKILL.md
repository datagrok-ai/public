---
name: custom-views
version: 0.1.0
description: |
  Register a workspace screen from a Datagrok package that isn't a
  `TableView` — your own data backing (notebook, study summary, designer
  canvas, dashboard), your own URL path, your own icon, and state that
  restores when a saved project reopens. Subclass `DG.ViewBase` and
  expose a function tagged `view` so the platform dispatches URL paths
  into your class and serializes state back on save.
  Use when asked to "register a main-area screen from my package",
  "open my plugin at a deep-linked URL like /foo/<id>", or "ship a
  dashboard that survives project save and reopen".
triggers:
  - register a main-area screen
  - deep-linkable dashboard from a package
  - non-table workspace tab
  - url-routable screen with restorable state
  - package-owned dashboard with icon
  - open my plugin at /foo/<id>
allowed-tools:
  - Read
  - Edit
  - Bash
harness-authored: true
---

# custom-views

## When to use

You need a workspace screen that isn't a `TableView` — a dashboard,
notebook host, study summary — that deep-links via URL, restores from a
saved project, and carries a package-owned icon. For a one-shot tab
with no URL/project needs, `grok.shell.newView` is enough — see Step 1.

## Prerequisites

- A package scaffold (`grok create <Name>`); code lives under `src/`.
- Standard `datagrok-api` imports (article omits these):
  ```typescript
  import * as grok from 'datagrok-api/grok';
  import * as DG   from 'datagrok-api/dg';
  import * as ui   from 'datagrok-api/ui';
  ```
- Familiarity with the `routing` skill if you'll wire URL state.

## Steps

1. **Decide which API tier you need.**
   `grok.shell.newView(name, children?, options?)` builds, names,
   registers, and docks a `View` in one call (`DG-FACT-213`,
   `js-api/src/shell.ts:234`). Fine for production, but the returned
   `View` does NOT dispatch URL paths or serialize state. Subclass
   `DG.ViewBase` only when you need URL routing, project restore, or a
   per-class icon (`DG-FACT-214`):
   ```typescript
   // Ad-hoc: one-liner, no routing/state — keep for prototypes & dialogs
   const v = grok.shell.newView('Quick View', [ui.divText('Hi!')]);
   ```
   Expected: a new tab "Quick View" appears with the children appended
   to `v.root`. Stop here if you don't need routing or state.

2. **Subclass `DG.ViewBase` with the standard constructor signature.**
   Use the `DG.` prefix — package scaffolds import `* as DG` from
   `datagrok-api/dg` (`DG-FACT-224`). The third constructor arg
   (`createHost`) defaults to `true`; package code MUST NOT pass
   `false`, that path is reserved for the internal `View` subclass that
   wraps a Dart handle (`DG-FACT-214`):
   ```typescript
   // src/views/dashboard-view.ts
   export class DashboardView extends DG.ViewBase {
     static readonly TYPE = 'Dashboard';
     static readonly PATH = '/dashboard';
     id: string | null = null;

     constructor(params: any = null, path: string = '') {
       super(params, path);  // do NOT pass createHost=false
       this.id = params?.id ?? path.replace(`${DashboardView.PATH}/`, '') ?? null;
       if (this.id) this.open(this.id);
     }
     // ...
   }
   ```
   Expected: `tsc` compiles; `new DashboardView({id: 'foo'}, '/dashboard/foo')`
   constructs without throwing.

3. **Override the contract surface — name, path, icon, state.**
   Defaults live in `js-api/src/views/view.ts:45-204`. Store identity
   on `this`, expose it via `path`, serialize via `saveStateMap` /
   `loadStateMap`, dispatch URLs via `handlePath` / `acceptsPath`
   (`DG-FACT-215`):
   ```typescript
   get type()    { return DashboardView.TYPE; }
   get name()    { return this.id ? `Dashboard: ${this.id}` : 'Dashboard'; }
   get path()    { return `${DashboardView.PATH}/${this.id ?? ''}`; }
   get helpUrl() { return '/help/develop/how-to/views/custom-views.md'; }

   getIcon(): HTMLElement {
     const i = document.createElement('i');
     i.className = 'grok-icon fal fa-th-large';
     return i;
   }

   saveStateMap()       { return {id: this.id}; }
   loadStateMap(s: any) { this.open(s['id']); }

   handlePath(p: string)         { this.open(p.replace(`${DashboardView.PATH}/`, '')); }
   acceptsPath(p: string): boolean { return p.startsWith(DashboardView.PATH); }

   open(id: string) {
     this.id = id;
     ui.empty(this.root);
     this.root.append(ui.divText(`Loading ${id}…`));
   }
   ```
   Expected: `view.path` returns `/dashboard/<id>`; `view.getIcon()`
   returns a non-null `HTMLElement`; `saveStateMap()` returns JSON.

4. **Register the view-producing function with `tags: view`.**
   The platform turns a function into a registered view by reading four
   header lines: `tags: view`, `input: map params`, `input: string path`,
   `output: view result` (`DG-FACT-216`). Without `tags: view` the
   function is callable but never dispatched on URL paths or wired into
   the navigation bar:
   ```typescript
   // src/package.ts
   //name: Dashboard
   //description: Opens the package's dashboard view
   //input: map params
   //input: string path
   //tags: view
   //output: view result
   export function dashboardView(params: any = null, path: string = '') {
     return new DashboardView(params, path);
   }
   ```
   Expected: after `grok link` regenerates `package.g.ts`, the new
   function appears with `//tags: view` and `//output: view result`
   intact. Notebooks (`packages/Notebooks/src/package.js:528-535`) is
   the canonical reference — only one production package ships this
   annotation today (`DG-FACT-217`).

5. **Build, publish, and exercise the deep-link URL.**
   ```bash
   webpack && grok publish
   ```
   Expected: upload succeeds; opening `https://<host>/dashboard/abc`
   instantiates `DashboardView` with `path='/dashboard/abc'`,
   `acceptsPath` returns `true`, `handlePath` runs, and the tab opens
   already pointed at `abc`.

## Common failure modes

- **Tab opens but URL doesn't change / deep-link doesn't restore.** You
  used `grok.shell.newView` — it returns a leaf `View` with no path
  dispatcher (`DG-FACT-213`). Move the body into a `DG.ViewBase`
  subclass and register it via `tags: view`.
- **Reopening a saved project shows an empty tab.** State lives in
  module-level variables instead of `saveStateMap()` / `loadStateMap()`.
  Move identity onto `this` and round-trip via the state map
  (`DG-FACT-215`).
- **Function compiles but never dispatches the path.** Missing
  `tags: view`, or `input: string path` was renamed. The platform reads
  the EXACT header tokens; rename and re-run `grok link`
  (`DG-FACT-216`).
- **Type error on `class … extends ViewBase`.** Default scaffold
  imports `* as DG`, so the canonical form is `DG.ViewBase`. Bare
  `ViewBase` only works if you `import {ViewBase} from 'datagrok-api/dg'`
  (`DG-FACT-224`).
- **Icon invisible in the tab strip.** `getIcon()` returned a string
  instead of an `HTMLElement`, or `setIcon(...)` was called from the
  constructor. Return an `<img>` / `<i>` from `getIcon()` (`DG-FACT-215`).

## Verification

- TypeScript build (`webpack` or `grok check`) exits `0`.
- `grok publish` lists the function with `tags: view`; the regenerated
  `src/package.g.ts` retains all four header lines.
- Opening the deep-linked URL (e.g. `/dashboard/abc`) opens the tab
  pointed at `abc` — constructor sees `path='/dashboard/abc'`.
- Add the open screen to a project, save, log out, reopen — the same
  state restores (`saveStateMap` / `loadStateMap` round-trip).
- The tab strip shows your icon (`getIcon()` returned a node).

## See also

- Source articles: `help/develop/how-to/views/custom-views.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-213` … `DG-FACT-217`, plus `DG-FACT-224` (`* as DG` import
  convention).
- Reference packages:
  - `packages/Notebooks/src/package.js:27-128, 528-535` — canonical
    `tags: view` registration plus full state + path round-trip.
  - `packages/Flow/src/funcflow-view.ts:24` — `ViewBase` subclass
    constructed by an `meta.role: app` function instead.
  - `packages/PreclinicalCase/src/views/study-summary-view.ts:25` —
    minimal `ViewBase` that builds `this.root` lazily.
- Related skills: `routing` — wires URL state for `meta.role: app`
  functions; pair with this skill when the screen is opened FROM an app.
