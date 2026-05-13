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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# custom-views

## When to use

You need a workspace screen that isn't a `TableView` — a dashboard,
notebook host, study summary — that deep-links via URL, restores from a
saved project, and carries a package-owned icon. For a one-shot tab
with no URL/project needs, `grok.shell.newView` is enough — see Step 1.

## Prerequisites

- Familiarity with the `routing` skill if you'll wire URL state.

## Steps

1. **Decide which API tier you need.**
   `grok.shell.newView` is the one-liner — fine for production, but
   the returned `View` does NOT dispatch URL paths or serialize state
   (see `DG-FACT-213`). Subclass `DG.ViewBase` only when you need
   routing, project restore, or a per-class icon (see `DG-FACT-214`):
   ```typescript
   // Ad-hoc: one-liner, no routing/state — keep for prototypes & dialogs
   const v = grok.shell.newView('Quick View', [ui.divText('Hi!')]);
   ```
   Stop here if you don't need routing or state.

2. **Subclass `DG.ViewBase` with the standard constructor signature.**
   Use the `DG.` prefix (`DG-FACT-224`). Never pass `createHost=false`
   from package code — that path is reserved for the internal `View`
   subclass (see `DG-FACT-214`):
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
   Store identity on `this`, expose it via `path`, serialize via
   `saveStateMap` / `loadStateMap`, dispatch URLs via `handlePath` /
   `acceptsPath` (see `DG-FACT-215` for the full overridable surface):
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
   The platform reads four exact header lines (see `DG-FACT-216`).
   Without `tags: view` the function is callable but never URL-dispatched:
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
   Notebooks (`packages/Notebooks/src/package.js:528-535`) is the only
   production reference today (see `DG-FACT-217`).

5. **Build, publish, and exercise the deep-link URL.**
   ```bash
   webpack && grok publish
   ```
   Opening `https://<host>/dashboard/abc` instantiates `DashboardView`
   with `path='/dashboard/abc'`; `handlePath` runs and the tab opens
   pointed at `abc`.

## Common failure modes

- URL doesn't change / deep-link doesn't restore — you used
  `grok.shell.newView` (no path dispatcher; `DG-FACT-213`). Subclass
  `DG.ViewBase` and register via `tags: view`.
- Reopened project shows empty tab — state lives in module-level vars.
  Move identity onto `this` and round-trip via state map (`DG-FACT-215`).
- Function compiles but path never dispatched — missing `tags: view` or
  renamed `input: string path` header (`DG-FACT-216`).
- Type error on `extends ViewBase` — use `DG.ViewBase` with the default
  `* as DG` import (`DG-FACT-224`).
- Icon invisible — `getIcon()` returned a string, not `HTMLElement`
  (`DG-FACT-215`).

## See also

- Source articles: `help/develop/how-to/views/custom-views.md`
- Knowledge: facts `DG-FACT-213`…`DG-FACT-217`, `DG-FACT-224`.
- Reference packages: `packages/Notebooks/src/package.js:27-128, 528-535`
  (canonical `tags: view`); `packages/Flow/src/funcflow-view.ts:24`;
  `packages/PreclinicalCase/src/views/study-summary-view.ts:25`.
- Related skills: `routing` — pair when the screen is opened FROM a
  `meta.role: app` function.
