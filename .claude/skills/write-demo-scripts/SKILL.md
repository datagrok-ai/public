---
name: write-demo-scripts
version: 0.1.0
description: |
  Wire a stepwise, narrated walkthrough of a package feature so it
  appears as an interactive entry under **Apps → Tutorials → Demo**.
  For plugin authors who want users to click one card and see a
  scripted sequence — open a dataset, run an analysis, add a viewer —
  with descriptions and (optional) auto-advance between steps.
  Use when asked to "ship a click-to-run feature walkthrough", "add
  an interactive showcase to the gallery", or "sequence a stepwise
  narrated example with timing between stages".
triggers:
  - interactive walkthrough in the gallery
  - click-to-run feature showcase
  - stepwise narrated example
  - guided package walkthrough
  - showcase a feature with timed steps
  - breadcrumb entry under Apps → Tutorials
allowed-tools:
  - Read
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# write-demo-scripts

## When to use

You want a stepwise, narrated walkthrough of a package feature to appear
under **Apps → Tutorials → Demo → \<Group\> → \<Your entry\>**, with
either auto-advance (timed) or manual ("Next step") progression. For a
static dashboard with no step list, use `meta.isDemoDashboard: true`
instead (out of scope).

## Prerequisites

- A package scaffold (`grok create <Name>`); the functions you want to
  showcase already callable from the package.
- Familiarity with package header annotations (`//name:`, `//meta.*:`)
  or the `@grok.decorators.func` form.

## Steps

1. **Create the walkthrough source file and import `DemoScript`.**
   Convention — `src/demo/<feature>.ts` (e.g.,
   `packages/Chem/src/demo/demo.ts:6`). The import path is fixed —
   the library declares no package root entry, so the full
   `src/demo-script` subpath is required (knowledge `DG-FACT-374`).
   ```bash
   mkdir -p src/demo && touch src/demo/my-feature.ts
   ```
   ```typescript
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
   ```
   Expected: `npx tsc --noEmit` resolves the import. If it can't, the
   dependency from Prerequisites didn't land in `node_modules/`.

2. **Instantiate `DemoScript` with name, description, mode, options.**
   Constructor (knowledge `DG-FACT-375`):
   `new DemoScript(name, description, isAutomatic = false,
   options?: {autoStartFirstStep?: boolean, path?: string})`.
   ```typescript
   export async function demoMyFeature(): Promise<void> {
     const script = new DemoScript(
       'My Feature',                            // name
       'Walks through My Feature end-to-end.',  // description
       false,                                   // isAutomatic — user clicks "Next"
       {autoStartFirstStep: true, path: 'MyArea/MyFeature'},
     );
     // steps added below
   }
   ```
   Pick `isAutomatic`:
   - `true` → auto-advance using each step's `delay` (default 2000 ms,
     knowledge `DG-FACT-377`). `autoStartFirstStep` is ignored.
   - `false` → manual; user clicks "Next step". `delay` values are
     parsed but inert (knowledge `DG-FACT-378`). Set
     `autoStartFirstStep: true` to start step 1 without a click.

   `options.path` becomes the in-view breadcrumb suffix and the
   `grok.shell.v.path` value, appended to `apps/Tutorials/Demo/`
   (knowledge `DG-FACT-379`). Spaces in `path` are replaced with `-`.

3. **Chain `.step(...)` calls and terminate with `.start()`.**
   Signature (knowledge `DG-FACT-376`): `step(name, async () => {...},
   {description?, delay?}): this`. The fluent return means steps chain.
   `func` must be `async` returning `void`; `await ...start()` is the
   trigger — without it nothing runs.
   ```typescript
   await script
     .step('Open dataset', async () => {
       grok.shell.addTableView(grok.data.demo.demog());
     }, {description: 'Load the demog dataset.', delay: 2000})
     .step('Show a scatter plot', async () => {
       grok.shell.tv.addViewer(DG.VIEWER.SCATTER_PLOT);
     }, {description: 'Add a scatter plot viewer.', delay: 2000})
     .step('Finish', async () => {
       grok.shell.info('Walkthrough complete.');
     })
     .start();
   ```
   Expected: `npx tsc --noEmit` clean.

4. **(Optional) wrap the chain in try/catch.**
   So a broken step doesn't orphan the view:
   ```typescript
   try { await script.step(/* ... */).start(); }
   catch (err: any) { grok.shell.error(`Walkthrough failed: ${err?.message ?? err}`); }
   ```

5. **Register the function in `src/package.ts` with `meta.demoPath`
   and `meta.isDemoScript: true`.**
   `meta.demoPath` is what makes the entry appear in the gallery;
   `meta.isDemoScript: true` renders the stepwise sidebar (knowledge
   `DG-FACT-380`). Write the value lowercase `true` — the article
   shows `True`, but every in-tree usage is lowercase (e.g.,
   `packages/Chem/src/package.g.ts:1059`); live execution wins.
   ```typescript
   // src/package.ts
   import {demoMyFeature} from './demo/my-feature';

   //name: Demo My Feature
   //description: Walks through My Feature end-to-end.
   //meta.demoPath: MyArea | My Feature
   //meta.isDemoScript: true
   export async function demoMyFeatureWrapper(): Promise<void> {
     await demoMyFeature();
   }
   ```
   The `@grok.decorators.func({meta: {demoPath: '...', isDemoScript:
   'true'}})` decorator form is equivalent — codegen emits the same
   header into `package.g.ts`. Reference —
   `packages/Chem/src/package.g.ts:1057-1063`.

6. **Build, publish, and open the entry.**
   ```bash
   webpack
   grok publish dev
   ```
   Expected: `src/package.g.ts` contains your `//meta.demoPath:` and
   `//meta.isDemoScript: true` headers. In Datagrok, navigate to
   **Apps → Tutorials → Demo → MyArea → My Feature** — the sidebar
   shows your step list with play/next controls.

## Common failure modes

- **Entry doesn't appear in the gallery.** `meta.demoPath` is missing
  or the function isn't exported (knowledge `DG-FACT-380`). Fix: add
  `//meta.demoPath: <Group> | <Sub>` and confirm `package.g.ts` has
  the line after rebuild.
- **Entry appears but shows a static "dashboard" page, not steps.**
  `meta.isDemoScript: true` is missing (or `meta.isDemoDashboard:
  true` is set, which overrides). Fix: add the lowercase
  `//meta.isDemoScript: true` annotation; remove any
  `meta.isDemoDashboard`.
- **Steps don't run; the dialog opens but the sidebar is empty.**
  You called `.step(...)` but forgot the terminating `.start()`, or
  awaited the constructor instead (knowledge `DG-FACT-376`). Fix:
  `await script.start()`.
- **Mode-option mismatch.** `autoStartFirstStep` is consulted only in
  manual mode; `delay` is consumed only by the automatic countdown
  (knowledge `DG-FACT-378`). Choose one mode; don't combine.
- **`Cannot find module '@datagrok-libraries/tutorials/...'`.** The
  dependency wasn't added or `npm i` wasn't rerun (knowledge
  `DG-FACT-374`). Fix: add the entry to `package.json` `dependencies`
  and reinstall.

## See also

- Source articles:
  - `help/develop/how-to/misc/write-demo-scripts.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-374..380` (import, constructor, fluent step API,
  delay semantics, `DEMO_PATH` breadcrumb, registration metadata).
- Related skills:
  - `publish-packages` — `webpack && grok publish dev` flow that
    propagates decorator-form metadata into `package.g.ts`.
  - `manipulate-viewers` — building blocks used inside step
    callbacks (`grok.shell.addTableView`, `tv.addViewer(...)`).
  - `write-tutorials` — fully interactive learning lessons (users
    drive the real UI with hints) rather than scripted playback.
