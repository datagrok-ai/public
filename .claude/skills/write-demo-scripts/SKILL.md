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

1. **Create the walkthrough source file and import `DemoScript`**
   (`DG-FACT-374` — full subpath required).
   ```bash
   mkdir -p src/demo && touch src/demo/my-feature.ts
   ```
   ```typescript
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
   ```

2. **Instantiate `DemoScript` with name, description, mode, options**
   (constructor at `DG-FACT-375`).
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
   - `true` → auto-advance per-step `delay` (default 2000 ms,
     `DG-FACT-377`). `autoStartFirstStep` is ignored.
   - `false` → manual ("Next step" button). `delay` values are parsed
     but inert (`DG-FACT-378`). Set `autoStartFirstStep: true` to start
     step 1 without a click.

   `options.path` becomes the breadcrumb suffix on `apps/Tutorials/Demo/`;
   spaces replaced with `-` (`DG-FACT-379`).

3. **Chain `.step(...)` calls and terminate with `.start()`** —
   fluent; `func` must be `async`; nothing runs without `await
   ...start()` (`DG-FACT-376`).
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

4. **(Optional) wrap the chain in try/catch.**
   ```typescript
   try { await script.step(/* ... */).start(); }
   catch (err: any) { grok.shell.error(`Walkthrough failed: ${err?.message ?? err}`); }
   ```

5. **Register the function in `src/package.ts` with `meta.demoPath`
   and `meta.isDemoScript: true`** (`DG-FACT-380`). Use lowercase
   `true` (article incorrectly shows `True`).
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
   `@grok.decorators.func({meta: {demoPath: '...', isDemoScript:
   'true'}})` is equivalent.

6. **Build, publish, and open the entry.**
   ```bash
   webpack
   grok publish dev
   ```
   Navigate to **Apps → Tutorials → Demo → MyArea → My Feature**.

## Common failure modes

- Entry not in gallery — `meta.demoPath` missing or function not exported (`DG-FACT-380`).
- Static dashboard instead of steps — missing `meta.isDemoScript: true` or `meta.isDemoDashboard: true` overrides.
- Sidebar empty — missing `.start()` call (`DG-FACT-376`).
- Mode-option mismatch — `autoStartFirstStep` is manual-only; `delay` is automatic-only (`DG-FACT-378`).
- `Cannot find module '@datagrok-libraries/tutorials/...'` — missing dep, run `npm i` (`DG-FACT-374`).

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
