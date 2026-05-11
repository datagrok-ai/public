---
name: write-demo-script
description: Write a Datagrok demo script that walks users through package functionality step-by-step in the Demo app
harness-authored: true
---

# write-demo-script

## When to use

You want a stepwise, narrated walkthrough of a package feature to appear
under **Apps → Tutorials → Demo → <Group> → <Your demo>**. Triggers:
"add a demo to the gallery," "make a tutorial that runs itself,"
"register a demo script," "guide the user through these viewers." For a
static dashboard (no step list), use `meta.isDemoDashboard: true`
instead — out of scope here (`packages/Chem/src/package.g.ts:1080-1095`).

## Prerequisites

- A package scaffold (`grok create <Name>`); functions you want to
  demonstrate already callable.
- Add the tutorials library to the package: in `package.json`,
  `"@datagrok-libraries/tutorials": "^<x>.<y>.<z>"` under
  `dependencies`, then `npm i` (knowledge `DG-FACT-374`).
- Familiarity with package header annotations (`//name:`, `//meta.*:`)
  or the `@grok.decorators.func` form.

## Steps

1. **Create the demo source file and import.**
   Convention — `src/demo/<feature>.ts`
   (`packages/Bio/src/demo/bio03-atomic-level.ts:1-75`,
   `packages/Chem/src/demo/demo.ts:1-90`). Import path fixed
   (knowledge `DG-FACT-374`).
   ```bash
   mkdir -p src/demo && touch src/demo/my-feature-demo.ts
   ```
   ```typescript
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
   ```
   Expected: `npx tsc --noEmit` resolves the import. If it can't, the
   dependency from Prerequisites didn't land in `node_modules/`.

2. **Instantiate `DemoScript` with name, description, mode, options.**
   Constructor signature (knowledge `DG-FACT-375`):
   `new DemoScript(name, description, isAutomatic = false,
   options?: {autoStartFirstStep?: boolean, path?: string})`.
   ```typescript
   export async function demoMyFeature(): Promise<void> {
     const script = new DemoScript(
       'My Feature',                            // name
       'Walks through My Feature end-to-end.',  // description
       false,                                   // isAutomatic — manual: user clicks "Next"
       {autoStartFirstStep: true, path: 'MyArea/MyFeature'},
     );
     // steps added below
   }
   ```
   Pick `isAutomatic`:
   - `true` → the script auto-advances using each step's `delay`
     (default 2000 ms, knowledge `DG-FACT-377`). `autoStartFirstStep`
     is silently ignored in this branch (drift
     `DG-FACT-DRIFT-DEMO-002`).
   - `false` → manual; user clicks "Next step." `delay` values are
     parsed but inert (knowledge `DG-FACT-378`). Set
     `autoStartFirstStep: true` to kick off step 1 without a click.

   `options.path` becomes the in-view breadcrumb suffix and the
   `grok.shell.v.path` value, appended to `apps/Tutorials/Demo/`
   (knowledge `DG-FACT-379`). Spaces in `path` are replaced with `-`.

3. **Chain `.step(...)` calls.**
   Signature (knowledge `DG-FACT-376`): `step(name, async () => {...},
   {description?, delay?}): this`. The fluent API means you keep
   appending. `func` must be `async` returning `void`.
   ```typescript
   await script
     .step('Open dataset', async () => {
       const df = grok.data.demo.demog();
       grok.shell.addTableView(df);
     }, {description: 'Load the demog dataset.', delay: 2000})
     .step('Show a scatter plot', async () => {
       grok.shell.tv.addViewer(DG.VIEWER.SCATTER_PLOT);
     }, {description: 'Add a scatter plot viewer.', delay: 2000})
     .step('Finish', async () => {
       grok.shell.info('Demo complete.');
     })
     .start();
   ```
   Expected: `npx tsc --noEmit` clean. The terminal `.start()` is
   mandatory — without `await ...start()` nothing executes
   (knowledge `DG-FACT-376`).

4. **(Optional) wrap in try/catch.**
   Production demos handle exceptions so a broken step doesn't orphan
   the view (`packages/Bio/src/demo/bio03-atomic-level.ts:31-74`):
   ```typescript
   try { await new DemoScript(/* ... */).step(/* ... */).start(); }
   catch (err: any) { grok.shell.error(`Demo failed: ${err?.message ?? err}`); }
   ```

5. **Register the demo as a package function.**
   Add a header-annotated wrapper in `src/package.ts` (or the
   decorator form below). `meta.demoPath` is what makes it appear in
   the Demo gallery; `meta.isDemoScript: true` renders the stepwise
   sidebar (knowledge `DG-FACT-380`). Write the boolean lowercase —
   the article shows `True` but every in-tree usage is `true` (drift
   `DG-FACT-DRIFT-DEMO-001`).
   ```typescript
   // src/package.ts
   import {demoMyFeature} from './demo/my-feature-demo';

   //name: Demo My Feature
   //description: Walks through My Feature end-to-end.
   //meta.demoPath: MyArea | My Feature
   //meta.isDemoScript: true
   export async function demoMyFeatureWrapper(): Promise<void> {
     await demoMyFeature();
   }
   ```
   Decorator equivalent (codegen regenerates `package.g.ts`):
   ```typescript
   @grok.decorators.func({
     name: 'Demo My Feature',
     description: 'Walks through My Feature end-to-end.',
     meta: {demoPath: 'MyArea | My Feature', isDemoScript: 'true'},
   })
   static async demoMyFeatureWrapper(): Promise<void> {
     await demoMyFeature();
   }
   ```
   Do **not** write `return new demoMyFeature();` — the article's
   example uses `new` on a plain function (drift
   `DG-FACT-DRIFT-DEMO-003`); production code calls or awaits the
   function directly. Reference —
   `packages/Chem/src/package.g.ts:1057-1063`.

6. **Build, publish, and open the demo.**
   ```bash
   webpack
   grok publish dev
   ```
   Expected: `src/package.g.ts` contains your `//meta.demoPath:` and
   `//meta.isDemoScript: true` headers. Then in Datagrok navigate to
   **Apps → Tutorials → Demo → MyArea → My Feature** — the sidebar
   shows your step list with play/next controls.

## Common failure modes

- **Demo doesn't appear in the gallery.** `meta.demoPath` is missing
  or the function isn't exported (knowledge `DG-FACT-380`). Fix: add
  `//meta.demoPath: <Group> | <Sub>` and confirm `package.g.ts` has
  the line after rebuild.
- **Demo appears but shows a static "dashboard" page, not steps.**
  `meta.isDemoScript: true` is missing (or `meta.isDemoDashboard:
  true` is set, which overrides). Fix: add the lowercase
  `//meta.isDemoScript: true` annotation; remove any
  `meta.isDemoDashboard`.
- **Steps don't run; the dialog opens but the sidebar is empty.**
  You called `.step(...)` but forgot the terminating `.start()`, or
  awaited the constructor instead (knowledge `DG-FACT-376`). Fix:
  `await script.start()` — note `start()` returns a `Promise`.
- **`isAutomatic: true` ignores my `autoStartFirstStep`.** Expected —
  the field is consulted only in the manual branch (drift
  `DG-FACT-DRIFT-DEMO-002`). Choose one mode; don't combine.
- **`delay` has no effect on a manual demo.** Expected — `delay` is
  used only by the automatic countdown (knowledge `DG-FACT-378`). To
  make a manual demo time itself, switch `isAutomatic` to `true`.
- **`Cannot find module '@datagrok-libraries/tutorials/...'`.** The
  dependency wasn't added or `npm i` wasn't rerun (knowledge
  `DG-FACT-374`). Fix: add the entry to `package.json` `dependencies`
  and reinstall.

## Verification

- After `grok publish dev`, open
  `apps/Tutorials/Demo/<your-demoPath-joined-with-/>` in the URL bar.
  The sidebar should show your steps in order; clicking **Next**
  (manual) or watching the countdown (automatic) executes each step
  and updates `grok.shell.v.path` to
  `apps/Tutorials/Demo/<path-with-spaces-replaced-by-dashes>`
  (knowledge `DG-FACT-379`).
- The Demo app's gallery card for your demo carries the
  `description` text from the package function header.

## See also

- Source articles:
  - `help/develop/how-to/misc/write-demo-scripts.md`
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts
    `DG-FACT-374..380` and drifts
    `DG-FACT-DRIFT-DEMO-001..003`.
- Related skills:
  - `publish-packages` — `webpack && grok publish dev` flow that
    propagates decorator-form metadata into `package.g.ts`.
  - `manipulate-viewers` — common building blocks used inside step
    callbacks (`grok.shell.addTableView`, `tv.addViewer(...)`).
