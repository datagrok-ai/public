---
name: write-tutorials
version: 0.1.0
description: |
  Build an interactive in-app lesson — a `Tutorial` subclass that drives
  the real Datagrok UI with hint blobs on real widgets and waits for
  user actions to advance. For plugin authors who want their package
  to teach a domain workflow (open this file, filter that column, add
  this viewer) and have it appear in **Apps → Tutorials** grouped
  under a learning path.
  Use when asked to "ship interactive in-app onboarding", "wire a
  guided lesson that hints at real widgets and waits for clicks", or
  "build a step-by-step learning path inside the plugin".
triggers:
  - interactive in-app onboarding
  - guided lesson with hints on real widgets
  - teach a workflow inside the plugin
  - step-by-step learning path
  - in-app guided walkthrough that waits for clicks
  - group lessons under a learning track
allowed-tools:
  - Read
  - Edit
  - Write
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# write-tutorials

## When to use

You want users to learn a domain workflow by *doing it* — hint
indicators decorate real Datagrok widgets, the lesson advances on real
user events (filter change, dialog open, viewer added), and a sequence
of such lessons groups under a track. For scripted playback users only
watch, use `write-demo-scripts` instead.

## Prerequisites

- The Tutorials package installed on the target server — it owns the
  **Apps → Tutorials** view that hosts external tutorials (knowledge
  `DG-FACT-437`).

## Steps

1. **Create the source file and subclass `Tutorial`.**
   Convention — `src/tutorials/<name>.ts`. The import path is the full
   subpath; the library has no root entry (knowledge `DG-FACT-434`).
   Override `demoTable` if you don't want `'demog.csv'` auto-loaded
   into `this.t!` BEFORE `_run` (knowledge `DG-FACT-436`). Override
   `prerequisites` to gate on packages/services (knowledge
   `DG-FACT-439`); the runner aborts with a toast if any fail.
   ```typescript
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   import {Tutorial, TutorialPrerequisites} from
     '@datagrok-libraries/tutorials/src/tutorial';

   export class MyFeatureTutorial extends Tutorial {
     get name(): string { return 'My Feature'; }
     get description(): string { return 'Walks through My Feature.'; }
     get icon(): string { return '📘'; }
     get steps(): number { return 3; }

     demoTable: string = '';  // '' skips auto-load; default is 'demog.csv'
     prerequisites: TutorialPrerequisites = {grokConnect: true};

     protected async _run(): Promise<void> { /* see step 2 */ }
   }
   ```
   Expected: `npx tsc --noEmit` clean. `steps` should equal the number
   of `action(...)` calls inside `_run` — it drives the progress bar.

2. **Fill `_run()` with `action(...)` calls keyed off real events.**
   `action(instructions, completed, hint?, description?)` advances
   only when `completed` (Observable or Promise) resolves (knowledge
   `DG-FACT-435`). The `completed` source is what couples the step to
   the user actually doing the thing — never a fixed timeout.
   ```typescript
   import {filter} from 'rxjs/operators';

   protected async _run(): Promise<void> {
     this.header.textContent = this.name;
     this.title('Filter the table');
     this.describe('Narrow the dataset by a single category.');

     await this.action(
       'Click on a value in the "DIS_POP" filter',
       this.t!.onFilterChanged.pipe(filter(() => {
         const fs = this.t!.rows.filters;
         return fs.length === 1 && fs.get(0).startsWith('DIS_POP:');
       })),
       null,
       'Click <b>AS</b> in the DIS_POP filter.',
     );

     await this.openPlot('scatter plot',
       (v) => v.type === DG.VIEWER.SCATTER_PLOT);
   }
   ```
   Prefer the higher-level helpers where they fit — `openPlot`,
   `openDialog`, `dlgInputAction`, `buttonClickAction` — they encode
   the right event source. Reference —
   `packages/Tutorials/src/tracks/eda/tutorials/filters.ts:65-79`.

3. **Register the tutorial in `src/package.ts`.**
   Required: `//tags: tutorial`, `//meta.name`, `//output: object`.
   Optional: `//meta.track`, `//meta.icon`, `//description`
   (knowledge `DG-FACT-437`). Without `meta.track` the tutorial lands
   under a track named after the package's friendly name; with it
   set to an existing track name, it joins that track.
   ```typescript
   import {MyFeatureTutorial} from './tutorials/my-feature';

   //tags: tutorial
   //meta.name: My Feature
   //meta.track: My Plugin Onboarding
   //meta.icon: images/my-feature-tutorial.png
   //description: Walks through My Feature end-to-end.
   //output: object
   export function myFeatureTutorial() {
     return new MyFeatureTutorial();
   }
   ```

4. **(Optional) Register a track for a custom `help-url`.**
   `//tags: track`, `//meta.name`, `//help-url`, `//output: object`.
   Caveat — the Tutorials package does NOT invoke the function body;
   it builds an empty Track from metadata (knowledge `DG-FACT-438`).
   Skip this step if you're fine with `meta.track` creating the track
   implicitly without a help link.
   ```typescript
   //tags: track
   //meta.name: My Plugin Onboarding
   //help-url: https://datagrok.ai/help/my-plugin
   //output: object
   export function myPluginOnboardingTrack() { return null as any; }
   ```

5. **Build, publish, and open the entry.**
   ```bash
   webpack
   grok publish dev
   ```
   Expected: `src/package.g.ts` carries the `tags: tutorial` /
   `meta.name:` header. In Datagrok open **Apps → Tutorials** — your
   track shows your tutorial card; clicking it walks the steps as
   users interact with the real UI.

## Common failure modes

- **`Cannot find module '@datagrok-libraries/tutorials/src/tutorial'`.**
  Dependency missing or `npm i` not rerun. The full subpath is
  mandatory — no root entry (knowledge `DG-FACT-434`).
- **Card never appears in Apps → Tutorials.** Missing `//tags:
  tutorial` or `//output: object`, function not exported, or the
  Tutorials package isn't installed on the server. Confirm
  `package.g.ts` carries the header after `webpack`.
- **Step never completes — the indicator stays unchecked.** The
  `completed` observable filter never matches. Log candidate events
  (`this.t!.onFilterChanged.subscribe(console.log)`) and align the
  predicate; prefer `buttonClickAction` over hand-rolled DOM polling
  (knowledge `DG-FACT-435`).
- **Lesson aborts with "Please install package …" or "Service not
  available".** A `prerequisites` entry isn't met; install/enable or
  drop the optional prereq (knowledge `DG-FACT-439`).
- **Tutorial sees `demog.csv` already loaded.** Base class
  auto-loads it (knowledge `DG-FACT-436`). Set `demoTable = ''`.
- **Track help link goes nowhere.** No `tags: track` function with
  `//help-url` is registered (knowledge `DG-FACT-438`); add one.

## See also

- Source articles:
  - `help/develop/how-to/misc/write-tutorials.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-434..439` (subclass contract, `action` signature,
  auto-loaded demo table, `tags: tutorial` registration,
  `tags: track` registration drift, `TutorialPrerequisites`).
- Related skills:
  - `write-demo-scripts` — scripted playback (users watch
    `DemoScript.step(...)`) versus this skill's interactive lesson.
  - `manipulate-viewers` — building blocks called inside `action`
    callbacks (`grok.shell.addTableView`, `tv.addViewer(...)`).
