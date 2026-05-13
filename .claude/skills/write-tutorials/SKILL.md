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

1. **Create the source file and subclass `Tutorial`.** Convention — `src/tutorials/<name>.ts`. Import the full subpath (see DG-FACT-434). Override `demoTable = ''` to skip the auto-loaded `'demog.csv'` (see DG-FACT-436). Override `prerequisites` to gate on packages/services (see DG-FACT-439).
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
   `steps` should equal the number of `action(...)` calls inside `_run` — it drives the progress bar.

2. **Fill `_run()` with `action(...)` calls keyed off real events.** `action(instructions, completed, hint?, description?)` advances when `completed` (Observable or Promise) resolves (see DG-FACT-435). Use higher-level helpers (`openPlot`, `openDialog`, `dlgInputAction`, `buttonClickAction`) where they fit.
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
   Reference: `packages/Tutorials/src/tracks/eda/tutorials/filters.ts:65-79`.

3. **Register the tutorial in `src/package.ts`.** Required: `//tags: tutorial`, `//meta.name`, `//output: object`. Optional: `//meta.track`, `//meta.icon`, `//description` (see DG-FACT-437).
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

4. **(Optional) Register a track for a custom `help-url`.** The Tutorials package does NOT invoke the function body — it builds an empty Track from metadata (see DG-FACT-438). Skip if you don't need a help link.
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
   Open **Apps → Tutorials** — your track shows the tutorial card; clicking it walks the steps as users interact with the real UI.

## Common failure modes

- **`Cannot find module '@datagrok-libraries/tutorials/src/tutorial'`.** Full subpath mandatory; re-run `npm i` (see DG-FACT-434).
- **Card never appears in Apps → Tutorials.** Missing `//tags: tutorial` or `//output: object`, function not exported, or Tutorials package not installed.
- **Step never completes.** `completed` observable predicate never matches — log candidate events and align (see DG-FACT-435).
- **Lesson aborts with "Please install package…".** Unmet `prerequisites`; install or drop the optional prereq (see DG-FACT-439).
- **Tutorial sees `demog.csv` already loaded.** Set `demoTable = ''` (see DG-FACT-436).
- **Track help link goes nowhere.** Add a `tags: track` function with `//help-url` (see DG-FACT-438).

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
