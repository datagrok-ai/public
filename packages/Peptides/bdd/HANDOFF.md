# Peptides → bdd translation: handoff (2026-09-11, afternoon)

The task: translate the 14 hand-written Playwright specs of this package (`playwright/*.test.ts`,
their `.md` under `public/playwright-public/Peptides/`) into Gherkin features on
`@datagrok-libraries/bdd`, following `public/.claude/skills/bdd-translate/SKILL.md`. Nothing is
committed. This file is where the previous session stopped: environment prepared, surveys done,
decisions taken, live DOM probed, **no feature and no core/package code written yet**. Read this,
then the three survey reports in `docs/`, then start at "Plan" below.

The lead's standing rules (from earlier rounds, still binding): nothing committed without the
lead's order; a wait/sleep/pixel scan/retry in a test means a missing name or signal in the core —
fix the core; never restart the user's `pub serve`; never kill node processes wholesale; ≤3 agents
at a time; probe the DOM before naming anything; a step written for one viewer that another could
want goes to the library's `bindings/tiers/viewers/widgets.ts`.

## 1. Environment state (what was done, what runs)

- **Stand.** At the start of the session datlas (:8082) and pub serve (:63343) were not running
  (host nginx on :8888 answered 502). They were started with the user's own script
  `C:\Users\rizhi\Desktop\GROK-CORE\dg-start.bat` (two `cmd /k` console windows:
  `dart ./bin/server.dart -c start_full` in `core/server/datlas` with
  `GROK_SCRIPTING_MANAGE_CONTAINERS=false`, and `pub serve --port=63343` in `core/client/xamgle`).
  Both were up and the 32 MB client bundle compiled when the session ended. Dev-key login works
  through nginx only: `POST http://localhost:8888/api/users/login/dev/admin` (the direct :8082
  call answers "Invalid session"). Peptides, Bio, Dendrogram, EDA-family packages are published
  on the stand (versions unknown; `Peptides` was published 2026-09-11 10:50 with the rest).
- **Package install.** `package.json` gained `"@datagrok-libraries/bdd": "file:../../libraries/bdd"`
  (devDependencies) and the script `"test:bdd": "grok-bdd run"`; `npm install` ran (943 packages,
  lock refreshed); `npx grok-bdd link` made the package's `@playwright/test` the library's copy
  (1.61.1 — Node 18 CI pin). `node_modules/@datagrok-libraries/bio` is the npm 6.0.2 (no
  automation surface in it — irrelevant: the WebLogo viewer used at runtime is the **Bio
  package's**, `packages/Bio/src/viewers/web-logo-viewer.ts`, published on the stand).
- **Scaffold.** `npx grok-bdd init` created `bdd/` (`bdd.config.json` now `{"tiers": ["viewers"]}`,
  `package.json`, `tsconfig.json`, `bindings/elements.ts`, `bindings/steps.ts`,
  `features/smoke.feature`, `generated/smoke.test.ts`), `.vscode/settings.json` (Bio/DiffStudio
  commit theirs too) and the `.gitignore` lines for `bdd/test-results/`, `bdd/e2e/`,
  `bdd/.auth.json`. `tsconfig.json` `exclude` now includes `"bdd"` (by hand — the file has
  comments). The scaffold's `bindings/steps.ts` "user opens the Peptides app" cannot work: the
  landing view is not tagged `app` (see §4).
- **Verified:** `npx grok-bdd run --reporter=list` → the smoke feature passes (login, 4.9 s);
  `npx webpack` builds the package unchanged (17.9 s, 2 warnings). So the toolchain works end to
  end: edit `src/` → `npx webpack` → `grok publish localhost --key admin --skip-check` (the Bio
  README's recipe) → `npx grok-bdd compile` → `npx grok-bdd run`.
- **Git (public submodule, branch master, nothing staged):** modified `package.json`,
  `package-lock.json`, `.gitignore`, `tsconfig.json`; untracked `.vscode/`, `bdd/`. The core
  checkout is clean apart from the pre-existing untracked `core/server/datlas/dev-key.txt`.

## 2. Findings from the live stand (probe `docs/probe-sar.mjs`, its output in `docs/probe-sar-output.md`)

- `System:DemoFiles/bio/peptides.csv` (repo file `data/demo/bio/peptides.csv`, NOT under
  `public/`): 647 rows, columns `ID, AlignedSequence, IC50`; `AlignedSequence` is detected as
  `Macromolecule`, `units=separator`, renderer `sequence`. Every sequence is 17 `-`-separated
  tokens `NH2 … COOH`.
- The column's context panel panes: `pane-Details | pane-Filter | pane-Actions | pane-Colors |
  pane-Style | pane-Settings | pane-Plots | pane-Advanced | pane-Peptides | pane-Bioinformatics`
  (`grok.shell.o = column` was used by the probe; the features should click
  `the "header AlignedSequence" area of grid` — untested whether that makes the column current).
- `pane-Peptides` content names: `div-section--Peptides`, `input-host-Activity`/`input-Activity`,
  `input-host-Scaling`/`select[input-Scaling]`, `input-host-Clusters`, `input-Generate-clusters`,
  `button[button-Launch-SAR]`, `viewer-Histogram`, and the MCL form (hidden until the gear
  `Adjust clustering parameters` is clicked): `input-Similarity-Threshold`, `input-Inflation-Factor`,
  `input-Max-Iterations`, `input-Min-Cluster-Size`, `select[input-Distance-Function]`,
  `select[input-Fingerprint-Type]`, `input-Gap-Open-Penalty`, `input-Gap-Extend-Penalty`,
  `input-Use-WebGPU`. The pane's WebLogo (`.bio-wl-host`, 2 canvases) has **no `viewer-` name**
  (`DG.Widget.find` on its host's parent does find the Bio WebLogo viewer, with `getWidgetStatus`).
- Launch SAR from the pane: model in 4.7 s, all four viewers (`Sequence Variability Map`,
  `Most Potent Residues`, `MCL`, `Logo Summary Table`) at **38 s** on this stand. Settings after
  launch: `activityColumnName: IC50`, `activityScaling: none`, `mclSettings: {threshold: 70,
  inflation: 1.4, maxIterations: 16, minClusterSize: 5, distanceF: Needlemann-Wunsch, …}`.
  Columns added: `Activity`, `"1"…"17"`, `EmbedX (MCL)`, `EmbedY (MCL)`, `Cluster (MCL)`,
  `Cluster size (MCL)`, `Connectivity (MCL)`.
- `tv.viewers` after launch: `Grid`, `Sequence Variability Map`, `Most Potent Residues`, `MCL`,
  `Logo Summary Table`, `Grid` (the docked "Selection" grid). **Every JS viewer answers
  `typeof getWidgetStatus === 'function'` and `typeof isRenderPending === 'boolean'`** — that is
  the JS API base (`js-api/src/viewer.ts:144, 341` → Dart host: an EMPTY status and a host flag
  that is always false), so `settle` returns at once and "the viewer has no readings" is what a
  feature would see today. `onRendered` is undefined on all of them. Each must override all three.
- DOM roots (`[name^="viewer-"]`, with containment): `viewer-Grid` (main), `viewer-MCL` ⊃
  `viewer-Scatter-plot`, `viewer-Sequence-Variability-Map` ⊃ `viewer-Grid`,
  `viewer-Most-Potent-Residues` ⊃ `viewer-Grid`, `viewer-Logo-Summary-Table` ⊃ `viewer-Grid`,
  `viewer-Grid` (Selection). So **five elements are named `viewer-Grid`** — the reserved `grid`
  phrase is ambiguous after a launch (Playwright strict mode). The JS viewers' own `root` carries
  no name (the name is on the dock host around it); the runtime's `findViewer` still resolves the
  named element to the JS viewer by containment.
- The SVM's inner grid status (a Dart `Grid`): areas `cell <r> of AAR`, `cell <r> of 1`…`17`,
  `header <col>`; values `text of cell <r> of AAR` = the monomer (`A`, `C`, `COOH`, …), rows 22,
  columns `AAR, 1..17`, current column `∑` (the hidden total-count column). Rows are sorted by
  monomer, table rows counted from 1. LST inner grid: `cell <r> of Cluster`, `cell <r> of Members`.
- Ribbon: the settings wrench is `i.grok-icon.fal.fa-wrench` with
  `aria-label="Peptides analysis settings"` and no `name` → the `icon` kind matches by aria:
  `"Peptides analysis settings" icon`. SVM header inputs: `input-Mutation-Cliffs`,
  `input-Invariant-Map` (bool inputs retyped to radio), `input-Search`.
- Dock tab titles after launch: `Table | MCL | Sequence Variability Map | Most Potent Residues |
  Logo Summary Table | Selection`.
- Library facts checked: `locateActionable` keeps several visible matches ambiguous (strict);
  `kind('icon')` matches `aria`/`name`; `kind('section')` aliases `pane` with dart name
  `pane-{q}`; `kind('viewer')` dart name `viewer-{q}`; the d4 accordion pane header ALREADY sets
  `aria-expanded` (`accordion.dart:351, 365`) — survey A's "needs-core" for expanding a pane is
  wrong, `user expands Peptides pane` should work as is; `expectCustomEvent(page, id, timeoutMs)`
  takes a budget (the Then step's default is 30 s).

## 3. Decisions taken (execute unless the lead says otherwise)

1. **Pass 0 is the job**: every Peptides viewer gets `getWidgetStatus()`, `get isRenderPending()`,
   `get onRendered()` (pattern: `packages/Bio/src/viewers/web-logo-viewer.ts:1003-1041,
   1396-1406`; a Subject fed by the inner Dart viewer's render event). Delegate geometry to the
   inner Dart viewer's status and re-key it by meaning:
   - `SARViewer` (both subclasses): `parts.canvas/overlay` = the inner grid's; `isRenderPending` =
     inner grid pending OR mutation cliffs still computing (`calculateMutationCliffs` promise);
     `onRendered` from `grid.onAfterDrawContent` (subscribe when the grid is created in
     `createViewerGrid`; `_viewerGrid` is nulled and recreated on many property changes — hook
     every creation); values `activity scaling`, `data source`, `rows shown`, `positions`,
     `monomers`, `selected monomer-positions` (`"2:A, 4:Q"`, the model's own join
     `model.ts:456-463`), `cliff cells` / `cliff pairs` (omit both keys while cliffs are null),
     `error` (the "Please, select…" text when shown).
   - `MonomerPosition` (Sequence Variability Map): areas `cell <monomer> at <position>` from the
     grid's `cell <r> of <pos>` with monomer = `text of cell <r> of AAR`; `monomer <m>`; values
     `mode` (`Mutation Cliffs` / `Invariant Map`), per visible cell `count of cell <m> at <p>`,
     `mean difference of cell …`, `cliffs of cell …` from `monomerPositionStats` / `mutationCliffs`.
   - `MostPotentResidues`: areas `position <p>` = the `Diff` cell of the row (the only clickable
     column), `monomer of position <p>`; values `monomer at <p>`, `mean difference at <p>`,
     `p-value at <p>`, `count at <p>`, `positions`.
   - `LogoSummaryTable`: areas `cluster <name>` (Cluster cell), `weblogo of cluster <name>`,
     `distribution of cluster <name>`; values `clusters`, `clusters shown`, `members total`,
     `Members of cluster <n>`, `Mean difference of cluster <n>`, `P-Value of cluster <n>` (`''`
     for `DG.FLOAT_NULL`), `clusters column`, `selected clusters`, `message`; pending includes the
     per-row WebLogo creation cache and the 200 ms invalidate timer.
   - `MutationCliffsViewer`: delegate to the inner line chart + `position`, `cliff rows`,
     `message`; pending = the 300 ms debounce timer + line chart pending.
   - `ClusterMaxActivityViewer`: delegate to the inner scatter plot + `cluster size threshold`,
     `activity threshold`, `activity target`; pending = `renderTimeout != null` + scatter pending.
   - `SequencePositionStatsViewer`: delegate to the inner box plot + `positions`.
   - Forward `immediateRendering` to the inner Dart viewer when the library arms the wrapper.
2. **Names** (in package code): the inner grids get `name="<Type>-grid"` (e.g.
   `Sequence-Variability-Map-grid`, not `viewer-…`, so they stop colliding with `grid`); the
   docked Selection grid `name="viewer-Selection"` (`widgets/peptides.ts:353-361`; verify Dart
   does not re-set it); the pane/dialog WebLogo root `name="viewer-WebLogo"`
   (`widgets/peptides.ts`, both `fromType('WebLogo')` sites) so `WebLogo viewer in Peptides pane`
   resolves; the manual-alignment textarea gets the caption `Sequence`
   (`widgets/manual-alignment.ts:17`).
3. **Launch/settle signal**: `grok.events.fireCustomEvent('peptides-sar-ready', {table})` at the
   end of `startAnalysis` (`widgets/peptides.ts`, after the Selection grid is docked) and after a
   settings apply finishes its async work (`model.ts` `set settings`: collect the promises of
   `addMCLClusters`/`addClusterMaxActivityViewer`/`addDendrogram`/… and fire when all settle).
   Features: `Given user listens for "peptides-sar-ready" custom event` (library) before OK, then
   the package step `Then the SAR analysis should be ready` = `expectCustomEvent(page,
   'peptides-sar-ready', 180000)` (MCL is ~35 s on 647 rows here, 40–105 s on dev).
4. **Package steps** (`bdd/bindings/steps.ts`): `the Peptides package is initialized`
   (`grok.functions.call('Peptides:initPeptides')`, returns nothing — Bio's pattern); `the SAR
   analysis should be ready`; `the SAR setting {string} should be {string}` (reads the `settings`
   tag JSON, dotted path, e.g. `mclSettings.threshold`); `user opens the Peptides landing view`
   (`grok.functions.call('Peptides:Peptides')` + `grok.shell.addView`; the func is not tagged
   `app`, registering it as one is the lead's call); `only rows with {string} at position {int}
   of {string} column should be selected` (split on the column's separator) — or use the `"2"`
   position columns after a launch with `only rows where "2" is "A" should be selected`.
   Dataset: `dataset('peptides', {path: 'System:DemoFiles/bio/peptides.csv'})` in `elements.ts`.
5. **Core seam for the WebLogo header glyphs** (the one core change): the main grid's header
   glyphs are drawn by Peptides (`utils/cell-renderer.ts:254, 388` fills `model.webLogoBounds`
   in canvas CSS px) and clicked on the grid overlay (`:403-445`), but the grid's Dart status
   cannot know them. Add status providers: Dart `Widget` keeps a list of JS functions
   (`Widget_AddStatusProvider` in `xamgle/lib/src/interop/grok_api.dart`, next to
   `Widget_GetWidgetStatus` at :424, which merges each provider's `hitAreas`/`values` into the
   status JSON), JS API `Widget.addStatusProvider(fn)` in `js-api/src/widgets.ts` (+ the
   `grok_api.g.ts` line by hand; rebuild the js-api bundle into `xamgle/web/js/api`). Peptides
   registers one on `analysisView.grid` in `model.ts updateGrid()` publishing `<monomer> at
   <position>` areas from `webLogoBounds` and `highlighted rows`. Peptides builds against npm
   `datagrok-api` 1.26.0, so call it as `(grid as any).addStatusProvider?.(…)`. Features then say
   `user clicks on the "A at 2" area of grid` (+ `holding Shift` / `holding Control`).
6. **Library fix** (`libraries/bdd/src/runtime/viewer-runtime.ts findViewer`): the element named
   is the viewer meant — check an exact root match, then `DG.Widget.find(el)` when `el` itself is
   a `[name^="viewer-"]` root with a status, and only then containment. Needed for `scatter plot
   viewer in MCL viewer`, `line chart viewer in Sequence Mutation Cliffs viewer`, `Histogram
   viewer in Peptides pane`; without it those resolve to the outer JS wrapper.
7. **Product bugs found by the surveys** (verify on the stand, then fix or `@known-failure`):
   `widgets/manual-alignment.ts:26-29` writes split index `i` into column `i.toString()` but the
   position columns are named `i+1` (`libraries/bio/src/utils/splitter.ts:16`) — Apply shifts the
   row's cells left by one and never writes column "1"; `model.ts` `closeViewer(DENDROGRAM)`
   can never find the dendrogram (it is a GridNeighbor, not a viewer) so the toggle-off is a
   no-op and `settings.ts:139` cannot disable the toggle; the model re-asserts
   `grok.shell.o = acc.root` on a 1.5 s timer (`model.ts:911-914`).
8. **Out of scope, say so in the README**: the Dendrogram toggle (needs a name on the Dendrogram
   package's neighbor root + the bug above), the LST per-row WebLogo cells (Bio viewers inside
   cells), the Sequence space toggle (no old spec covers it), the demo-gallery route.
9. **Not restored** from the old specs: everything listed in the surveys' §4 — sleeps, synthetic
   `dispatchEvent`s, coordinate clicks, `getImageData` scans, `childCount > 5`, `lastError`
   regexes (`lastError` is a Promise: `String()` gives `[object Promise]`), private-state reads,
   `model.addClusterMaxActivityViewer()` as a "UI" path, `_doExportMutationCliffs([idCol])`
   instead of driving the dialog, "restored OR re-initialized" either/or claims.

## 4. Proposed features (≈11, all on `peptides`; group by launched analysis, `@journey`)

| feature | from old spec(s) | scenarios (short) |
|---|---|---|
| `panel/peptides-pane.feature` | info-panels, peptides | column renders as sequence; Details + Peptides panes; Details text; SAR parameters with IC50 default; scaling change rebuilds the histogram (`axis max` reading, root is replaced — fresh locate); WebLogo glyph click selects the rows (`T at 5` = 630 by naive split — verify numbering) |
| `sar/from-panel.feature` | sar | Launch SAR adds SVM+MPR; clustering adds MCL, `Cluster (MCL)`, LST; wrench → settings with MCL defaults (70/1.4); threshold 50 re-clusters; Invariant Map mode (`mode` reading + repaint); cell click selects rows + Mutation Cliffs pairs / Distribution panes (`present`); Distribution "Positions" split; back to Mutation Cliffs |
| `sar/from-top-menu.feature` | peptide-space, sar-viewer-lifecycle | Analyze Peptides dialog defaults; OK launches; settings panes General/Viewers/MCL; inflation 2.5 re-runs MCL (`scatter plot viewer in MCL viewer should be painted`); Active peptide selection toggles on / off with SVM+MPR+MCL staying |
| `sar/weblogo-selection.feature` | collaborative-selection | glyph click → `only rows where "2" is "A"` (299), Selection Sources `2:A`, Distribution/Selection panes; Shift adds `Q at 4` (314); Control removes it; clearing empties the panes |
| `sar/tooltips.feature` | monomer-position-hover-tooltip (first 100 rows) | IM cell tooltip shows `14 (`; another cell replaces it (`74 (`); MC cliff-cell tooltip shows `Pairs count`; header glyph hover; leaving hides tooltip + `highlighted rows` 0; settings round-trip keeps tooltips; LST cluster tooltip |
| `sar/mutation-cliffs.feature` | mutation-cliffs-compute-pipeline (first 200 rows) | `cliff cells`/`cliff pairs` ≥ 1, `count of cell A at 2` = 59; MC mode repaints; Sequence Mutation Cliffs viewer at position 2; export; LST members sum = 200, p-value in [0,1]; cluster click selects members |
| `sar/similarity-threshold.feature` | sar-similarity-threshold-matrix (first 200 rows, Outline 10/50/75/90) | gear reveals the threshold input; setting recorded (`the SAR setting "mclSettings.threshold"`); viewers attach; header ≥100 px and painted; IM cell click selects |
| `sar/export.feature` | export-invariant-map, export-mutation-cliffs | Invariant Map: 18 cols × 22 rows, `"1"` of NH2 = 647, `"2"` of A = 299; same from MPR; Mutation Cliffs: 6 columns, 6253 rows, `Mutation` = `Seq 1#Seq 2`, semtypes; extra column ID → 8 columns (drive the `ui.input.columns` picker — probe its dialog) |
| `sar/manual-alignment.feature` | manual-alignment | Monomer cell → Manual Alignment pane with the sequence; Apply rewrites the sequence and cells `"1"…"3"` (the off-by-one); Reset re-binds; Reset discards an unsaved edit; a glyph click still selects |
| `sar/project-round-trip.feature` | sar-save-reopen | save as project (library step, cleaned at feature end), close all, open: viewers, `settings` tag, 17 position columns, header height; the reopened analysis answers a cell click (one true claim about the selection, not either/or) |
| `entry/demo-dashboard.feature`, `entry/landing.feature` | peptide-sar-demo-dashboard | demo builds `Simple peptides` (647 rows, 15 positions, tags PT/fasta/SEQ.MSA, scaling `-lg`, MCL threshold 94); landing view: three demo buttons, Simple 647/fasta, Complex 540/`MSA` separator, HELM 334/`HELM` helm (all three tables are named `Peptides`, views `PeptidesView` — claim through the current table) |

Fixture counts (naive `-` split, token 1 = NH2; **verify on the stand whether column "1" holds
NH2 or the splitter drops terminals**): pos 5 `T` 630; pos 13 `N` 643; pos 3 `A` 535; `A` at 2
= 299 (first 100 rows 14, first 200 rows 59); `Q` at 4 = 116; `A@2 ∪ Q@4` = 314; `N` at 4 = 474
(first 100: 74). Invariant map rows: 22 monomers (`A,C,COOH,D,E,F,G,H,I,K,L,M,N,NH2,P,Q,R,S,T,V,W,Y`).

## 5. Plan (order that worked in the viewer rounds)

1. Core seam (§3.5) + library fix (§3.6): edit, `npm run build` in `libraries/bdd`, rebuild the
   js-api bundle, touch the Dart files and ping `http://localhost:8888/login.dart.js_1.part.js`
   (a recompile never takes more than 4 min; verify by grepping a string literal you added out
   of `login.dart.js_1.part.js`, not by etag).
2. Package: statuses (§3.1), names (§3.2), event (§3.3), the manual-alignment fix (§3.7 — after
   confirming it). `npx webpack`, `grok publish localhost --key admin --skip-check`. Re-run
   `docs/probe-sar.mjs` (extend it) to read the new statuses and answer the open questions in the
   surveys' §5: position numbering, which cells carry cliffs on 100/200 rows, cluster count on
   200 rows, the `ui.input.columns` picker DOM, how `readValue` reads a Dart column input
   (`editor of Activity input … should have text "IC50"` is Bio's phrase), the tooltip root class,
   the selection-border hex, whether `header AlignedSequence` click makes the column current.
3. Bindings (`bdd/bindings/elements.ts`, `steps.ts`, one file per subject as Bio does), then the
   features in §4, `npx grok-bdd compile`, `npx grok-bdd run --reporter=list generated/<folder>`;
   green twice, then `PLAYWRIGHT_WORKERS=2` (4 workers collapse this stand), then headed once.
4. Pass 2: one read-only reviewer per old-spec/feature pair (prompt in the skill), consolidate,
   fix; then the package README (Bio's as the template: what each folder claims, how to run,
   timings from a JSON-reporter run), the library CLAUDE.md facts, the surveys' stale-md notes,
   and CHANGELOG `## v.next` lines for the package changes.
5. Leave the old `playwright/` specs and md files in place unless the lead orders otherwise
   (Bio's precedent).

## 6. Files in this folder

- `docs/survey-launch-panels.md` — spec group A (peptides, sar, info-panels, peptide-space,
  sar-viewer-lifecycle). One correction: the accordion already has `aria-expanded`.
- `docs/survey-selection-tooltips-pipeline.md` — group B (collaborative-selection,
  monomer-position-hover-tooltip, mutation-cliffs-compute-pipeline, sar-similarity-threshold-matrix).
- `docs/survey-export-alignment-project-demo.md` — group C (export-invariant-map,
  export-mutation-cliffs, manual-alignment, sar-save-reopen, peptide-sar-demo-dashboard).
- `docs/probe-sar.mjs` — the live probe (Playwright via the library's copy, the storage state
  the smoke run wrote to `bdd/e2e/.auth.json`); `docs/probe-sar-output.md` — what it printed.
- `features/smoke.feature`, `generated/smoke.test.ts` — the scaffold's smoke test (passes;
  delete when the real features exist). `bindings/*.ts` — scaffold placeholders to replace.
