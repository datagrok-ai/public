---
name: bdd-translate
description: Translate hand-written Playwright specs (or any existing UI tests) into Gherkin features on @datagrok-libraries/bdd, then prove the features test what they claim by backward-matching them against the originals with independent reviewers, and fix the gaps in the core, the library and the features
when-to-use: When the user asks to migrate, translate, port or "convert to bdd" existing Playwright/TestTrack specs, to review whether bdd features actually test what they claim, or to compare features against the tests they replaced
context: fork
effort: high
argument-hint: "<old spec folder or files> [<features folder>]"
---

# Translating tests into features, and proving them

A feature file is only as honest as the bindings behind its phrases. Translating a spec is three
passes, and the second and third are where the value is: the first produces green tests, the
second finds which of them would stay green with the behaviour broken, the third turns those into
real checks — usually by adding a name or a signal to the core rather than a wait to the test.

The lead's rules that hold throughout: **a wait, sleep, pixel scan or retry in a test means a
missing signal or name in the core — fix the core**; nothing is reimplemented that the platform
already has; nothing is committed without the lead's order.

## Pass 1 — translate

1. Read the library guide first: `public/libraries/bdd/README.md` (vocabulary, tiers, "How a
   feature runs") and `public/libraries/bdd/CLAUDE.md` (the rules that must not regress and the
   facts that cost a run each). `npx grok-bdd list-steps` in the package prints every phrase.
2. Map each old spec to one feature. Scenarios that share a dataset and a viewer become one
   `@journey` (Background once, scenarios in order, each puts back what it changed). Features that
   share a subject go in one folder (`features/viewers/box-plot/`) — one page per folder.
3. For every old assertion, write the step that says the same thing. When the old spec probed
   pixels, clicked candidate offsets, sniffed colors or slept, ask what the platform knows that it
   did not say — a hit area (`getWidgetStatus().hitAreas`), a reading (`getWidgetStatus().values`),
   a name on a control, an event — and add that to the core with one line in the viewer's
   `CLAUDE.md` "Automation surface" section. Probe the live DOM before naming anything.
4. `npx grok-bdd compile`, run, iterate until green. Keep the old spec, its helpers and its md in
   place; the lead decides what to delete.

## Pass 2 — backward-match with independent reviewers

Green is not done. Spawn **one read-only reviewer agent per old-spec/feature pair**, in parallel,
with the prompt below. Its report is the deliverable of this pass; consolidate the reports into
systemic findings (one library or core fix serves every feature) and per-feature restorations.

The reviewer prompt (fill the paths):

```
You are a read-only reviewer. Do not edit, create or run anything.
Pair: NEW feature <path>, its generated spec <path> (shows which binding each phrase calls),
OLD spec <path> with its helpers and md.
Bindings to read in full for every phrase the generated spec imports: the tier's steps.ts, the
platform and common bindings, src/runtime/viewers.ts (the in-page runtime), harness.ts,
locate.ts, args.ts, the package's own bindings. The core file that provides a hit area or a name
when you need to know what it points at.
1. Table every assertion of the OLD spec: covered / weakened (how) / dropped (what is lost) /
   replaced (only if at least as strong), with the new feature line.
2. For every Then/And step of the NEW feature trace phrase → binding → in-page code and decide
   whether it can pass with the behaviour broken. Think in mutations: the property ignored, the
   render never happened, the baseline taken after the change, a stale hit area or tooltip, a
   table with the right shape and blank cells, a one-pixel change. Flag tautologies, checks of
   values the step itself wrote, "exists" where "shown" is claimed, and wording that promises
   more than the code verifies.
3. For every When/Given step: does it perform the action it names on the element it names, and
   would the following Then notice if it silently did nothing?
4. Flag any sleep, fixed wait, cap used as a wait, or retry loop in the bindings used.
5. List what in the OLD spec was itself vacuous or a hack, so nobody restores it.
Report: verdict in two sentences; the coverage table; findings most serious first with feature
line, binding file:line, claimed vs checked, a concrete false-pass scenario, a one-sentence fix.
```

What the box plot review found, as the checklist for the next translation:

- **"Repainted" meant one pixel.** A change detector is not a shape check. Where a shape is
  claimed, pair it with `more/less ink`, a per-area `the "X" area … should have more ink than
  before`, a `should have a "stats" area`, a color in an area; a pure chrome toggle gets `should
  have repainted by at least N pixels`.
- **A hidden element keeps its last text.** The tooltip is one element, hidden between hovers;
  `contain text` passed on the previous hover's text. Fixed in the core (`Tooltip.hide()` empties
  it) and in the binding (a tooltip is matched visible-only).
- **Negative checks read the first tick.** "Same range as before", "not repainted", "coloring
  survives" succeeded before a reset that landed a tick later. Negative checks read once the
  viewer is quiet — after the frame and after whatever it says is pending.
- **A silent cap is a wait.** A settle that gives up after 300 ms hides a late repaint, which then
  satisfies the next step. The viewer now says whether a render is pending
  (`viewer.isRenderPending`); the settle waits exactly when it has to and never otherwise.
- **Margins.** A highlight check that accepts one more orange pixel passes on the hover halo; the
  margin scales with the selection (`highlightMargin` in viewers.ts). A palette count that admits
  axis grey and label black proves nothing; compare areas.
- **Round trips of the property bag.** `sets X` then `X should be` proves the look kept the value,
  not that the viewer applied it: add paint evidence, read `viewer.dataFrame` (not the Table
  property) for a rebind, read a hit area or a reading that only exists when the work was done.
- **Shape without content.** A result table with the right name, row count and columns can be
  all blanks: `should have no missing values in "<computed column>" column`.
- **Floors owned by nobody.** Error and balloon floors ran from login, so an earlier scenario's
  error failed a later one; each journey scenario now clears them at its start.
- **A typo passes.** "All rows where RACE is Blakc are selected" matched 0 of 0; the positive
  data steps fail when no row matches.
- **Precondition of a fixture.** A scenario about empty categories must assert the fixture has
  blanks before it asserts what the viewer does with them.

## Pass 3 — fix, in this order

1. Put the decisions to the lead first, one question per systemic finding, options with the
   recommended one first. The lead chooses; do not fix on your own judgment where the reading
   changes the work.
2. Core signals and names (the viewer, the tooltip, `WidgetStatus`), with the JS API line
   (`grok_api.dart` + `grok_api.g.ts` by hand + the `js-api` getter, rebuild the bundle), the
   analyzer, and a line in the viewer's `CLAUDE.md`.
3. Library runtime and bindings (`viewers.ts`, the tier's `steps.ts`, `platform/data.ts`,
   `harness.ts`), `npm run build`, the unit tests.
4. The features, `grok-bdd compile`, the run; a failure is read as evidence first (the fixture,
   the platform) and a phrase second.
5. The docs: the library README vocabulary block and CLAUDE.md rule, the package README, the
   core CLAUDE.md, the memory file. Numbers (scenarios, seconds) from the last green run.

## Where things are

- Library: `public/libraries/bdd` — `src/runtime/viewers.ts` (in-page `__bdd`), `bindings/tiers/viewers/steps.ts`,
  `bindings/platform/data.ts`, `src/runtime/harness.ts` (`feature`, `journey`, floors).
- The first translation, as the worked example: `public/packages/UsageAnalysis/bdd/features/viewers/box-plot/`
  from `public/packages/UsageAnalysis/files/TestTrack/Viewers/BoxPlot/`.
- Core automation surfaces: `core/client/d4/lib/src/viewers/<viewer>/CLAUDE.md` "Automation surface",
  `core/client/d4/CLAUDE.md` (AppEvents, `isRenderPending`, `WidgetStatus.values`).
- Run: `cd <package>/bdd && node ../../../libraries/bdd/bin/grok-bdd.js run --reporter=list`
  (junction-linked checkout) or `npx grok-bdd run`; per-step timings with
  `PLAYWRIGHT_JSON_OUTPUT_NAME=run.json … --reporter=list,json`.
