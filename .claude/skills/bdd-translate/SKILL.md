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
passes: the first produces green tests, the second finds which of them would stay green with the
behaviour broken, the third turns those into real checks — usually by adding a name or a signal to
the core rather than a wait to the test.

Rules that hold throughout: **a wait, sleep, pixel scan or retry in a test means a missing signal
or name in the core — fix the core**; nothing is reimplemented that the platform already has;
nothing is committed without the lead's order. The vocabulary and the invariants are in
`public/libraries/bdd/README.md` and `CLAUDE.md` (`npx grok-bdd list-steps` prints every phrase);
a viewer's areas and readings are in `core/client/d4/lib/src/viewers/<viewer>/CLAUDE.md`.

## Pass 0 — is there anything to test against?

Look at the viewer's `getWidgetStatus()` first. **A viewer that has none is the whole job**; the
features are the easy half. Write the status, get it into the served bundle, and only then write a
feature — a feature run against a stale bundle reports "the viewer has no readings", which reads
exactly like a wrong phrase.

With agents, three at a time:

1. **Survey** — one read-only agent per viewer: a per-scenario inventory of the old spec marked
   `direct` / `needs-core` / `manual` / `vacuous`; the fixture numbers verified against the data
   files; the core signals to add as "old hack → the reading or area that replaces it → where in
   the code it comes from"; the step phrases marked `exists` / `new`.
2. **Core** — one implementer per viewer. Status code goes in a part file next to the viewer as a
   private feature class (added to the library's `part` list); finish with the viewer's
   `CLAUDE.md` "Automation surface": every area and value, what it reads from, and when it is
   absent.
3. **Features** — one agent per viewer, that section as the authoritative list. Every new step
   goes in the agent's own package binding file; promote the generic ones yourself afterwards.

Tell agents to push back: a brief is a hypothesis, and the best returns came from an implementer
refusing one (a shortcut the viewer does not bind, a reading that would report 0 for every group).

## Pass 1 — translate

1. Map each old spec to one feature. Scenarios that share a dataset and a viewer become one
   `@journey` — Background once, scenarios in order, each putting back what it changed.
2. For every old assertion, write the step that says the same thing. Where the old spec probed
   pixels, clicked candidate offsets, sniffed colours or slept, ask what the platform knows that it
   did not say — a hit area, a reading, a name on a control, an event — and add that to the core.
   Probe the live DOM before naming anything.
3. `npx grok-bdd compile`, run, iterate until green **twice**, then on four workers, then headed
   once (a GPU-rasterized canvas repaints differently after enough pixel readbacks).

## Pass 2 — backward-match with independent reviewers

Green is not done. One read-only reviewer agent per old-spec/feature pair, in parallel:

```
You are a read-only reviewer. Do not edit, create or run anything.
Pair: NEW feature <path>, its generated spec <path> (which binding each phrase calls), OLD spec
<path> with its helpers and md. Read in full every binding the spec imports, the library's
runtime (src/runtime/viewer-runtime.ts, viewers.ts, viewer-pixels.ts, harness.ts, locate.ts) and
the core file behind any hit area or reading you need to understand.
1. Table every assertion of the OLD spec: covered / weakened (how) / dropped (what is lost) /
   replaced (only if at least as strong), with the new feature line.
2. For every Then of the NEW feature trace phrase → binding → in-page code and decide whether it
   can pass with the behaviour broken. Think in mutations: the property ignored, the render never
   happened, the baseline taken after the change, a stale hit area or tooltip, a table with the
   right shape and blank cells, a one-pixel change. Flag tautologies, checks of values the step
   itself wrote, "exists" where "shown" is claimed, wording that promises more than is verified.
3. For every When/Given: does it do what it names on the element it names, and would the Then
   after it notice if it silently did nothing?
4. Flag any sleep, fixed wait, cap used as a wait, or retry loop in the bindings used.
5. List what in the OLD spec was itself vacuous or a hack, so nobody restores it.
Report: verdict in two sentences; the coverage table; findings most serious first with feature
line, binding file:line, claimed vs checked, a concrete false-pass scenario, a one-sentence fix.
```

Consolidate the reports into systemic findings (one library or core fix serves every feature)
and per-feature restorations. A reviewer's finding is a hypothesis until the suite runs it.

## Pass 3 — fix, in this order

1. The decisions to the lead first, one question per systemic finding, recommended option first.
2. Core signals and names (with the JS API line where one is needed, the analyzer, the viewer's
   `CLAUDE.md`).
3. Library runtime and bindings, `npm run build`, `npm run test:unit`.
4. The features, `grok-bdd compile`, the run; a failure is evidence first and a phrase second.
5. The docs: the library README and CLAUDE.md, the package README. Numbers from the last green run,
   timed with the JSON reporter.

## What a green feature still gets wrong

Each line is a way a passing test was found to prove nothing.

- **A property set, then the same property read back** proves the look kept a value, not that the
  viewer applied it: pair it with a reading or an area that only moves when the work was done.
- **"Repainted" is one pixel** and ink means "pixels moved", not what they now mean: equal bars
  swapping places keep every colour count. Ink is unusable for anything that widens an axis.
  Where a shape is claimed, pair it with per-area ink, a colour in an area, or a signature.
- **A canvas liveness probe** ("a `<canvas>` still exists") is not an assertion; the viewer's own
  `error` reading is.
- **Counting what the test itself wrote**; **a negative that never had a positive** (four canvas
  clicks that may miss every node); **a trivially true comparison** (state both numbers); **a
  guard that passes when its subject is missing** (a rect asserted inside an `if (input)`).
- **A silent cap is a wait**: a settle that gives up after 300 ms hides a late repaint, which then
  satisfies the next step. Check what the viewer's pending flags actually cover before blaming a
  flake on load (found this way: an overlay invalidation, a refresh timer, a legend's settle timer,
  a 3D plot's label loads, a zoom animation).
- **Negative checks read the first tick**: "same as before" and "not repainted" must read once the
  viewer is quiet.
- **A hidden element keeps its last text and its last geometry**: a status reports an area only
  while the thing is drawn; the tooltip is what is visible now.
- **A picture that cannot be read** (WebGL, an inner viewer) gets a signature hashed in one task,
  and a signature that "differs" must have exactly one cause.
- **Shape without content**: a result table with the right columns can be all blanks. **A typo
  passes**: "rows where RACE is Blakc" matched 0 of 0 — the data steps fail on an empty match.
- **Count the fixture against the data file**; an old `md` goes stale (three "open bugs" were all
  fixed) — verify before writing a `@known-failure`; a `test.fail(true, …)` spec proves nothing.
- **Write patch scripts with the Write tool, never a heredoc** — the Bash tool collapses `\\`.

## Where a step belongs

A step written for one viewer is usually not about that viewer. It belongs in
`bindings/tiers/viewers/widgets.ts` and takes `{widget}` when another viewer could want it: the
same body under two names, a gesture rather than a fact (a lasso, a drag between widgets, the
viewer's own menu), a property every viewer has, a convention several viewers publish (`fields`
and `<COL> of row <r>`, `columns`/`rows`). It stays in the package when it reads a status key
shape only that viewer publishes (`bar <category>`, `share of "<slice>"`, `record of card <n>`),
performs a gesture only that viewer's DOM has (a sketch designer, an axis-label reorder), or
checks an arithmetic only it can (an aggregation against `groupBy`, a correlation against
`DG.Stats`). When promoting, keep what the specific one did (`hitArea(…, beforeChange: true)`
both polls for the area and takes the baseline) and check the phrase the library already
resolves.

A JS/TS viewer needs the same three signals as a Dart one (`getWidgetStatus`, `isRenderPending`,
`onRendered`), and its surface reaches a running stand only when the **package** is rebuilt and
published — for a viewer that lives in a library, the package's `node_modules/@datagrok-libraries/<lib>`
must point at the checkout. Establish this before writing features.

Worked examples: `public/packages/UsageAnalysis/bdd/features/viewers/` (the box plot was the
first; `viewer-chrome.feature` is the outline over what every viewer shares). Run with
`npx grok-bdd run --reporter=list generated/<folder>` from the package.
