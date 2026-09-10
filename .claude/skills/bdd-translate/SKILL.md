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

The platform facts these passes rest on — what a viewer reports, how the menus and hover behave,
what `isRenderPending` covers, which fixture has which counts — live in
`public/libraries/bdd/CLAUDE.md`. This file is the process.

## Pass 0 — is there anything to test against?

Look at the viewer's `getWidgetStatus()` before you plan anything. **A viewer that has none is the
whole job**, and the features are the easy half. Write the status first, get it into the served
bundle, and only then write a feature — a feature run against a stale bundle reports "the viewer
has no readings", which reads exactly like a wrong phrase and will cost you an afternoon.

The order that works, with agents:

1. **Survey** — one read-only agent per viewer or small group. Ask for: a per-scenario inventory of
   the old spec marked `direct` / `needs-core` / `manual` / `vacuous`; the fixture numbers verified
   against the data files; **the core signals to add**, as a table of "old hack → the reading or
   area that replaces it → where in the code it comes from"; and the step phrases, each marked
   `exists` / `exists-but-hardcoded` / `new`.
2. **Core** — one implementer per viewer, the survey as its brief. Status code goes in a part file
   next to the viewer as a private feature class; the core keeps the render flags and one
   delegating line. A new file must be added to the viewer library's `part` list or Dart will not
   compile it. Finish with the viewer's `CLAUDE.md` "Automation surface" section: every area and
   every value, what it reads from, **and when it is absent**. That section is what the feature
   writers work from, so it has to be complete.
3. **Features** — one agent per viewer, the `CLAUDE.md` as the authoritative list.

Three agents at a time. Tell every feature agent to put **all** its new steps in its own package
binding file and to touch neither the library nor the core, or three agents edit `steps.ts` at
once. Promote the generic ones yourself, afterwards, in one pass.

**Tell agents to push back.** The best returns in the last round came from implementers refusing a
brief: one would not declare an Escape shortcut because the viewer binds no keys at all; one found
the `aria-checked` change had already been made that morning and checked the git history instead of
duplicating it; one found the reading I specified would have reported 0 for every group because the
tree files rows on leaves only. A brief is a hypothesis.

## Pass 1 — translate

1. Read `public/libraries/bdd/README.md` (vocabulary, tiers, how a feature runs) and
   `public/libraries/bdd/CLAUDE.md` (the rules that must not regress, and the platform facts that
   cost a run each). `npx grok-bdd list-steps` prints every phrase.
2. Map each old spec to one feature. Scenarios that share a dataset and a viewer become one
   `@journey` — Background once, scenarios in order, each putting back what it changed. Features
   that share a subject go in one folder; one page per folder.
3. For every old assertion, write the step that says the same thing. Where the old spec probed
   pixels, clicked candidate offsets, sniffed colours or slept, ask what the platform knows that it
   did not say — a hit area, a reading, a name on a control, an event — and add that to the core.
   Probe the live DOM before naming anything.
4. `npx grok-bdd compile`, run, iterate until green **twice**.

## Pass 2 — backward-match with independent reviewers

Green is not done. Spawn **one read-only reviewer agent per old-spec/feature pair**, in parallel.
Its report is the deliverable; consolidate the reports into systemic findings (one library or core
fix serves every feature) and per-feature restorations.

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

## Pass 3 — fix, in this order

1. Put the decisions to the lead first, one question per systemic finding, options with the
   recommended one first.
2. Core signals and names, with the JS API line where one is needed, the analyzer, and the viewer's
   `CLAUDE.md`.
3. Library runtime and bindings, `npm run build`, the unit tests.
4. The features, `grok-bdd compile`, the run; a failure is read as evidence first and a phrase
   second.
5. The docs: the library README and CLAUDE.md, the package README, the memory file. Numbers from
   the last green run.

## What a green feature still gets wrong

The checklist, from five rounds of reviews. Each line is a way a passing test was found to prove
nothing.

**Claims that are not claims**
- **A property set, then the same property read back.** That proves the look kept a value, not that
  the viewer applied it. Pair it with a reading or an area that only moves when the work was done.
- **"Repainted" meant one pixel**, and an ink threshold means "a lot of pixels moved", not what
  they now mean. Equal bars swapping places keep every colour count. Where a shape is claimed, pair
  it with per-area ink, a colour in an area, or a signature hash.
- **Ink is unusable for anything that widens an axis**: turning SPC on *reduces* painted pixels,
  because the axis stretches and the series compresses by more than the overlay adds.
- **A canvas liveness probe.** "A `<canvas>` element still exists under the viewer root" was an
  entire round's idea of an assertion. The viewer's own `error` is the honest replacement.
- **Counting what the test itself wrote.** One spec wrote two formula lines as JSON, then parsed
  that same JSON and counted 2. Another set a `<select>`'s value and asserted the `<select>` had it.
- **A negative that never had a positive.** Four canvas clicks that may miss every node prove
  nothing when nothing gets selected; a claim about order or membership needs the thing to have
  been there to lose.
- **A trivially true comparison**: "the nested leaf's rows are fewer than the largest group's" is
  480 vs 896. State both numbers.
- **A guard that passes when its subject is missing** — a rect asserted ±1 px inside an `if (input)`.

**Waits, and what they were hiding**
- A silent cap is a wait: a settle that gives up after 300 ms hides a late repaint, which then
  satisfies the next step. The viewer says whether a render is pending; the settle waits exactly
  when it has to.
- Every helper file full of settle loops is a missing signal. One 106-line helper existed entirely
  because a viewer deferred its first refresh with a bare `Timer.run`, which no flag covered; its
  own comment measured the cost at 0.9–1.2 s per spec.
- Check what the viewer's pending flags actually cover before blaming a flake on load. Found this
  way: an overlay invalidation, a bar chart's refresh timer, a legend's settle timer, a 3D plot's
  label loads, a zoom animation walking the viewport over ten 50 ms frames.
- `VIEWER_RENDERED` can fire before the data lands — one viewer fired it synchronously ahead of the
  future that put the aggregation on the grid, so every settle in the tier returned early.
- Negative checks read the first tick. "Same as before", "not repainted" must read once the viewer
  is quiet.

**Evidence**
- A hidden element keeps its last text and its last geometry. A status reports an area only while
  the thing is drawn; "how many does it hold" is read off the look, "how many did it draw" off the
  frame.
- A picture that cannot be read (WebGL, an inner viewer, a cell) gets a signature: sampled and
  hashed in one task. A signature that "differs" must have exactly one cause — fonts fetched per
  scene, labels added asynchronously and an auto-rotating camera each faked a repaint until the
  core stopped doing them on a test page.
- Shape without content: a result table with the right name, row count and columns can be all
  blanks.
- A typo passes: "rows where RACE is Blakc are selected" matched 0 of 0.
- Assert the fixture in the Background, and **count it against the data file**. The old specs were
  written against a bigger table; a number carried over unchecked is the first thing a
  green-looking draft gets wrong.
- An old `md` goes stale. One documented three open bugs that had all been fixed in the code —
  verify before writing a `@known-failure`.
- A spec marked `test.fail(true, …)` proves nothing at all. Check for it before trusting a coverage
  claim.

**Process**
- **Green twice on a quiet machine is not green.** Run on four workers more than once, and headed
  once — Chrome moves a 2D canvas from GPU to CPU after enough `getImageData` readbacks and the
  first paint after the move differs in antialiasing, which headless never shows.
- A reviewer's finding is a hypothesis until the suite runs it. Three "fixes" in one round asserted
  things the product does not do.
- `@known-failure` is how an open bug is translated: the scenario is written honestly, its failure
  is expected, and its **passing** fails the test ("the bug is fixed, remove the tag"). Never soften
  an assertion to keep a suite green.
- Every journey scenario owns its errors and balloons.
- Time the suite with the JSON reporter before calling it done — every Gherkin step is a Playwright
  step with a duration.
- **Write patch scripts with the Write tool, never a heredoc.** The Bash tool collapses `\\` inside
  heredocs, so a `\\s` in a regex becomes `\s`, the replacement silently fails to match, and you
  lose twenty minutes to "why didn't this apply".

## Recurring bindings — promote as soon as a second viewer wants one

A step written for one viewer is usually not about that viewer. The test: **could another viewer
want this phrase?** If yes, it belongs in `bindings/tiers/viewers/widgets.ts` (the tier's second
half, for exactly this) and takes `{widget}`, not a named viewer.

The signs, all of which produced duplicates before anyone noticed:

- **The same body under two names.** `user clicks on empty plot space of box plot viewer` and
  `… of bar chart viewer` were the same function twice. So were a scatter plot's and a PC plot's
  range-slider drags, and a scatter plot's and a pie chart's column selectors.
- **A phrase that already takes `{widget}` but lives in a package.** Its home is wrong, nothing
  else is.
- **A gesture, not a fact.** A lasso, a drag between two widgets, a right-click on the viewer's own
  menu, a click on a plain dialog checkbox — none of these are any viewer's business.
- **A property every viewer has.** Description Position, the title, the row source.
- **A convention several viewers publish.** Once two viewers report `fields` and `<COL> of row <r>`,
  the card steps are shared vocabulary; once two report `columns`/`rows`, so is
  `the cells of {widget} should be N wide and M tall`.

What genuinely stays in a package: a gesture only that viewer's DOM has (a sketch designer, an
axis-label reorder), or an arithmetic only it can check (an aggregation against an independent
`groupBy`, a correlation against `DG.Stats`).

Two things to watch when promoting. A generic implementation must keep what the specific one did —
one promotion silently dropped `hitArea(…, beforeChange: true)`, which both polls for the area and
takes the canvas baseline, and the viewer lost its "lower than before" baseline and the frame its
slider needed to lay itself out. And check the phrase the library already resolves: a status that
reports `x scroll min handle` when the library looks for `range min handle "x"` takes a DOM
fallback instead.

## Package viewers deploy differently

A JS/TS viewer needs the same three signals as a Dart one (`getWidgetStatus`, `isRenderPending`,
`onRendered`), but its surface reaches a running stand only when the **package** is rebuilt and
published — and a feature run against a stale package fails exactly like a mistyped phrase. Worse
when the viewer lives in a library and is registered by a package: the package's
`node_modules/@datagrok-libraries/<lib>` must point at the checkout, or it compiles the registry
copy and the new keys silently do not appear. Establish this before writing features, not after.

## Where things are

- Library: `public/libraries/bdd` — `src/runtime/viewers.ts` (in-page `__bdd`),
  `bindings/tiers/viewers/steps.ts` and `widgets.ts`, `bindings/platform/`, `src/runtime/harness.ts`.
- Worked examples: `public/packages/UsageAnalysis/bdd/features/viewers/` — the box plot was the
  first, the density plot the most recent.
- Core automation surfaces: `core/client/d4/lib/src/viewers/<viewer>/CLAUDE.md` "Automation
  surface"; `core/client/d4/CLAUDE.md` for `isRenderPending` and `WidgetStatus`.
- Run: `cd <package>/bdd && node ../../../libraries/bdd/bin/grok-bdd.js run --reporter=list`, or
  Playwright directly with `BDD_ROOT=$(pwd -W)` and the library's `dist/playwright.config.js`.
