# @datagrok-libraries/bdd — behavioral automation (Gherkin → Playwright)

Feature files bound to the u2/platform vocabulary, compiled deterministically into Playwright specs
that packages own. Design record: `core/docs/features/ui2/automation/BRAINSTORM.md` (rulings by the
lead: standard Gherkin; committed, drift-gated codegen; u2-centred; overridable composition;
generic kinds; reserved platform names; packages own their tests; one page per worker; a global
CLI; NOT wired into `grok test`; §8 there is the current state, gaps and next steps). The contract
it consumes and the u2 findings: `core/docs/features/ui2/AUTOMATION.md`; the lane map:
`core/docs/features/ui2/TESTING.md`.
Authoring guide: `README.md` here.

## Layout

```
bin/grok-bdd.js    launcher: prefers the package-local install (its specs import that copy), else this one; runs dist/src/cli.js
src/project.ts     a project = a dir with features/ (a package's bdd/, or the library root); config bdd.config.json {tiers}; binding sources
src/discover.ts    imports binding modules (library: this build's bindings/; project: .ts through tsx) and maps step fn → export name
src/gherkin.ts     @cucumber/gherkin → flat model (And/But resolved, outlines expanded, rule backgrounds folded)
src/registry.ts    Given/When/Then/Step, element(), context(), alias(), kind(), dataset(), defineParameterType(); one global registry
src/nouns.ts       phrase (+ context) → NounRef (pure; shared by compiler and runtime); parts of elements AND kinds
src/match.ts       cucumber-expressions matching + specificity
src/compile.ts     FeatureModel → *.test.ts (feature(test) session, test() per scenario — or one test with soft-step scenarios for a @journey feature, test.step per Gherkin step, enter/leave)
src/init.ts        `grok-bdd init`: scaffold bdd/ in a package (never overwrites; merges .vscode settings, .gitignore, package.json)
src/cli.ts         init | compile [--check] | lint | list-steps | run [playwright args]  (run = check + Playwright with playwright.config)
src/runtime/       args (el/ds/enter/leave), locate, gestures, assertions, harness (feature session + resetShell + error floor), viewers (in-page `window.__bdd`: immediate rendering, render settles, properties by caption, hit areas, canvas ink, context menus), global-setup, index
src/index.ts       what package bindings import from '@datagrok-libraries/bdd'
bindings/common/   parameter-types, kinds (ALL u2 data-u2 kinds + Dart conventions), steps, session — base, always loaded
bindings/platform/ the shell: elements (reserved names), datasets, steps (open a dataset, switch views, projects), data (the current table's selection, filter, rows, columns, colors, workspace tables) — base, always loaded
bindings/tiers/<t>/ opt-in tiers (`viewers`: add/configure viewers, properties by caption, context menus and hit areas, canvas ink, events, the error floor); a project names them in bdd.config.json
features/ generated/ bdd.config.json   the library's own project (platform smoke feature, tier viewers)
playwright.config.ts  the one config every project runs with (BDD_ROOT → testDir/outputDir/storageState under the project)
dist/              `npm run build` output (tsc, ESM, .d.ts, source maps) — what `exports` and the launcher use; gitignored
tests/             node:test via tsx (nouns, compile, project, init)
```

A package project: `<pkg>/bdd/{package.json {"type":"module"}, bdd.config.json, features/, bindings/,
generated/}`; sample `packages/U2Demo/bdd`; the first production project is `packages/UsageAnalysis/bdd`
(`features/viewers/box-plot/*.feature`: six journeys, 50 scenarios, the whole of the six TestTrack
box plot specs under `files/TestTrack/Viewers/BoxPlot/` and their helpers, 36 s for all six on one page,
reviewed against the specs they replaced on 2026-09-07 (late) — the `/bdd-translate` skill is that process).
The second is `packages/Bio/bdd` (2026-09-08): eighteen journeys, 101 scenarios, 1.9 min on one worker,
from the TestTrack Bio specs, the package's `playwright/` specs and their md files — the first package
translation, and the one that added the top menu, the command signal, the column and function steps
and the JS-viewer conventions below.

## We test our own platform, not a black box

The one principle every step here follows, and the reason the box plot feature runs in 14 s where
the hand-written spec it replaced took 49 s for one test: **when a test waits, something is
missing, and the fix goes into the core, not into the test.**

- A viewer that repaints on a debounce or an animation frame is waiting for a human. Tests do not
  need that: the runtime sets `viewer.immediateRendering` on every viewer the page holds or adds
  (`onViewerAdded`), so a property set repaints on the next task and `onViewerRendered` is the
  signal. No `waitForTimeout`, no polling the canvas on a timer.
- A menu that builds on a timer says so: `grok.events.onContextMenuShown` / `onContextMenuClosed`
  (added to the core for this). Arm the event, do the gesture, await the event.
- A region a test wants to right-click is a fact the viewer knows: `getWidgetStatus().hitAreas`
  (box plot: view, x axis, y axis, stats, p value, group comparison, color scale, marker). A test
  that scans pixel offsets until a menu looks right is guessing at what the viewer could tell it.
- A state a test wants to read must be in the DOM as a state: grayed menu items and gated
  property rows carry `aria-disabled` (added), not only opacity 0.5.
- Data is a lever too: every paint of a marker viewer costs one draw per row. `demog-1000` is the
  stratified 1000-row demog for viewer features; a feature that does not assert on row counts has
  no business painting 5850 markers seventy times.
- Selector names are the platform's contract: `camelCaseToCss` (prop_gen) stopped emitting
  `marker--size` for "Marker Size" — the name is `div-column-combobox-marker-size` now; the old
  double-dash names in the TestTrack branch specs are stale.

When a feature needs a wait, a sleep, a pixel scan or a retry loop, stop and find the missing
signal or name in the core (`d4` viewers, `xamgle` shell, the js-api), add it there with the
owner's blessing, and write the step against it. Duct tape in the test is a bug report that
nobody filed.

## Rules that must not regress

- **Everything a spec imports from the library goes through `dist/` by package subpath**
  (`@datagrok-libraries/bdd/runtime`, `…/bindings/common/steps`); project bindings are imported
  relative (`../bindings/steps.js` → Playwright maps to `.ts`). One registry instance per process is
  what makes context/element registration visible to the runtime — never mix `src/` and `dist/`
  imports in one run. The launcher resolves the package-local copy first for the same reason.
- **ESM everywhere.** The library is `"type": "module"`; a package's `bdd/package.json` marks its
  specs and bindings as ESM too — a CJS spec would `require()` the ESM dist and fail on Node 20.
  Project `.ts` bindings load at compile time through `tsx/esm/api` `register()`; at run time
  Playwright transpiles them.
- **One `@playwright/test` per run.** Specs import it from the package, the library's runtime from
  its own real path; a second copy fails with "Requiring @playwright/test second time". Until the
  library is on npm, packages depend on it by path (`"file:../../libraries/bdd"`, written by
  `init` from a checkout; npm symlinks the directory, `npm ci` needs no registry, `.bin/grok-bdd`
  appears — a `^version` on the unpublished package broke `npm ci` with a 404) and `grok-bdd link`
  (cli.ts) links the package's `node_modules/@playwright/test` to the library's copy (the package's
  own goes to `node_modules/.bdd-link-backup`, `--undo` or an `npm ci` restores it; the library
  link too when npm did not make it). On a registry install the peer dependency is shared. The
  harness receives `test` from the spec and never imports it.
- **One page per worker** (`src/runtime/harness.ts`). It was one page per feature *folder* until
  Playwright's file-by-file hand-off made that meaningless under four workers — consecutive files
  of a worker are never from the same folder, so every feature booted its own page (20 s of login
  each). `feature(test)` registers `afterEach` (leave the context, `resetShell`) and `afterAll`
  (the feature's `atFeatureEnd` cleanups); the page is created inside the first test
  (`session.page(browser)`) so Playwright merges the project's context options, and kept in the
  module-level `shared`, so `user is logged in` only resets the shell (0.1 s instead of ~4 s).
  **Feature isolation is the reset shell plus the feature's own cleanups, not the tab.** Every test
  still gets its own trace (Playwright starts a chunk per existing context at each test start); a
  failed test restarts the worker, which drops the page. NEVER leave several Datagrok pages open in
  one browser — six live clients made every step 2–3× slower. The generated spec calls `test()`
  itself so reports point at the spec line, and every step is `session.step(line, title, fn)`, a
  Playwright step whose `location` is the feature line.
- **A throw inside an `expect.poll` callback ends the poll** — Playwright calls it outside its own
  try (`invokePollMatcher`: `const value = await actual()`), so the first read that finds no area
  fails the step instead of retrying it. A callback that reads something the viewer may be between
  layouts of returns `false` and keeps the reason in a variable the `catch` after the poll throws
  (`expectAreaGrew`, `areaBiggerThanArea`). A single non-polling `hitAreas()` read behind a claim is
  the same bug without the poll: the word cloud reported no areas at all for one frame and
  `taller than` had already failed. An in-page read counts too — `should be painted` polled
  `__bdd.ink`, which throws "has no canvas" for a viewer between two layouts, and gave up on the
  first one; it catches, keeps the reason and tells it if the poll really does run out.
- **Every check in the library imports `expect` from `src/runtime/patience.js`**, never from
  `@playwright/test`. It is the same `expect` behind a proxy, and while a `@known-failure` scenario
  runs it is `configure`d down to `KNOWN_FAILURE_MS` (3 s) — such a scenario states a defect the
  product has, and spending the suite's 15 s expect budget on each of its assertions was 100 s of a
  run proving what the tag already says (the heat map's journey 21.8 s → 10.1 s, pc plot's
  transformations 21 s → 6.8 s). A check that names its own `{timeout}` wraps it in `pollMs(…)`,
  which the same switch narrows; a check that must wait longer than that inside such a scenario is
  the one that says so.
- **A gesture aims at where the viewer has finished putting the thing.** `hitArea(…, beforeChange)`
  settles the viewer before it reads the area and takes the baseline: a network diagram still
  running its physics moves the node between the read and the click, and the click then selects
  nothing (`network-diagram-selection`, `node "F"` → 0 rows). A viewer that reports nothing pending
  answers at once, so the settle costs a roundtrip; one that never stops pending is the following
  claim's problem, not the gesture's. **And a finished render is not the end of the moving**: a
  title just set shows an auto-sizing text area, whose resize observer moves the pivot's inner
  grid down a row on a later frame, and the right-click meant for `grid header None med(AGE)`
  opened the *cell* menu one row below (one in twelve runs of `pivot-table-persistence`, four in
  the full suite). So before a right-click the anchor's box — what every area is relative to —
  must agree on two consecutive frames (`stableArea`, in `menuPoint`): one frame of cost, and a
  cheap read. It is not on every gesture: a first cut put it in `hitArea(…, beforeChange)` too,
  reading the whole widget status on each of three frames, and that alone cost the suite about
  80 s of test time (762 → 845 s) for a move seen only under a right-click. And it is a defence,
  not the cure: the title
  grows through a **100 ms polling** resizer (`handleResize`, `d4/src/utils/utils.dart`), later
  than any frame count, and the product now sizes the title in the same pass as the property write
  (`refreshTitle`, `viewer_base.dart`). A late mover the library cannot see is a product fix,
  every time.
- **Typed text is verified before it is committed** (`gestures.typeVerified`). Control+A, the text,
  then the editor is read back and retyped until it holds exactly the text: a keystroke that
  creates or rebuilds the editor lands at an unpredictable moment (the column picker's `"EXS"`),
  and an editor the widget closes under the typing keeps only what came after — the grid's cell
  editor was typed `"51"` and committed `"1"`. `typeInto` also leaves an editor that already has
  the focus unclicked: a click is what would blur a cell editor, and a blurred cell editor commits
  and closes. An editor whose value cannot be read back (a contenteditable) is typed into once, and
  so is one that refuses the text (read-only, disabled) — refusing it is what the step after such
  a typing claims. The verification is what found the product cause of all three losses: a table
  view focuses its grid a second after it opens (`table_view.dart`), whatever the user is typing
  by then; a temporary hook on `Element.remove` in the runtime caught the editor's `onChange`
  (a blur) destroying it, with the stack. That timer now leaves an editable element alone.
- **An area typed into must own the focus** (`typeIntoArea`). The histogram's range input took
  the click and not the focus once in twenty runs (the same timer), and `"18"` then opened a
  cell editor on the grid, unseen, while the filter stayed where it was. The click is repeated at
  the area's current place until an editor inside the viewer is focused, the editor is pinned by
  a `data-bdd-editor` mark of its own — a `:focus` locator waits on nothing once the focus has
  moved — the text goes in verified, and the failure names what has the focus instead.
- **Every column picker goes through `gestures.pickInColumnGrid`** — one place, because the Dart
  `ColumnComboBox` cost this suite a full day of intermittent failures across the histogram, the
  scatter plot, the filter panel and the pivot table, and each of its three facts had to be learned
  separately. (1) Its search box does not exist until a letter is typed **at the selector**, and the
  selector does not always keep the focus its own mouse-down gave it, so the letter is pressed on
  that element (`selector.press(name[0])`), not on the page. (2) The letter that opens the box lands
  in it at an unpredictable moment — before anything typed afterwards, or after all of it: `"SEX"`
  came back `"EXS"` and `"RACE"` as `"RACER"` — so the name goes in over whatever is there until the
  box holds it and nothing else. A name that arrives mangled matches no column, `currentColumnName`
  resolves to −1, the row does not move, nothing is announced and **the picker stays open having
  taken nothing**, showing the row the pointer previewed — which is how
  `user picks "HEIGHT" in the "y" column selector` came back `WEIGHT`. (3) Enter is pressed **on the
  box**, because the grid moves the focus while it filters, and the popup's disappearance is the
  step's own outcome check: swallowing it leaves the reading it was for to time out three steps
  later. Where the pointer may go afterwards is the caller's business — an on-viewer selector needs
  it off the popup (a row it rests on is previewed onto the selector), the filter panel's needs it
  on the panel (that picker lives in a header shown only while the panel is hovered).
- **A baseline is taken on a viewer that has finished rendering.** `user sets properties of …` and
  `user resizes …` used to snapshot at whatever instant the write started, so a viewer caught
  between two layouts — a word cloud reports no areas at all there — left the "than before" claim
  after it with nothing to compare with (`before undefined`). Both quiet the viewer first
  (`writeProperties`, `resize` in `src/runtime/viewers.ts`); it costs a frame and removes a whole
  class of "the area was not there" failures.
- **The size a step asks for is held** (`resize`, `src/runtime/viewers.ts`). The dock manager sizes
  the element it hosts, so an inline size written while it is still laying a freshly docked viewer
  out is gone on the next pass — and the feature then reads a viewer at whatever width the dock
  gave it (the statistics journey's `med` and `stdev` columns were simply off screen, three
  scenarios failing on readings that were never going to arrive). The step waits for the box to be
  the same for two frames before writing, and a `MutationObserver` puts the size back if a later
  layout pass takes it away, until `user restores the size of …` or the feature ends.
- **An in-page step returns nothing it does not read** (2026-09-08): `page.evaluate` serializes
  its return value to Node, and a platform object can be enormous — Bio's init step returned the
  SeqHelper, which holds the RDKit module and its 16 MB WASM heap, so every call cost 10 s of
  base64 (`typedArrayToBase64` in Playwright's serializer, seen in a CDP profile) — 312 s of an
  8.6 min serial run, 615 s under four workers. Await inside the evaluate and return `undefined`;
  keep results in the page (`__bddLastResult`) and read them with a small expression.
- **A failure is the feature line, the step, one sentence, and what was there instead**
  (`src/runtime/failure.ts`, `harness.ts` `session.step`, `locate.ts` `explain`; 2026-09-07, the
  lead: "the error PW gave me is very uninformative … this should be true everywhere"). `reasonOf`
  strips Playwright's API prefix, an evaluate's in-page stack, ANSI, the ms timeout phrasing and
  the `waiting for locator(…)` selector line; the rest of a call log and a matcher's own text
  stay. On a wait failure (`isWaitFailure`) the report adds `explain(page)`: the last phrase
  `refOf` parsed — its scope `not open`, or the visible elements of its kind by their labels
  (`data-u2-name`/aria-label/title/first text line), so the author sees the name to write. A
  `StepFailure` has no frames except for a programming error (`TypeError` & co., own frames only);
  its one frame is the feature line, so Playwright prints the Gherkin snippet. Errors thrown by our
  own code must already be sentences that name the alternatives (`has no "x" area; it has: …`,
  `has no "x" property; nearest: …`, `no "x" in the menu; it shows: …`) — a new throw follows that
  shape. `journeyFailure` lists the failed scenarios in the same shape, no `toEqual` diff.
  `report()` in cli.ts prints notes only for `lint` or `--verbose`: 134 of them buried the failure.
- **`@journey` = one test, the Background once, scenarios as soft steps** (`journey(test, n)` in
  harness.ts; codegen `emitJourney`): a failing scenario is recorded, the next runs, `finish()`
  fails the test with the list — the hand-written specs' `softStep`. The budget is the per-test
  timeout × the scenario count. Scenarios of a journey restore what they change (the feature's
  contract, not the harness's). Introduced 2026-09-07 because per-scenario Background + `resetShell`
  cost ~1.9 s × 13 on the box plot — two thirds of its 38 s.
- **Roundtrips are the cost, not the page**: an idle evaluate is ~1 ms, but every Playwright action
  carries the trace's per-action work, so a step is as slow as its number of actions. Hence:
  `locator.evaluate` waits for its element (no `waitFor` before it); the canvas baseline is taken
  inside the same in-page call as the change (`writeProperties`, `resize`, `menuPoint`,
  `findArea(…, beforeChange)`); `installViewerRuntime` remembers the page (forgotten on main-frame
  navigation); `locateActionable` is `filter({visible: true})` with no count; `expectVisible` and
  `expectEnabled` are one query each. A step should be locate + one action.
- **Traces keep no DOM snapshots and no per-action screenshots** (`playwright.config.ts`:
  `{mode: 'retain-on-failure', snapshots: false, screenshots: false}`): serializing the shell's
  DOM around every action was ~43% of the box plot's scenario time (11.4 s → 6.5 s measured with
  `--trace off`), and the screenshot per action another ~12 s over six features (2026-09-07, the
  lead: "we only need error screenshots"). Actions, console and network stay in the trace, the
  failure screenshot is the base config's `screenshot: 'only-on-failure'`; `grok-bdd run --trace
  on` records everything, `--video on` a video.
- **Hit areas are awaited like elements** (`viewers.ts` `hitArea`): a viewer that renders twice on
  a change (the box plot after a category switch) can report no `marker` between the two paints;
  the lookup polls up to 5 s and fails naming the areas it does report.
- **Names are global and unique; platform names are reserved.** `element()` throws on a name or
  alias already registered; a context's `element()` throws on a global name too. App vocabulary goes
  on a `context(name, def)`.
- **Context switching is compile-time tracked and runtime explicit.** Step meta `enters: '<ctx>'`
  makes the compiler validate later phrases against that context and emit `enter(page, '<ctx>')`;
  the runtime keeps the current context per page; `resetShell` leaves it. Every runtime parse goes
  through `refOf(page, target)` (locate.ts) — a parse without the page loses the context.
- **Noun resolution order** (`src/nouns.ts`): whole phrase registered (context names first, then
  global) → ordinal → split at the first scope word outside quotes (inner within outer, recursive;
  `of` names a part of a registered element or of a kind) → registered element, else generic kind by
  suffix with EVERY matching suffix kept, longest first.
- **Context-first lookup**: a context-local name or a generic kind without an explicit scope is
  searched inside the context root first, then on the whole page (portaled dialogs, notifications).
- **Kind qualifier strategies** run in the kind's `match` order, first hit wins, union (`.or`) when
  nothing matches (negative assertions still get a locator). Kinds cover every u2 `data-u2` value
  (`grep -rhoE "dataset\.u2 = '[^']+'" libraries/u2/src`) AND Dart `name=` conventions. Plain u2
  `button()`, toolbar buttons and tab headers carry NO `data-u2` — those kinds match by tag/role +
  text. Qualifiers are lowercased; text matching is case-insensitive (`exactText()` regex).
- **Owner edge**: nothing inside the outer → retry inside `[data-u2-owner="<outer's name>"]`.
- **Playwright scopes inner selectors to the element**: a `labelSelector`, a part, a `has:` filter
  is evaluated from the outer element, so `.u2demo-status > span:first-child` never matches the
  first span OF a `.u2demo-status` — write `span:first-child`. Cost a full failed run (2026-09-03).
- **`page.evaluate` must return nothing DOM/Dart-shaped**: returning `grok.shell.addView(...)`
  fails with "Cannot serialize result: object reference chain is too long" — wrap in a block and
  return undefined.
- **States** (`src/runtime/assertions.ts`): `selected` = `loc.and('[aria-selected=true], [aria-pressed=true],
  [aria-checked=true], [aria-current]')`; `expanded`/`collapsed` read the element's own
  `aria-expanded` or the first one inside (section/accordion/category headers); `visible`/`hidden`
  over several matches = any/none (`filter({visible: true})`) — a bare `toBeHidden` on a multi-match
  locator is a strict-mode error. `expectCount` prefers `tbody tr` (no header row).
- **Gestures** (`src/runtime/gestures.ts`): `editorOf` = the element itself when editable, else
  the first of `EDITOR` (inputs, selects, textareas, contenteditables), else the
  `[data-u2-part="editor"]` trigger (icon/function/columns pickers), else the element (a slider
  handle, a list). `select` = native `<select>` first; else click the editor, ArrowDown when it is
  a `role=combobox` (u2 comboboxes/typeaheads open on a keystroke, not a click), then the option
  by whole text → primary-text part (`OPTION_LABEL`) → `title`/`aria-label` (icon cells) → substring.
  `setExpanded` clicks the twistie of a tree row (a row click only selects), else the
  `aria-expanded` control. `fillIn` recognises switches through `role=switch` inside the editor.
- **Step specificity**: more literal characters first, then fewer parameters; ties are ambiguous
  errors. The literal text decides so that `user right-clicks on the {string} area of {element}`
  beats `user right-clicks (on ){element}` (which would otherwise swallow the whole phrase into
  `{element}`); `{key}` has no spaces so `user presses {key} in {element}` cannot lose to
  `user presses {key}`. Changed 2026-09-06; U2Demo's 25 specs compiled identically.
- **A phrase says what it needs** (2026-09-07, the lead's ruling): a step that reads or changes a
  viewer takes `{widget}` (`parameter-types.ts`: an element phrase ending in "viewer" or "widget",
  compiled like `{element}`), and the property steps spell it out — `user sets {string} property of
  {widget} to {string}`, `{string} property of {widget} should be {string}`. `{string} of {element}
  should be {string}` is too general to own: it would swallow every later "X of Y should be Z"
  phrase. Context-menu and pointer steps take any element because they work on any element.
- **Gestures act on the visible match** (`locateActionable` = `filter({visible: true})`, ordinals
  exempt). The Dart context menu mirrors every property under a zero-size "Properties..." submenu,
  a closed view leaves its viewers behind — both duplicate the labels a phrase names. Several
  visible matches stay a strict-mode error. `enabled`/`disabled` do NOT filter: `expectEnabled`
  evaluates all matches in-page and takes the visible ones when there are any, else all — a
  property row in the context panel has a zero-size box under `simpleMode` and still says
  `aria-disabled`.
- **Labels are found first, items second** (`locate.ts` `byLabel`): a kind whose `labelSelector`
  is `:scope > …` resolves as label → `xpath=parent::*` ∩ kind selector. The `has:` form over a
  315-item Dart popup costs ~35 ms per query, the parent form ~2 ms (probe 2026-09-07).
- **Menu items match their own label** (`:scope > .u2-menu-label, :scope > .d4-menu-item-label`,
  label before text): a group item contains its children's labels (u2 nests the submenu inside
  the item, Dart its children), so a `has:` over descendants makes "As CSV" match the Export group
  too, and "Markers" the group, the "Properties..." group and its nested mirror. Text-first was
  no better: the exact text of "Markers" is only the hidden mirror's.
- **Hover is two pointer events, and never sleeps**: leave the element to its left on the same
  line, land on its centre in one move, then check it is still there (a view still docking moves
  out from under the pointer). Leaving upwards would cross the neighbouring row and close the
  submenu. Every pointer event costs a frame (~16 ms), and hover-driven layout is synchronous, so
  there is nothing to wait for — except two browser behaviours that are not sleeps:
  **pointer moves are delivered frame-aligned**, so a DOM read issued right after `mouse.move` can
  run before the move's handler (a step reading geometry the hover reveals waits for it to be
  visible first — this was a 2-in-45 flake of the value-axis zoom); and **moves queued while the
  main thread is busy are coalesced**, so the leave-then-enter pair can collapse into one move that
  enters nothing when the pointer already rested inside. `hover` waits in-page for the element's own
  `mouseenter` and repeats the pair when it did not come.
- **A step ends when the platform is done, not when the DOM shows** (`platform/steps.ts`
  `openDataset`): opening a table starts semantic-type detection in the background (package
  detectors — Chem's SMILES detector took 300–500 ms on spgi-100), which used to land on
  whichever step came next. The step waits for the platform's own `ddt-semantic-type-detected`
  event for that data frame (identity by the Dart handle — the same file opened twice is two
  tables), and outlasts the table view's startup timer (`table_view.dart`, 1000 ms after the grid
  is created: focus the grid and make row 0 current) by making row 0 current itself, so the timer's
  `currentRow == -1` check skips the reset. The focus at one second still happens: a step typing
  into an input within the first second of a feature would lose it.
- **Headed and headless rasterize canvases differently.** Chrome rasterizes a 2D canvas on the GPU
  in a headed window and moves it to the CPU after enough `getImageData` readbacks (the harness
  reads pixels on every check); the first paint after the move differs in antialiasing — 2023 px of
  "repaint" on identical geometry — while headless is software throughout and never shows it. The
  config launches with `--disable-accelerated-2d-canvas`. This is why "should not have repainted"
  names the region that changed (the box in CSS px and the hit areas it touches): a failure that
  names its region is the reproduction when the machine that fails is not yours.
- **The top menu is driven by name, and a command's end is its function call's end**
  (`src/runtime/menus.ts`, `bindings/platform/commands.ts`). The Dart menu bar names its items
  `div-Bio---Analyze---Sequence-Space...` (spaces to dashes, `---` between levels), so a path needs
  no label scan. Facts that cost a run each: the bar's group opens on `mouseenter` and closes on
  `mouseleave`, so the library's leave-left `hover` CLOSES it when the item sits at the popup's left
  edge — a vertical group is entered with two moves *within* the item (`enterGroup`); a narrow bar
  folds its groups into a "more" group (the viewport is 1920; a probe at 1280 folded everything); a
  package's group appears only once the table has a column it applies to, a moment after detection,
  so the pick waits up to 5 s; Escape does not close the bar's group, moving to the page centre
  does. The command is the function call whose `Func.topMenu` equals the picked path (added to the
  core for this), watched through `onBeforeRunAction`/`onAfterRunAction`. A function-editor dialog
  keeps the call running until the work is done; a hand-made dialog returns at once and the work
  shows up as columns, which the new-column steps poll for. The pick also baselines every viewer and
  remembers the table's columns.
- **A Dart dialog button says "disabled" with a class only** — `CommandBar.setButtonActive` now
  sets `aria-disabled` too (core `ui.dart`), so `OK button in Identity dialog should be
  disabled` reads a state, not a class.
- **The Dart column selector** (`.d4-column-selector`, `ColumnComboBox`): a `mousedown` on it
  opens a `ColumnGrid` popup (a real grid, `.d4-column-grid`), typing opens the grid's search
  box, and Enter there takes the name TYPED as the column (`popup.currentColumnName =
  searchInput.value` — a partial name selects nothing). `gestures.select` does exactly that for
  it and waits for the selector to read the name. A click opened nothing in the probe, and the
  popup is not a `.d4-menu-popup`.
- **A package's init is the platform's hold on its calls**: Bio initializes on its first call
  (RDKit, the monomer libraries; 8.6 s on a fresh page) and the dialog of the first command
  waits that long; `grok.events.onPackageLoaded` fires at login for every package, not after
  init, and `Package.load()` hung in the probe. So a package binds its own readiness step
  (`the Bio package is initialized` = a call of `Bio:getSeqHelper`, which the platform holds
  until `initBio` has run). A generic "package initialized" signal in the core is still wanted.
- **A JS viewer announces itself**: `onRendered` (the runtime's `arm` subscribes to it before
  the host's `onViewerRendered`, which never fires for a JS viewer), a `get isRenderPending()`
  true from the render request to the paint (WebLogo: `_renderPending` around its debounced
  syncer; the search viewers: a counter around `renderPromise`), and `getWidgetStatus()` with
  the canvas under `parts` and hit areas in CSS px of it (WebLogo's monomer bounds are device
  px — divided by dpr, over the last laid-out range only: stale bounds of positions scrolled
  out keep old coordinates). JsViewer's Dart-backed `getWidgetStatus`/`isRenderPending` are
  overridden by the class's own members.
- **The grid reports itself** (core `grid_core.dart` `getWidgetStatus`, 2026-09-08): every
  visible cell, row header and column header as a hit area (`cell 3 of fasta`, rows as the table
  counts them from 1, `grid2table` applied), `cell type of <column>` for each visible column,
  `rows shown`, `rows`, `columns shown`, `current row`, `current column`; `{widget}` accepts
  `grid`. So a renderer claim is the grid's own word (`the "cell type of HELM string" reading of
  grid should be "helm"`) plus the pixels of the cell (`painted in at least 3 colors` — hues
  grouped by `near`, greys aside, a text cell has none), a cell click is a real click on the
  area, and the cell context menu is picked by path (`Copy > helm`). Values are reported for the
  visible column range only: a wide table (38 monomer columns) may show columns 9–38 after a
  transform, and `cell type of 1` is then "no such reading" — pick a fixture that fits.
- **A package's word about work finished off-screen is a custom platform event**
  (`grok.events.fireCustomEvent`; `src/runtime/events.ts`, `bindings/platform/events.ts`):
  listened for by id, claimed with `should have fired` (30 s, the read zeroes the count). Bio
  fires `bio-monomer-lib-loaded` after every monomer library load, with the sources loaded,
  through the bio library's `monomer-works/lib-events.ts` (`onMonomerLibLoaded`), so any package
  or test can wait for a reload instead of polling the library — and Bio now awaits the library
  update before reporting the load complete (`updateLibs` was fire-and-forget).
- **A worker pool spawns for the job, not for the machine** (ml `DistanceMatrixService`,
  2026-09-08): the diversity search's "importScripts … failed to load" errors under four workers
  were `net::ERR_ABORTED` — 30 workers spawned per `hardwareConcurrency`, 6 given work, all 30
  terminated when the 6 answered, so 24 died mid-import and Chrome reported each as an uncaught
  NetworkError on the page; alone, 5 of 30 aborted too. Now `min(threads, jobs)` workers are
  spawned when the job is known. The diagnosis: Playwright's trace records no worker network;
  a standalone page with `page.on('requestfailed')` showed the aborts, and 48 concurrent fetches
  of the chunk answered 200 in ~100 ms, so the server was never the cause.
- **A u2 name is one token**: the locator tries a phrase without its spaces and with dashes
  (`canonicalaas`, `canonical-aas`), never the spaced form — a package stamping `data-u2-name`
  from a display name writes the dashed form (Bio's collection cards, `New-Collection`).
- **Linking libraries into a package duplicates their dependencies**: `npm link` of
  `@datagrok-libraries/bio` and `ml` into Bio built with 20 "separate declarations of a private
  property" errors — each library's own `node_modules` held its own `datagrok-api` and `utils`.
  Junctions from those copies to Bio's (`New-Item -ItemType Junction`) made the types one
  again; `npm link` also restores the package's own Playwright, so `grok-bdd link` runs again
  after it.
- **A bdd page runs in simple mode, so the view tabs are hidden**: `user is logged in` sets
  `grok.shell.windows.simpleMode = true`, and the `view-handle: …` elements of every view are then
  present but hidden. U2Demo's "clicks on "U2 Demo" view" passed for days only because the rogue
  `Datagrokdsmf` autostart flipped simple mode off mid-run; with that package gone (2026-09-07) it
  failed, and became `user switches to the "U2 Demo" view` (`platform/steps.ts` `switchView`,
  `grok.shell.v = view` by name). A step that clicks a view tab is not a step on a bdd page.
- **Package autostarts run 3 s after the app started, on whatever step is running then**
  (`func_sync.dart`; found 2026-09-07 when a stand package's autostart forced `simpleMode = false`
  mid-journey). The core now exposes `grok.shell.autostartsCompleted` (a promise, added for
  this); `user is logged in` does NOT await it — the lead's call, it would cost up to 3 s per
  feature — so a feature that depends on what an autostart sets up (a package's top menu, a
  registered editor) awaits it in its own step.
- **The other five box plot specs became five journeys (2026-09-07 evening, 37 scenarios, 50 s;
  all six run in 63 s)** and each hack they carried became a core name or signal:
  `getWidgetStatus().hitAreas` gained `category <label>` / `<label> values` per category (the
  pointer-select spec used to click candidate label-band offsets until one selected a whole
  category), `p value of <group>` and `<effect> effect` under group comparison (the click that
  opens the stats in the context panel), `p value` maps to the comparison's overall result when
  the t-test box is gone; the bar chart got `getWidgetStatus` with `bar <category>`; the on-chart
  comparison selects are `input-host-method` / `control-group` / `adjustment` / `baseline` (the
  spec used to find them by their option values), the Simpson cue is `icon-simpson-warning`;
  balloons fire `AppEvents.BALLOON_SHOWN` (`d4-balloon-shown`, args type/message — the balloons
  helper was a MutationObserver over `body`), read by `no error or warning balloon should have
  been shown`, cleared at login. The data steps (`platform/data.ts`: selection, filter, rows,
  calculated columns, column colors, workspace tables) snapshot every viewer (`baselineAll`)
  before they act. `properties of {widget} should be:` is the ladder's read-back. A project
  round-trip saves through dapi with the view state of `saveLayout({saveWithData: true})` (a plain
  `getInfo()` drops the viewport) and registers its deletion with `atFeatureEnd` (harness.ts).
  `should not have repainted` is the canvas after one animation frame and one task: a mouse-over
  runs a render pass that draws the same picture, so the render count is reported, not asserted.
  Facts that cost a run each: the range slider's `max-handle` is the BOTTOM handle on an inverted
  axis (drag whichever handle is higher), and the slider lays out only once the pointer entered
  the axis strip; a viewer's title-bar close icon is `name="Close"`; `expectValue` on a `<select>`
  is the option's text, not its value; the bare p-value hover DOES show the test name (the old
  spec hovered the icon slot instead); the `T` key is `root.onKeyPress` on the viewer, so click
  the plot before pressing; the ANCOVA table's control row has no p-value (do not assert
  completeness there); `demog-1000`'s auto-picked category is DIS_POP.
- **The bar chart and the 3D scatter plot became eight journeys (2026-09-08)** — seven bar chart
  features under `features/viewers/bar-chart/` from the seven TestTrack bar chart specs, one 3D
  scatter plot feature — and what each old hack became: the bar chart reports `values` (`rows
  shown`, `bars`, `stack segments`, `clipped bars`) and gates its axes on their show flags, and
  the package's `bar-chart.ts` reads the order and the lengths of the bars and a spot no bar
  covers from the `bar <category>` hit areas (the old specs scanned the canvas for the bar green
  and clicked fixed fractions); the 3D scatter plot cannot give pixels (a WebGL canvas: no 2D
  context, no preserved buffer — `pixels()` returns an empty bitmap for it), so it reports a
  `scene signature` (the frame hashed after a render, in one task), the camera and a `point` hit
  area (the marker nearest the camera, projected), read through the generic reading steps
  (`the {string} reading of {widget} should be / differ from before / be the same as before /
  lower / higher`). `repainted` counts pixels that differ from the snapshot's bitmap now, not the
  color-histogram distance: a reorder of equal bars (Bar Sort Order) keeps the histogram and is a
  repaint all the same (the old spec had noted it as a "fault guard only"); a bitmap that changed
  size falls back to the histogram plus the area difference. The legend's 100 ms settle timer
  relays out the canvas after a legend shows or hides and was invisible to `isRenderPending` —
  "Relative Values alone is inert" saw 39000 px of that relayout as a repaint; the getter now
  includes it. `pickMenuPath` enters a group item from its left (`openGroup`): a pointer already
  resting on the item from a hover before the menu was reopened moved nowhere and opened
  nothing. Facts that cost a run each: the on-chart column selectors are named by the property
  they BIND (`div-column-combobox-split`, not the caption — `ColumnComboBox.bind` re-annotates),
  so the bar chart's are `Split` / `Value` / `Stack column input` and the 3D plot's `X` / `Y` /
  `Z` / `Color`; a selector's text has no space after the caption (`X:AGE`); the bar chart's Row
  Source defaults to Filtered, so under On Click = Filter a bar click leaves only the clicked bar
  and the Filtered Rows overlay (dark-blue outline `#0000A0`) needs Row Source = All; a fully
  negative category has no stack segment under Relative Values (its share is drawn outside
  0..1); the Bottom legend slot yields to the value selector on a short chart (assert the
  property, not the side); "Show Labels: never" removes white glyphs from inside colored bars,
  so ink goes UP, not down — a chrome toggle takes `repainted by at least N pixels`; the 3D
  plot's camera auto-rotates from the first scene until the first mouse-down or look change.
  **The backward-match round on the eight (same day, eight reviewers)** turned into core and
  library changes rather than feature patches: the bar chart reports its overlay shares as
  `selected <category>` / `filtered <category>` hit areas (an orange pixel anywhere and a blue
  outline that was there before the filter proved nothing), its `x axis` area is a frame fact (it
  echoed the property), its own debounced refresh counts as render-pending (settles resolved
  before the refresh ran and lived on timer FIFO); the 3D plot fetched its label font per scene
  and added the labels asynchronously (a signature that "differed" for any rebuild, font-XHR
  race included) — the font is now cached and pending label loads are render-pending — its
  auto-rotation is off under `immediateRendering`, and it reports `current row` and
  `highlighted rows`; a Dart menu kept its static `_lastMouseMove` across menus, so a move over
  the next menu at the same client point opened no group (`hide()` clears it); `painted` reads
  ink without moving the snapshot; `user saves the layout of the current table view to the
  server` (dapi + `atFeatureEnd` delete) is the honest round-trip where the old spec had one;
  `{string} column should have no color coding / be color-coded categorically` observe a data
  step that had no Then; the bar chart's `hang below`, `lie one under another` and the Alt-drag
  `zooms into the categories` steps replace claims the geometry never checked; every journey
  scenario ends on `no errors should have been logged` because a scenario owns its floor.
  **A settle ends on the last render the viewer announces, not the first** (found with a render
  timeline after a Stack removal: renders at 19, 131 and 243 ms with `isRenderPending` true
  throughout — the legend re-anchors in up to four passes 100 ms apart — while `settle` had
  returned at 20 ms, so the later passes landed on the next step's baseline and "Relative Values
  alone is inert" saw 78000 px). `settle` now resolves only once a render has landed AND the
  viewer says nothing is pending.
- **What a step costs, measured (2026-09-07 evening)**: the Playwright floor on this machine is
  2.4 ms for `page.evaluate`, 6 ms for `locator.evaluate`, ~2 ms for `expect.poll`, and
  `test.step`/`session.step` add nothing measurable; the trace (retain-on-failure, no snapshots,
  no screenshots) adds nothing measurable either. So a Gherkin step is the viewer's real render
  plus ~10 ms: property sets p50 43 ms (the settle waits for the paint of 1000 markers), property
  read-backs p50 13 ms, and the page is near idle after a settle (63 ms busy in a 250 ms window).
  The six box plot features: 510 steps, 33 s wall, of which one 4 s shell boot (0.8 s to first
  byte, 3.2 s of client boot on the 27 MB dev bundle, 0.2 s reset) and ~4 s of Playwright and
  global-setup start-up. Profile a slow step with a CDP `Profiler` (scratchpad `probe/profile.mjs`)
  rather than guessing: the one outlier, a numeric marker-color column at 600–700 ms, was
  `Marker.draw`'s sprite cache — a miss drew into the 16000-pixel cache canvas mid-frame and
  every following `drawImage` re-uploaded it; fixed in the core (`marker.dart` `_Cache.add`: a
  miss draws directly and the sprite lands in the cache in a microtask after the frame), 597 →
  31 ms, for every marker viewer with a continuous coloring.
- **The error floor** (`harness.ts` `watchErrors`/`takeErrors`): console errors and page errors
  from page open; `user is logged in` clears what the stand logs while booting, `resetShell`
  clears the teardown's, `journey.scenario` clears errors and balloons at each scenario's start
  (a scenario owns its floor; before 2026-09-07 (late) an earlier scenario's error failed a later one's
  check), `no errors should have been logged` reads and clears. `Failed to load resource` is not
  collected — a help page the stand does not serve is logged by the browser, not raised by the
  platform.
- **Viewer settles are armed before the change** (`writeProperties`, `resize`): the repaint an
  action causes lands on the next task, so a settle subscribed after the action already missed it
  and burns its cap (that was 1.5 s per resize).
- **A settle ends when the viewer says nothing is pending, not when a cap runs out**
  (`viewers.ts` `settle`/`quiet`, 2026-09-07 (late)): `viewer.isRenderPending` (core
  `ViewerBase.isRenderPending`: a `debounced` timer armed, a resize the poll has not handled,
  `_invalidateRequested`) is polled every task until false or a render lands; a property that
  paints nothing returns at once, a repaint is waited for as long as it takes, one pending for
  10 s is a platform failure and is thrown. The 300 ms cap remains only for a viewer without the
  signal (a JS viewer). The old silent cap let a late repaint from step N satisfy step N+1's
  "repainted". The data steps (`platform/data.ts` `changeTable`) baseline every viewer, act, then
  `settleAll` — so the next step's baseline is the state after the change. Found while wiring it:
  the resize path is a 100 ms polling `Resizer`, invisible to `_invalidateRequested` — hence
  `Resizer.isResizePending` and `parent` on the Resizer.
- **"Before" is before the last change, never after the last check** (2026-09-07 (late)): the
  than-before family (`repainted`, `less/more ink`, area ink, highlight, value range, color scale)
  does NOT move the snapshot when it passes; every change (property set, menu pick, gesture on a
  hit area, data step, `takes a snapshot`) takes it. Before, `repainted` re-snapshotted after
  passing, so `Then repainted / And the "M values" area more ink than before` compared the
  after-state with itself (0 px difference, a false failure). The snapshot holds the histogram,
  the ink of every hit area, the value range, the color scale's range and the viewer's `values`.
- **What a check means, after the review of 2026-09-07 (late)** (six reviewers backward-matched the six
  box plot features against the specs they replaced; findings in the `/bdd-translate` skill):
  `repainted` is a change detector — one pixel — and stays one; a claim about a shape says it
  (`the "stats" area … more ink than before`, `should have a "stats" area`, `should contain the
  color`, `areas … painted in different colors`, `should show fewer rows than before` from the
  `rows shown` reading); a chrome toggle takes `repainted by at least N pixels`. Highlight
  checks carry a margin (`highlightMargin`: max(200, 2 × selected rows) × dpr², capped at a quarter
  of the view) — the mouse-over halo alone moves the hue count by hundreds. Colors compare by hue
  (`near`: ±20°, greys by lightness) because markers are drawn with alpha, so `#00FF00` lands as
  `#A6FFA6`. `should be bound to table` reads `viewer.dataFrame`, not the Table property. The
  positive data steps fail on a category no row has (a typo matched "0 of 0"). Negative checks
  (`should not have repainted`, `the same value range as before`) read after `quiet`.
- **The tooltip is one element, hidden between hovers** (2026-09-07 (late)): `contain text` on it passed
  on the previous hover's text. Core: `Tooltip.hide()` empties `content`; binding: `expectText`
  on the `tooltip` kind matches visible elements only. A viewer's `values`
  (`WidgetStatus.values`, JS `IWidgetStatus.values`) are named readings a test compares before and
  after — the box plot reports `rows shown` and `color scale min`/`max` (the range the scale
  labels, the filtered rows' when they narrow it).
- **Hit areas say what is shown**: the box plot reports `x axis`/`y axis` only while
  `showCategoryAxis`/`showValueAxis` are on (the box existed, hidden, before), and the control
  band as `control band` (pooled) or `control band <stratum>` for every stratum of a matched
  band (adjacent strata share one band and it answers to each name). The covariate selector
  (`Adjust by`) is hidden by design with two category levels (`covariateSelectorVisible` needs
  `_singleCat`) — do not assert it visible in a two-level scenario.
- **Viewer event subscriptions have a lifetime** (`viewers.ts` `listen`/`unlisten`/`forget`): one
  subscription per viewer and event, replaced by a repeated "listens for", ended by the "should have
  fired" read or by `grok.events.onViewerClosed` (which also drops the render stamp). Before
  2026-09-07 every "listens for" stacked another subscription that lived as long as the viewer.
- **The second viewer round (2026-09-08: grid, scatter plot, histogram, form, PowerGrid's forms,
  filter panel, legend; design record `core/docs/features/ui2/automation/VIEWERS_ROUND_2026-09.md`)**
  added to the runtime before any feature was written: hit areas are anchored at `parts.canvas`,
  else `parts.root`, else the viewer's root (`anchorOf`), so a DOM viewer reports areas and
  readings and only the pixel steps need a canvas; every pixel reading composites the viewer's
  `overlay` part over its canvas when both have the same size (`pixelsOf` — the scatter plot
  draws regression lines, labels and stats on its overlay, the grid its selection and current
  cell; before, `repainted` and `ink` were blind to them); the snapshot keeps every hit area's
  rectangle (`taller/wider than before`) and the legend's mode, slot, keys and size; readings can
  be remembered by viewer type (`rememberValue`, the value-range pair generalised); the legend is
  read from the `data-legend-*` attributes it already publishes and its items are a kind, never
  hit areas (a tooltip-hosted legend leaves the viewer root — `legendRoot` looks in `.d4-tooltip`
  too); a legend item click and a splitter drag snapshot first and `settle` after, since the
  common click does neither. Kinds `legend item` and `filter card`, the `filter panel` element
  (accepted by `{widget}`) and the states `partially checked` / `invalid` came with it. The
  ownership rule of the round: one widget reports everything it contains — the filter panel's
  categories, checkboxes and handles are ITS areas and readings, its inner grids are named
  `filter-grid`, not `viewer-Grid`, so the reserved `grid` stays unique.
- **Keys pressed "in" a viewer go to the viewer** (`gestures.ts` `editorOf`, 2026-09-09): `focuses on`
  and `presses {key} in` resolve an element's editor — itself when editable, else the first input
  inside — and for a viewer that first input is a field, not an editor: the form's arrow-key row
  walk never fired because the arrows went into a read-only field input. A target matching
  `[name^="viewer-"], .d4-viewer` is now its own editor (its root, which the form makes focusable
  with `tabIndex = 0`). The probe that found it: the same key through `locator.press` on the root
  moved the row; through the library's step it did not.
- **Global setup fetches the client before a browser does** (`global-setup.ts` `warmClient`,
  2026-09-08): a dev stand serves the Dart client from `pub serve`, which recompiles the whole
  bundle after any source edit — the first request then blocks for minutes, and every run started
  in that window died in global setup with `waiting for locator('[name="Browse"]')` after 60 s.
  One plain `fetch` of `login.dart.js` and its deferred part pays that wait once, outside the
  browser and with no cap, and costs milliseconds when nothing is compiling. It is not a wait for
  a signal the platform owes us: the compile is the dev environment's, and the alternative is a
  60 s cap that fails a good suite while someone edits Dart.
- **A menu group needs a signal, not a hover** (2026-09-09): `openGroup` used to move the pointer
  over the group and return, so the next lookup read the parent menu, where the wanted item is
  present but hidden — the click then timed out on an element the error message listed as being
  there. It now waits until the item it was asked for is visible, and tries **every** item that
  carries the group's label: a viewer menu holds two groups called "Annotations" (the axis one and
  the viewer's own), and `.first()` picked whichever came first in the DOM. A group the menu
  renders inline, or one already open, returns at once.
- **A hidden thing keeps the geometry of the frame that drew it** (2026-09-09): an annotation
  region hidden by `showViewerAnnotationRegions = false` still had `polygonScreen` from its last
  render, so its hit area would have kept answering. A status reports an area only while the thing
  is drawn; and a count of what the look *holds* is read off the look, not off the parsed list the
  viewer keeps for drawing, which the same flag empties.
- **`@known-failure`** (2026-09-09): a scenario tagged with it is expected to fail — its failure
  does not fail the test, and its **passing** does, with "the bug it describes is fixed, so the tag
  has to go". An open bug is translated honestly and the suite still says something when it is
  fixed; nothing is softened to stay green.
- **The overlay was invisible to every settle** (2026-09-09): `CanvasViewerMixin.invalidateOverlay`
  sets `_invalidateOverlayRequested`, which `ViewerBase.isRenderPending` did not read — so a
  current-row, hover-line or row-group repaint was pending work no step waited for (the scatter
  plot worked around it with a private flag). It is in the base getter now.
- **Codegen** emits names, never selectors; `\n` line endings, no timestamps, EOL-normalized drift
  check; orphaned generated files are removed on compile and reported on `--check`.
- **Session readiness**: `user is logged in` skips navigation when the page is already in the shell
  (chained scenarios), calls `grok.shell.closeAll()` and waits for the Home view to be current
  (`grok.shell.v.type === 'datagrok'`) — closeAll re-adds Home asynchronously.
- Do not reinvent what `@datagrok-libraries/test/src/playwright/*` has (storage-state login,
  `openTableFromFile`, base config): import it (compiled `.js` + `.d.ts`, `.js` suffix).

## Environment facts that cost time to learn

- Dev setup on this checkout: the library is junction-linked into U2Demo
  (`packages/U2Demo/node_modules/@datagrok-libraries/bdd` → `libraries/bdd`, and
  `node_modules/@playwright/test` → the library's copy); create junctions with PowerShell
  `New-Item -ItemType Junction` — Git Bash mangles `mklink /J`. Run `npm run build` in the library
  after engine or binding changes; the U2Demo project is driven with
  `node ../../libraries/bdd/bin/grok-bdd.js <command>` from `packages/U2Demo` (a global install /
  `npm link` gives plain `grok-bdd`).
- Stand: `http://localhost:8888` (nginx over datlas `:8082`), admin/admin, developer key `admin`.
  When every login "hangs", the hand-started Postgres container (`affectionate_einstein`, no restart
  policy) died with a Docker restart: `docker start affectionate_einstein`.
- To publish U2Demo: `cd packages/U2Demo && npx webpack && grok publish localhost --key admin
  --skip-check`. Webpack and `grok check` ignore `bdd/`; U2Demo un-ignores `src/raw.d.ts`.
- U2Demo/u2 builds: `@datagrok-libraries/u2` is NOT on npm; `datagrok-api` in `libraries/u2` and
  `packages/U2Demo` node_modules must be declaration-only copies of `public/js-api`; `grok api`
  rewrites `package-api.ts`/`package.g.ts` with LF endings (content-identical; `git checkout --`).
- In the U2Demo app view the platform `.d4-toolbox` is present and visible but empty; the app's
  navigation is a u2 splitter panel of its own. U2Demo's `bdd/bindings/demo.ts` makes the
  sub-demo content pane (`.u2demo-content`) the `U2 Demo` context so page trees/lists win over
  the navigation tree (`demo navigation` by name), registers a package kind `readout` for the
  pages' "name = value" lines, and opens a sub-demo with `user opens the "<label>" demo page`
  (cold page: `page.goto`; warm page: `grok.functions.call('U2Demo:u2DemoApp', {path})` +
  `grok.shell.addView`, so tables opened by earlier steps survive — a `goto` reloads the client;
  the shared helper is `openSubDemo(page, route, ready)`, which the MSA workbench step uses too).
  One feature per sub-demo lives under `bdd/features/demo/<group>/`.
- **Live DOM outline probe**: before writing phrases for a page, dump what the engine will see.
  Recipe (scratch dir): `createRequire('<libraries/bdd>/package.json')('@playwright/test').chromium`,
  `newContext({storageState: '<pkg>/bdd/e2e/.auth.json'})`, goto the page, then `page.evaluate` a
  walker over the root printing tag, classes, `data-u2*`, `role`, `aria-*`, `name`/`placeholder`,
  input values and own text per element (mark `offsetParent === null` as hidden); for popups dump
  `[data-u2="menu"], [data-u2="dialog"], [data-u2="tooltip"], [data-u2="notify"], [data-u2="tour"],
  .d4-dialog, .d4-balloon` after the opening action. Pass the walker as a string expression that
  CALLS the function (`page.evaluate(\`(${fn})(${JSON.stringify(args)})\`)`) — a bare function
  string is evaluated, not invoked. u2 popups are portaled under `body > .u2-overlay` with
  `data-u2-owner` = the nearest NAMED ancestor of the trigger (the app root, "U2 Demo", when the
  input has no name), so the owner edge rarely helps inside the demo — the whole-page fallback does.
- Facts the probe established (2026-09-03): u2 ChoiceInput is a native `<select>` with a leading
  empty option; BoolInput `switch` is `span[role=switch][aria-checked]` (no input); NumberInput's
  editor part wraps `span > input[type=text]`; TextInput `search` wraps the input in the editor
  part; icon/function/columns pickers are `[data-u2-part=editor][role=button][aria-haspopup]`;
  Combobox has no input root (`[data-u2=combobox] > input[role=combobox]`) and opens on input or
  ArrowDown; MultiSelect's field is `[role=combobox]` opening a `.u2-multi-select-popup` with
  `[role=option]` rows and a `Select all` `[role=button]`; ButtonGroup items are `button`s
  (`role=radio` + `aria-checked` in single-toggle mode, `aria-pressed` in multi); toolbar toggles
  carry `aria-pressed`; menu items are `.u2-menu-item[role=menuitem]` with `.u2-menu-label` and
  `.u2-menu-shortcut`, `aria-disabled`, submenu `aria-haspopup`; menu-bar items are
  `button.u2-menu-bar-item[role=menuitem]`; breadcrumbs are `button.u2-breadcrumbs-item` +
  `span.u2-breadcrumbs-current[aria-current=page]`; accordion panes `.u2-accordion-pane` with a
  `[role=button][aria-expanded]` header and `.u2-accordion-title`; Section header
  `.u2-section-header[role=button][aria-expanded]`; Wizard steps `li.u2-wizard-step` with
  `.u2-wizard-title` and `aria-current=step`, buttons BACK/NEXT/FINISH; Tour overlay
  `[data-u2=tour]` + sibling `.u2-tour-popup[role=dialog]` with SKIP/NEXT/DONE and a "1 / 4"
  counter, `onFinish` gets `done`/`skipped`; notify balloons `.u2-notify[role=status|alert]` with
  aria-label Close/Copy icon buttons; BasicTable is a real `<table>` (`tbody tr[aria-selected]`);
  VirtualGrid/VirtualList are `[role=listbox]` with `[role=option]` cells/rows (`.u2-list-row`,
  `.u2-grid-cell[title]`); VirtualTree rows `[role=treeitem][aria-expanded]` with
  `.u2-tree-twistie`/`.u2-tree-label` (the twistie toggles, the row selects); PropertyGrid rows
  `.u2-propgrid-row` (`.u2-propgrid-name` + inline inputs) under `.u2-propgrid-category[role=button]
  [aria-expanded]`; Card `[data-u2=card]` with `.u2-card-title`, `role=button` when clickable,
  `aria-pressed` when selectable; StatCard `.u2-stat-label`/`.u2-stat-value`; ProgressBar
  `[role=progressbar]` + `.u2-progress-percent` + `.u2-progress-description`; MessageInput editor
  `.u2-msg-editor[contenteditable][role=textbox]` + `button.u2-msg-send`; the Dart `ui.dialog`
  is `.d4-dialog[role=dialog][name="dialog-<title>"]` with `.d4-dialog-title` and
  `button.ui-btn[name="button-OK"]`; a bridged u2 input inside it keeps `data-u2="text-input"` on
  the `.ui-input-root` with a `.ui-input-label`; the shell status bar is `.layout-status-bar` with
  the view's panels in `.d4-view-status-panel`; a platform "Debugging packages" `.d4-balloon.warning`
  is up after every reload (text checks on `notification` are any-of, so it does not interfere).
- Node prints "Importing JSON modules is experimental" whenever tsx registers; the launcher
  filters that one warning.
- **VS Code Cucumber extension** (`cucumberopen.cucumber-official`). It parses
  `defineParameterType({name, regexp: /…/})` from the glue, so `cucumber.parameterTypes` holds only
  `state` (whose regexp is built from a list in code); listing the others again yields "already a
  parameter type". Settings live at the core root (gitignored), in `packages/U2Demo/.vscode` and
  here. **A settings JSON written through a Bash heredoc lost the `\\` of `\\w`**: VS Code's lenient
  parser kept `[w -]+?`, so `{viewer}` steps showed as undefined while the compiler resolved them —
  write JSON with the Write tool and validate with `json.load`. To see what the extension resolves,
  build an `ExpressionBuilder` over `@cucumber/language-service`'s `WasmParserAdapter` (the Node
  adapter has no TypeScript grammar) and match against it.
- Bash tool here: cwd persists across calls and very long heredocs fail to parse — use absolute
  paths and the Write tool for whole files.
- **u2 finding, not a harness bug (2026-09-03, unreported)**: in a `funcForm` a number field
  typed key by key starting with "-" ends up EMPTY — the transient "-" parses to null, the form
  writes null into the FuncCall, and the param's `onChanged` echo writes null back into the input
  after the next keystroke landed (probe: `scratchpad/probe/numinput.mjs`; "7" sticks, "-7" does
  not, `fill('-5')` does). A human loses the sign the same way. The U2Demo features use positive
  numbers meanwhile; do not add a `fill` fallback to `typeInto` to hide it.
- Two phrase traps met while writing the U2Demo features: a qualifier starting with an ordinal
  word (`First name input` → "first" + `name input`) needs quotes (`"First name" input`), and a
  label equal to a kind name (`Columns input`) now resolves to the labelled input first (nouns.ts
  puts the whole-kind reading last).
- **Viewer facts true of every viewer** (a specific viewer's areas and readings are in its own
  `core/client/d4/lib/src/viewers/<viewer>/CLAUDE.md` "Automation surface"): `Property.caption` in
  the JS API is the raw name (`showInsideValues`) unless `@Prop(name:)` set one, and the runtime
  resolves a caption, a name and a `…ColumnName` by normalizing both sides — two properties sharing
  a caption (a plot's `Min`/`Max` for X and Y) must be named, not captioned. Hit areas come back in
  **canvas coordinates** and the runtime adds the canvas's client rect. The on-canvas column
  selectors (`div-column-combobox-*`) are hover-revealed — `visibility: hidden` until the pointer is
  over the viewer, and again after a property change — while `display: none` is a disabled selector
  and stays hidden. Property rows are `tr.property-grid-item[name="prop-marker-type"]` with
  `aria-disabled` while `dependsOn` gates them; menu items are `role=menuitem`,
  `name="div-Misc---Show-Inside-Values"`, own label a direct `.d4-menu-item-label`, and a checked
  one carries `aria-checked`.
- `user filters rows where ...` writes the table's filter bitset directly (`df.filter.init`). That
  write is transient: anything that calls `dataFrame.rows.requestFilter` recomputes the filter from
  the registered filters and drops it. A histogram does that on **every context-menu pick**
  (`_refresh()` while `filteringEnabled`, its default), so a menu pick after such a step silently
  restores all rows. When a scenario must hold a filter across viewer interaction, filter through a
  filter card (`user adds a categorical filter on ... keeping ...`), which is a registered filter.
- `{widget} should have repainted` measures the whole canvas even when the phrase names an area —
  use `the {string} area of {widget} should have repainted` for one area's own pixels. A hover
  highlight can be gone by the time a later step reads it (the tooltip lands under the pointer), so
  put the repaint checks right after the gesture and the tooltip text after them.
- **Reviews found these, and the fixes are load-bearing** (the process itself is `/bdd-translate`):
  a per-area repaint claim measured the whole canvas; `remembered value range` compared top and
  bottom only, so a Reset View that lost the X window passed; the tooltip readers matched a hidden
  tooltip's leftover text; `expectLegendPlacedAsBefore` ignored size; and `the filter panel should
  have N filters` called `getFiltersGroup`, which CREATES a panel when the view has none, so
  "0 filters" quietly resurrected one. Two steps exist because a claim had no honest phrase:
  `the open tableview should have {int} {viewer} viewer(s)` (an "added"/"closed" claim is vacuous on
  a view that already owns one of that type) and `no|all rows where {string} is {string} should pass
  the filter` (a row COUNT cannot tell a card that filters from a card that stopped).
  Write the fix, run it, and keep what the product says.
- Two core defects the review pass surfaced: `grid.cell2screen` returns null for a cell scrolled out,
  so the grid's status emitted a null rect and the filter panel's status (which reads the inner grid's
  areas) crashed with `TypeError … reading 'get$left'` whenever a card was re-added — both guard now;
  and the filter panel's counter host started visible with an empty count until the first
  `_refreshCounts`, so an empty panel showed an empty badge (`htmlSetVisible(false)` at creation).
- A grid's selected columns are a set: `columns "..." should be selected` compares regardless of
  order. `the "current row" reading of grid` is 1-based, 0 meaning none, and Escape clears the row
  and column selection but leaves the current cell where it was.
- Probe scripts for these live in the scratchpad (`probe/boxplot.mjs`, `selectors.mjs`,
  `submenu.mjs`): dev-key token via `POST /users/login/dev/admin`, cookie + localStorage `auth`.

### Platform behaviour these features had to learn (2026-09-09)

- `RowSet.SelectedOrCurrent` is the selection **when it has anything**, and the current row only
  otherwise — not the union. `MouseOverGroup` is `dataFrame.highlight`, which a hover elsewhere
  leaves standing, so a scenario that wants "empty" must state what it hovered last.
- A viewer's `rows shown` is `combinedFilter.trueCount` everywhere except the viewers that mean
  "what the frame drew": the scatter plot (markers past the filter, drawable on the current axes,
  inside the viewport), the PC plot under a transformation (the aggregated frame), and the density
  plot (the rows the last binning pass counted — on demog-1000 with Y=HEIGHT that is 872, not 1000,
  because a blank is never binned).
- `demog-1000` has 128 rows with a blank HEIGHT, so a scatter plot on HEIGHT draws 872, not 1000.
  Pick columns with no blanks (AGE, WEIGHT) when the claim is about row sources rather than blanks.
- The grid reports a cell's colour as **lowercase** `#rrggbb` (`htmlColor` → `toRadixString(16)`).
- A linked colour coding is two tags — `.color-coding-type = Linked` and
  `.%color-coding-linked-column-name` — and "Apply to: Text" is the *grid column's*
  `isTextColorCoded`, not a column tag. Switching a colouring off rewrites only the type tag, so
  the colours it stored survive and come back when the type is switched on again.
- `grok.shell.t` is the current view's table, not the viewer's: after a viewer is rebound to
  another table, a column step still reads the view's table.
- PowerPack's Formula Lines dialog opens 500 ms after a region is drawn, from a bare `Timer`; the
  region itself is in the look synchronously, so the claims about it need no dialog. Closing the
  dialog with OK keeps it.

### Platform behaviour the third viewer round had to learn (2026-09-10)

- **A zoom is animated, and it used to win.** `CanvasViewportMixin.zoom` walks the viewport over
  ten 50 ms frames; a viewport written from outside during that window — Reset View, a slider, a
  bound property — was silently overwritten by the frames still to come, so Reset View immediately
  after an Alt-drag did nothing. The animation now stops when it sees the viewport moved under it,
  and `isRenderPending` reports a running zoom (it did not). Both fixes are in
  `viewer_base/canvas_viewport_mixin.dart` and `viewer_base/viewer_base.dart`.
- **A hidden group label is still a match.** The Dart menu mirrors every property under a zero-size
  `Properties...` submenu, so `Columns`, `Tooltip`, `Correlation Type` each occur twice in one
  popup. `openGroup` waits on the first **visible** candidate; before that it waited on
  `.first()`, which was whichever came first in the DOM.
- **`painted in at least {int} colors` groups by hue**, so a linear colour scale (five shades of one
  blue) counts as one. Correct for markers, useless as a "this is colour-coded" precondition — use
  `should be painted` or compare two areas.
- **A package viewer's surface is not on the stand until the package is published**, and when the
  viewer lives in a library (the Forms viewer is in `libraries/utils`, registered by PowerGrid) the
  package's `node_modules/@datagrok-libraries/<lib>` must point at the checkout or it compiles the
  registry copy. The failure looks exactly like a mistyped phrase.
- **Viewer facts worth knowing before writing a claim.** The density plot recolours its bins rather
  than drawing a selection overlay, so it shows no selection highlight; `bins = 1` draws **three**
  hexagon bins, because hexagon rows interlock (`getCountsIndex` adds `ceil(y / 2)`) — one bin is
  one bin only under rectangles. `descriptionPosition` defaults to **Top**. A calendar at 1920×1080
  lays out about 12 weeks, so a claim about all 541 days of demog-1000 is unreachable — pin an early
  date or claim relatively. The heat map's `row height` depends on the docked layout (0.86–0.98
  measured), so bound it rather than pinning it.
- **`DensityPlotViewer` had no `viewport` in the JS API** although the Dart viewer has one and
  `ScatterPlotViewer` and `BoxPlot` both publish it — every value-range step read `undefined`.
  Added to `public/js-api/src/viewer.ts`; a js-api rebuild lands in `xamgle/web/js/api` and needs no
  pub serve restart.

### What the flakiness tail turned out to be (2026-09-10)

Six full runs landed at 118-120 of 121 with a different feature failing each time, and every one of
them passed in isolation. Three causes, all found by reading what the failing run recorded rather
than by re-running:

- **The column picker does not take a column that is already current.** `ColumnComboBox`'s search
  box sets `popup.currentColumnName` on Enter and the popup closes from
  `dfColumns.onCurrentRowChanged` — a frame whose current row does not move announces nothing, so
  Enter on the name already highlighted left the picker open and took nothing. The pivot table's
  journey re-adds the column it has just removed, which is exactly that case. Fixed in
  `d4/lib/src/common/column_combo_box.dart` (the commit is now a method both paths call). The
  binding no longer swallows the picker's failure to close: a picker still open after Enter took no
  column, and the reading it was for timed out 15 s later three steps away.
- **A dock layout pass takes an inline size back** — see the resize rule above.
- **A viewer that empties its root before it can lay out has no picture and no message.** The word
  cloud's `render` emptied the root and then read `this.root.parentElement!.clientWidth`; on a host
  the dock had not sized yet that throws, leaving no canvas, no error and `isRenderPending` stuck
  true, which the status reports as a viewer with no areas at all. It now keeps the picture that is
  up and stays pending until a host it can measure comes back.

## Conventions

TS strict, 2-space, single quotes, `.js` suffixes on relative imports (NodeNext), comments near
zero. Nothing is committed without the lead's order; stage by explicit path — the tree carries
other sessions' work.
