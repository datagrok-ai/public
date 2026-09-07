# @datagrok-libraries/bdd — behavioral automation (Gherkin → Playwright)

Feature files bound to the u2/platform vocabulary, compiled deterministically into Playwright specs
that packages own. Design record: `core/docs/features/ui2/automation/BRAINSTORM.md` (rulings by the
lead: standard Gherkin; committed, drift-gated codegen; u2-centred; overridable composition;
generic kinds; reserved platform names; packages own their tests; one page per feature folder; a global
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
box plot specs under `files/TestTrack/Viewers/BoxPlot/` and their helpers, 36 s for all six on one page).

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
- **One page per feature folder** (`src/runtime/harness.ts`, the lead's rule 2026-09-07: "one
  folder, one tab"): `feature(test)` registers `afterEach` (leave the context, `resetShell`) and
  `afterAll` (the feature's `atFeatureEnd` cleanups); the page is created inside the first test
  (`session.page(browser)`), so Playwright merges the project's context options (storage state,
  viewport), and kept in the module-level `shared` with its folder — the folder's next feature
  file finds it and `user is logged in` only resets the shell (0.1 s instead of the ~4 s boot);
  a feature from another folder closes that context first. Playwright 1.62 starts a trace chunk
  on every existing context at each test start (`ArtifactsRecorder.willStartTest` →
  `didCreateBrowserContext`), so every test still gets its own trace; a failed test restarts the
  worker, which drops the page. NEVER leave several Datagrok pages open in one browser: six live
  clients made every step 2–3× slower (measured 2026-09-07, 124 s for the suite). The generated
  spec calls `test()` itself so reports point at the spec line, not the harness; every step is
  `session.step(line, title, fn)` (`feature(test, "features/x.feature", import.meta.url)`), a
  Playwright step whose `location` is the feature line.
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
- **Hover is two pointer events, and never sleeps**: the gesture leaves the element to its left on
  the same line, lands on its centre in one move, then checks that the element is still where it
  was (a view still docking moves it out from under the pointer). Every pointer event costs a
  frame (~16 ms; a six-step move was 100 ms), and the animation-frame wait it had was another two
  frames for nothing: hover-driven layout is synchronous. Leaving upwards would cross the
  neighbouring row and close the submenu the item sits in. The stepped move existed because a
  Dart menu group did not open on the first `mousemove` when another group's submenu state was
  stale — fixed in the core (`menu.dart` `_initItem`: the move that closes a sibling's submenu
  opens this one; `hide()` clears `_expandedItem`), 2026-09-07. Scrolling into view only when
  the box is outside the viewport.
- **A step ends when the platform is done, not when the DOM shows** (`platform/steps.ts`
  `openDataset`): opening a table starts semantic-type detection in the background (package
  detectors — Chem's SMILES detector took 300–500 ms on spgi-100), which used to land on
  whichever step came next (a 500 ms "Table" read, a 500 ms property set). The step now waits
  for the platform's own `ddt-semantic-type-detected` event for that data frame (identity by the
  Dart handle — the same file opened twice is two tables).
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
- **The error floor** (`harness.ts` `watchErrors`/`takeErrors`): console errors and page errors
  from page open; `user is logged in` clears what the stand logs while booting, `resetShell`
  clears the teardown's, `no errors should have been logged` reads and clears. `Failed to load
  resource` is not collected — a help page the stand does not serve is logged by the browser, not
  raised by the platform.
- **Viewer settles are armed before the change** (`writeProperties`, `resize`): the repaint an
  action causes lands on the next task, so a settle subscribed after the action already missed it
  and burns its cap (that was 1.5 s per resize).
- **Viewer event subscriptions have a lifetime** (`viewers.ts` `listen`/`unlisten`/`forget`): one
  subscription per viewer and event, replaced by a repeated "listens for", ended by the "should have
  fired" read or by `grok.events.onViewerClosed` (which also drops the render stamp). Before
  2026-09-07 every "listens for" stacked another subscription that lived as long as the viewer.
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
- **VS Code Cucumber extension** (`cucumberopen.cucumber-official` 1.11.0, language-server 1.7.0;
  the lead also has `alexkrechik.cucumberautocomplete`, unconfigured). It parses
  `defineParameterType({name, regexp: /…/})` calls from the glue, so `cucumber.parameterTypes` in
  settings holds only `state` (its regexp is built from a list in code); listing the others again
  yields "already a parameter type" errors. Settings files live at the core root (gitignored), in
  `packages/U2Demo/.vscode` and here. **A settings JSON written through a Bash heredoc lost the
  `\\` of `\\w`** (2026-09-03): VS Code's lenient parser kept `[w -]+?`, so `{viewer}` and
  `{dataset}` steps showed as undefined ("scatter plot viewer should be added …") while the compiler
  resolved them — write JSON with the Write tool and validate with `json.load`. To see exactly what
  the extension resolves: `npm i @cucumber/language-service` in a scratch dir, `new
  WasmParserAdapter('<pkg>/dist')` (the Node adapter has no TypeScript grammar), sources
  `{languageName: 'tsx', uri, content}`, `new ExpressionBuilder(adapter).build(sources,
  parameterTypesFromSettings)`, then `expressionLinks[].expression.match(stepText)`.
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
- **Viewer facts (box plot, 2026-09-06)**: `Property.caption` in the JS API is the raw name
  (`showInsideValues`) unless the Dart `@Prop(name:)` set one (`Min`, `Adjust By`); the viewer
  runtime resolves "Show Inside Values", "Category 1" (`category1ColumnName`), "Marker Size Column"
  (`markerSizeColumnName`, not `markerSize`) by normalizing both sides. Hit areas come back in
  canvas coordinates (`getWidgetStatus().hitAreas`, `Rect.toMap` → `{x, y, width, height}`), the
  runtime adds the canvas's client rect. The on-canvas column selectors (`div-column-combobox-*`)
  are hover-revealed: `visibility: hidden` until the pointer is over the viewer and again after a
  property change — hover the viewer before asserting them visible; `display: none` (a disabled
  selector) is hidden regardless. The property panel's rows are `tr.property-grid-item[name="prop-
  marker-type"][aria-label="Marker Type"]` with `aria-disabled` while `dependsOn` gates them. The
  vertical range slider is `svg[type="range-slider"][name="y-slider"]`, `max-handle` at the top.
  Dart context-menu items: `role=menuitem`, `name="div-Misc---Show-Inside-Values"`, own label a
  direct `.d4-menu-item-label`; the stats-region menu names drop the "Statistics" prefix.
- Probe scripts for these live in the scratchpad (`probe/boxplot.mjs`, `selectors.mjs`,
  `submenu.mjs`): dev-key token via `POST /users/login/dev/admin`, cookie + localStorage `auth`.

## Conventions

TS strict, 2-space, single quotes, `.js` suffixes on relative imports (NodeNext), comments near
zero. Nothing is committed without the lead's order; stage by explicit path — the tree carries
other sessions' work.
