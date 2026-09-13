# @datagrok-libraries/bdd — behavioral automation (Gherkin → Playwright)

Feature files bound to the u2/platform vocabulary, compiled deterministically into Playwright specs
that packages own. Authoring guide: `README.md`. Design record and rulings:
`core/docs/features/ui2/automation/BRAINSTORM.md`; the round-by-round history and the postmortems:
`core/docs/features/ui2/automation/VIEWERS_ROUND_2026-09.md`; the u2 selector contract:
`core/docs/features/ui2/AUTOMATION.md`. A viewer's areas and readings are in its own
`core/client/d4/lib/src/viewers/<viewer>/CLAUDE.md` "Automation surface". Translating old specs is
the `/bdd-translate` skill.

## Layout

```
bin/grok-bdd.js         launcher: the package-local install first, else this one; runs dist/src/cli.js
src/project.ts          a project = a dir with features/ (a package's bdd/, or the library root); bdd.config.json {tiers}
src/discover.ts         imports binding modules (library: this build's bindings/; project: .ts through tsx), maps step fn → export
src/gherkin.ts          @cucumber/gherkin → flat model (And/But resolved, outlines expanded, rule backgrounds folded)
src/registry.ts         Given/When/Then, element(), context(), kind(), dataset(), defineParameterType(); one global registry
src/nouns.ts            phrase (+ context) → NounRef (pure; shared by compiler and runtime)
src/match.ts            cucumber-expressions matching + specificity
src/compile.ts          FeatureModel → *.test.ts (feature(test) session, test() per scenario, @journey = one test)
src/states.ts           the {state} list, shared by the assertions, the parameter type and init's VS Code settings
src/init.ts, src/cli.ts init | compile [--check] | lint | list-steps | run [playwright args] | link
src/runtime/            args, locate, gestures, assertions, harness (session, journey, resetShell, error floor),
                        viewer-runtime (in-page window.__bdd), viewers (readers over it), viewer-pixels,
                        viewer-menus, viewer-legend, menus (top menu), events, functions, patience, failure
bindings/common/        parameter-types, kinds (every u2 data-u2 kind + Dart conventions), steps, session — always loaded
bindings/platform/      the shell: elements, datasets, steps, data, columns, commands, functions, events — always loaded
bindings/tiers/viewers/ opt-in: steps (properties, menus, areas, pixels, legend, events, floor), widgets (shared per-viewer steps)
tests/                  node:test via tsx: nouns, compile, project, init, failure, locate (Chromium over a static page)
playwright.config.ts    the one config every project runs with (BDD_ROOT → testDir/outputDir/storageState; 4 workers)
```

## The principle

We test our own platform, not a black box. **When a test would wait, sleep, scan pixels or retry,
a signal or a name is missing in the core, and the fix goes there** (`d4` viewers, `xamgle`, the
js-api) — never a `waitForTimeout` in a step. What was added that way: `viewer.immediateRendering`
and `isRenderPending`, `onContextMenuShown/Closed`, `getWidgetStatus().hitAreas/values/parts`,
`aria-disabled` on menu items, property rows and dialog buttons, `Func.topMenu`,
`d4-balloon-shown`, `grok.shell.autostartsCompleted`, `Resizer.isResizePending`, `data-legend-*`.

## Invariants — what must not regress

- **One registry, through `dist/`.** Specs import the library by package subpath, project bindings
  by relative path; never mix `src/` and `dist/` in one run. ESM everywhere. One `@playwright/test`
  per run: packages depend by path (`file:../../libraries/bdd`) and `grok-bdd link` makes the
  package's Playwright the library's copy (redo after `npm ci`).
- **One page per worker** (`harness.ts`): `feature(test)` reuses the worker's page, `afterEach`
  resets the shell (Escape for dialogs and menus, `ui.tooltip.hide`, notices removed, `closeAll`,
  Home current), `afterAll` runs the feature's `atFeatureEnd` cleanups. Never open several
  Datagrok pages in one browser.
- **`user is logged in` only resets when the page is in the shell**; it sets `simpleMode` (view
  tabs hidden — switch views by name), clears the error and balloon floors, installs the in-page
  runtime. Package autostarts land 3 s after boot; a feature that needs one awaits
  `the package autostarts have completed`. The first shell load of a page waits 180 s and warns
  past 30 s: a starved `pub serve` hands out the 32 MB bundle in a minute, and that is a delay
  once per page, not the feature's failure.
- **A row test is checked before the rows are scanned** (`viewer-runtime.ts` `checkTest`): a
  value the column does not hold (a typo) fails naming the values it has, a range over a text
  column fails; an empty result is legal for a filter or a selection (a range that keeps no row
  empties a chart on purpose), and a claim about rows none of which exist fails.
- **A dataset is read once per page and every feature gets a clone** (`platform/steps.ts`
  `openTable`): the clone keeps the semantic types, so the platform's detection on it skips the
  typed columns; the step still ends on `ddt-semantic-type-detected` for that frame, and makes
  row 0 current itself so the table view's 1 s timer does not repaint every viewer mid-step.
- **`expect` comes from `src/runtime/patience.js`** in every check: a `@known-failure` scenario
  narrows it to 3 s, and a check that names its own timeout wraps it in `pollMs`.
- **A throw inside an `expect.poll` callback ends the poll.** A read the viewer may be between
  layouts of returns `false` and keeps the reason for the failure message.
- **Every check is one sentence naming the alternatives** (`has no "x" area; it has: …`); a
  `StepFailure` carries the feature line, the step, the reason and what the page shows instead
  (`failure.ts`, `locate.ts` `explain`). Notes print only for `lint` / `--verbose`.
- **`@journey`**: one test, Background once, scenarios as soft steps that restore what they
  change, each owning its error and balloon floor; budget = the test timeout + 20 s per scenario.
- **Snapshots are lazy** (`viewer-runtime.ts` `snapshot`): a change keeps the bitmap, the areas,
  the readings, the range, the scale and the legend; the histogram and the per-area ink are
  computed on the first "than before" read. Every property set, menu pick, area gesture, resize
  and data step takes one; **no check moves it**, so several can follow one change.
- **A settle ends when the viewer says nothing is pending** (`isRenderPending`) and a render has
  landed — through every pass it announces; 10 s pending is a platform failure; a viewer without
  the signal falls back to a 300 ms cap. The resizer includes the first layout: otherwise a grid
  can report ready with its interactive overlay still 300×150. Settles are armed before the change. Negative checks
  (`not repainted`, `same range`, `same reading`) read after `quiet`.
- **A gesture aims where the viewer has finished putting the thing**: `hitArea(…, beforeChange)`
  settles first; `menuPoint` also waits for the anchor's box to hold for two frames; hit areas are
  polled for up to 5 s like elements. A context menu opens on `{element}` (a tree node too):
  `menuPoint` skips the viewer waits for a non-viewer.
- **Typed text is verified and retyped** (`typeVerified`); an editor that already has the focus
  is not clicked; a hit area typed into must end up owning the focus (`typeIntoArea`).
- **Platform keys are normalized in the shared gesture helpers.** `Control` / `Ctrl` becomes
  `ControlOrMeta`, including held modifiers for drags and legends; typing and clearing also
  select all with the platform modifier. `Delete` / `Del` follows `d4/shortcuts.dart`
  (Backspace on macOS, Delete elsewhere). Physical keys are `ControlLeft` / `ControlRight`,
  `ForwardDelete` and `Backspace`; the grid's custom current-cell copy requires `ControlLeft+Shift+C`.
- **Every Dart column picker goes through `pickInColumnGrid`**: the first letter is pressed on the
  selector, the name retyped until the box holds it, Enter pressed on the box, and a popup still
  open afterwards is the failure.
- **Baselines and resizes wait for a finished viewer**; a size a step asks for is held by a
  `MutationObserver` until `restores the size` or the feature ends.
- **An in-page step returns nothing it does not read** — a platform object serializes for
  seconds (Bio's SeqHelper, 10 s per call).
- **Roundtrips are the cost**: a step is locate + one action; `onViewer`/`evaluate` bundle the
  install check with the call; traces keep no DOM snapshots and no per-action screenshots.
- **Names are global and unique; platform names are reserved**; app names live on a
  `context()`. Context switching is compile-time tracked (`enters`) and runtime explicit
  (`enter(page, …)`); every runtime parse goes through `refOf(page, target)`.
- **Noun resolution**: whole registered phrase → ordinal → split at the first scope word (`of`
  names a part) → registered element, else kind by every matching suffix, longest first. Kinds
  try their `match` strategies in order; a scope that is not on the page is not waited for.
- **Gestures act on the visible match** (`locateActionable`); `enabled`/`disabled` read every
  match and prefer the visible ones; `visible`/`hidden` over several matches = any/none.
- **Labels are found first, items second** (`byLabel`, `:scope >` labels); menu items match
  their own label, not their children's.
- **The context panel renders the current object (`grok.shell.o`) and nothing else**
  (`property_panel.dart` on `onCurrentObjectChanged`; `grok.shell.windows.showContextPanel`
  shows it). The setter drops a change to the object already current, one within 2 s of a
  property edit and one while the object is frozen. An explicit click in a top-level grid releases
  the property-edit guard, so a cell clicked just after expanding a pane still becomes current:
  `the context panel is open` is the Given, `the context panel should show "X"` names the object
  before any pane is read, and a failure inside the panel reports the current object and the
  panes in the DOM but not shown (`explain`). **A context pane that counts its items is hidden
  while the count is 0** (`accordion.css`, `.grok-prop-panel .d4-accordion-pane[d4-info="0"]`),
  and the count arrives asynchronously: Activity on a space created a second ago is `present`,
  not `visible`, until the server has logged the creation.
- **Escape goes to the topmost dialog** (`press`): the dialog closes on a keydown inside its own
  root, and the focus is not reliably there (the grid's 1 s timer, a menu that just closed).
- **A gesture is dispatched once; the target is decided before it** (`pickMenuPath`): a click
  whose handler rebuilds a viewer synchronously (the tile viewer, 0.9 s idle for 1000 rows, past
  3 s under four workers) outlives a short cap with its work done, and a fallback click undoes
  the toggle. Never `click().catch(() => otherClick())`.
- **Hover is two pointer events and never sleeps**; it waits in-page for the element's own
  `mouseenter` and repeats the pair when a coalesced move swallowed it.
- **Step specificity**: more literal text wins, then fewer parameters; a tie is a compile error.
  Viewer steps take `{widget}` (a phrase ending in viewer/widget, `grid`, `filter panel`).
- **Playwright scopes inner selectors to the element**: a `labelSelector`, a part or a `has:`
  filter is evaluated from the outer element (`span:first-child`, not `.x > span:first-child`).
- **Headless and headed rasterize canvases differently**; the config launches with
  `--disable-accelerated-2d-canvas`, so a repaint check reads the same pixels either way.
- **Codegen emits names, never selectors**; `\n` endings, no timestamps; orphans removed on
  compile and reported by `--check`.

## Facts that cost a run each

- Dart names: viewers `viewer-<Type>`, inputs `input-host-<Caption>`, sections
  `div-section--<Name>`, on-canvas selectors `div-column-combobox-<bound property>` (hover-revealed,
  `X:AGE` with no space), menu items `div-Group---Item` with `aria-checked`, the close icon
  `name="Close"`, `camelCaseToCss` single-dashed. `Property.caption` is the raw name unless
  `@Prop(name:)` set one; two properties sharing a caption must be named.
- Menus: a Dart group opens on the first pointer move; the popup mirrors every property under a
  zero-size "Properties..." group, so labels occur twice — `openGroup` waits on the first visible
  candidate and tries every one; the top menu bar folds into a "more" group under 1920 px, its
  vertical groups are entered with two moves inside the item, Escape does not close it.
- Filters: `user filters rows where …` writes the filter bitset, and anything that calls
  `requestFilter` (a histogram on every menu pick) recomputes it — hold a filter across viewer
  interaction through a filter card. `getFiltersGroup` creates a panel when there is none.
- Readings: `rows shown` is `combinedFilter.trueCount` except where the viewer means "what the
  frame drew" (scatter plot, PC plot under a transformation, density plot); demog-1000 has 128
  blank HEIGHTs, so a plot on HEIGHT draws 872; its auto-picked category is DIS_POP; a calendar at
  1920×1080 shows about 12 weeks; a hierarchical node is `partially checked` when its children
  disagree. The grid reports colours lowercase; `current row` is 1-based, 0 for none.
- `RowSet.SelectedOrCurrent` is the selection when it has anything, else the current row;
  `MouseOverGroup` is `dataFrame.highlight`, which a hover elsewhere leaves standing.
- `painted in at least N colors` groups by hue: a linear scale is one colour. `repainted`
  measures the whole canvas; `the "x" area … should have repainted` one area. A hover highlight
  can be gone by the time a later step reads it: repaint checks right after the gesture.
- A JS viewer joins by `getWidgetStatus()` (canvas under `parts`, `hitAreas` in CSS px,
  `values`), `get isRenderPending()` and `onRendered`; a package viewer's surface reaches the
  stand only when the package is republished, and a library viewer only when the package's
  `node_modules/@datagrok-libraries/<lib>` points at the checkout.
- The u2 side: `ChoiceInput` is a native `<select>`; comboboxes open on a keystroke; a tree row
  click selects and the twistie toggles; popups are portaled under `.u2-overlay` with
  `data-u2-owner` = the nearest named ancestor; plain `button()`, toolbar buttons and tab headers
  carry no `data-u2`. A `funcForm` number field loses a typed leading "-" (unreported u2 bug).
- A u2 name is one token: the locator tries the phrase without spaces and with dashes.
- The Dart column picker ("Select columns...") is a grid viewer: `text of cell N of __name` names
  a row's column, `cell N of x` is its checkbox, its Search input filters without renumbering. A
  Dart property grid category (`tr.property-grid-category`) has no aria state, only its icon
  (`property-grid-icon-minus` open, `-plus` folded), which `readExpanded` reads. The 2 s drop of the
  Invariants (`AppEvents.propertyEdited`) reaches across features: a settings click on a new viewer
  right after another feature edited a property leaves the panel on the old one.
- A Dart choice input's phrase can resolve to its `<select>` itself; `select` handles both. The Share
  dialog of an entity that is not a project (a model) fetches the entity's project after it opens and
  its OK throws "Not initialized" before that: wait for the owner's grant row ("Full access").
- A viewer outside a table view (a function view's docked chart, a facet's small multiples) is
  reached through `DG.Widget.find(root)`; the function view's tabs are dock-spawn-ts handles in a
  shadow root (`.dockspan-tab-handle`, a CSS locator pierces it), and the viewers of its other
  tabs stay in the DOM with no rectangle — a claim names the tab it reads.
- A compute form's parameter switch is a Dart `SwitchInput` (`role="switch"`, `aria-checked`)
  that is not inside the input it governs: the sensitivity form puts it in the input's host, the
  fitting form before it as a sibling — `switchOf` looks in the element, then back over the
  siblings. Switching a parameter on replaces its input with a min and a max.
- The platform's script view is CodeMirror 5 (`.CodeMirror`), the packages' editors CodeMirror 6
  (`.cm-editor`); a document is not an input value, the text goes in at the caret. The Model Hub
  gallery is on the page and empty for seconds after `Compute2:modelCatalog` returns: wait for
  cards, not the element.

## Environment

- Stand: `http://localhost:8888` (nginx over datlas `:8082`), dev key `admin`; a dev stand's
  `pub serve` recompiles after a Dart edit and `global-setup.ts` fetches the client once before a
  browser does. When every login hangs, the hand-started Postgres container died: `docker start
  affectionate_einstein`.
- `grok-bdd run` defaults to 4 workers (`PLAYWRIGHT_WORKERS` overrides); the installed
  `@datagrok-libraries/test` base config is older than the checkout's and says one.
- A sharing feature shares with `DATAGROK_SHARING_LOGIN` or, unset, with the `bddsecond` user
  `global-setup.ts` creates through `POST /public/v1/users` with the dev-key token (a missing
  login answers 200 with an `ApiError` body; users cannot be deleted, so it stays).
- `pub serve` degrades under a run: the direct port (`:63343`) has served the bundle in 46 s
  while nginx at `:8888` answered from cache in 3 s; a run started while it is starved fails every
  feature at the shell load. `curl -o /dev/null -w '%{time_total}' localhost:63343/login.dart.js_1.part.js`
  under a few seconds first.
- Junctions, not `mklink /J`, from Git Bash (`New-Item -ItemType Junction`); `grok-bdd link` again
  after any `npm link`/`npm ci` in a package.
- Bash tool: cwd persists across calls, long heredocs fail — write files with the Write tool.
- Node 18 is what the libraries CI runs, and `@playwright/test` 1.62 exits at load below Node 20:
  the library pins `~1.61` (the last that runs on 18) until CI moves on; `npm test` there is the
  build, the drift check of the library's own project and the unit tests, with the locator tests
  skipping themselves where no Chromium is installed.

## Conventions

TS strict, 2-space, single quotes, `.js` suffixes on relative imports (NodeNext), comments only
for a non-obvious why. Every step is an exported const with a one-line `description` where the
phrase alone does not say what is checked. A step for one viewer that a second viewer could want
goes to `bindings/tiers/viewers/widgets.ts` and takes `{widget}`; a step that reads a status key
shape only that viewer publishes stays in the package. Nothing is committed without the lead's
order; stage by explicit path.
