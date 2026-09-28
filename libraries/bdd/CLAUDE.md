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
src/init.ts, src/cli.ts init | compile [--check] | lint | list-steps | run [playwright args] | guide <features> | link
src/runtime/            args, locate, gestures, assertions, harness (session, journey, resetShell, error floor),
                        viewer-runtime (in-page window.__bdd), viewers (readers over it), viewer-pixels,
                        viewer-menus, viewer-legend, menus (top menu), events, functions, patience, failure,
                        guide (BDD_GUIDE: per-step screenshots, located element, menu stops (hop), the page's own pointer events per stop → steps.json; full shell)
tool/guide-render.py    steps.json → guide.mp4 / step-NN.png / steps.md / audit.png (+ --gif: guide.gif, guide-thumb.png)
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
`d4-balloon-shown`, `grok.shell.autostartsCompleted`, `Resizer.isResizePending`, `data-legend-*`,
`DG.Widget.addStatusProvider` (a package's own areas and readings on a native widget, e.g. the
WebLogo glyphs Peptides draws in grid headers), the annotation regions' `region "<title>" title`
hit areas and `title strip top` / `title strip right` / `region titles shown` readings
(`AnnotationRegionsMixin.addStatus`).

## What never becomes a feature — hard rule

A feature tests what a user does and sees in the browser. Two kinds of test never go in one — not
when translating, not when looking for gaps, not when a TestTrack case asks for it (lead's ruling,
2026-09-23):

- **Anything that runs a server-side script or compute container, or reaches an outside web
  service**: Python/R/Julia scripts, Jupyter kernels, Docker containers, third-party lookups.
  Examples: Chem Curate, Mutate, IUPAC Name, Generate Conformers, Butina, Synthon Search, the 3D
  Structure and Gasteiger panes (Python); Descriptors and Map Identifiers (the chem-chem
  container); the Identifiers pane (UniChem and PubChem lookups); Bio Molecules to HELM and its 3D
  embedding. The outcome depends on the stand's kernel, containers and network (a cold kernel after
  a restart hangs past any budget; a late error from a lookup fails the next scenario's "no errors"),
  so the suite reports the environment, never the UI. Test the script with a package test. Check the
  code path, not the menu name: a JS-looking command or pane can call
  `grok.functions.call('<Pkg>:<PythonScript>')`; the package's `scripts/` folder lists the scripts.
  A pane builds when expanded, and the expanded state persists in `localStorage` for the worker's
  page, so a scenario that expands a server-backed pane makes it build in later scenarios too.
- **Anything with nothing UI-specific**: a function called with arguments and its result checked,
  a server outcome no UI shows. That is a package test (`src/tests/`) or an `ApiTests` test.

One exception, by the lead's ruling: the Scaffold Tree features stay. The viewer's tree comes from
the Python `GenerateScaffoldTree`, but what they test is the viewer's own UI (checking, colouring,
filtering, editing and removing nodes), which no package test reaches. Similar things should stay/be translated as well, as long as they actually test ui.

A TestTrack case marked `target_layer: manual-only` or `apitest` is never translated. In a
`playwright` case, a scenario of either kind is skipped, and the feature description says so in
one line. A gap hunt counts these as covered elsewhere, not as gaps.

## Everything a feature puts on the server goes — hard rule

A feature leaves the stand as it found it (lead's ruling, 2026-09-24). Whatever it adds or changes
on the server is removed or restored at feature end (`atFeatureEnd`), and swept again when the
feature starts, since a killed run never reached its end; the cleanup reads the server back to
prove it is gone. That covers:

- entities it creates: projects, tables, layouts, queries, scripts, connections, spaces, groups,
  files it uploads, rows or tables it writes into a database;
- what the UI makes on the side: the layout and project a query or script save writes, the chat a
  Chats pane post creates, the grant a Share dialog adds;
- changes to things the feature does not own: a connection's identifiers configuration, a catalog's
  comment, a shared connection's parameters, the account's settings — put back as they were.

A feature that cannot undo what it does (a user cannot be deleted) works on a fixed fixture it makes
once and reuses, never one per run. A change nothing can undo does not go in a feature.

## Invariants — what must not regress

- **One registry, through `dist/`.** Specs import the library by package subpath, project bindings
  by relative path; never mix `src/` and `dist/` in one run. ESM everywhere. One `@playwright/test`
  per run: a package of the pnpm workspace depends on `workspace:^` and resolves the library's copy
  with nothing to link; one outside it depends by path (`file:…`) and `grok-bdd link` makes its
  Playwright the library's copy (redo after `npm ci`). A package's `grok-bdd` loads `dist/`, so it
  stops while a source is newer than its build (`staleBuild` in cli.ts): an "unknown kind" or "no step
  matches" right after a pull is a stale build, and a kind is defined in `bindings/`, not `src/`.
- **One page per worker** (`harness.ts`): `feature(test)` reuses the worker's page, `afterEach`
  resets the shell (first waiting up to 60 s for the command the scenario armed and up to 25 s until the task bar has no progress entry — an
  analysis a scenario left running reopens its closed table and makes it current in the next feature;
  a menu command's `onAfterRunAction` can come before its work ends — then Escape for dialogs and
  menus, `ui.tooltip.hide`, notices removed, `closeAll`, Home current), `afterAll` runs all the feature's `atFeatureEnd` cleanups and fails if any fails. Never open several
  Datagrok pages in one browser. The renderer never gives back what a feature took (about 2 GB a
  minute of work, outside the JS heap), so a page older than `BDD_PAGE_MAX_MIN` (2) minutes is
  closed with its context after its feature and the next feature opens a new one (a new page of
  the same context shares the old renderer process). A journey
  scenario that fails closes its dialogs and menus (Escape) before the next scenario starts.
- **`user is logged in` only resets when the page is in the shell**; it sets `simpleMode` (view
  tabs hidden — switch views by name), clears the error and balloon floors, installs the in-page
  runtime. It waits for the PowerPack Home widgets to finish loading first, so what a widget logs
  (GROK-20891) stays out of the scenarios. Package autostarts land 3 s after boot; a feature that needs one awaits
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
  match and prefer the visible ones; `visible`/`hidden` over several matches = any/none. An
  ordinal counts the visible matches: the Home page keeps its own viewers in the DOM, hidden,
  under a table view, so "first line chart viewer" must not land on one of them.
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
- **An area hover settles the viewer after moving the pointer.** Grid cell tooltip requests use
  the tracked debounce, including in nested correlation grids. Negative tooltip text checks count
  visible matching tooltips; a hidden or absent tooltip has no displayed text.
- **Step specificity**: more literal text wins, then fewer parameters; a tie is a compile error.
  Viewer steps take `{widget}` (a phrase ending in viewer/widget, `grid`, `filter panel`).
- **Playwright scopes inner selectors to the element**: a `labelSelector`, a part or a `has:`
  filter is evaluated from the outer element (`span:first-child`, not `.x > span:first-child`).
- **Headless and headed rasterize canvases differently**; the config launches with
  `--disable-accelerated-2d-canvas`, so a repaint check reads the same pixels either way.
- **Codegen emits names, never selectors**; `\n` endings, no timestamps; orphans removed on
  compile and reported by `--check`.

- **Server fixture names may include `{run}` and `{time}`**: the compiler resolves strings, element
  phrases, tables and doc strings through the feature session. One UUID and one start time per
  feature instance keep workers and repeated runs independent; never generate either at compile time.
- **Spaces cleanup verifies IDs against every page of the root listing.** Spaces smart filters
  can return an empty list for an existing ID, so a filtered result cannot prove deletion. Match
  the captured IDs locally, delete by exact ID, and retain unrelated roots. Include a fixture's
  parent root in its cleanup names because the listing does not include child spaces.
  After setup cleanup, refresh an open Browse tree: API deletion leaves cached nodes behind,
  so recreating the same name otherwise targets a stale node or resolves to two nodes.
- **A killed run never reaches its feature-end cleanup**: a fixture named with `{run}` or `{time}`
  is also swept by family — the same name with any run suffix, older than an hour — whenever a
  `no … named` or `a … named` step runs, and by the project save (`isStaleFixture`,
  platform/steps.ts). Fixed names (the spaces, most projects) are swept by their exact name.
- **Groups and roles are cleaned like spaces** (the complete listing, never a name or ID filter —
  `grok.dapi.groups.filter('name = …')` missed a group that existed). Their global permissions are
  revoked before `grok.dapi.groups.delete`, which refuses a role that holds one (GROK-20904). Never
  delete a group as an entity: that leaves its grants behind with no group, and the Global
  Permissions pane of every role shows "error" (GROK-20901).
  **A group's chats go first**, found by the group's id (`/api/chats/with_groups?ids=`, as the
  client's Chat does) and deleted through `DELETE /api/chats/{id}` (the JS API has no chats): Chat on
  a group makes a private chat in a hidden group. Deleting that hidden group, or the group, first
  leaves a chat the server can no longer delete and that throws in every user profile's chat
  listing (`forum.dart` `_refreshChats`) for that account — it happened once on dev.
- **A translated stack trace is not a second error**: the platform logs "… Look below, ID = X" and,
  seconds later, "Stack trace X"; the floor joins it to its error, or drops it once reported.
- **Model cards are not completion signals.** The Train Model preview reports `aria-busy` before
  debounce/queued training and `aria-invalid` for unavailable or failed results. `model preview
  should be ready` requires the latest completed training, predictions, charts and history.
- **Nested viewers resolve by their own root.** A scatter plot inside a JS viewer must not resolve
  to the enclosing viewer merely because that viewer contains its element.


## Facts that cost a run each

- Dart names: viewers `viewer-<Type>`, inputs `input-host-<Caption>`, sections
  `div-section--<Name>`, on-canvas selectors `div-column-combobox-<bound property>` (hover-revealed,
  `X:AGE` with no space), menu items `div-Group---Item` with `aria-checked`, the close icon
  `name="Close"`, `camelCaseToCss` single-dashed. `Property.caption` is the raw name unless
  `@Prop(name:)` set one; two properties sharing a caption must be named. `annotate` prefixes the
  element's tag (`div-`, `span-`, `icon-`…) and turns `:_; *\[]{}|` into dashes — the `dart` match
  tries that form; a sketch box's name sits on its `.d4-host`, not on `.d4-sketch-item`.
- A context menu of a non-viewer element is opened at a point of its visible part that
  `elementFromPoint` gives back to it: a tree row runs past its panel's edge, and a panel docked a
  moment ago (the console) covers part of a gallery. `scrollIntoView` is the last resort — it moves
  overflow-hidden ancestors too, and the page with them.
- An entity saved through the JS API gets its first Activity entry late, and a counting pane hides
  at 0: claim its count (`should count at least`), never the pane's presence right after the save.
- A guide takes the pointer from the page, not from `page.mouse`: a locator's `click()`, `hover()`,
  `dragTo()` never pass through it, and a recorder that wrapped only the mouse had 60 of 80 clicks
  with no place (drawn unmarked, at the element's centre). `guide.ts` logs trusted pointer events
  in the page (capture listeners on every document) and drains them into the stop of the step they
  belong to; `page.mouse` is wrapped only to picture a drag.
- Package tools that decorate a view on `onViewAdded` (DevTools' Signature Editor icon) attach in
  the package's autostart; a view opened before it runs used to miss them (fixed in DevTools
  2026-09-24 by decorating the open views too).
- Menus: a Dart group opens on the first pointer move; the popup mirrors every property under a
  zero-size "Properties..." group, so labels occur twice — `openGroup` waits on the first visible
  candidate and tries every one; the top menu bar folds into a "more" group under 1920 px, its
  vertical groups are entered with two moves inside the item, Escape does not close it. The bar
  rebuilds when a package's entries arrive, with no signal that it is done: a pick can find its
  leaf and then click a collapsed group, so `pickTopMenu` makes the whole pick once more from the
  bar (a core "menu settled" signal would remove that).
- The Dart property grid's choice editor is lazy: its `<select>` enters the value cell only once
  that cell is clicked (`select` clicks it first). The Save project dialog's name field is a bare
  `<input>` (aria-label "Name"), which `text input` reaches. Typing into a column picker's search
  box used to toggle the scatter plot's regression line on every "r" (the R shortcut listened on
  the plot's root) — fixed 2026-09-21 in `regression_line.dart`; the box plot's T (p-value) and the
  line chart's R had the same hole until 2026-09-24 (picking "HEIGHT" hid the p-value). A
  single-key shortcut on a viewer's root ignores keys whose target is an input or a text area. A
  column selector named after a property with a space is reached without it (`"Category 1"` is
  `div-column-combobox-category1`).
- Filters: `user filters rows where …` writes the filter bitset, and anything that calls
  `requestFilter` (a histogram on every menu pick) recomputes it — hold a filter across viewer
  interaction through a filter card. `getFiltersGroup` creates a panel when there is none. A click
  on a categorical card's category name applies "only this one" from a 1 ms debounce that reads the
  card's current row when it fires (`grid_filter_base.dart`); Chrome runs queued input before timers,
  so a checkbox click on the same card right after it can move that row first and the card keeps
  the wrong category — claim the count after the name click before the next click on that card.
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
  `should show a/no selection highlight` reads the canvas as it is, so a viewer no gesture touched
  (a histogram that shows a selection made in the WebLogo) needs no snapshot; `more`/`less …
  than before` compare with the snapshot.
- An area phrase can name a part of a hit area: `left edge of x axis` (a 16 px strip along that
  edge — the axis away from the column selector in its middle), `top left corner of region Older`
  (a 16 px square), `overlap of region Tall and region Heavy` (the rectangle two areas share), and
  they nest (`left edge of overlap of …`). Resolved in-page by `edgeOf`, through `findArea` and
  `quietAreaRects`, so every gesture step, `should have a … area`, the menu steps and the checks
  that compare areas with each other or with a remembered place take them; the ink readings and
  `taller/wider than before` (`areaRectChange`) do not.
- A marker under the pointer takes precedence over an annotation region: the scatter plot hit-tests
  its regions only while no marker is hovered (`navigation.dart`), so a region gesture aims at a
  marker-free point (half-integer AGE on demog, small `markerDefaultSize`) and a title click is
  preceded by a hover. A right-click over a region opens the region's own menu (Edit… | Show
  Annotation Regions | Title Font), not the viewer's; `mouseOverRegions` is not cleared by hiding
  the regions, only by the pointer leaving the viewer. The histogram hit-tests a region only where
  no bin is under the pointer and its bins slider sits at the top centre of the plot; the scatter
  plot's Color and Size selectors cover the top-right corner; the bar chart hit-tests between the
  bars, and a band on its aggregated axis selects the rows of the bars whose value lies in the band,
  not rows by their own value. A line chart binds a region to its aggregated caption
  (`y: "avg(WEIGHT)"`).
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
  Invariants (`AppEvents.propertyEdited`) reaches across features: an implicit current-object change
  on a new viewer right after another feature edited a property leaves the panel on the old one;
  a viewer's "Properties..." command forces the change (since 2026-09-23).
- Users, groups, roles: a login takes `[a-z0-9._-]` only (`grok_user.dart` `validateLogin`). A user
  cannot be deleted: a feature takes the `bddviewed` fixture user to look at, or `bddmanaged` to
  join, disable and favorite (the `@serial` features, never at the same time), both made once per
  stand as `bddsecond` is; only users-create adds users, the two it tests. A group
  saved through the JS API without a friendly name is listed by `camelCaseToWords(name)`
  ("BDD-probe" shows as "BD D-probe"), so `a group named` sets both. A role is a group the JS API
  cannot flag, so a feature makes one in the New Role dialog. The users search is fuzzy (a login
  brings up every login sharing its letters), and a gallery counter reads `N`, `N of M` (M is the
  list, only N rendered) or `shown / total`, with `...` before it knows — an item outside a search
  is no claim, the counter is. A view-mode icon says it is current with `d4-current`, which
  `selected` reads inside the gallery toolbar only (elsewhere it marks the current card); the gallery itself carries `mode="Brief|Card|Grid"`.
  Sort list keeps the chosen order in the browser (`OrderMenu.savePersonal`), so a feature that sorts
  ends on Default; "Apply for all users" writes the order for everyone — never a step. While the
  gallery reloads after an order pick the counter reads "..." and the first item is empty, which is
  not "another item". A user's personal group (friendly name = login) is not listed by the Groups
  view; `user.tag()` throws on a User, so no feature can tag one.
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

## Claims that cannot fail — what every audit finds again

Each of these passed green while the thing it named was broken (audits of 2026-09-21 and
2026-09-23). Read every Then of a new feature against this list before running it.

- **A page function sees only its own source.** A Node-side helper named inside `page.evaluate` /
  `waitForFunction` is a `ReferenceError` in the page, and a `.catch` around the call turns the
  check into a no-op: the shell reset's task-bar wait never waited from the day it was written.
  Pass the function itself, or its source as text (`` `(${fn})(…)` ``). Under `tsx` (the unit
  tests, a package's bindings) a named function inside one — `const add = (e) => …` — is wrapped
  in a `__name` call the page does not have: keep inner functions anonymous, passed inline.
- **Absence is not a value.** A throw ends an `expect.poll`, so a reading or a column a computation
  adds late fails the claim at once — the reading steps return `MissingReading` and the column
  polls the missing column's text instead; neither may ever satisfy a negative.
- **A negative, a zero or "unchanged" right after an async gesture reads the state before it.**
  Two search types in a row that both keep 0 rows, a header that echoes the property just set while
  the cards re-render 200 ms later, "no new column" read once while the analysis runs: pair every
  such claim with one the gesture must change (a remembered reading that must differ, an end signal
  first — a balloon, a column, `… should have finished updating`).
- **`the top menu command should have completed` is the menu function's call.** A package command
  whose function only builds and shows its own dialog has ended before OK; the step now fails after
  that OK (`waitCommand`), and the claim is what the OK produces.
- **Coarse where the number is known.** `fewer than N`, `at least 1`, `lower than before` all pass
  on 0; "split differently" on label text passes for the same partition renumbered. Write the exact
  count the description states — measured in a run, never guessed.
- **Echoes.** A property, tag or header read back right after the step wrote it, a setting-derived
  count (`regions shown`, `formula lines` = active, not drawn): claim what the frame drew (its hit
  area) or what the change caused downstream.
- **Name resolution lands on journey leftovers.** A second `mutations` table, an `R1` column from
  the first run: close what a scenario made, or claim the new name (`R1_1`).
- **Text and pixel matchers are wide.** `contain text` is case-insensitive textContent over the
  whole element, hidden children included; a list reading's `contain` is membership after a comma
  split (an item with a comma is refused); white is the blank the colour steps skip (refused), grey
  gridlines are ink (`painted` skips a 2 px border) — a cell claim needs a hue or a reading.
- **Fixtures that cannot tell the semantics apart**: Shift-adding a region nested in the one already
  selected gives the same rows as a replace; pick a pair whose union differs from either.
- **State the harness or an earlier feature left**, stated as the product's: `openTable`'s current
  cell, the context panel or toolbox open, a custom-event count (reset by `listens for`), a per-account
  setting toggled through the UI (pin it with a Given that restores it at feature end).
- **Package bindings take `expect` and `pollMs` from `@datagrok-libraries/bdd/runtime`**, never from
  `@playwright/test`: the `@known-failure` narrowing and `BDD_EXPECT_TIMEOUT` go through them.
- **Titles and descriptions that promise more than the Thens claim** ("…and dropping the tree
  removes it", "builds a ball-and-stick view" checked as "shows no error") — trim the text or add
  the claim.

## Environment

- Stand: `http://localhost:8888` (nginx over datlas `:8082`), dev key `admin`; a dev stand's
  `pub serve` recompiles after a Dart edit and `global-setup.ts` fetches the client once before a
  browser does. When every login hangs, the hand-started Postgres container died: `docker start
  affectionate_einstein`.
- `grok-bdd run` defaults to 4 workers (`PLAYWRIGHT_WORKERS` overrides); the installed
  `@datagrok-libraries/test` base config is older than the checkout's and says one.
- Budgets a slow stand may raise: `BDD_EXPECT_TIMEOUT` (every check, 15 s), `BDD_COMMAND_TIMEOUT`
  (the columns a top-menu command adds, 180 s). `BDD_FRESH_PAGE` reloads the shell before every
  feature instead of resetting it: a feature that fails only after another one, and passes with
  the variable set, is failing on state the other left behind.
- Every run leaves its JSON report in the project's `test-results/report.json`, with the stand and
  the machine in `config.metadata` (a `--reporter` on the command line gets `json` added). The run
  history in `packages/UsageAnalysis/bdd/history` records such reports **only when the user asks
  for it** — never as part of a run, a review or a fix (`node history.mjs record --note …`, then
  `html`; see that package's bdd README).
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
