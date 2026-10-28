# @datagrok-libraries/u2

Next-gen Datagrok UI library: design tokens, a headless signal-based component core, and a
standalone component gallery. Compact, IDE-class visual language by default.

Architecture (see the design docs in the core repo, `core/docs/features/ui2/`):

- **L0 tokens** — CSS custom properties, single source of truth for both old CSS and u2.
- **L1 headless core** — TypeScript + signals; state, behavior, a11y, validation.
- **L2 spec & registry** — JSON-serializable UI, machine-readable component manifest.
- **Platform layer** (`src/dg`) — the bridge onto Datagrok's views, inputs and events.

L0–L2 are **platform-free**: their only `datagrok-api` import is `datagrok-api/u2core`, the shared
signal core. Platform integration lives in `src/dg`. Enforced by `tests/boundaries.test.js`.

`vendor/` holds frozen copies of permissively-licensed upstream code; see
[vendor/LICENSES.md](vendor/LICENSES.md).

## Layout

```
src/core/          signals, Scope, Component/Control, Input base, elements, async primitives
src/components/    platform-free controls, by role — the same vocabulary as the registry categories
  inputs/          everything with a value the form and the Dart bridge read and write (+ combobox, typeahead,
                   choice chip)
  forms/           form, property-grid — edit an object through rows of inputs
  containers/      page, splitter, tabs, accordion, dialog, wizard — arrange children, carry no value
  actions/         buttons, button-group, row/context actions
  navigation/      menu, menu-bar, toolbar, breadcrumbs
  collections/     virtualized list, tree, table, catalog table, master-detail
  display/         icon, badge, progress, async-view, app-header
src/spec/          dg-ui/1 document: renderer, registry, patch engine, binding paths
src/sources/       the data sources a spec can declare (platform-free, behind `backends`)
src/dg/            the platform layer — the only code that imports datagrok-api
  inputs/          platform-bound inputs (file, files, user, molecule, meta, rsa, column and table pickers) and the two Dart bridges
  forms/           object-form, property-editor, the Editors registry
  entities/        entity chips/cards, column rendering, dapi-backed sources
  shell/           view-host (ViewHost: a platform view as an owner), widget-host, platformContext, the sources'
                   backends, collection-selection, view-sync, table-view-chrome
  designer/        the visual designer
css/               one sheet per component, flat — installed once per page (src/styles.ts), never imported
```

Domain tables (EMS) are in [`@datagrok-libraries/ems`](../ems/README.md).

Import from the entry points only, never from a path inside the package (`…/u2/src/…`):
`@datagrok-libraries/u2` (platform-free), `@datagrok-libraries/u2/dg` (the platform layer),
`@datagrok-libraries/u2/filter` and `@datagrok-libraries/u2/dg/filter` (the filter model and its platform side), and
`@datagrok-libraries/u2/testing` (the testing kit).

## Install

Packages in this repo take it from the pnpm workspace:

```json
"dependencies": {"@datagrok-libraries/u2": "workspace:^"}
```

It is not on npm yet. `datagrok-api` is a peer dependency: it carries `datagrok-api/u2core`, the
signal core, which the platform serves once per page as `DG.U2` (every package build externalizes
it through `build-config`). The two entry points most code takes:

```ts
import {Form, TextInput, ChoiceInput, Dialog} from '@datagrok-libraries/u2';  // platform-free (L0–L2)
import {ViewHost, ObjectForm} from '@datagrok-libraries/u2/dg';  // platform layer
```

**No CSS to import.** The first u2 control a page loads installs u2's stylesheet once, as a single
`<style data-u2-styles data-revision="N">`, placed after the platform's stylesheets and before every
package's, so a package's own CSS wins a tie with u2's. Every package carries a copy; a newer revision replaces an
older one in place, so the page always runs the newest u2 CSS whatever order packages load in.
Changing a sheet in `css/` needs `npm run styles -- --bump`, which mints the next revision in
`css/styles.lock.json` — commit it with the CSS (a build keeps the committed revision; CI fails on a
stale lock).
Only the opt-in sheets are importable: `css/icons.css` (Font Awesome for hosts outside the
platform) and `css/theme-dark.css`. A standalone page that links the sheets itself opts out with
`<meta name="u2-styles" content="linked">` (the galleries do).

## Rules that matter

The rules u2 code keeps. Each names the API that makes it easy and what holds it: a lint rule `check-consumer`
runs, a usage error the engine raises, a finding of the operability probe, a type error or a ratchet count
(`consumers-baseline.json`). [AGENTS.md](AGENTS.md) has the same rules, each as a line done right and a line done
wrong; both lists are generated from `src/core/rules.ts` (`npm run manifest`).

<!-- rules: generated from src/core/rules.ts by `npm run manifest`; edit them there -->
- **Import from the entry points only.** `@datagrok-libraries/u2` (platform-free), `@datagrok-libraries/u2/dg` (the
  platform layer), `/filter` and `/dg/filter` (the filter model and its platform side) and `/testing` (the testing
  kit); never a path inside the package (`…/u2/src/…`, `…/u2/dist/…`). There is no CSS to import: the first u2
  control a page loads installs u2's stylesheet. Held by: lint `u2/import-entry`; ratchet `deepImports`,
  `deepDgIndex`.
- **Build inside an owner, as a function.** Element factories and `each()` need an owner building. A piece of UI with
  a life of its own (shown and closed, rebuilt, in a dialog, addressed by tests) is a function component,
  `component((self) => { …; return {root, …api}; }, {name})`: its builder owns everything it makes, and all of it is
  disposed with the control; `self` has the owner verbs (`self.effect`, `self.subscribe`, `self.setInterval`,
  `self.task`). A piece that is only part of its parent is a plain function that builds it, called inside the
  parent's builder — it is the parent's. After an `await` nothing is building:
  `self.task(async (t) => … t.build(…))`. A tree that nothing else owns (an app's home that outlives its views)
  starts on purpose as `Scope.root('<what>')`, and whoever starts it disposes it. Consumer code does not subclass
  u2's controls. A control's own members (`name`, `root`, `scope`, `visible`, `effect`, `own`, `on`, `dispose`, …)
  are not names for a part: name it for what it holds (`nameInput`); u2 reports the collision (`reserved-key`). The
  first sample: docs/samples/first-control.ts. Held by: lint `u2/build-in-owner`, `u2/no-ui-in-class`,
  `u2/no-subclass`, `u2/no-unowned-root`; usage error `owner-missing`, `owner-disposed`, `reserved-key`; ratchet
  `subclassU2`, `classBuildsUi`, `newScope`, `runWithUndefined`.
- **An input's value is `input.value`.** Read inside an effect or a `computed` it is followed, and assigning it
  writes it; `input.peek()` reads it untracked, and `input.signal` is the signal for wiring (`bind: input.signal`,
  `span(input.signal)` — `span(input.value)` shows a snapshot). A validator may read other inputs: it re-judges
  whenever they change. Held by: lint `u2/no-snapshot-in-live-slot`; API convention `value-signal`; ratchet
  `valueValue`, `valuePeek`.
- **Register on the thing that lives as long.** Inside an `each()` row or an effect run, put effects and cleanups on
  a control built in that row (`input.effect(...)`), not on the long-lived control around it — that one would gain
  one more per rebuild. With u2 debug on (`localStorage['u2.debug'] = '1'`, and in every bdd run with the leak audit)
  u2 warns when the same call site registers on a long-lived control from two different rows or runs. Held by: debug
  warning (`tests/rebuilt-scope-warning.test.js`).
- **Act on the user with `onInput`.** `onChanged` also hears writes from code (it matches the platform's
  `onChanged`). "The user picked a library, load it" is `onInput`.
- **Load a record with `reset(value)`; follow a model with `follow(read)`.** Neither fires `onChanged`; `reset` also
  clears the touched state and any unreadable text. `Input.system()` is the input's own pruning, and it is *not*
  silent.
- **Save goes by `validity`; what a field shows goes by `error`.** An input has two verdicts. `validity` judges the
  value from the start (what Save goes by — disable it and say why in its tooltip); `error` is what the input shows,
  only once it was edited, left or validated (what a tab marker or a status line goes by). A `Form` has both.
- **A tab marks itself: give it a `validity`, or make its content a `Form`.** Before reaching for tabs in a narrow
  panel, consider showing what matters and folding the rest into a `Section` with a `hint` (what it holds), a
  `validity` (an error is never hidden by the fold) and a `persistKey` (the user's open or close is remembered).
- **Say why something is disabled.** Every off takes the reason the user reads, a string or a signal of one (`null`
  while it is on): `disabled: 'Pick a file first'` on an input, a button, an action, a menu item, a ribbon button;
  `dialog.okDisabled(reason)`. There is no `enabled` and no `disabled: true`. The reason is the tooltip; a dialog's
  OK also names the form's first failing field. A disabled u2 button is `aria-disabled` and keeps its tab stop, so
  keyboard focus shows the reason too; a press on it does nothing. Held by: type `ButtonOptions.disabled`,
  `InputOptions.disabled`; API convention `no-enabled`, `disabled-reason`; lint `u2/reason-text`; operability probe
  `disabled-without-reason`; ratchet `enabledOpt`, `okEnabled`.
- **An action's `run` returns a promise or nothing.** `run: () => save()` hands back what `save` answers; a value
  handed back by mistake (a `DG.Row`) is a type error, so wrap a call that answers one in braces. A `run` or a click
  handler that throws or rejects goes to its owner's `onError`, then up its owners, then to `Scope.onUnhandledError`,
  which logs it to the console: the user sees nothing. Say what failed where the work is (`notify.error(…)`), or once
  for a whole component (`self.onError = (e) => notify.error(…)`); `confirmDialog` and a dialog's `onOK` keep their
  own failure in the dialog. Held by: type `Action.run`.
- **Forms that must line up share one `LabelColumn`.** The tabs of one editor, a header over a body: each `Form`
  takes the same `labelColumn`; `labelWidth` fixes a column (one-letter captions).
- **Widths are options, not CSS.** `fill: true` stretches an editor across its row; `width` (px or a CSS length) sets
  it outright. No `!important` needed. A form is 700px at most: `new Form({maxWidth: 820})` sets its own. Several
  inputs under one caption (Colors: Line, Background, Text) are `form.addGroup(caption, inputs)`. A form's button row
  is on the right (`form.addButtons(build, {align: 'start' | 'editors'})` moves it). Spacing too:
  `divH([cancel, save], {gap: 'm'})` on the token scale (`xs` 2px … `xxl` 24px), never `style.gap`. Held by: ratchet
  `cssImportant`, `cssInternalLines`.
- **Lay a view out with `Page`; a list beside its details is a `MasterDetail`.**
  `new Page({header, sections, footer, maxWidth})`: the header stays put, the body alone scrolls and the action row
  stays pinned under it, its buttons on the right and a `statusLine` on their left; the forms on the page share one
  label column. A section is a heading with its actions over its content (`{title, actions, content, maxWidth}`).
  `new MasterDetail({items, keyOf, row, detail, disabled})` puts the rows — a check, a hint — beside the current
  one's detail, built when first shown and kept, off with its reason; the arrow keys move the current row and Space
  checks it. Spacing and heights are the `--u2-page-*` and `--u2-master-detail-height` properties, set on an element
  of your own. Held by: ratchet `cssLines`, `cssInternalLines`.
- **Style your own classes, never u2's.** A package's CSS styles its own classes and sets u2's `--u2-*` custom
  properties on its own elements — the ones each kind lists as `styling` in manifest.json (README › Styling); it
  never targets `.u2-*` classes or `[data-u2…]` attributes, and code never looks u2's parts up by class name. Tests
  reach a control through the testing kit (`byName`, `byKind`, `state`). Held by: lint `u2/no-internal-selectors`;
  ratchet `tsInternalSelectors`, `cssInternalLines`.
- **A value shown as an element is a `DynamicInput` row; the rail takes buttons.** Inside a form, a value drawn as an
  element (a structure drawing, a chart) is a `DynamicInput`: a row of the form, in its label column, that grows to
  the element. The rail beside an editor takes buttons and links (`addOptions(iconButton(…))`). Read-only details in
  a pane of their own are a `tableFromMap`: muted keys, wrapping values (docs/recipes/forms-and-details.md). Held by:
  type `Input.addOptions`; operability probe `oversized-option`.
- **Listen with `on(el, type, handler)`, follow with `effect(fn)`, on the elements a row builds.** Both go with the
  row or the card that built the element. `self.effect` on the component around the row would tie them, and the
  element, to the long-lived control. Held by: choke point `handlers`; debug warning
  (`tests/rebuilt-scope-warning.test.js`).
- **A docked editor panel is laid out by its form: columns, a status line and folds.**
  `form.addColumns(inputs, {minWidth})` puts short fields two to a row, each editor filling its column, one to a row
  when narrow. `statusLine(problem, otherwise)` beside Save says what it does or why it cannot yet; a flush
  `Accordion` pane with a `hint` folds the rest; a `ChoiceChip` (each item's `description` says what it stands for)
  changes a small choice in two clicks. The recipe: docs/recipes/forms-and-details.md. `new Form({quiet: true})`
  drops the underlines at rest: for an inspector, not a form people fill in.
- **Dock with `ViewHost.of(view).dock(control, place)`.** Closing the panel — its ✕, or `dockManager.close` from code
  — disposes the control, and so does the view's close; disposing the control closes the panel. Chrome that is not a
  view (a context-panel widget, a pane) takes `host(control)`, which the platform's close of that chrome disposes.
- **Find a view's control from its root.** The platform hands out a new wrapper object for a view each time
  (`grok.shell.views`), so nothing stored on one is there the next time; reach the content with
  `Control.forElement(view.root.querySelector('.my-view'))`.
- **A dialog comes with CANCEL; one whose work applies as it is done is `closeOnly`.** The first `onOK` / `addButton`
  adds CANCEL, as in every platform dialog. A dialog whose work applies as it is done (a tracker, a viewer of
  progress) is `new Dialog(title, {closeOnly: true})`: one CLOSE, never a CANCEL; ✕ and Escape close it the same way.
- **Personal settings apply at once; shared and admin settings stage behind SAVE N CHANGES.** A form of the user's
  own settings (`grok.userSettings`, a view's look) writes and keeps each valid edit as it is made, with no Save
  button: `propertyForm(props, settings, {onChanged: save})` (docs/samples/property-form.ts). Settings shared with a
  team or set by an admin keep the edits in the form and save them together with one button that counts them ("SAVE 3
  CHANGES"), off with its reason while there is nothing to save or a field is invalid.
- **A dialog opened from a row outlives the row.** Shown, a dialog or a `confirmDialog` leaves the row or effect run
  that built it, so a rebuild does not close it; the content's own handlers and live text stay with whoever built
  them, so build that content with `dialog.add(() => …)`. A `keepAlive` dialog belongs to its builder: make it in the
  control that reopens it, not in a row.
- **A balloon's content is the balloon's.** `notify.info(() => div(['Deleted. ', link('Undo', undo)]), {onClose})`
  builds the content owned by the balloon and ends it when the balloon goes; an element built in a row's handler
  would keep its handlers on the row.
- **An app home starts with `AppHeader`.** The icon, what the app does as bullets (text, an element or a link each)
  and Learn more: the same header in every app, built in one call (docs/recipes/shell-integration.md). It follows its
  own width, so build it under the view's host — in the `ViewHost.create` builder, or in `ViewHost.of(view).build(…)`
  on a view you made — and it goes with the view.
- **Own a platform subscription with the owner's `subscribe`.** `self.subscribe(grok.events.onViewRemoved, handler)`
  — or `host.subscribe(…)` on a view's host — ends with that owner, runs the handler batched and contained, and shows
  in the leak audit. The free `subscribe` from u2 takes whoever is building: an effect run or an `each()` row, so a
  subscription made there goes when the run or the row does. After an `await` nothing is building: subscribe on the
  owner, and start timers on it (`self.setInterval`). Held by: lint `u2/no-raw-subscribe`; usage error
  `source-missing`; API convention `owner-positional`; ratchet `rawSubscribe`, `freeSubscribe`.
- **Code living in a platform view is owned by its `ViewHost`.** `ViewHost.of(tv)` is one host per view, whichever
  wrapper you hold, disposed when the shell removes the view: build under it (`host.build(…)`; after an `await`,
  `host.task(async (t) => … t.build(…))`), subscribe with `host.subscribe(stream, handler)`, read `host.signals` (the
  current row, selection, columns and the `grid` as signals; an effect over `grid` re-subscribes per grid, no
  timers), put the ribbon group up with `host.ribbon({actions, buttons})`, show its dialogs with `host.show(dialog)`
  and dock with `host.dock`. A new view is `ViewHost.create({name, content})` (docs/recipes/shell-integration.md).
  Held by: usage error `view-not-ready`; ratchet `viewZoo`.
- **Ask with `confirmDialog()` before what cannot be undone.** CANCEL and the action as a verb; `onConfirm` keeps the
  dialog busy until the work is done.
- **Rows the user edits are a `RecordsInput`; a pick with an order is a `ChecklistInput`.** A template's fields or
  stages, the functions to run in order — never an `ItemsGrid` or a hidden columns input
  (docs/recipes/forms-and-details.md).
- **A list of records is a `CatalogTable`.** Columns (`value` to sort and group by, `render` to draw), groupings, row
  `actions` and a `persist` adapter; it owns the order, the groups and the chosen columns, and `chrome()` gives the
  heading its columns, group and refresh icons. Never a `ui.table` with header listeners and injected rows
  (docs/recipes/catalog-table.md).
- **A column choice is `columnInput`.** Picking a column of a table is u2/dg's `columnInput(label, df, …)` (or
  `Platform.listedColumns(df, filter)` for a list of your own): they list what the platform's pickers list, hidden
  `~` columns left out. Never a `ChoiceInput` over `df.columns`. Held by: lint `u2/no-raw-column-list`; choke point
  `column listing`.
- **A table choice is `tableInput`.** `tableInput(label)` (u2/dg) offers the open tables the platform lists and holds
  one of them; `Platform.tableNames()` is that list for a picker of your own. A table with no name is never offered:
  a choice list holding one reports `item-null`, and a table input given a table it does not offer reports
  `table-not-offered`. Held by: usage error `item-null`, `table-not-offered`; choke point `table names`.
- **One text for a value: `Format`.** `Format.value(v)` is the text u2 shows for a value — a date in the user's
  platform date format, read in the browser's zone; `Format.date(d, {format})` takes a platform pattern; a Date on
  midnight UTC that stands for a day reads in UTC with `{utcDates: true}` (docs/recipes/forms-and-details.md). Held
  by: choke point `formatter`.
- **Comment only what the code cannot say.** Doc blocks go on the API a caller reads; no file-header narratives and
  no doc blocks on private or protected members; a one-line `//` says why where the code cannot. Held by: lint
  `u2/comments`; ratchet `fileHeaders`, `privateDocs`.
<!-- /rules -->

## Styling

A package styles its own classes and tunes u2 with the custom properties each kind lists as `styling` in
`manifest.json`, set on an element of its own — never with a rule over `.u2-*` classes or `[data-u2…]` attributes:

```css
.my-home { --u2-page-max-width: 960px; --u2-page-section-indent: 14px; }
.my-compute { --u2-master-detail-height: 320px; --u2-row-selected-bg: var(--dg-selection-bg); }
```

| Property | Kind | What it sets |
|---|---|---|
| `--u2-page-max-width` | page | the body's widest; its forms go as wide |
| `--u2-page-padding-x` | page | the space left and right of the header, the body and the action row (24px) |
| `--u2-page-padding-y` | page | the space above the header and the body (12px) |
| `--u2-page-gap` | page | the space between the sections (24px) |
| `--u2-page-section-indent` | page | how far a section's content sits in from its heading (0) |
| `--u2-form-max-width` | form | the widest a form gets (700px; as wide as the page in a `Page`) |
| `--u2-form-label-width` | form | the label column, which the form sets as it measures (`labelWidth` fixes it) |
| `--u2-form-section-gap` | form | the space between two sections of the form (0) |
| `--u2-master-detail-height` | master-detail | its height; it fills its host by default |
| `--u2-row-base-bg`, `--u2-row-hover-bg`, `--u2-row-hover-stripe`, `--u2-row-selected-bg`, `--u2-section-header-bg` | table, catalog-table (`--u2-row-selected-bg`: master-detail's current row) | a row at rest, under the pointer, selected; a group's header |

## Testing

`@datagrok-libraries/u2/testing` is the kit for unit and package tests, the same in `node:test` and in
`@datagrok-libraries/test`:

```ts
import {byName, expectOperable, mount, u2test, user} from '@datagrok-libraries/u2/testing';
import {TextInput, iconButton} from '@datagrok-libraries/u2';

test('the pencil turns editing on', u2test(() => {
  const name = mount(() => new TextInput({label: 'Name', value: 'KHL-3', readOnly: true}));
  name.addOptions(iconButton('pencil-alt', () => name.readOnly = false, {tooltip: 'Edit Name'}));
  expectOperable(name);
  user.click(byName(name, 'Edit Name'));
  user.type(name, 'KHL-4', {replace: true});
}));
```

- `u2test(body)` runs the body under a scope of the test's own and fails the test when it leaves a scope, a timer
  or a global listener live, when u2 reports a usage error meanwhile, or when an error reaches
  `Scope.onUnhandledError`. `mount(build)` builds under it and attaches; both go when the test ends.
- A usage error throws (`U2UsageError`) in u2's node runner, in bdd, in `grok test` and in dev (u2 debug on); in
  production it is reported once per call site. `U2_STRICT=report` (Node, bdd, `grok test`) or
  `localStorage['u2.strict'] = 'report'` lists them instead.
- `user.click`, `type`, `press` and `pick` act as a user does, and refuse what a user could not do — a covered or
  inert element, a disabled one (saying its reason), one that takes no pointer. A row action revealed on hover is
  reached with the keyboard (`user.press(button, 'Enter')`): no script can hover.
- `expectOperable(el)` fails, listing each finding, when anything offered under `el` cannot be operated as offered.
- `byName(within, name)` finds an element by its automation name or the name a user reads; `byKind(within, kind)`
  the controls of a kind; `state(el)` the owning control's state (`DG.U2.state`).
- `DG.U2.Control.forElement(el)` reaches the live control from any element; `getState()` reads its
  state, `getItems()` its items, `data-u2-name` / `data-u2-part` / `data-u2-item` address it.
- Labels align in a microtask after a form is built, before the first paint; `form.layoutNow()`
  aligns at once for a test that measures in the same task. A resize re-aligns in the next frame.
- BDD (`libraries/bdd`) drives u2 controls through generic phrases — a new screen usually needs
  no bindings.

## Development

```bash
pnpm install                 # from public/
npm run build                # styles, tsc to dist/, manifest and token checks
npm test                     # compiles to .test-build/ (never into src/) and runs tests/*.test.js
npm run styles -- --bump     # after a CSS change: the next stylesheet revision, to commit with it
npm run manifest             # regenerates manifest.json and the rules lists here and in AGENTS.md (from dist)
npm run manifest -- --from .test-build   # the same from what npm test compiled
npm run gallery              # serves the standalone gallery over that build and prints its URL
```

The reference is this README, [`docs/recipes/`](docs/recipes/README.md), [`docs/samples/`](docs/samples/README.md)
(the "Build with u2" samples, whole code) and `manifest.json` (each component's `usage`, each kind's `whenToUse` and
`notFor`); [`AGENTS.md`](AGENTS.md) is the digest for coding agents.
`tests/docs-drift.test.js` fails when they name an API that does not exist; `tests/recipes.test.js` compiles every
TypeScript block of this README and the recipes and every sample against u2's sources, and runs the live ones
under the testing kit.
`core/docs/features/ui2/` holds the design history — plans and reviews may name APIs that changed.
