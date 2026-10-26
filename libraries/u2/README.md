# @datagrok-libraries/u2

Next-gen Datagrok UI library: design tokens, a headless signal-based component core, and a
standalone component gallery. Compact, IDE-class visual language by default.

Architecture (see the design docs in the core repo, `core/docs/features/ui2/`):

- **L0 tokens** — CSS custom properties, single source of truth for both old CSS and u2.
- **L1 headless core** — TypeScript + signals; state, behavior, a11y, validation.
- **L2 spec & registry** — JSON-serializable UI, machine-readable component manifest.
- **L3 delivery** — fluent `u2.*` API, custom-element tags, generated framework wrappers.

L0–L2 are **platform-free**: their only `datagrok-api` import is `datagrok-api/u2core`, the shared
signal core. Platform integration lives in `src/dg`. Enforced by `tests/boundaries.test.js`.

`vendor/` holds frozen copies of permissively-licensed upstream code; see
[vendor/LICENSES.md](vendor/LICENSES.md).

## Layout

```
src/core/          signals, Scope, Component/Control, Input base, elements, async primitives
src/components/    platform-free controls, by role — the same vocabulary as the registry categories
  inputs/          everything with a value the form and the Dart bridge read and write (+ combobox, typeahead)
  forms/           form, property-grid — edit an object through rows of inputs
  containers/      splitter, tabs, accordion, dialog, wizard — arrange children, carry no value
  actions/         buttons, button-group, row/context actions
  navigation/      menu, menu-bar, toolbar, breadcrumbs
  collections/     virtualized list, tree, table
  display/         icon, badge, progress, async-view
src/spec/          dg-ui/1 document: renderer, registry, patch engine, binding paths
src/sources/       the data sources a spec can declare (platform-free, behind `backends`)
src/dg/            the platform layer — the only code that imports datagrok-api
  inputs/          platform-bound inputs (file, files, user, molecule, meta, rsa, column and table pickers) and the two Dart bridges
  forms/           object-form, property-editor, the Editors registry
  entities/        entity chips/cards, column rendering, dapi-backed sources
  shell/           app-view, widget-host, platformContext, the sources' backends, collection-selection, view-sync
  designer/        the visual designer
css/               one sheet per component, flat — installed once per page (src/styles.ts), never imported
```

Domain tables (EMS) are in [`@datagrok-libraries/ems`](../ems/README.md).

The root `src/index.ts` is the single public surface; import from it (or from `src/dg/index.ts`),
never from a component's file path.

## Install

Packages in this repo take it from the pnpm workspace:

```json
"dependencies": {"@datagrok-libraries/u2": "workspace:^"}
```

It is not on npm yet. `datagrok-api` is a peer dependency: it carries `datagrok-api/u2core`, the
signal core, which the platform serves once per page as `DG.U2` (every package build externalizes
it through `build-config`). Two entry points:

```ts
import {Form, TextInput, ChoiceInput, Dialog} from '@datagrok-libraries/u2';  // platform-free (L0–L2)
import {host, appView, ObjectForm} from '@datagrok-libraries/u2/src/dg/index.js';  // platform layer
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

What the first platform view built on u2 (Bio's monomer editor, 2026-10) got wrong before it got
right. Each rule names the API that makes it easy.

- **Build inside an owner.** Element factories and `each()` need an ambient scope: build in a
  control's constructor under `this.runInScope(() => ...)`. A component built there is disposed
  with the control.
- **Register on the thing that lives as long.** Inside an `each()` row or an effect run, put
  effects and cleanups on a control built in that row (`input.effect(...)`), not on the long-lived
  control around it — that one would gain one more per rebuild. With u2 debug on
  (`localStorage['u2.debug'] = '1'`, and in every bdd run with the leak audit) u2 warns when the same
  call site registers on a long-lived control from two different rows or runs.
- **Act on the user with `onInput`.** `onChanged` also hears writes from code (it matches the
  platform's `onChanged`). "The user picked a library, load it" is `onInput`.
- **Load a record with `reset(value)`; follow a model with `follow(read)`.** Neither fires
  `onChanged`; `reset` also clears the touched state and any unreadable text. `Input.system()` is
  the input's own pruning, and it is *not* silent.
- **Two verdicts.** `validity` judges the value from the start (what Save goes by — disable it and
  say why in its tooltip); `error` is what the input shows, only once it was edited, left or
  validated (what a tab marker or a status line goes by). A `Form` has both.
- **Tabs mark themselves.** Give a tab a `validity` signal, or make its content a `Form`. Before
  reaching for tabs in a narrow panel, consider showing what matters and folding the rest into a
  `Section` with a `hint` (what it holds), a `validity` (an error is never hidden by the fold) and a
  `persistKey` (the user's open or close is remembered).
- **Say why something is disabled.** `button(text, fn, {disabled: reason})` and
  `dialog.okDisabled(reason)` take the reason (a string or a signal of one) and show it as the
  tooltip; a dialog's OK also names the form's first failing field. A disabled u2 button is
  `aria-disabled` and keeps its tab stop, so keyboard focus shows the reason too; a press on it does
  nothing.
- **One label column across forms.** Forms that must line up — the tabs of one editor, a header
  over a body — share a `LabelColumn`; `labelWidth` fixes a column (one-letter captions).
- **Widths are options, not CSS.** `fill: true` stretches an editor across its row; `width`
  (px or a CSS length) sets it outright. No `!important` needed. Several inputs under one caption
  (Colors: Line, Background, Text) are `form.addGroup(caption, inputs)`.
- **Listen with `on(el, type, handler)`, follow with `effect(fn)`** on elements you build in a row
  or a card: both go with them. `this.scope.listen` or `this.effect` would tie them, and the element,
  to the long-lived control.
- **A docked editor panel.** `form.addColumns(inputs, {minWidth})` puts short fields two to a row,
  each editor filling its column, one to a row when narrow. `statusLine(problem, otherwise)` beside
  Save says what it does or why it cannot yet; a flush `Accordion` pane with a `hint` folds the rest; a
  chip that opens a `Menu` (an item's `detail` says what it stands for) changes a small choice in two
  clicks. The recipe: docs/recipes/forms-and-details.md. `new Form({quiet: true})` drops the underlines
  at rest: for an inspector, not a form people fill in.
- **Dock with `host(control)`.** It disposes the control when the platform closes its panel. A
  panel closed from code with `dockManager.close` skips that: check `root.isConnected`.
- **Find a view's control from its root.** The platform hands out a new wrapper object for a view
  each time (`grok.shell.views`), so nothing stored on one is there the next time; reach the content with
  `Control.forElement(view.root.querySelector('.my-view'))`.
- **Dialogs come with CANCEL.** The first `onOK` / `addButton` adds it, as in every platform
  dialog.

## Testing

- `DG.U2.Control.forElement(el)` reaches the live control from any element; `getState()` reads its
  state, `getItems()` its items, `data-u2-name` / `data-u2-part` / `data-u2-item` address it.
- Labels align in a microtask after a form is built, before the first paint; `form.layoutNow()`
  aligns at once for a test that measures in the same task. A resize re-aligns in the next frame.
- `Scope.liveCount` before and after a test proves it leaked nothing.
- BDD (`libraries/bdd`) drives u2 controls through generic phrases — a new screen usually needs
  no bindings.

## Development

```bash
pnpm install                 # from public/
npm run build                # styles, tsc to dist/, manifest and token checks
npm test                     # compiles in place and runs tests/*.test.js
npm run styles -- --bump     # after a CSS change: the next stylesheet revision, to commit with it
npm run manifest             # regenerates manifest.json from the registry
npm run gallery              # serves the standalone gallery at /gallery/
```

The reference is this README, [`docs/recipes/`](docs/recipes/README.md) and `manifest.json` (each
component's `usage`); `tests/docs-drift.test.js` fails when they name an API that does not exist.
`core/docs/features/ui2/` holds the design history — plans and reviews may name APIs that changed.
