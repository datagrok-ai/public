# Recipe: shell integration

A u2 app runs inside the Datagrok shell, which already provides per-view chrome. Use it —
don't rebuild it inside the view content.

| App part | Where it goes | How |
|---|---|---|
| Main view controls (filter box, refresh, mode switch, menu bar) | The view ribbon | the `ribbon` a `ViewHost.create` builder answers |
| Overall state ("86 of 120 reports", "Saving…") | The per-view status bar | the `status` a `ViewHost.create` builder answers |
| Details of the currently selected object | The context panel | `grok.shell.o = entity` |
| The view itself | Opens from Browse on single click | app function with `//output: view result` that returns the view |
| The app home's header (icon, what the app does, Learn more) | The top of the home view | an `AppHeader` built in the `ViewHost.create` builder |
| A package's actions over a platform TableView (a campaign's grid) | The view's ribbon, pinned at its end | `ViewHost.of(tv).ribbon({actions, buttons})` (u2/dg) |

A view is an owner through its `ViewHost` (u2/dg): `ViewHost.of(view)` for a view the platform made,
`ViewHost.create(…)` for a new one. Whatever is built under it goes when the shell removes the view.

## The shape

```ts
import {ViewHost} from '@datagrok-libraries/u2/dg';

//name: Reports Browser
//tags: app
//meta.browsePath: Dev
//output: view result
export function reportsBrowserApp(): DG.ViewBase {
  return ViewHost.create({name: 'Reports', content: () => {
    const {content, ribbon, status} = buildReportsBrowser();   // your builder returns the parts
    return {root: content, ribbon, status};
  }}).view;
}
```

`ViewHost.create` runs the builder under the view's host, so the content and the chrome are
built owned. It hosts the content as the view's root, released through the platform's kill
channel when the view closes (every effect goes with it), routes `ribbon` groups to
`view.setRibbonPanels`, and puts `status` — a `Signal<string>` for a live text panel, or
elements — into the shell's per-view status bar. The chrome is disposed together with the
content (the shell owns the ribbon/status containers, so the view's kill-walk alone would never
reach them). `host.view` is the view.

- A search `TextInput({search: true, inline: true})` is a fine ribbon citizen.
- Inline toolbars are still right for *panels within* a view — the ribbon is for the view's
  main controls.
- Empty / loading / error-with-retry states belong in the content area (the user acts
  there); the status bar carries the one-line summary.

## The app home

An app's home view opens with the same header across Datagrok: the app's big icon beside what
the app does, one bullet each, and Learn more. Then what the user starts from (a list, a form). The
whole home — `ViewHost.create` building the `AppHeader` (`icon`, `items`, `learnMore`) and a `CatalogTable`
under it — is the sample [docs/samples/app-home.ts](../samples/app-home.ts); the header alone, with every
option, is [docs/samples/app-header.ts](../samples/app-header.ts). Both compile and the header runs as a test.

`title` and `subtitle` put the app's name beside a smaller icon; `bottomLine` rules it off from
what follows; `learnMore: {onClick}` opens the app's own help instead of a page. The header
follows its own width, so it needs an owner: build it in the `ViewHost.create` builder, or, on a view
you made yourself (`DG.View.create()`), in `ViewHost.of(view).build(…)`, and it goes when the view closes.

## An app that lives in a platform TableView

The grid is the platform's; the view's host owns the app's session, its ribbon group, its
subscriptions and its dialogs, all ended by the view's close:

```ts
import {subscribe} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

const host = ViewHost.of(tv);                                // one per view, whichever wrapper you hold
const session = host.build(() => campaignSession(tv));      // built under the view: gone when the shell removes tv
host.effect(() => {
  const grid = host.signals.grid.value;                      // a layout swaps the grid: this re-runs,
  if (grid !== null)                                         // and the last run's subscriptions are gone
    subscribe(grid.onCellValueEdited, (cell) => session.recalc(cell));
});
host.subscribe(grok.events.onFileImportRequest, (args) => session.importFile(args));
host.ribbon({
  hideSave: true, name: 'hit-design-actions',
  breadcrumbs: {parent: landing, onClick: () => session.backToLanding()},
  actions: [
    {icon: 'sync', tooltip: 'Refresh Campaign', run: () => session.refresh(),
      hint: session.newChanges, hintDetail: 'New changes'},
    {icon: 'users', tooltip: 'Edit campaign permissions', text: 'Edit campaign permissions...',
      run: () => session.editPermissions(), visible: session.canEdit},
  ],
  buttons: [{text: 'Submit...', run: () => session.submit(), cls: 'hit-design-submit-button'}],
  status: session.statusText,
});
host.show(trackerDialog, {modal: true});                     // hidden while another view is current
host.dock(detailsPanel(), 'right', {title: 'Details'});      // gone when its panel closes, or with the view
host.task(async (t) => {
  const template = await loadTemplate();                     // nothing is building after an await:
  t.build(() => session.apply(template));                    // back in under the view; skipped once it is gone
});
```

- `host.signals` holds `isCurrent`, `grid`, `dataFrame`, `currentRowIdx` (a write moves it),
  `selection`, `filter`, `columns` and `rowCount`, followed through the platform's events.
- The group sits after a rule at the ribbon's end and stays visible while the ribbon overflows;
  below the width it needs, the icons fold into one "…" menu whose items read as the tooltips,
  and the buttons stay. An action's `hint` pulses its icon (folded: a dot on "…").
- Code after an `await` has no builder around it: run it in `host.task(async (t) => …)` and build
  with `t.build(…)`; start timers on `host.scope`.
- A docked panel is `host.dock(control, place)`: closing the panel (its ✕, or `dockManager.close`
  from code) disposes the control, and disposing the control closes the panel.
- The shell renders the crumb (`tv.parentView` + `tv.showParentCrumb`); a press on it makes the
  parent current, and `breadcrumbs.onClick` adds what the app does besides (close the view).
  `hideSave` is `tv.showSaveButton = false`. Never restyle `.d4-ribbon-name` from package CSS.

## Routing

The `path` a `ViewHost.create` builder answers writes the view's URL. A JS app owns the whole
string — the platform prefixes nothing — so **include the app root** (`/apps/<Package>/<App>`).
Answer a `computed` and every in-app navigation round-trips to the address bar; the shell replaces
rather than pushes, so it costs no history entries. After U2Demo's `package.ts`:

```ts
const APP_PATH = '/apps/U2demo/U2Demo';

export function u2DemoApp(path?: string): DG.ViewBase {
  return ViewHost.create({name: 'U2 Demo', content: () => {
    const {content, shell, status} = buildDemo({initial: leafForPath(path)?.id});
    return {root: content, status, ribbon: [[bar], demoRibbon(shell)],
      path: computed(() => APP_PATH + pathOf(shell.current.value))};
  }}).view;
}
```

The inbound direction is the platform's: give the app function a
`//input: string path { meta.url: true; optional: true }` parameter and the tail of the URL
arrives in it — resolve it to your initial state.

## Anti-pattern

```ts
// DON'T: chrome rebuilt inside the content — a second toolbar row above the data,
// a hand-made status strip below it, details locked in a pane the shell can't reuse
const bar = new Toolbar().add(search.root).addButton('Refresh', () => refresh());
const view = grok.shell.newView('Reports');
view.root.append(divV([bar, split, statusStrip]));
```

Standalone (gallery, no-shell prototypes) the same builder still works: render the ribbon
parts into an inline `Toolbar` and the status signal into a footer line yourself.
