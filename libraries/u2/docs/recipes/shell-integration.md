# Recipe: shell integration

A u2 app runs inside the Datagrok shell, which already provides per-view chrome. Use it —
don't rebuild it inside the view content.

| App part | Where it goes | How |
|---|---|---|
| Main view controls (filter box, refresh, mode switch, menu bar) | The view ribbon | `appView({ribbon})` |
| Overall state ("86 of 120 reports", "Saving…") | The per-view status bar | `appView({status})` |
| Details of the currently selected object | The context panel | `grok.shell.o = entity` |
| The view itself | Opens from Browse on single click | app function with `//output: view result` that returns the view |

## The shape

```ts
import {appView} from '@datagrok-libraries/u2/dg';

//name: Reports Browser
//tags: app
//meta.browsePath: Dev
//output: view result
export function reportsBrowserApp(): DG.ViewBase {
  const {content, ribbon, status} = buildReportsBrowser();   // your builder returns the parts
  return appView({name: 'Reports', content, ribbon, status});
}
```

`appView` hosts `content` as a `DG.Widget` (view close disposes every effect), routes
`ribbon` groups to `view.setRibbonPanels`, and puts `status` — a `Signal<string>` for a live
text panel, or elements — into the shell's per-view status bar. Chrome components passed to
it are disposed together with the content (the shell owns the ribbon/status containers, so
the view's kill-walk alone would never reach them).

- A search `TextInput({search: true, inline: true})` is a fine ribbon citizen.
- Inline toolbars are still right for *panels within* a view — the ribbon is for the view's
  main controls.
- Empty / loading / error-with-retry states belong in the content area (the user acts
  there); the status bar carries the one-line summary.

## Anti-pattern

```ts
// DON'T: chrome rebuilt inside the content — a second toolbar row above the data,
// a hand-made status strip below it, details locked in a pane the shell can't reuse
const bar = new Toolbar(); bar.add(search.root); bar.addButton('Refresh', ...);
const view = grok.shell.newView('Reports');
view.root.append(divV([bar, split, statusStrip]).root);
```

Standalone (gallery, no-shell prototypes) the same builder still works: render the ribbon
parts into an inline `Toolbar` and the status signal into a footer line yourself.
