# Current view in the status bar

A small status-bar badge shows the name of the current view, `grok.shell.v?.name`, or **No view** when there is
none, and keeps up as the user switches views (the shell fires `grok.events.onCurrentViewChanged`).

The plugin calls:

```ts
export function currentViewBadge(): Control
```
