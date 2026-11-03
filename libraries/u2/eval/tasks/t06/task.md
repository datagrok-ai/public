# Rename campaigns

A context-panel pane lists the user's campaigns for quick renaming. `load()` fetches them from the server, which
takes a moment, so the pane shows that it is loading meanwhile. Then it shows one editable row per campaign: a
name field captioned with the campaign's code. Under the rows a line counts, as the user types, how many names
differ from the loaded ones: `0 renamed`, `1 renamed`, and so on.

The plugin calls:

```ts
export function campaignsPanel(load: () => Promise<{code: string, name: string}[]>): Control
```
