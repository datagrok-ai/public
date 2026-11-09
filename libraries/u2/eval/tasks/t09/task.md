# Campaign actions on the ribbon

A campaign opens in a platform TableView, `tv`. Put the campaign's actions on the view's ribbon, as icons with
tooltips:

- **Add new row** (a plus) appends a row to the campaign's table: `tv.dataFrame.rows.addNew(null, true)`.
- **Refresh Campaign** (a sync icon) reloads the campaign: `session.refresh()`.

The plugin calls this once the view is open:

```ts
export function campaignRibbon(tv: DG.TableView, session: {refresh(): Promise<void>}): void
```
