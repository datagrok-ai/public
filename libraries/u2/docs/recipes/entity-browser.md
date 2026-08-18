# Recipe: entity browser (master-detail)

A virtualized list of server objects on the left, details for the selected one on the right,
main controls in the ribbon. Reference implementation: the U2Demo **Reports Browser**
(`packages/U2Demo/src/reports-browser.ts`).

## The shape

```ts
const query = signal('');
const pager = dapiPager<DG.UserReport>(() => grok.dapi.reports.include('reporter,assignee'), {
  pageSize: 30, order: 'createdOn', desc: true,
  filter: () => query.value.trim() || undefined,       // dapi smart-filter string
});

const list = new VirtualList({itemHeight: 44, keyOf: (r) => r.id, render: reportRow});
list.setItems(pager.items);                            // near-bottom scroll → pager.loadMore()

const current = computed(() => pager.items.value[list.selectedIndex.value] ?? null);
// selection → details pane AND the shell's context panel:
scope.effect(() => { if (current.value) grok.shell.o = current.value; });
```

The pieces, each written once and reused:

- **`dapiPager`** (u2/dg) owns paging, ordering, the smart-filter string, debounce and
  cancellation. Never hand-roll `page(n)` loops or spinners.
- **`VirtualList`** renders only visible rows — item count does not matter. Row rendering:
  see the [list-item-rendering](list-item-rendering.md) recipe.
- **Loading / empty / error-with-retry** live in the content area under the list; the
  one-line summary (`${loaded} of ${total}`) is a computed signal sent to the status bar.
- **Details**: a key/value pane (see [forms-and-details](forms-and-details.md)) and/or the
  context panel via `grok.shell.o`. Rebuild details under a per-selection `Scope` so the
  previous selection's components are disposed.
- **Chrome**: search input and refresh go to the ribbon via `appView` (see
  [shell-integration](shell-integration.md)).

## Anti-patterns

- Fetching everything up front and filtering client-side — use the dapi filter + pager.
- A non-virtualized `divV` of row elements — dies at thousands of rows.
- Re-rendering the whole list on selection — selection is a signal; rows stay put.
