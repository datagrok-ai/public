# Recipe: a catalog of records (CatalogTable)

The records a person comes back to — campaigns, profiles, templates, saved queries — listed whole
under an app's home: a real table whose columns fit their content, sorted by a header gesture,
grouped into foldable sections with counts, its columns chosen from a menu, each row with its
actions, and every choice kept. The reference is the U2Demo page Collections › Catalog table
(`packages/U2Demo/src/pages/catalog-table.ts`).

Pick it when the whole list fits the page (tens to a few hundred rows: every row is in the DOM,
real table semantics). Thousands of rows belong in `DataTable` (virtualized); one line or a card
per item in `VirtualList` ([list-item-rendering](list-item-rendering.md)).

## Composition

```
Continue Campaign  [👁] [≡] [⟳]          ← h1 + table.chrome(), 2px apart
  Name ▲            Code   Created   …  Status      [👥][🗑]   ← the name pinned first, actions pinned last
  ˅ Kinase hit-to-lead  11                                    ← a group: chevron, label, count
    PI3Kδ quick wins KHL-12 2026/10/19 … ● In Progress
  37 campaigns · 3 view only                                  ← footer
```

```ts
const campaigns = new AsyncValue(() => loadCampaigns());          // first load: a skeleton
const table = new CatalogTable<Campaign>({
  name: 'campaigns', label: 'Campaigns',
  items: campaigns,                                                // rows stay while a refresh runs
  keyOf: (c) => c.code,
  columns: computed(() => [
    {name: 'Name', value: (c) => c.friendlyName ?? c.code, pinned: true, link: 'Continue Campaign', maxWidth: '320px'},
    {name: 'Code', value: (c) => c.code},
    {name: 'Created', value: (c) => c.createDate},                 // 2026/3/30 sorts as a date
    {name: 'Molecules', value: (c) => c.rowCount, align: 'right'},
    {name: 'Status', value: (c) => c.status, render: (c) => divH([dot(variantOf(c.status)), span(c.status)], {gap: 's'})},
    ...fieldsOf(campaigns.value).map((f): CatalogColumn<Campaign> => ({name: `campaignFields.${f}`, header: f,
      value: (c) => c.campaignFields?.[f], menuGroup: 'Campaign Fields', menuDetail: templatesWith(f)})),
  ]),
  groupings: [{name: 'Template', key: (c) => c.template.key, groupLabel: (key, first) => first.template.name},
    {name: 'Status', key: (c) => c.status}],
  defaultColumns: ['Code', 'Created', 'Molecules', 'Status'],
  defaultOrder: (a, b) => b.modified - a.modified,                // newest first while nothing is sorted
  actions: (c) => [
    {text: 'Manage campaign permissions', icon: 'users', requires: 'edit', run: () => share(c)},
    {text: 'Delete campaign', icon: 'trash-alt', requires: 'edit', destructive: true, run: () => remove(c)},
  ],
  access: (c) => editRights(c),                                    // 'edit' | 'view' | 'pending', or a signal of one
  actionsVisible: 'always',
  viewOnly: {icon: 'eye', tooltip: 'View only: you can open this campaign but not change it'},
  onActivate: (c) => open(c),
  empty: 'No campaigns yet. Start the first one under New Campaign.',
  footer: (shown) => `${shown.length} campaigns`,
  persist: 'my-app-campaigns',
});
host.append(divH([h1('Continue Campaign'), ...table.chrome()]), table.root);
```

## What each part does

- **Columns** are live (`LiveOption`): fields that arrive with the items become columns when they do,
  and a column name kept from an earlier session shows once its column exists. `value` is what sorts
  and groups; `render` only draws. The first `pinned` column stays put while the table scrolls sideways
  in a narrow host, as does the actions column at the end; a 10px edge shades what runs under it.
- **Sort**: a double click on a header (`sortGesture: 'click'` for a single one; Enter or Space on a
  focused header either way) cycles ascending → descending → the default order, with the arrow and
  `aria-sort`. `type: 'auto'` reads the whole column: numbers when every value is one, dates when every
  value is a `Date` or date text (2026/3/30, 2026-03-30), else text in natural order (KHL-2 before
  KHL-12 — a code is never taken for a date). Empty values go last both ways and ties keep the default
  order. A column's `compare`, or the table's `sorter` (the whole list, for an app whose order probes
  the column its own way), replace it. `CatalogTable.sortItems` is the same comparator for code outside the
  table. The order is one string, as every u2 table keeps it (`SortOrder`): `'Code'`, `'!Created'`, `''` for the
  default order; a kept `{column, ascending}` still reads.
- **Groups** are computed over every item, start open, show their count and fold on a click; folds are
  kept per grouping. Order: first appearance in the current order, the empty group (`emptyGroupLabel`)
  last; a grouping's `compare` orders them its way. One group shows no header unless
  `singleGroupHeader`.
- **Menus**: `table.chrome()` gives the columns, group and refresh icons, each opening its menu under
  itself (tooltips through `CatalogChromeOptions`). A column's or grouping's `menuGroup` makes a
  submenu, `menuDetail` the grey text after an item. Columns turned on are appended, so the shown order
  is the order of turning on; the columns menu stays open while they are toggled (a `keepOpen` menu item).
  `columnsMenu()` and `groupMenu()` build the same menus for a ribbon.
- **Row actions**: the icon-bearing `Action`s show in the last cell — on the hovered row by default,
  on every row with `actionsVisible: 'always'`; a right click (or the context-menu key) opens all of
  them as a menu. `requires` hides an action where the row's access is not 'edit'; 'pending' keeps
  their room empty, so nothing shifts when the rights arrive.
- **States**: the first load shows `skeleton()`, no items the `empty` text, a failed first load the
  reason with Retry. A refresh keeps the rows (`aria-busy` on the table, the refresh icon spins); a
  failed one says why in a balloon.
- **Code and assistants** write `table.state` (`{columns, sort, groupBy, collapsed}`) or call
  `sortBy`, `groupBy`, `showColumn`, `foldGroup`. Every change is saved through `persist`.

## What a hovered row does

Every u2 table marks the row under the pointer, as the platform grid marks its own; nothing else changes on hover,
so nothing moves and nothing is reachable only by hovering:

| Table | A hovered row |
|---|---|
| `CatalogTable`, `BasicTable` | the grid's mouse-over tint (`--dg-grid-mouse-over-row`) and its stripe at the row's start; a selected row keeps the selection colour, its pinned cells too. A catalog's icon actions show on the hovered row (by default) — they keep their room while hidden and are reached by Tab or a right click |
| `DataTable`, `VirtualList`, `VirtualTree` | the subtle hover fill (`--dg-hover-subtle`); a `rowActions` block in a list row shows while the row is hovered or a button in it is focused |
| `tableFromMap` | no tint: a detail pane is read, not acted on row by row |

The colours are tokens (`--u2-row-hover-bg` and `--u2-row-hover-stripe`, the grid's own by default), never app
rules over u2's classes. A test cannot hover (no script can produce `:hover`): it reaches a hover-revealed action
with the keyboard, `user.press(button, 'Enter')`.

## Keeping the choices where an app already keeps them

`persist` takes a key (the whole state as JSON under `u2.catalog.<key>`) or an adapter, so an app
that has its own storage keeps its keys and formats:

```ts
const persist: CatalogPersistence = {
  load: () => ({
    columns: readColumns(localStorage.getItem('MyColumns')),
    sort: fromLegacySort(localStorage.getItem('MySorting')),
    groupBy: localStorage.getItem('MyGrouping') ?? null,
  }),
  save: (s) => {
    localStorage.setItem('MyColumns', JSON.stringify(s.columns));
    localStorage.setItem('MyGrouping', s.groupBy ?? 'None');
  },
};
```

`load` runs once; anything of the wrong shape is ignored. `save` runs after each change, not on load.

## Tests and automation

The control is kind `catalog-table` (`data-u2-name` from `name`). Its items are the group headers
(`section:<key>`, `expanded`) and the rows (by `keyOf`, labelled by the first column), a row of a folded
group reported but not rendered; `revealItem` unfolds its group. `extra` reads `sort` ('Molecules desc'),
`groupBy`, `columns` and `hidden` (comma lists), `collapsed`, `groups` ('In Progress (8), Submitted (3)'),
`rows` and `total`. Header cells are kind `column-header`, so a feature double-clicks
`"Molecules" column header in campaigns catalog table`.

## Anti-patterns

- A `ui.table` with listeners added on its header cells and accordions injected as rows — every sort,
  group and toggle rebuilds the table and loses the rest.
- A cap on the rows shown ("Show all") — the catalog is the whole list.
- Sorting the items yourself and handing them over: the table owns the order, so its state block, its
  persistence and its arrow agree.
- Icon colours, row heights or a status dot's gap in app CSS — the table carries them.
