# Recipe: rendering complex items in a list

How to render rich items — reports, mail, issues, runs — in a `VirtualList`, the way Slack
or Outlook Web render theirs. An item is a small layout, not a line of text. The patterns
here are u2 primitives; the reference implementation is the U2Demo Reports Browser
(`packages/U2Demo/src/reports-browser.ts`) and the gallery's "Rich rows" demo (`#list`).

## Row anatomy

```
┌────────────────────────────────────────────────┐
│ Title, ellipsized…                [👤 J. Smith] │   primary line
│ Aug 11 ·                              [🗨 ✓ ⋯]  │   meta line; actions reveal on hover
└────────────────────────────────────────────────┘
```

- **Primary line**: the item's identity (subject / description), ellipsized, never wrapped —
  rows keep a fixed height for virtualization.
- **Meta line**: secondary facts — timestamp, status dot, counts — in
  `--dg-font-size-small` and `--dg-text-color-light`.
- **People render as markup, not plain text**: icon + name via the entity renderer
  (`chip(user)` from u2/dg — handler markup, tooltip, click → context panel for free).

## Short keyed lists outside a VirtualList: `each` / `show`

For a handful of items rendered in place (chips, access rows, a popup's options), `each(host,
items, keyOf, build, {update?, empty?})` keeps one element per key: built once in its own scope,
patched or moved in place, disposed when its key leaves — focus and caret survive an edit.
`show(host, when, build, otherwise?)` is the conditional sibling. Both need an owner (call them
inside a `component((self) => …)` builder or `owner.build(…)`). Cards in a grid take `onClick(e, card)` — the event and the
card — and `selectMode: 'toggle' | 'select'`.

## Actions: one list, two surfaces

Declare the item's actions once (`Action[]`), then:

```ts
const actionsOf = (r: Report): Action[] => [
  {text: 'Copy description', icon: 'copy', run: () => copy(r.description)},
  ...(r.jiraTicket ? [{text: `Open ${r.jiraTicket}`, icon: 'external-link-alt',
    run: () => {
      window.open(ticketUrl(r));                     // a run returns a promise or nothing, not the window
    }}] : []),
  {text: 'Copy id', run: () => copy(r.id)},          // no icon — menu only
];

const list = new VirtualList<Report>({
  itemHeight: 44,
  contextActions: actionsOf,        // right-click: selects the row, opens the FULL list
  render: (r) => divV([
    primaryLine(r),
    divH([timestamp(r.createdOn), rowActions(actionsOf(r))]),   // hover: icon shortcuts
  ]),
});
```

- **`rowActions(actions)`** renders the icon-bearing subset as a hover-revealed block
  (`.u2-row-actions`, css/list.css). It hides by opacity, not display, so the buttons keep
  their space (no layout shift) and stay tabbable — tabbing into one reveals the block for
  keyboard users. Each button carries the action's `text` as tooltip and `aria-label`. While hidden the block
  takes no pointer, so a test reaches its buttons the way the keyboard does (`user.press(button, 'Enter')`): no
  script can hover. Any element carrying `u2-row` reveals a block inside it on hover too.
- **`contextActions`** on `VirtualList` wires right-click: the row is selected and the full
  list opens as a `Menu` at the cursor (`actionsMenu(actions)` does the same standalone).
  Hover shortcuts are never the only path.
- For platform entities, the chip's context menu already carries the object's registered
  actions — don't build a parallel list for the entity itself; item actions are for the
  row's own verbs.

## Timestamps: `timestamp(date)`

Compact visual, full truth on hover: `timestamp(r.createdOn)` spells the date out in the locale
(`Aug 11, 2026`, with the time of day where the value carries one) and the full date-time in the tooltip,
styled small and light (`.u2-timestamp`). It accepts `Date`, epoch numbers, strings, or
anything with `toDate()` — dayjs included; invalid dates render empty. It reads the date in the browser's zone; a
day stored as midnight UTC reads in UTC with `timestamp(day, undefined, {utcDates: true})` (see
[forms-and-details](forms-and-details.md#dates-and-the-text-of-a-value)).

## Adaptive rendering: `u2-adaptive` + `u2-p*`

A row must degrade gracefully as its pane narrows (dock panes get thin). Mark the list host
`u2-adaptive` and rank the hideable parts:

```ts
list.root.classList.add('u2-adaptive');
// in the row: the assignee name is the first thing to go, its icon the last
span(`@${r.assignee}`, 'u2-p2');
```

`u2-p2` hides below 420px of container width, `u2-p1` below 340px (css/adaptive.css).
Container queries cannot read custom properties, so the breakpoints are fixed — an app
needing different thresholds writes its own `@container` rules against the same classes.
The primary line owns the flexible space (`min-width: 0` + ellipsis); secondary items are
`flex-shrink: 0` until their breakpoint hides them.

## Anti-patterns

- Actions always visible on every row — noise at 40 rows; hover + right-click is the pattern.
- Full date-times printed in every row — that is what the tooltip is for.
- Wrapping titles — rows must keep a fixed height for virtualization.
- Rendering the user as plain text when an entity renderer exists.
- Hand-rolling the hover CSS or the context-menu listener — `rowActions` / `contextActions`
  already carry the accessible behavior.
