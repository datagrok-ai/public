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

## Actions: one list, two surfaces

Declare the item's actions once (`Action[]`), then:

```ts
const actionsOf = (r: Report): Action[] => [
  {name: 'Copy description', icon: 'copy', run: () => copy(r.description)},
  ...(r.jiraTicket ? [{name: `Open ${r.jiraTicket}`, icon: 'external-link-alt',
    run: () => window.open(ticketUrl(r))}] : []),
  {name: 'Copy id', run: () => copy(r.id)},          // no icon — menu only
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
  keyboard users. Each button carries the action name as tooltip and `aria-label`.
- **`contextActions`** on `VirtualList` wires right-click: the row is selected and the full
  list opens as a `Menu` at the cursor (`actionsMenu(actions)` does the same standalone).
  Hover shortcuts are never the only path.
- For platform entities, the chip's context menu already carries the object's registered
  actions — don't build a parallel list for the entity itself; item actions are for the
  row's own verbs.

## Timestamps: `timestamp(date)`

Compact visual, full truth on hover: `timestamp(r.createdOn)` renders the locale short date
(`Aug 11`; the year appears only when not current) with the full date-time in the tooltip,
styled small and light (`.u2-timestamp`). It accepts `Date`, epoch numbers, strings, or
anything with `toDate()` — dayjs included; invalid dates render empty.

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
