# Recipe: rendering complex items in a list

How to render rich items — reports, mail, issues, runs — in a `VirtualList`, the way Slack
or Outlook Web render theirs. An item is a small layout, not a line of text.

## Row anatomy

```
┌────────────────────────────────────────────────┐
│ Title, ellipsized…                [👤 J. Smith] │   primary line
│ Aug 11 ·                              [🗨 ✓ ⋯]  │   meta line; actions replace it on hover
└────────────────────────────────────────────────┘
```

- **Primary line**: the item's identity (subject / description), ellipsized, never wrapped.
- **Meta line**: secondary facts — timestamp, status dot, counts — in
  `--dg-font-size-small` and `--dg-text-color-light`.
- **People render as markup, not plain text**: icon + name via the entity renderer
  (`chip(user)` from u2/dg — handler markup, tooltip, click → context panel for free).

## Hover-activated actions

Show the 2–4 most common actions as an icon block that appears on row hover (top-right,
over the meta line), like Slack message actions or Outlook's archive/delete:

```ts
function reportRow(r: DG.UserReport): HTMLElement {
  const actions = divH([
    iconButton('check', () => resolve(r), {tooltip: 'Resolve'}),
    iconButton('external-link', () => openTicket(r), {tooltip: 'Open ticket'}),
  ], 'u2-row-actions');
  return divV([primaryLine(r), divH([metaLine(r), actions], 'row-meta')], 'row');
}
```

```css
.row .u2-row-actions { display: none; }
.row:hover .u2-row-actions { display: flex; gap: var(--dg-space-s); }
```

- Hover actions are shortcuts, never the only path: **right-click shows ALL actions** in a
  context menu (`Menu`). For platform entities, the chip's context menu already carries the
  object's registered actions — don't build a parallel list.
- Icons only, each with a tooltip; no text buttons inside rows.

## Timestamps

Compact visual, full truth on hover: render the short form (`Aug 11`; add the year when not
current), and put the full date-time in the tooltip/`title`. Small font, light color. Right
side of a line, so columns of dates align across rows.

```ts
const el = span(d.format(d.year() === now.year() ? 'MMM D' : 'MMM D, YYYY'), 'row-date');
el.title = d.format('MMM D, YYYY HH:mm');
```

## Adaptive rendering

A row must degrade gracefully as its pane narrows (dock panes get thin). Rank the parts and
hide from least important down — e.g. user *name* first, then the user icon — using
container queries on the list host:

```css
.report-list { container-type: inline-size; }
@container (max-width: 420px) { .row .row-user-name { display: none; } }
@container (max-width: 340px) { .row .row-user     { display: none; } }
```

Never let secondary items squeeze the primary line into wrapping; the primary line owns the
flexible space (`min-width: 0` + ellipsis), secondary items are `flex-shrink: 0` until their
breakpoint hides them.

## Anti-patterns

- Actions always visible on every row — noise at 40 rows; hover + right-click is the pattern.
- Full date-times printed in every row — that is what the tooltip is for.
- `white-space: normal` titles — rows must keep a fixed height for virtualization.
- Rendering the user as plain text when an entity renderer exists.
