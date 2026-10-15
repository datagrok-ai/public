# Chord

A circos chord diagram in svg: one arc per category of the From and To columns, one chord per
From–To pair the aggregation found. Registered by Charts as `Chord`, so a feature reaches it as
**`chord viewer`**.

| File | Purpose |
|---|---|
| `chord-viewer.ts` | The viewer: column pick, aggregation, ordering, the circos layout, the automation signals |
| `chord-status.ts` | The areas and readings `getWidgetStatus` reports |
| `utils.ts` | The layout configuration (shared by every chord viewer) and the topological sort |

## Automation surface (`getWidgetStatus`)

There is no canvas: `parts.root` is the svg (the viewer root when there is none), and every area is
in px of the svg's box. Everything is read from the svg the last render built.

| Area | What |
|---|---|
| `category <name>` | a 4 px square inside the category's arc, at mid-radius |
| `chord <from> -> <to>` | a 4 px square inside the chord, just in from the circle at one of its ends |

An area is reported only at a point inside the shape (`isPointInFill`) that `elementFromPoint` also
resolves to it. The arc radii are measured on this svg's own arcs: the layout configuration in
`utils.ts` is one object every chord viewer writes. An empty value is `(empty)`. With From = To no
chords are drawn.

| Reading | What |
|---|---|
| `categories` | category arcs drawn |
| `chords` | chords drawn |
| `rows shown` | `filter.trueCount` — the viewer's own combined filter |

## `isRenderPending` and `onRendered`

Both come from `RenderSignals` (`utils/utils.ts`). `render` rebuilds the svg synchronously and fires
`onRendered` on its way out, on every path. The selection, size and filter subscriptions are
debounced 50 ms; each holds the flag from the first event of a burst until its own render. A filter
change renders twice — at once through `onSourceRowsChanged`, and 50 ms later through the debounce —
and the flag stays up until the second.
