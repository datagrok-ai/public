# Word cloud

An echarts viewer (`echarts` + the `echarts-wordcloud` extension) that counts one string column and
draws each distinct value at a size proportional to how many rows carry it. Registered by Charts as
`Word cloud`, so a feature reaches it as **`word cloud viewer`**.

| File | Purpose |
|---|---|
| `word-cloud-viewer.ts` | The viewer: column pick, counting, the echarts option, the three automation signals |
| `word-cloud-status.ts` | The areas and readings `getWidgetStatus` reports |

## What it counts

`render` (`:181-185`) walks **`this.filter`** — the viewer's own combined filter, not the whole
table — and counts the column's values into a `Map`. That map is the viewer's record of the frame:
it is what `rows of word` reports, and it is why filtering must move the counts. A word's *name* is
unchanged by a filter that leaves at least one row of it, so a claim about filtering that compares
names only proves nothing.

The column is picked in `onTableAttached` (`:107`): the first string column with more than 1 and at
most `MAX_UNIQUE_CATEGORIES_NUMBER = 500` distinct values. On **demog-1000** that is RACE —
Caucasian 896, Other 62, Black 27, Asian 15; SEX is F 553 / M 447; USUBJID has 1000 distinct values
and trips the 500 gate.

## The two message states

`render` returns early, **before creating the chart**, when there is no string column at all
(`:156-159`) or the chosen column exceeds the gate (`:161-164`). It writes the message into the root
as a `.d4-viewer-error` div and leaves the previous chart standing — `this.chart` still holds the
last cloud's geometry, and the root was never emptied.

The status therefore reports `error` and **nothing else about a cloud**: no `word` area, no `view`
area, no `canvas` part, no `words` / `word names` / `rows of word` readings. A hidden thing must not
keep the geometry of the frame that drew it. `column` and `rows shown` are still reported — they are
the viewer's state, not the picture's.

## Automation surface (`getWidgetStatus`)

`parts`: `root`, and `canvas` — the echarts canvas, absent in a message state. Hit areas are in CSS
px of that canvas, which is where zrender's own coordinates are.

| Area | What it is, and when it is absent |
|---|---|
| `view` | The chart canvas box. Absent in a message state. |
| `word "<name>"` | The axis-aligned box of one laid-out word: the series' graphic element's bounding rect through its computed transform, so rotation is included. **Absent for a word the layout dropped** — the graphic element is the only record the layout leaves, and a word it could not place has none. Also absent for every word in a message state. |

| Reading | Reads from |
|---|---|
| `column` | `wordColumnName`. Reported in every state. |
| `rows shown` | `filter.trueCount` — the population the counts are over. Reported in every state. |
| `words` | How many words the layout **drew**. |
| `word names` | Those words' names, comma-joined, in series order. |
| `rows of word "<name>"` | The viewer's own count for that word over its filter. Reported for **every counted word**, drawn or not — a word can have a count and no area, which is exactly the `drawOutOfBound` / no-room case. |
| `font` | `<fontWeight> <fontFamily>` as it reached the first drawn word — the applied font, not the property. Absent when nothing was drawn. |

`shortcuts` and `events` are empty: the viewer binds no keys and fires no custom event.

## `isRenderPending` and `onRendered`

The flag goes up when a change **arrives** — `render` sets it, and so does each raw data stream
before its 50 ms debounce (`:112-118`), so a settle right after a filter change does not read the
previous frame.

It comes down one macrotask plus one animation frame after `setOption`. That is not arbitrary:
`WordCloudSeries` does not declare `layoutAnimation`, the layout helper copies the undefined value
over its own `true` default (`echarts-wordcloud/src/layout.js:229, 250-256`), and the falsy branch
(`:1263-1270`) lays **every** word out inside the single macrotask that `setOption` queued; zrender
paints them on the next frame. echarts' own `finished` event is not the signal — the wordcloud
layout is not part of echarts' scheduler, so `finished` fires on the empty chart before the layout
tick has run, and again after it.

`render` disposes and re-creates the chart on every pass (`:195-198`), so nothing survives a render
except the counts and the error the viewer records itself.

## Writing features against this viewer

**Point at a word by name, never by pixel.** `word "Caucasian"` is exact; scanning the canvas in
buckets and hovering candidates is what these areas replace.

**Text-size and rotation claims are box claims.** Equal `minTextSize` and `maxTextSize` means the
`word` areas are the same height; a rotation range shows up as the boxes' aspect. Comparing whole-
canvas ink before and after proves nothing — three different layouts of the same words always
differ.

**Filtering claims compare counts.** `rows of word "Caucasian"` falls; `word names` may not change
at all.

**The 500-category gate is a message, not an empty cloud.** Assert `error` contains
`500 or fewer unique categories` and that the viewer reports no `word "…"` area.

## Rebuilding after a change here

Charts is a published package: `npm run build` (which runs `grok api`, `grok check --soft` and
webpack) then `grok publish <alias>` — or `npm run build-charts` for both. The features see the new
keys only after that publish; nothing in the client or the Dart core has to be rebuilt.
