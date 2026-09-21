# Viewers in `@datagrok-libraries/utils`

One viewer lives here: the **Forms viewer** (`forms-viewer.ts`), a `DG.JsViewer` that lays a table
out as a row of cards, one field per column, with a shared label column down the left. It is
implemented in this library but **registered by PowerGrid** (`packages/PowerGrid/src/package.ts:404-410`)
under the name `Forms`, so on a stand it is `[name="viewer-Forms"]` and a feature reaches it as
**`forms viewer`** — never `form viewer`, which is the Dart `Form` viewer (`viewer-Form`), a
different viewer with a different automation surface.

| File | Purpose |
|---|---|
| `forms-viewer.ts` | The viewer: cards, header labels, sorting, pinning, the three automation signals |
| `forms-viewer-status.ts` | The areas and readings `getWidgetStatus` reports |

## What a card is

`render` (`forms-viewer.ts:555-584`) rebuilds the header and hands the virtual view a card count of
`indexes.length + (showCurrentRow ? 1 : 0) + (showMouseOverRow ? 1 : 0)`. The **leading positions are
fixed**: position 0 is the current row when Show Current Row is on, position 1 (or 0) is the
mouse-over row when Show Mouse Over Row is on (`currentRowPos` / `mouseOverPos`, `:373-374`), and the
records follow. `renderForm` builds an empty card for a row of `-1`, which is why **`card 2` is an
empty card whenever nothing is hovered** — it is not missing, it is blank.

`indexes` holds the selected rows intersected with the filter (`:563-569`), sorted by the grid's sort
or by `sortBy` when set, minus any pinned rows (`:570-571`). Pinned cards live in their own pane
(`pinnedFormsDiv`) above the list, and that pane is `display: none` while nothing is pinned
(`:614`).

**The DOM order of the cards is not their position order.** `virtualView.refreshItem` re-appends the
card it rebuilt, so a card that has been refreshed sits last in the DOM. `listCards`
(`forms-viewer-status.ts:10-17`) sorts by `left`, then `top`, which is the order the user reads —
so **`card 1` is always the leading card**, whatever the DOM says.

## Automation surface (`getWidgetStatus`)

`parts`: `root`, `header` (`columnHeadersDiv`), `list` (the virtual view's root), `pinned`
(`pinnedFormsDiv`). There is no canvas part, so every hit area is in CSS px of the viewer's **root**.

`hitAreas`:

| Area | What it is, and when it is absent |
|---|---|
| `card <n>` | The n-th card by position, `n` from 1. Only cards the virtual view currently has in the DOM — scrolling changes which. |
| `current card` / `mouse-over card` | Aliases for the leading positions, present only while `showCurrentRow` / `showMouseOverRow`. |
| `pinned card <n>` | The n-th card of the pinned pane, `n` from 1. Absent while nothing is pinned. |
| `field <COL> of <label>` | One `[column]` element inside that card, for every card label above. Absent for a column that is not a field, and for a card whose row is `-1` (nothing is built for it). |
| `label <COL>` | The header label. |
| `remove <COL>` | The `✕` icon in the header. The label container toggles its `visibility` on hover (`:364-365`), which does not change its box — the area is reported whenever the header is laid out, hovered or not, and clicking it removes the field. |
| `sort indicator <COL>` | The `↑` / `↓` glyph, appended only to the column the viewer sorts by (`:338-340`). |

`values`:

| Reading | Reads from |
|---|---|
| `cards` | How many cards the virtual view has laid out. |
| `records shown` | The cards showing a real row **plus every pinned card** — a pinned row is a record on screen, so pinning one and leaving the list untouched moves this reading. |
| `pinned records` | Pinned cards only. |
| `fields shown` / `fields` | `fieldsColumnNames` — the **configured** field set, in order. |
| `header labels` | The ordered text of the labels the header actually drew. Distinct from `fields`: a field survives in `fieldsColumnNames` until the viewer prunes it, and `renderHeader` skips a column the built form has no element for (`:328`). This is the reading that makes "a dropped or `~`-renamed column prunes the field" honest. |
| `<COL> of <label>` | The field's text: the input's value, `true`/`false` for a checkbox, the literal `canvas` for a renderer-drawn field. |
| `width of` / `height of <COL> of <label>` | The field's box, rounded — what a Renderer Size claim compares. |
| `background of <COL> of <label>` | The computed background as `#RRGGBB`, `''` when transparent. |
| `align of <COL> of <label>` | `getComputedStyle(field).textAlign`. The viewer writes `textAlign` only for a `center` / `right` grid column style (`:472-473`), so a left-aligned column reads whatever the stylesheet gives (`start` or `left`), not `left` by the viewer's doing. |
| `font of <COL> of <label>` | The inline `font` shorthand on the field: the viewer's `font` property (`:443`), replaced by the grid column's `contentCellStyle.font` when it has one (`:475-476`). `''` for a canvas field — a renderer paints its own text. |
| `field kind of <COL>` | `canvas` or `input` — which branch of `buildForm` drew it. Recorded once, from the first card that has the field. |
| `record of card <n>` | The table row the card shows, **1-based**, `''` when the card is blank (no current row, nothing hovered). This is what a sort-mirror claim compares: an expected row order, instead of re-deriving `getSortedOrder` in the test. |
| `card kind of card <n>` | `current`, `mouse-over` or `record`; `card kind of pinned card <n>` is `pinned`. Derived from the position, which is what the `d4-multi-form-form-indicator-*` classes are derived from too. |
| `pinned pane shown` | `pinnedFormsDiv.style.display !== 'none'` (`:614`). |
| `pinned values` / `pinned by` | `pinnedRowValues` / `pinnedRowColumnNames`, comma-joined — the **persisted** identity of a pinned row, the pair a layout carries and `resolvePinnedRows` (`:637-653`) maps back to a row index. |
| `sort column` / `sort direction` | The first of `getSortByColumns()` and `↑` / `↓`; both `''` when nothing sorts. |
| `current record` / `mouse-over record` | The dataframe's current / mouse-over row, 1-based, `''` when none. |

`shortcuts` and `events` are empty, and `error` is always `null` — see below.

`isRenderPending` / `onRendered` (`forms-viewer.ts:59-72`) cover both places a render is deferred:
the debounced data subscriptions (`:214-234`) and the 110 ms the virtual view takes to lay cards out
while its root is still out of the document (`:586-600`).

## Writing features against this viewer

**The viewer has no error state.** An empty field set draws zero labels and zero fields and shows
nothing at all — no balloon, no banner, no message. The 20-column cap is silent the same way: the
source still declares `COLS_LIMIT_EXCEEDED_WARNING` (`:10`) but nothing reads it, the message having
been removed deliberately. So a feature asserts `fields shown` is `0` and `header labels` is empty,
never that the viewer reports an error.

**Read a card's fields, never the grid's.** `<COL> of card 1` is the text the card shows, which is
the viewer's own formatting (`numberFormat`, `:444-449`) — it is not the grid cell's text and is not
meant to be.

**A blank leading card is normal.** With Show Current Row and Show Mouse Over Row both on and nothing
hovered, `card 2` exists, has a box and has no field values. Assert `card kind of card 2` is
`mouse-over` and `record of card 2` is `''`; do not assert the card is absent.

**Pinning is by value, not by index.** `pinnedRowValues` / `pinnedRowColumnNames` are what a layout
persists, and a non-unique value warns on pinning (`:620-621`) because it will not survive the round
trip. A layout claim compares `pinned by` and `pinned values`, not row numbers.

## Rebuilding after a change here

This library is consumed as **TypeScript sources** (`@datagrok-libraries/utils/src/viewers/...`), so
a change here reaches a stand only through the packages that import it: rebuild and republish
**PowerGrid** (`npm run build && grok publish <alias>`), with `node_modules/@datagrok-libraries/utils`
pointing at this checkout (`npm link` or a junction). Publishing this library to npm is only needed
for consumers that install it from the registry.
