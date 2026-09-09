---
feature: markup
target_layer: playwright
boot_lane: server
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - markup-spec.ts
related_bugs:
  - id: GROK-20637
    status: fixed
    fixed_in: 1.28.0
---

# Markup (Playwright)

All scenarios start with:

1. Close all
2. Open **System:DemoFiles/demog.csv**
3. Add **Markup** from **Toolbox > Viewers**

This viewer renders into the DOM rather than a canvas, so every assertion reads
the rendered markup itself — the heading elements, the list items, the `<pre>`
block — instead of counting pixels. The content is edited the two ways a user
can: the **Edit content...** item in the viewer's context menu, and the
**Content** editor in the Context Panel property grid.

## Add the viewer

1. Click the **Markup** icon in **Toolbox > Viewers**
2. The viewer opens on the default Markdown sample, already rendered: a heading,
   a four-item list and several links — not the Markdown source

## Edit content

1. Right-click the viewer and choose **Edit content...**
2. The dialog opens pre-filled with the content currently on screen
3. Replace it with a Markdown heading and a two-item list referencing columns:

   ```markdown
   # Demographics

   * Age: ${AGE}
   * Sex: ${SEX}
   ```

4. Click **OK** — the viewer renders the new heading and both list items

## Column references

1. Click a cell in the Grid to make a row current
2. The `${AGE}` and `${SEX}` references render that row's values
3. Press the down arrow to move to the next row — the rendered values follow it
4. A reference to a column that does not exist is left exactly as written

## Markup engine expressions

1. Set the content to `Rows: #{t.rowCount} Selected: #{t.selection.trueCount}`
2. Both expressions render as numbers, and the selected count starts at 0
3. Click **Select all** on the toolbar — the selected count becomes the row count
4. Click **Select none** — it goes back to 0

## Markup Enabled

1. With the content set to `Rows: #{t.rowCount}`, the expression renders as a number
2. Uncheck **Markup Enabled** in **Misc** — the expression is left as written
   (GROK-20637, fixed in 1.28.0)
3. Check it again

## Interpretation mode

With the content set to a single Markdown heading (`# Heading probe`), the four
modes are told apart by whether an `<h1>` is produced:

1. **Auto** (default) — the content has no leading `<`, so it is read as Markdown
   and the heading is rendered
2. **None** — nothing is interpreted; the source is shown as preformatted text
3. **Html** — the content is inserted as HTML, where a hash is just a hash, so no
   heading appears
4. **Markup** — Markdown again, whatever the content looks like

With the content set to `<b>bold probe</b> plain probe`, **None** shows the
source with its angle brackets (GROK-20637, fixed in 1.28.0).

## Bold

1. In **Markup** mode set the content to `**md bold** plain tail`
2. A `<strong>` is produced, and *md bold* is visibly heavier than the tail
   (GROK-20637, fixed in 1.28.0 — the spec reads the computed weight rather than
   trusting the element's presence)

## Title

1. Set **Title** in **Context Panel > Description** to *Patient card*
2. Check **Show Title** — the viewer title bar shows the text
3. Uncheck it — the title bar is hidden again

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Bugs this scenario asserts against

**GROK-20637** (fixed in 1.28.0, verified on dev 2026-09-03) covered three defects
that are now asserted hard: **Markup Enabled** was never read, **Mode = None** did
not escape the content, and `.grok-help`'s relative `font-weight: lighter` cancelled
the browser's `bolder` on `strong` so emphasis landed on the ordinary weight.
