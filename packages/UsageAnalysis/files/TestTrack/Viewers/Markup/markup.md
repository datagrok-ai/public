---
feature: markup
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - markup-spec.ts
related_bugs:
  - id: GROK-20637
    status: open
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
   (GROK-20637 — it stays resolved; the property is never read)
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
source with its angle brackets (GROK-20637 — the tag is parsed instead, so the
markers disappear).

## Bold

1. In **Markup** mode set the content to `**md bold** plain tail`
2. A `<strong>` is produced, and *md bold* is visibly heavier than the tail
   (GROK-20637 — it is not; `.grok-help` uses `font-weight: lighter`, and the
   browser's relative `bolder` on `strong` cancels it back to ordinary weight)

## Title

1. Set **Title** in **Context Panel > Description** to *Patient card*
2. Check **Show Title** — the viewer title bar shows the text
3. Uncheck it — the title bar is hidden again

## Closing the viewer

1. Click **Close** on the viewer title bar — the viewer is gone

## Open bugs guarded by this scenario

All three are parts of **GROK-20637**. Each is asserted for the DESIRED behaviour
and wrapped in `knownOpenBug()`, so the step is green while the bug reproduces
and goes loud once it is fixed.

* **Markup Enabled does nothing.** The property is declared in
  `markup_viewer_look.dart` and never read anywhere in `d4` — the content always
  goes through the Markup engine. Switching it off leaves `#{t.rowCount}` and
  friends rendering as numbers.
* **Mode = None does not escape the content.** It wraps the text in `<pre>`
  without escaping it, so the browser parses `<b>bold probe</b>` and the markers
  vanish instead of being shown as source. Markdown is genuinely not interpreted,
  which is what the spec asserts hard.
* **Bold is never rendered bold.** `.grok-help` sets the relative
  `font-weight: lighter`, which computes to 100; the browser's own
  `b, strong { font-weight: bolder }` is relative too and resolves to 400 against
  it. The two cancel, so emphasis lands on the document's ordinary weight — the
  spec reads the computed weight rather than trusting the element's presence.

## Manual scenarios (not automated)

> Manual

1. Add the viewer from the **Add viewer** dialog rather than the toolbox icon
2. Embed an external page and check that it works inside the platform:

   ```html
   <iframe fremeborder="0" id="iframe_opkomst" src="https://dirkmjk.nl/files/articles/2016/opkomst/en.html"
   width="100%" height="100%">
   </iframe>
   ```

3. Add a **Markup view** from the **+** menu (a view, not a viewer), check the
   sample text on it, and edit that text from the Context Panel
4. **Stretch** — switch it on and off and watch how the content is laid out
5. Edit the values of a row through `<input data-field="COLUMN">` fields placed
   in the content — typing into them writes back into the current row

---
{
  "order": 26,
  "datasets": ["System:DemoFiles/demog.csv"]
}
