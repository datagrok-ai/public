# Markup (Playwright) — Run Results

**Date**: 2026-08-06
**URL**: https://dev.datagrok.ai (1.28.0)
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Markup from the Viewers toolbox | PASS | Default Markdown sample arrives rendered: `h1`, 4 list items, 4+ links |
| 2 | Edit content... replaces what the viewer renders | PASS | Dialog pre-filled with the current content; new heading and list rendered after **OK** |
| 3 | `${COLUMN}` renders the current row values | PASS | Values follow a grid click and then an arrow-down |
| 4 | An unknown `${COLUMN}` is left as written | PASS | `${NO_SUCH_COLUMN}` stays literal |
| 5 | The Markup engine renders table expressions live | PASS | `#{t.rowCount}` → 5850; **Select all** drives `#{t.selection.trueCount}` to 5850 and **Select none** back to 0 |
| 6 | Markup Enabled decides whether expressions are evaluated | PASS | Checkbox reads back `false`, the expression stays resolved — `[KNOWN_BUG_REPRODUCED:GROK-20637]` |
| 7 | Mode decides how the content is interpreted | PASS | `# Heading probe`: Auto → `h1`, None → `pre` with the source, Html → no heading, Markup → `h1` |
| 8 | Mode = None shows HTML as source | PASS | `<pre>` produced, but the `<b>` survives in it — `[KNOWN_BUG_REPRODUCED:GROK-20637]` |
| 9 | Bold content is rendered bold | PASS | `strong` produced, computed weight 400 — `[KNOWN_BUG_REPRODUCED:GROK-20637]` |
| 10 | Title and Show Title drive the viewer title bar | PASS | Title text shown, then hidden |
| 11 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 30.0s / 27.4s on two consecutive runs (was 15.6s before the three guarded steps, which spend their timeouts waiting for a change that never comes) |

## Notes

* This viewer renders into the DOM (`.grok-help` inside `.grok-help-host`), not a
  canvas, so every assertion reads elements — `h1`, `li`, `pre`, `a` — and no pixel
  helper is involved. That is also why it is the fastest viewer spec in the folder.
* The content is edited both ways a user can, and the two paths open **different**
  dialogs:
  * **Edit content...** in the viewer's context menu (`.d4-menu-popup
    [name="div-Edit-content..."]`) opens an *Edit* dialog with a plain
    `textarea.ui-input-editor`;
  * the **Content** row in the property grid opens through its ellipsis button
    (`button.property-grid-ellipsis-editor-ellipsis`) an editor dialog titled
    *Content*, backed by **CodeMirror** — there is no textarea to fill, so the text
    goes in as real keystrokes. Keep such content on one line: the markdown editor
    continues lists and indentation by itself.
* Right-click the **bottom** of the viewer. The default sample is full of links,
  and a right-click on one of them opens the browser's own context menu instead.
* Acting on the table replaces the Context Panel, so the spec re-opens the viewer's
  properties by probing for a row only this viewer has
  (`tr[name="prop-markup-enabled"]`) and checking it **attached** rather than
  visible — the row lives in a category that may be collapsed.
* The current row has to be set from the UI first: nothing is current on a freshly
  opened table, and the viewer substitutes only for a current row. Clicking the
  leftmost strip of the grid does not do it — that is the row header, and it
  selects rows instead.

## Open bugs — GROK-20637, guarded by `knownOpenBug()`

| What | Detail |
|---|---|
| Markup Enabled does nothing | Declared in `markup_viewer_look.dart` and never read anywhere in `d4`; the content always goes through the Markup engine, so `#{t.rowCount}` keeps rendering with the option off |
| Mode = None does not escape | It wraps the content in `<pre>` without escaping, so the browser parses `<b>bold probe</b>` and the markers disappear instead of being shown as source. Markdown genuinely is not interpreted, which step 7 asserts hard |
| Bold is never bold | `.grok-help` sets the relative `font-weight: lighter` (= 100) and the UA's `b, strong { font-weight: bolder }` resolves to 400 against it, so emphasis lands on ordinary weight. Measured, not eyeballed: the `<b>` element exists in the DOM while the text looks perfectly plain — which is why step 9 reads `getComputedStyle` |

## Not automated

* External `<iframe>` embedding, the **Markup view** from the **+** menu, the
  **Stretch** property and write-back input fields — kept as a manual list in
  `markup.md`.
