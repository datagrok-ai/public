---
feature: formsviewer
realizes_atlas:
  - formsviewer.cp.interactions-and-row-binding
realizes:
  - viewers.form
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs: []
realized_as:
  - formsviewer-interactions-and-row-binding-spec.ts
expected_results:
  - anchor: "Scenario 1"
    expectation: >-
      After clicking a non-current card, df.currentRowIdx equals that card's
      row; a field value of the clicked card matches what df.get reports for
      that (row, column) pair.
  - anchor: "Scenario 2"
    expectation: >-
      After clicking a field inside a card, df.currentCell.column.name equals
      the clicked field's column attribute and df.currentCell.rowIndex equals
      the card's row.
  - anchor: "Scenario 3"
    expectation: >-
      After clicking a column label in the Forms viewer header pane,
      df.currentCol.name equals that label's column name.
  - anchor: "Scenario 4a"
    expectation: >-
      Starting from a selected card, after Ctrl+clicking it once
      df.selection.get(row) is false for that card's row; after a second
      Ctrl+click on the same card df.selection.get(row) is true again — a
      toggled round-trip on the BitSet.
  - anchor: "Scenario 4b"
    expectation: >-
      After Shift+clicking a card at index N (0-based in the full row set), the
      selected index set equals exactly the rows 0 through N inclusive; rows
      after N are not selected.
  - anchor: "Scenario 4c"
    expectation: >-
      After Ctrl+Shift+clicking a card at index N, the rows 0 through N
      inclusive are deselected (their BitSet bits are cleared); any rows outside
      that range that were previously selected remain selected.
  - anchor: "Scenario 5a"
    expectation: >-
      While the pointer hovers over a card, df.mouseOverRowIdx equals that
      card's row; after the real pointer is moved away to a neutral point the
      mouse-over binding releases (df.mouseOverRowIdx is -1) and the mouse-over
      card (position 1) holds no field value, while the card set is otherwise
      unchanged (position 0 stays the current-row card).
  - anchor: "Scenario 5b"
    expectation: >-
      While a grid row is hovered, the mouse-over card (the card immediately
      after the current-row card) fills with that row's field values; the field
      values go from empty (or a prior row's values) to the hovered row's
      values. Moving the pointer off the grid row returns the mouse-over card
      fields to empty.
  - anchor: "Scenario 6a"
    expectation: >-
      After turning Show Mouse Over Row OFF and waiting for the debounced
      re-render, the mouse-over position card is absent from the card area; the
      total card count drops by exactly one compared to the pre-toggle count;
      hovering a grid row no longer fills any card with that row's values.
  - anchor: "Scenario 6a restore"
    expectation: >-
      After turning Show Mouse Over Row back ON and waiting for the debounced
      re-render, the total card count returns to the pre-toggle baseCount.
  - anchor: "Scenario 6b"
    expectation: >-
      After turning Show Current Row OFF and waiting for the debounced
      re-render, the current-row card and its green current-row indicator are
      absent from the card area.
---

# Forms Viewer — Mouse Interactions and Row Binding

## Setup

1. Close all open views and dialogs.
2. Open **System:DemoFiles/demog.csv** so the table view is active.
3. Add the **Forms** viewer from **Toolbox > Viewers** (click the Forms icon in
   the viewer toolbox on the right side of the screen).
4. Verify that **Show Selected Rows**, **Show Current Row**, and **Show Mouse
   Over Row** are all ON in the Forms viewer property panel (Context Panel >
   Misc). These are the defaults; confirm rather than toggle blindly.
5. In the grid, click row 5 to make it the current row and select a handful of
   rows (Ctrl+click rows 2, 5, 10, 15) so the ordinary card area has content to
   interact with. Wait for the Forms viewer to re-render (allow for the 50 ms
   debounce settle).

## Scenarios

### Scenario 1: Card click sets current row

Steps:
1. Identify one of the selected-row cards that is NOT the card
   carrying the green current-row indicator. Note the row this card belongs to
   by reading the value of the USUBJID field on that card.
2. Click anywhere on the body of that card (not a field's label — click the
   card's background area or any visible field value).
3. After the re-render settles, confirm the grid highlights the row that card
   belonged to as the new current row (the row highlight moves in the grid).
4. Also verify that the SEX field value shown on that card matches the value
   displayed in the grid for the same row and column.

Expected:
- After clicking a non-current card, the grid's current-row highlight moves to
  that card's row; the field values shown on the card match the grid's values
  for that same row.

### Scenario 2: Field click sets current cell

Steps:
1. Inside any ordinary card, click directly on the value area of the **AGE**
   field (the data entry area for AGE inside the card body). Do not click the
   label in the header pane — click the value area in the card itself.
2. After the click settles, confirm that the AGE column is highlighted as the
   current column in the grid and that the current row matches the card's row
   (the grid's current-cell indicator lands on that row's AGE cell).

Expected:
- After clicking a field inside a card, the grid's current cell moves to that
  field's column ("AGE") and the card's row.

### Scenario 3: Column label click sets current column

Steps:
1. In the Forms viewer header pane (the left column listing field names), click
   the label for the **HEIGHT** column (the text "HEIGHT" in the header area,
   NOT the field value inside a card).
2. After the click settles, confirm that the HEIGHT column is highlighted as the
   current column in the grid.

Expected:
- After clicking a column label in the Forms viewer header pane, the HEIGHT
  column becomes the current column in the grid.

### Scenario 4: Modifier-click selection round-trips

Steps:
1. Establish a clean entry state: clear the current selection (click an empty
   area in the grid with no modifier key, or press Escape if a selection
   shortcut is available), then click row 20 in the GRID with no modifier key so
   that row 20 becomes the CURRENT row. Wait for re-render.
2. Now Ctrl+click row 20 in the GRID to select it as well. Keep this order —
   current row first, selection second — so that making the row current cannot
   undo the selection you are about to test. Then work with the leading card that
   carries the green current-row indicator, not with the ordinary card: an
   ordinary card vanishes the moment its row is deselected, so there would be
   nothing left to click a second time, while the current-row card stays in place
   whatever the selection does — and a Ctrl+click on it toggles the selection just
   the same. Ctrl+click the current-row card and wait for the re-render.
3. Confirm that row 20 is now deselected — the grid no longer shows it as selected
   and no ordinary card for row 20 remains (the Ctrl+click toggled it off since it
   was already selected). This is the toggle half of the round-trip. Label this
   **Scenario 4a** first half.
4. Ctrl+click the SAME current-row card again. Confirm row 20 is selected again —
   its ordinary card reappears in the ordinary-card area. This confirms the toggle
   is a round-trip. Label this **Scenario 4a** complete.
5. Reset selection: clear it and Ctrl+click rows 3, 7, 12, 18, 22 in the grid
   (five rows) so there are enough ordinary cards. Wait for re-render.
6. In the Forms viewer, Shift+click the card for row 12 (the card corresponding to
   grid row 12 — the third of the five selected rows when they are read in row
   order: 3, 7, 12, 18, 22). As in Scenario 4a, actuate this on the stable leading
   current-row card: make row 12 the current row first so its card sits at the fixed
   leading position, then Shift+click that leading card. A Shift+click on it drives
   the same 0..N range-select handler as clicking row 12's own ordinary card would,
   but the leading current-row card stays put across the re-render, whereas an
   ordinary card can shift or detach mid-gesture. Wait for re-render.
7. Confirm that rows 0 through 12 inclusive are now selected (the grid shows
   those rows highlighted) and that rows above 12 (rows 18 and 22) are no
   longer selected. Label this **Scenario 4b**.
8. Do NOT reset the selection here — carry over the state Scenario 4b left behind
   (rows 0 through 12 selected). Confirm before clicking that rows 9, 10, 11 and
   12 are still selected; they are the rows above the click that must survive it.
9. In the Forms viewer, Ctrl+Shift+click the card for row 8. As in Scenarios 4a and
   4b, actuate this on the stable leading current-row card: make row 8 the current row
   first so its card sits at the fixed leading position, then Ctrl+Shift+click that
   leading card — it drives the same 0..N range-clear handler as clicking row 8's own
   ordinary card would, while staying put across the re-render. Wait for re-render.
10. Confirm that rows 0 through 8 inclusive are deselected (they no longer
    appear highlighted in the grid). Confirm that rows 9 through 12, which were
    selected and are above row 8, remain selected. Label this **Scenario 4c**.

Expected:
- After Ctrl+clicking a card once, that card's row is deselected and disappears
  from the ordinary-card area; after a second Ctrl+click, it is selected again
  — a toggled round-trip.
- After Shift+clicking a card at row N, the grid shows rows 0 through N
  inclusive as selected; rows above N are deselected.
- After Ctrl+Shift+clicking a card at row N, rows 0 through N inclusive are
  deselected; any rows above N that were previously selected remain selected.

### Scenario 5: Hover sets mouse-over row and fills the mouse-over card

Steps:
1. Ensure Show Mouse Over Row is ON and at least one ordinary card is visible.
   Move the pointer away from all cards first so no card shows a hover state.
2. Move the pointer over one of the selected-row cards (not the
   leading current-row card). Keep the pointer still and wait for the re-render
   (50 ms debounce settle).
3. Confirm that the card now shows a hover highlight and that the grid row
   corresponding to this card is indicated as the mouse-over row.
4. Move the pointer off the card entirely (into a blank area of the screen
   below the Forms viewer panel, off every card and grid row). Wait for
   re-render.
5. Confirm the mouse-over binding is released — the grid shows no mouse-over row
   and the mouse-over card (the card right after the current-row card) is empty.
   The card set is otherwise unchanged: the current-row card stays in place.
   Label this **Scenario 5a**.
6. Reset: move the pointer away from all cards and confirm no hover state is
   active.
7. Move the pointer over a row in the GRID (not the Forms viewer). Keep it
   still.
8. In the Forms viewer, identify the mouse-over position card — the card
   immediately after the current-row card. Read the value of the USUBJID field
   on that card.
9. Confirm that the USUBJID value shown in the mouse-over card matches the
   USUBJID value in the grid for the row currently being hovered.
10. Move the pointer off the grid row. Wait for re-render. Confirm the
    mouse-over card's USUBJID field is now empty (read the field value, not
    just whether the card element is present — the card element stays in place
    but its fields clear when no row is hovered).
    Label this **Scenario 5b**.
    Run steps 7-10 with the mouse exactly as written; the automated run reaches
    the same state by another route, described in the Automation notes.

Expected:
- While the pointer hovers over a card, the grid marks the corresponding row as
  the mouse-over row. When the pointer is moved away to a neutral point, the
  mouse-over binding is released — the grid shows no mouse-over row and the
  mouse-over card empties — while the rest of the card set, including the
  current-row card, stays in place.
- While a grid row is hovered, the mouse-over card (the card immediately after
  the current-row card) fills with that row's field values; the field values go
  from empty (or a prior row's values) to the hovered row's values. Moving the
  pointer off the grid row returns the mouse-over card fields to empty.

### Scenario 6: showMouseOverRow and showCurrentRow toggles

Steps:
1. Ensure both Show Mouse Over Row and Show Current Row are ON. Count the total
   number of card panels currently visible in the Forms viewer.
   Record this count as `baseCount`.
2. In the Forms viewer property panel (Context Panel > Misc), set **Show Mouse
   Over Row** to OFF. Wait for the debounced re-render (50 ms).
3. Count the card panels again. The count must equal `baseCount - 1`.
4. Hover over a grid row. Confirm no card fills with that row's values (the
   mouse-over position no longer exists).
   Label this **Scenario 6a**.
5. Set **Show Mouse Over Row** back to ON. Wait for re-render. Confirm the
   count returns to `baseCount`. Label this **Scenario 6a restore**.
6. In the Forms viewer property panel, set **Show Current Row** to OFF. Wait
   for the debounced re-render.
7. Confirm the current-row card is gone: no card in the area carries the green
   current-row indicator any more.
8. Do NOT judge this by counting cards. With Show Mouse Over Row still ON and
   Show Current Row OFF, this build renders ZERO ordinary cards, so there is no
   "one card fewer" to observe and no "first card" to read — the disappearance
   of the green current-row indicator is the whole signal here.
   Label this **Scenario 6b**.
9. Set **Show Current Row** back to ON to restore defaults.

Expected:
- After turning Show Mouse Over Row OFF and waiting for the debounced
  re-render, the mouse-over position card is absent from the card area; the
  total card count drops by exactly one compared to the pre-toggle count;
  hovering a grid row no longer fills any card with that row's values.
- After turning Show Mouse Over Row back ON and waiting for the debounced
  re-render, the total card count returns to the pre-toggle baseCount.
- After turning Show Current Row OFF and waiting for the debounced re-render,
  the current-row card and its green current-row indicator are absent from the
  card area.

## Automation notes

- **Debounce settle**: all selection, filter, metadata and sort changes are
  debounced by 50 ms in the Forms viewer source
  (DG.debounce(stream, 50), forms-viewer.ts#L196). Every assert taken after
  a selection change, a hover event, or a property toggle needs an explicit
  settle wait.
- **Card counting**: cards are `.d4-multi-form-form` elements. The leading two
  positions (current-row at index 0 when showCurrentRow is on, mouse-over next)
  are built unconditionally even when no row is current or hovered — their
  fields are empty but the card div exists. Subtract the leading offset before
  comparing card count against an expected selected-row count.
- **Gestures made ON THE GRID are driven through the dataframe, by section
  standard — this is a declared actuation substitution, not a workaround.** The
  grid renders on a canvas and exposes no per-row DOM target, so a hover over a
  grid row cannot be actuated as a DOM gesture. Wherever a step says "hover a grid
  row" (Scenario 5b steps 7-10, Scenario 6 step 4), the spec writes
  `df.mouseOverRowIdx` — the exact binding a real grid hover produces — and then
  asserts the DOWNSTREAM observable (the mouse-over card's field values), which is
  a different channel from the actuation. A human running these steps by hand does
  hover with the mouse and sees the same result. Gestures made ON A CARD stay real
  DOM gestures (`.hover()`, `.click({modifiers})`) — the cards are ordinary
  elements — so this substitution applies to the grid side only.
- **Mouse-over card field is the signal, not the element**: the mouse-over card
  element always exists (built unconditionally), so an element-presence or
  element-visibility assert is false-green. Assert the FIELD VALUE inside it:
  the USUBJID INPUT's value goes from empty to the hovered row's USUBJID when
  a grid row is hovered.
- **Pointer-away release is a REAL-gesture assertion — never zero the binding
  from JS** (operator live measurement, demog, 2026-08-11). Moving the real
  pointer to a neutral point off every card and grid row does NOT change the card
  set: position 0 stays the current-row card (green indicator), the mouse-over
  card at position 1 empties, and `df.mouseOverRowIdx` settles to -1 — measured
  stable over six move-away rounds across two neutral points (just below the host
  and the window corner), with no flicker. The historical `Expected -1 /
  Received 5` failure was a SETUP bug, not an unstable assertion: the setup wrote
  `df.mouseOverRowIdx = -1` from JS while the physical pointer still rested on a
  card, and the next re-render re-fired that card's `onmouseenter` and re-bound
  the row under the pointer (5 = the row that card carried). The fix is to move
  the pointer for real, which Scenario 5a now does — hover a card (binding tracks
  its row), then `parkPointerBelowViewer`, then assert `df.mouseOverRowIdx === -1`
  AND the position-1 card is empty. When a JS-driven read (Scenario 5b) needs the
  physical pointer out of the way, `parkPointerBelowViewer` removes it as a
  confounder and the hovered state there is driven by a JS write, not read back on
  the strength of the move.
- **Shift+click selection range is 0..N** (source-derived): the
  onCardMouseDown handler computes
  `df.selection.setAll(false); df.selection.handleClick(row, e)` which sets
  the range from row 0 to the clicked row. Asserting "rows 0..N selected" is
  the correct claim, not "count of rows equals N+1" (which passes by accident
  if some other rows happened to be selected).
- **Ctrl+Shift+click deselects 0..N**: the same handler with ctrl+shift results
  in clearing the 0..row range. Only rows strictly above N that were previously
  selected should remain selected.
- **Scoping**: the Forms viewer's own root carries no `name` attribute; use the
  host `[name="viewer-Forms"]` as the scope and exclude `.temp` — `renderHeader`
  appends a throwaway row-0 form to `document.body` with class `temp`, and
  during that window a global `.d4-multi-form-form` / `[column=…]` query matches
  nodes outside the viewer. Ordinary cards are
  `[name="viewer-Forms"] #vlist .d4-multi-form-form:not(.temp)`.
- **Field elements**: ordinary columns render as INPUT elements carrying
  `column="<ColumnName>"` inside `.d4-multi-form-form`. Use that selector to
  read field values; `INPUT[column="USUBJID"]` scoped to the target card.
- **No automation tokens in step prose**: all `df.*`, `page.*`, and selector
  references above are in this Automation notes section only. Step prose uses
  UI-action language.

### Spec must keep (converged — do not regress)

- **Scoping**: cards addressed as `[name="viewer-Forms"] #vlist
  .d4-multi-form-form:not(.temp)`; the viewer root has no `name` (use the host
  `[name="viewer-Forms"]` as scope).
- **Leading offset 2**: within the ordinary card set the current-row card
  (position 0) and the always-built mouse-over card (position 1) occupy the two
  leading positions, ahead of the selected-row cards; the empty mouse-over card
  contributes `null` field values by position.
- **Readiness barrier must be DISCRIMINATING, never card-count alone.** Two
  different states can produce the same number of cards (e.g. a fresh setup and
  the exit state of a predecessor that changed the selection), so a `count() == N`
  barrier is satisfied instantly by the STALE set before the 50 ms debounced
  re-render builds the fresh one, and the next positional read samples stale rows
  (this caused a Gate B FLAKY: 5a read row 9 from 4c's leftover set). `establishSetup`
  polls the ordinary cards' USUBJID sequence BY POSITION (current row + selected
  rows ascending; the volatile mouse-over position is excluded) against the
  sequence computed live from the dataframe, so the barrier only clears once the
  DOM matches THIS setup's rows — not merely their count.
- **Every read after an action polls for the settled state** (50 ms debounce +,
  for sort, two nested `setTimeout`s); no fixed sleeps as the assert channel.
- **Click/selection channels asserted from the dataframe** (`currentRowIdx`,
  `currentCell`, `currentCol.name`, `selection.get`), a different channel than
  the actuating gesture; `page.evaluate` bodies return primitives / JSON only.
- **Scenario 5a asserts both halves**: real `.hover()` on a card, then poll
  `df.mouseOverRowIdx` equals that card's row (guarded against a stale-state
  false pass by a pointer-free `!= targetRow` read first); then a REAL pointer
  move-away (`parkPointerBelowViewer`) and poll that `df.mouseOverRowIdx` is -1
  AND the position-1 mouse-over card is empty. Never zero the binding from JS in
  this step — it re-binds on the next re-render (see the pointer-away note above).
- **Scenario 5b mouse-over fill** asserts the mouse-over card's USUBJID field
  value (empty → hovered row → empty), never mere card presence; it parks the
  physical pointer off the cards first only to stop it confounding the JS-driven
  `df.mouseOverRowIdx` (a park-then-read whose read is JS-driven, not a
  pointer-away binding assertion).
- **Scenario 6b** asserts the current-row card and its green indicator vanish
  (count 0); the card-count / leading-offset arithmetic is deliberately not
  asserted, because with Show Current Row OFF and Show Mouse Over Row still ON
  this build renders ZERO ordinary cards (measured on demog: 6 → 0), so a
  "one card fewer" count would be a false RED.
