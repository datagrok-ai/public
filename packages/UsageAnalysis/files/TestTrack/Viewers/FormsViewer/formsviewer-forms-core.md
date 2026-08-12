---
feature: formsviewer
realizes_atlas:
  - formsviewer.cp.forms-core
  - formsviewer.int.sort-mirrors-grid
  - formsviewer.int.pinned-rows-persist-by-value
  - formsviewer.int.selection-intersects-filter
  - formsviewer.edge.pin-non-unique-value-warns
  - formsviewer.edge.pinned-row-absent-from-ordinary-cards
realizes:
  - viewers.form
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-20026
    status: fixed
  - id: GROK-19340
    status: fixed
  - id: GROK-20027
    status: fixed
  - id: GROK-20380
    status: open
realized_as:
  - formsviewer-forms-core-spec.ts
expected_results:
  - anchor: "Step 1"
    expectation: >-
      The Forms viewer mounts and its default field set equals the first 20
      visible columns of demog.csv, with no column whose name starts with '~'
      among the header labels, and no warning banner or balloon visible.
  - anchor: "Step 2"
    expectation: >-
      A field element with attribute column="HEIGHT" (or the first visible
      column) inside the current-row card shows a value that matches what the
      grid renders for that cell; after moving the current row to row 77, the
      field shows the value for row 77.
  - anchor: "Step 3"
    expectation: >-
      With showSelectedRows on and three rows selected, exactly three extra
      cards appear in the scrollable card area beyond the two leading cards;
      each extra card's field values identify the selected row it belongs to.
  - anchor: "Step 4"
    expectation: >-
      After applying a filter that excludes one of the selected rows, the
      selected-row cards (the positions after the leading offset) equal
      selection.clone().and(filter) computed at runtime — the filtered-out row
      produces no card, and the remaining selected rows still do.
  - anchor: "Step 5a"
    expectation: >-
      After sorting the grid by HEIGHT, the card order mirrors the grid sort and
      the HEIGHT header label carries a direction indicator (up or down arrow).
  - anchor: "Step 5b"
    expectation: >-
      After setting sortByColumnName to a viewer-specific column, the card order
      follows that column while the grid's own sort is unchanged; the indicator
      appears on the viewer-set column's header label, not on the grid-sort
      column's label.
  - anchor: "Step 5c"
    expectation: >-
      With useGridSort toggled OFF, the grid sort is NOT mirrored — the card
      order stops following the grid sort. GROK-20380 is open, so the product
      still mirrors: the strict assertion of the desired behaviour is kept
      unchanged and wrapped in knownOpenBug('GROK-20380'), which absorbs the
      failure and reports the defect as reproduced. If the product is fixed the
      wrapper throws loudly, which is the signal to unwrap it — the check is
      never softened and never skipped.
  - anchor: "Step 5d"
    expectation: >-
      Three successive double-clicks on the header label of the viewer-sort
      column produce three DIFFERENT indicator states, one of which is "no
      indicator" (the other two being the two arrow directions). The ORDER of
      the three states is not claimed: the incoming sort direction is not a
      controlled or registered setting, so which state the cycle starts from
      depends on the entry state. From the descending entry state, which the
      step establishes deliberately, double-clicking a DIFFERENT column's label
      leaves the indicator on the sort column; this claim is conditional on that
      entry state, because the direction-branched handler clears the sort from
      the ascending state regardless of which label is double-clicked.
  - anchor: "Step 6a"
    expectation: >-
      After right-clicking a selected-row card and choosing Pin Row, that card
      moves into the pinned pane (outside the scrollable area), and the row no
      longer produces an ordinary card even though it is still selected.
  - anchor: "Step 6b"
    expectation: >-
      After right-clicking the pinned card and choosing Unpin Row, the row
      returns to the ordinary card set and the pinned pane is hidden (display:
      none) once the pinned count reaches zero.
  - anchor: "Step 6c"
    expectation: >-
      Pinning a card through a field whose value is NON-unique in its column (an
      earlier row carries the same value) raises a warning balloon with the
      exact text "You have pinned a non-unique value. It won't be applied from
      the layout."; the pin still happens for the session. The pre-existing pin
      through the unique USUBJID field raises no such warning and is left
      intact.
  - anchor: "Step 7a"
    expectation: >-
      Before saving, the field set is moved to a deliberately NON-DEFAULT one —
      non-default in both composition and order — and that is what makes the
      field-set round-trip discriminating rather than a default-vs-default match.
      After saving a layout and re-applying it over a deliberately corrupted
      view (the Forms viewer closed and a foreign Histogram viewer added), the
      Forms viewer is present with its accumulated registered settings restored:
      the field set matches the pre-save non-default field set, and the header
      label that carried the sort direction indicator before saving still carries
      an indicator (the identity of the label, not the direction arrow, is the
      assertion). Applying the saved layout also removes the foreign Histogram
      viewer that the layout does not contain — a layout apply drops a widget
      that is not part of it.
  - anchor: "Step 7b"
    expectation: >-
      The pinned row is present in the pinned pane after the layout round-trip
      and its key field value equals the value recorded before saving; the row
      index is NOT used as the assertion (the pin is restored by value, not by
      index).
  - anchor: "Step 7c"
    expectation: >-
      After a project save (via the ribbon Save button), Close All, and project
      reopen, the same field-set and pinned-row assertions hold across the
      session boundary, and the header label that carried the sort direction
      indicator before saving still carries an indicator (the identity of the
      label, not the direction arrow, is the assertion).
---

# Forms Viewer — Core Ladder (p0)

## Setup

1. Close all open views and dialogs.
2. Open **System:DemoFiles/demog.csv** so the table view is active.
3. Add the **Forms** viewer from **Toolbox > Viewers** (click the Forms icon in
   the viewer toolbox on the right side of the screen).

## Scenarios

### Scenario 1: Default field set on mount

Steps:
1. Inspect the header pane of the Forms viewer (the left column listing field names).
   Confirm the header labels are exactly the first 20 visible (non-`~`) columns of
   demog.csv, in dataframe order — not merely "at most 20" but the same set, in the
   same order (demog has 11 columns, so all 11 appear).
   Confirm that no label starts with the character `~`. Note: demog carries no
   `~`-prefixed columns, so on this dataset this check is degenerate — it is a guard
   against such a column ever leaking into the field set, not live coverage. The
   substantive coverage of `~`-column exclusion lives in the fields-lifecycle scenario
   of this section (which renames a column to a `~`-prefixed name); do not read this
   step as covering that edge.
   Confirm that no warning banner or tooltip balloon is visible anywhere in the
   Forms viewer panel.

Expected:
- The Forms viewer mounts and its default field set equals the first 20 visible
  columns of demog.csv, with no column whose name starts with '~' among the
  header labels, and no warning banner or balloon visible.

### Scenario 2: Current-row field value and row tracking

Steps:
1. Note which row is current in the grid (highlighted in green) and read the value
   the grid shows for the HEIGHT column on that row.
2. In the Forms viewer, read the value displayed in the HEIGHT field of the card
   that carries the green current-row indicator.
3. In the grid, click row 77 to make it current (if row 77 is already current,
   click any other row first, then row 77).
4. After the Forms viewer re-renders (allow for the 50 ms debounce settle), read
   the HEIGHT field of the current-row card again.

Expected:
- A field element with attribute column="HEIGHT" (or the first visible column)
  inside the current-row card shows a value that matches what the grid renders
  for that cell; after moving the current row to row 77, the field shows the
  value for row 77.

### Scenario 3: Show Selected Rows — cards per selected row

Steps:
1. Verify that **Show Selected Rows** is on in the Forms viewer's property panel
   (Context Panel > Misc > Show Selected Rows; it ships enabled, but confirm
   rather than toggle blindly).
2. In the grid, Ctrl+click three distinct rows to select them (choose rows that
   are all visible and unfiltered).
3. Wait for the Forms viewer to re-render (50 ms debounce settle plus grid
   interaction lag).
4. Count the cards in the scrollable card area that are beyond the two leading
   cards (current-row position and mouse-over position).
5. For each extra card, read one field value and match it against the row the
   grid reports as selected.

Expected:
- With showSelectedRows on and three rows selected, exactly three extra cards
  appear in the scrollable card area beyond the two leading cards; each extra
  card's field values identify the selected row it belongs to.

### Scenario 4: Selection intersects filter — filtered-out rows get no card

Steps:
1. Starting from the state after Scenario 3 (three rows selected), open the
   filter panel (Filter icon in the toolbar or Shift+Ctrl+F).
2. Apply a filter on the SEX column that excludes at least one of the three
   selected rows (for example, if you selected rows with both M and F values,
   filter to show only one sex).
3. Wait for the Forms viewer to re-render after the filter change (50 ms debounce).
4. Compute the expected card set: starting from the three selected rows, remove
   those that the active filter excludes. Confirm the number of selected-row
   cards (the positions after the leading offset) matches this intersection.
5. Confirm that a row that is selected but filtered out produces no card among
   the selected-row cards.
6. Remove the filter (clear the filter panel) to restore the full dataset for
   subsequent steps.

Expected:
- After applying a filter that excludes one of the selected rows, the
  selected-row cards (the positions after the leading offset) equal
  selection.clone().and(filter) computed at runtime — the filtered-out row
  produces no card, and the remaining selected rows still do.

### Scenario 5: Sort — grid sort mirroring, viewer-specific sort, useGridSort toggle, double-click cycle

Steps:
1. In the grid, click the HEIGHT column header to sort by HEIGHT ascending. Wait
   for the Forms viewer to re-render.
2. Read the order of the first three selected-row cards and confirm it matches the
   sorted grid order (lowest HEIGHT first). Confirm the HEIGHT header label in
   the Forms viewer carries a direction indicator (up or down arrow beside the
   label text). Label this Step 5a.
3. In the Forms viewer property panel (Context Panel > Data), set **Sort By** to
   **WEIGHT** (the registered property name is sortByColumnName; address it as
   such in any property-panel interaction). Wait for the re-render.
4. Confirm the card order now follows WEIGHT values (the row with the lowest or
   highest WEIGHT value appears first depending on sortAscending state), while
   the HEIGHT direction indicator in the Forms viewer header is gone and the
   WEIGHT label now carries the indicator. Confirm the grid's own HEIGHT sort
   indicator is unchanged. Label this Step 5b.
5. Clear **Sort By** so the viewer has no sort of its own and the grid's HEIGHT
   sort is the only ordering authority. Wait until the card order actually follows
   the grid's HEIGHT sort. Now set **Use Grid Sort** to OFF and wait for the
   re-render. Read the card order again: it must NO LONGER match the grid's HEIGHT
   sort order. IMPORTANT: if the product continues to mirror the grid sort after
   Use Grid Sort is turned off, this step must FAIL as a hard assertion failure —
   do not weaken the check or skip it. This is the deliberate known-red state for
   GROK-20380. Label this Step 5c.
6. Reset: re-enable **Use Grid Sort** and make sure **Sort By** is still empty.
   Now set **Sort By** to AGE and confirm the AGE direction indicator appears. Now double-click the AGE
   header label in the Forms viewer. Record the indicator state after the
   double-click. Double-click again and record. Double-click a third time and
   record. The three recorded states must all be different and one of them must be
   "no indicator"; do not expect a fixed order, because the direction the cycle
   starts from depends on the direction the sort was already in, which is not a
   setting you can control.
   Also: the claim "double-clicking a DIFFERENT column's label leaves the indicator
   on the sort column" holds ONLY from the descending entry state, and the step must
   establish that state on purpose. The viewer's double-click handler branches only on
   whether a sort column is set and on its direction — it does NOT look at which label
   was clicked; from the ascending state a double-click on ANY label (even a different
   column) CLEARS the sort and the indicator falls onto the grid-sort column. So first
   drive the sort to the descending state deterministically (clear Sort By, then a
   single double-click on AGE, which from a cleared state always lands descending),
   confirm the indicator is on AGE, then double-click a DIFFERENT column header (e.g.
   HEIGHT) and confirm the indicator stays on the AGE label — from descending the
   handler only flips AGE to ascending, it never reassigns the column to the clicked
   label nor clears it.
   Label this Step 5d.

Expected:
- After sorting the grid by HEIGHT, the card order mirrors the grid sort and the
  HEIGHT header label carries a direction indicator (up or down arrow).
- After setting sortByColumnName to a viewer-specific column, the card order
  follows that column while the grid's own sort is unchanged; the indicator
  appears on the viewer-set column's header label, not on the grid-sort column's label.
- With useGridSort toggled OFF, the grid sort is NOT mirrored — the card order
  stops following the grid sort. GROK-20380 is open, so the product still
  mirrors: the strict assertion is kept unchanged and wrapped in
  knownOpenBug('GROK-20380'), which absorbs the failure and reports the defect as
  reproduced; a fix makes the wrapper throw loudly, which is the signal to unwrap
  it. The check is never softened and never skipped.
- Three successive double-clicks on the header label of the viewer-sort column
  produce three DIFFERENT indicator states, one of which is "no indicator". The
  order of the three states is not claimed — the incoming sort direction is not a
  controlled setting, so the cycle's starting point depends on the entry state.
  From the descending entry state, which the step establishes deliberately,
  double-clicking a different column's label leaves the indicator on the sort column;
  this claim is conditional on that entry state, because the direction-branched
  handler clears the sort from the ascending state regardless of which label is
  double-clicked.

### Scenario 6: Row pinning — pin, verify exclusion from ordinary set, unpin round-trip

Steps:
1. Ensure at least one row is selected so ordinary cards exist. Clear any active
   filter. Make sure Show Selected Rows is on.
2. Right-click the **USUBJID** field of one of the ordinary (selected-row) cards —
   not the leading current-row or mouse-over card. Right-click that field, not the
   card's empty area: the pin remembers the field the right-click landed on, and
   USUBJID is the field whose value is unique. The context menu must offer the item
   "Pin Row". Click "Pin Row".
3. Wait for the Forms viewer to re-render (50 ms debounce).
4. Confirm the card appears in the pinned pane: a container outside the scrollable
   card area that carries a "Pinned row" indicator on the card. Note the value of
   the **USUBJID** field — call this the "anchor value". The anchor field must be
   one whose value is UNIQUE in the table (USUBJID is; HEIGHT is not): the pin is
   remembered by value, so a repeated value would resolve back to the first row
   that carries it, not to the row you pinned.
5. Confirm the row no longer produces an ordinary card: count the ordinary cards
   and compare with the pre-pin count; the pinned row's card is absent from the
   ordinary set even though the row is still selected. (Turn Show Mouse Over Row
   OFF or move the pointer away from the pinned card before counting, to avoid the
   mouse-over position adding it back to the visible count.)
   Label this Step 6a.
6. Right-click the pinned card. The context menu must offer "Unpin Row". Click it.
7. Wait for re-render. Confirm the row produces an ordinary card again (its card
   reappears among the selected-row cards). Confirm the pinned pane is no longer
   visible (its display collapses to none once the pinned count reaches zero).
   Label this Step 6b.
8. Re-pin the same row (repeat steps 2-4) so the peak in Scenario 7 persists a
   pinned row.
9. Additionally pin a DIFFERENT ordinary card by right-clicking its SEX field — a
   value that repeats earlier in the column (not the first M or first F). A warning
   balloon must appear with the exact text "You have pinned a non-unique value. It
   won't be applied from the layout." The pin still happens for the session. Then
   Unpin that SEX-pinned row again, leaving only the unique USUBJID pin from the
   earlier step for Scenario 7.

Expected:
- After right-clicking a selected-row card and choosing Pin Row, that card moves
  into the pinned pane (outside the scrollable area), and the row no longer
  produces an ordinary card even though it is still selected.
- After right-clicking the pinned card and choosing Unpin Row, the row returns
  to the ordinary card set and the pinned pane is hidden (display: none) once
  the pinned count reaches zero.
- Step 6c: pinning through a NON-unique field value raises the warning balloon
  with that exact text and still pins for the session; the pre-existing unique
  USUBJID pin raises no warning and is left intact.

### Scenario 7: Persistence peak — layout and project round-trip

Steps:
1. First move the field set OFF its default: pick a NON-DEFAULT set that differs
   from the default both in composition (a proper subset — but keep USUBJID, the
   pinned-row anchor field, and the field that carries the sort indicator) and in
   order (reorder it so the drawn sequence is not the default column order). This
   is not decoration: the field-set round-trip checked in step 5 is only
   discriminating because the saved set is non-default — a default-vs-default
   comparison would pass even on a build that ignored the saved field set and
   rebuilt the default one. Then record the current accumulated state: note the
   list of field names shown in the header pane, note which header label carries
   the sort direction indicator (if any), and note the anchor value of the pinned
   row — the **USUBJID** value read in Scenario 6, step 4. Use USUBJID and not a
   field like HEIGHT: the pin is restored by value, so a non-unique value would
   resolve back to the first row carrying it instead of the pinned row.
2. Save the layout via **View | Layout | Save to Gallery** (keyboard shortcut
   Ctrl+S). The save is SILENT — no dialog appears and no confirmation balloon is
   shown; the layout is named automatically from the view name. Do NOT wait for a
   dialog that will not come. Confirm the layout now exists — any route that
   identifies the saved layout is acceptable; the gallery is one such route, not
   a requirement.
3. Deliberately edit the current layout: close the Forms viewer (click Close on
   its title bar) and add a different viewer (e.g. Histogram from the Toolbox).
   The table view should now show the Histogram but NOT the Forms viewer.
4. Re-apply the layout saved in step 2 to the current view. What matters is that
   the SAME saved layout is applied — the gallery entry is one way to reach it,
   not the thing under test.
5. After re-apply, confirm the Forms viewer is present (the Histogram added in
   step 3 is gone). Confirm the field set matches the recorded field names from
   step 1. Confirm the header label that carried the sort indicator in step 1
   still carries an indicator (assert the IDENTITY of the label, not the direction
   of the arrow — the direction is out of scope per the known GROK-20666 defect
   and must NOT be asserted here). Confirm the pinned pane is visible and the
   pinned card's anchor field shows the same anchor value recorded in step 1.
   Label the field-set and sort-indicator half **Step 7a** and the pinned-row
   half **Step 7b**.
6. Clean up: delete the probe layout saved in step 2, even if the step above
   failed.
7. Return to the table view with the restored layout still active. Save the
   entire view as a project: click the ribbon **Save** button (the floppy-disk
   icon in the main toolbar — this is the Save that captures the view layout;
   do not use File > Export). Enter a distinct project name (e.g. "FormsViewer
   Core Test Project"). Save.
8. Close all views: close the table view and any other open views.
9. Reopen the saved project from the Projects browser (Browse > Projects, find
   the entry, double-click).
10. After the project opens, confirm the same assertions hold across the session
    boundary: the Forms viewer is present with the recorded field set, the sort
    column's header label carries an indicator (identity, not direction), and the
    pinned card's anchor field value equals the pre-save value.
    Label this Step 7c.
11. In a finally block (runs even on failure): delete the probe project from the
    Projects browser.

Expected:
- Before saving, the field set is moved to a deliberately NON-DEFAULT one (a
  reordered proper subset), which is what makes the field-set round-trip
  discriminating rather than a default-vs-default match. After saving a layout via
  View | Layout | Save to Gallery and re-applying it, the Forms viewer is present
  with its accumulated registered settings restored: the field set matches the
  pre-save non-default field set, and the header label that carried the sort
  direction indicator before saving still carries an indicator (the identity of the
  label, not the direction arrow, is the assertion). Applying the saved layout also
  removes the foreign Histogram viewer that the layout does not contain — a layout
  apply drops a widget that is not part of it.
- The pinned row is present in the pinned pane after the layout round-trip and
  its key field value equals the value recorded before saving; the row index is
  NOT used as the assertion (the pin is restored by value, not by index).
- After a project save (via the ribbon Save button), Close All, and project reopen,
  the same field-set and pinned-row assertions hold across the session boundary, and
  the header label that carried the sort direction indicator before saving still
  carries an indicator (the identity of the label, not the direction arrow, is the
  assertion).

## Automation notes

- **Debounce settle**: all selection, filter, metadata and sort changes are debounced
  by 50 ms in the Forms viewer source (DG.debounce(stream, 50), forms-viewer.ts#L196).
  The grid-sort path has an additional setTimeout inside the subscription handler
  (forms-viewer.ts#L217-224). Every assert taken after one of these changes needs
  an explicit settle wait beyond the ordinary one.
- **Card counting**: cards are `.d4-multi-form-form` elements. The leading two
  positions (current-row at index 0 when showCurrentRow is on, mouse-over next)
  are built unconditionally even when no row is current or hovered — their fields
  are empty but the card div exists. Subtract the leading offset before comparing
  card count against a selected-row expectation.
- **Card count is a debounce-settle, not a virtual-view clip — for SMALL selections**
  (live-diagnosed 2026-08-11). With 3 rows selected the ordinary count settles
  `2 → 5` within one debounce cycle (≈250–500 ms); polling with a generous
  timeout is sufficient and reliable. The alternative worry — that
  `ui.virtualView` materialises only the on-screen cards, making a DOM count
  unreliable — does NOT bite here: all five cards were present in a 362 px-tall
  `#vlist` even though each card is ≈300 px, so `virtualView` renders the whole
  small set regardless of fit. This holds only because the probe keeps the
  selection tiny (3 rows). A large selection WOULD be clipped and a DOM count of
  it would be unreliable — keep the selected set to a handful.
- **Pinned pane**: `.d4-multi-form-pinned-forms`. Its display collapses to `none`
  when empty; check `display !== 'none'` rather than element presence.
- **Sort indicator**: the direction indicator lives on the header label element as
  a child with class `.d4-multi-form-column-sort-indicator`. Assert its presence
  on the expected label and absence on all others; do NOT assert the arrow text
  direction (↑ / ↓) for the persistence check (GROK-20666 exception).
- **GROK-20380 is guarded by `knownOpenBug`**: the useGridSort OFF assertion in
  Scenario 5c keeps the REAL strict `expect` on the desired behaviour, wrapped in
  `knownOpenBug('GROK-20380', ...)`. This is not a softening — the same hard assert
  runs every time; the wrapper only makes it self-flipping (an expected failing
  reproduction is logged and green while the defect is open, and an unexpected pass
  throws loudly once it is fixed). Do NOT weaken the assertion inside the wrapper,
  and do NOT remove the wrapper while GROK-20380 is open.
- **Layout save is silent**: View | Layout | Save to Gallery (Ctrl+S) produces no
  dialog and no confirmation balloon. The step verifies the save by locating the
  saved layout itself — never by waiting for a dialog.
- **Sort direction is OUT OF SCOPE for persistence**: GROK-20666 (open) means the
  direction cannot be correctly asserted after a round-trip. Assert only the
  identity (which label carries the indicator), not the direction. This exclusion
  is deliberate and permanent until GROK-20666 is fixed.
- **Read a card only after its row is current** (root cause of the original
  Gate B red, live-confirmed 2026-08-11): `buildForm(row = -1)` returns one empty
  `div` per field with NO `column` attribute, so the current-row card (position 0)
  has zero `[column="…"]` inputs until a row is made current. A freshly opened
  table can sit at `currentRowIdx = -1`, and `grid.cell('HEIGHT', -1).cell.valueString`
  is `""` — together that produces the exact `Expected: "" Received: null`
  signature. Step 2 therefore sets `currentRowIdx` FIRST, then polls the card;
  no wait alone can populate a card whose row is -1.
- **Pin records the column under the right-click** (live-confirmed 2026-08-11):
  the pinned pair is `(e.target.closest('[column]'), stringified value)`. A
  right-click on the card centre lands on the field at the card's vertical
  midpoint — HEIGHT on demog's 11-field card, a NON-UNIQUE float — so
  `pinnedRowValues` records e.g. `"175.35"` and the by-value restore resolves to
  the wrong row (or none). The spec right-clicks the **USUBJID** field
  (`cardContextMenu(..., 'USUBJID')`) so the recorded key is unique, matches the
  anchor read from that same field, and survives the layout round-trip. Do not
  revert to a centre right-click.
- **Scoping**: the Forms viewer's own root carries no `name` attribute; use the
  host `[name="viewer-Forms"]` as the scope and exclude `.temp` — `renderHeader`
  appends a throwaway row-0 form to `document.body` with class `temp`, and
  during that window a global `.d4-multi-form-form` / `[column=…]` query matches
  nodes outside the viewer. Ordinary cards are
  `[name="viewer-Forms"] #vlist .d4-multi-form-form:not(.temp)`.

- **Declared actuation substitutions (section standard, not workarounds).** Several
  prose-prescribed gestures are driven through the dataframe / JS API rather than the
  literal UI control, because the target is canvas-drawn (no per-row DOM handle) or not
  drivable headless. Each still asserts the DOWNSTREAM observable — a channel distinct
  from the actuation — and a human running the steps by hand reaches the same state:
  - "click the HEIGHT column header to sort" → `grid.sort(['HEIGHT'], [true])` (Step 5a).
  - "Ctrl+click three rows to select" → `df.selection.set(...)` on the grid's canvas
    rows, which expose no per-row DOM target (Steps 3, 6a, 6b, 6c).
  - "click row 77 to make it current" → `df.currentRowIdx = 77` (Step 2).
  - "open the Filter panel and apply / clear a filter" → the categorical filter is
    applied and reset through the filter group (Step 4, `v.applyCategoricalFilter` /
    `v.resetFilters`).
  - "close the viewer from its title bar" and "add Histogram from the Toolbox" →
    `forms.close()` and `tv.addViewer('Histogram')` (Step 7a corruption setup).
  - "reopen the project from the Projects browser" → `project.open()` by id (Step 7c,
    `projects.reopenAndAssertProvenance`).
  What stays a REAL UI gesture: the layout Save (`View | Layout | Save to Gallery` via
  `v.driveTopMenuLeaf`, Step 7a), the project Save (ribbon Save button via
  `projects.saveProjectViaUI`, Step 7c), the Pin / Unpin context menus, and the Use Grid
  Sort and Show Mouse Over Row toggles (property panel via `v.setPropertyGridCheckbox`).
  The Sort By column substitution is declared separately below.

### Spec must keep (survive future re-authorings)

- **Step 2**: set `currentRowIdx` before reading the current-row card; poll the
  card field, never read once. Reason above.
- **Pin actuation targets the USUBJID field**, not the card centre (Steps 6a, 6b
  re-pin). Reason above. The anchor is read from USUBJID and must equal the
  recorded `pinnedRowValues` entry.
- **Card/field/label queries are scoped to `[name="viewer-Forms"]` and exclude
  `.temp`** (the header builds a throwaway `.temp` form on `document.body`).
- **Leading offset is positional**: within the ordinary card set the current-row
  card (position 0 when showCurrentRow, identified by its green indicator — the
  `CURRENT` selector) and the mouse-over card (next) occupy the two leading
  positions, ahead of the selected-row cards, and are NOT filtered against the
  pinned set. Slice by the offset FIRST, then drop empties.
- **Expected values are computed at runtime**, not from `col.getString`: the current
  card value is compared to `grid.cell(col, row).cell.valueString`; the ordinary
  set to `selection.clone().and(filter)`; the sorted order to the grid/viewer
  sort. Keep those reference channels.
- **Field set after a round-trip is verified by the DRAWN LABELS, not the property.**
  The baseline (`drawnLabelNames`, the sequence of `.d4-multi-form-column-name`
  texts) is captured BEFORE the save and compared after restore in BOTH Step 7a
  (layout) and Step 7c (project). Reading `vw.props.fieldsColumnNames` on both sides
  is prop-echo false green: the refdoc documents a live failure it cannot see (10 of
  11 labels drawn while the property still listed 11). The property MAY be checked
  additionally, but NEVER instead of the drawn labels.
- **Step 1 is a real set equality**, not `length > 0 && <= 20`: the drawn label
  sequence must `toEqual` the first 20 visible (non-`~`) dataframe columns computed
  on the fly. And the "silent cap" check counts balloons of ANY severity
  (`.d4-balloon`), not just `.d4-balloon.error` — a warning cap notice would be
  `.d4-balloon.warning`.
- **Step 4 recomputes selection∩filter INSIDE the poll**, every iteration, together
  with the DOM tail in one `evaluate`, and REQUIRES that the filter actually excluded
  a selected row (`inter.length < selectedCount` — an invariant, NOT a hardcoded
  count). Capturing the expected set once right after the filter write races the 50 ms
  re-render and can degrade the check to "cards == selection" without proving the
  filtered-out row drops.
- **Step 5a compares the card order to `df.getSortedOrder([col],[asc], selection∩filter)`**
  computed at runtime, read in the same `evaluate` as the DOM tail — not monotonicity
  of three HEIGHT values (passes by luck ~1/6 on three rows).
- **Step 5b compares the card order to `df.getSortedOrder(['WEIGHT'],[asc], selection∩filter)`**
  the same way (direction `asc` taken from the live WEIGHT indicator arrow, `nonEmpty`
  vacuum guard) — NOT monotonicity of three WEIGHT values. One standard for order
  assertions across the file (5a, 5b, 5c all use `getSortedOrder`). Its two other
  assertions — indicator sits on the WEIGHT label, and the grid's own sort stays on
  HEIGHT — are strict and prove the sort-authority handover; keep them verbatim.
- **Step 6a gates on ORDER stability** (`waitForOrderStable` — two equal consecutive
  USUBJID-sequence reads) before reading the anchor, because the card COUNT settles
  before the ORDER after a sort clear; a count-only gate reads the anchor off a
  pre-reorder card while the right-click hits the post-reorder card at the same index
  (the Gate B flake, byte-identical 2/3). It also re-confirms the pinned row is STILL
  selected after pinning.
- **Step 3 confirms the `showSelectedRows` default by reading it BEFORE any write** —
  setting it true then reading it back is a tautology.
- **The project part of the persistence peak (Step 7c) goes ONLY through
  `saveProjectViaUI`** (the real ribbon Save button) — only the UI save captures
  the view layout (viewers + their settings). Substituting `DG.Project.create()`
  + `addChild` + `grok.dapi.projects.save` persists an empty project that reopens
  as a bare Grid, so the step passes without proving anything: that is false
  green, not a persistence check. Reopen via `reopenAndAssertProvenance` and clean
  up via `deleteProjectWithCleanup` in `finally`.
- **Every `page.evaluate` returns primitives or JSON** (avoids the "object
  reference chain is too long" serialization failure) — never return a DG object.
- **Sort By is a DECLARED ACTUATION SUBSTITUTION.** The steps prescribe setting it in
  the property panel, but every spec sets it as `setOptions({sortByColumnName})`
  (`props.sortBy` throws — that is the registered name). The panel's own column
  selector `[name="div-column-combobox-sort--by"]` (`.d4-column-selector`) is NOT
  driven by any spec in this section, so the refdoc flow `task-set-the-viewers-sort-column`
  stays unexercised. This is an open declared debt recorded in the section chain's
  `gate_f_verdict`, not a claim that the panel route works.
- **Probe cleanup (layout, project) runs in a `finally`** so it happens even when
  an assertion above it fails.
- **Step 5c (useGridSort OFF, GROK-20380) is wrapped in
  `knownOpenBug('GROK-20380', () => { expect(stillMirrorsGridSort).toBe(false); })`**
  with the REAL strict assertion inside. `stillMirrorsGridSort` is the exact-sequence
  equality of the card order with `df.getSortedOrder(grid.sortByColumns, grid.sortTypes,
  selection∩filter)` — NOT monotonicity of three HEIGHT values, which could pass by
  accident on the (fixed-state) unsorted set and then the wrapper would NOT throw on
  fix-day, silently losing the self-flip signal exactly when it is the only thing that
  matters. The wrapper is not a softening: the same hard `expect` on the desired
  behaviour (Use Grid Sort off must stop mirroring the grid sort) runs every time.
  While the defect is open the reproduction fails and is logged green; the day it is
  fixed the wrapper throws loudly (`[KNOWN_BUG_FIXED:GROK-20380]`). Do NOT remove the
  wrapper and do NOT weaken the assertion inside it. When the defect is fixed, replace
  the wrapper with a plain hard `expect` and set `related_bugs[].status: fixed` (+ `fixed_in`).
- **Persistence asserts sort-indicator IDENTITY, not arrow direction** (GROK-20666
  keeps the direction out of scope).
