---
feature: formsviewer
realizes_atlas:
  - formsviewer.cp.fields-lifecycle-and-number-format
  - formsviewer.int.number-format-vs-grid
  - formsviewer.edge.empty-fields-column-names
  - formsviewer.edge.tilde-columns-excluded
  - formsviewer.edge.over-20-columns-capped-silently
  - formsviewer.edge.renamed-column-follows-rename
  - formsviewer.edge.number-format-float-only
realizes:
  - viewers.form
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-20367
    status: fixed
  - id: GROK-20027
    status: fixed
  - id: GROK-14962
    status: fixed
realized_as:
  - formsviewer-fields-lifecycle-and-number-format-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      Header labels appear in exactly the order the columns were picked in the
      Fields dialog (RACE, AGE, SEX) — not the table's column order — and the
      field elements inside each card follow that same sequence. No extra labels
      or fields appear.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      After removing one field via its close icon, that field's label and every
      [column="<removed>"] element are absent from all cards; the remaining
      labels and fields retain their original sequence without gaps.
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After the column is dropped from the dataframe, no label or field element
      for the removed column name appears in the form; the remaining fields keep
      their order and no console error is raised.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      The renamed column keeps rendering; its header label shows the new column
      name immediately (or after the next re-render — confirmed by live recon);
      no console error is raised. The field element carries the new column name
      in its column attribute.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      With a table carrying more than 20 visible columns, the Forms viewer
      renders exactly 20 field labels in the header pane and 20 field elements in
      the current-row card. No warning balloon, no warning banner, and no message
      of any kind appears — the cap is silent.
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      No label and no field element whose name starts with '~' appears in the
      viewer, either in the initial auto-selected set or after renaming a column
      to a '~'-prefixed name.
  - anchor: "Scenario 1 Step 8b"
    expectation: >-
      With every field removed (an empty field set), the viewer draws zero
      header labels and zero [column] field elements, raises no console error,
      and shows no error balloon — the empty-list build is a silent no-op.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      For a FLOAT column with no explicit format tag, the text rendered in the
      form field for the current row equals the grid column's declared display
      format applied to the raw value (DG.format(raw, gridColumn.format)) — the
      reference is computed independently of the viewer's own col.getString(row)
      path, not read back from it.
  - anchor: "Scenario 2 Step 6"
    expectation: >-
      After switching numberFormat to an explicit mask ('3 digits after comma'),
      the same FLOAT field's rendered text changes to the mask's shape (exactly
      three decimal places) and now differs from the grid cell's text. The
      underlying dataframe value is unchanged.
  - anchor: "Scenario 2 Step 7"
    expectation: >-
      For an INT column and a string column, the rendered field text is identical
      before and after the numberFormat change to an explicit mask; named number
      formats apply to FLOAT columns only. demog carries no QNUM-typed
      (qualified-number) column, so the integer and string channels are the two
      distinct non-FLOAT channels this check can exercise.
---

# Forms Viewer — Field Lifecycle and Number Format

## Setup

1. Open the demo file `System:DemoFiles/demog.csv`.
2. Add a computed FLOAT column with no explicit format tag — in the console or
   via the Functions panel, compute `HEIGHT * 1.0` and name the result `COMPUTED_H`
   (this gives a FLOAT column that carries no format tag, which is the entry state
   required for the GROK-20367 regression guard).
3. Add the Forms viewer to the table view from the Viewers toolbox.
4. Set the current row in the grid to row 0 (click the first data row) so the
   current-row card is populated.

## Scenarios

### Scenario 1: Field set lifecycle — order, removal, column changes, caps

Steps:

1. Open the Forms viewer's property panel (right-click the viewer header or use
   the gear icon) and locate the **Fields** property. Click its "..." button to
   open the Select Columns dialog.
2. Clear the current selection entirely, then pick exactly three columns in the
   following non-table order: **RACE**, then **AGE**, then **SEX**. Confirm with
   OK.
3. Read the sequence of header labels in the left pane of the Forms viewer and
   the sequence of field elements inside the current-row card.
4. In the viewer's header pane, click the close icon on the **AGE** label to
   remove that field from the field set. Wait for the viewer to re-render.
5. From the grid's top toolbar or via the column header context menu, drop the
   **RACE** column from the dataframe itself. Wait for the onColumnsRemoved
   re-render (debounced 50 ms).
6. Right-click the **SEX** column header in the grid and choose the rename
   option; rename the column to **GENDER** using the rename dialog. Wait for the
   re-render. (The rename MUST be performed through the grid's column-header
   context menu with a real mouse gesture — assigning `col.name` through the JS
   API is forbidden here and produces different behaviour.)
7. Open any table with MORE THAN TWENTY visible columns, add the Forms viewer to
   it, and read the field count in the left header pane. Any such table will do —
   the claim is about the column-count cap, not about a particular dataset. If no
   wide table is at hand, make one: open demog.csv in a new view and add computed
   columns to it (e.g. `HEIGHT * 1.0` repeated under names `C01`, `C02`, …) until
   more than twenty visible columns exist.
8. In the table from Step 7, rename any visible column to start with a `~`
   prefix (e.g. `~SERVICE`) through the grid's column-header context menu. Wait
   for the re-render, then read the header labels.
9. Remove every field from the field set (open the Fields dialog and click None, or
   clear the field list). Wait for the re-render. Label this Step 8b.

Expected:

- Scenario 1 Step 3: header labels appear in the order RACE, AGE, SEX (the picked order,
  not the table's column order — GROK-14962); field elements inside each card
  follow that same sequence; no extra labels or fields appear.
- Scenario 1 Step 4: the AGE label and every `[column="AGE"]` element are absent from all
  cards; RACE and SEX retain their positions without gaps.
- Scenario 1 Step 5: no label and no field element for "RACE" appears; the remaining field
  (GENDER, still named SEX at this point) keeps rendering; no console error is
  raised.
- Scenario 1 Step 6: the field now labelled GENDER continues to render; its header label
  shows "GENDER" (the new column name); no console error is raised.
- Scenario 1 Step 7: exactly 20 field labels appear in the header pane; no warning balloon,
  no warning banner, and no message of any kind is shown — the 20-column cap is
  silent (the source still declares a COLS_LIMIT_EXCEEDED_WARNING constant, but
  nothing reads it: the message was removed deliberately).
- Scenario 1 Step 8: no label and no `[column]` element whose name starts with `~` appears
  in the viewer; the renamed column is absent from the field set (GROK-20027).
- Scenario 1 Step 8b: with an empty field set the viewer draws zero header labels and zero
  [column] value elements, raises no console error, and shows no error balloon —
  the empty-list build is a silent no-op.

### Scenario 2: Number format — regression guard and float-only scope (GROK-20367)

Steps:

1. Open demog.csv again and re-create the COMPUTED_H computed column — the
   previous step closed the environment, so the table is opened fresh and the
   column rebuilt. Add the Forms viewer, confirm Number Format reads its default
   "Same as grid" without changing it, and set the current row to row 0.
2. Read the text of the `[column="COMPUTED_H"]` field element inside the
   current-row card. This is the **form field text**.
3. Compute the reference text for COMPUTED_H at row 0 as the grid column's declared
   display format applied to the raw value; do not read it from the drawn grid cell. This is the **grid cell text**.
4. Compare the two texts.
5. In the Forms viewer's property panel, change **Number Format** from "Same as
   grid" to "3 digits after comma". Wait for the re-render. (Three digits, not
   two: the grid already renders this no-format-tag column to 2 decimals, so a
   2-digit mask would COINCIDE with the grid's own rendering and the "differs
   from the grid" half of the expectation could never hold. Do not "simplify"
   this back to 2 digits.)
6. Read the form field text for COMPUTED_H again and compare it with the grid
   cell text for the same row.
7. Now read the form field text for the **AGE** column (an INT column) before
   and after the format change, and the **SEX** column (a string column) before
   and after. These are two distinct non-FLOAT channels — integer and string.
   demog carries no QNUM-typed (qualified-number) column, which would be a third
   channel; the two present channels already show that the named format leaves
   non-FLOAT fields untouched.

Expected:

- Scenario 2 Step 4: the form field text for COMPUTED_H at row 0 equals the grid column's
  declared display format applied to the raw value (GROK-20367 guard — the
  reference is DG.format(raw, gridColumn.format), computed independently of the
  viewer's own col.getString(row) path, not read back from it).
- Scenario 2 Step 6: the form field text for COMPUTED_H now shows exactly three decimal
  places (the explicit mask shape) and differs from the grid cell text; the
  underlying dataframe value is unchanged.
- Scenario 2 Step 7: the INT and string fields render identically before and
  after the format change — named numberFormat applies only to FLOAT columns.
  demog has no QNUM-typed (qualified-number) column, so integer and string are
  the two distinct non-FLOAT channels this check can exercise.

## Automation notes

- **Header close icon is visibility:hidden, revealed only by a per-node inline
  handler**: columnLabelContainer.onmouseenter sets closeIcon.style.visibility =
  'visible' (forms-viewer.ts renderHeader). renderHeader() starts with
  ui.empty(this.columnHeadersDiv) and rebuilds every container and icon from
  scratch, so any re-render landing after a hover discards the inline visibility.
- **Column rename / remove actuate ONLY through the grid column-header UI with a
  real right-click**; assigning col.name from JS raises semType-null errors and
  loses the field (actuation-path artifact, not a viewer defect).
- **Removing a FIELD via the header close icon** is the viewer's own affordance
  (it filters fieldsColumnNames), distinct from dropping a COLUMN from the dataframe.
- **50 ms debounce on all viewer subscriptions**: after any column-set change,
  poll for the settled state, never sleep a fixed time.
- **Number-format reference is the grid column's declared format applied to the raw value** (DG.format(raw, gridColumn.format)), never
  col.getString(row) — the viewer itself goes through col.getString, so comparing
  against it stays green on a broken build.
- **Number Format is changed through the property-panel choice editor**: its VALUE
  element (not the row) creates the `<select>`, which stays in edit mode after the
  first use, so a later change reuses the open select. The row lives under Misc.

### Spec must keep (survive future re-authorings)

- **Header close-icon removal must use a REAL-mouse bundle, rebuild-proof.** Re-query
  the label container each iteration (renderHeader() may have rebuilt it), then
  page.mouse.move onto the label centre (fires onmouseenter → reveals the icon) and
  immediately page.mouse.click the icon centre — NO intervening wait and NO toBeVisible.
  Retry the whole move→click bundle until the field leaves fieldsColumnNames (the source
  of truth; a transient empty header mid-rebuild must not report a false "removed"). The
  move→click window is ~ms vs the viewer's 50 ms subscription debounce, so it usually
  wins on the first try; the retry only covers a rebuild landing in that window. A
  hover→toBeVisible→click with the visibility assert as a separate step straddled the
  rebuild and went Gate B FLAKY — that is a lost-state race, NOT a timeout, so do not
  "fix" it by lengthening waits. Do NOT substitute a synthetic dispatchEvent(mouseenter)
  + el.click(): it exercises only the handler, not that a user can reach the icon. The
  icon keeps a real ~8x13 px box while visibility:hidden (verified live), so its centre
  is a valid target; if it ever collapses to zero size the helper throws rather than
  silently falling back to synthetic input.
- **Readiness barriers must not rely on label names/counts alone when a node
  rebuild can produce a same-name, same-count generation** — gate on the source of
  truth (fieldsColumnNames) or on two consecutive agreeing reads.
- **Picked field ORDER is the GROK-14962 signal** — exact sequence equality on the
  drawn labels AND the card [column] elements, never set or subset checks.
- **Column-lifecycle actions assert a zero console-error floor and no error balloon**;
  the number-format guard asserts against the grid column's declared format applied to the raw value, never col.getString(row).
