---
feature: stickymeta
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [add-and-edit-metadata-on-cell, materialise-sticky-columns-in-grid, edit-values-in-context-panel-then-indicator-fills, materialise-properties-as-grid-columns]
realizes: []
realized_as:
  - 02-add-and-edit.test.ts
related_bugs: []
---

# 02 — Add and edit metadata on a dataframe (cell, sticky column, batch)

Covers applying sticky metadata to data points in a grid: a single cell, a whole-dataset
"sticky column", and a multi-row batch edit.

**Preconditions**

- A metadata schema matching molecules exists (entity type `semtype=molecule`, with properties
  `rating`/int, `notes`/string, `verified`/bool, `review_date`/datetime). The autotest creates
  its own `PW_SM_Schema_<suffix>` in `beforeAll` and deletes it in `afterAll` — it does not rely
  on a pre-seeded schema.
- `System:DemoFiles/chem/SPGI.csv` opens with a `Structure` column detected as `Molecule`
  (3624 rows).

---

## Test 2.1 — Add metadata to a single cell

**Steps**

1. Open `SPGI.csv` and wait for chem rendering (the `Structure` cells render molecules).
2. Right-click a `Structure` cell (row 0) → **Sticky meta** → **Edit for current cell...**.
   A "Sticky meta" dialog opens with the schema's section.
3. Fill: `rating` = 5, `notes` = "test note", `verified` = ON. (`review_date` may auto-populate.)
4. Click the schema's **Save** button in the dialog.
5. Close the dialog. Re-open it on the same cell (right-click → Sticky meta → Edit for current
   cell...).

**Expected**

- A **blue marker** (dark-blue dot, top-right of the cell) appears on the cell after save,
  indicating it now has metadata. A hover tooltip shows the values.
- On re-open, the saved values are retained: `rating` = 5, `notes` = "test note",
  `verified` = true.

> Note: the **Context Panel > Sticky meta** pane, when the current object is a *cell*, does not
> render editable inputs — the edit path for a single cell is the **right-click → Edit for current
> cell...** dialog. The Context Panel pane is for *column*-level actions (see Test 2.2).

---

## Test 2.2 — Create a sticky column

**Steps**

1. With `SPGI.csv` open, select the `Structure` column (so the Context Panel targets the column).
2. Open the Context Panel **Sticky meta** pane. It shows one section per matching schema.
3. In the schema's section, click **+** "Add all properties as columns" (or the per-property **+**).
4. Wait for the sticky columns to appear (schema matching is asynchronous).
5. Sort ascending by the `rating` sticky column.
6. Remove the `rating` column from the view, then re-add it from the Sticky meta pane.

**Expected**

- Sticky columns are added for the schema properties (`rating`, `notes`, `verified`,
  `review_date`). Each sticky column header carries a **blue circle** marker.
- Sticky columns behave like regular columns — sorting works.
- Removing a sticky column from the view does **not** delete its metadata: after re-adding,
  the previously entered values reappear (e.g. the cell edited in Test 2.1 still shows
  `rating` = 5).

---

## Test 2.3 — Batch edit metadata on multiple rows

**Steps**

1. Select several rows (e.g. shift-click or via the row header — three rows).
2. Right-click a `Structure` cell in the selection → **Sticky meta** → **Edit for all
   properties**. (The **Rows** selector defaults to **Selected** when multiple rows are selected.)
3. Set `verified` = ON and `notes` = "batch note". Apply.

**Expected**

- The metadata is applied to **all** selected rows — each selected row shows
  `verified` = true and `notes` = "batch note" (verify via the sticky columns or per-cell tooltip).
- No selected row is skipped.

---

## Cleanup

- Delete the `PW_SM_Schema_<suffix>` schema and its entity type (removes all values written
  during this scenario).
- Close all views (`grok.shell.closeAll()`).
