---
feature: stickymeta
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: []
realized_as:
  - 03-persistence-copy-delete.test.ts
related_bugs: []
---

# 03 — Persistence across copy / clone / reload, and deleting metadata

Verifies that sticky metadata is bound to the object (not the in-memory dataframe) and survives
copy/clone/save/reopen/reload/relogin, and that clearing values removes the marker.

**Preconditions**

- A molecule-matching schema (`PW_SM_Schema_<suffix>`) exists with `rating`/int and
  `notes`/string properties (autotest creates and later deletes it).
- `System:DemoFiles/SPGI.csv` open, with metadata seeded on a few molecules (e.g. `rating` and
  `notes` on three rows). Seeding may use the UI (Test 2.1 flow) or, in the autotest, the JS API
  for speed.

> Scope note: clone / new-view / refresh persistence is UI-observable. Save-as-project,
> move-to-space, and binary export/import are driven through the **JS API** in the autotest
> (there is no pure-UI gesture for them that fits a DOM assertion); each such step is annotated
> with a `SCOPE NOTE:` in code.

---

## Test 3.1 — Metadata survives copy / clone / new view

**Steps**

1. With the seeded table open, **clone** the table (right-click view tab → Clone, or the grid
   "Clone" action) and open the clone.
2. Open a **new view** on the same dataframe.
3. For each, inspect the seeded molecules.

**Expected**

- The **blue marker** and the metadata values are present on the same molecules in the cloned
  table and the new view — no metadata is lost.

---

## Test 3.2 — Metadata survives save-as-project and reopen

**Steps**

1. Save the seeded table as a **project**.
2. Close everything, then **reopen** the project.
3. Inspect the seeded molecules.

**Expected**

- After reopen, the seeded molecules still show the same metadata (marker + values).

---

## Test 3.3 — Metadata survives page refresh and relogin

**Steps**

1. Seed metadata on a molecule.
2. **Refresh** the browser tab; re-open the table and inspect.
3. **Log out and log back in**; re-open the table and inspect.

**Expected**

- Metadata remains intact in both cases — nothing is lost across refresh or a fresh session.

---

## Test 3.4 — Delete metadata values and verify removal

**Steps**

1. On a cell that has metadata, open **right-click > Sticky meta > Edit for current cell...**.
2. Clear the `rating` and `notes` values. Save.
3. Refresh the page and/or reopen the table; inspect the same cell.

**Expected**

- The values are removed. If the cell has no remaining sticky meta, the **blue marker disappears**.
- After reload, the metadata is still absent (the deletion persisted).

> Platform note (from recon): the JS `setAllValues` API silently ignores a `null` for an `int`
> property — clearing a numeric value via the API requires a sentinel (e.g. `0`) or deleting the
> schema. The **UI** "Edit for current cell" dialog is the path that clears a value cleanly; the
> autotest prefers the UI for the clear step and documents the API limitation where relevant.

---

## Cleanup

- Delete the project created in Test 3.2.
- Delete the `PW_SM_Schema_<suffix>` schema and its entity type.
- Close all views.
