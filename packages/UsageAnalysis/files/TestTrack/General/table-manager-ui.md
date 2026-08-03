---
feature: general
target_layer: manual-only
coverage_type: regression
produced_from: split
original_path: public/packages/UsageAnalysis/files/TestTrack/General/table-manager.md
split_date: 2026-06-16
related_bugs: []
manual_only_reason: |
  The Table Manager is a d4 canvas grid. Row selection (single and Shift+range
  multi-select), the "Open as table" right-click action, and the "Show > All"
  attribute-column toggle are all driven by pixel coordinates on the canvas with
  no stable DOM/API hook. Opening the manager and verifying it lists the open
  tables IS automated in table-manager-spec.ts; these canvas/context-menu
  interactions remain manual.
---

# Table Manager — canvas grid interactions

Manual companion to `table-manager.md`. Opening **View | Tables** (Alt+T) and
verifying the manager lists all open tables is automated by
`table-manager-spec.ts`. The interactions below operate on the canvas grid and
its context menus and are verified manually.

## Preconditions

1. Load 3-4 datasets (e.g. drag in a couple of demo files).
2. Open the Table Manager via **View | Tables** (or `Alt + T`).

## Steps

1. **Single-row navigation + Context Panel.**
   - Click a dataset row in the Table Manager — a tab with the selected table
     opens / activates.
   - The selected dataset is reflected on the Context Panel.

2. **"Open as table".**
   - Right-click a row → **Open as table**.
   - A new tab opens containing a table built from the Table Manager's dataset
     metadata, and a new row appears in the Table Manager.

3. **Shift multi-select.**
   - With `Shift` held, select all added tables in the Table Manager.
   - The context menu shows a submenu whose actions apply to all selected tables
     (its title reads "N tables", e.g. "3 tables").
   - The Context Panel shows actions related to all selected tables.

4. **"Show" submenu — attribute columns.**
   - Right-click the Table Manager view → expand the **Show** submenu (lists the
     table attributes that can be added as columns, plus an **All** item).
   - Click **All**: all table attributes are added as columns; **All** and each
     attribute show a ✓ (cells are empty where an attribute is not available for
     a given table).
   - Click **All** again: the attribute columns are removed and the ✓ marks
     clear.

---
{
  "order": 6
}
