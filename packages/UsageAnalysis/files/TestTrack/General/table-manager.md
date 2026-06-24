### Table Manager

**Table Manager**, available via **View | Tables**, contains a list of currently open tables. Use it to navigate between tables, select them, or perform batch actions. It also allows to view metadata on multiple tables in
a tabular format. 

The implementation is based on the **grid**, so many of the grid's features apply.

- Precondition: load 3-4 random datasets.

1. Open "Table" dialog via **View | Tables** (or ```Alt + T```):
* A list of opened datasets opens.
* The Table Manager lists every currently open table.

> Automated by `table-manager-spec.ts` (opening the manager via View | Tables and
> verifying it lists all open tables).
>
> The canvas-grid interactions — single-row navigation + Context Panel,
> "Open as table", Shift multi-select ("N tables" submenu), and the "Show > All"
> attribute-column toggle — are manual: see `table-manager-ui.md`.

---
{
  "order": 6
}
