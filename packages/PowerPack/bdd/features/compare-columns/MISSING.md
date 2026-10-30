# Compare Columns feature — what is still missing

## Compare Columns, Filter Builder

Feature: `compare-columns.feature` (30 Oct 2026, local stand).

### Questions (Compare Columns)

- **The second table and second column selectors have no label.** `compare-columns.ts` builds them with
  `ui.input.choice('', …)`, so they carry no caption, `name` or aria-label; the feature reaches them by position
  (`second choice input`, `fourth choice input`), and a screen reader hears two unnamed combo boxes. Should they
  take an accessible name (e.g. "Second table", "Second column")?
- **No help page** covers Compare Columns (only `_old_release-history.md` mentions it). The feature claims what
  the code does.
- **The function's name on the stand.** The source registers it as `compareColumns` (decorator `name`), but the
  local stand lists it as `PowerPack:_compareColumns`, and `PowerPack:compareColumns` answers "Unable to get project
  asset". The feature calls `_compareColumns`; once the menu entry is back (Findings) it should use the menu.

### Findings 30 Oct 2026 (Compare Columns)

- **Data > Compare Columns... is missing from the table view's Data menu.** Steps: open any table; open the Data
  menu. Measured: the menu lists Aggregate Rows... through Categorize..., Calculate MSR..., Curves; no Compare
  Columns.... Expected: the item the function declares (`'top-menu': 'Data | Compare Columns...'`, PowerPack
  `src/package.ts`). Cause: the table view's menu adds a package function only when it takes a dataframe input
  (core `xamgle/lib/src/views/table_view.dart` `fillMenu`, the `f.inputParams.any((p) => p.propertyType ==
  Types.DATA_FRAME)` filter), while `compareColumns` takes none and declares `requiresTable: true` instead (#4226,
  which `Check.forFunc` honours for enabling but the menu filter does not). Curves' Data to Curves
  (`requiresTableView`, no inputs) is missing from Data > Curves the same way. The feature's first scenario is a
  `@known-failure` on this.
- **ADD COLUMN over columns of two tables gives a column that is false in every row.** Steps: open a table cmp
  (a = 1, 2, 3, 4) and a table twin (c = 1, 2, 9, 9); open Compare Columns (Run function `PowerPack:_compareColumns`);
  Tables cmp and twin, Columns a and c; the dialog says Matched: 2; click ADD COLUMN. Measured: the info balloon
  "Column was added to cmp" and an "a == c" column in cmp that is false in all four rows. Expected: true in rows 1
  and 2, as the dialog's own count says. Cause: `compare-columns.ts` adds the calculated column
  `'${a} == ${c}'` to the first table, where no column c exists (a probe `addNewCalculated('a == c', '${a} == ${c}')`
  on cmp gives four falses and no error). The feature's last scenario is a `@known-failure` on this.
