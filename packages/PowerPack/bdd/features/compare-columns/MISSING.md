# Compare Columns — what the feature does not claim

## Resolved

- core#1246 (Compare Columns missing from the table view's Data menu): fixed in core (`table_view.dart` `fillMenu`
  takes functions with `requiresTable`); the feature picks it from the menu.
- public#4475 (ADD COLUMN over columns of two tables gave a column false in every row): PowerPack
  `src/compare-columns.ts` adds the match as values.

## Won't do

- A column compared with itself, or per column type: one rule (`getString` equality of every row).
- The counts and the row numbers of the length mismatch: text of the same message.
