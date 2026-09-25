---
title: "Change type"
description: Convert a column to a different data type, choosing a format and handling non-convertible values.
keywords:
  - change column type
  - convert column data type
  - cast column
  - column format conversion
  - retype column
---

Use this dialog to change the data type of one or several columns. To open it, right-click a column header and
select **Change Type...**. To undo the change, press <kbd>Ctrl+Z</kbd>.

## Parameters

| Parameter       | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| Columns         | Columns to convert                                                          |
| Type            | Current type, read-only. Shows `(mixed)` when the columns have different types |
| New type        | Type to convert to. The list offers only the types that all selected columns can be converted to |
| Format          | Format to parse the values with, such as a date format                      |
| Not convertible | Value to use for the values that can't be converted                         |

Under the parameters, the dialog shows how many values are **Convertible**, **Not convertible**, and **Empty**. Next
to each count, use the icons to go to the previous or next such row, filter the rows, or select them.

See also:

* [Column](../datagrok/concepts/table.md#column)
