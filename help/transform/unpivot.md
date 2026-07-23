---
title: "Unpivoting"
description: Transform data from a short/wide format to a tall/skinny format using the Unpivot function.
keywords:
  - unpivot
  - melt
  - wide to long
  - reshape data
  - stack columns to rows
---

Unpivoting transforms data from a short/wide to a tall/skinny format.

![unpivot](unpivot.gif)

## Scripting

```
Unpivot(table, copyColumns, mergeColumns, categoryColumnName, valueColumnName, result)

Inputs:
  dataframe table
  string_list copyColumns
  string_list mergeColumns
  string categoryColumnName +
  string valueColumnName +
Output:
  dataframe result
```

Example

```
Unpivot("scores", ["student"], ["math", "english", "history", "science"])
```

See also:

* [Aggregation functions](functions/aggregation-functions.md)
