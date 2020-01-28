<!-- TITLE: Unpivoting -->
<!-- SUBTITLE: -->

# Unpivoting

Unpivoting transforms data from a short/wide to a tall/skinny format.

![Aggregation](../uploads/gifs/aggregate.gif "Aggregation")

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
* [Aggregation functions](../features/aggregation-functions.md)

