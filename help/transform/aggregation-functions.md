<!-- TITLE: Aggregation Functions -->
<!-- SUBTITLE: -->

# Aggregation Functions

The same set of aggregation functions is used across the whole platform in different context 
(["Aggregate Rows" dialog](aggregate-rows.md), transformations, [markup](../features/markup.md)) 

## Common
| Function | Description|
|----------|------------|
| first    | first (any) value                              |
| count    | total number of all values (including missing) |
| #values  | number of non-missing values                   |
| #unique  | number of unique values                        |
| #nulls   | number of missing values                       |


## Numerical
| Function | Description|
|----------|------------|
| min      |            |
| max      |            |
| sum      |            |
| med      |            |
| avg      |            |
| stdev    |            |
| variance |            |
| skew     | skewness   |
| kurt     | kurtosis   |
| q1       | first quartile, the median of the lower half |
| q2       | second quartile, same as med           |
| q3       | third quartile, the median of the upper half |

## Textual
| Function | Description|
|----------|------------|
| min      |            |
| max      |            |
| sum      |            |
| med      |            |
| avg      |            |
| stdev    |            |
| variance |            |
| skew     | skewness   |
| kurt     | kurtosis   |
| q1       | first quartile, the median of the lower half |
| q2       | second quartile, same as med           |
| q3       | third quartile, the median of the upper half |

## User-defined functions

Datagrok API provides a way to register custom aggregation functions; once registered, they will
automatically appear in UI across the whole platform. 

Markup syntax for referring to aggregations: `t.stats.avg(AGE)`

See also:
  * [Aggregate rows](aggregate-rows.md)
 