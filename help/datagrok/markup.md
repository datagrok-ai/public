<!-- TITLE: Markup -->
<!-- SUBTITLE: -->

# Markup

A powerful mechanism that allows to embed platform-specific visual elements right into the HTML document. Simply embed
the expression like that: `#{expression}`

## Expressions

| Description         | Expression                                                 | Output                     |
|---------------------|------------------------------------------------------------|----------------------------|
| Any entity          | x.User.196900d0-f9e0-11e7-c473-6fef20f86714."Andrew"       | Andrew                     |
| Function            | x.Func.MissingValuesImputation                             | MissingValuesImputation    |
| Aggregations        | t.aggr.sum($AGE)                                           | 267266.00                  |
| Date                | date                                                       | 2018-03-07 22:53:01.423700 |
| Table name          | t.name                                                     | demog                      |
| Table tag           | t.tags\[source.file\]                                      | demog.xlsx                 |
| Row count           | t.rowCount                                                 | 5850                       |
| Selected rows count | t.selection.trueCount                                      | 55                         |
| Filtered rows count | t.filter.trueCount                                         | 5850                       |
| Current filters     | t.rows.filters                                             |                            |
| Current row         | t.currentRow                                               | 12                         |
| Current column      | t.currentCol                                               | USUBJID                    |
| Current cell        | t.currentCell                                              | X0273T21000400001          |
| Current row value   | t.row\[USUBJID\]                                           | X0273T21000200001          |
| Color               | color(AGE)                                                 | #e5e5ff                    |
| Formula             | formula(${AGE} * 2)                                        | 116                        |
| Value editor        | t.editor\[USUBJID\]                                        |                            |
| Statistics          | t.stats.avg(AGE)                                           | 45.69                      |
| Chart               | chart{"type":"Histogram","look":{"valueColumnName":"AGE"}} |                            |

### Any entity

Shows any platform entity.

Format:

```
x.<type>.<id>."<name>"
```

### Function

Shows function.

Format:

```
x.Func.<name>.<option>
```

Options:

* **run** - runs action by click.

### Aggregations

Shows aggregated value of defined column.

Format:

```
t.aggr<aggregation>($<column name>)
```

| Aggregations | Description               |
|--------------|---------------------------|
| count        | Number of rows            |
| #values      | Number of not null values |
| #unique      | Number of unique values   |
| #nulls       | Number of nulls           |
| min          | Minimum                   |
| max          | Maximum                   |
| sum          | Sum                       |
| med          | Median                    |
| avg          | Average                   |
| stdev        | Standard deviation        |
| variance     | Variance                  |
| skew         | Skewness                  |
| kurt         | Kurtosis                  |
| q1           | First quartile            |
| q2           | Second quartile           |
| q3           | Third quartile            |

### Date

Shows current date and time.

Format:

```
date
```

### Table name

Shows table name.

Format:

```
t.name
```

### Table tag

Shows table tag.

Format:

```
t.tags[<tag name>]
```

### Row count

Shows number of rows.

Format:

```
t.rowCount
```

### Selected rows count

Shows number of selected rows.

Format:

```
t.selection.trueCount
```

### Filtered rows count

Shows number of filtered rows.

Format:

```
t.filter.trueCount
```

### Current filters

Shows list of selected filters.

Format:

```
t.rows.filters
```

### Current row

Shows current row index.

Format:

```
t.currentRow
```

### Current column

Shows current column index.

Format:

```
t.currentCol
```

### Current cell

Shows current cell value.

Format:

```
t.currentCell
```

### Current row value

Shows current row value in selected column.

Format:

```
t.row[<column name>]
```

### Color

Shows column's color in HEX format (32 bits).

Format:

```
color(<column name>)
```

### Formula

Shows result of defined formula, see formulas format [there](../transform/add-new-column.md).

Format:

```
formula(<formula>)
```

### Value editor

Shows text input of specified column in selected row.

Format:

```
t.editor[<column name>]
```

### Statistics

Shows aggregated value of defined column.

Format:

```
t.stats.<statistics>(<column name>)
```

| Statistics   | Description               |
|--------------|---------------------------|
| count        | Number of rows            |
| #values      | Number of not null values |
| #unique      | Number of unique values   |
| #nulls       | Number of nulls           |
| min          | Minimum                   |
| max          | Maximum                   |
| sum          | Sum                       |
| med          | Median                    |
| avg          | Average                   |
| stdev        | Standard deviation        |
| variance     | Variance                  |
| skew         | Skewness                  |
| kurt         | Kurtosis                  |
| q1           | First quartile            |
| q2           | Second quartile           |
| q3           | Third quartile            |

### Chart

Chart markup allows to embed any of [Viewers](../visualize/viewers.md) into page.

Format:

```
chart{"type":"<chart type>","look":{<chart parameters>}}
```

Chart types:

* Scatter plot
* Line chart
* Bar chart
* Density plot
* Tree map
* Matrix plot
* Histogram
* Filters
* 3D scatter plot
* Grid
* Column viewer
* Calendar
* Trellis plot
* Box plot
* Pie chart
* Heat map
* PC Plot
* Statistics
* Correlation plot
* Map
* Form designer
* Web viewer
* Markup Viewer
* Word cloud
* Network diagram
* Protein viewer
* Globe
* Card
* Viewer host
* R Viewer

See also:

* [Markup viewer](../visualize/viewers/markup.md)
