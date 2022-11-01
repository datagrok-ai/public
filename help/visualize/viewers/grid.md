<!-- TITLE: Grid -->
<!-- SUBTITLE: -->

# Grid

## Selection

|                                 |                                        |
|---------------------------------|----------------------------------------|
| Shift+Mouse Drag                | Select rows                            |
| Ctrl+Shift+Mouse Drag           | Deselect rows                          |
| Mouse Drag row headers          | Select rows                            |
| Shift+drag column headers       | Select columns                         |
| Ctrl+click column headers       | Select columns                         |
| Ctrl+Shift+click column headers | Deselect columns                       |
| (Ctrl+) Shift + ↑↓              | (Un)select rows                        |
| (Ctrl+) Shift + ←→              | (Un)select columns                     |
| (Ctrl+) Shift + mouse-drag      | (Un)select rows                        |
| (Ctrl+) Shift + ENTER           | (Un)Select rows with the current value |
| Ctrl + Shift + Home             | Select rows above current              |
| Ctrl + Shift + End              | Select rows below current              |

## Navigation

|                       |                      |
|-----------------------|----------------------|
| Up, Down, Left, Right | Navigate             |
| Page Up, Page Down    | Navigate             |
| Ctrl+Home             | Jump to first row    |
| Ctrl+End              | Jump to last row     |
| Home                  | Jump to first column |
| End                   | Jump to last column  |

## Sorting

|                           |                |
|---------------------------|----------------|
| Double-click column header | Sort          |
| Menu\| Current Column\| Sort      | Sort |

## Editing

|              |                         |
|--------------|-------------------------|
| Double-click | Edit cell               |
| Ctrl+C       | Copy cell               |
| Ctrl+V       | Paste into cell         |
| Ctrl+click   | Invert column selection |
| ≡            | Open column filter      |
| Enter or the "+" icon at the last row | Add a row (requires `Allow Edit` set to true) |
| Edit \| Add Rows... | Add a specific number of rows at a specified position |

## Resizing and reordering

|                           |                      |
|---------------------------|----------------------|
| Drag column header        | Reorder columns      |
| Drag column header border | Resize columns       |
| Drag row header border    | Copy cell            |
| Menu -> Column sizing     | Batch sizing options |

## Formatting

|                                                |                         |
|------------------------------------------------|-------------------------|
| Right-click column header \|  Format           |  Change datetime format |
| Right-click cell \|  Current column \|  Format |  Change datetime format |
| Column Properties (F2) \|  Tags \|  format     |  Change datetime format |

![Date and number formatting](grid-formatting.gif "Date and number formatting")

## Color coding

### Grid color coding

|              |                         |
|--------------|-------------------------|
| Menu \|  Color coding \|  On/Off  | Turn color-coding on/off for all columns  |
| Menu \|  Color coding \|  Color scheme  | Select a palette  |

![Color-coding](grid-color-coding.gif "Color-coding")

### Column color coding

Color coding can be defined on the column level. Color coding types are
suggested according to the column type: `Categorical` applies to categorical
columns (`string` and `bool`); `Conditional` and `Linear` apply to numeric types
(linear color coding additionally includes `datetime`).

|                                           |                                    |
|-------------------------------------------|------------------------------------|
| Menu \|  Color coding \| Off              | Turn off column color coding       |
| Menu \|  Color coding \| Categorical      | Turn on categorical color coding   |
| Menu \|  Color coding \| Linear           | Turn on linear color coding        |
| Menu \|  Color coding \| Conditional      | Turn on conditional color coding   |
| Menu \|  Color coding \| Edit...          | Edit column color coding           |
| Menu \|  Color coding \| Pick Up Coloring | Clone coloring settings            |
| Menu \|  Color coding \| Apply Coloring   | Apply copied coloring settings     |

Color coding configuration can be copied from one column to another via commands
`Pick Up Coloring` and `Apply Coloring`. Application can work on multiple
columns. Both standard and custom color coding are copied. The `Off` setting can
get picked up as well (when applied, it will turn off the coloring on a column
in question). However, except for the `Off` option, copied settings cannot be
applied if you try to transfer them from a numeric column to a categorical one
or vice versa. This means that the column type is always taken into account. The
`Apply Coloring` command is disabled if nothing has been picked up. Settings are
remembered for a viewer instance (currently, they are not preserved through
layout serialization).

## Row summary columns

Summary columns is a way to visualize multiple values numerical across the row. This feature is useful for quick visual
profiling of values. In the following picture, each inline viewer visualizes the values of five numerical columns, which
allows for quick visual comparison between rows.

![Summary columns](../../uploads/viewers/grid-summary-columns.png "Summary columns")

The following summary column types are available:

* Sparkline
* Bar Chart
* Radar
* Pie Bar Chart
* Markup

To add a summary column: **Menu | Add | Summary Columns**

### Forms

An HTML (or Markdown) template that renders row values can be embedded in each row.

![Forms](../../uploads/viewers/grid-form.png "Forms")

To add a default form: **Menu | Add | Forms | Default**
To add a custom form: **Menu | Add | Forms | Custom...**

### Custom cell renderers

![Molecules](../../uploads/viewers/grid-molecules.png "Molecule renderer")

### Current rows

Rows in a grid can not only be selected or filtered, in addition to that, the grid keeps track of a current row and
highlights it in green. This indication is a neat and lightweight way to update information related to the current value
and lets users explore and compare rows with ease.

To make a row current, simply click on it, or navigate up and down the grid using the cursor up and down keys. Info
panels in the property panel get synchronized with the current cell.

It is also integrated into Datagrok's visualizations and cheminformatics functionality, e.g., similarity search, so as
you move from one row to another you immediately see where the row values belong on the chart or which molecules have
the most similar structure to the reference. This also works the other way around: by first clicking on a visual
element, you will see the row it represents in the grid.

![Current rows](../current-rows-2.gif "Current rows")

## Videos

[![Grid](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=2971s)

See also:

* [Viewers](../viewers.md)
* [Table View](../../datagrok/table-view.md)
* [JS API: Grid](https://public.datagrok.ai/js/samples/ui/viewers/types/grid)
