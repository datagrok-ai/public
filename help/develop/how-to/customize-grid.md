<!-- TITLE: Grid Customization -->

# Grid Customization

Whenever users open a spreadsheet on the platform, they see it presented in a [grid](../../visualize/viewers/grid.md). This view can be altered in many ways, both from the user interface and programmatically. The aspects subject to change include:

* [formatting](#formatting) numbers and dates
* [ordering](#ordering) rows and columns
* controlling the [visibility](#column-visibility) of columns
* changing cell [size](#resizing)
* encoding categories with [colors](#color-coding) and many more

## Formatting

The Datagrok platform automatically sets the most appropriate format for numbers and dates in grid columns of a particular dataset. However, the way such columns look may be changed. The global formatting settings allow users to specify their default formats:

![](../../uploads/navigation/user-settings-formatting.png "Settings | Format")

In addition to that, users may change the format directly in [column properties](../../visualize/viewers/grid.md#formatting) of an open table. This is done by defining the value of `format` tag. The syntax and standard formats are described in this [article](../../discover/tags.md#format) and given for reference in the respective [code snippet](https://public.datagrok.ai/js/samples/grid/data-format). Here is a brief example to illustrate how the representation of data can be changed:

```javascript
let view = grok.shell.addTableView(grok.data.demo.demog());

view.grid.col('height').format = 'scientific';
view.grid.col('weight').format = '#.0000';
view.grid.col('started').format = 'dd.MM.yyyy';
```

If a format is specified on a column level, it takes precedence over global settings. Another thing to keep in mind is that the tag applies only to numeric and datetime columns. Columns of other data types will ignore the `format` tag. As the actual values remain unchanged, sorting and filtering produce the same result regardless of selected data representation.

## Column Visibility

It is possible to hide parts of data in a grid without actually removing them (see an [example](https://public.datagrok.ai/js/samples/grid/hide-columns)). To achieve this, either explicitly specify the columns you intend to show:

```javascript
let data = grok.data.demo.demog();
let view = grok.shell.addTableView(data);
view.grid.columns.setVisible(['age', 'sex', 'race']);
```

or hide the particular columns by adding `~` prefix to their names before creating a table view:

```javascript
let data = grok.data.demo.demog();
data.columns.byName('age').name = '~age';
let view = grok.shell.addTableView(data);
```

The approach largely depends on which columns are easier to list, but not only. Thus, to make the prefixed column visible again, it is not enough to invoke `setVisible(['~columnName'])`, you have to invert the previous change first. Besides, such renaming affects all table views derived from the data.

See also:
  * [Grid](../../visualize/viewers/grid.md)
  * [Table View](../../overview/table-view.md)
  * [JavaScript API: Grid](https://datagrok.ai/js-api/Grid)
