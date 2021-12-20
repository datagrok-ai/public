<!-- TITLE: Show formula lines -->

# How to show formula lines

Dataframe and viewer can contain information about lines and bands.

These figures are used by some viewers to draw additional lines on the charts. These can be reference lines, highlighting different areas of the chart and data, etc.

![Example of formula lines](../../uploads/viewers/formula-lines-example.png)

Lines information is stored in a special storage in a dataframe or viewer. The viewer automatically reads storages when it connects to the dataframe.

You can create and modify information about lines and bands by changing the `.formula-lines` dataframe tag or by changing the `formulaLines` property of the viewer. The content of these storages is a JSON string.

There is a more convenient ways to create lines:

- method `dataframe.meta.addFormulaLine()` - for creating and saving a line in a dataframe
- method `viewer.meta.addFormulaLine()` - for creating and saving a line in a viewer

To completely remove lines from dataframe or viewer use methods `dataframe.meta.removeFormulaLines()` or `viewer.meta.removeFormulaLines()`. Calling these methods without parameters will delete all lines from the storage. You can also remove only specific lines by listing their IDs, for example: `viewer.meta.removeFormulaLines('123', 'abc')`. If a complete removal is not necessary, then you can simply hide the unnecessary line using the `visible` attribute.

Lines saved in the dataframe will be displayed on all viewers in the same way. Lines saved in the viewer are displayed only in this viewer and do not affect other viewers.

An example of creating and displaying a line in this way:

```javascript
let demog = grok.data.demo.demog(100);

// Add line to dataframe:
demog.meta.addFormulaLine({
  title: 'Parabola',
  formula: '${height} = 180 + 0.01 * ${weight} * ${weight} - 1.5 * ${weight}',
  zindex: -30,
  color: "#FFA500",
  width: 2,
  visible: true,
});

// Add another line to dataframe:
demog.meta.addFormulaLine({
  id: 'MyLine',
  formula: '${height} = 200'
});

// Remove line with id = 'MyLine':
demog.meta.removeFormulaLines('MyLine');

let view = grok.shell.addTableView(demog);

let plot = view.scatterPlot({
  x: 'weight',
  y: 'height',
  showDataframeFormulaLines: true,    // Hide or show all lines stored in the dataframe.
  showViewerFormulaLines: true        // Hide or show all lines stored in the viewer.
});

// Add line to viewer:
plot.meta.addFormulaLine({
  formula: '${weight} = 150',
  color: "#ff0000",
  width: 10
});

```

A similar methods is used to create bands - `dataframe.meta.addFormulaBand()` or `viewer.meta.addFormulaBand()`. Most of the parameters for lines and bands are the same. But there are also some parameters that are specific for lines and bands. See them in the description of the parameters for lines and bands.

More examples of creating lines and bands can be found [here](https://public.datagrok.ai/js/samples/data-frame/metadata/formula-lines).

## Line parameters

Method to create a line: `dataframe.meta.addFormulaLine(parameters)` or `viewer.meta.addFormulaLine(parameters)`

Only one parameter ("formula") is required. All other parameters have their default values.

| Parameter       | Type       | Example              | Default              | Description                                                          |
|-----------------|------------|-------------------------------|-------------------------------|----------------------------------------------------------------------|
| `id`           | string     | '123'              | Empty string              | Line ID. Used to completely remove unnecessary line.            |
| `title`           | string     | 'Reference line'              | Empty string              | Short name of the line used when displaying the tooltip.             |
| `description`     | string     | 'Normal distribution of data' | Empty string | Detailed description of the line used when displaying the tooltip.   |
| `formula`        | string     | '${height} = 2.2 * ${weight}' | Required parameter that must be specified | Formula for line. There should be one column to the left of the "=". And any formula using the second column on the right side. The formula uses syntax and formulas similar to the [Add New Column](../../transform/add-new-column.md) form.   |
| `color`           | string     | '#FF0000'                     | '#838383' (dark gray color)         | Line color in HEX format.   |
| `zindex`          | integer    | 25                            | 100                            | The "depth" of the line along the Z axis. The higher the number - the higher the line is located, overlapping other lines with a lower zindex value. The viewer's chart itself has zindex = 0. Values less than zero lead to the placement of lines under the chart. Values greater than zero cause lines to be placed on top of the chart. Lines with the same depth value are displayed in the order in which they were created.   |
| `opacity`           | float     | 0.7                     | 1.0                     | Opacity is a number in the range [0..1], where 0 is completely invisible, 1 is completely opaque.   |
| `visible`           | boolean     | false                     | true                     | Indicates whether a line is displayed or hidden.   |
| `min`           | float     | 50                     | No minimum limit                     | Line boundaries along the value axis. In this example, the line will be drawn for a "Weight" greater than 50 kg.   |
| `max`           | float     | 300                     | No maximum limit                      | Line boundaries along the value axis. In this example, the line will be drawn for a "Weight" less than 300 kg.   |
| `width`           | float     | 3                     | 1                     | Line width in pixels.   |
| `spline`           | float     | 0.5                     | 0.9                     | Smoothness of curve line in range [0..1], where 0 - no smoothing, 1 - max smoothing.   |
| `style`           | string     | 'dashed'                     | 'solid'                     | Line style. Possible styles: 'solid', 'dotted', 'dashed', 'longdash', 'dotdash'.   |

## Band parameters

Method to create a band: `dataframe.meta.addFormulaBand(parameters)` or `viewer.meta.addFormulaBand(parameters)`

Only 3 parameters ("formula", "column" and "column2") are required. All other parameters have their default values.

| Parameter       | Type       | Example              | Default              | Description                                                          |
|-----------------|------------|-------------------------------|-------------------------------|----------------------------------------------------------------------|
| `id`           | string     | '123'              | Empty string              | Line ID. Used to completely remove unnecessary line.            |
| `title`           | string     | 'Clipping range'              | Empty string              | Short name of the band used when displaying the tooltip.             |
| `description`     | string     | 'Ignored range of data'                | Empty string | Detailed description of the band used when displaying the tooltip.   |
| `formula`        | string     | '< 40' | Required parameter that must be specified | Band boundary formula. The formula can contain expressions of the form: "< 200", "> 50", "in(18, 60)", "in(q1, q3)", etc. The numbers are specified in the units of the column (in this case in centimeters).    |
| `color`           | string     | '#00FF00'                     | '#F0F0F0' (light gray color)         | Band color in HEX format.   |
| `zindex`          | integer    | -10                            | 100                            | The "depth" of the band along the Z axis. The higher the number - the higher the band is located, overlapping other lines and bands with a lower zindex value. The viewer's chart itself has zindex = 0. Values less than zero lead to the placement of lines and bands under the chart. Values greater than zero cause bands to be placed on top of the chart. Lines with the same depth value are displayed in the order in which they were created.   |
| `opacity`           | float     | 0.7                     | 1.0                     | Opacity is a number in the range [0..1], where 0 is completely invisible, 1 is completely opaque.   |
| `visible`           | boolean     | false                     | true                     | Indicates whether a band is displayed or hidden.   |
| `min`           | float     | 50                     | No minimum limit                     | Band boundaries along the value axis. In this example, the band will be drawn for a "Weight" greater than 50 kg.   |
| `max`           | float     | 300                     | No maximum limit                      | Band boundaries along the value axis. In this example, the band will be drawn for a "Weight" less than 300 kg.   |
| `column`           | string     | ${weight}                     | Required parameter that must be specified                     | Column for which the band is set. In this example, the formula means that ${weight} < 40. |
| `column2`           | string     | ${height}                     | Required parameter that must be specified                     | Second column for which the band will be drawn.   |

See also:

- [Examples of using formula lines](https://public.datagrok.ai/js/samples/data-frame/metadata/formula-lines)
- [Adding new columns](help/transform/add-new-column)
- [Math functions](https://datagrok.ai/help/transform/functions/math-functions)
- [Operators](https://datagrok.ai/help/transform/functions/operators)
- [Constants](help/transform/functions/constants)
- [Statistical functions](https://datagrok.ai/help/transform/functions/stats-functions)
- [Conversion functions](https://datagrok.ai/help/transform/functions/conversion-functions)
- [Date and time functions](https://datagrok.ai/help/transform/functions/datetime-functions)
- [Text functions](https://datagrok.ai/help/transform/functions/text-functions)
