# Datagrok API Reference: Inputs and Viewers

For interactive scientific applications.

## 1. Inputs (`ui.input.*`)

All inputs return an object inheriting from [`InputBase`](https://datagrok.ai/api/js/dg/classes/InputBase).

UI Documentation: [Datagrok UI](../../help/develop/advanced/ui.md) | Namespace [`ui.input`](https://datagrok.ai/api/js/ui/namespaces/input/)

### 1.1. Input Catalog

| Method | Value Type | Description | Docs |
|---|---|---|---|
| `ui.input.int(label, options?)` | `number` | Integer | [int()](https://datagrok.ai/api/js/ui/namespaces/input/functions/int) |
| `ui.input.float(label, options?)` | `number` | Floating-point number | [float()](https://datagrok.ai/api/js/ui/namespaces/input/functions/float) |
| `ui.input.string(label, options?)` | `string` | String | [string()](https://datagrok.ai/api/js/ui/namespaces/input/functions/string) |
| `ui.input.bool(label, options?)` | `boolean` | Checkbox (true/false) | [bool()](https://datagrok.ai/api/js/ui/namespaces/input/functions/bool) |
| `ui.input.toggle(label, options?)` | `boolean` | Toggle switch | [toggle()](https://datagrok.ai/api/js/ui/namespaces/input/functions/toggle) |
| `ui.input.choice(label, options?)` | `string` | Selection from a list (dropdown) | [choice()](https://datagrok.ai/api/js/ui/namespaces/input/functions/choice) |
| `ui.input.multiChoice(label, options?)` | `string[]` | Multiple selection | [multiChoice()](https://datagrok.ai/api/js/ui/namespaces/input/functions/multiChoice) |
| `ui.input.dateTime(label, options?)` | `DateTime` | Date and time | [dateTime()](https://datagrok.ai/api/js/ui/namespaces/input/functions/dateTime) |
| `ui.input.textArea(label, options?)` | `string` | Multiline text | [textArea()](https://datagrok.ai/api/js/ui/namespaces/input/functions/textArea) |
| `ui.input.search(label, options?)` | `string` | String with search icon, Esc to clear | [search()](https://datagrok.ai/api/js/ui/namespaces/input/functions/search) |
| `ui.input.column(label, options?)` | `DG.Column` | Column selection | [column()](https://datagrok.ai/api/js/ui/namespaces/input/functions/column) |
| `ui.input.columnList(label, options?)` | `DG.Column[]` | Multiple column selection | [columnList()](https://datagrok.ai/api/js/ui/namespaces/input/functions/columnList) |
| `ui.input.markdown(label, options?)` | `string` | Markdown editor | [markdown()](https://datagrok.ai/api/js/ui/namespaces/input/functions/markdown) |

### 1.2. Common Options (`options`)

| Option | Type | Description |
|---|---|---|
| `value` | matches the type | Initial value |
| `items` | `string[]` | List of choices (for `choice`, `multiChoice`) |
| `nullable` | `boolean` | Whether the value can be `null` |
| `min` | `number` | Minimum value (for numeric types) |
| `max` | `number` | Maximum value (for numeric types) |
| `step` | `number` | Step (for numeric types) |
| `placeholder` | `string` | Hint text in an empty field |
| `icon` | `string \| HTMLElement` | Icon in the field (for string types) |
| `clearIcon` | `boolean` | Clear icon (for string types) |
| `escClears` | `boolean` | Esc clears the field (for string types) |

### 1.3. Properties and Methods of [`InputBase`](https://datagrok.ai/api/js/dg/classes/InputBase)

#### Value and State

| Property/Method | Description |
|---|---|
| `.value` | Current value (get/set) |
| `.enabled` | Input enabled state (get/set) |
| `.root` | Root `HTMLElement` (div with label and input) |
| `.input` | `HTMLElement` of the input field itself |
| `.captionLabel` | `HTMLElement` of the label |

#### Tooltips

```typescript
input.setTooltip('Parameter description');
```

#### Validators

```typescript
// Adding a validator: returns null (valid) or an error string
input.addValidator((value) => {
  if (value < 0) return 'Value must be non-negative';
  return null;
});

// Checking validity
const isValid = input.validate(); // true / false
```

#### Events

| Event | Description |
|---|---|
| `.onChanged` | Value changed (by user or programmatically). Subscribe: `.onChanged.subscribe(callback)` |
| `.onInput` | Value changed by user. Subscribe: `.onInput.subscribe(callback)` |
| `.fireChanged()` | Programmatically trigger the `onChanged` event |
| `.fireInput()` | Programmatically trigger the `onInput` event |

### 1.4. Binding Inputs ([`ui.bindInputs`](https://datagrok.ai/api/js/ui/functions/bindInputs))

```typescript
// Combining subscriptions from multiple inputs
const subs: rxjs.Subscription[] = ui.bindInputs([input1, input2, input3]);
```

### 1.5. Grouping Inputs into a Form ([`ui.inputs`](https://datagrok.ai/api/js/ui/functions/inputs) | [UI: Forms](../../help/develop/advanced/ui.md#forms))

```typescript
// Vertical form
const form = ui.inputs([
  ui.input.string('Name'),
  ui.input.int('Age'),
  ui.buttonsInput([
    ui.bigButton('Apply'),
    ui.button('Cancel'),
  ]),
]);

// Also: ui.form([...]), ui.narrowForm([...]), ui.wideForm([...])
```

## 2. Buttons and Icons

| Method | Description | Docs |
|---|---|---|
| `ui.button(text, onClick, tooltip?)` | Standard button | [button()](https://datagrok.ai/api/js/ui/functions/button) |
| `ui.bigButton(text, onClick, tooltip?)` | Accent button (for the primary action) | [bigButton()](https://datagrok.ai/api/js/ui/functions/bigButton) |
| `ui.iconFA(name, onClick, tooltip?)` | FontAwesome icon as a button | [iconFA()](https://datagrok.ai/api/js/ui/functions/iconFA) |
| `ui.iconFAB(name, onClick, tooltip?)` | FontAwesome icon (blue) | [iconFAB()](https://datagrok.ai/api/js/ui/functions/iconFAB) |
| `ui.iconSvg(svgContent, onClick, tooltip?)` | SVG icon as a button | [iconSvg()](https://datagrok.ai/api/js/ui/functions/iconSvg) |

## 3. Tooltips ([`ui.tooltip`](https://datagrok.ai/api/js/ui/classes/Tooltip))

```typescript
// Binding a tooltip to any HTMLElement
ui.tooltip.bind(element, 'Tooltip text');

// Showing a tooltip programmatically
ui.tooltip.show('Text', x, y);

// Hiding a tooltip
ui.tooltip.hide();

// Showing a tooltip for a group of table rows
ui.tooltip.showRowGroup(dataFrame, predicate, x, y);
```

## 4. Viewers ([`DG.Viewer`](https://datagrok.ai/api/js/dg/classes/JsViewer))

Documentation: [Viewers](../../help/visualize/viewers/) | [Viewer API](https://datagrok.ai/api/js/dg/classes/JsViewer) | [UI: Viewers](../../help/develop/advanced/ui.md#viewers)

### 4.1. Standard Viewer Catalog

| Factory Method | TableView Method | Description | Docs |
|---|---|---|---|
| `DG.Viewer.barChart(df, options?)` | `view.barChart(options?)` | Bar chart | [Bar Chart](../../help/visualize/viewers/bar-chart.md) |
| `DG.Viewer.boxPlot(df, options?)` | `view.boxPlot(options?)` | Box plot | [Box Plot](../../help/visualize/viewers/box-plot.md) |
| `DG.Viewer.calendar(df, options?)` | `view.calendar(options?)` | Calendar | [Calendar](../../help/visualize/viewers/calendar.md) |
| `DG.Viewer.correlationPlot(df, options?)` | `view.corrPlot(options?)` | Correlation matrix | [Correlation Plot](../../help/visualize/viewers/correlation-plot.md) |
| `DG.Viewer.densityPlot(df, options?)` | `view.densityPlot(options?)` | Point density | [Density Plot](../../help/visualize/viewers/density-plot.md) |
| `DG.Viewer.filters(df, options?)` | `view.filters(options?)` | Filter set | [Filters](../../help/visualize/viewers/filters.md) |
| `DG.Viewer.form(df, options?)` | `view.form(options?)` | Form (single row) | [Form](../../help/visualize/viewers/form.md) |
| `DG.Viewer.grid(df, options?)` | `view.grid` | Table grid | [Grid](../../help/visualize/viewers/grid.md) |
| `DG.Viewer.heatMap(df, options?)` | `view.heatMap(options?)` | Heat map | [Heat Map](../../help/visualize/viewers/heat-map.md) |
| `DG.Viewer.histogram(df, options?)` | `view.histogram(options?)` | Histogram | [Histogram](../../help/visualize/viewers/histogram.md) |
| `DG.Viewer.lineChart(df, options?)` | `view.lineChart(options?)` | Line chart | [Line Chart](../../help/visualize/viewers/line-chart.md) |
| `DG.Viewer.markup(df, options?)` | `view.markup(options?)` | HTML/Markdown | [Markup](../../help/visualize/viewers/markup.md) |
| `DG.Viewer.matrixPlot(df, options?)` | `view.matrixPlot(options?)` | Matrix of plots | [Matrix Plot](../../help/visualize/viewers/matrix-plot.md) |
| `DG.Viewer.network(df, options?)` | `view.networkDiagram(options?)` | Network diagram | [Network Diagram](../../help/visualize/viewers/network-diagram.md) |
| `DG.Viewer.pcPlot(df, options?)` | `view.pcPlot(options?)` | Parallel coordinates | [PC Plot](../../help/visualize/viewers/pc-plot.md) |
| `DG.Viewer.pieChart(df, options?)` | — | Pie chart | [Pie Chart](../../help/visualize/viewers/pie-chart.md) |
| `DG.Viewer.scatterPlot(df, options?)` | `view.scatterPlot(options?)` | Scatter plot | [Scatter Plot](../../help/visualize/viewers/scatter-plot.md) |
| `DG.Viewer.scatterPlot3d(df, options?)` | `view.scatterPlot3d(options?)` | 3D scatter plot | [3D Scatter Plot](../../help/visualize/viewers/3d-scatter-plot.md) |
| `DG.Viewer.statistics(df, options?)` | `view.statistics(options?)` | Descriptive statistics | [Statistics](../../help/visualize/viewers/statistics.md) |
| `DG.Viewer.tile(df, options?)` | `view.tileViewer(options?)` | Tile view | [Tile Viewer](../../help/visualize/viewers/tile-viewer.md) |
| `DG.Viewer.treeMap(df, options?)` | `view.treeMap(options?)` | Tree map | [Tree Map](../../help/visualize/viewers/tree-map.md) |
| `DG.Viewer.trellisPlot(df, options?)` | — | Facet grid | [Trellis Plot](../../help/visualize/viewers/trellis-plot.md) |
| `DG.Viewer.wordCloud(df, options?)` | — | Word cloud | [Word Cloud](../../help/visualize/viewers/word-cloud.md) |

Additional viewers (require data with coordinates):

| Viewer | Description | Docs |
|---|---|---|
| `view.googleMap(options?)` | Google Maps with data overlay | [Google Map](../../help/visualize/viewers/google-map.md) |
| `DG.Viewer.fromType('Globe', df)` | 3D globe | [Globe](../../help/visualize/viewers/globe.md) |
| `view.shapeMap(options?)` | Region map | [Shape Map](../../help/visualize/viewers/shape-map.md) |

### 4.2. Creating by Type

```typescript
// Creating a viewer by string type
const viewer = DG.Viewer.fromType('Scatter plot', dataFrame);
```

### 4.3. Configuring Options

```typescript
// At creation time
const plot = view.scatterPlot({
  x: 'height',
  y: 'weight',
  size: 'age',
  color: 'race',
});

// After creation
plot.setOptions({
  showRegressionLine: true,
  markerType: 'square',
});
```

### 4.4. Docking to TableView ([UI: Docking](../../help/develop/advanced/ui.md#docking))

```typescript
const view = grok.shell.addTableView(df);

// Docking a viewer
const chart = DG.Viewer.lineChart(df);
view.dockManager.dock(chart, 'right', null, 'Line Chart');

// Docking an arbitrary element
const div = ui.div([/* content */]);
const node = view.dockManager.dock(div, 'down', null, 'Panel', 0.3);

// Docking types: 'left', 'right', 'top', 'down', 'fill'
// Last parameter is the dock ratio (0..1)
```

## 5. Notifications ([`grok.shell`](https://datagrok.ai/api/js/dg/classes/Shell))

```typescript
grok.shell.info('Informational message');
grok.shell.warning('Warning');
grok.shell.error('Error message');
```

## 6. Dialogs ([`ui.dialog`](https://datagrok.ai/api/js/ui/functions/dialog) | [UI: Dialogs](../../help/develop/advanced/ui.md#dialogs))

```typescript
// Standard dialog
ui.dialog('Title')
  .add(ui.inputs([
    ui.input.float('Parameter 1', {value: 1.0}),
    ui.input.float('Parameter 2', {value: 2.0}),
  ]))
  .onOK(() => { /* handling */ })
  .show();

// Modal dialog
ui.dialog('Title')
  .add(/* content */)
  .onOK(() => { /* handling */ })
  .showModal();
```

## 7. Subscriptions and Cleanup

```typescript
// Subscribing to an event
const sub = input.onChanged.subscribe((value) => {
  // handling
});

// Unsubscribing
sub.unsubscribe();

// For viewers
viewer.sub(eventId, callback);        // registers a subscription
viewer.registerCleanup(cleanupFunc);   // will be called on close
```

## 8. Progress Bar

```typescript
const pi = DG.TaskBarProgressIndicator.create('Task description...');
pi.update(50, 'Progress 50%');
// ...
pi.close();
```

## 9. Layouts and Containers ([UI: Layouts](../../help/develop/advanced/ui.md#layouts))

| Method | Description | Docs |
|---|---|---|
| `ui.div([...])` | Container | [div()](https://datagrok.ai/api/js/ui/functions/div) |
| `ui.divH([...])` | Horizontal flex container | [divH()](https://datagrok.ai/api/js/ui/functions/divH) |
| `ui.divV([...])` | Vertical flex container | [divV()](https://datagrok.ai/api/js/ui/functions/divV) |
| `ui.panel([...])` | Panel with padding | [panel()](https://datagrok.ai/api/js/ui/functions/panel) |
| `ui.box(element)` | Fixed-size container | [box()](https://datagrok.ai/api/js/ui/functions/box) |
| `ui.splitH([...])` | Horizontal splitter (resizable) | [splitH()](https://datagrok.ai/api/js/ui/functions/splitH) |
| `ui.splitV([...])` | Vertical splitter (resizable) | [splitV()](https://datagrok.ai/api/js/ui/functions/splitV) |
| `ui.tabControl({...})` | Tabs | [tabControl()](https://datagrok.ai/api/js/ui/functions/tabControl) |
| `ui.accordion()` | Accordion | [accordion()](https://datagrok.ai/api/js/ui/functions/accordion) |
