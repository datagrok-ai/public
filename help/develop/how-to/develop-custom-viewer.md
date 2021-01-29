<!-- TITLE: Develop Custom Viewers -->
<!-- SUBTITLE: -->

# Custom Viewers

Developers can extend Datagrok with special visual components bound to data, which are called
[viewers](../../visualize/viewers.md).There are two major alternatives for developing viewers on
Datagrok. The first one is JavaScript-based development, which lets you create interactive viewers
via [Datagrok JavaScript API](../js-api.md). The second option would be utilizing visualizations
available for popular programming languages, such as Python, R, or Julia. This implementation uses
[scripting](../scripting.md) internally, so the code runs on the server, which slightly affects
interactivity. Nonetheless, both options are applicable and support the essential functionality,
such as data filtering and selection.

Typically, development starts with [package](../develop.md#packages) creation. Packages are
convenient units for distributing content within the platform. Using them you can extend Datagrok
with your widgets, applications, plugins, and more. Besides, scripting viewers must be part of a
package in order to run their code in a separate environment.

Table of contents
  * [JavaScript-Based Viewers](#javascript-based-viewers)
    * [External Dependencies](#external-dependencies)
    * [Properties](#properties)
    * [Events](#events)
  * [Scripting Viewers](#scripting-viewers)
  * [Registering Viewers](#registering-viewers)
  * [Examples](#examples)

## JavaScript-Based Viewers

Let's create an `Awesome` viewer for your package. This can be done in two simple steps. First, we
need to define a subclass of `JsViewer` in a separate file:

```javascript
export class AwesomeViewer extends DG.JsViewer {
  /* AwesomeViewer contents */
}
```

The naming convention assumes that you would add a `Viewer` postfix for such a class. Once this
class is added, you can refer to it in the main JavaScript file of your package: 

```javascript
import {AwesomeViewer} from './awesome-viewer.js'

//name: Awesome
//description: Creates an awesome viewer
//tags: viewer
//output: viewer result
export function awesome() {
  return new AwesomeViewer();
}
```

The annotated function above [registers](#registering-viewers) our viewer and makes it available on
the platform.

Now we will start adding new functionality to this template. If you have
[datagrok-tools](https://www.npmjs.com/package/datagrok-tools) installed and want to follow along,
you can obtain the starter code with these commands:

```
grok create AwesomePackage
grok add viewer AwesomeViewer
```

In case you are new to package development workflow, [this article](../develop.md) will get you
started. Also, visit our JavaScript API Samples, there you can run this [code
snippet](https://public.datagrok.ai/js/samples/functions/custom-viewers/viewers) to get a simple
working example right away and see how it functions within the platform. Finally, if you would like
to explore the methods discussed in the article on your own, jump right to our [JavaScript API
documentation](https://datagrok.ai/js-api/JsViewer).

### External Dependencies

Just like regular `npm` packages, Datagrok packages may have dependencies. Typically you would list
them in your [package.json](../develop.md#package.json) file and get them locally via `npm install`.
For instance, we will use a JavaScript library called [D3](https://d3js.org/) to create a bar chart.
In this case, the `dependencies` field will look as follows:

```json
"dependencies": {
  "datagrok-api": "latest",
  "d3": "^6.2.0"
}
```

Some libraries, such as [datagrok-api](https://www.npmjs.com/package/datagrok-api), are marked as
external modules in a package's [webpack configuration](../develop.md#webpack.config.js). This means
that they will not be included in a bundle file of your package. The platform provides them in its
environment, so if you use such a library, it will be taken from there.

### Properties

First and foremost, let's take a look at how we can define properties of a viewer. They will include
all the parameters you want users to be able to tweak from the UI in the [property
panel](../../overview/navigation.md#properties), be it data to display, a color scheme, or some
numeric values. We need to set a couple of properties for our bar chart:

```javascript
import {axisBottom, axisLeft, scaleBand, scaleLinear, select} from 'd3';

export class AwesomeViewer extends DG.JsViewer {
  constructor() {
    super();

    // Register properties and define fields initialized to properties' default values
    // Properties that represent columns should end with the 'ColumnName' postfix
    this.splitColumnName = this.string('splitColumnName', 'site');
    this.valueColumnName = this.int('valueColumnName', 'age');

    // Provide a list of acceptable values in the field `choices`
    this.valueAggrType = this.string('valueAggrType', 'avg', { choices: ['avg', 'count', 'sum'] });
    this.color = this.string('color', 'steelblue', { choices: ['darkcyan', 'seagreen', 'steelblue'] });
    this.initialized = false;
  }
}
```

In the class `constructor`, you can create properties of the following data types:

  * integer number `this.int(propertyName[, defaultValue, options])`
  * floating point number `this.float(propertyName[, defaultValue, options])`
  * string `this.string(propertyName[, defaultValue, options])`
  * string array `this.stringList(propertyName[, defaultValue, options])`
  * boolean value `this.bool(propertyName[, defaultValue, options])`
  * datetime `this.dateTime(propertyName[, defaultValue, options])`

Follow the naming conventions for JavaScript variables (`lowerCamelCase` will do). Certainly, the
names will be nicely displayed in the property panel: capitalized, with spaces between words (e.g.,
`property` becomes `Property` and `propertyName` becomes `Property Name`). But for them to appear
so, pay attention to the names you use in your code. Consistency in naming also ensures that
property values are updated correctly. When you set a value of a property, either as default or
within the class body as a result of assignment, the specified value appears next to the
corresponding property in UI. Then changing these values is up to users, and only correctly named
properties are synchronized.

Properties are divided into several main groups depending on their nature. Each group has a separate
tab at the property panel. Here is a complete list of these groups and formal criteria determining
whether a property belongs to one of them:

  * `Data`: data-related properties have the `ColumnName` postfix
  * `Colors`: properties changing a color of visualization end with `color`
  * `Axes`: settings controlling the axes contain the word `axis`
  * `Legend`: properties modifying a legend start with the `legend` prefix
  * `Margins`: properties belonging to this group have `margin` in their name
  * `Markers`: properties containing `marker` appear in this tab
  * `Description`: properties should be named `title` or `description`
  * `Misc`: the rest of properties are placed under this tab

Some properties, such as `stringList` and `dateTime`, can only be used to represent columns and,
therefore, require the `ColumnName` postfix. Let's proceed to our example. After defining the
properties, we specify additional chart settings, which will not appear in the property panel, and
initialize the viewer. The next step is to retrieve data for our visualization. This is typically
done in the `onTableAttached` method where you can test if a viewer is applicable to the currently
open table, select the most relevant columns, and so on. However, as we already provided the default
column names, we use this method only to subscribe to events modifying the displayed data (selection
and filtering) and events that have impact on the visualization's size.

```javascript
export class AwesomeViewer extends DG.JsViewer {
  constructor() {...}
  
  // Additional chart settings
  init() {
    this.margin = { top: 10, right: 30, bottom: 30, left: 90 };
    this.xScale = scaleLinear();
    this.yScale = scaleBand();
    this.data = [];
    this.initialized = true;
  }

  // Stream subscriptions
  onTableAttached() {
    this.init();

    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));

    this.render();
  }

  // Cancel subscriptions when the viewer is detached
  detach() {
    this.subs.forEach(sub => sub.unsubscribe());
  }
}
```

The properties of a viewer are updated automatically, however, the exact behavior can be adjusted or
completely overridden in the `onPropertyChanged` method. For example, the following code will show a
message for users that choose a non-categorical, or rather non-string, column:

```javascript
class AwesomeViewer extends DG.JsViewer {
  constructor {...}
  init {...}
  onTableAttached {...}
  detach {...}

  // Override to handle property changes
  onPropertyChanged(property) {
    super.onPropertyChanged(property);
    if (this.initialized) {
      if (property.name === 'splitColumnName' &&
          this.dataFrame.getCol(this.splitColumnName).type !== property.propertyType) {
        grok.shell.info('Wrong property type');
        return;
      }
      this.render();
    }
  }
}
```

### Events

Datagrok makes working with browser events easy. For a start, you can make viewers more interactive
by showing custom tooltips. Their content may be of any kind (see an
[example](https://dev.datagrok.ai/js/samples/ui/tooltips/tooltips)). For instance, there are
[special tooltips](https://dev.datagrok.ai/js/samples/ui/tooltips/row-group-tooltips) for elements
that represent multiple rows in a dataframe. With them, users will immediately see the rows that
match pre-defined criteria and will be able to select them on click (pay attention to the
`handleClick` function in the [example](https://dev.datagrok.ai/js/samples/ui/tooltips/row-group-tooltips)).
Plus, other open viewers will highlight the elements that represent the respective row group as you move the mouse
pointer over the current viewer (you don't need to configure anything, read more about Datagrok's
efficient visualizations [here](../../visualize/viewers.md)). This also works the other way around:
you can show only a selected portion of the data in your viewer. To do that, you need to know which
rows the user selected so that you can narrow the data set for rendering accordingly. And this is
fairly simple: the method `dataFrame.filter.getSelectedIndexes()` gives you exactly the data you need.

## Scripting Viewers

## Registering Viewers

Tagging scripts or functions as `viewers` registers them within the platform. Registering a viewer
makes it available in the top menu and enables common viewer operations, such as cloning, docking,
embedding, and switching to full screen mode. This also means that users can persist this viewer as
part of a [project](../../overview/project.md).

![](leaflet-menu.png "Viewer Menu")

## Examples

You can find more inspiring examples in our [public repository](https://github.com/datagrok-ai/public):

  * JavaScript-based viewers:
    * [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts): constructs graphs of various types using the [Echarts](http://echarts.apache.org) framework
    * [Leaflet](https://github.com/datagrok-ai/public/tree/master/packages/Leaflet): integrates with the [Leaflet](https://leafletjs.com/) library to build interactive maps
    * [Sunburst](https://github.com/datagrok-ai/public/tree/master/packages/Sunburst): uses the [D3](https://d3js.org/) library for a sunburst chart
    * [Viewers](https://github.com/datagrok-ai/public/tree/master/packages/Viewers): showcases creating JavaScript viewers using various visualization libraries
  * Scripting viewers (R, Python, Julia):
    * [ChaRPy](https://github.com/datagrok-ai/public/tree/master/packages/ChaRPy): translates a Datagrok viewer to Python and R code using scripting viewers for the respective programming languages
    * [JuliaScripts](https://github.com/datagrok-ai/public/tree/master/packages/JuliaScripts): demonstrates the scripting functionality for Julia
    * [PythonScripts](https://github.com/datagrok-ai/public/tree/master/packages/PythonScripts): contains Python scripts for a number of applications, ranging from plotting to computer vision
    * [RScripts](https://github.com/datagrok-ai/public/tree/master/packages/RScripts): offers a variety of R visualizations

    Most of these scripts are also available by the `viewers` tag in the script browser: [https://public.datagrok.ai/scripts?q=%23viewers](https://public.datagrok.ai/scripts?q=%23viewers).

## Videos

<iframe width="560" height="315" src="https://www.youtube.com/embed/JaJgxtHAb98?start=202" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

See also:

  * [Datagrok JavaScript API](../js-api.md)
  * [JavaScript API Samples](https://public.datagrok.ai/js/samples/functions/custom-viewers/viewers)
  * [JavaScript Development](../develop.md)
  * [Viewers](../../visualize/viewers.md)
  * [Custom Viewers](../js-api.md#custom-viewers)
  * [Scripting Viewers](../../visualize/viewers/scripting-viewer.md)
