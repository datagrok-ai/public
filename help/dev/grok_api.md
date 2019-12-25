<!-- TITLE: JavaScript API -->
<!-- SUBTITLE: -->

# JavaScript API

[Grok API](api/grok_api.js.html) allows to control all aspects of the Datagrok platform.

* Data manipulation
* Views
* Pre-defined viewers
* Custom viewers
* Registering functions

## Data manipulation

Use [DataFrame](api/DataFrame.html), [Column](api/Column.html), [ColumnList](api/ColumnList.html), 
and [Row](api/Row.html) classes for table manipulation.

```javascript
demog = gr.testData('demog', 5000);
demog.cols.remove('sex');
foo = demog.cols.addNew('foo', 'int');
demog.rows.removeAt(1, 3);
demog.rows.insertAt(2, 2);
demog.rows.addNew(['Spiderman', 'studyX', 'NYC', 32, 'Spider', 'Net', new Date(2020), 180, 80, 666]);
demog.rows.addNew().subj = 'Iron Man';

// alternative ways of setting values
foo.set(1, 777);
demog.set('age', 1, 44);
``` 

Each [DataFrame](api/DataFrame.html) is associated with two [bitsets](api/BitSet.html): selection and filter.

## Views

Control [views](../features/navigation.md) via the following methods:

```javascript
gr.addTableView(table);
```

Dock an arbitrary visual element in a platform:

```javascript
let e = document.createElement('DIV');
e.innerText = 'This element has been created in JavaScript';
gr.dockElement(e, 'JS', 'left', 0.5);
```

## Pre-defined viewers

[Viewers](../viewers/viewers.md) are very important components of the Datagrok platform. Grok API 
exposes functionality for manipulating pre-defined viewers 
(such as [scatter plot](../viewers/scatter-plot.md) or [histogram](../viewers/histogram.md)), as
well as for developing custom viewers.

Add a new viewer and set up its properties:

```javascript
view = gr.addTableView(gr.testData('demog', 5000));
hist = view.addViewer('histogram');
hist.options({'valueColumnName': 'weight'});
```

## Custom viewers

Extend [JsViewer](api/JsViewer.html) class to develop viewers that become first-class citizens in 
the Grok platform. Once a viewer is registered, you can do the following:

* Add viewer using `Add | JsDemoViewer`, or from the toolbar 'viewers' popup
* Persist viewer as part of the project
* Common viewer operations under the "Viewer" popup menu, such as cloning, embedding, etc

The following code defines a new viewer by subclassing JsViewer. The viewer listens to changes
of filter and selection in the attached table, and updates numbers of filtered/selected rows accordingly.

```javascript
class JsDemoViewer extends JsViewer {
    onFrameAttached(dataFrameHandle) {
        this.dataFrame = new DataFrame(dataFrameHandle);
        this.dataFrame.selection.onChanged(() => this.render());
        this.dataFrame.filter.onChanged(() => this.render());

        this.render();
    }

    render() {
        this.root.innerHTML =
            `${this.dataFrame.toString()}<br>
            Selected: ${this.dataFrame.selection.trueCount}<br>
            Filtered: ${this.dataFrame.filter.trueCount}`;
    }
}

gr.registerViewer('JsDemoViewer', 'JavaScript-based viewer', () => new JsDemoViewer());
```

## Registering functions

Pretty much anything in Grok is a [function](../entities/function.md), it is a concept that
connects together [scripts](../features/scripting.md) written in different languages, predictive models, statistical
functions, query transformations, data flows, and many other features.

The following code registers a "jsConcat" function that becomes a first-class
citizen in the platform (i.e., it can be used from console, gets registered
in help, there could be an optional audit trail associated with the invocations, etc)

```javascript
gr.functions.register({
    signature: 'String jsConcat(int foo, int bar)',
    run: (foo, bar) => `${foo}_${bar}`});
```

Internally, JavaScript-based application are also functions that are annotated accordingly.


## Custom file handlers

If custom file format support is required, just add function to [application](app_development.md) with 
"file-handler" tag in Grok. Input can be string or list of bytes, output is list of 
[tables](../entities/table.md). Extensions are specified in "meta.ext" option and separated with comma. 

```js
//input: string content
//output: list tables
//tags: file-handler
//meta.ext: fasta
fastaFileHandler(content) {
    ...
    return tables;
}
```
