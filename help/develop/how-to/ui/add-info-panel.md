---
title: "Add an info panel"
---

[Info panels](../../../datagrok/navigation/panels/info-panels.md) are a powerful tool for bringing
new context-specific data to the sight. You can inform users about an object
they see through these panels, which is why they have such a name. New details
typically appear along with the rest of the information in the [context panel](../../../datagrok/navigation/panels/panels.md#context-panel) and [visibility conditions](#visibility-conditions) will be re-evaluated whenever the object
changes.

## Development

Info panels are added as part of a [package](../../develop.md). There are two ways
of developing them for Datagrok: either as panel
[scripts](../../../compute/scripting/scripting.mdx) or as JavaScript panel
[functions](../../../datagrok/concepts/functions/functions.md). Panel scripts can be written
in any language supported by the platform (the full list of supported languages
is available [here](../../../compute/compute.md#functions-and-cross-language-support)). In this
case, the main difference between the two implementations pertains to where the
code is executed. Panel functions defined in the package entry point will run on
the client side, whereas panel scripts get executed on the server.

### Visibility conditions

#### Semantic types

Sometimes it is desirable to show an info panel only for data of a specific
[semantic type](../../../govern/catalog/semantic-types.md). To make use of detectors
available out of the box, simply specify a relevant semantic type either from a
script or from a panel function written in JavaScript.

```python
# name: string length
# language: python
# tags: panel
# input: string s {semType: text}
# output: int length
# condition: true

length = len(s)
```

In the above example, the input parameter has `{semType: Text}`. So when you
have a table open and go to a cell in a column with the semantic type `Text`,
you will see this panel. You can use other semantic types in a similar way, for
example, set `{semType: Molecule}` to display properties of various chemical
structures. See the full list of semantic types
[here](../../../govern/catalog/semantic-types.md#automatic-semantic-type-detection).
Since the type is common to all values in a column, it is often convenient to
check it in the panel condition like this:

```Grok Script
condition: columnName.semType == "Molecule"
```

Our [JavaScript API](../../packages/js-api.md) provides the means to override the types
automatically detected by the platform. Refer to this [code
snippet](https://public.datagrok.ai/js/samples/data-frame/semantic-type-detection)
as an example.

It is possible to define your own semantic types. To apply a custom detector,
first define a function that determines the corresponding semantic type in your
[detectors.js](../../develop.md#package-structure) file. Then you can write a
function for auto-detection<!--, as done in our demo package [Pedometer](https://github.com/datagrok-ai/labs/tree/master/packages/Pedometer) -->.

### Scripts

To create a panel script, you should tag it as `panel` and specify conditions
for the panel to be shown in the `condition` header parameter:

```python
# name: DetectCats
# description: Detects cats on image
# language: python
# input: file file
# output: bool hascats
# tags: demo, files, panel, ml, opencv
# condition: file.isfile && file.size < 1e6 && file.path.contains("/cats/") && (file.name.endsWith("jpg") || file.name.endswith("jpeg"))

import cv2

image = cv2.imread(file)
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
detector = cv2.CascadeClassifier(cv2.data.haarcascades + "haarcascade_frontalcatface.xml")
hasCats = len(detector.detectMultiScale(gray, scaleFactor=1.3, minNeighbors=3, minSize=(75, 75))) != 0
```

Regardless of a script's language, conditions are written in [Grok script](../../under-the-hood/grok-script.md) syntax.

### Functions

A different approach is used to introduce an info panel in a JavaScript file. In
addition to the `panel` tag indication and `condition`, the function should be
properly annotated to return a widget. A simplified example is shown below:

```javascript
//name: Translation
//tags: panel, widgets
//input: file file
//output: widget result
//condition: isTextFile(file)
export function translationPanel(file) {
    return new DG.Widget(ui.divText("Lost in Translation"));
}
```

### Semantic value

There is a special type of input `semantic_value`. It is most commonly used to
preserve information about the value's context or its representation in cell.
The following code demonstrates how to get the column that contains the value.

```javascript
//name: get_column
//tags: panel, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export function valueWidget(value) {
    return value ? new DG.Widget(ui.divText('column - ' + value.cell.column.name)) : new DG.Widget(ui.divText('value is empty'));
}
```

## Examples

### Visibility

#### User condition

You can control access to specific types of info panes based on user attributes
such as roles, groups, etc. 

<details>
<summary>Code snippet</summary>

To set the user condition in a script, use the `user` variable like this:

```
# condition: user.name == "john doe" || user.name == "jack smith"
# condition: user.hasrole("chemist")
# condition: user.inteam("high-throughput screening")
```

</details>

#### Dataset condition

Certain info panes are only relevant when applied to data retrieved from
specific data sources. You can define the dataset condition to specify the data
source for the table. 

<details>
<summary>Code snippet</summary>

To set the dataset condition in a script, use the `table` variable like this:

```
# condition: table.gettag("database") == "northwind"
```

</details>

#### Context condition

Info panes accept only one input parameter, which can be a column, a table, a table
cell, or any other [object](../../../datagrok/concepts/objects.md). A condition
may check against that object using the parameter name ("x" in a
sample code snippet below).

<details>
<summary>Code snippet</summary>

```
#input: column x
#condition: x.isnumerical && x.name == "f3" && x.stats.missingvaluecount > 0
```

</details>

### Using predefined visualizations

Oftentimes, it is beneficial to show users an interactive plot, pre-customized
based on the structure of the table that is currently open.

See the following info pane (viewer-scatter.grok) in action by opening
(project:demog). It creates a [Scatterplot](../../../visualize/viewers/scatter-plot.md), sets the axes to the predefined
columns, and adds a regression line.

<details>
<summary>Code snippet</summary>

```
#name: Scatter plot
#description: Panel that contains an interactive Scatter plot
#language: grok
#tags: panel
#input: dataframe table
#condition: table.name == "demog" && table.columns.containsAll(["height", "weight", "age", "sex"])
#output: viewer plot

plot = table.ScatterPlot("height", "weight", "age", "sex")
plot.showRegressionLine = true
```
</details>

### Digital signal processing

Depending on the nature of the data, our platform makes assumptions regarding the ways users would want to analyze it.
For instance, whenever a scientist working in the digital signal processing domain opens a dataset that contains digital
signals, it is likely that they would want to easily see certain features from the frequency domain
(common actions would be analyzing results of the Fourier transform, doing an auto-correlation, or checking out
spectrogram or scalogram). Typically, when such a need arises, the scientist would fire up Matlab, load that dataset,
then either write a script or use Matlab's DSP toolbox (along the way, they would need to enter metadata that is often
not included in the dataset, such as sampling rate).

Within the Datagrok platform, a simple info pane can be developed that would understand when the data is a digital
signal, perform all of the above-mentioned operations in the background (utilizing the metadata stored within the
dataset), and push the results to the user. The result of the script below is a "Spectrogram" pane that would get shown
in the **Context Panel** on the right when user clicks on a column with the digital signal.

See the following info pane (spectrogram-panel.grok) in action by opening (project:eeg)

```
#name: Spectrogram info panel
#description: Panel that contains graphics produced by the R script
#language: grok
#tags: panel,dsp
#input: column signal {type:numerical}
#output: graphics pic
#condition: "F3" == signal.name

pic = Spectrogram("eeg", signal, 256.0, 1024, 0.1, true)
```

### Actions

In addition to providing additional information (such as data, graphics, interactive viewers, etc), it is also possible
to add commands that would do something with the current context. For instance, you might want to send an email to a
user, or update a record in the database.

```
#name: Transaction review panel
#description: Actions available for the credit card transaction
#language: grok
#tags: panel
#input: row activity
#output: string actions {action: markup}
#condition: activity.table.name = "credit card transactions"

actions = "#{button("Flag as suspicious", "http.Post(myserver, row.transactionId)")}"
```


### Predicting molecule solubility

The following pane calculates different molecular properties of a given chemical structure. It appears whenever user
clicks on a structure.

`#{x.ChemScripts:SolubilityPrediction}`

```
#name: Solubility prediction
#description: Predicts solubility by molecule descriptors ("Ipc", "MolWt", "NumValenceElectrons", "MolLogP", "LabuteASA", "TPSA", "HeavyAtomCount", "NumhAcceptors", "NumHDonors", "NumRotatableBonds", "RingCount")
#language: grok
#tags: panel, prediction, chem
#condition: smiles.semtype == "Molecule"
#input: dataframe table
#input: column smiles {semtype: Molecule} [Column with molecules, in smiles format]
#output: dataframe predictions {action: join(table)}
featureNames = ["TPSA", "Ipc", "NumHAcceptors", "NumHDonors", "LabuteASA", "RingCount", "MolWt", "NumValenceElectrons", "HeavyAtomCount", "MolLogP", "NumRotatableBonds"]
ChemDescriptors(table, smiles, featureNames)
MissingValuesImputation(table, featureNames, featureNames, 5)
ApplyModel(Demo:PredictSolubility, table, showProgress=false)
predictions = ExtractColumns(table, ["outcome"])
```


## See also:

* [Info panels](../../../datagrok/navigation/panels/info-panels.md)
* [Datagrok JavaScript API](../../packages/js-api.md)
* [JavaScript API
  Samples](https://public.datagrok.ai/js/samples/functions/info-panels/info-panels)
* [JavaScript development](../../develop.md)
* [Scripting](../../../compute/scripting/scripting.mdx)
* [Functions](../../../datagrok/concepts/functions/functions.md)
* [Semantic types](../../../govern/catalog/semantic-types.md)
