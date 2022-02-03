<!-- TITLE: Add an info panel -->
<!-- SUBTITLE: -->

# Info panels

[Info panels](../../discover/info-panels.md) are a powerful tool for bringing new context-specific data to the sight.
You can inform users about an object they see through these panels, which is why they have such a name. New details
typically appear along with the rest of the information in the [property panel](../../overview/navigation.md#properties)
and [visibility conditions](#visibility-conditions) will be re-evaluated whenever the object changes.

## Development

Info panels are added as part of a [package](../develop.md). There are two ways of developing them for Datagrok: either
as panel [scripts](../../compute/scripting.md) or as JavaScript panel [functions](../../overview/functions/function.md).
Panel scripts can be written in any language supported by the platform (the full list of supported languages is
available [here](../../compute/scripting.md#supported-languages)). In this case, the main difference between the two
implementations pertains to where the code is executed. Panel functions defined in the package entry point will run on
the client side, whereas panel scripts get executed on the server.

### Visibility conditions

#### Semantic types

Sometimes it is desirable to show an info panel only for data of a
specific [semantic type](../../discover/semantic-types.md). To make use of detectors available out of the box, simply
specify a relevant semantic type either from a script or from a panel function written in JavaScript.

```python
# name: string length
# language: python
# tags: panel
# input: string s {semtype: text}
# output: int length
# condition: true

length = len(s)
```

In the above example, the input parameter has `{semType: Text}`. So when you have a table open and go to a cell in a
column with the semantic type `Text`, you will see this panel. You can use other semantic types in a similar way, for
example, set `{semType: Molecule}` to display properties of various chemical structures. See the full list of semantic
types [here](../../discover/semantic-types.md#automatic-semantic-type-detection). Since the type is common to all values
in a column, it is often convenient to check it in the panel condition like this:

```
condition: columnName.semType == "Molecule"
```

Our [JavaScript API](../js-api.md) provides the means to override the types automatically detected by the platform.
Refer to this [code snippet](https://public.datagrok.ai/js/samples/data-frame/semantic-type-detection) as an example.

It is possible to define your own semantic types. To apply a custom detector, first define a function that determines
the corresponding semantic type in your [detectors.js](../develop.md#package-structure) file. Then you can write a
function for auto-detection, as done in our demo
package [Pedometer](https://github.com/datagrok-ai/public/tree/master/packages/Pedometer).

### Scripts

To create a panel script, you should tag it as `panel` and specify conditions for the panel to be shown in
the `condition` header parameter:

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

Regardless of a script's language, conditions are written in [Grok script](../../overview/grok-script.md) syntax.

### Functions

A different approach is used to introduce an info panel in a JavaScript file. In addition to the `panel` tag indication
and `condition`, the function should be properly annotated to return a widget. A simplified example is shown below:

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

## Domain examples

See also:

* [Info panels](../../discover/info-panels.md)
* [Datagrok JavaScript API](../js-api.md)
* [JavaScript API Samples](https://public.datagrok.ai/js/samples/functions/info-panels/info-panels)
* [JavaScript development](../develop.md)
* [Scripting](../../compute/scripting.md)
* [Functions](../../overview/functions/function.md)
* [Semantic types](../../discover/semantic-types.md)
