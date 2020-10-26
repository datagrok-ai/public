<!-- TITLE: Add an Info Panel -->
<!-- SUBTITLE: -->

# Info Panels

[Info panels](../../discover/info-panels.md) are a powerful tool for bringing new context-specific data to the sight. You can inform users about an object they see through these panels, which is why they have such a name. New details typically appear along with the rest of the information in the [property panel](../../overview/navigation.md#properties) and [visibility conditions](#visibility-conditions) will be re-evaluated whenever the object changes.

## Development

Info panels are added as part of a [package](../develop.md). There are two ways of developing them for Datagrok: either as panel [scripts](../scripting.md) or as JavaScript panel [functions](../overview/functions/function.md). Panel scripts can be written in any language supported by the platform (the full list of supported languages is available [here](../scripting.md#supported-languages)). In this case, the main difference between the two implementations pertains to where the code is executed. Panel functions defined in the package entry point will run on the client side, whereas panel scripts get executed on the server.

### Visibility Conditions

### Scripts

To create a panel script, you should tag it as `panel` and specify conditions for the panel to be shown in the `condition` header parameter:

```python
#name: Detect Cats
#description: Detects cats on image
#language: python
#input: file file
#output: bool hasCats
#tags: demo, files, panel, ml, opencv
#condition: file.isFile && file.size < 1e6 && file.path.contains("/cats/") && (file.name.endsWith("jpg") || file.name.endsWith("jpeg"))

import cv2

image = cv2.imread(file)
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
detector = cv2.CascadeClassifier(cv2.data.haarcascades + "haarcascade_frontalcatface.xml")
hasCats = len(detector.detectMultiScale(gray, scaleFactor=1.3, minNeighbors=3, minSize=(75, 75))) != 0
```

Regardless of a script's language, conditions are written in [Grok script](../../overview/grok-script.md).

### Functions

A different approach is used to introduce an info panel in a JavaScript file. In addition to the `panel` tag indication and `condition`, the function should be properly annotated to return a widget. A simplified example is shown below:

```javascript
//name: Translation
//tags: panel, widgets
//input: file file
//output: widget result
//condition: detectTextFile(file)
export function translationPanel(file) {
    return new DG.Widget(ui.divText("Lost in Translation"));
}
```

## Domain Examples

See also:

  * [Info Panels](../../discover/info-panels.md)
  * [Datagrok JavaScript API](../js-api.md)
  * [JavaScript API Samples](https://public.datagrok.ai/js/samples/functions/info-panels/info-panels)
  * [JavaScript Development](../develop.md)
  * [Scripting](../scripting.md)
  * [Functions](../overview/functions/function.md)
