FileSchemasDemo is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform.
It demonstrates the capabilities for file browsing and integration with the metadata.

Under `/files`, we have different pictures in two folders, "cars" and "animals".

Under `/schemas`, we have defined the following three schemas:

  * [pictures.yaml](schemas/auto.yaml) — common picture properties (width, height)
  * [auto.yaml](schemas/auto.yaml) — automobile properties (make, model, hasPlates)
  * [mammals.yaml](schemas/auto.yaml) — mammal properties (color, isCat)

In [files/meta.yaml](files/meta.yaml), we assign schemas to the corresponding folders. Since an object inherits all properties up its hierarchy, it means that all files under /files might have properties associated with the "pictures" schema.

_Note: Schema creation and associating folders with schemas could be done via the UI as well._

```yaml
meta:
  /:
    - picture
  /animals:
    - mammals
  /cars:
    - auto
``` 

Once a folder is associated with the schema, file properties will show up in the property panel (they could also be edited directly if you have proper privileges). It is also possible to filter by property values in the object browser – this applies not only to files, but also to connections, queries, datasets, and other objects.

## Scripting

```python
#name: Detect Car Plate Numbers
#description: Detect car plate numbers
#language: python
#input: file file { condition: entity:domain == "auto" }
#output: bool hasNumbers { set: file.auto:hasTags }
#tags: demo, files, panel, ml, opencv
#condition: file.isFile && file.size < 1e6 && (file.name.endsWith("jpg") || file.name.endsWith("jpeg"))

import cv2

image = cv2.imread(file)
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
detector = cv2.CascadeClassifier(cv2.data.haarcascades + "haarcascade_russian_plate_number.xml")
hasNumbers = len(detector.detectMultiScale(gray, scaleFactor=1.3, minNeighbors=3, minSize=(100, 25))) != 0
```

Let's take a closer look at the following two lines:

```python
#input: file file { condition: file.schema == auto }
#output: bool hasNumbers { set: file.auto:hasTags }
```

The first one imposes a condition on the input parameter that the object (file) has to have the "auto" schema associated with it. This means that the info panel will only be shown when user clicks on images in the "cars" folder.

See also:

  * Community: [File Management](https://community.datagrok.ai/t/new-feature-file-share-browser/17/6)
  * [Files](https://datagrok.ai/help/access/connectors/files)
  * [Scripting](https://datagrok.ai/help/develop/scripting)