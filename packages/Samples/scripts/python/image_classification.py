#name: Image Classification
#description: Image classification based of Efficientnet model
#language: python
#input: file file
#output: map classes [Detected classes with probabilities]
#tags: demo, panel, files, efficientnet
#condition: file.isFile && file.size < 1e6 && (file.name.endsWith("jpg") || file.name.endsWith("jpeg") || file.name.endsWith("png"))
#help-url: https://keras.io/api/applications/efficientnet/

import numpy as np
from skimage.io import imread
from skimage.transform import resize
from keras.applications.efficientnet import EfficientNetB0, preprocess_input, decode_predictions

image = imread(file)
model = EfficientNetB0(weights='imagenet')
image_size = model.input_shape[1]
_image = resize(image, (image_size, image_size), preserve_range=True, anti_aliasing=True)
_image = preprocess_input(_image)
_image = np.expand_dims(_image, 0)
predicted = model.predict(_image)
predicted = decode_predictions(predicted)[0]

classes = {}
for p in predicted:
    classes[p[1]] = float(p[2])
