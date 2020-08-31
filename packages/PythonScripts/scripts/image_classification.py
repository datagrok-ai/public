#name: Image Classification
#description: Image classification based of Efficientnet model
#language: python
#input: file file
#output: map classes [Detected classes with probabilities]
#tags: demo, panel, files, efficientnet
#condition: file.isFile && file.size < 1e6 && (file.name.endsWith("jpg") || file.name.endsWith("jpeg") || file.name.endsWith("png"))
#help-url: https://github.com/qubvel/efficientnet

import numpy as np
from skimage.io import imread
from efficientnet.keras import EfficientNetB0
from keras.applications.imagenet_utils import decode_predictions
from efficientnet.keras import center_crop_and_resize, preprocess_input

image = imread(file)
model = EfficientNetB0(weights='imagenet')
image_size = model.input_shape[1]
_image = center_crop_and_resize(image, image_size=image_size)
_image = preprocess_input(_image)
_image = np.expand_dims(_image, 0)
predicted = model.predict(_image)
predicted = decode_predictions(predicted)[0]

classes = {}
for p in predicted:
    classes[p[1]] = float(p[2])
