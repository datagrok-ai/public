#name: Image Pixel Count
#description: image pixel count
#language: python
#environment: my_test_env
#input: file fileInput
#output: int count

import cv2
import numpy as np

img = cv2.imread(fileInput)
count = np.sum(img)
