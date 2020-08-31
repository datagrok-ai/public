#name: Cell Imaging Segmentation
#description: Cell imaging segmentation based on Watershed algorithm
#language: python
#input: file file
#output: int cells [Number of detected cells segments]
#output: graphics segmented
#tags: demo, panel, files, opencv
#condition: file.isFile && file.size < 1e6 && file.path.contains("/cells/") && (file.name.endsWith("jpg") || file.name.endsWith("jpeg") || file.name.endsWith("png"))
#help-url: https://docs.opencv.org/3.4/d3/db4/tutorial_py_watershed.html

import cv2 as cv
import numpy as np
import matplotlib.pyplot as plt

# Gray color image
img = cv.imread(file)
img = cv.cvtColor(img, cv.COLOR_RGB2BGR)
orig_img = img.copy()
gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)
_, thresh = cv.threshold(gray, 0, 255, cv.THRESH_BINARY_INV + cv.THRESH_OTSU)

# Noise removal
kernel = np.ones((3, 3), dtype=np.uint8)
opening = cv.morphologyEx(thresh, cv.MORPH_OPEN, kernel, iterations=2)

# Sure background and foreground areas
sure_bg = cv.dilate(opening, kernel, iterations=3)
dist_transform = cv.distanceTransform(opening, cv.DIST_L2, 5)
_, sure_fg = cv.threshold(dist_transform, np.median(dist_transform), 255, 0)

# Finding unknown region
sure_fg = np.uint8(sure_fg)
unknown = cv.subtract(sure_bg, sure_fg)

# Marker segments
_, markers = cv.connectedComponents(sure_fg)
markers = markers + 1
markers[unknown == 255] = 0
markers = cv.watershed(img, markers)
img[markers == -1] = [255, 0, 255]

# Count cells segments
gray = cv.cvtColor(gray, cv.COLOR_GRAY2BGR)
gray[markers == -1] = [255, 0, 255]
th = cv.inRange(gray, (255, 0, 255), (255, 0, 255)).astype(np.uint8)
th = cv.bitwise_not(th)
_, _, stats, _ = cv.connectedComponentsWithStats(th, connectivity=4)
areas = stats[:, cv.CC_STAT_AREA]
th = np.median(areas) * 1.6
cells = np.sum(areas < th)

# Plot beautified image
markers[0, :] = 0; markers[:, 0] = 0; markers[-1, :] = 0; markers[:, -1] = 0
orig_img[markers == -1] = [255, 0, 255]
_, m2 = cv.threshold(markers.astype(np.uint8), 0, 255, cv.THRESH_BINARY | cv.THRESH_OTSU)
contours, hierarchy = cv.findContours(m2, cv.RETR_LIST, cv.CHAIN_APPROX_NONE)
for c in contours:
    cv.drawContours(orig_img, c, -1, (255, 0, 255), 3)
segmented = plt.imshow(orig_img)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.axis('off')
plt.show()
