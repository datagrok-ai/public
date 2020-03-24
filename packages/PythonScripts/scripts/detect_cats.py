#name: Detect Cats
#description: Detects cats on image
#language: python
#input: file filePath
#output: bool hasCats
#tag: panel, widget
#condition: x.size < 1e6 && x.ext == "jpg"

import cv2

image = cv2.imread(filePath)
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
detector = cv2.CascadeClassifier(cv2.data.haarcascades + "haarcascade_frontalcatface.xml")
hasCats = len(detector.detectMultiScale(gray, scaleFactor=1.3, minNeighbors=3, minSize=(75, 75))) != 0
