#name: Detect Car Plate Numbers
#description: Detect car plate numbers
#language: python
#input: file file
#output: bool hasNumbers
#tags: demo, files, panel, ml, opencv
#condition: file.size < 1e6 && (file.name.endsWith("jpg") || file.name.endsWith("jpeg"))

import cv2

image = cv2.imread(file)
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
detector = cv2.CascadeClassifier(cv2.data.haarcascades + "haarcascade_russian_plate_number.xml")
hasNumbers = len(detector.detectMultiScale(gray, scaleFactor=1.3, minNeighbors=3, minSize=(100, 25))) != 0
