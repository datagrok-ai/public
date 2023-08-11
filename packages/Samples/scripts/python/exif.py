#name: EXIF
#description: Gets EXIF data from PNG, JPEG, WEBP
#language: python
#input: file file
#output: map exif
#tags: demo, files, panel, images, metadata-extractor
#condition: file.isFile && file.size < 1e6 && (file.name.endsWith("jpg") || file.name.endsWith("jpeg") || file.name.endsWith("png") || file.name.endsWith("webp"))

from PIL import Image, ExifTags

img = Image.open(file)
info = img.getexif()
exif = {}
if info:
    info = dict(info)
    for key in info.keys():
        if key in ExifTags.TAGS:
            try:
                json.dumps({key: info[key]})
                exif[ExifTags.TAGS[key]] = info[key]
            except:
                pass
else:
    info = img.info
    for key in info.keys():
        try:
            json.dumps({key: info[key]})
            exif[key] = info[key]
        except:
            pass
