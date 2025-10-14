#name: Text Extractor
#description: Extracts text from a file
#language: python
#tags: nlp, panel, extractor
#input: file file [A file that contains text]
#input: string extension
#output: string text {semType: text} [Extracted text]
#condition: file.isFile && 0 < file.size && file.size < 1e6 && supportedExtExtract(file.name)

import cld3
import pycountry
import textract
import cleantext

params = {'filename': file, 'extension': extension}

if file.endswith(('gif', 'jpg', 'jpeg', 'png', 'tiff', 'tif')):
    text = textract.process(**params).decode().strip()
    # The default setting for OCR
    code = "eng"

    # Detect up to 3 most frequent languages in the text
    for lang in cld3.get_frequent_languages(text, 3):
        # Look up the name
        langcode = lang.language
        if len(langcode) == 2:
            language = pycountry.languages.get(alpha_2=langcode).name
        elif len(langcode) == 3:
            language = pycountry.languages.get(alpha_3=langcode).name
        else:
            # No language detected
            break

        # Update code for multiple languages scenario, e.g., "eng+fra".
        # Tesseract might subsequently use only some of the given codes
        # instead of the full list, which helps to extract non-English texts,
        # i.e., the default "eng" will be ignored for texts written in French.
        code += "+" + pycountry.languages.get(name=language).alpha_3

    params.update({'method': 'tesseract', 'language': code})

# Extract text
text = textract.process(**params).decode()

text = cleantext.clean(text, lower=False, to_ascii=False, no_line_breaks=True)
