#name: Language Detection
#description: Detect the language in which a text is written
#language: python
#input: file file {semType: text} [A text to analyze]
#output: string language {semType: lang} [Detected language]
#output: string alpha_2 {semType: langCode} [ISO 639-3 language code]
#output: string alpha_3 {semType: langCode} [ISO 639-3 language code]
#tags: nlp
#reference: https://github.com/google/cld3
#condition: file.isFile && file.size < 1e6 && supportedExt(file.name)

import cld3
import pycountry
import textract

params = {'filename': file, 'extension': file[file.rfind('.', 0, -10) : -10]}

# Extract text
text = textract.process(**params).decode().strip()

# with open(file) as f:
#     text = " ".join(f.read().splitlines())

language = "Undetermined"
alpha_2 = "—"
alpha_3 = "—"

# Detect language
detected = cld3.get_language(text)

# Look up its name
if detected:
    langcode = detected.language
    if len(langcode) == 2:
        lang = pycountry.languages.get(alpha_2=langcode)
    elif len(langcode) == 3:
        lang = pycountry.languages.get(alpha_3=langcode)
    language = lang.name
    alpha_2 = lang.alpha_2
    alpha_3 = lang.alpha_3
