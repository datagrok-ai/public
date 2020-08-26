#name: Language Detection
#description: Detect the language in which a text is written
#language: python
#input: string text {semType: text} [A text to analyze]
#output: string language {semType: lang} [Detected language]
#tags: nlp, panel

import cld3
import pycountry


language = "Undetermined"

# Detect language
detected = cld3.get_language(text)

# Look up its name
if detected:
    langcode = detected.language
    if len(langcode) == 2:
        language = pycountry.languages.get(alpha_2=langcode).name
    elif len(langcode) == 3:
        language = pycountry.languages.get(alpha_3=langcode).name
