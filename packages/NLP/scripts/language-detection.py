#name: Language Detector
#description: Detect the language in which a text is written
#language: python
#input: string text {semType: text} [A text to analyze]
#output: string language {semType: lang} [Detected language]
#output: string alpha_2 {semType: langCode} [ISO 639-3 language code]
#output: string alpha_3 {semType: langCode} [ISO 639-3 language code]
#meta.role: panel
#reference: https://github.com/google/cld3

import cld3
import pycountry

def detect_language():
    detected = cld3.get_language(text)
    undetermined = ("Undetermined", "un", "und")
    if not detected: return undetermined
    langcode = detected.language
    if langcode == "iw": langcode = "he"
    if len(langcode) == 2:
        lang = pycountry.languages.get(alpha_2=langcode)
    elif len(langcode) == 3:
        lang = pycountry.languages.get(alpha_3=langcode)
    if not lang or not detected.is_reliable or detected.probability < 0.95:
        return undetermined
    return (lang.name, lang.alpha_2, lang.alpha_3)

language, alpha_2, alpha_3 = detect_language()
