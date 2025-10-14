#name: Text Statistics
#description: Check how easy your text is to understand
#language: python
#tags: nlp, panel
#input: file file {semType: text} [Given text]
#output: string language {semType: lang} [Detected language]
#output: string words [Number of words]
#output: string sentences [Number of sentences]
#output: string long_words [Number of long words (7+ letters)]
#output: string readability_test [FRES test for English and LIX index for other languages]
#output: double score [Readability score]
#output: string edu_level [U.S. school grade level]
#output: string comment [Comment on the complexity of a text]
#condition: file.isFile && file.size < 1e6 && supportedExt(file.name)
#reference: https://en.wikipedia.org/wiki/Flesch%E2%80%93Kincaid_readability_tests#Flesch_reading_ease
#reference: https://en.wikipedia.org/wiki/Lix_(readability_test)

import cld3
import pycountry
import nltk
from nltk.corpus import cmudict
from nltk.tokenize import sent_tokenize, word_tokenize
import textract


nltk.download('cmudict')
nltk.download('punkt')
nltk.download('punkt_tab')

# Pronunciation dictionary
PROND = cmudict.dict()
# Lowercase vowel letters
VOWELS = "aeiou"


def count_syllables(word):
    """
    Counts the number of syllables as the total of pronounced vowels per word.
    In Carnegie Mellon University Pronouncing Dictionary, vowels carry a special
    stress marker ranging from 0 to 2, which allows for their identification.
    For non-dictionary words, it is simply the count of vowel letters per word.
    """
    if word in PROND:
        return len([phoneme for phoneme in PROND[word][0] if phoneme[-1].isdigit()])
    return len([char for char in word if char in VOWELS])


def fres(n_syllables, n_words, n_sentences):
    """
    The Flesch Reading Ease (FRES) scores a text based on the average number of
    words per sentence and the average number of syllables per word.
    """
    return 206.835 - 1.015 * n_words / n_sentences - 84.6 * n_syllables / n_words


def lix(n_words, n_sentences, n_long_words):
    """
    Läsbarhetsindex (LIX) measures readability given the average number of words
    per sentence and the percentage of long words (7+ letters) in the text.
    """
    return n_words / n_sentences + 100 * n_long_words / n_words


def detect_language(text):
    """Detects the language of a text and returns its name."""
    detected = cld3.get_language(text)
    if detected:
        langcode = detected.language
        if len(langcode) == 2:
            return pycountry.languages.get(alpha_2=langcode).name
        elif len(langcode) == 3:
            return pycountry.languages.get(alpha_3=langcode).name
    return "Undetermined"


def get_counts(text, language):
    """
    Preprocesses text and returns the number of syllables, (long) words, and
    sentences. If language is not supported, uses the default setting (English)
    """
    try:
        sentences = sent_tokenize(text, language=language.lower())
        words = [word.lower() for word in word_tokenize(text,
                 language=language.lower()) if word.isalpha()]
    except LookupError:
        sentences = sent_tokenize(text)
        words = [word.lower() for word in word_tokenize(text) if word.isalpha()]
    long_words = [word for word in words if len(word) > 6]
    n_syllables = sum(count_syllables(word) for word in words)
    n_words = len(words)
    n_sentences = len(sentences)
    n_long_words = len(long_words)
    return n_syllables, n_words, n_sentences, n_long_words


def comment(score, test):
    """Returns a comment on the complexity of text based on FRES and LIX."""
    if (score > 90 and test == 'FRES') or (score < 25 and test == 'LIX'):
        return ("Very easy to read", "5th grade")
    elif (score > 80 and test == 'FRES') or (score < 30 and test == 'LIX'):
        return ("Easy to read", "6th grade")
    elif (score > 70 and test == 'FRES') or (score < 35 and test == 'LIX'):
        return ("Fairly easy to read", "7th grade")
    elif (score > 60 and test == 'FRES') or (score < 45 and test == 'LIX'):
        return ("Plain English", "8th to 9th grade")
    elif (score > 50 and test == 'FRES') or (score < 55 and test == 'LIX'):
        return ("Fairly difficult to read", "10th to 12th grade")
    elif (score > 30 and test == 'FRES') or (score < 60 and test == 'LIX'):
        return ("Difficult to read", "College")
    return ("Very difficult to read", "College graduate")


# Extract text
params = {'filename': file, 'extension': file[file.rfind('.', 0, -10) : -10]}
text = textract.process(**params).decode().strip()

# Detect the language
language = detect_language(text)

# Calculate counts
syllables, words, sentences, long_words = get_counts(text, language)

# Compute the score
if language == 'English':
    readability_test = 'Flesch Reading-Ease'
    score = round(fres(syllables, words, sentences))
else:
    readability_test = 'LIX (Läsbarhetsindex)'
    score = round(lix(words, sentences, long_words))

# Comment on complexity
comment, edu_level = comment(score, readability_test)
