#name: Readability
#description: Check how easy your text is to understand
#language: python
#tags: nlp, panel
#input: string text {semType: text} [Given text]
#output: double FRES [Flesch Reading Ease Score]
#output: string comment [Comment on the complexity of a text]

import nltk
from nltk.corpus import cmudict
from nltk.tokenize import sent_tokenize, word_tokenize


nltk.download('cmudict')
nltk.download('punkt')

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


def get_counts(text):
    """Preprocesses text and returns the number of syllables, words, and sentences."""
    sentences = sent_tokenize(text)
    words = [word.lower() for word in word_tokenize(text) if word.isalpha()]
    n_syllables = sum(count_syllables(word) for word in words)
    n_words = len(words)
    n_sentences = len(sentences)
    return n_syllables, n_words, n_sentences


def comment(score):
    """Returns a comment on the complexity of text based on FRES."""
    if score > 90:
        return "Very easy to read"
    elif score > 80:
        return "Easy to read"
    elif score > 70:
        return "Fairly easy to read"
    elif score > 60:
        return "Plain English"
    elif score > 50:
        return "Fairly difficult to read"
    elif score > 30:
        return "Difficult to read"
    return "Very difficult to read"


# Compute the score
FRES = round(fres(*get_counts(text)))
# Comment on complexity
comment = comment(FRES)
