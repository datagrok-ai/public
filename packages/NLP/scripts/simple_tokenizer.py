#name: Tokenizer
#description: Segment a text into tokens (words, numbers, punctuation)
#language: python
#input: string text [A text to tokenize]
#output: string tokens [Comma-separated tokens in the order in which they appear in the text]
#tags: nlp

import re


TOKEN_REGEX = re.compile(r'[a-z]+|[+-]?\d*[.,]?\d+|\S', re.I)

def tokenize(text, regex=TOKEN_REGEX):
    return regex.findall(text)

tokens = ','.join(tokenize(text))
