#name: Sentiment Analysis
#description: Detect polar sentiments in a text
#language: python
#input: string text {semType: text} [Statement where sentiments are expressed]
#output: string sentiment [Detected sentiment]
#output: double polarity [Polarity of expression]
#tags: nlp, panel
#environment: flair_sentiment

from flair.data import Sentence
from flair.models import TextClassifier


# Download a model for sentiment analysis
classifier = TextClassifier.load('en-sentiment')

# Predict a label for the given text
sentence = Sentence(text)
classifier.predict(sentence)
label = sentence.labels[0]

# Interpret the results
sentiment = label.value
polarity = label.score
