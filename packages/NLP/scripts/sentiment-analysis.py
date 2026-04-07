#name: Sentiment Analysis
#description: Detect polar sentiments in a text
#language: python
#input: string text {semType: text} [Statement where sentiments are expressed]
#output: string sentiment [Detected sentiment]
#output: double polarity [Polarity of expression]
#meta.role: panel
#reference: https://en.wikipedia.org/wiki/Sentiment_analysis

from transformers import pipeline


# Download a model for sentiment analysis
classifier = pipeline("sentiment-analysis")

# Predict a label for the given text
label = classifier(text)[0]

# Interpret the results
sentiment = label["label"]
polarity = label["score"]
