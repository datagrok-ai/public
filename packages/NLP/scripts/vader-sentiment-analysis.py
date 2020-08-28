#name: Rule-Based Sentiment Analysis
#description: Detect polar sentiments in a text with the VADER lexicon
#language: python
#input: string text {semType: text} [Statement where sentiments are expressed]
#output: string sentiment [Detected sentiment]
#output: double polarity [Polarity of expression on the scale from -1 to 1]
#tags: nlp, panel

import nltk
from nltk.sentiment.vader import SentimentIntensityAnalyzer
from nltk.tokenize import sent_tokenize


# Download sentence tokenizer and sentiment lexicon
nltk.download('punkt')
nltk.download('vader_lexicon')

# Segment the text into sentences
sentences = sent_tokenize(text)

# Set up the text's average sentiment
polarity = 0.0

# Analyze the text at the sentence level
sent_analyzer = SentimentIntensityAnalyzer()
for sentence in sentences:
    scores = sent_analyzer.polarity_scores(sentence)
    polarity += scores["compound"]

# Normalize the polarity
polarity = round(polarity / len(sentences), 4)

if polarity >= 0.05:
    sentiment = "POSITIVE"
elif -0.05 < polarity < 0.05:
    sentiment = "NEUTRAL"
else:
    sentiment = "NEGATIVE"
