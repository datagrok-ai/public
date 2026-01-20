#name: Rule-Based Sentiment Analysis
#description: Detect polar sentiments in a text with the VADER lexicon
#language: python
#input: dataframe data [Table with text data]
#input: column col {type: string} [Name of a text column where sentiments are expressed]
#output: dataframe polarity {action: join(data)} [Column with distribution of polarity on the scale from -1 to 1]
#output: dataframe sentiment {action: join(data)} [Column with the polarity interpretation]
#sample: tweets.csv
#reference: https://github.com/cjhutto/vaderSentiment

import nltk
from nltk.sentiment.vader import SentimentIntensityAnalyzer
from nltk.tokenize import sent_tokenize


# Download sentence tokenizer and sentiment lexicon
nltk.download('punkt')
nltk.download('vader_lexicon')


def get_polarity_score(text):
    """Analyzes `text` at the sentence level and returns its overall polarity."""
    sentences = sent_tokenize(text)
    if text and not sentences:
        sentences = [text]
    polarity = 0.0
    sent_analyzer = SentimentIntensityAnalyzer()
    for sentence in sentences:
        scores = sent_analyzer.polarity_scores(sentence)
        polarity += scores["compound"]
    return round(polarity / len(sentences), 4)

def get_sentiment(polarity_score):
    """Returns a string representation for `polarity_score`"""
    if polarity_score >= 0.05:
        return "POSITIVE"
    elif -0.05 < polarity_score < 0.05:
        return "NEUTRAL"
    return "NEGATIVE"

# Create output dataframes
polarity = data[col].apply(get_polarity_score).to_frame(name='polarity')
sentiment = polarity['polarity'].apply(get_sentiment).to_frame(name='sentiment')
