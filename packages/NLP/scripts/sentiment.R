#name: Sentiment Classification
#description: Sentiment classification of emotions and polarity
#reference: https://en.wikipedia.org/wiki/Sentiment_analysis
#language: r
#input: dataframe data [Input data table]
#input: column col {semType: text} [Name of a column that contains text data]
#output: dataframe emotion {action:join(data)} [Column with distribution of emotions]
#output: dataframe polarity {action:join(data)} [Column with distribution of polarity]


require(sentiment)

# Select input column
data = data[[col]]

# Classify emotion, get emotion best fit and substitute NA's by "unknown"
class_emo = classify_emotion(data, algorithm="bayes", prior=1.0)
emotion = class_emo[,7]
emotion[is.na(emotion)] = "unknown"
emotion = data.frame(emotion=emotion, stringsAsFactors=FALSE)

# Classify polarity and get polarity best fit
class.pol = classify_polarity(data, algorithm="bayes")
polarity = data.frame(polarity=class.pol[,4], stringsAsFactors=FALSE)
