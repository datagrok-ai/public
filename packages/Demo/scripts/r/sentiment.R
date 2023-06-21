#name: Sentiment Classification
#description: Sentiment classification of emotions and polarity
#reference: https://en.wikipedia.org/wiki/Sentiment_analysis
#language: r
#tags: demo
#sample: tweets.csv
#input: dataframe data [Input data table]
#input: column col {type:text} [Name of column contains text data]
#output: dataframe emotion {action:join(data)} [Column with distribution of emotions]
#output: dataframe polarity {action:join(data)} [Column with distribution of polarity]
#output: graphics plotEmo [Plot with distribution of emotions]
#output: graphics plotPol [Plot with distribution of polarity]

require(sentiment)
require(ggplot2)
require(RColorBrewer)

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

# Plot distribution of emotions
plotEmo =
  ggplot(emotion, aes(x=emotion)) +
    geom_bar(aes(y=..count.., fill=emotion)) +
    scale_fill_brewer(palette="Dark2") +
    labs(x="Emotion categories", y="Number of text lines") +
    ggtitle("Sentiment analysis\n(by emotion)") +
    theme(title=element_text(size=11),
          plot.title=element_text(size=11))

# Plot distribution of polarity
plotPol =
  ggplot(polarity, aes(x=polarity)) +
    geom_bar(aes(y=..count.., fill=polarity)) +
    scale_fill_brewer(palette="RdGy") +
    labs(x="Polarity categories", y="Number of text lines") +
    ggtitle("Sentiment analysis\n(by polarity)") +
    theme(title=element_text(size=11),
          plot.title=element_text(size=11))

print(plotEmo)
print(plotPol)
