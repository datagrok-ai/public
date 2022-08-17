#name: LDA
#description: Linear Discriminant Analysis
#reference: https://en.wikipedia.org/wiki/Linear_discriminant_analysis
#language: r
#tags: demo, analysis, prediction, lda
#sample: winequality.csv
#input: dataframe table [Input data set table]
#input: column predict {type:categorical;allowNulls:false} [Name of column contains values to predict]
#input: column_list features {type:numerical;allowNulls:false} [Input data set table columns]
#input: double perc = 0.5 [Percentage of data that goes to training]
#output: dataframe predicted {action:join(table)} [Predicted column]
#output: graphics ldaModel

require(MASS)

table = table[, c(predict, unlist(features))]
train <- sample(1:nrow(table), floor(nrow(table) * perc))
formula <- as.formula(paste(predict, " ~ ."))
model <- lda(formula, table, subset=train)
prediction <- predict(model, newdata=table)
predicted <- data.frame(prediction$class)

plot(model)
