#name: LMER
#description: Mixed Effects Model (LMER)
#reference: https://en.wikipedia.org/wiki/Mixed_model
#language: r
#tags: demo, analysis, prediction, lda
#sample: winequality.csv
#input: dataframe table [Input data set table]
#input: column_list features {type:numerical;allowNulls:false} [Input data set table columns]
#input: column random {type:numerical;allowNulls:false} [Name of column that will be interped as random effects]
#input: column predict {type:numerical;allowNulls:false} [Name of column contains values to predict]
#input: double perc = 0.5 [Percentage of data that goes to training]
#output: dataframe predicted {action:join(table)} [Predicted column]
#output: graphics ldaModel

require(lme4)

table = table[, c(predict, random, unlist(features))]
features = features[features != predict]
train <- sample(1:nrow(table), floor(nrow(table) * perc))
formula <- as.formula(paste(predict, ' ~ ', paste(features, collapse=' * '), " + (1|", random, ')'))
model <- lmer(formula, data=table, subset=train, REML=FALSE)
predicted <- data.frame(predict(model, newdata=table, allow.new.levels=TRUE))
colnames(predicted) <- c("predicted")

ldaModel <- plot(model)
print(ldaModel)
