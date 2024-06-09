#name: R Empty Dataframe
#description: df input/output
#language: r
#output: dataframe resultDf

resultDf <- data.frame(Date=as.Date(character()),
                             File=character(),
                             User=character(),
                             stringsAsFactors=FALSE)
