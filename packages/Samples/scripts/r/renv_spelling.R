#name: renvSpellingExample
#description: Simple Spell Checking
#language: r
#output: dataframe packagesBefore
#output: dataframe packagesAfter
#output: bool isCorrect1
#output: bool isCorrect2
#output: bool isCorrect3

packagesBefore = as.data.frame(installed.packages()[,c(1,3:4)])
renv::init()
renv::install("hunspell@3.0.1")
packagesAfter = as.data.frame(installed.packages()[,c(1,3:4)])
library(hunspell)
words <- c("deployment", "servise", "analysis")
correct <- hunspell_check(words)
isCorrect1 <- correct[1]
isCorrect2 <- correct[2]
isCorrect3 <- correct[3]
