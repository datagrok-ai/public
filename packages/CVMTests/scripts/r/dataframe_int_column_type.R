#name: R Int Column
#language: r
#output: dataframe resultInBound
#output: dataframe resultOutBound


col1 <- as.integer(c(1000, 10000, 100000))
resultInBound <- data.frame(col1)
col1 <- c(1000, 10000, 2147483648)
resultOutBound <- data.frame(col1)