#name: R Env String Reverse
#description: reverse a string using stringi from a custom env
#language: r
#environment: channels: [conda-forge], dependencies: [r-stringi]
#input: string input_string
#output: string reversed

library(stringi)
reversed <- stri_reverse(input_string)
