#name: R Env File String Reverse
#description: reverse a string using stringi from a yaml env file
#language: r
#environment: r_test_env
#meta.timeout: 600
#input: string input_string
#output: string reversed

library(stringi)
reversed <- stri_reverse(input_string)
