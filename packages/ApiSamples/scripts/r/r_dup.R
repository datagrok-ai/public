#name: RDup
#description: Duplicates a string in R
#language: r
#tags: test, selenium
#input: string s = "abc"
#output: string res
#test: RDup("abc") == "abcabc"

res <- paste(s, s, sep='')
