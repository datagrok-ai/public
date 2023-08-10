#name: RDup
#description: Duplicates a string in R
#language: r
#tags: test, selenium
#input: string s
#output: string res
#test: RDup("abc") == "abcabc"

res <- paste(s, s, sep='')
